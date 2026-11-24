# Licensed under the Apache License, Version 2.0.
"""Shared bulk RNA-seq ``(qc, cleanup, normalize)`` algebra, NNLS mixture, & Beta quantiles (#778, #780).

Provides dataset-agnostic primitives reused across :mod:`oncoref.gtex`, non-GTEx human
normal bulk panels, non-human species :class:`~oncoref.species_expression.ExpressionDataset`
atlases, and bulk tumor cohorts:

1. **Operation & alias expansion** (:func:`expand_ops`): maps ``None``, ``"all"`` (or other
   named aliases in ``aliases``), ``"+"``/``","``-delimited strings, or iterables of atomic
   operation names to a canonical ``frozenset[str]``.
2. **Thread/process-safe spec-keyed Parquet caching** (:func:`read_or_compute_cached_parquet`):
   deterministic SHA-256 cache entries with atomic ``os.replace`` under :class:`filelock.FileLock`.
3. **``NaN``-safe coordinate-descent NNLS tissue-archetype mixture solver**
   (:func:`solve_nonnegative_archetype_mixture`): pure ``numpy`` (no ``scipy``).
4. **Exact ``scipy``-free purity-weighted Harrell-Davis Beta quantiles**
   (:func:`purity_weighted_harrell_davis_quantile`): Lentz continued-fraction regularized
   incomplete beta ``I_x(a, b)`` smooth across all ``n_eff`` and tail quantiles ``q``.
"""

from __future__ import annotations

import contextlib
import json
import math
import os
import re
import tempfile
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path

import numpy as np
import pandas as pd
from filelock import FileLock


def expand_ops(
    value: str | bool | Iterable[str] | None,
    *,
    allowed: tuple[str, ...],
    aliases: dict[str, frozenset[str]],
    param_name: str,
    domain_label: str = "bulk",
) -> frozenset[str]:
    """Expand ``None``, ``"all"`` (or other aliases), or any list/set of operations to a ``frozenset``."""
    if value is None or value is False:
        return frozenset()
    if value is True:
        return aliases.get("all", frozenset(allowed))
    if isinstance(value, Mapping):
        raise TypeError(
            f"{param_name} must be None, an alias ({sorted(aliases)}), or a list/set of operations; "
            f"got dict/mapping"
        )
    if isinstance(value, str):
        raw_tokens = [p.strip().lower() for p in re.split(r"[+,]", value) if p.strip()]
        if not raw_tokens:
            return frozenset()
    elif isinstance(value, Iterable):
        raw_tokens = []
        for item in value:
            if not isinstance(item, str):
                raise TypeError(
                    f"{param_name} elements must be operation strings; got {type(item).__name__}"
                )
            for p in re.split(r"[+,]", item):
                tok = p.strip().lower()
                if tok:
                    raw_tokens.append(tok)
    else:
        raise TypeError(
            f"{param_name} must be None, an alias ({sorted(aliases)}), or an iterable of operations; "
            f"got {type(value).__name__}"
        )

    expanded: set[str] = set()
    unknown: list[str] = []
    allowed_set = set(allowed)
    for tok in raw_tokens:
        if tok in aliases:
            expanded.update(aliases[tok])
        elif tok in allowed_set:
            expanded.add(tok)
        else:
            unknown.append(tok)
    if unknown:
        raise ValueError(
            f"Unknown {domain_label} {param_name} step(s) {sorted(set(unknown))}; "
            f"choose from {allowed}, aliases {sorted(aliases)}, or None"
        )
    return frozenset(expanded)


def read_or_compute_cached_parquet(
    cache_dir: Path,
    key: str,
    manifest_payload: dict,
    builder_fn: Callable[[], pd.DataFrame],
    *,
    use_cache: bool = True,
    spec_attr_name: str = "spec_manifest",
) -> pd.DataFrame:
    """Read or atomically materialize a spec-keyed Parquet DataFrame under ``cache_dir``."""
    if not use_cache:
        df = builder_fn()
        df.attrs[spec_attr_name] = manifest_payload
        return df

    pq_path = cache_dir / f"{key}.parquet"
    meta_path = cache_dir / f"{key}.json"
    lock_path = cache_dir / f".{key}.lock"

    def _try_read_existing() -> pd.DataFrame | None:
        if pq_path.is_file() and meta_path.is_file():
            with contextlib.suppress(Exception):
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                if meta.get("spec_key") == key:
                    df_cached = pd.read_parquet(pq_path)
                    df_cached.attrs[spec_attr_name] = meta
                    df_cached.attrs["cache_path"] = str(pq_path)
                    return df_cached
        return None

    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        lock_ctx = FileLock(str(lock_path))
        lock_ctx.acquire()
    except OSError:
        existing = _try_read_existing()
        if existing is not None:
            return existing
        df = builder_fn()
        df.attrs[spec_attr_name] = manifest_payload
        return df

    try:
        existing = _try_read_existing()
        if existing is not None:
            return existing

        df = builder_fn()
        with contextlib.suppress(Exception):
            fd, tmp_pq = tempfile.mkstemp(dir=cache_dir, prefix=f".{key}-", suffix=".parquet.part")
            os.close(fd)
            try:
                df.to_parquet(tmp_pq, index=False)
                os.replace(tmp_pq, pq_path)
                meta_out = {**manifest_payload, "rows": len(df), "columns": list(df.columns)}
                meta_path.write_text(
                    json.dumps(meta_out, indent=2, sort_keys=True) + "\n", encoding="utf-8"
                )
            finally:
                Path(tmp_pq).unlink(missing_ok=True)
        df.attrs[spec_attr_name] = manifest_payload
        if pq_path.is_file():
            df.attrs["cache_path"] = str(pq_path)
        return df
    finally:
        with contextlib.suppress(Exception):
            lock_ctx.release()


def solve_nonnegative_archetype_mixture(
    basis_tpm: pd.DataFrame,
    sample_tpm: pd.DataFrame,
    *,
    l1_penalty: float = 1e-3,
    max_iter: int = 80,
    tol: float = 1e-6,
) -> pd.DataFrame:
    """Solve joint non-negative least squares ``Y_tilde ~= A_tilde @ W`` (``W >= 0``).

    Row-normalizes each marker gene ``g`` by ``max(max_k A_{g,k}, 10.0)`` so ultra-abundant
    secretory transcripts and moderate parenchymal transcription factors contribute on a
    balanced ``[0, 1]`` scale, then solves the L1-regularized convex NNLS problem across all
    samples simultaneously via vectorized coordinate descent in pure ``numpy``. ``NaN``/``inf``
    observations are sanitized to ``0.0`` so missing values never poison ``gram`` or ``rhs``.
    """
    common = basis_tpm.index.intersection(sample_tpm.index)
    arch_names = list(basis_tpm.columns)
    samp_names = list(sample_tpm.columns)
    if common.empty or not arch_names or not samp_names:
        return pd.DataFrame(0.0, index=arch_names, columns=samp_names)

    a_arr = np.nan_to_num(
        basis_tpm.loc[common].to_numpy(dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0
    )
    y_arr = np.nan_to_num(
        sample_tpm.loc[common].to_numpy(dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0
    )
    a_raw = np.fmax(a_arr, 0.0)
    y_raw = np.fmax(y_arr, 0.0)
    row_scale = np.maximum(a_raw.max(axis=1, keepdims=True), 10.0)
    a_til = a_raw / row_scale
    y_til = y_raw / row_scale

    gram = a_til.T @ a_til  # K x K
    rhs = a_til.T @ y_til  # K x N
    diag = np.diag(gram).copy()
    k_dim = len(arch_names)
    w = np.zeros((k_dim, len(samp_names)), dtype=np.float64)
    active_ks = [k for k in range(k_dim) if diag[k] > 1e-12]
    if not active_ks:
        return pd.DataFrame(w, index=arch_names, columns=samp_names)

    for _ in range(max_iter):
        max_delta = 0.0
        for k in active_ks:
            old_wk = w[k].copy()
            numer = rhs[k] - (gram[k] @ w - diag[k] * old_wk) - l1_penalty
            new_wk = np.fmax(0.0, numer / diag[k])
            w[k] = new_wk
            diff = float(np.max(np.abs(new_wk - old_wk)))
            if diff > max_delta:
                max_delta = diff
        if max_delta < tol:
            break

    return pd.DataFrame(w, index=arch_names, columns=samp_names)


def _betacf_scalar(a: float, b: float, x: float) -> float:
    """Modified Lentz continued fraction for the incomplete beta function."""
    fpmin = 1e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, 201):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-14:
            break
    return h


def regularized_betainc(a: float, b: float, x_arr: np.ndarray) -> np.ndarray:
    """Exact scipy-free regularized incomplete beta ``I_x(a, b)`` over ``x_arr in [0, 1]``."""
    out = np.empty_like(x_arr, dtype=np.float64)
    log_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    sym_pivot = (a + 1.0) / (a + b + 2.0)
    for idx, x_raw in enumerate(x_arr):
        x = float(x_raw)
        if math.isnan(x):
            out[idx] = float("nan")
        elif x <= 0.0:
            out[idx] = 0.0
        elif x >= 1.0:
            out[idx] = 1.0
        else:
            log_bt = a * math.log(x) + b * math.log1p(-x) - log_beta
            bt = math.exp(log_bt) if log_bt > -740.0 else 0.0
            if bt == 0.0:
                out[idx] = 0.0 if x < sym_pivot else 1.0
            elif x < sym_pivot:
                out[idx] = min(1.0, max(0.0, bt * _betacf_scalar(a, b, x) / a))
            else:
                out[idx] = min(1.0, max(0.0, 1.0 - bt * _betacf_scalar(b, a, 1.0 - x) / b))
    return out


def purity_weighted_harrell_davis_quantile(
    values: np.ndarray | pd.Series,
    weights: np.ndarray | pd.Series | None = None,
    quantile: float = 0.95,
) -> float:
    """Compute exact purity-weighted Harrell-Davis Beta quantile without ``scipy``.

    Evaluates ``w_i = I_{U_i}(a, b) - I_{U_{i-1}}(a, b)`` with ``a = q*(n_eff + 1)`` and
    ``b = (1 - q)*(n_eff + 1)`` across cumulative purity-weight intervals ``[U_{i-1}, U_i]``.
    """
    if not (0.0 < quantile < 1.0):
        raise ValueError(f"quantile must be in (0, 1), got {quantile!r}")
    x = np.asarray(values, dtype=np.float64)
    w = np.ones_like(x) if weights is None else np.asarray(weights, dtype=np.float64)
    valid = np.isfinite(x) & np.isfinite(w) & (w > 0.0)
    x = x[valid]
    w = w[valid]
    if x.size == 0:
        return float("nan")
    if x.size == 1:
        return float(x[0])
    order = np.argsort(x, kind="mergesort")
    x = x[order]
    w = w[order]
    w_norm = w / w.sum()
    n_eff = float((w.sum() ** 2) / np.dot(w, w))
    a = quantile * (n_eff + 1.0)
    b = (1.0 - quantile) * (n_eff + 1.0)
    edges = np.concatenate(([0.0], np.cumsum(w_norm)))
    edges[-1] = 1.0
    cdf = regularized_betainc(a, b, edges)
    hd_weights = np.fmax(0.0, np.diff(cdf))
    hd_sum = float(hd_weights.sum())
    if hd_sum <= 0.0 or not math.isfinite(hd_sum):
        cum = np.cumsum(w_norm) - 0.5 * w_norm
        return float(np.interp(quantile, cum, x))
    hd_weights /= hd_sum
    return float(np.dot(hd_weights, x))


__all__ = [
    "expand_ops",
    "purity_weighted_harrell_davis_quantile",
    "read_or_compute_cached_parquet",
    "regularized_betainc",
    "solve_nonnegative_archetype_mixture",
]

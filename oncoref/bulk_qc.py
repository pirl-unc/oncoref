# Licensed under the Apache License, Version 2.0.
"""Shared bulk RNA-seq ``(qc, cleanup, normalize)`` algebra, NNLS mixture, & Beta quantiles (#778, #780).

Provides dataset-agnostic primitives reused across :mod:`oncoref.gtex`, non-GTEx human
normal bulk panels (:mod:`oncoref.human_protein_atlas`), non-human species
:class:`~oncoref.species_expression.ExpressionDataset` atlases, and bulk tumor cohorts
(:mod:`oncoref.source_matrices` / :mod:`oncoref.library_qc`):

1. **Operation & alias expansion** (:func:`expand_ops`, :func:`resolve_normalize_option`,
   :class:`BulkProcessingSpec`): maps ``None``, ``"all"`` (or other named aliases in
   ``aliases``), ``"+"``/``","``-delimited strings, or iterables of atomic operation names to a
   canonical ``frozenset[str]`` with identical equality, ``hash()``, and SHA-256 Parquet cache keys.
2. **Thread/process-safe spec-keyed Parquet caching** (:func:`read_or_compute_cached_parquet`):
   deterministic SHA-256 cache entries with atomic ``os.replace`` under :class:`filelock.FileLock`.
3. **Stoichiometric splash masking, lane-bleed subtraction, & NNLS mixture adjustment**
   (:func:`estimate_foreign_panel_fractions`, :func:`mask_panel_splash_cells`,
   :func:`subtract_panel_index_hopping`, :func:`deconvolve_panel_admixture`,
   :func:`solve_nonnegative_archetype_mixture`): pure ``numpy``/``pandas`` (no ``scipy``).
4. **Exact ``scipy``-free purity-weighted Harrell-Davis Beta quantiles**
   (:func:`purity_weighted_harrell_davis_quantile`): Lentz continued-fraction regularized
   incomplete beta ``I_x(a, b)`` smooth across all ``n_eff`` and tail quantiles ``q``.
5. **Species-native 16/9/75 compartment scale renormalization**
   (:func:`species_compartment_clean_scale`): applies conserved ribosomal/technical/biological
   budget pinning without cross-species gene-ID canonicalization.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import math
import os
import re
import tempfile
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from filelock import FileLock

NORMALIZE_ALIASES: dict[str | None, str | None] = {
    None: None,
    "none": None,
    "raw": None,
    "native": None,
    "clean_tpm": "clean_tpm",
    "tpm_clean": "clean_tpm",
}

#: Canonical mapping from GTEx anchor ``source_site`` to somatic/germline marker panel names.
GTEX_ANCHOR_SITE_TO_PANEL: dict[str, str] = {
    "adrenal_gland": "adrenal_gland",
    "brain_cerebellar_hemisphere": "brain",
    "brain_cerebellum": "brain",
    "brain_nucleus_accumbens_basal_ganglia": "brain",
    "brain_spinal_cord_cervical_c-1": "brain",
    "colon_transverse": "stomach",
    "esophagus_mucosa": "squamous_epithelium",
    "heart_atrial_appendage": "striated_muscle",
    "heart_left_ventricle": "striated_muscle",
    "kidney_cortex": "kidney",
    "kidney_medulla": "kidney",
    "liver": "liver",
    "lung": "lung",
    "minor_salivary_gland": "salivary_gland",
    "muscle_skeletal": "striated_muscle",
    "ovary": "ovary",
    "pancreas": "pancreas",
    "pituitary": "pituitary",
    "prostate": "prostate",
    "skin_not_sun_exposed_suprapubic": "squamous_epithelium",
    "skin_sun_exposed_lower_leg": "squamous_epithelium",
    "small_intestine_terminal_ileum": "stomach",
    "stomach": "stomach",
    "testis": "testis",
    "thyroid": "thyroid",
}

REPRODUCTIVE_ORGAN_PANELS: frozenset[str] = frozenset(
    {
        "testis",
        "ovary",
        "placenta",
        "embryo",
        "epididymis",
        "fallopian_tube",
        "uterus",
        "cervix",
        "vagina",
        "germ_cell",
        "trophoblast",
    }
)

#: Conserved mammalian/vertebrate 1:1 organ-lineage anchor symbols used to derive
#: within-dataset organ reference profiles without cross-species expression borrowing.
VERTEBRATE_ORGAN_ANCHOR_SYMBOLS: dict[str, tuple[str, ...]] = {
    "liver": ("ALB", "FGG", "FGA", "FGB", "APOA2", "CRP", "HP", "ORM1"),
    "pancreas": ("INS", "INS1", "INS2", "PRSS1", "PRSS2", "CPA1", "CLPS", "CELA3A"),
    "lung": ("SFTPC", "SFTPA1", "SFTPA2", "SFTPB", "SCGB1A1", "AGER"),
    "kidney": ("UMOD", "SLC12A1", "AQP2", "SLC22A8", "NPHS2"),
    "thyroid": ("TG", "TPO", "IYD", "TSHR"),
    "adrenal_gland": ("CYP11B1", "CYP11B2", "CYP21A2", "STAR"),
    "pituitary": ("GH1", "PRL", "POMC", "TSHB", "LHB", "CGA"),
    "brain": ("GFAP", "MBP", "MOBP", "MOG", "GABRA6", "SNAP25", "SYP"),
    "striated_muscle": ("CKM", "ACTA1", "MB", "NPPA", "MYL7", "TNNC2", "MYH7"),
    "stomach": ("PGA5", "PGC", "GKN1", "ATP4B", "LIPF", "CBLIF"),
    "salivary_gland": ("MUC7", "STATH", "BPIFB2", "ZG16B"),
    "squamous_epithelium": ("KRT1", "KRT10", "KRT13", "KRT4", "SPRR3", "CRNN"),
    "testis": ("PRM1", "PRM2", "TNP1", "SMCP", "LELP1", "PHF7", "AKAP4"),
    "ovary": ("ZP3", "ZP2", "FIGLA", "GDF9", "BMP15", "FOXL2"),
    "placenta": ("PLAC1", "PSG1", "CYP19A1", "GCM1"),
}

#: Stoichiometric estimation anchors (TPM/nTPM scale) in pure human platelets and
#: peripheral reticulocytes for FACS/MACS immune-cell sort carryover estimation.
#: Avoids ``NRGN`` (endogenous in classical monocytes) and ``PF4``/``SLC4A1`` in the
#: q25 estimator (``PF4`` splits reads with paralog ``PF4V1`` in Kallisto; ``SLC4A1``
#: decays in enucleated reticulocytes), while subtracting the full panel in
#: :data:`SORTED_IMMUNE_CARRYOVER_PROFILES`.
SORTED_IMMUNE_CARRYOVER_PANELS: dict[str, dict[str, float]] = {
    "platelets": {
        "ENSG00000163736": 8500.0,  # PPBP
        "ENSG00000169704": 420.0,  # GP9
        "ENSG00000005961": 180.0,  # ITGA2B
        "ENSG00000101162": 350.0,  # TUBB1
    },
    "erythroid": {
        "ENSG00000244734": 45000.0,  # HBB
        "ENSG00000188536": 26000.0,  # HBA2
        "ENSG00000206172": 5500.0,  # HBA1
        "ENSG00000223609": 45.0,  # HBD
    },
}

SORTED_IMMUNE_CARRYOVER_PROFILES: dict[str, dict[str, float]] = {
    "platelets": {
        "ENSG00000163736": 8500.0,  # PPBP
        "ENSG00000163737": 6200.0,  # PF4
        "ENSG00000169704": 420.0,  # GP9
        "ENSG00000005961": 180.0,  # ITGA2B
        "ENSG00000101162": 350.0,  # TUBB1
        "ENSG00000185245": 120.0,  # GP1BA
    },
    "erythroid": {
        "ENSG00000244734": 45000.0,  # HBB
        "ENSG00000188536": 26000.0,  # HBA2
        "ENSG00000206172": 28000.0,  # HBA1
        "ENSG00000223609": 220.0,  # HBD
        "ENSG00000158578": 850.0,  # ALAS2
        "ENSG00000004939": 60.0,  # SLC4A1
    },
}

_RIBOSOMAL_SYMBOL_RE = re.compile(r"^(?:RPL|RPS)\d+[A-Z0-9]*$", re.IGNORECASE)
_TECHNICAL_SYMBOL_RE = re.compile(
    r"^(?:MT-|mt-|MTRNR|mt-Rnr|RNA5S|RNA18S|RNA28S|RN7SL|RN7SK|7SK)", re.UNICODE
)


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


def resolve_normalize_option(
    value: str | None,
    *,
    param_name: str = "normalize",
    domain_label: str = "bulk",
) -> str | None:
    """Validate and canonicalize ``normalize=None | 'clean_tpm'`` across bulk loaders."""
    if value is None:
        return None
    if not isinstance(value, str):
        raise TypeError(
            f"{domain_label} {param_name} must be None or 'clean_tpm'; got {type(value).__name__}"
        )
    key = value.strip().lower()
    if key not in NORMALIZE_ALIASES:
        raise ValueError(
            f"Unknown {domain_label} {param_name}={value!r}; choose None (native scale) or 'clean_tpm'"
        )
    return NORMALIZE_ALIASES[key]


@dataclass(frozen=True)
class BulkProcessingSpec:
    """Immutable specification of ``(qc, cleanup, normalize)`` for any bulk RNA dataset."""

    domain: str
    source_id: str
    source_version: str
    qc: frozenset[str] = frozenset()
    cleanup: frozenset[str] = frozenset()
    normalize: str | None = None
    extra: tuple[tuple[str, object], ...] = ()

    def __post_init__(self) -> None:
        if isinstance(self.qc, str):
            raise TypeError(
                "BulkProcessingSpec.qc must be an expanded frozenset/iterable, not bare str"
            )
        if isinstance(self.cleanup, str):
            raise TypeError(
                "BulkProcessingSpec.cleanup must be an expanded frozenset/iterable, not bare str"
            )
        object.__setattr__(self, "qc", frozenset(str(x).strip().lower() for x in self.qc))
        object.__setattr__(self, "cleanup", frozenset(str(x).strip().lower() for x in self.cleanup))
        object.__setattr__(
            self,
            "normalize",
            resolve_normalize_option(self.normalize, domain_label=self.domain),
        )
        object.__setattr__(self, "extra", tuple((str(k), v) for k, v in self.extra))

    @property
    def is_noop(self) -> bool:
        return not self.qc and not self.cleanup and self.normalize is None

    def _payload(self) -> dict[str, object]:
        return {
            "schema": "oncoref.bulk-processing-spec.v1",
            "domain": self.domain,
            "source_id": self.source_id,
            "source_version": self.source_version,
            "qc": sorted(self.qc),
            "cleanup": sorted(self.cleanup),
            "normalize": self.normalize,
            "extra": [[k, v] for k, v in self.extra],
        }

    def cache_key(self) -> str:
        raw = json.dumps(self._payload(), sort_keys=True, separators=(",", ":")).encode("utf-8")
        digest = hashlib.sha256(raw).hexdigest()[:20]
        safe_src = re.sub(
            r"[^A-Za-z0-9_.-]+", "_", f"{self.domain}-{self.source_id}-{self.source_version}"
        )
        return f"{safe_src}-{digest}"

    def manifest(self) -> dict[str, object]:
        key = self.cache_key()
        return {**self._payload(), "spec_key": key}


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
        prior = df.attrs.get(spec_attr_name)
        merged = {**(prior if isinstance(prior, dict) else {}), **manifest_payload}
        df.attrs[spec_attr_name] = merged
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
        prior = df.attrs.get(spec_attr_name)
        df.attrs[spec_attr_name] = {
            **(prior if isinstance(prior, dict) else {}),
            **manifest_payload,
        }
        return df

    try:
        existing = _try_read_existing()
        if existing is not None:
            return existing

        df = builder_fn()
        prior = df.attrs.get(spec_attr_name)
        merged_meta = {
            **(prior if isinstance(prior, dict) else {}),
            **manifest_payload,
            "rows": len(df),
            "columns": list(df.columns),
        }
        with contextlib.suppress(Exception):
            fd, tmp_pq = tempfile.mkstemp(dir=cache_dir, prefix=f".{key}-", suffix=".parquet.part")
            os.close(fd)
            try:
                df.to_parquet(tmp_pq, index=False)
                os.replace(tmp_pq, pq_path)
                meta_path.write_text(
                    json.dumps(merged_meta, indent=2, sort_keys=True) + "\n", encoding="utf-8"
                )
            finally:
                Path(tmp_pq).unlink(missing_ok=True)
        df.attrs[spec_attr_name] = merged_meta
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


def estimate_foreign_panel_fractions(
    values: pd.DataFrame,
    panels: Mapping[str, pd.Series],
    *,
    sample_native_panels: Mapping[str, frozenset[str]] | None = None,
    protected_panels: frozenset[str] = frozenset(),
    quantile: float = 0.25,
    min_present: int = 2,
) -> pd.DataFrame:
    """Estimate foreign-lineage bleed fraction ``rho_{p, s}`` per ``(panel, sample)``.

    Uses a conservative lower order statistic (``quantile=0.25`` across ``>= min_present``
    anchor genes) of ``observed / reference`` so single ectopic genes never trigger a
    stoichiometric multi-gene contamination call. Native panels for a sample and any
    ``protected_panels`` are identically zeroed.
    """
    cols = list(values.columns)
    panel_names = sorted(panels)
    if not cols or not panel_names or values.empty:
        return pd.DataFrame(0.0, index=panel_names, columns=cols)

    dedup = values.groupby(level=0).sum(min_count=1)
    out = np.zeros((len(panel_names), len(cols)), dtype=np.float64)
    for idx, pname in enumerate(panel_names):
        if pname in protected_panels:
            continue
        ref = panels[pname]
        ref = ref[ref > 0.0]
        present = ref.index.intersection(dedup.index)
        eff_min = min(min_present, len(ref))
        if len(present) < max(1, eff_min):
            continue
        sub = dedup.loc[present].to_numpy(dtype=np.float64)
        denom = ref.loc[present].to_numpy(dtype=np.float64)[:, None]
        ratios = np.fmax(0.0, np.nan_to_num(sub / denom, nan=0.0, posinf=0.0, neginf=0.0))
        if ratios.shape[0] == 1:
            est = ratios[0]
        else:
            est = np.quantile(ratios, quantile, axis=0)
        if sample_native_panels:
            for c_idx, cname in enumerate(cols):
                if pname in sample_native_panels.get(cname, frozenset()):
                    est[c_idx] = 0.0
        out[idx] = np.clip(est, 0.0, 1.0)
    return pd.DataFrame(out, index=panel_names, columns=cols)


def _build_profile_matrix(
    row_index: pd.Index,
    panels: Mapping[str, pd.Series],
    full_profiles: pd.DataFrame | None = None,
) -> pd.DataFrame:
    panel_names = sorted(panels)
    if full_profiles is not None and not full_profiles.empty:
        prof = full_profiles.reindex(index=row_index, columns=panel_names).fillna(0.0)
        return prof.clip(lower=0.0)
    data = {}
    for pname in panel_names:
        s = panels[pname].groupby(level=0).max()
        data[pname] = s.reindex(row_index).fillna(0.0).clip(lower=0.0)
    return pd.DataFrame(data, index=row_index, columns=panel_names)


def mask_panel_splash_cells(
    values: pd.DataFrame,
    panels: Mapping[str, pd.Series],
    *,
    sample_native_panels: Mapping[str, frozenset[str]] | None = None,
    protected_panels: frozenset[str] = frozenset(),
    full_profiles: pd.DataFrame | None = None,
    min_bleed_fraction: float = 0.005,
    min_bleed_value: float = 0.15,
    min_bleed_share: float = 0.20,
) -> pd.DataFrame:
    """``NaN``-mask ``(feature, sample)`` cells dominated by ``>= 0.5%`` foreign-panel splash.

    Surviving unmasked observations and ``0.0`` zeros remain byte-for-byte untouched.
    """
    out = values.copy()
    if out.empty or not panels:
        return out
    fracs = estimate_foreign_panel_fractions(
        out,
        panels,
        sample_native_panels=sample_native_panels,
        protected_panels=protected_panels,
    )
    active = fracs.where(fracs >= min_bleed_fraction, 0.0)
    if float(active.to_numpy().max(initial=0.0)) <= 0.0:
        return out
    prof = _build_profile_matrix(out.index, panels, full_profiles=full_profiles)
    bleed_est = prof.to_numpy(dtype=np.float64) @ active.to_numpy(dtype=np.float64)
    arr = out.to_numpy(dtype=np.float64, copy=True)
    mask = (
        np.isfinite(arr)
        & (arr > 0.0)
        & (bleed_est >= min_bleed_value)
        & (bleed_est >= min_bleed_share * arr)
    )
    if mask.any():
        arr[mask] = np.nan
        return pd.DataFrame(arr, index=out.index, columns=out.columns)
    return out


def subtract_panel_index_hopping(
    values: pd.DataFrame,
    panels: Mapping[str, pd.Series],
    *,
    sample_native_panels: Mapping[str, frozenset[str]] | None = None,
    protected_panels: frozenset[str] = frozenset(),
    protected_features: frozenset[str] = frozenset(),
    full_profiles: pd.DataFrame | None = None,
    min_hop_fraction: float = 2e-5,
    max_hop_fraction: float = 5e-3,
) -> pd.DataFrame:
    """Subtract stoichiometric patterned-flowcell lane bleed (``2e-5 .. 5e-3``), preserving ``NaN``s."""
    out = values.copy()
    if out.empty or not panels:
        return out
    fracs = estimate_foreign_panel_fractions(
        out,
        panels,
        sample_native_panels=sample_native_panels,
        protected_panels=protected_panels,
    )
    hop = fracs.where(fracs >= min_hop_fraction, 0.0).clip(upper=max_hop_fraction)
    if float(hop.to_numpy().max(initial=0.0)) <= 0.0:
        return out
    prof = _build_profile_matrix(out.index, panels, full_profiles=full_profiles)
    if protected_features:
        keep_zero = prof.index.isin(protected_features)
        if keep_zero.any():
            prof.loc[keep_zero, :] = 0.0
    bleed = prof.to_numpy(dtype=np.float64) @ hop.to_numpy(dtype=np.float64)
    arr = out.to_numpy(dtype=np.float64, copy=True)
    finite = np.isfinite(arr)
    arr[finite] = np.fmax(0.0, arr[finite] - bleed[finite])
    return pd.DataFrame(arr, index=out.index, columns=out.columns)


def deconvolve_panel_admixture(
    values: pd.DataFrame,
    panels: Mapping[str, pd.Series],
    *,
    sample_native_panels: Mapping[str, frozenset[str]] | None = None,
    protected_panels: frozenset[str] = frozenset(),
    protected_features: frozenset[str] = frozenset(),
    full_profiles: pd.DataFrame | None = None,
    min_shadow_weight: float = 0.005,
    max_subtract_fraction: float = 0.85,
) -> pd.DataFrame:
    """Fit NNLS foreign-archetype weights and subtract off-lineage background admixture."""
    out = values.copy()
    if out.empty or not panels:
        return out
    panel_names = sorted(panels)
    basis_genes = sorted({g for p in panel_names for g in panels[p].index if g in out.index})
    if not basis_genes:
        return out
    basis = _build_profile_matrix(pd.Index(basis_genes), panels, full_profiles=full_profiles)
    dedup_vals = out.loc[out.index.isin(basis_genes)].groupby(level=0).sum(min_count=1)
    weights = solve_nonnegative_archetype_mixture(basis, dedup_vals)
    w_arr = weights.to_numpy(dtype=np.float64, copy=True)
    for r_idx, pname in enumerate(weights.index):
        if pname in protected_panels:
            w_arr[r_idx, :] = 0.0
            continue
        if sample_native_panels:
            for c_idx, cname in enumerate(weights.columns):
                if pname in sample_native_panels.get(cname, frozenset()):
                    w_arr[r_idx, c_idx] = 0.0
    w_arr[w_arr < min_shadow_weight] = 0.0
    if float(w_arr.max(initial=0.0)) <= 0.0:
        return out
    prof = _build_profile_matrix(out.index, panels, full_profiles=full_profiles).reindex(
        columns=weights.index, fill_value=0.0
    )
    if protected_features:
        keep_zero = prof.index.isin(protected_features)
        if keep_zero.any():
            prof.loc[keep_zero, :] = 0.0
    shadow = prof.to_numpy(dtype=np.float64) @ w_arr
    arr = out.to_numpy(dtype=np.float64, copy=True)
    finite = np.isfinite(arr)
    sub = np.minimum(shadow[finite], max_subtract_fraction * np.fmax(arr[finite], 0.0))
    arr[finite] = np.fmax(0.0, arr[finite] - sub)
    return pd.DataFrame(arr, index=out.index, columns=out.columns)


def species_compartment_clean_scale(
    values: pd.DataFrame,
    *,
    symbols: pd.Series | None = None,
    target_total: float = 1_000_000.0,
) -> pd.DataFrame:
    """Species-native 16% ribosomal / 9% technical / 75% biological scale renormalization.

    Used for non-human bulk matrices (e.g. canine/murine TPM or TMM-CPM) so conserved
    ribosomal and mtDNA/rRNA abundance drift is removed without mapping to human ENSG IDs.
    Missing ``NaN`` entries remain ``NaN``; any absent compartment is recorded in
    ``result.attrs['unfilled_compartments']``.
    """
    if values.empty:
        return values.copy()
    sym_series = (
        pd.Series(values.index.astype(str), index=values.index)
        if symbols is None
        else symbols.reindex(values.index).fillna("").astype(str)
    )
    ribo_mask = sym_series.map(lambda s: bool(_RIBOSOMAL_SYMBOL_RE.match(s.strip()))).to_numpy()
    tech_mask = (
        ~ribo_mask
        & sym_series.map(lambda s: bool(_TECHNICAL_SYMBOL_RE.match(s.strip()))).to_numpy()
    )
    if not ribo_mask.any() and not tech_mask.any():
        raise ValueError(
            "species_compartment_clean_scale requires identifiable ribosomal or technical features "
            "(e.g. a GeneSymbol column); none matched"
        )
    bio_mask = ~(ribo_mask | tech_mask)

    arr = np.array(values.to_numpy(dtype=np.float64, na_value=np.nan), copy=True)
    out_arr = np.array(arr, copy=True)
    unfilled: dict[str, tuple[str, ...]] = {}
    cols = list(values.columns)

    for comp_name, mask, frac in (
        ("ribosomal_protein", ribo_mask, 0.16),
        ("technical", tech_mask, 0.09),
        ("biological", bio_mask, 0.75),
    ):
        if not mask.any():
            unfilled[comp_name] = tuple(str(c) for c in cols)
            continue
        sub = arr[mask]
        totals = np.nansum(sub, axis=0)
        pos = totals > 0.0
        if not pos.all():
            unfilled[comp_name] = tuple(str(cols[i]) for i in range(len(cols)) if not pos[i])
        scale = np.zeros_like(totals)
        np.divide(frac * target_total, totals, out=scale, where=pos)
        out_arr[mask] = sub * scale[None, :]

    res = pd.DataFrame(out_arr, index=values.index, columns=values.columns)
    if unfilled:
        res.attrs["unfilled_compartments"] = unfilled
    return res


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
    for m in range(1, 1001):
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
    "GTEX_ANCHOR_SITE_TO_PANEL",
    "NORMALIZE_ALIASES",
    "REPRODUCTIVE_ORGAN_PANELS",
    "SORTED_IMMUNE_CARRYOVER_PANELS",
    "SORTED_IMMUNE_CARRYOVER_PROFILES",
    "VERTEBRATE_ORGAN_ANCHOR_SYMBOLS",
    "BulkProcessingSpec",
    "deconvolve_panel_admixture",
    "estimate_foreign_panel_fractions",
    "expand_ops",
    "mask_panel_splash_cells",
    "purity_weighted_harrell_davis_quantile",
    "read_or_compute_cached_parquet",
    "regularized_betainc",
    "resolve_normalize_option",
    "solve_nonnegative_archetype_mixture",
    "species_compartment_clean_scale",
    "subtract_panel_index_hopping",
]

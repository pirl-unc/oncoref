# Licensed under the Apache License, Version 2.0.
"""Unified GTEx normal-tissue RNA interface across versions, providers, QC, cleanup, and normalization (#778).

Separates three orthogonal operations into one central, cacheable :class:`GtexSpec`:

1. **``qc`` — Reject invalid columns or contaminated ``(gene, sample)`` cells**
   (*never alters surviving numbers; drops bad columns or masks splashed cells as ``NaN``*):
   - ``"cell_lines"``: drop *in vitro* ``Cells_EBV-transformed_lymphocytes`` & ``Cells_Cultured_fibroblasts``.
   - ``"swaps"``: drop globally quarantined foreign-organ tube swaps & cross-sex chimeras.
   - ``"splash"``: mask (``NaN``) positive sample cells materially corrupted by ``>= 0.5%``
     non-boundary foreign splash or biologically impossible opposite-sex bleed (never masks
     ``0.0`` zeros or high endogenous expression where foreign bleed is ``< 20%`` of observed TPM).
   Pass ``qc=None`` (untouched), ``qc="all"`` (alias for :data:`ALL_QC`), or any list/set of operations.

2. **``cleanup`` — Improve continuous TPM & donor-quantile estimates on surviving observations:**
   - ``"index_hopping"``: subtract stoichiometric patterned-flowcell multiplexed pool bleed
     (``2e-5 .. 5e-3`` using conservative 25th-percentile anchor ratios).
   - ``"deconvolution"``: apply joint 52-shard NNLS shadow-source suppression, sample
     ``purity_weights``, and purity-weighted Harrell-Davis Beta donor quantiles.
   Pass ``cleanup=None`` (no model subtraction/reweighting), ``cleanup="all"`` (alias for
   :data:`ALL_CLEANUP`), or any list/set of operations.

3. **``normalize`` — Cross-cohort expression scale renormalization (applied after ``qc`` & ``cleanup``):**
   - ``normalize=None`` (native ``TPM`` / ``nTPM``) or ``normalize="clean_tpm"`` (16/9/75
     multi-compartment renormalization via :func:`oncoref.normalization.clean_tpm`).
"""

from __future__ import annotations

import contextlib
import functools
import gzip
import hashlib
import json
import os
import re
import tempfile
import warnings
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from filelock import FileLock

from . import gtex_deconvolution, normal_tissues, normalization, reference_data

SOURCE_ID = "gtex_tissue_rna"
SHARD_SOURCE_ID = "gtex_individual_tissue_tpm"
HPA_TISSUE_GTEX_SOURCE_ID = "rna_tissue_gtex"
HPA_BRAIN_GTEX_SOURCE_ID = "rna_brain_gtex"

_UNSET = object()

#: Atomic sample/value rejection rules (drop columns or mask contaminated cells as NaN).
ALL_QC: tuple[str, ...] = (
    "cell_lines",
    "swaps",
    "splash",
)

#: Named aliases that expand directly to sets of atomic :data:`ALL_QC` operations.
QC_ALIASES: dict[str, frozenset[str]] = {
    "all": frozenset(ALL_QC),
    "none": frozenset(),
    "raw": frozenset(),
}

#: Atomic quantitative TPM / quantile estimation steps on surviving observations.
ALL_CLEANUP: tuple[str, ...] = (
    "index_hopping",
    "deconvolution",
)

#: Named aliases that expand directly to sets of atomic :data:`ALL_CLEANUP` operations.
CLEANUP_ALIASES: dict[str, frozenset[str]] = {
    "all": frozenset(ALL_CLEANUP),
    "none": frozenset(),
    "raw": frozenset(),
}

CELL_LINE_COLUMNS: frozenset[str] = frozenset(
    {
        "Cells_Cultured_fibroblasts",
        "Cells_EBV-transformed_lymphocytes",
        "Cells - Cultured fibroblasts",
        "Cells - Transformed fibroblasts",
        "Cells - EBV-transformed lymphocytes",
    }
)

GTEX_V11_SUBCOMPARTMENT_COLUMNS: frozenset[str] = frozenset(
    {
        "Colon_Transverse_Mixed_Cell",
        "Colon_Transverse_Mucosa",
        "Colon_Transverse_Muscularis",
        "Liver_Hepatocyte",
        "Liver_Mixed_Cell",
        "Liver_Portal_Tract",
        "Pancreas_Acini",
        "Pancreas_Islets",
        "Pancreas_Mixed_Cell",
        "Small_Intestine_Terminal_Ileum_Lymphoid_Aggregate",
        "Small_Intestine_Terminal_Ileum_Mixed_Cell",
        "Stomach_Mixed_Cell",
        "Stomach_Mucosa",
        "Stomach_Muscularis",
    }
)

HPA_UNIQUE_GTEX_TISSUES: frozenset[str] = frozenset({"retina"})

_HPA_TISSUE_TO_GTEX_SITE: dict[str, str] = {
    "adipose tissue": "adipose_subcutaneous",
    "adrenal gland": "adrenal_gland",
    "amygdala": "brain_amygdala",
    "basal ganglia": "brain_putamen_basal_ganglia",
    "breast": "breast_mammary_tissue",
    "caudate": "brain_caudate_basal_ganglia",
    "cerebellum": "brain_cerebellum",
    "cerebral cortex": "brain_cortex",
    "cervix": "cervix_ectocervix",
    "colon": "colon_transverse",
    "endometrium": "uterus",
    "esophagus": "esophagus_mucosa",
    "fallopian tube": "fallopian_tube",
    "heart muscle": "heart_left_ventricle",
    "hippocampal formation": "brain_hippocampus",
    "hippocampus": "brain_hippocampus",
    "hypothalamus": "brain_hypothalamus",
    "kidney": "kidney_cortex",
    "liver": "liver",
    "lung": "lung",
    "midbrain": "brain_substantia_nigra",
    "nucleus accumbens": "brain_nucleus_accumbens_basal_ganglia",
    "ovary": "ovary",
    "pancreas": "pancreas",
    "pituitary gland": "pituitary",
    "prostate": "prostate",
    "putamen": "brain_putamen_basal_ganglia",
    "retina": "retina",
    "salivary gland": "minor_salivary_gland",
    "skeletal muscle": "muscle_skeletal",
    "skin": "skin_sun_exposed_lower_leg",
    "small intestine": "small_intestine_terminal_ileum",
    "spinal cord": "brain_spinal_cord_cervical_c-1",
    "spleen": "spleen",
    "stomach": "stomach",
    "substantia nigra": "brain_substantia_nigra",
    "testis": "testis",
    "thyroid gland": "thyroid",
    "urinary bladder": "bladder",
    "vagina": "vagina",
}


def _expand_ops(
    value: str | bool | Iterable[str] | None,
    *,
    allowed: tuple[str, ...],
    aliases: dict[str, frozenset[str]],
    param_name: str,
) -> frozenset[str]:
    """Expand ``None``, ``"all"`` (or other aliases), or any list/set of operations to a ``frozenset``."""
    if value is None or value is False:
        return frozenset()
    if value is True:
        return aliases["all"]
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
            f"Unknown GTEx {param_name} step(s) {sorted(set(unknown))}; "
            f"choose from {allowed}, aliases {sorted(aliases)}, or None"
        )
    return frozenset(expanded)


def expand_qc_ops(value: str | bool | Iterable[str] | None) -> frozenset[str]:
    """Expand a ``qc`` alias (e.g. ``'all'``) or list/set of operations into canonical :data:`ALL_QC` members."""
    return _expand_ops(value, allowed=ALL_QC, aliases=QC_ALIASES, param_name="qc")


def expand_cleanup_ops(value: str | bool | Iterable[str] | None) -> frozenset[str]:
    """Expand a ``cleanup`` alias (e.g. ``'all'``) or list/set of operations into canonical :data:`ALL_CLEANUP` members."""
    return _expand_ops(value, allowed=ALL_CLEANUP, aliases=CLEANUP_ALIASES, param_name="cleanup")


@dataclass(frozen=True)
class GtexSpec:
    """Composable GTEx sample/value ``qc``, quantitative ``cleanup``, and ``normalize`` spec (#778)."""

    qc: frozenset[str] | str | bool | Iterable[str] | None = frozenset()
    cleanup: frozenset[str] | str | bool | Iterable[str] | None = frozenset()
    normalize: str | None = None
    sex: str | None = None
    include_subcompartments: bool = True
    hpa_unique_tissues_only: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "qc", expand_qc_ops(self.qc))
        object.__setattr__(self, "cleanup", expand_cleanup_ops(self.cleanup))
        norm_val = (
            None
            if self.normalize is None or str(self.normalize).strip().lower() == "none"
            else str(self.normalize).strip().lower()
        )
        if norm_val is not None and norm_val != "clean_tpm":
            raise ValueError(f"normalize must be None or 'clean_tpm'; got {self.normalize!r}")
        object.__setattr__(self, "normalize", norm_val)
        if self.sex is not None and self.sex not in {"XX", "XY"}:
            raise ValueError(f"sex must be None, 'XX', or 'XY'; got {self.sex!r}")

    @property
    def drop_cell_lines(self) -> bool:
        return "cell_lines" in self.qc

    @property
    def quarantine_swaps(self) -> bool:
        return "swaps" in self.qc

    @property
    def reject_splash(self) -> bool:
        return "splash" in self.qc

    @property
    def subtract_index_hopping(self) -> bool:
        return "index_hopping" in self.cleanup

    @property
    def nnls_deconvolution(self) -> bool:
        return "deconvolution" in self.cleanup

    @property
    def clean_tpm(self) -> bool:
        return self.normalize == "clean_tpm"

    @property
    def is_raw_default(self) -> bool:
        return (
            not self.qc
            and not self.cleanup
            and self.normalize is None
            and self.sex is None
            and self.include_subcompartments
            and not self.hpa_unique_tissues_only
        )

    def to_dict(self) -> dict:
        return {
            "qc": [x for x in ALL_QC if x in self.qc],
            "cleanup": [x for x in ALL_CLEANUP if x in self.cleanup],
            "normalize": self.normalize,
            "sex": self.sex,
            "include_subcompartments": self.include_subcompartments,
            "hpa_unique_tissues_only": self.hpa_unique_tissues_only,
        }


def resolve_spec(
    *,
    qc: GtexSpec | str | bool | Iterable[str] | None = None,
    cleanup: object = _UNSET,
    normalize: object = _UNSET,
    sex: object = _UNSET,
    include_subcompartments: object = _UNSET,
    hpa_unique_tissues_only: object = _UNSET,
) -> GtexSpec:
    """Construct and validate a :class:`GtexSpec` from ``qc``, ``cleanup``, and ``normalize``."""
    if isinstance(qc, GtexSpec):
        return GtexSpec(
            qc=qc.qc,
            cleanup=qc.cleanup if cleanup is _UNSET else cleanup,  # type: ignore[arg-type]
            normalize=qc.normalize if normalize is _UNSET else normalize,  # type: ignore[arg-type]
            sex=qc.sex if sex is _UNSET else sex,  # type: ignore[arg-type]
            include_subcompartments=(
                qc.include_subcompartments
                if include_subcompartments is _UNSET
                else bool(include_subcompartments)
            ),
            hpa_unique_tissues_only=(
                qc.hpa_unique_tissues_only
                if hpa_unique_tissues_only is _UNSET
                else bool(hpa_unique_tissues_only)
            ),
        )
    return GtexSpec(
        qc=qc,
        cleanup=None if cleanup is _UNSET else cleanup,  # type: ignore[arg-type]
        normalize=None if normalize is _UNSET else normalize,  # type: ignore[arg-type]
        sex=None if sex is _UNSET else sex,  # type: ignore[arg-type]
        include_subcompartments=(
            True if include_subcompartments is _UNSET else bool(include_subcompartments)
        ),
        hpa_unique_tissues_only=(
            False if hpa_unique_tissues_only is _UNSET else bool(hpa_unique_tissues_only)
        ),
    )


def versions(provider: str = "direct", *, granularity: str = "medians") -> tuple[str, ...]:
    """Registered concrete GTEx versions for ``provider`` (``'direct'`` or ``'hpa'``)."""
    prov = provider.strip().lower()
    gran = granularity.strip().lower()
    if prov == "direct":
        if gran == "medians":
            return reference_data.versions(SOURCE_ID)
        if gran == "samples":
            pins = reference_data.REFERENCE_SOURCES[SHARD_SOURCE_ID]["pins"]
            rels = sorted({pk.split("-", 1)[0] for pk in pins}, reverse=True)
            return tuple(rels)
        raise ValueError(f"Unknown granularity {granularity!r}; expected 'medians' or 'samples'")
    if prov == "hpa":
        return reference_data.versions(HPA_TISSUE_GTEX_SOURCE_ID)
    raise ValueError(f"Unknown GTEx provider {provider!r}; expected 'direct' or 'hpa'")


def _primary_shard_sites() -> tuple[str, ...]:
    pins = sorted(reference_data.REFERENCE_SOURCES[SHARD_SOURCE_ID]["pins"])
    return tuple(pk.removeprefix("v11-") for pk in pins)


def normalize_column_to_site(column: str) -> str:
    """Normalize a GTEx v8/v10/v11 or HPA-GTEx column name to a canonical GTEx site key."""
    raw = str(column).strip()
    low = raw.lower()
    if low in _HPA_TISSUE_TO_GTEX_SITE:
        return _HPA_TISSUE_TO_GTEX_SITE[low]
    norm = (
        re.sub(r"[ ()-]+", "_", raw.replace("c-1", "cMINUS1"))
        .strip("_")
        .replace("cMINUS1", "c-1")
        .lower()
        .removeprefix("v11-")
    )
    primary = set(_primary_shard_sites())
    if norm in primary:
        return norm
    if norm == "small_intestine_terminal_ileum_lymphoid_aggregate":
        return "small_intestine_terminal_ileum"
    for suffix in (
        "_mixed_cell",
        "_mucosa",
        "_muscularis",
        "_hepatocyte",
        "_portal_tract",
        "_acini",
        "_islets",
    ):
        if norm.endswith(suffix) and norm[: -len(suffix)] in primary:
            return norm[: -len(suffix)]
    return norm


def _is_boundary_columns(src_col: str, tgt_col: str) -> bool:
    s_site = normalize_column_to_site(src_col)
    t_site = normalize_column_to_site(tgt_col)
    if gtex_deconvolution.is_anatomical_boundary(s_site, t_site):
        return True
    return bool(
        str(tgt_col).lower().endswith("lymphoid_aggregate")
        and gtex_deconvolution.organ_lineage_family(s_site) == "lymphoid_hematopoietic"
    )


def derived_cache_dir() -> Path:
    """Directory holding materialized GTEx Parquet tables (read-only; creates no dirs)."""
    return reference_data.cache_base_dir() / "derived" / "gtex"


def spec_cache_key(
    *,
    provider: str,
    version: str,
    granularity: str,
    stat: str,
    value_column: str,
    spec: GtexSpec,
    upstream_sha256: dict[str, str],
) -> tuple[str, dict]:
    """Compute a deterministic cache key and provenance manifest for a GTEx query spec."""
    payload = {
        "schema": "oncoref.gtex-derived-spec.v4",
        "policy_id": gtex_deconvolution.GTEX_DECONVOLUTION_POLICY_ID,
        "provider": provider,
        "version": version,
        "granularity": granularity,
        "stat": stat,
        "value_column": value_column,
        "spec": spec.to_dict(),
        "upstream_sha256": dict(sorted(upstream_sha256.items())),
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]
    slug = f"{provider}-{version}-{granularity}-{stat}-{value_column.lower()}"
    if spec.sex:
        slug += f"-{spec.sex.lower()}"
    slug = re.sub(r"[^a-zA-Z0-9_.-]+", "_", slug)
    key = f"{slug}-{digest}"
    payload["spec_key"] = key
    return key, payload


def cached_specs() -> pd.DataFrame:
    """Read-only inventory of materialized derived GTEx tables in ``derived/gtex``."""
    cdir = derived_cache_dir()
    cols = [
        "spec_key",
        "provider",
        "version",
        "granularity",
        "stat",
        "value_column",
        "qc",
        "cleanup",
        "normalize",
        "bytes",
        "path",
    ]
    if not cdir.is_dir():
        return pd.DataFrame(columns=cols)
    rows = []
    for meta_path in sorted(cdir.glob("*.json")):
        pq_path = meta_path.with_suffix(".parquet")
        if not pq_path.is_file():
            continue
        with contextlib.suppress(OSError, json.JSONDecodeError):
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            sp = meta.get("spec", {})
            rows.append(
                {
                    "spec_key": meta.get("spec_key", meta_path.stem),
                    "provider": meta.get("provider", ""),
                    "version": meta.get("version", ""),
                    "granularity": meta.get("granularity", ""),
                    "stat": meta.get("stat", ""),
                    "value_column": meta.get("value_column", ""),
                    "qc": ";".join(sp.get("qc", ())),
                    "cleanup": ";".join(sp.get("cleanup", ())),
                    "normalize": sp.get("normalize") or "none",
                    "bytes": pq_path.stat().st_size,
                    "path": str(pq_path),
                }
            )
    return pd.DataFrame(rows, columns=cols)


def clear_derived_cache(spec_key: str | None = None) -> int:
    """Remove one or all materialized GTEx derived Parquet caches; returns files removed."""
    cdir = derived_cache_dir()
    if not cdir.is_dir():
        return 0
    removed = 0
    targets = (
        [cdir / f"{spec_key}.parquet", cdir / f"{spec_key}.json"]
        if spec_key
        else list(cdir.iterdir())
    )
    for p in targets:
        if p.is_file() and not p.name.endswith(".lock"):
            with contextlib.suppress(OSError):
                p.unlink()
                removed += 1
    return removed


def _read_or_compute_cached_df(
    key: str,
    manifest_payload: dict,
    builder_fn,
    *,
    use_cache: bool,
) -> pd.DataFrame:
    if not use_cache:
        df = builder_fn()
        df.attrs["gtex_spec"] = manifest_payload
        return df

    cdir = derived_cache_dir()
    cdir.mkdir(parents=True, exist_ok=True)
    pq_path = cdir / f"{key}.parquet"
    meta_path = cdir / f"{key}.json"
    lock_path = cdir / f".{key}.lock"

    with FileLock(str(lock_path)):
        if pq_path.is_file() and meta_path.is_file():
            with contextlib.suppress(Exception):
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                if meta.get("spec_key") == key:
                    df = pd.read_parquet(pq_path)
                    df.attrs["gtex_spec"] = meta
                    df.attrs["cache_path"] = str(pq_path)
                    return df

        df = builder_fn()
        with contextlib.suppress(Exception):
            fd, tmp_pq = tempfile.mkstemp(dir=cdir, prefix=f".{key}-", suffix=".parquet.part")
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
        df.attrs["gtex_spec"] = manifest_payload
        if pq_path.is_file():
            df.attrs["cache_path"] = str(pq_path)
        return df


def apply_summary_spec(
    df: pd.DataFrame,
    *,
    spec: GtexSpec,
    provider: str = "direct",
    value_column: str = "TPM",
) -> pd.DataFrame:
    """Pass any GTEx summary matrix (direct ``v8``/``v10``/``v11`` or HPA ``v23``) through ``spec``."""
    if spec.sex is not None:
        raise ValueError(
            "Donor sex filtering (sex='XX' or 'XY') requires per-sample GTEx shards "
            "(provider='direct', version='v11')"
        )
    out = df.copy()
    id_cols = [c for c in ("Name", "Description") if c in out.columns]
    val_cols = [c for c in out.columns if c not in id_cols]

    if spec.drop_cell_lines:
        val_cols = [c for c in val_cols if c not in CELL_LINE_COLUMNS]
    if not spec.include_subcompartments:
        val_cols = [c for c in val_cols if c not in GTEX_V11_SUBCOMPARTMENT_COLUMNS]

    out = out[id_cols + val_cols].copy()
    if not val_cols:
        return out

    mat = out[val_cols].to_numpy(dtype=np.float64, copy=True)
    orig_nan = np.isnan(mat)

    # On pre-aggregated summary matrices, > 0.5% splash bleed is cleaned when either splash QC
    # or deconvolution cleanup is requested, while <= 0.5% lane hopping is subtracted when
    # index_hopping cleanup is requested.
    clean_summary_splash = spec.reject_splash or spec.nnls_deconvolution
    if spec.subtract_index_hopping or clean_summary_splash:
        anchors = gtex_deconvolution.gtex_deconvolution_anchors()
        unver = out["Name"].astype(str).str.split(".").str[0]
        first_pos: dict[str, int] = {}
        for pos, gid in enumerate(unver):
            if gid not in first_pos:
                first_pos[gid] = pos

        mat_clean = np.fmax(np.nan_to_num(mat, nan=0.0, posinf=0.0, neginf=0.0), 0.0)
        col_sites = [normalize_column_to_site(c) for c in val_cols]
        site_to_col_idxs: dict[str, list[int]] = {}
        for j, st in enumerate(col_sites):
            site_to_col_idxs.setdefault(st, []).append(j)

        nnls_w_mat: pd.DataFrame | None = None
        if spec.nnls_deconvolution:
            anchor_present = [g for g in anchors["gene_id"].unique() if g in first_pos]
            if len(anchor_present) >= 6:
                a_pos = [first_pos[g] for g in anchor_present]
                sub_df = pd.DataFrame(mat_clean[a_pos, :], index=anchor_present, columns=val_cols)
                basis_dict = {
                    st: sub_df[[val_cols[i] for i in idxs]].max(axis=1)
                    for st, idxs in site_to_col_idxs.items()
                }
                nnls_w_mat = gtex_deconvolution.solve_nonnegative_archetype_mixture(
                    pd.DataFrame(basis_dict, index=anchor_present),
                    sub_df,
                    l1_penalty=5e-4,
                )

        # Collect candidate bleed contributors per (target_col_j, source_lineage_family) so
        # sibling archetypes (e.g. 6 brain shards or 5 squamous shards) never double-subtract!
        best_by_col_fam: dict[tuple[int, str], tuple[float, float, np.ndarray]] = {}
        for src_site, grp in anchors.groupby("source_site", sort=True):
            if src_site not in site_to_col_idxs:
                continue
            gids = [g for g in grp["gene_id"] if g in first_pos]
            if len(gids) < 3:
                continue
            src_col_idxs = site_to_col_idxs[src_site]
            src_profile = mat_clean[:, src_col_idxs].max(axis=1)
            a_rows = np.array([first_pos[g] for g in gids], dtype=int)
            src_anchor_vals = src_profile[a_rows]
            valid_a = src_anchor_vals >= 20.0
            if int(valid_a.sum()) < 3:
                continue
            a_rows = a_rows[valid_a]
            src_anchor_vals = src_anchor_vals[valid_a]
            excl_sets = [
                {x for x in str(ex).split(";") if x}
                for ex, keep in zip(grp["excluded_target_sites"], valid_a)
                if keep
            ]
            src_fam = gtex_deconvolution.organ_lineage_family(str(src_site))

            for j, tgt_col in enumerate(val_cols):
                tgt_site = col_sites[j]
                if _is_boundary_columns(src_site, tgt_col):
                    continue
                keep_t = np.array([tgt_site not in ex for ex in excl_sets], dtype=bool)
                if int(keep_t.sum()) < 3:
                    continue
                use_rows = a_rows[keep_t]
                use_src = src_anchor_vals[keep_t]
                obs_a = mat[use_rows, j]
                det = np.isfinite(obs_a) & (obs_a >= 0.05)
                if int(det.sum()) < 3:
                    continue
                ratios = obs_a[det] / use_src[det]
                med_r = float(np.median(ratios))
                if med_r < 2e-5 or med_r > 0.25:
                    continue
                log_r = np.log10(np.maximum(ratios, 1e-12))
                q75, q25 = np.percentile(log_r, [75.0, 25.0])
                if float(q75 - q25) > 0.65:
                    continue
                if (
                    spec.nnls_deconvolution
                    and nnls_w_mat is not None
                    and med_r >= 5e-3
                    and src_site in nnls_w_mat.index
                    and float(nnls_w_mat.at[src_site, tgt_col]) < 0.15 * med_r
                ):
                    continue

                is_hop = med_r <= 5e-3
                if (is_hop and not spec.subtract_index_hopping) or (
                    not is_hop and not clean_summary_splash
                ):
                    continue

                coeff = float(np.percentile(ratios, 25.0)) if is_hop else med_r
                prev = best_by_col_fam.get((j, src_fam))
                if prev is None or med_r > prev[0]:
                    best_by_col_fam[(j, src_fam)] = (med_r, coeff, src_profile)

        expected_bleed = np.zeros_like(mat_clean)
        for (j, _sfam), (_med_r, coeff, src_profile) in best_by_col_fam.items():
            expected_bleed[:, j] += coeff * src_profile

        substantial = (
            ~orig_nan
            & (mat_clean > 0.0)
            & (expected_bleed > 0.0)
            & (expected_bleed >= 0.20 * mat_clean)
        )
        mat = np.where(substantial, np.fmax(0.0, mat_clean - expected_bleed), mat)
        mat[orig_nan] = np.nan
        out[val_cols] = mat

    if provider.strip().lower() == "hpa" and spec.hpa_unique_tissues_only:
        keep_unique = [c for c in val_cols if c.lower() in HPA_UNIQUE_GTEX_TISSUES]
        out = out[id_cols + keep_unique].copy()
        val_cols = keep_unique
        orig_nan = (
            np.isnan(out[val_cols].to_numpy(dtype=np.float64))
            if val_cols
            else np.zeros((0, 0), dtype=bool)
        )

    if spec.clean_tpm and val_cols:
        if value_column.upper() != "TPM":
            raise ValueError("normalize='clean_tpm' requires linear 'TPM', not 'nTPM'")
        gene_tbl = pd.DataFrame({"Ensembl_Gene_ID": out["Name"].astype(str).str.split(".").str[0]})
        filled = out[val_cols].fillna(0.0)
        norm_arr = normalization.clean_tpm(
            filled, gene_tbl, on_empty_compartment="record"
        ).to_numpy(dtype=np.float64, copy=True)
        norm_arr[orig_nan] = np.nan
        out[val_cols] = norm_arr

    return out


def _load_hpa_gtex_wide(version: str, *, value_column: str = "TPM") -> tuple[pd.DataFrame, dict]:
    if value_column not in {"TPM", "nTPM", "pTPM"}:
        raise ValueError(f"value_column must be 'TPM', 'nTPM', or 'pTPM'; got {value_column!r}")
    v_t = reference_data.resolve_version(HPA_TISSUE_GTEX_SOURCE_ID, version)
    v_b = reference_data.resolve_version(HPA_BRAIN_GTEX_SOURCE_ID, version)
    p_t = reference_data.ensure(HPA_TISSUE_GTEX_SOURCE_ID, v_t)
    p_b = reference_data.ensure(HPA_BRAIN_GTEX_SOURCE_ID, v_b)
    df_t = pd.read_csv(p_t, sep="\t")
    df_b = pd.read_csv(p_b, sep="\t").rename(columns={"Brain region": "Tissue"})

    t_existing = set(df_t["Tissue"].unique())
    df_b_extra = df_b[~df_b["Tissue"].isin(t_existing)]
    combined = pd.concat([df_t, df_b_extra], ignore_index=True)

    wide = (
        combined.pivot_table(
            index=["Gene", "Gene name"],
            columns="Tissue",
            values=value_column,
            aggfunc="first",
        )
        .reset_index()
        .rename(columns={"Gene": "Name", "Gene name": "Description"})
    )
    wide.columns.name = None
    prov_t = reference_data.provenance(HPA_TISSUE_GTEX_SOURCE_ID, v_t)
    prov_b = reference_data.provenance(HPA_BRAIN_GTEX_SOURCE_ID, v_b)
    up_sha = {
        f"{HPA_TISSUE_GTEX_SOURCE_ID}:{v_t}": str(prov_t.get("sha256") or ""),
        f"{HPA_BRAIN_GTEX_SOURCE_ID}:{v_b}": str(prov_b.get("sha256") or ""),
    }
    return wide, up_sha


def tissue_rna(
    version: str | None = None,
    *,
    provider: str = "direct",
    qc: GtexSpec | str | bool | Iterable[str] | None = None,
    cleanup: object = _UNSET,
    normalize: object = _UNSET,
    value_column: str = "TPM",
    include_subcompartments: object = _UNSET,
    hpa_unique_tissues_only: object = _UNSET,
    use_cache: bool = True,
) -> pd.DataFrame:
    """Return a GTEx gene x tissue summary matrix across any version, provider, ``qc``, and ``cleanup``."""
    prov = provider.strip().lower()
    spec = resolve_spec(
        qc=qc,
        cleanup=cleanup,
        normalize=normalize,
        include_subcompartments=include_subcompartments,
        hpa_unique_tissues_only=hpa_unique_tissues_only,
    )
    if spec.sex is not None:
        raise ValueError(
            "Donor sex filtering (sex='XX' or 'XY') requires per-sample GTEx shards "
            "(provider='direct', version='v11')"
        )

    if prov == "direct":
        concrete = reference_data.resolve_version(SOURCE_ID, version)
        if spec.is_raw_default:
            path = reference_data.ensure(SOURCE_ID, concrete)
            result = pd.read_csv(path, sep="\t", skiprows=2)
            result.attrs.update(reference_data.provenance(SOURCE_ID, concrete))
            result.attrs["unit"] = "TPM"
            result.attrs["statistic"] = "donor median"
            result.attrs["provider"] = "direct"
            result.attrs["spec"] = spec.to_dict()
            return result

        prov_rec = reference_data.provenance(SOURCE_ID, concrete)
        up_sha = {f"{SOURCE_ID}:{concrete}": str(prov_rec.get("sha256") or "")}
        key, manifest_payload = spec_cache_key(
            provider="direct",
            version=concrete,
            granularity="medians",
            stat="median",
            value_column="TPM",
            spec=spec,
            upstream_sha256=up_sha,
        )

        def _build_direct() -> pd.DataFrame:
            p = reference_data.ensure(SOURCE_ID, concrete)
            raw_df = pd.read_csv(p, sep="\t", skiprows=2)
            return apply_summary_spec(raw_df, spec=spec, provider="direct", value_column="TPM")

        result = _read_or_compute_cached_df(
            key, manifest_payload, _build_direct, use_cache=use_cache
        )
        result.attrs.update(prov_rec)
        result.attrs["unit"] = "clean_TPM" if spec.clean_tpm else "TPM"
        result.attrs["statistic"] = "donor median"
        result.attrs["provider"] = "direct"
        result.attrs["spec"] = spec.to_dict()
        return result

    if prov == "hpa":
        concrete = reference_data.resolve_version(HPA_TISSUE_GTEX_SOURCE_ID, version)
        prov_t = reference_data.provenance(HPA_TISSUE_GTEX_SOURCE_ID, concrete)
        prov_b = reference_data.provenance(HPA_BRAIN_GTEX_SOURCE_ID, concrete)
        up_sha = {
            f"{HPA_TISSUE_GTEX_SOURCE_ID}:{concrete}": str(prov_t.get("sha256") or ""),
            f"{HPA_BRAIN_GTEX_SOURCE_ID}:{concrete}": str(prov_b.get("sha256") or ""),
        }
        key, manifest_payload = spec_cache_key(
            provider="hpa",
            version=concrete,
            granularity="means",
            stat="mean",
            value_column=value_column,
            spec=spec,
            upstream_sha256=up_sha,
        )

        def _build_hpa() -> pd.DataFrame:
            wide, _ = _load_hpa_gtex_wide(concrete, value_column=value_column)
            if spec.is_raw_default:
                return wide
            return apply_summary_spec(wide, spec=spec, provider="hpa", value_column=value_column)

        result = _read_or_compute_cached_df(key, manifest_payload, _build_hpa, use_cache=use_cache)
        result.attrs["version"] = concrete
        result.attrs["unit"] = "clean_TPM" if spec.clean_tpm else value_column
        result.attrs["statistic"] = "tissue arithmetic mean (GTEx v8 via HPA)"
        result.attrs["provider"] = "hpa"
        result.attrs["spec"] = spec.to_dict()
        return result

    raise ValueError(f"Unknown GTEx provider {provider!r}; expected 'direct' or 'hpa'")


def donor_sex() -> pd.DataFrame:
    """Return splash-robust chromosomal sex (``XX`` vs ``XY``) for all 946 GTEx v11 donors."""
    return gtex_deconvolution.gtex_donor_sex()


def coverage_grid(*, category: str | None = None, missing_from_gtex: bool = False) -> pd.DataFrame:
    """Return the comprehensive human normal-tissue x source coverage grid (#778)."""
    return normal_tissues.normal_tissue_source_grid(
        category=category, missing_from_gtex=missing_from_gtex
    )


def tissues(
    version: str | None = None,
    *,
    provider: str = "direct",
    include_subcompartments: bool = True,
    include_cell_lines: bool = False,
) -> pd.DataFrame:
    """List available GTEx tissues/shards with sample/donor counts and lineage families."""
    prov = provider.strip().lower()
    if prov == "hpa":
        concrete = reference_data.resolve_version(HPA_TISSUE_GTEX_SOURCE_ID, version)
        rows = [
            {
                "tissue": hpa_t,
                "canonical_site": site,
                "lineage_family": gtex_deconvolution.organ_lineage_family(site),
                "kind": (
                    "hpa_unique_eyegex" if hpa_t in HPA_UNIQUE_GTEX_TISSUES else "hpa_gtex_v8_mean"
                ),
                "provider": "hpa",
                "version": concrete,
            }
            for hpa_t, site in sorted(_HPA_TISSUE_TO_GTEX_SITE.items())
        ]
        return pd.DataFrame(rows)

    if prov != "direct":
        raise ValueError(f"Unknown GTEx provider {provider!r}; expected 'direct' or 'hpa'")
    concrete = reference_data.resolve_version(SOURCE_ID, version)
    tsum = gtex_deconvolution.gtex_deconvolution_tissue_summary().set_index("site")
    rows = []
    for site in _primary_shard_sites():
        rec = tsum.loc[site] if site in tsum.index else None
        rows.append(
            {
                "tissue": site,
                "canonical_site": site,
                "lineage_family": gtex_deconvolution.organ_lineage_family(site),
                "kind": "primary_shard",
                "has_individual_samples": concrete == "v11",
                "total_samples": (
                    int(rec["total_samples"]) if rec is not None and concrete == "v11" else None
                ),
                "total_donors": (
                    int(rec["total_donors"]) if rec is not None and concrete == "v11" else None
                ),
                "provider": "direct",
                "version": concrete,
            }
        )
    if concrete == "v11" and include_subcompartments:
        for col in sorted(GTEX_V11_SUBCOMPARTMENT_COLUMNS):
            st = normalize_column_to_site(col)
            rows.append(
                {
                    "tissue": col,
                    "canonical_site": st,
                    "lineage_family": gtex_deconvolution.organ_lineage_family(st),
                    "kind": "v11_subcompartment_median",
                    "has_individual_samples": False,
                    "total_samples": None,
                    "total_donors": None,
                    "provider": "direct",
                    "version": concrete,
                }
            )
    if include_cell_lines:
        for col in ("Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes"):
            st = normalize_column_to_site(col)
            rows.append(
                {
                    "tissue": col,
                    "canonical_site": st,
                    "lineage_family": gtex_deconvolution.organ_lineage_family(st),
                    "kind": "transformed_or_cultured_cells",
                    "has_individual_samples": False,
                    "total_samples": None,
                    "total_donors": None,
                    "provider": "direct",
                    "version": concrete,
                }
            )
    return pd.DataFrame(rows)


def _resolve_shard_pin(tissue: str, version: str = "v11") -> tuple[str, str]:
    site = normalize_column_to_site(tissue)
    pin_key = f"{version}-{site}"
    pins = reference_data.REFERENCE_SOURCES[SHARD_SOURCE_ID]["pins"]
    if pin_key not in pins:
        avail = ", ".join(_primary_shard_sites())
        raise ValueError(f"Unknown GTEx {version} tissue shard {tissue!r}; available: {avail}")
    return site, pin_key


@functools.lru_cache(maxsize=2)
def _cached_v11_primary_medians(path_str: str, mtime_ns: int, size_bytes: int) -> pd.DataFrame:
    """Parse and cache unversioned primary-shard medians from ``gtex_tissue_rna:v11`` once."""
    del mtime_ns, size_bytes
    raw = pd.read_csv(path_str, sep="\t", skiprows=2)
    raw["_gid"] = raw["Name"].astype(str).str.split(".").str[0]
    raw = raw.drop_duplicates("_gid").set_index("_gid")
    site_cols: dict[str, pd.Series] = {}
    for c in raw.columns:
        if (
            c in ("Name", "Description")
            or c in GTEX_V11_SUBCOMPARTMENT_COLUMNS
            or c in CELL_LINE_COLUMNS
        ):
            continue
        site_cols[normalize_column_to_site(c)] = pd.to_numeric(raw[c], errors="coerce").fillna(0.0)
    return pd.DataFrame(site_cols, index=raw.index)


def _v11_primary_medians() -> pd.DataFrame:
    p = reference_data.ensure(SOURCE_ID, "v11")
    st = p.stat()
    return _cached_v11_primary_medians(str(p.resolve()), st.st_mtime_ns, st.st_size)


def _stream_gtex_shard(
    path: Path,
    *,
    genes: set[str] | None,
) -> pd.DataFrame:
    """Stream a GTEx ``.gct.gz`` shard (full or gene-subsetted)."""
    want_all = genes is None
    want_set = {str(g).strip() for g in genes} if genes is not None else set()
    want_bytes = {g.encode("ascii") for g in want_set}

    names: list[str] = []
    descs: list[str] = []
    vecs: list[np.ndarray] = []

    with gzip.open(path, "rb") as fh:
        fh.readline()
        fh.readline()
        header = fh.readline().decode("utf-8").rstrip("\n").split("\t")
        samps = header[2:]
        for line in fh:
            rn_b, _, rest = line.partition(b"\t")
            unver_b = rn_b.split(b".", 1)[0]
            desc_b, _, vals_b = rest.partition(b"\t")
            if not want_all and not (
                rn_b in want_bytes or unver_b in want_bytes or desc_b in want_bytes
            ):
                continue
            names.append(rn_b.decode("ascii"))
            descs.append(desc_b.decode("ascii"))
            vecs.append(np.asarray(vals_b.split(), dtype=np.float64))

    mat = np.vstack(vecs) if vecs else np.zeros((0, len(samps)), dtype=np.float64)
    df = pd.DataFrame(mat, columns=samps)
    df.insert(0, "Description", descs)
    df.insert(0, "Name", names)
    return df


def sample_tpm(
    tissue: str,
    *,
    genes: list[str] | set[str] | tuple[str, ...] | None = None,
    version: str = "v11",
    qc: GtexSpec | str | bool | Iterable[str] | None = None,
    cleanup: object = _UNSET,
    normalize: object = _UNSET,
    sex: object = _UNSET,
) -> pd.DataFrame:
    """Load per-sample TPM for a GTEx primary shard in strict ``qc -> cleanup -> normalize`` order.

    1. ``qc`` (sample/value rejection; never alters surviving observed numbers):
       - ``"swaps"``: drops globally quarantined foreign-organ swap / chimera columns.
       - ``"splash"``: masks (sets to ``NaN``, never fake ``0.0``!) positive cells in samples
         carrying ``>= 0.5%`` non-boundary foreign splash or biologically impossible opposite-sex
         bleed where expected bleed is ``>= 0.15 TPM`` and ``>= 20%`` of observed TPM.
    2. ``cleanup`` (quantitative TPM / purity estimation on surviving native-TPM observations):
       - ``"index_hopping"``: subtracts stoichiometric multiplexed pool hopping (``2e-5 .. 5e-3``)
         via conservative 25th-percentile anchor ratios (deduplicated per source lineage family).
       - ``"deconvolution"``: attaches sample ``purity_weights`` in ``df.attrs['purity_weights']``.
    3. ``normalize=None | "clean_tpm"`` (applied after ``qc`` and ``cleanup``):
       - Runs exact :func:`oncoref.normalization.clean_tpm` across all transcriptome compartments,
         preserving ``NaN`` masks before subsetting ``genes``.
    """
    spec = resolve_spec(qc=qc, cleanup=cleanup, normalize=normalize, sex=sex)
    site, pin_key = _resolve_shard_pin(tissue, version=version)
    shard_path = reference_data.ensure(SHARD_SOURCE_ID, pin_key)
    want_genes = {str(g).strip() for g in genes} if genes is not None else None
    need_live_anchors = spec.subtract_index_hopping or spec.reject_splash
    anchors_df = (
        gtex_deconvolution.gtex_deconvolution_anchors() if need_live_anchors else pd.DataFrame()
    )

    # Stream full shard when clean_tpm is requested so normalization.clean_tpm runs on the
    # complete cleaned transcriptome with zero inlined compartment math (#778 review fix).
    if spec.clean_tpm or want_genes is None:
        stream_genes = None
    elif need_live_anchors:
        stream_genes = want_genes | set(anchors_df["gene_id"])
    else:
        stream_genes = want_genes

    df_full = _stream_gtex_shard(shard_path, genes=stream_genes)
    samp_cols = [c for c in df_full.columns if c not in ("Name", "Description")]

    dsex_df = donor_sex().set_index("donor_id")["chromosomal_sex"].to_dict()
    swaps_set = (
        set(gtex_deconvolution.gtex_quarantined_swaps()["sample_id"])
        if spec.quarantine_swaps
        else set()
    )

    keep_samps: list[str] = []
    keep_samp_pos: list[int] = []
    for pos_j, s in enumerate(samp_cols):
        if s in swaps_set:
            continue
        if spec.sex is not None:
            did = gtex_deconvolution.parse_gtex_sample_barcode(s).donor_id
            if dsex_df.get(did) != spec.sex:
                continue
        keep_samps.append(s)
        keep_samp_pos.append(pos_j)

    anchor_obs_by_gid: dict[str, np.ndarray] = {}
    if need_live_anchors and len(df_full) > 0 and keep_samps:
        unver_all = df_full["Name"].astype(str).str.split(".").str[0].tolist()
        mat_keep_all = df_full[keep_samps].to_numpy(dtype=np.float64)
        a_set = set(anchors_df["gene_id"])
        for r_i, gid in enumerate(unver_all):
            if gid in a_set and gid not in anchor_obs_by_gid:
                anchor_obs_by_gid[gid] = mat_keep_all[r_i]

    # When clean_tpm is False and we streamed extra anchor genes, subset rows before vectorized work.
    if not spec.clean_tpm and want_genes is not None and stream_genes is not want_genes:
        unver_s = df_full["Name"].astype(str).str.split(".").str[0]
        keep_r = (
            df_full["Name"].isin(want_genes)
            | unver_s.isin(want_genes)
            | df_full["Description"].isin(want_genes)
        )
        df = df_full.loc[keep_r, ["Name", "Description", *keep_samps]].reset_index(drop=True)
    else:
        df = df_full[["Name", "Description", *keep_samps]].copy()

    purity_by_samp = dict.fromkeys(keep_samps, 1.0)
    if (
        (spec.reject_splash or spec.subtract_index_hopping or spec.nnls_deconvolution)
        and keep_samps
        and len(df) > 0
    ):
        ev_all = gtex_deconvolution.gtex_sample_splash_events(include_boundary=False)
        ev_site = ev_all[ev_all["target_site"].eq(site) & ev_all["sample_id"].isin(keep_samps)]
        med_primary = _v11_primary_medians()

        unver_rows = df["Name"].astype(str).str.split(".").str[0].tolist()
        vals = df[keep_samps].to_numpy(dtype=np.float64, copy=True)
        samp_idx = {s: j for j, s in enumerate(keep_samps)}
        src_med_vecs: dict[str, np.ndarray] = {}

        def _get_src_vec(ssrc: str) -> np.ndarray | None:
            if ssrc in src_med_vecs:
                return src_med_vecs[ssrc]
            if ssrc not in med_primary.columns:
                return None
            v = med_primary[ssrc].reindex(unver_rows).fillna(0.0).to_numpy(dtype=np.float64)
            src_med_vecs[ssrc] = v
            return v

        splashed_pairs: set[tuple[int, str]] = set()
        for sid, grp in ev_site.groupby("sample_id"):
            j = samp_idx[str(sid)]
            by_fam: dict[str, dict] = {}
            for rec in grp.to_dict("records"):
                if rec["regime"] == "tube_swap":
                    continue
                sfam = gtex_deconvolution.organ_lineage_family(str(rec["source_site"]))
                if sfam not in by_fam or float(rec["alpha_hat"]) > float(by_fam[sfam]["alpha_hat"]):
                    by_fam[sfam] = rec

            tot_alpha = sum(float(e["alpha_hat"]) for e in by_fam.values())
            purity_by_samp[str(sid)] = max(0.05, 1.0 - min(0.95, tot_alpha))

            for sfam, e in by_fam.items():
                ssrc = str(e["source_site"])
                a_hat = float(e["alpha_hat"])
                splashed_pairs.add((j, sfam))
                if (
                    spec.reject_splash
                    and a_hat >= 0.005
                    and e["regime"] in {"plate_pipetting_splash", "same_donor_carryover"}
                ):
                    src_vec = _get_src_vec(ssrc)
                    if src_vec is None:
                        continue
                    bleed_est = a_hat * src_vec
                    splash_mask = (
                        np.isfinite(vals[:, j])
                        & (vals[:, j] > 0.0)
                        & (bleed_est >= 0.15)
                        & (bleed_est >= 0.20 * vals[:, j])
                    )
                    vals[splash_mask, j] = np.nan

        if need_live_anchors and anchor_obs_by_gid:
            target_fam = gtex_deconvolution.organ_lineage_family(site)
            # Deduplicate live <0.5% hopping candidates per (sample_j, source_lineage_family)
            best_hop: dict[tuple[int, str], tuple[float, float, np.ndarray, bool]] = {}
            for src_site, grp in anchors_df.groupby("source_site", sort=True):
                src_fam = gtex_deconvolution.organ_lineage_family(str(src_site))
                if src_site == site or src_fam == target_fam:
                    continue
                if gtex_deconvolution.is_anatomical_boundary(str(src_site), site):
                    continue
                keep_a = [
                    g in anchor_obs_by_gid and site not in {x for x in str(ex).split(";") if x}
                    for g, ex in zip(grp["gene_id"], grp["excluded_target_sites"])
                ]
                sub_grp = grp.loc[keep_a]
                if len(sub_grp) < 3:
                    continue
                src_vec = _get_src_vec(str(src_site))
                if src_vec is None:
                    continue
                a_mat = np.vstack([anchor_obs_by_gid[g] for g in sub_grp["gene_id"]])
                a_meds = np.maximum(sub_grp["source_median_tpm"].to_numpy(dtype=np.float64), 1.0)
                det = np.isfinite(a_mat) & (a_mat >= 0.08)
                cand_js = np.where(det.sum(axis=0) >= 3)[0]
                for j_raw in cand_js:
                    j = int(j_raw)
                    if (j, src_fam) in splashed_pairs:
                        continue
                    dm = det[:, j]
                    r_vec = a_mat[dm, j] / a_meds[dm]
                    a_hat = float(np.median(r_vec))
                    if not (2e-5 <= a_hat < 5e-3):
                        continue
                    lr = np.log10(np.maximum(r_vec, 1e-12))
                    q75, q25 = np.percentile(lr, [75.0, 25.0])
                    if float(q75 - q25) > 0.65:
                        continue
                    a_lo = float(min(a_hat, np.percentile(r_vec, 25.0)))
                    samp_dsex = dsex_df.get(
                        gtex_deconvolution.parse_gtex_sample_barcode(keep_samps[j]).donor_id,
                        "unassessed",
                    )
                    cross_sex = (
                        samp_dsex == "XX"
                        and str(src_site) in gtex_deconvolution.MALE_SPECIFIC_SITES
                    ) or (
                        samp_dsex == "XY"
                        and str(src_site) in gtex_deconvolution.FEMALE_SPECIFIC_SITES
                    )
                    prev = best_hop.get((j, src_fam))
                    if prev is None or a_hat > prev[0]:
                        best_hop[(j, src_fam)] = (a_hat, a_lo, src_vec, cross_sex)

            for (j, _sfam), (a_hat, a_lo, src_vec, cross_sex) in best_hop.items():
                if spec.reject_splash and cross_sex:
                    bleed_est = a_hat * src_vec
                    splash_mask = (
                        np.isfinite(vals[:, j])
                        & (vals[:, j] > 0.0)
                        & (bleed_est >= 0.10)
                        & (bleed_est >= 0.20 * vals[:, j])
                    )
                    vals[splash_mask, j] = np.nan
                    continue
                if spec.subtract_index_hopping:
                    bleed_est = a_lo * src_vec
                    sub_mask = (
                        np.isfinite(vals[:, j])
                        & (vals[:, j] > 0.0)
                        & (bleed_est > 0.0)
                        & (bleed_est >= 0.25 * vals[:, j])
                    )
                    vals[sub_mask, j] = np.fmax(0.0, vals[sub_mask, j] - bleed_est[sub_mask])

        df[keep_samps] = vals

    if spec.clean_tpm and keep_samps and len(df) > 0:
        cur_vals = df[keep_samps].to_numpy(dtype=np.float64)
        nan_mask = np.isnan(cur_vals)
        gene_tbl = pd.DataFrame({"Ensembl_Gene_ID": df["Name"].astype(str).str.split(".").str[0]})
        clean_arr = normalization.clean_tpm(
            df[keep_samps].fillna(0.0), gene_tbl, on_empty_compartment="record"
        ).to_numpy(dtype=np.float64, copy=True)
        clean_arr[nan_mask] = np.nan
        df[keep_samps] = clean_arr

    if spec.clean_tpm and want_genes is not None:
        unver_s = df["Name"].astype(str).str.split(".").str[0]
        keep_r = (
            df["Name"].isin(want_genes)
            | unver_s.isin(want_genes)
            | df["Description"].isin(want_genes)
        )
        df = df.loc[keep_r].reset_index(drop=True)

    df.attrs["tissue"] = site
    df.attrs["version"] = version
    df.attrs["unit"] = "clean_TPM" if spec.clean_tpm else "TPM"
    df.attrs["spec"] = spec.to_dict()
    df.attrs["purity_weights"] = purity_by_samp
    df.attrs["donor_sex"] = {
        s: dsex_df.get(gtex_deconvolution.parse_gtex_sample_barcode(s).donor_id, "unassessed")
        for s in keep_samps
    }
    return df


def _match_summary_columns(mat_columns: Iterable[str], requested: list[str]) -> list[str]:
    val_cols = [c for c in mat_columns if c not in ("Name", "Description")]
    by_exact = {c: c for c in val_cols}
    by_lower = {c.lower(): c for c in val_cols}
    by_site: dict[str, list[str]] = {}
    for c in val_cols:
        by_site.setdefault(normalize_column_to_site(c), []).append(c)

    out: list[str] = []
    seen: set[str] = set()
    for req in requested:
        r = str(req).strip()
        if r in by_exact:
            matched = [by_exact[r]]
        elif r.lower() in by_lower:
            matched = [by_lower[r.lower()]]
        else:
            st = normalize_column_to_site(r)
            if st in by_site:
                matched = by_site[st]
            else:
                raise ValueError(f"Unknown GTEx tissue {req!r}; available: {val_cols}")
        for m in matched:
            if m not in seen:
                seen.add(m)
                out.append(m)
    return out


def tissue_summary(
    tissues_select: str | list[str] | tuple[str, ...] | None = None,
    *,
    genes: list[str] | set[str] | tuple[str, ...] | None = None,
    version: str | None = None,
    provider: str = "direct",
    stat: str = "p95",
    qc: object = _UNSET,
    cleanup: object = _UNSET,
    normalize: object = _UNSET,
    sex: object = _UNSET,
    value_column: str = "TPM",
    include_subcompartments: object = _UNSET,
    hpa_unique_tissues_only: object = _UNSET,
    use_cache: bool = True,
) -> pd.DataFrame:
    """Summarize GTEx donor expression across one or all tissues (#778).

    - ``stat``: ``"median"`` | ``"p95"`` | ``"p99"`` | ``"max"`` | ``"host_exposure"``
    - ``qc="all"`` (default when ``qc`` is not a :class:`GtexSpec`): rejects cell lines,
      quarantined swaps, and splashed cells; pass ``qc=None`` for no sample/value rejection.
    - ``cleanup="all"`` (default when ``qc`` is not a :class:`GtexSpec`): applies index-hopping
      subtraction and NNLS purity-weighted Harrell-Davis quantiles; pass ``cleanup=None`` for
      unadjusted TPMs on QC-passing observations.
    - ``normalize=None | "clean_tpm"``: optional 16/9/75 renormalization.
    - ``sex=None | "XX" | "XY"``: optional donor chromosomal sex filter on ``v11`` shards.
    """
    prov = provider.strip().lower()
    if version is None or version == "latest":
        version = "v23" if prov == "hpa" else "v11"

    if isinstance(qc, GtexSpec):
        spec = resolve_spec(
            qc=qc,
            cleanup=cleanup,
            normalize=normalize,
            sex=sex,
            include_subcompartments=include_subcompartments,
            hpa_unique_tissues_only=hpa_unique_tissues_only,
        )
    else:
        spec = resolve_spec(
            qc="all" if qc is _UNSET else qc,  # type: ignore[arg-type]
            cleanup="all" if cleanup is _UNSET else cleanup,
            normalize=normalize,
            sex=sex,
            include_subcompartments=include_subcompartments,
            hpa_unique_tissues_only=hpa_unique_tissues_only,
        )

    stat_norm = stat.strip().lower()
    valid_stats = {"median", "p95", "p99", "max", "host_exposure"}
    if stat_norm not in valid_stats:
        raise ValueError(f"Unknown stat {stat!r}; choose from {sorted(valid_stats)}")

    if prov == "hpa" or (prov == "direct" and version != "v11"):
        if prov == "direct" and stat_norm != "median":
            raise ValueError(
                f"GTEx {version} pins donor medians only; per-sample {stat!r} requires version='v11'"
            )
        mat = tissue_rna(
            version=version,
            provider=prov,
            qc=spec,
            value_column=value_column,
            use_cache=use_cache,
        )
        if genes is not None:
            wset = {str(g).strip() for g in genes}
            unver = mat["Name"].astype(str).str.split(".").str[0]
            mat = mat[
                mat["Name"].isin(wset) | unver.isin(wset) | mat["Description"].isin(wset)
            ].reset_index(drop=True)
        if tissues_select is not None:
            req = [tissues_select] if isinstance(tissues_select, str) else list(tissues_select)
            keep_c = _match_summary_columns(mat.columns, req) if req else []
            mat = mat[["Name", "Description", *keep_c]].copy()
        return mat

    if tissues_select is None:
        target_sites = list(_primary_shard_sites())
    elif isinstance(tissues_select, str):
        target_sites = [_resolve_shard_pin(tissues_select, version=version)[0]]
    else:
        target_sites = [_resolve_shard_pin(t, version=version)[0] for t in tissues_select]

    if not target_sites:
        empty = pd.DataFrame(columns=["Name", "Description"])
        empty.attrs["provider"] = "direct"
        empty.attrs["version"] = version
        empty.attrs["statistic"] = stat_norm
        empty.attrs["unit"] = "clean_TPM" if spec.clean_tpm else "TPM"
        empty.attrs["spec"] = spec.to_dict()
        return empty

    can_cache_full = use_cache and genes is None and tissues_select is None

    def _compute_summary() -> pd.DataFrame:
        base_meta: pd.DataFrame | None = None
        site_cols: dict[str, np.ndarray] = {}
        q_map = {"median": 0.50, "p95": 0.95, "p99": 0.99}

        for site in target_sites:
            sdf = sample_tpm(site, genes=genes, version=version, qc=spec)
            if base_meta is None:
                base_meta = sdf[["Name", "Description"]].copy()
            s_cols = [c for c in sdf.columns if c not in ("Name", "Description")]
            if not s_cols or len(sdf) == 0:
                site_cols[site] = np.full(len(sdf), np.nan, dtype=np.float64)
                continue

            vals = sdf[s_cols].to_numpy(dtype=np.float64)
            pw_map = sdf.attrs.get("purity_weights", {})
            p_wts = np.array([float(pw_map.get(s, 1.0)) for s in s_cols], dtype=np.float64)

            donors = [gtex_deconvolution.parse_gtex_sample_barcode(s).donor_id for s in s_cols]
            unique_donors = list(dict.fromkeys(donors))
            if len(unique_donors) != len(s_cols):
                d_idx_map: dict[str, list[int]] = {}
                for j_i, d_id in enumerate(donors):
                    d_idx_map.setdefault(d_id, []).append(j_i)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    d_vals = np.vstack(
                        [np.nanmedian(vals[:, idxs], axis=1) for idxs in d_idx_map.values()]
                    ).T
                d_wts = np.array([p_wts[idxs].mean() for idxs in d_idx_map.values()])
            else:
                d_vals = vals
                d_wts = p_wts

            out_vec = np.empty(len(sdf), dtype=np.float64)
            for r_i in range(len(sdf)):
                row_v = d_vals[r_i]
                fin = np.isfinite(row_v)
                if not fin.any():
                    out_vec[r_i] = np.nan
                    continue
                rv = row_v[fin]
                rw = d_wts[fin]
                if stat_norm == "max":
                    out_vec[r_i] = float(rv.max())
                elif stat_norm == "host_exposure":
                    sv = np.sort(rv)
                    nd = sv.size
                    if nd >= 50:
                        p99 = float(np.percentile(sv, 99.0))
                        out_vec[r_i] = min(
                            float(sv[-1]), max(p99, float(sv[-2])) if nd < 100 else p99
                        )
                    elif nd >= 30:
                        out_vec[r_i] = float(sv[-2])
                    else:
                        out_vec[r_i] = float(sv[-1])
                elif spec.nnls_deconvolution:
                    out_vec[r_i] = gtex_deconvolution.purity_weighted_harrell_davis_quantile(
                        rv, rw, quantile=q_map[stat_norm]
                    )
                else:
                    out_vec[r_i] = float(np.percentile(rv, 100.0 * q_map[stat_norm]))
            site_cols[site] = out_vec

        assert base_meta is not None
        return pd.concat([base_meta, pd.DataFrame(site_cols, index=base_meta.index)], axis=1)

    pins = reference_data.REFERENCE_SOURCES[SHARD_SOURCE_ID]["pins"]
    up_sha = {f"{SHARD_SOURCE_ID}:v11-{s}": pins[f"v11-{s}"]["sha256"] for s in target_sites}
    key, manifest_payload = spec_cache_key(
        provider="direct",
        version=version,
        granularity="samples",
        stat=stat_norm,
        value_column="TPM",
        spec=spec,
        upstream_sha256=up_sha,
    )
    out = _read_or_compute_cached_df(
        key, manifest_payload, _compute_summary, use_cache=can_cache_full
    )
    out.attrs["provider"] = "direct"
    out.attrs["version"] = version
    out.attrs["statistic"] = stat_norm
    out.attrs["unit"] = "clean_TPM" if spec.clean_tpm else "TPM"
    out.attrs["spec"] = spec.to_dict()
    return out


__all__ = [
    "ALL_CLEANUP",
    "ALL_QC",
    "CELL_LINE_COLUMNS",
    "CLEANUP_ALIASES",
    "GTEX_V11_SUBCOMPARTMENT_COLUMNS",
    "HPA_BRAIN_GTEX_SOURCE_ID",
    "HPA_TISSUE_GTEX_SOURCE_ID",
    "HPA_UNIQUE_GTEX_TISSUES",
    "QC_ALIASES",
    "SHARD_SOURCE_ID",
    "SOURCE_ID",
    "GtexSpec",
    "apply_summary_spec",
    "cached_specs",
    "clear_derived_cache",
    "coverage_grid",
    "derived_cache_dir",
    "donor_sex",
    "expand_cleanup_ops",
    "expand_qc_ops",
    "normalize_column_to_site",
    "resolve_spec",
    "sample_tpm",
    "spec_cache_key",
    "tissue_rna",
    "tissue_summary",
    "tissues",
    "versions",
]

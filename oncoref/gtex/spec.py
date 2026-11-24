# Licensed under the Apache License, Version 2.0.
"""Composable ``(qc, cleanup, normalize)`` specification and derived cache keys for GTEx (#778)."""

from __future__ import annotations

import contextlib
import hashlib
import json
import re
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .. import bulk_qc, reference_data
from . import deconvolution as gtex_deconvolution

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


def expand_qc_ops(value: str | bool | Iterable[str] | None) -> frozenset[str]:
    """Expand a ``qc`` alias (e.g. ``'all'``) or list/set of operations into canonical :data:`ALL_QC` members."""
    return bulk_qc.expand_ops(
        value, allowed=ALL_QC, aliases=QC_ALIASES, param_name="qc", domain_label="GTEx"
    )


def expand_cleanup_ops(value: str | bool | Iterable[str] | None) -> frozenset[str]:
    """Expand a ``cleanup`` alias (e.g. ``'all'``) or list/set of operations into canonical :data:`ALL_CLEANUP` members."""
    return bulk_qc.expand_ops(
        value,
        allowed=ALL_CLEANUP,
        aliases=CLEANUP_ALIASES,
        param_name="cleanup",
        domain_label="GTEx",
    )


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
        norm_val = bulk_qc.resolve_normalize_option(self.normalize, domain_label="GTEx")
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


def primary_shard_sites() -> tuple[str, ...]:
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
    primary = set(primary_shard_sites())
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


def derived_cache_dir() -> Path:
    """Directory holding materialized GTEx Parquet tables (read-only; creates no dirs)."""
    return reference_data.cache_base_dir() / "derived" / "gtex"


def _effective_derived_cache_dir() -> Path:
    pkg = sys.modules.get("oncoref.gtex")
    if pkg is not None and hasattr(pkg, "derived_cache_dir"):
        return pkg.derived_cache_dir()
    return derived_cache_dir()


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
    cdir = _effective_derived_cache_dir()
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
    cdir = _effective_derived_cache_dir()
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
    "cached_specs",
    "clear_derived_cache",
    "derived_cache_dir",
    "expand_cleanup_ops",
    "expand_qc_ops",
    "normalize_column_to_site",
    "primary_shard_sites",
    "resolve_spec",
    "spec_cache_key",
]

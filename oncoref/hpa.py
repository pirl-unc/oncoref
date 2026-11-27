# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Loaders over the Human Protein Atlas normal-tissue reference data.

These fetch (once, cached) and parse the HPA tables registered in
:mod:`oncoref.reference_data`:

- ``hpa_rna_consensus`` — per-tissue RNA nTPM (``Gene``, ``Gene name``,
  ``Tissue``, ``nTPM``).
- ``hpa_normal_tissue`` — IHC protein detection (``Gene``, ``Gene name``,
  ``Tissue``, ``Cell type``, ``Level``, ``Reliability``).
- ``hpa_single_cell`` — single-cell-type RNA expression (``Gene``, ``Gene name``,
  ``Cell type``, ``nCPM`` in v25.1+ / ``nTPM`` in v23).

This is the normal-tissue expression evidence behind the cancer-testis-antigen
tissue-restriction definition and protein-level / single-cell comparisons.
"""

from __future__ import annotations

import contextlib
import os
import tempfile
from dataclasses import dataclass
from functools import cache, lru_cache
from pathlib import Path

import pandas as pd

from . import bulk_qc, reference_data
from .cta_tissues import SAFETY_TISSUE_GROUPS
from .load_dataset import _register_derived_cache, get_data

ALL_HPA_BULK_QC: tuple[str, ...] = ("duplicate_gtex_tissues", "splash")
HPA_BULK_QC_ALIASES: dict[str, frozenset[str]] = {
    "all": frozenset(ALL_HPA_BULK_QC),
    "none": frozenset(),
    "raw": frozenset(),
}

ALL_HPA_BULK_CLEANUP: tuple[str, ...] = ("index_hopping", "deconvolution")
HPA_BULK_CLEANUP_ALIASES: dict[str, frozenset[str]] = {
    "all": frozenset(ALL_HPA_BULK_CLEANUP),
    "none": frozenset(),
    "raw": frozenset(),
}

#: HPA tissue labels that duplicate GTEx v11's 52 direct donor shards (excluding EyeGEx retina).
HPA_DUPLICATE_GTEX_TISSUES: frozenset[str] = frozenset(
    {
        "adipose tissue",
        "adrenal gland",
        "amygdala",
        "basal ganglia",
        "blood vessel",
        "breast",
        "caudate",
        "cerebellum",
        "cerebral cortex",
        "cervix",
        "colon",
        "endometrium",
        "esophagus",
        "fallopian tube",
        "heart muscle",
        "hippocampal formation",
        "hippocampus",
        "hypothalamus",
        "kidney",
        "liver",
        "lung",
        "midbrain",
        "nucleus accumbens",
        "ovary",
        "pancreas",
        "pituitary gland",
        "prostate",
        "putamen",
        "salivary gland",
        "skeletal muscle",
        "skin",
        "small intestine",
        "spinal cord",
        "spleen",
        "stomach",
        "substantia nigra",
        "testis",
        "thyroid gland",
        "urinary bladder",
        "vagina",
    }
)

_HPA_TISSUE_NATIVE_PANELS: dict[str, frozenset[str]] = {
    "adrenal gland": frozenset({"adrenal_gland"}),
    "amygdala": frozenset({"brain"}),
    "basal ganglia": frozenset({"brain"}),
    "caudate": frozenset({"brain"}),
    "cerebellum": frozenset({"brain"}),
    "cerebral cortex": frozenset({"brain"}),
    "choroid plexus": frozenset({"brain"}),
    "hippocampal formation": frozenset({"brain"}),
    "hippocampus": frozenset({"brain"}),
    "hypothalamus": frozenset({"brain"}),
    "midbrain": frozenset({"brain"}),
    "nucleus accumbens": frozenset({"brain"}),
    "putamen": frozenset({"brain"}),
    "substantia nigra": frozenset({"brain"}),
    "spinal cord": frozenset({"brain"}),
    "retina": frozenset({"brain"}),
    "esophagus": frozenset({"squamous_epithelium"}),
    "skin": frozenset({"squamous_epithelium"}),
    "tonsil": frozenset({"squamous_epithelium"}),
    "tongue": frozenset({"squamous_epithelium", "striated_muscle"}),
    "cervix": frozenset({"squamous_epithelium"}),
    "vagina": frozenset({"squamous_epithelium"}),
    "heart muscle": frozenset({"striated_muscle"}),
    "skeletal muscle": frozenset({"striated_muscle"}),
    "kidney": frozenset({"kidney"}),
    "liver": frozenset({"liver"}),
    "gallbladder": frozenset({"liver"}),
    "lung": frozenset({"lung"}),
    "pancreas": frozenset({"pancreas"}),
    "pituitary gland": frozenset({"pituitary"}),
    "salivary gland": frozenset({"salivary_gland"}),
    "stomach": frozenset({"stomach"}),
    "duodenum": frozenset({"stomach"}),
    "small intestine": frozenset({"stomach"}),
    "colon": frozenset({"stomach"}),
    "rectum": frozenset({"stomach"}),
    "appendix": frozenset({"stomach"}),
    "testis": frozenset({"testis"}),
    "epididymis": frozenset({"testis"}),
    "thyroid gland": frozenset({"thyroid"}),
    "parathyroid gland": frozenset({"thyroid"}),
}


def _resolve_hpa_native_panels(
    tissue_label: str, *, brain_subregions: frozenset[str] = frozenset()
) -> frozenset[str]:
    low = tissue_label.strip().lower()
    if low in _HPA_TISSUE_NATIVE_PANELS:
        return _HPA_TISSUE_NATIVE_PANELS[low]
    if low in brain_subregions:
        if "pituitary" in low or "hypophys" in low:
            return frozenset({"brain", "pituitary"})
        return frozenset({"brain"})
    inferred: set[str] = set()
    if (
        any(
            tok in low
            for tok in (
                "cortex",
                "gyrus",
                "nucleus",
                "thalam",
                "callosum",
                "cerebell",
                "medulla",
                "pons",
                "tegment",
                "white matter",
                "olfactory",
                "optic",
                "ventricle",
                "brain",
                "spinal",
            )
        )
        and "kidney" not in low
        and "renal" not in low
        and "adrenal" not in low
        and "heart" not in low
    ):
        inferred.add("brain")
    if "pituitary" in low or "hypophys" in low:
        inferred.add("pituitary")
    if ("muscle" in low and "smooth" not in low) or "heart" in low:
        inferred.add("striated_muscle")
    return frozenset(inferred)


SAFETY_TISSUE_MAPPING_KINDS: frozenset[str] = frozenset(
    {"exact", "reviewed_equivalent", "reviewed_substructure", "unavailable"}
)
SAFETY_TISSUE_COVERAGE_LEVELS: frozenset[str] = frozenset({"complete", "partial", "unavailable"})
SAFETY_TISSUE_COVERAGE_STATES: frozenset[str] = frozenset({"complete", "partial", "unavailable"})

_SAFETY_TISSUE_MAPPING_COLUMNS = (
    "source_name",
    "source_version",
    "safety_group",
    "requested_tissue",
    "source_tissue",
    "mapping_kind",
    "coverage_level",
    "mapping_reference",
    "notes",
)


class SafetyTissueResolutionError(ValueError):
    """A safety-tissue group cannot be resolved under the requested contract."""


@dataclass(frozen=True)
class SafetyTissueMapping:
    """One conceptual safety tissue mapped to source-native tissue labels."""

    requested_tissue: str
    source_tissues: tuple[str, ...]
    mapping_kind: str
    coverage_level: str
    mapping_references: tuple[str, ...]
    notes: tuple[str, ...]

    @property
    def available(self) -> bool:
        """Whether at least one source-native tissue represents this concept."""
        return bool(self.source_tissues)

    @property
    def complete(self) -> bool:
        """Whether the source fully represents this conceptual tissue."""
        return self.coverage_level == "complete"


@dataclass(frozen=True)
class SafetyTissueResolution:
    """Deterministic source-aware resolution of a safety-tissue group."""

    safety_group: str
    source_name: str
    source_version: str
    source_url: str
    mappings: tuple[SafetyTissueMapping, ...]
    coverage_state: str

    @property
    def requested_tissues(self) -> tuple[str, ...]:
        """Conceptual tissues configured for the safety group."""
        return tuple(mapping.requested_tissue for mapping in self.mappings)

    @property
    def source_tissues(self) -> tuple[str, ...]:
        """Deduplicated source-native labels safe to use for exact matching."""
        return tuple(
            sorted({tissue for mapping in self.mappings for tissue in mapping.source_tissues})
        )

    @property
    def fully_covered_tissues(self) -> tuple[str, ...]:
        """Conceptual tissues with complete source coverage."""
        return tuple(
            mapping.requested_tissue
            for mapping in self.mappings
            if mapping.coverage_level == "complete"
        )

    @property
    def partially_covered_tissues(self) -> tuple[str, ...]:
        """Conceptual tissues represented only by reviewed substructures."""
        return tuple(
            mapping.requested_tissue
            for mapping in self.mappings
            if mapping.coverage_level == "partial"
        )

    @property
    def unavailable_tissues(self) -> tuple[str, ...]:
        """Conceptual tissues with no source-native label in this source version."""
        return tuple(
            mapping.requested_tissue
            for mapping in self.mappings
            if mapping.coverage_level == "unavailable"
        )

    @property
    def complete(self) -> bool:
        """Whether every requested conceptual tissue is completely covered."""
        return self.coverage_state == "complete"


def _read_hpa(name: str, version: str | None = None) -> pd.DataFrame:
    """Load an HPA table, caching a columnar **parquet** copy next to the raw TSV
    so repeated/cold reads are fast and compact (HPA TSVs are large — single-cell
    is tens of MB; parquet is a few-fold smaller and skips re-parsing). The parquet
    cache is best-effort and regenerated if the TSV is re-downloaded."""
    concrete = reference_data.resolve_version(name, version)
    tsv = reference_data.ensure(name, concrete)
    parquet = tsv.with_suffix(".parquet")
    # Older writers could expose an incomplete cache to another reader. A
    # damaged optional cache must never prevent reading the authoritative TSV.
    with contextlib.suppress(Exception):
        if parquet.exists() and parquet.stat().st_mtime >= tsv.stat().st_mtime:
            result = pd.read_parquet(parquet)
            result.attrs.update(reference_data.provenance(name, concrete))
            return result
    df = pd.read_csv(tsv, sep="\t")
    # parquet cache is an optimization, not required — never fail a read over it.
    with contextlib.suppress(Exception):
        fd, temporary = tempfile.mkstemp(dir=parquet.parent, suffix=".parquet.part")
        os.close(fd)
        try:
            df.to_parquet(temporary, index=False)
            os.replace(temporary, parquet)
        finally:
            Path(temporary).unlink(missing_ok=True)
    df.attrs.update(reference_data.provenance(name, concrete))
    return df


@cache
def _hpa_rna_consensus_for_version(version: str) -> pd.DataFrame:
    return _read_hpa("hpa_rna_consensus", version)


def _human_marker_panels_and_profiles() -> tuple[dict[str, pd.Series], pd.DataFrame]:
    panels_df = get_data("tissue-marker-panels.csv")
    panels: dict[str, pd.Series] = {}
    for pname, grp in panels_df.groupby("panel", sort=True):
        s = grp.set_index("gene_id")["reference_tpm"].astype(float)
        panels[str(pname)] = s[~s.index.duplicated(keep="first")]
    anchors_df = get_data("gtex-deconvolution-anchors")
    # Ensure adrenal_gland (and ovary/prostate) anchor panels are populated alongside tissue-marker-panels.csv
    for extra_site in ("adrenal_gland", "ovary", "prostate"):
        if extra_site not in panels:
            sub_a = anchors_df[anchors_df["source_site"].eq(extra_site)]
            if not sub_a.empty:
                s = sub_a.set_index("gene_id")["source_median_tpm"].astype(float)
                panels[extra_site] = s[~s.index.duplicated(keep="first")]

    prof_rows: dict[str, dict[str, float]] = {}
    for pname, ser in panels.items():
        for gid, val in ser.items():
            prof_rows.setdefault(str(gid), {})[pname] = float(val)
    for row in anchors_df.itertuples(index=False):
        mapped = bulk_qc.GTEX_ANCHOR_SITE_TO_PANEL.get(str(row.source_site))
        if mapped and mapped in panels:
            prof_rows.setdefault(str(row.gene_id), {})[mapped] = max(
                prof_rows.get(str(row.gene_id), {}).get(mapped, 0.0),
                float(row.source_median_tpm),
            )
    prof_df = pd.DataFrame.from_dict(prof_rows, orient="index").fillna(0.0)
    return panels, prof_df


def _process_hpa_long_bulk_table(
    raw_df: pd.DataFrame,
    *,
    source_name: str,
    version: str,
    tissue_col: str,
    value_cols: tuple[str, ...],
    qc_set: frozenset[str],
    cleanup_set: frozenset[str],
    norm: str | None,
    brain_subregions: frozenset[str] = frozenset(),
) -> pd.DataFrame:
    import numpy as np

    from .normalization import clean_tpm

    df = raw_df.copy()
    if "duplicate_gtex_tissues" in qc_set:
        df = df[
            ~df[tissue_col].astype(str).str.strip().str.lower().isin(HPA_DUPLICATE_GTEX_TISSUES)
        ].reset_index(drop=True)
    if df.empty or (not (qc_set - {"duplicate_gtex_tissues"}) and not cleanup_set and norm is None):
        return df

    if norm == "clean_tpm" and "TPM" not in df.columns:
        raise ValueError(
            f"normalize='clean_tpm' requires linear 'TPM', not non-linear consensus 'nTPM' ({source_name})"
        )

    if (
        source_name == "hpa_rna_consensus"
        and cleanup_set
        and "duplicate_gtex_tissues" not in qc_set
        and version == "v23"
        and reference_data.is_cached("rna_tissue_gtex", "v23")
    ):
        from . import gtex

        gtex_cons = gtex.tissue_rna(
            version="v23",
            source="hpa_consensus",
            qc=tuple(sorted(qc_set & {"splash"})) or None,
            cleanup=tuple(sorted(cleanup_set)),
            normalize=None,
            use_cache=False,
        )
        # Update ONLY overlapping tissues present in raw_df so HPA-unique surgical tissues stay intact!
        gtex_cons_idx = gtex_cons.set_index(["Gene", tissue_col])["nTPM"]
        key_idx = pd.MultiIndex.from_frame(df[["Gene", tissue_col]])
        repl = gtex_cons_idx.reindex(key_idx).to_numpy()
        has_repl = ~np.isnan(repl)
        df = df.copy()
        df.loc[has_repl, "nTPM"] = repl[has_repl]

    df["_unversioned_gid"] = (
        df["Gene"].astype(str).str.replace(r"\.\d+(?=_PAR_Y$|$)", "", regex=True)
    )
    all_gids = pd.Index(pd.unique(df["_unversioned_gid"]), name="Gene")
    all_tissues = pd.Index(pd.unique(df[tissue_col].astype(str)), name=tissue_col)
    gene_names = (
        df.drop_duplicates("_unversioned_gid").set_index("_unversioned_gid")["Gene name"].to_dict()
    )
    panels, prof_df = _human_marker_panels_and_profiles()
    native_map = {
        str(col): _resolve_hpa_native_panels(str(col), brain_subregions=brain_subregions)
        for col in all_tissues
    }

    present_vcols = [c for c in value_cols if c in df.columns]
    raw_w_map: dict[str, pd.DataFrame] = {}
    for vcol in present_vcols:
        df[vcol] = pd.to_numeric(df[vcol], errors="coerce")
        piv = df.pivot_table(
            index="_unversioned_gid",
            columns=tissue_col,
            values=vcol,
            aggfunc="max",
            dropna=False,
        ).reindex(index=all_gids, columns=all_tissues)
        piv.columns = [str(c) for c in piv.columns]
        raw_w_map[vcol] = piv

    primary_col = present_vcols[0]
    splash_nan_mask: np.ndarray | None = None
    if "splash" in qc_set:
        prim_before = raw_w_map[primary_col]
        prim_masked = bulk_qc.mask_panel_splash_cells(
            prim_before, panels, sample_native_panels=native_map, full_profiles=prof_df
        )
        splash_nan_mask = (
            (prim_before.to_numpy(dtype=float) > 0.0)
            & prim_before.notna().to_numpy()
            & prim_masked.isna().to_numpy()
        )

    w_map: dict[str, pd.DataFrame] = {}
    for vcol in present_vcols:
        piv = raw_w_map[vcol].copy()
        if splash_nan_mask is not None and splash_nan_mask.any():
            arr = piv.to_numpy(dtype=float, copy=True)
            arr[splash_nan_mask] = np.nan
            piv = pd.DataFrame(arr, index=piv.index, columns=piv.columns)
        if "index_hopping" in cleanup_set:
            piv = bulk_qc.subtract_panel_index_hopping(
                piv, panels, sample_native_panels=native_map, full_profiles=prof_df
            )
        if "deconvolution" in cleanup_set:
            piv = bulk_qc.deconvolve_panel_admixture(
                piv, panels, sample_native_panels=native_map, full_profiles=prof_df
            )
        w_map[vcol] = piv

    if norm == "clean_tpm":
        gene_tbl = pd.DataFrame(
            {
                "Ensembl_Gene_ID": all_gids,
                "Symbol": [gene_names.get(str(g), str(g)) for g in all_gids],
            },
            index=all_gids,
        )
        prim_pre = w_map["TPM"]
        prim_clean = clean_tpm(prim_pre, gene_table=gene_tbl, on_empty_compartment="record")
        w_map["TPM"] = prim_clean
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio_scale = np.where(
                prim_pre.to_numpy(dtype=float) > 0.0,
                prim_clean.to_numpy(dtype=float) / prim_pre.to_numpy(dtype=float),
                1.0,
            )
        for vcol in present_vcols:
            if vcol == "TPM":
                continue
            scaled = w_map[vcol].to_numpy(dtype=float) * ratio_scale
            w_map[vcol] = pd.DataFrame(scaled, index=all_gids, columns=w_map[vcol].columns)

    r_idx = all_gids.get_indexer(df["_unversioned_gid"])
    c_idx = pd.Index([str(c) for c in all_tissues]).get_indexer(df[tissue_col].astype(str))
    dup_mask = df.duplicated(["_unversioned_gid", tissue_col], keep=False).to_numpy()
    for vcol in present_vcols:
        clean_arr = w_map[vcol].to_numpy(dtype=float)[r_idx, c_idx]
        if dup_mask.any():
            raw_sum_piv = df.pivot_table(
                index="_unversioned_gid",
                columns=tissue_col,
                values=vcol,
                aggfunc="sum",
                dropna=False,
            ).reindex(index=all_gids, columns=all_tissues)
            raw_sum_arr = raw_sum_piv.to_numpy(dtype=float)[r_idx, c_idx]
            orig_vals = df[vcol].to_numpy(dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                prop = np.where(raw_sum_arr > 0.0, orig_vals / raw_sum_arr, 0.0)
            clean_arr = np.where(dup_mask, clean_arr * prop, clean_arr)
        df[vcol] = clean_arr
    return df.drop(columns=["_unversioned_gid"]).reset_index(drop=True)


def hpa_rna_consensus(
    version: str | None = None,
    *,
    qc=None,
    cleanup=None,
    normalize: str | None = None,
    use_cache: bool = True,
) -> pd.DataFrame:
    """HPA tissue nTPM; supports ``(qc=, cleanup=, normalize=)`` with ``'all'`` aliases (#780)."""
    concrete = reference_data.resolve_version("hpa_rna_consensus", version)
    qc_set = bulk_qc.expand_ops(
        qc,
        allowed=ALL_HPA_BULK_QC,
        aliases=HPA_BULK_QC_ALIASES,
        param_name="qc",
        domain_label="HPA bulk",
    )
    cleanup_set = bulk_qc.expand_ops(
        cleanup,
        allowed=ALL_HPA_BULK_CLEANUP,
        aliases=HPA_BULK_CLEANUP_ALIASES,
        param_name="cleanup",
        domain_label="HPA bulk",
    )
    norm = bulk_qc.resolve_normalize_option(normalize, domain_label="HPA bulk")
    if not qc_set and not cleanup_set and norm is None:
        return _hpa_rna_consensus_for_version(concrete)

    spec = bulk_qc.BulkProcessingSpec(
        domain="hpa_bulk",
        source_id="hpa_rna_consensus",
        source_version=concrete,
        qc=qc_set,
        cleanup=cleanup_set,
        normalize=norm,
    )
    cache_dir = reference_data.local_path("hpa_rna_consensus", concrete).parent / "bulk_processed"
    return bulk_qc.read_or_compute_cached_parquet(
        cache_dir,
        spec.cache_key(),
        spec.manifest(),
        lambda: _process_hpa_long_bulk_table(
            _hpa_rna_consensus_for_version(concrete),
            source_name="hpa_rna_consensus",
            version=concrete,
            tissue_col="Tissue",
            value_cols=("nTPM",),
            qc_set=qc_set,
            cleanup_set=cleanup_set,
            norm=norm,
        ),
        use_cache=use_cache,
        spec_attr_name="bulk_spec",
    )


def hpa_rna_tissue(
    version: str | None = None,
    *,
    include_brain: bool = False,
    qc=None,
    cleanup=None,
    normalize: str | None = None,
    use_cache: bool = True,
) -> pd.DataFrame:
    """HPA surgical normal-tissue bulk RNA (``rna_tissue_hpa``), with ``(qc=, cleanup=, normalize=)``."""
    concrete = reference_data.resolve_version("rna_tissue_hpa", version)
    qc_set = bulk_qc.expand_ops(
        qc,
        allowed=ALL_HPA_BULK_QC,
        aliases=HPA_BULK_QC_ALIASES,
        param_name="qc",
        domain_label="HPA surgical bulk",
    )
    cleanup_set = bulk_qc.expand_ops(
        cleanup,
        allowed=ALL_HPA_BULK_CLEANUP,
        aliases=HPA_BULK_CLEANUP_ALIASES,
        param_name="cleanup",
        domain_label="HPA surgical bulk",
    )
    norm = bulk_qc.resolve_normalize_option(normalize, domain_label="HPA surgical bulk")

    def _raw() -> tuple[pd.DataFrame, frozenset[str]]:
        base = _read_hpa("rna_tissue_hpa", concrete).copy()
        brain_subs: frozenset[str] = frozenset()
        if include_brain:
            brain_ver = reference_data.resolve_version("rna_brain_hpa", version)
            brain = _read_hpa("rna_brain_hpa", brain_ver).copy()
            if "Subregion" in brain.columns and "Tissue" not in brain.columns:
                brain = brain.rename(columns={"Subregion": "Tissue"})
            brain_subs = frozenset(brain["Tissue"].dropna().astype(str).str.strip().str.lower())
            keep_cols = [c for c in base.columns if c in brain.columns]
            base = (
                pd.concat([base[keep_cols], brain[keep_cols]], ignore_index=True)
                .drop_duplicates(["Gene", "Tissue"], keep="first")
                .reset_index(drop=True)
            )
        return base, brain_subs

    if not qc_set and not cleanup_set and norm is None:
        return _raw()[0]

    spec = bulk_qc.BulkProcessingSpec(
        domain="hpa_surgical_bulk",
        source_id="rna_tissue_hpa",
        source_version=concrete,
        qc=qc_set,
        cleanup=cleanup_set,
        normalize=norm,
        extra=(("include_brain", bool(include_brain)),),
    )
    cache_dir = reference_data.local_path("rna_tissue_hpa", concrete).parent / "bulk_processed"

    def _build_surgical() -> pd.DataFrame:
        raw_tbl, b_subs = _raw()
        return _process_hpa_long_bulk_table(
            raw_tbl,
            source_name="rna_tissue_hpa",
            version=concrete,
            tissue_col="Tissue",
            value_cols=("TPM", "pTPM", "nTPM"),
            qc_set=qc_set,
            cleanup_set=cleanup_set,
            norm=norm,
            brain_subregions=b_subs,
        )

    return bulk_qc.read_or_compute_cached_parquet(
        cache_dir,
        spec.cache_key(),
        spec.manifest(),
        _build_surgical,
        use_cache=use_cache,
        spec_attr_name="bulk_spec",
    )


@cache
def _hpa_normal_tissue_for_version(version: str) -> pd.DataFrame:
    """Cache one HPA IHC table under its concrete source-version key."""
    return _read_hpa("hpa_normal_tissue", version)


def hpa_normal_tissue(version: str | None = None) -> pd.DataFrame:
    """HPA IHC protein detection per tissue/cell type for one source version."""
    concrete_version = reference_data.resolve_version("hpa_normal_tissue", version)
    return _hpa_normal_tissue_for_version(concrete_version)


@cache
def _hpa_single_cell_for_version(version: str) -> pd.DataFrame:
    return _read_hpa("hpa_single_cell", version)


def _single_cell_value_col(df: pd.DataFrame) -> str:
    """Resolve the HPA single-cell expression column ('nCPM' in v25.1+, 'nTPM' in v23)."""
    for col in ("nCPM", "nTPM"):
        if col in df.columns:
            return col
    raise KeyError(
        f"HPA single-cell table lacks an 'nCPM' or 'nTPM' column; got {list(df.columns)}"
    )


def hpa_single_cell(version: str | None = None) -> pd.DataFrame:
    """HPA cell-type RNA (``nCPM`` in v25.1+, ``nTPM`` in v23); newest release unless pinned."""
    return _hpa_single_cell_for_version(reference_data.resolve_version("hpa_single_cell", version))


@lru_cache(maxsize=4)
def _hpa_cell_type_expression_for_version(version: str) -> pd.DataFrame:
    df = hpa_single_cell(version)
    if df.empty:
        out = pd.DataFrame(columns=["Ensembl_Gene_ID", "Symbol"])
        out.attrs.update(df.attrs)
        return out
    value_col = _single_cell_value_col(df)
    wide = df.pivot_table(
        index=["Gene", "Gene name"],
        columns="Cell type",
        values=value_col,
        aggfunc="sum",
        fill_value=0.0,
    )
    wide = wide.reset_index().rename(columns={"Gene": "Ensembl_Gene_ID", "Gene name": "Symbol"})
    wide.columns.name = None
    wide.attrs.update(df.attrs, value_column=value_col)
    return wide


def hpa_cell_type_expression(version: str | None = None) -> pd.DataFrame:
    """Wide per-gene HPA single-cell-type expression matrix (``nCPM`` in v25.1+, ``nTPM`` in v23).

    Returns ``Ensembl_Gene_ID``, ``Symbol``, then one numeric column per HPA
    single-cell type. This is the analysis-facing companion to the raw long-form
    :func:`hpa_single_cell` table.
    """
    concrete = reference_data.resolve_version("hpa_single_cell", version)
    cached = _hpa_cell_type_expression_for_version(concrete)
    out = cached.copy()
    out.attrs.update(cached.attrs)
    return out


hpa_cell_type_expression.cache_clear = _hpa_cell_type_expression_for_version.cache_clear  # type: ignore[attr-defined]
hpa_cell_type_expression.cache_info = _hpa_cell_type_expression_for_version.cache_info  # type: ignore[attr-defined]
_register_derived_cache(_hpa_cell_type_expression_for_version.cache_clear)


@lru_cache(maxsize=1)
def _safety_tissue_mapping_table() -> pd.DataFrame:
    table = get_data("hpa-safety-tissue-map", copy=False).fillna("")
    missing = set(_SAFETY_TISSUE_MAPPING_COLUMNS) - set(table.columns)
    if missing:
        raise SafetyTissueResolutionError(
            f"HPA safety-tissue mapping table is missing columns: {sorted(missing)}"
        )
    table = table.loc[:, _SAFETY_TISSUE_MAPPING_COLUMNS].astype(str)

    invalid_kinds = set(table["mapping_kind"]) - SAFETY_TISSUE_MAPPING_KINDS
    invalid_levels = set(table["coverage_level"]) - SAFETY_TISSUE_COVERAGE_LEVELS
    if invalid_kinds or invalid_levels:
        raise SafetyTissueResolutionError(
            "invalid HPA safety-tissue mapping vocabulary: "
            f"mapping_kind={sorted(invalid_kinds)}, coverage_level={sorted(invalid_levels)}"
        )

    duplicate_rows = table.duplicated(
        ["source_name", "source_version", "safety_group", "requested_tissue", "source_tissue"]
    )
    if duplicate_rows.any():
        raise SafetyTissueResolutionError("duplicate HPA safety-tissue mapping rows")

    unavailable = table["coverage_level"] == "unavailable"
    invalid_availability = (unavailable & table["source_tissue"].ne("")) | (
        ~unavailable & table["source_tissue"].eq("")
    )
    if invalid_availability.any():
        raise SafetyTissueResolutionError(
            "unavailable mappings must omit source_tissue and available mappings must provide it"
        )

    for (source_name, source_version), source_rows in table.groupby(
        ["source_name", "source_version"], sort=True
    ):
        mapped_groups = set(source_rows["safety_group"])
        configured_groups = set(SAFETY_TISSUE_GROUPS)
        if mapped_groups != configured_groups:
            raise SafetyTissueResolutionError(
                f"{source_name!r} {source_version!r} maps safety groups "
                f"{sorted(mapped_groups)}, expected {sorted(configured_groups)}"
            )
        for group, configured_tissues in SAFETY_TISSUE_GROUPS.items():
            mapped_tissues = set(
                source_rows.loc[source_rows["safety_group"] == group, "requested_tissue"]
            )
            if mapped_tissues != configured_tissues:
                raise SafetyTissueResolutionError(
                    f"{source_name!r} {source_version!r} group {group!r} maps "
                    f"{sorted(mapped_tissues)}, expected {sorted(configured_tissues)}"
                )
    return table


_register_derived_cache(_safety_tissue_mapping_table.cache_clear)


def safety_tissue_mapping_table() -> pd.DataFrame:
    """Reviewed conceptual-to-source-native safety-tissue mappings.

    The returned table is a defensive copy. Use :func:`resolve_safety_tissue_group`
    for a validated, immutable selection contract.
    """
    return _safety_tissue_mapping_table().copy()


def resolve_safety_tissue_group(
    safety_group: str,
    *,
    source_name: str = "hpa_normal_tissue",
    source_version: str | None = None,
    require_complete: bool = True,
) -> SafetyTissueResolution:
    """Resolve a conceptual safety group to exact source-native tissue labels.

    Resolution is pinned to a concrete ``(source_name, source_version)`` and
    includes reviewed aliases/substructures and unavailable conceptual tissues.
    By default incomplete coverage raises instead of silently weakening a safety
    filter. Set ``require_complete=False`` only when the returned partial and
    unavailable coverage has been explicitly reviewed by the caller.
    """
    group = str(safety_group).strip().casefold()
    if group not in SAFETY_TISSUE_GROUPS:
        raise SafetyTissueResolutionError(
            f"unknown safety tissue group {safety_group!r}; "
            f"available: {', '.join(sorted(SAFETY_TISSUE_GROUPS))}"
        )
    try:
        version = reference_data.resolve_version(source_name, source_version)
    except reference_data.ReferenceDataError as exc:
        raise SafetyTissueResolutionError(str(exc)) from None

    table = _safety_tissue_mapping_table()
    rows = table[
        (table["source_name"] == source_name)
        & (table["source_version"] == version)
        & (table["safety_group"] == group)
    ]
    if rows.empty:
        supported = sorted(
            {
                (str(name), str(source_version))
                for name, source_version in zip(table["source_name"], table["source_version"])
            }
        )
        supported_text = ", ".join(f"{name}@{ver}" for name, ver in supported)
        raise SafetyTissueResolutionError(
            f"no safety-tissue mapping for {source_name!r} {version!r}; supported: {supported_text}"
        )

    mappings = []
    for requested_tissue, tissue_rows in rows.groupby("requested_tissue", sort=True):
        kinds = set(tissue_rows["mapping_kind"])
        levels = set(tissue_rows["coverage_level"])
        if len(kinds) != 1 or len(levels) != 1:
            raise SafetyTissueResolutionError(
                f"inconsistent mapping rows for {group!r}/{requested_tissue!r}"
            )
        mappings.append(
            SafetyTissueMapping(
                requested_tissue=str(requested_tissue),
                source_tissues=tuple(sorted(t for t in tissue_rows["source_tissue"] if t)),
                mapping_kind=next(iter(kinds)),
                coverage_level=next(iter(levels)),
                mapping_references=tuple(
                    sorted({ref for ref in tissue_rows["mapping_reference"] if ref})
                ),
                notes=tuple(sorted({note for note in tissue_rows["notes"] if note})),
            )
        )

    coverage_levels = {mapping.coverage_level for mapping in mappings}
    if coverage_levels == {"complete"}:
        coverage_state = "complete"
    elif coverage_levels == {"unavailable"}:
        coverage_state = "unavailable"
    else:
        coverage_state = "partial"
    if coverage_state not in SAFETY_TISSUE_COVERAGE_STATES:
        raise SafetyTissueResolutionError(f"invalid coverage state {coverage_state!r}")

    resolution = SafetyTissueResolution(
        safety_group=group,
        source_name=source_name,
        source_version=version,
        source_url=reference_data.provenance(source_name, version)["url"],
        mappings=tuple(mappings),
        coverage_state=coverage_state,
    )
    if require_complete and not resolution.complete:
        raise SafetyTissueResolutionError(
            f"safety tissue group {group!r} is {coverage_state} for "
            f"{source_name}@{version}; partial={list(resolution.partially_covered_tissues)}, "
            f"unavailable={list(resolution.unavailable_tissues)}. Pass "
            "require_complete=False only after reviewing this coverage."
        )
    return resolution


@cache
def _hpa_normal_tissue_labels_for_version(version: str) -> tuple[str, ...]:
    labels = _hpa_normal_tissue_for_version(version)["Tissue"].dropna().astype(str)
    return tuple(sorted({label.strip() for label in labels if label.strip()}))


def hpa_normal_tissue_labels(version: str | None = None) -> tuple[str, ...]:
    """Exact HPA IHC tissue labels available in one concrete source version."""
    concrete_version = reference_data.resolve_version("hpa_normal_tissue", version)
    return _hpa_normal_tissue_labels_for_version(concrete_version)


def resolve_hpa_normal_tissue_label(label: str, *, version: str | None = None) -> str:
    """Return the exact source-native HPA IHC label or reject an unavailable label.

    Resolving the label independently of a gene query distinguishes a recognized
    tissue with no observation for that gene from a label absent from the source.
    """
    concrete_version = reference_data.resolve_version("hpa_normal_tissue", version)
    native_labels = hpa_normal_tissue_labels(concrete_version)
    normalized = str(label).strip().casefold()
    matches = [
        native_label for native_label in native_labels if native_label.casefold() == normalized
    ]
    if len(matches) == 1:
        return matches[0]
    raise SafetyTissueResolutionError(
        f"{label!r} is not an HPA IHC tissue label in hpa_normal_tissue@{concrete_version}"
    )


def _strip_version(gene_id: str) -> str:
    return str(gene_id).split(".")[0]


def gene_tissue_ntpm(gene_id: str) -> dict[str, float]:
    """``{tissue (lowercased): nTPM}`` of normal RNA expression for one gene
    (unversioned Ensembl ID)."""
    gid = _strip_version(gene_id)
    df = hpa_rna_consensus()
    sub = df[df["Gene"].astype(str).map(_strip_version) == gid]
    return {str(t).strip().lower(): float(v) for t, v in zip(sub["Tissue"], sub["nTPM"])}


def gene_cell_type_ntpm(gene_id: str, *, version: str | None = None) -> dict[str, float]:
    """``{cell_type (lowercased): value}`` of HPA single-cell RNA for one gene (``nCPM`` in v25.1+, ``nTPM`` in v23)."""
    gid = _strip_version(gene_id)
    df = hpa_single_cell(version)
    if df.empty:
        return {}
    value_col = _single_cell_value_col(df)
    sub = df[df["Gene"].astype(str).map(_strip_version) == gid]
    return {str(t).strip().lower(): float(v) for t, v in zip(sub["Cell type"], sub[value_col])}


def gene_protein_tissues(gene_id: str, *, levels=("Low", "Medium", "High")) -> set[str]:
    """Tissues (lowercased) where the gene has detected IHC protein at one of
    ``levels`` (default Low/Medium/High)."""
    gid = _strip_version(gene_id)
    df = hpa_normal_tissue()
    sub = df[df["Gene"].astype(str).map(_strip_version) == gid]
    sub = sub[sub["Level"].astype(str).isin(levels)]
    return {str(t).strip().lower() for t in sub["Tissue"]}

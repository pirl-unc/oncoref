# Licensed under the Apache License, Version 2.0.
"""Human Protein Atlas normal-tissue, single-cell and cancer evidence.

These are the existing HPA loaders and caches, exposed together under the full
dataset name. Native units, assay limitations and historical imports are retained,
while bulk & tissue-cluster RNA loaders (:func:`rna_consensus`, :func:`rna_tissue_hpa`,
:func:`immune_cells`, :func:`single_cell_tissues`) share the unified
``(qc=, cleanup=, normalize=)`` + ``"all"`` alias contract (#779, #780).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import bulk_qc, reference_data
from .hpa import (
    ALL_HPA_BULK_CLEANUP,
    ALL_HPA_BULK_QC,
    HPA_BULK_CLEANUP_ALIASES,
    HPA_BULK_QC_ALIASES,
    HPA_DUPLICATE_GTEX_TISSUES,
    SAFETY_TISSUE_COVERAGE_LEVELS,
    SAFETY_TISSUE_COVERAGE_STATES,
    SAFETY_TISSUE_MAPPING_KINDS,
    SafetyTissueMapping,
    SafetyTissueResolution,
    SafetyTissueResolutionError,
    gene_cell_type_ntpm,
    gene_protein_tissues,
    gene_tissue_ntpm,
    resolve_safety_tissue_group,
    safety_tissue_mapping_table,
)
from .hpa import (
    hpa_cell_type_expression as cell_type_expression,
)
from .hpa import (
    hpa_normal_tissue as normal_tissue,
)
from .hpa import (
    hpa_normal_tissue_labels as normal_tissue_labels,
)
from .hpa import (
    hpa_rna_consensus as rna_consensus,
)
from .hpa import (
    hpa_rna_tissue as rna_tissue_hpa,
)
from .hpa import (
    hpa_single_cell as single_cell,
)
from .hpa import (
    resolve_hpa_normal_tissue_label as resolve_normal_tissue_label,
)
from .hpa_cancer import (
    hpa_cancer_assay_limitations as cancer_assay_limitations,
)
from .hpa_cancer import (
    hpa_cancer_crosswalk as cancer_crosswalk,
)
from .hpa_cancer import (
    hpa_cancer_ihc_prevalence as cancer_ihc_prevalence,
)
from .hpa_cancer import (
    hpa_cancer_rna_cohorts as cancer_rna_cohorts,
)
from .hpa_cancer import (
    hpa_cancer_rna_ihc_comparison as cancer_rna_ihc_comparison,
)
from .hpa_cancer import (
    hpa_cancer_rna_prevalence as cancer_rna_prevalence,
)
from .hpa_cancer import (
    hpa_cancer_sources as cancer_sources,
)
from .normal_cell_data import _gene_filter, observations, table

ALL_IMMUNE_BULK_QC: tuple[str, ...] = ("ex_vivo_stimulated", "low_purity_sort")
IMMUNE_BULK_QC_ALIASES: dict[str, frozenset[str]] = {
    "all": frozenset(ALL_IMMUNE_BULK_QC),
    "none": frozenset(),
    "raw": frozenset(),
}

ALL_IMMUNE_BULK_CLEANUP: tuple[str, ...] = ("sort_carryover", "index_hopping")
IMMUNE_BULK_CLEANUP_ALIASES: dict[str, frozenset[str]] = {
    "all": frozenset(ALL_IMMUNE_BULK_CLEANUP),
    "none": frozenset(),
    "raw": frozenset(),
}

ALL_SC_TISSUE_QC: tuple[str, ...] = ("misannotated_erythroid",)
SC_TISSUE_QC_ALIASES: dict[str, frozenset[str]] = {
    "all": frozenset(ALL_SC_TISSUE_QC),
    "none": frozenset(),
    "raw": frozenset(),
}

ALL_SC_TISSUE_CLEANUP: tuple[str, ...] = ("ambient_soup",)
SC_TISSUE_CLEANUP_ALIASES: dict[str, frozenset[str]] = {
    "all": frozenset(ALL_SC_TISSUE_CLEANUP),
    "none": frozenset(),
    "raw": frozenset(),
}

_IMMUNE_HOP_PANELS: dict[str, dict[str, float]] = {
    "b_lineage": {
        "ENSG00000177455": 650.0,  # CD19
        "ENSG00000156738": 1400.0,  # MS4A1
        "ENSG00000105369": 950.0,  # CD79A
        "ENSG00000007312": 900.0,  # CD79B
    },
    "t_lineage": {
        "ENSG00000167286": 1200.0,  # CD3D
        "ENSG00000198851": 1100.0,  # CD3E
        "ENSG00000153563": 800.0,  # CD8A
    },
    "nk_lineage": {
        "ENSG00000149294": 450.0,  # NCAM1
        "ENSG00000134539": 900.0,  # KLRD1
    },
    "myeloid_granulocyte": {
        "ENSG00000143546": 8000.0,  # S100A8
        "ENSG00000163220": 9500.0,  # S100A9
        "ENSG00000119535": 900.0,  # CSF3R
    },
}

_UNFRACTIONATED_IMMUNE_MIXTURES: frozenset[str] = frozenset(
    {"total pbmc", "pbmc", "whole blood", "unfractionated pbmc"}
)

_SOUP_RECEIVER_KEYWORDS: tuple[str, ...] = (
    "macrophage",
    "t-cell",
    "b-cell",
    "plasma cell",
    "nk-cell",
    "monocyte",
    "granulocyte",
    "dendritic",
    "mast cell",
    "fibroblast",
    "endothelial",
    "adipocyte",
    "smooth muscle",
    "peritubular",
    "hofbauer",
)


def _native_immune_panels(cell_label: str) -> frozenset[str]:
    low = cell_label.strip().lower()
    native: set[str] = set()
    if "b-cell" in low or "b cell" in low or "plasmablast" in low:
        native.add("b_lineage")
    if "t-cell" in low or "t cell" in low or "t-reg" in low or "gdt" in low or "mait" in low:
        native.add("t_lineage")
    if "nk-cell" in low or "nk cell" in low:
        native.add("nk_lineage")
    if any(
        tok in low for tok in ("monocyte", "neutrophil", "basophil", "eosinophil", "myeloid", " dc")
    ) or low.endswith("dc"):
        native.add("myeloid_granulocyte")
    if "pbmc" in low or "progenitor" in low:
        native.update({"b_lineage", "t_lineage", "nk_lineage", "myeloid_granulocyte"})
    return frozenset(native)


def _is_ex_vivo_stimulated_label(label: str) -> bool:
    low = label.strip().lower()
    return any(tok in low for tok in ("activated", "stimulated", "_stim", "pma", "lps"))


def _process_immune_bulk_table(
    raw_df: pd.DataFrame,
    *,
    qc_set: frozenset[str],
    cleanup_set: frozenset[str],
    norm: str | None,
) -> pd.DataFrame:
    from .normalization import clean_tpm

    df = raw_df.copy()
    cell_col = "Immune cell" if "Immune cell" in df.columns else "Cell type"
    val_cols = [c for c in ("TPM", "pTPM", "nTPM") if c in df.columns]
    if df.empty or not val_cols:
        return df

    if "ex_vivo_stimulated" in qc_set:
        df = df[~df[cell_col].astype(str).map(_is_ex_vivo_stimulated_label)].reset_index(drop=True)
    if "low_purity_sort" in qc_set:
        # Reject unfractionated mixture pools (e.g. 'total PBMC') and gross >= 4% platelet
        # sort failures; never drop FACS-sorted leukocyte lineages over co-isolated RBC/globin
        # transcripts (e.g. Monaco B-cells), which are handled continuously by sort_carryover.
        df = df[
            ~df[cell_col].astype(str).str.strip().str.lower().isin(_UNFRACTIONATED_IMMUNE_MIXTURES)
        ].reset_index(drop=True)
    if df.empty:
        return df

    if norm == "clean_tpm" and "TPM" not in val_cols:
        raise ValueError(
            "normalize='clean_tpm' requires linear 'TPM', not non-linear 'nTPM' in immune bulk"
        )

    df["_unversioned_gid"] = (
        df["Gene"].astype(str).str.replace(r"\.\d+(?=_PAR_Y$|$)", "", regex=True)
    )
    all_gids = pd.Index(pd.unique(df["_unversioned_gid"]), name="Gene")
    keep_cols = list(pd.unique(df[cell_col].astype(str)))
    gene_names = (
        df.drop_duplicates("_unversioned_gid").set_index("_unversioned_gid")["Gene name"].to_dict()
        if "Gene name" in df.columns
        else {}
    )
    carryover_panels = {
        k: pd.Series(v, dtype=float) for k, v in bulk_qc.SORTED_IMMUNE_CARRYOVER_PANELS.items()
    }
    carryover_profiles = {
        k: pd.Series(v, dtype=float) for k, v in bulk_qc.SORTED_IMMUNE_CARRYOVER_PROFILES.items()
    }
    hop_panels = {k: pd.Series(v, dtype=float) for k, v in _IMMUNE_HOP_PANELS.items()}
    # Progenitors possess native megakaryocyte-erythroid lineage potential (ITGA2B/ALAS2)
    carryover_native = {
        c: frozenset({"platelets", "erythroid"}) if "progenitor" in c.lower() else frozenset()
        for c in keep_cols
    }

    primary_vcol = val_cols[0]
    raw_w: dict[str, pd.DataFrame] = {}
    for vcol in val_cols:
        df[vcol] = pd.to_numeric(df[vcol], errors="coerce")
        piv = df.pivot_table(
            index="_unversioned_gid",
            columns=cell_col,
            values=vcol,
            aggfunc="max",
            dropna=False,
        ).reindex(index=all_gids, columns=keep_cols)
        piv.columns = [str(c) for c in piv.columns]
        raw_w[vcol] = piv

    if "low_purity_sort" in qc_set and keep_cols:
        plt_fracs = bulk_qc.estimate_foreign_panel_fractions(
            raw_w[primary_vcol],
            {"platelets": carryover_panels["platelets"]},
            sample_native_panels=carryover_native,
        )
        keep_cols = [c for c in keep_cols if float(plt_fracs.loc["platelets", c]) < 0.04]
        if not keep_cols:
            out_cols = [
                c for c in ("Gene", "Gene name", cell_col, *val_cols) if c in raw_df.columns
            ]
            return raw_df.loc[:, out_cols].iloc[0:0].copy()
        for vcol in val_cols:
            raw_w[vcol] = raw_w[vcol].loc[:, keep_cols]
        carryover_native = {c: carryover_native[c] for c in keep_cols}

    w_map: dict[str, pd.DataFrame] = {}
    for vcol in val_cols:
        piv = raw_w[vcol].copy()
        if "sort_carryover" in cleanup_set and not piv.empty:
            c_fracs = bulk_qc.estimate_foreign_panel_fractions(
                piv, carryover_panels, sample_native_panels=carryover_native
            ).clip(upper=0.60)
            prof_c = (
                pd.DataFrame(carryover_profiles)
                .reindex(index=piv.index, columns=c_fracs.index)
                .fillna(0.0)
            )
            sub_c = prof_c.to_numpy(dtype=float) @ c_fracs.to_numpy(dtype=float)
            arr = piv.to_numpy(dtype=float, copy=True)
            fin = np.isfinite(arr)
            arr[fin] = np.fmax(0.0, arr[fin] - sub_c[fin])
            piv = pd.DataFrame(arr, index=piv.index, columns=piv.columns)
        if "index_hopping" in cleanup_set and not piv.empty:
            native_map = {c: _native_immune_panels(c) for c in piv.columns}
            piv = bulk_qc.subtract_panel_index_hopping(
                piv, hop_panels, sample_native_panels=native_map
            )
        w_map[vcol] = piv

    if norm == "clean_tpm" and keep_cols:
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
        for vcol in val_cols:
            if vcol == "TPM":
                continue
            scaled = w_map[vcol].to_numpy(dtype=float) * ratio_scale
            w_map[vcol] = pd.DataFrame(scaled, index=all_gids, columns=w_map[vcol].columns)

    df = df[df[cell_col].astype(str).isin(set(keep_cols))].reset_index(drop=True)
    r_idx = all_gids.get_indexer(df["_unversioned_gid"])
    c_idx = pd.Index([str(c) for c in keep_cols]).get_indexer(df[cell_col].astype(str))
    dup_mask = df.duplicated(["_unversioned_gid", cell_col], keep=False).to_numpy()
    for vcol in val_cols:
        clean_arr = w_map[vcol].to_numpy(dtype=float)[r_idx, c_idx]
        if dup_mask.any():
            raw_sum_piv = df.pivot_table(
                index="_unversioned_gid",
                columns=cell_col,
                values=vcol,
                aggfunc="sum",
                dropna=False,
            ).reindex(index=all_gids, columns=keep_cols)
            raw_sum_arr = raw_sum_piv.to_numpy(dtype=float)[r_idx, c_idx]
            orig_vals = df[vcol].to_numpy(dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                prop = np.where(raw_sum_arr > 0.0, orig_vals / raw_sum_arr, 0.0)
            clean_arr = np.where(dup_mask, clean_arr * prop, clean_arr)
        df[vcol] = clean_arr
    return df.drop(columns=["_unversioned_gid"]).reset_index(drop=True)


def _process_single_cell_tissue_table(
    raw_df: pd.DataFrame,
    *,
    qc_set: frozenset[str],
    cleanup_set: frozenset[str],
    soup_alpha: float = 0.004,
    dominator_min_ntpm: float = 50000.0,
) -> pd.DataFrame:
    """Apply cluster QC and droplet ambient-soup subtraction on ``hpa_single_cell_tissue`` (#779, #780).

    - ``qc='misannotated_erythroid'``: drops non-erythroid clusters whose hemoglobin
      transcripts (``HBB + HBA1 + HBA2``) exceed ``>= 100,000 nTPM`` (e.g. ``breast c-24``
      annotated as ``t-cells`` at ``HBB = 597,120 nTPM``).
    - ``cleanup='ambient_soup'``: within each tissue dissociation, identifies hyper-abundant
      parenchymal secretory dominators (``peak nTPM >= 50,000``) and subtracts proportional
      cell-free soup bleed (``soup_alpha * tissue_peak_ntpm``, default ``0.4%``) from
      co-captured immune, stromal, endothelial, and adipocyte clusters in that same tissue.
    """
    df = raw_df.copy()
    if df.empty:
        return df
    df["nTPM"] = pd.to_numeric(df["nTPM"], errors="coerce")
    gids = df["Gene"].astype(str).str.replace(r"\.\d+(?=_PAR_Y$|$)", "", regex=True)

    if "misannotated_erythroid" in qc_set:
        hb_rows = df[gids.isin({"ENSG00000244734", "ENSG00000206172", "ENSG00000188536"})]
        if not hb_rows.empty:
            hb_sum = hb_rows.groupby(["Tissue", "Cluster", "Cell type"])["nTPM"].sum()
            bad_clusters = {
                (t, cl)
                for (t, cl, ct), total in hb_sum.items()
                if float(total) >= 100000.0 and "erythroid" not in str(ct).lower()
            }
            if bad_clusters:
                keep_mask = [
                    (t, cl) not in bad_clusters for t, cl in zip(df["Tissue"], df["Cluster"])
                ]
                df = df.loc[keep_mask].reset_index(drop=True)

    if "ambient_soup" in cleanup_set and not df.empty:
        peaks = df.groupby(["Gene", "Tissue"])["nTPM"].max()
        dom_peaks = peaks[peaks >= dominator_min_ntpm]
        if not dom_peaks.empty:
            dom_dict = {(str(g), str(t)): float(v) for (g, t), v in dom_peaks.items()}
            is_receiver = (
                df["Cell type"]
                .astype(str)
                .str.lower()
                .map(lambda ct: any(k in ct for k in _SOUP_RECEIVER_KEYWORDS))
                .to_numpy()
            )
            vals = df["nTPM"].to_numpy(dtype=float, copy=True)
            for idx, (g, t, recv) in enumerate(zip(df["Gene"], df["Tissue"], is_receiver)):
                if not recv:
                    continue
                pk = dom_dict.get((str(g), str(t)))
                if pk is not None and np.isfinite(vals[idx]) and vals[idx] < 0.05 * pk:
                    vals[idx] = max(0.0, vals[idx] - soup_alpha * pk)
            df["nTPM"] = vals

    return df


def primary_ihc(*, genes=None, version=None):
    """Primary antibody observations, including uncertain profiles omitted by curated IHC."""
    return table("hpa_antibody_xml", genes=genes, version=version)


def immune_cells(
    study: str = "hpa",
    *,
    genes=None,
    version=None,
    qc=None,
    cleanup=None,
    normalize: str | None = None,
    use_cache: bool = True,
):
    """Native HPA, Monaco or DICE/Schmiedel immune RNA; supports ``(qc=, cleanup=, normalize=)`` (#780)."""
    names = {
        "hpa": "hpa_immune_cell",
        "monaco": "hpa_immune_monaco",
        "dice": "hpa_immune_schmiedel",
        "schmiedel": "hpa_immune_schmiedel",
    }
    if study not in names:
        raise ValueError(f"Unknown immune study: {study}")
    src_name = names[study]
    qc_set = bulk_qc.expand_ops(
        qc,
        allowed=ALL_IMMUNE_BULK_QC,
        aliases=IMMUNE_BULK_QC_ALIASES,
        param_name="qc",
        domain_label="immune bulk",
    )
    cleanup_set = bulk_qc.expand_ops(
        cleanup,
        allowed=ALL_IMMUNE_BULK_CLEANUP,
        aliases=IMMUNE_BULK_CLEANUP_ALIASES,
        param_name="cleanup",
        domain_label="immune bulk",
    )
    norm = bulk_qc.resolve_normalize_option(normalize, domain_label="immune bulk")
    if not qc_set and not cleanup_set and norm is None:
        return table(src_name, genes=genes, version=version)

    concrete = reference_data.resolve_version(src_name, version)
    spec = bulk_qc.BulkProcessingSpec(
        domain="immune_bulk",
        source_id=src_name,
        source_version=concrete,
        qc=qc_set,
        cleanup=cleanup_set,
        normalize=norm,
    )
    cache_dir = reference_data.local_path(src_name, concrete).parent / "bulk_processed"
    processed = bulk_qc.read_or_compute_cached_parquet(
        cache_dir,
        spec.cache_key(),
        spec.manifest(),
        lambda: _process_immune_bulk_table(
            table(src_name, genes=None, version=concrete),
            qc_set=qc_set,
            cleanup_set=cleanup_set,
            norm=norm,
        ),
        use_cache=use_cache,
        spec_attr_name="bulk_spec",
    )
    if genes is not None:
        filtered = _gene_filter(processed, genes).reset_index(drop=True)
        filtered.attrs.update(processed.attrs)
        return filtered
    return processed


def single_cell_tissues(
    *,
    genes=None,
    version=None,
    qc=None,
    cleanup=None,
    normalize: str | None = None,
    use_cache: bool = True,
):
    """Native single-cell RNA with tissue/cluster identity, read counts, & optional ambient-soup QC."""
    qc_set = bulk_qc.expand_ops(
        qc,
        allowed=ALL_SC_TISSUE_QC,
        aliases=SC_TISSUE_QC_ALIASES,
        param_name="qc",
        domain_label="HPA tissue scRNA",
    )
    cleanup_set = bulk_qc.expand_ops(
        cleanup,
        allowed=ALL_SC_TISSUE_CLEANUP,
        aliases=SC_TISSUE_CLEANUP_ALIASES,
        param_name="cleanup",
        domain_label="HPA tissue scRNA",
    )
    norm = bulk_qc.resolve_normalize_option(normalize, domain_label="HPA tissue scRNA")
    if norm is not None:
        raise ValueError(
            "normalize='clean_tpm' requires linear 'TPM', not cluster 'nTPM' (hpa_single_cell_tissue)"
        )
    if not qc_set and not cleanup_set:
        return table("hpa_single_cell_tissue", genes=genes, version=version)

    concrete = reference_data.resolve_version("hpa_single_cell_tissue", version)
    spec = bulk_qc.BulkProcessingSpec(
        domain="hpa_sc_tissue",
        source_id="hpa_single_cell_tissue",
        source_version=concrete,
        qc=qc_set,
        cleanup=cleanup_set,
        normalize=None,
    )
    cache_dir = (
        reference_data.local_path("hpa_single_cell_tissue", concrete).parent / "sc_processed"
    )
    processed = bulk_qc.read_or_compute_cached_parquet(
        cache_dir,
        spec.cache_key(),
        spec.manifest(),
        lambda: _process_single_cell_tissue_table(
            table("hpa_single_cell_tissue", genes=None, version=concrete),
            qc_set=qc_set,
            cleanup_set=cleanup_set,
        ),
        use_cache=use_cache,
        spec_attr_name="sc_spec",
    )
    if genes is not None:
        filtered = _gene_filter(processed, genes).reset_index(drop=True)
        filtered.attrs.update(processed.attrs)
        return filtered
    return processed


__all__ = [
    "ALL_HPA_BULK_CLEANUP",
    "ALL_HPA_BULK_QC",
    "ALL_IMMUNE_BULK_CLEANUP",
    "ALL_IMMUNE_BULK_QC",
    "ALL_SC_TISSUE_CLEANUP",
    "ALL_SC_TISSUE_QC",
    "HPA_BULK_CLEANUP_ALIASES",
    "HPA_BULK_QC_ALIASES",
    "HPA_DUPLICATE_GTEX_TISSUES",
    "IMMUNE_BULK_CLEANUP_ALIASES",
    "IMMUNE_BULK_QC_ALIASES",
    "SAFETY_TISSUE_COVERAGE_LEVELS",
    "SAFETY_TISSUE_COVERAGE_STATES",
    "SAFETY_TISSUE_MAPPING_KINDS",
    "SC_TISSUE_CLEANUP_ALIASES",
    "SC_TISSUE_QC_ALIASES",
    "SafetyTissueMapping",
    "SafetyTissueResolution",
    "SafetyTissueResolutionError",
    "cancer_assay_limitations",
    "cancer_crosswalk",
    "cancer_ihc_prevalence",
    "cancer_rna_cohorts",
    "cancer_rna_ihc_comparison",
    "cancer_rna_prevalence",
    "cancer_sources",
    "cell_type_expression",
    "gene_cell_type_ntpm",
    "gene_protein_tissues",
    "gene_tissue_ntpm",
    "immune_cells",
    "normal_tissue",
    "normal_tissue_labels",
    "observations",
    "primary_ihc",
    "resolve_normal_tissue_label",
    "resolve_safety_tissue_group",
    "rna_consensus",
    "rna_tissue_hpa",
    "safety_tissue_mapping_table",
    "single_cell",
    "single_cell_tissues",
    "table",
]

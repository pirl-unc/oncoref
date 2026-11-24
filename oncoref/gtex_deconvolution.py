# Licensed under the Apache License, Version 2.0.
"""Holistic GTEx multi-tissue deconvolution and purity-weighted donor quantiles (#778).

Models cross-sample contamination across multiplexed GTEx v11 tissue shards across four
physical regimes—patterned-flowcell index hopping, 96-well plate / extraction-batch
pipetting splash (``SM-XXXX`` prep plates, ``SMNABTCH``/``SMGEBTCH`` batches, and
opposite-sex ``XX``/``XY`` impossibility), prosector/same-donor carryover
(``GTEX-XXXXX`` donors), and gross biopsy tube swaps—while preserving histological
lineage sharing, intravascular perfusion, surgical/anatomical boundary co-sampling
(e.g. prostate transcripts in ``bladder``, epidermis in subcutaneous adipose,
gastrocnemius in ``nerve_tibial``, epicardial coronary/aortic wall off myocardium),
and genuine single-gene ectopic expression.

Core-install safe: uses only Python stdlib + ``numpy`` + ``pandas`` (no ``scipy``).
"""

from __future__ import annotations

import copy
import functools
import hashlib
import json
import math
from dataclasses import dataclass
from importlib.resources import files

import numpy as np
import pandas as pd

from .load_dataset import get_data

GTEX_DECONVOLUTION_POLICY_ID = "oncoref.gtex-holistic-deconvolution.2026-11-24.3"

# Male-specific Y-chromosome sentinel loci (RPS4Y1, DDX3Y, UTY, KDM5D) for donor/sample
# chromosomal sex inference (XX vs XY).
Y_CHROMOSOME_SENTINEL_GENES: tuple[str, ...] = (
    "ENSG00000129824",
    "ENSG00000067048",
    "ENSG00000183878",
    "ENSG00000012817",
)

MALE_SPECIFIC_SITES: frozenset[str] = frozenset({"testis", "prostate"})
FEMALE_SPECIFIC_SITES: frozenset[str] = frozenset(
    {
        "ovary",
        "uterus",
        "vagina",
        "cervix_ectocervix",
        "cervix_endocervix",
        "fallopian_tube",
    }
)

# Histological / developmental lineage families: shared markers within a family reflect
# shared cell types (e.g. stratified squamous mucosa, GI enterocytes, sarcomeric muscle,
# CNS/neuroendocrine parenchyma, or mammary adipose pad) rather than foreign lab splash.
ORGAN_LINEAGE_FAMILIES: dict[str, str] = {
    "esophagus_mucosa": "squamous_mucosa",
    "vagina": "squamous_mucosa",
    "cervix_ectocervix": "squamous_mucosa",
    "cervix_endocervix": "squamous_mucosa",
    "minor_salivary_gland": "squamous_mucosa",
    "skin_sun_exposed_lower_leg": "squamous_mucosa",
    "skin_not_sun_exposed_suprapubic": "squamous_mucosa",
    "stomach": "gi_mucosa",
    "small_intestine_terminal_ileum": "gi_mucosa",
    "colon_transverse": "gi_mucosa",
    "colon_sigmoid": "gi_mucosa",
    "esophagus_gastroesophageal_junction": "gi_mucosa",
    "esophagus_muscularis": "gi_mucosa",
    "muscle_skeletal": "striated_cardiac_muscle",
    "heart_atrial_appendage": "striated_cardiac_muscle",
    "heart_left_ventricle": "striated_cardiac_muscle",
    "adipose_subcutaneous": "adipose_mammary",
    "adipose_visceral_omentum": "adipose_mammary",
    "breast_mammary_tissue": "adipose_mammary",
    "artery_aorta": "vascular_stroma",
    "artery_coronary": "vascular_stroma",
    "artery_tibial": "vascular_stroma",
    "uterus": "vascular_stroma",
    "fallopian_tube": "vascular_stroma",
    "bladder": "vascular_stroma",
    "cells_cultured_fibroblasts": "vascular_stroma",
    "whole_blood": "lymphoid_hematopoietic",
    "spleen": "lymphoid_hematopoietic",
    "cells_ebv_transformed_lymphocytes": "lymphoid_hematopoietic",
    "kidney_cortex": "renal_parenchyma",
    "kidney_medulla": "renal_parenchyma",
    "pituitary": "cns_neuroendocrine",
    "nerve_tibial": "cns_neuroendocrine",
}

# Symmetric anatomical-dissection / surgical-incision boundary pairs. Signals between
# these pairs (e.g. prostate/seminal-vesicle luminal carryover in urinary bladder,
# epidermis/muscle in a lower-leg incision, gastrocnemius around tibial nerve,
# epicardial coronary artery / ascending aorta off myocardium & lung mediastinum,
# pectoralis under breast, or cardia mucosa at the GEJ/stomach border) are recorded as
# ``anatomical_boundary`` and NEVER subtracted from somatic safety ``N``.
_RAW_BOUNDARY_PAIRS: tuple[tuple[str, str], ...] = (
    ("prostate", "bladder"),
    ("muscle_skeletal", "thyroid"),
    ("heart_atrial_appendage", "thyroid"),
    ("heart_left_ventricle", "thyroid"),
    ("lung", "thyroid"),
    ("esophagus_mucosa", "thyroid"),
    ("esophagus_muscularis", "thyroid"),
    ("muscle_skeletal", "prostate"),
    ("heart_left_ventricle", "prostate"),
    ("muscle_skeletal", "esophagus_muscularis"),
    ("muscle_skeletal", "esophagus_gastroesophageal_junction"),
    ("muscle_skeletal", "minor_salivary_gland"),
    ("muscle_skeletal", "nerve_tibial"),
    ("muscle_skeletal", "artery_tibial"),
    ("muscle_skeletal", "breast_mammary_tissue"),
    ("muscle_skeletal", "adipose_subcutaneous"),
    ("skin_sun_exposed_lower_leg", "muscle_skeletal"),
    ("skin_not_sun_exposed_suprapubic", "muscle_skeletal"),
    ("skin_sun_exposed_lower_leg", "adipose_subcutaneous"),
    ("skin_not_sun_exposed_suprapubic", "adipose_subcutaneous"),
    ("skin_sun_exposed_lower_leg", "nerve_tibial"),
    ("skin_sun_exposed_lower_leg", "artery_tibial"),
    ("esophagus_mucosa", "stomach"),
    ("esophagus_mucosa", "esophagus_gastroesophageal_junction"),
    ("esophagus_mucosa", "esophagus_muscularis"),
    ("adrenal_gland", "kidney_cortex"),
    ("adrenal_gland", "kidney_medulla"),
    ("adrenal_gland", "adipose_visceral_omentum"),
    ("pancreas", "stomach"),
    ("pancreas", "small_intestine_terminal_ileum"),
    ("pancreas", "spleen"),
    ("adipose_visceral_omentum", "pancreas"),
    ("adipose_visceral_omentum", "stomach"),
    ("adipose_visceral_omentum", "colon_transverse"),
    ("adipose_visceral_omentum", "colon_sigmoid"),
    ("adipose_visceral_omentum", "small_intestine_terminal_ileum"),
    ("adipose_visceral_omentum", "spleen"),
    ("heart_atrial_appendage", "artery_coronary"),
    ("heart_left_ventricle", "artery_coronary"),
    ("heart_atrial_appendage", "artery_aorta"),
    ("heart_left_ventricle", "artery_aorta"),
    ("lung", "artery_aorta"),
    ("lung", "artery_coronary"),
    ("lung", "heart_atrial_appendage"),
    ("lung", "heart_left_ventricle"),
    ("lung", "esophagus_muscularis"),
    ("liver", "stomach"),
    ("liver", "adrenal_gland"),
    ("liver", "kidney_cortex"),
    ("liver", "spleen"),
)
ANATOMICAL_BOUNDARY_PAIRS: frozenset[tuple[str, str]] = frozenset(
    pair for a, b in _RAW_BOUNDARY_PAIRS for pair in ((a, b), (b, a))
)


@dataclass(frozen=True)
class GtexSampleBarcode:
    """Parsed donor and Broad 96-well ``SM-XXXX`` prep-plate provenance from a GTEx ID."""

    sample_id: str
    donor_id: str
    site_code: str
    sm_plate_prefix: str
    sm_well_id: str


def organ_lineage_family(site: str) -> str:
    """Return the histological lineage family for a normalized GTEx shard site name."""
    norm = site.strip().lower().removeprefix("v11-")
    if norm.startswith("brain_"):
        return "cns_neuroendocrine"
    return ORGAN_LINEAGE_FAMILIES.get(norm, norm)


def is_anatomical_boundary(source_site: str, target_site: str) -> bool:
    """Return True if ``source_site -> target_site`` reflects shared histology, boundary, or blood."""
    s = source_site.strip().lower().removeprefix("v11-")
    t = target_site.strip().lower().removeprefix("v11-")
    if s == t or organ_lineage_family(s) == organ_lineage_family(t):
        return True
    # Unperfused post-mortem biopsies naturally contain intravascular blood.
    if s == "whole_blood" and not t.startswith("cells_"):
        return True
    return (s, t) in ANATOMICAL_BOUNDARY_PAIRS


def parse_gtex_sample_barcode(sample_id: str) -> GtexSampleBarcode:
    """Extract donor ID (``GTEX-XXXXX``) and 4-char ``SM-XXXX`` plate prefix from a GTEx ID."""
    raw = str(sample_id).strip()
    parts = raw.split("-")
    if len(parts) >= 2 and parts[0] == "GTEX":
        donor = f"{parts[0]}-{parts[1]}"
        site_code = parts[2] if len(parts) >= 3 and parts[2] != "SM" else ""
    else:
        donor = parts[0] if parts else raw
        site_code = ""
    sm_plate = ""
    sm_well = ""
    if "SM" in parts:
        idx = parts.index("SM")
        if idx + 1 < len(parts):
            sm_well = parts[idx + 1]
            sm_plate = sm_well[:4]
    return GtexSampleBarcode(
        sample_id=raw,
        donor_id=donor,
        site_code=site_code,
        sm_plate_prefix=sm_plate,
        sm_well_id=sm_well,
    )


def infer_chromosomal_sex(
    tpm_matrix: pd.DataFrame,
    *,
    y_genes: tuple[str, ...] = Y_CHROMOSOME_SENTINEL_GENES,
    xy_sum_threshold: float = 20.0,
    xx_sum_threshold: float = 2.0,
) -> tuple[pd.Series, pd.Series]:
    """Infer sample-level and splash-robust donor-level chromosomal sex (``XX`` vs ``XY``).

    Matches both unversioned (``ENSG00000129824``) and GENCODE-versioned
    (``ENSG00000129824.16``) Ensembl IDs, leaving all-``NaN`` samples/donors ``"unassessed"``.
    """
    sample_ids = [str(c) for c in tpm_matrix.columns]
    y_set = {str(g).split(".")[0] for g in y_genes}
    unver_idx = pd.Index([str(idx).split(".")[0] for idx in tpm_matrix.index])
    match_mask = unver_idx.isin(y_set)
    if not bool(match_mask.any()) or not sample_ids:
        empty_s = pd.Series("unassessed", index=sample_ids, dtype=object)
        donors = sorted({parse_gtex_sample_barcode(s).donor_id for s in sample_ids})
        return empty_s, pd.Series("unassessed", index=donors, dtype=object)

    sub = tpm_matrix.loc[match_mask].apply(pd.to_numeric, errors="coerce")
    has_obs = sub.notna().any(axis=0).to_numpy()
    y_sum_vals = sub.sum(axis=0, min_count=1).to_numpy(dtype=np.float64)
    y_sum = pd.Series(y_sum_vals, index=sample_ids)

    sample_sex = pd.Series(
        np.where(
            ~has_obs,
            "unassessed",
            np.where(
                y_sum_vals >= xy_sum_threshold,
                "XY",
                np.where(y_sum_vals <= xx_sum_threshold, "XX", "ambiguous"),
            ),
        ),
        index=sample_ids,
        dtype=object,
    )
    donor_ids = [parse_gtex_sample_barcode(s).donor_id for s in sample_ids]
    donor_y_med = y_sum.groupby(donor_ids, sort=True).median()
    d_vals = donor_y_med.to_numpy(dtype=np.float64)
    d_fin = np.isfinite(d_vals)
    donor_sex = pd.Series(
        np.where(
            ~d_fin,
            "unassessed",
            np.where(
                d_vals >= xy_sum_threshold,
                "XY",
                np.where(d_vals <= xx_sum_threshold, "XX", "ambiguous"),
            ),
        ),
        index=donor_y_med.index,
        dtype=object,
    )
    return sample_sex, donor_sex


def select_lineage_disjoint_anchors(
    shard_medians: pd.DataFrame,
    shard_p95s: pd.DataFrame,
    *,
    shard_sample_counts: dict[str, int] | None = None,
    symbols: dict[str, str] | None = None,
    min_source_median_tpm: float = 75.0,
    min_episodic_source_p95_tpm: float = 250.0,
    max_off_lineage_p95_tpm: float = 5.0,
    max_off_lineage_ratio: float = 0.01,
    max_off_lineage_shard_fraction: float = 0.12,
    hard_max_off_lineage_ratio: float = 0.10,
    target_exclusion_ratio: float = 0.0025,
    min_p95_shard_samples: int = 20,
    max_anchors_per_source: int = 10,
    min_anchors_per_source: int = 3,
) -> pd.DataFrame:
    """Select lineage-disjoint anchor genes across GTEx shards using medians and biopsy P95s.

    Requires ``p95 <= max(max_off_lineage_p95_tpm, max_off_lineage_ratio * ref_scale)`` across
    at least ``1 - max_off_lineage_shard_fraction`` (default ``>=88%``) of well-powered
    off-lineage shards while rejecting any marker with ``>10%`` of source expression in any
    off-lineage shard. Crucially, **every** off-lineage shard whose endogenous ``p95`` exceeds
    ``max(0.5, target_exclusion_ratio * src_med)`` (below the ``0.5%`` splash threshold) is
    recorded in ``excluded_target_sites`` so moderate-abundance anchors (e.g. ``CACNG7`` in
    testis/cerebellum or ``TNXA`` readthrough in adrenal stroma) never trigger false-positive
    splash calls in tissues with low-level endogenous expression.
    """
    common_genes = shard_medians.index.intersection(shard_p95s.index)
    common_sites = [c for c in shard_medians.columns if c in shard_p95s.columns]
    med = shard_medians.loc[common_genes, common_sites].astype(float)
    p95 = shard_p95s.loc[common_genes, common_sites].astype(float)
    sym_map = symbols or {}

    cols = [
        "gene_id",
        "symbol",
        "source_site",
        "source_family",
        "source_median_tpm",
        "source_p95_tpm",
        "max_off_lineage_p95_tpm",
        "excluded_target_sites",
    ]
    rows: list[dict] = []
    for gid in common_genes:
        m_vec = med.loc[gid]
        src_site = str(m_vec.idxmax())
        src_med = float(m_vec[src_site])
        src_p95 = float(p95.at[gid, src_site])
        if src_med < min_source_median_tpm and not (
            src_med >= 25.0 and src_p95 >= min_episodic_source_p95_tpm
        ):
            continue
        src_fam = organ_lineage_family(src_site)
        off_sites = [s for s in common_sites if not is_anatomical_boundary(src_site, s)]
        if not off_sites:
            continue

        p95_sites = (
            [
                s
                for s in off_sites
                if shard_sample_counts.get(s, min_p95_shard_samples) >= min_p95_shard_samples
            ]
            if shard_sample_counts is not None
            else off_sites
        )
        eval_sites = p95_sites if p95_sites else off_sites
        ref_scale = max(src_med, 0.25 * src_p95)
        p95_ceil = max(max_off_lineage_p95_tpm, max_off_lineage_ratio * ref_scale)
        hard_ceil = max(25.0, hard_max_off_lineage_ratio * ref_scale)

        off_p95_vec = p95.loc[gid, eval_sites]
        off_p95_max = float(off_p95_vec.max())
        if off_p95_max > hard_ceil:
            continue

        broad_exceed = [str(s) for s in eval_sites if float(off_p95_vec[s]) > p95_ceil]
        if len(broad_exceed) > math.floor(max_off_lineage_shard_fraction * len(eval_sites)):
            continue

        clean_broad = [s for s in off_sites if s not in broad_exceed]
        if clean_broad:
            off_med_ceil = max(2.0, 0.005 * ref_scale)
            if float(m_vec[clean_broad].max()) > off_med_ceil:
                continue

        # Exclude any target shard where endogenous P95 or median could mimic a >= 0.5% splash!
        per_site_p95_ceil = max(0.5, target_exclusion_ratio * max(src_med, 1.0))
        per_site_med_ceil = max(0.20, 0.001 * max(src_med, 1.0))
        excluded_targets = sorted(
            {
                str(s)
                for s in off_sites
                if float(p95.at[gid, s]) > per_site_p95_ceil or float(m_vec[s]) > per_site_med_ceil
            }
        )
        # Still require the anchor to remain usable across >= 70% of off-lineage shards
        if len(excluded_targets) > math.floor(0.30 * len(off_sites)):
            continue

        kept_eval = [s for s in eval_sites if s not in excluded_targets]
        clean_p95_max = float(p95.loc[gid, kept_eval].max()) if kept_eval else off_p95_max
        rows.append(
            {
                "gene_id": str(gid),
                "symbol": str(sym_map.get(str(gid), str(gid))),
                "source_site": src_site,
                "source_family": src_fam,
                "source_median_tpm": src_med,
                "source_p95_tpm": src_p95,
                "max_off_lineage_p95_tpm": clean_p95_max,
                "excluded_target_sites": ";".join(excluded_targets),
            }
        )
    if not rows:
        return pd.DataFrame(columns=cols)
    cand = pd.DataFrame(rows, columns=cols).sort_values(
        ["source_site", "source_median_tpm", "source_p95_tpm", "gene_id"],
        ascending=[True, False, False, True],
    )
    capped = cand.groupby("source_site", group_keys=False).head(max_anchors_per_source)
    counts = capped["source_site"].value_counts()
    valid_sources = set(counts[counts >= min_anchors_per_source].index)
    return capped[capped["source_site"].isin(valid_sources)].reset_index(drop=True)


def solve_nonnegative_archetype_mixture(
    basis_tpm: pd.DataFrame,
    sample_tpm: pd.DataFrame,
    *,
    l1_penalty: float = 1e-3,
    max_iter: int = 80,
    tol: float = 1e-6,
) -> pd.DataFrame:
    """Solve joint non-negative least squares ``Y_tilde ~= A_tilde @ W`` (``W >= 0``).

    ``NaN``/``inf`` safe: sanitizes missing observations to ``0.0`` via ``np.nan_to_num`` +
    ``np.fmax`` so a missing cell never poisons the Gram matrix or right-hand side.
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


def select_nnls_basis_genes(
    shard_medians: pd.DataFrame,
    anchors: pd.DataFrame | None = None,
    *,
    top_per_shard: int = 8,
    min_peak_tpm: float = 15.0,
    min_fold_over_off_lineage_p90: float = 3.0,
) -> list[str]:
    """Select a balanced multi-tissue marker gene set covering all 52 GTEx shards for NNLS.

    Scores each shard's candidate markers by ``med[site] / max(off_lineage_p90, 1.0)`` so
    pan-epithelial genes highly expressed in unrelated organs (e.g. ``lung``) are not selected
    as basis markers for smaller shards (e.g. ``ovary``).
    """
    chosen: set[str] = set()
    if anchors is not None and not anchors.empty:
        chosen.update(str(g) for g in anchors["gene_id"].unique() if g in shard_medians.index)

    med = shard_medians.astype(float).fillna(0.0)
    all_cols = list(med.columns)
    for site in all_cols:
        sfam = organ_lineage_family(str(site))
        off_cols = [c for c in all_cols if organ_lineage_family(str(c)) != sfam]
        if off_cols:
            off_p90 = np.percentile(med[off_cols].to_numpy(dtype=np.float64), 90.0, axis=1)
            denom = np.maximum(off_p90, 1.0)
        else:
            denom = np.maximum(med.median(axis=1).to_numpy(dtype=np.float64), 1.0)
        spec_series = pd.Series(med[site].to_numpy(dtype=np.float64) / denom, index=med.index)
        eligible = med.index[
            med[site].ge(min_peak_tpm) & spec_series.ge(min_fold_over_off_lineage_p90)
        ]
        if len(eligible) == 0:
            eligible = med.index[med[site].ge(min_peak_tpm)]
        top_genes = spec_series.loc[eligible].nlargest(top_per_shard).index
        chosen.update(str(g) for g in top_genes)
    return sorted(chosen)


def deconvolve_gtex_samples(
    anchor_tpm: pd.DataFrame,
    sample_sites: pd.Series,
    anchors: pd.DataFrame,
    *,
    basis_medians: pd.DataFrame | None = None,
    basis_sample_tpm: pd.DataFrame | None = None,
    donor_sex: pd.Series | dict[str, str] | None = None,
    sample_sex: pd.Series | dict[str, str] | None = None,
    cross_donor_batch_pairs: set[tuple[str, str]] | None = None,
    min_sentinels: int = 3,
    min_marker_tpm: float = 0.25,
    min_hop_fraction: float = 2e-5,
    splash_threshold: float = 5e-3,
    swap_threshold: float = 0.20,
    gross_chimera_threshold: float = 0.50,
    max_log10_iqr: float = 0.60,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Deconvolve foreign organ contamination across GTEx libraries (#778)."""
    if min_sentinels < 3:
        raise ValueError("min_sentinels must be >= 3 to prevent single-gene masking")
    sample_ids = [str(s) for s in anchor_tpm.columns]
    sites_norm = pd.Series(
        [str(sample_sites[s]).strip().lower().removeprefix("v11-") for s in sample_ids],
        index=sample_ids,
    )
    barcodes = {s: parse_gtex_sample_barcode(s) for s in sample_ids}

    donor_sex_map: dict[str, str] = {}
    if donor_sex is not None:
        donor_sex_map = {str(k): str(v) for k, v in dict(donor_sex).items()}
    sample_sex_map: dict[str, str] = {}
    if sample_sex is not None:
        sample_sex_map = {str(k): str(v) for k, v in dict(sample_sex).items()}
    batch_pairs = cross_donor_batch_pairs or set()

    # Key SM-XXXX plate donor occupancy by EXACT shard site (not 15-shard lineage family!)
    # so e.g. nerve_tibial on a plate never falsely flags pituitary as co-plated.
    plate_site_to_donors: dict[tuple[str, str], set[str]] = {}
    donor_to_sites: dict[str, set[str]] = {}
    plate_to_samples: dict[str, list[str]] = {}
    for s in sample_ids:
        bc = barcodes[s]
        st = sites_norm[s]
        donor_to_sites.setdefault(bc.donor_id, set()).add(st)
        if bc.sm_plate_prefix:
            plate_site_to_donors.setdefault((bc.sm_plate_prefix, st), set()).add(bc.donor_id)
            plate_to_samples.setdefault(bc.sm_plate_prefix, []).append(s)

    anchors_in_mat = anchors[anchors["gene_id"].isin(anchor_tpm.index)].copy()
    if "excluded_target_sites" not in anchors_in_mat.columns:
        anchors_in_mat["excluded_target_sites"] = ""
    anchors_in_mat["excluded_target_sites"] = anchors_in_mat["excluded_target_sites"].fillna("")
    mat = anchor_tpm.reindex(anchors_in_mat["gene_id"].unique()).fillna(0.0).astype(float)

    if basis_medians is not None and basis_sample_tpm is not None:
        nnls_basis = basis_medians.copy()
        nnls_basis.columns = [
            str(c).strip().lower().removeprefix("v11-") for c in nnls_basis.columns
        ]
        nnls_obs = basis_sample_tpm.reindex(columns=sample_ids).fillna(0.0)
    else:
        emp_cols: dict[str, pd.Series] = {}
        for st, s_idxs in sites_norm.groupby(sites_norm, sort=True).groups.items():
            emp_cols[str(st)] = mat[list(s_idxs)].median(axis=1)
        for src, grp in anchors_in_mat.groupby("source_site", sort=True):
            if src not in emp_cols:
                vec = pd.Series(0.0, index=mat.index)
                vec.loc[grp["gene_id"].to_numpy()] = grp["source_median_tpm"].to_numpy(dtype=float)
                emp_cols[str(src)] = vec
        nnls_basis = pd.DataFrame(emp_cols, index=mat.index)
        nnls_obs = mat

    nnls_basis_with_pool = nnls_basis.copy()
    nnls_basis_with_pool["__multiplex_pool__"] = nnls_obs.mean(axis=1)
    nnls_weights = solve_nonnegative_archetype_mixture(
        nnls_basis_with_pool,
        nnls_obs,
        l1_penalty=1e-3,
    )
    basis_arch_sites = [c for c in nnls_basis.columns if c != "__multiplex_pool__"]
    arch_has_signal = {c: bool(float(nnls_basis[c].max()) >= 10.0) for c in basis_arch_sites}
    fam_to_arch_cols: dict[str, list[str]] = {}
    for c in basis_arch_sites:
        if arch_has_signal[c]:
            fam_to_arch_cols.setdefault(organ_lineage_family(c), []).append(c)

    native_retention = pd.Series(1.0, index=sample_ids, dtype=float)
    has_native_basis: dict[str, bool] = {}
    for site, s_idxs in sites_norm.groupby(sites_norm, sort=True).groups.items():
        fam_cols = fam_to_arch_cols.get(organ_lineage_family(str(site)), [])
        has_native_basis[str(site)] = len(fam_cols) > 0
        if fam_cols:
            native_retention.loc[list(s_idxs)] = nnls_weights.loc[fam_cols, list(s_idxs)].sum(
                axis=0
            )

    anchors_by_src = {
        src: grp.set_index("gene_id")
        for src, grp in anchors_in_mat.groupby("source_site", sort=True)
    }

    global_pool = mat.mean(axis=1)
    plate_hop_beta = pd.Series(0.0, index=sample_ids, dtype=float)
    fams_by_anchor = anchors_in_mat.drop_duplicates("gene_id").set_index("gene_id")["source_family"]
    for s in sample_ids:
        t_site = sites_norm[s]
        t_fam = organ_lineage_family(t_site)
        bc = barcodes[s]
        peer_samps = [p for p in plate_to_samples.get(bc.sm_plate_prefix, ()) if p != s]
        pool_vec = mat[peer_samps].mean(axis=1) if len(peer_samps) >= 4 else global_pool
        off_genes = [
            g
            for g in mat.index
            if fams_by_anchor[g] != t_fam and pool_vec[g] >= 50.0 and mat.at[g, s] >= 0.05
        ]
        if len(off_genes) >= min_sentinels:
            ratios_pool = (mat.loc[off_genes, s] / pool_vec.loc[off_genes]).to_numpy()
            med_b = float(np.median(ratios_pool))
            if 1e-5 <= med_b <= 5e-3:
                plate_hop_beta[s] = med_b

    event_rows: list[dict] = []
    for target_site, target_samps_idx in sites_norm.groupby(sites_norm, sort=True).groups.items():
        target_samps = list(target_samps_idx)
        target_fam = organ_lineage_family(target_site)
        sub_target = mat[target_samps]
        for src_site, src_grp_full in anchors_by_src.items():
            if src_site == target_site or organ_lineage_family(src_site) == target_fam:
                continue
            keep_mask = [
                target_site not in {x for x in str(ex).split(";") if x}
                for ex in src_grp_full["excluded_target_sites"]
            ]
            src_grp = src_grp_full.loc[keep_mask]
            if len(src_grp) < min_sentinels:
                continue

            vals = sub_target.loc[src_grp.index].to_numpy()  # M x N_t
            meds = np.maximum(src_grp["source_median_tpm"].to_numpy()[:, None], 1.0)
            syms = src_grp["symbol"].to_numpy()
            active = vals >= min_marker_tpm
            n_det = active.sum(axis=0)
            cand_js = np.where(n_det >= min_sentinels)[0]
            boundary = is_anatomical_boundary(src_site, target_site)

            for j in cand_js:
                sid = target_samps[j]
                det_mask = active[:, j]
                r_vec = vals[det_mask, j] / meds[det_mask, 0]
                stoich_med = float(np.median(r_vec))
                if stoich_med < min_hop_fraction:
                    continue
                log_r = np.log10(np.maximum(r_vec, 1e-12))
                q75, q25 = np.percentile(log_r, [75.0, 25.0])
                iqr = float(q75 - q25)
                if iqr > max_log10_iqr:
                    continue

                # Gate NNLS shadow-source suppression at splash_threshold (0.005), NOT 0.01!
                nnls_w = (
                    float(nnls_weights.at[src_site, sid])
                    if src_site in nnls_weights.index
                    else stoich_med
                )
                if stoich_med >= splash_threshold and src_site in nnls_weights.index:
                    if nnls_w < 0.15 * stoich_med:
                        continue
                    alpha_hat = float(math.sqrt(stoich_med * nnls_w))
                else:
                    alpha_hat = stoich_med
                if alpha_hat < min_hop_fraction:
                    continue

                alpha_lo = float(min(alpha_hat, np.percentile(r_vec, 25.0)))
                bc = barcodes[sid]
                plate_donors = (
                    plate_site_to_donors.get((bc.sm_plate_prefix, src_site), set())
                    if bc.sm_plate_prefix
                    else set()
                )
                co_plated_other_donor = bool(any(d != bc.donor_id for d in plate_donors))
                co_plated_same_donor = bc.donor_id in plate_donors
                same_donor_sequenced = src_site in donor_to_sites.get(bc.donor_id, ())

                d_sex = donor_sex_map.get(bc.donor_id, "unassessed")
                s_sex = sample_sex_map.get(sid, "unassessed")
                cross_sex_impossible = (d_sex == "XX" and src_site in MALE_SPECIFIC_SITES) or (
                    d_sex == "XY" and src_site in FEMALE_SPECIFIC_SITES
                )
                sex_discordant_sample = {d_sex, s_sex} == {"XX", "XY"}
                shared_batch_other_donor = (sid, src_site) in batch_pairs

                nat_ret = float(native_retention[sid])
                has_nat = has_native_basis.get(target_site, False)

                if not boundary and (
                    (has_nat and alpha_hat >= swap_threshold and nat_ret < 0.25)
                    or alpha_hat >= gross_chimera_threshold
                    or (sex_discordant_sample and alpha_hat >= swap_threshold)
                ):
                    regime = "tube_swap"
                elif boundary and not cross_sex_impossible:
                    regime = "anatomical_boundary"
                elif alpha_hat >= splash_threshold and (
                    cross_sex_impossible
                    or shared_batch_other_donor
                    or (co_plated_other_donor and not co_plated_same_donor)
                ):
                    regime = "plate_pipetting_splash"
                elif alpha_hat >= splash_threshold:
                    regime = "same_donor_carryover"
                else:
                    regime = "flowcell_index_hopping"

                top_syms = ";".join(sorted(str(x) for x in syms[det_mask]))
                event_rows.append(
                    {
                        "sample_id": sid,
                        "donor_id": bc.donor_id,
                        "donor_sex": d_sex,
                        "sample_sex": s_sex,
                        "sm_plate_prefix": bc.sm_plate_prefix,
                        "target_site": target_site,
                        "source_site": src_site,
                        "regime": regime,
                        "alpha_hat": alpha_hat,
                        "nnls_weight": nnls_w,
                        "alpha_lower_q25": alpha_lo,
                        "log10_iqr": iqr,
                        "anchors_detected": int(det_mask.sum()),
                        "anchor_symbols": top_syms,
                        "co_plated_on_sm_batch": co_plated_other_donor,
                        "cross_donor_batch_matched": shared_batch_other_donor,
                        "cross_sex_impossible": cross_sex_impossible,
                        "same_donor_organ_present": same_donor_sequenced
                        and not cross_sex_impossible,
                        "native_retention": nat_ret,
                    }
                )

    events_cols = [
        "sample_id",
        "donor_id",
        "donor_sex",
        "sample_sex",
        "sm_plate_prefix",
        "target_site",
        "source_site",
        "regime",
        "alpha_hat",
        "nnls_weight",
        "alpha_lower_q25",
        "log10_iqr",
        "anchors_detected",
        "anchor_symbols",
        "co_plated_on_sm_batch",
        "cross_donor_batch_matched",
        "cross_sex_impossible",
        "same_donor_organ_present",
        "native_retention",
    ]
    events_df = pd.DataFrame(event_rows, columns=events_cols)
    tech_regimes = {
        "tube_swap",
        "plate_pipetting_splash",
        "same_donor_carryover",
        "flowcell_index_hopping",
    }
    tech_alpha_by_sample: dict[str, float] = {}
    swap_samples: set[str] = set()
    if not events_df.empty:
        for sid, grp in events_df.groupby("sample_id"):
            if grp["regime"].eq("tube_swap").any():
                swap_samples.add(str(sid))
            tech_sub = grp[grp["regime"].isin(tech_regimes)]
            if not tech_sub.empty:
                fam_max = float(
                    tech_sub.assign(sfam=tech_sub["source_site"].map(organ_lineage_family))
                    .groupby("sfam")["alpha_hat"]
                    .max()
                    .sum()
                )
                tech_alpha_by_sample[str(sid)] = float(min(1.0, fam_max))

    summary_rows = [
        {
            "sample_id": s,
            "donor_id": barcodes[s].donor_id,
            "donor_sex": donor_sex_map.get(barcodes[s].donor_id, "unassessed"),
            "sample_sex": sample_sex_map.get(s, "unassessed"),
            "sm_plate_prefix": barcodes[s].sm_plate_prefix,
            "site": sites_norm[s],
            "native_retention": float(native_retention[s]),
            "plate_hop_beta": float(plate_hop_beta[s]),
            "technical_bleed_total": float(tech_alpha_by_sample.get(s, 0.0)),
            "purity_weight": (
                0.0
                if s in swap_samples
                else float(max(0.0, 1.0 - min(1.0, tech_alpha_by_sample.get(s, 0.0))))
            ),
            "quarantined_swap": s in swap_samples,
        }
        for s in sample_ids
    ]
    return events_df, pd.DataFrame(summary_rows)


def _betacf_scalar(a: float, b: float, x: float) -> float:
    """Modified Lentz continued fraction for the incomplete beta function (Numerical Recipes)."""
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


def _regularized_betainc(a: float, b: float, x_arr: np.ndarray) -> np.ndarray:
    """Exact scipy-free regularized incomplete beta ``I_x(a, b)`` over ``x_arr in [0, 1]``."""
    out = np.empty_like(x_arr, dtype=np.float64)
    log_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    sym_pivot = (a + 1.0) / (a + b + 2.0)
    for idx, x_raw in enumerate(x_arr):
        x = float(x_raw)
        if x <= 0.0:
            out[idx] = 0.0
        elif x >= 1.0:
            out[idx] = 1.0
        else:
            log_bt = a * math.log(x) + b * math.log1p(-x) - log_beta
            bt = math.exp(log_bt) if log_bt > -740.0 else 0.0
            if bt == 0.0:
                out[idx] = 0.0 if x < sym_pivot else 1.0
            elif x < sym_pivot:
                out[idx] = bt * _betacf_scalar(a, b, x) / a
            else:
                out[idx] = 1.0 - bt * _betacf_scalar(b, a, 1.0 - x) / b
    return np.clip(out, 0.0, 1.0)


def purity_weighted_harrell_davis_quantile(
    values: np.ndarray | pd.Series,
    weights: np.ndarray | pd.Series | None = None,
    quantile: float = 0.95,
) -> float:
    """Compute exact purity-weighted Harrell-Davis Beta quantile without ``scipy``.

    Evaluates ``w_i = I_{U_i}(a, b) - I_{U_{i-1}}(a, b)`` with ``a = q*(n_eff + 1)`` and
    ``b = (1 - q)*(n_eff + 1)`` via the exact Lentz continued-fraction regularized
    incomplete beta function across all ``n_eff`` and ``q`` (no jump discontinuity and exact
    endpoint behavior at ``q = 0.99``).
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
    cdf = _regularized_betainc(a, b, edges)
    hd_weights = np.fmax(0.0, np.diff(cdf))
    hd_sum = float(hd_weights.sum())
    if hd_sum <= 0.0 or not math.isfinite(hd_sum):
        cum = np.cumsum(w_norm) - 0.5 * w_norm
        return float(np.interp(quantile, cum, x))
    hd_weights /= hd_sum
    return float(np.dot(hd_weights, x))


def summarize_deconvolved_donor_tpm(
    values: pd.Series,
    *,
    site: str,
    events: pd.DataFrame | None = None,
    source_medians: dict[str, float] | None = None,
    source_p95s: dict[str, float] | None = None,
    orthogonal_supported: bool = False,
    mode: str = "normal_p95",
    quantile: float = 0.95,
    heavy_splash_alpha: float = 0.005,
) -> dict[str, float | int]:
    """Summarize gene TPM across GTEx biopsies after swap quarantine & bleed deconvolution."""
    if mode not in {"normal_p95", "host_exposure"}:
        raise ValueError(f"Unknown mode {mode!r}; expected 'normal_p95' or 'host_exposure'")
    s = pd.to_numeric(values, errors="coerce").dropna()
    if s.empty:
        return {
            "value": float("nan"),
            "raw_p95": float("nan"),
            "raw_max": float("nan"),
            "clean_donors": 0,
            "dropped_samples": 0,
        }

    barcodes = {str(sid): parse_gtex_sample_barcode(str(sid)) for sid in s.index}
    raw_donors = s.groupby([barcodes[str(sid)].donor_id for sid in s.index], sort=False).median()
    raw_p95 = float(np.percentile(raw_donors.to_numpy(dtype=float), 95.0))
    raw_max = float(s.max())

    src_med = source_medians or {}
    src_p95 = source_p95s or src_med
    norm_site = site.strip().lower().removeprefix("v11-")

    ev_by_samp: dict[str, list[dict]] = {}
    if events is not None and not events.empty:
        sub_ev = events[events["target_site"].eq(norm_site)]
        for rec in sub_ev.to_dict("records"):
            ev_by_samp.setdefault(str(rec["sample_id"]), []).append(rec)

    adj_vals: dict[str, float] = {}
    samp_weights: dict[str, float] = {}
    dropped = 0
    tech_regimes = {"plate_pipetting_splash", "same_donor_carryover", "flowcell_index_hopping"}

    for sid_raw, y_raw in s.items():
        sid = str(sid_raw)
        y = float(y_raw)
        evs = ev_by_samp.get(sid, ())
        if any(e["regime"] == "tube_swap" for e in evs):
            dropped += 1
            continue
        if orthogonal_supported or not evs:
            adj_vals[sid] = y
            samp_weights[sid] = 1.0
            continue

        tech_evs = [e for e in evs if e["regime"] in tech_regimes]
        if not tech_evs:
            adj_vals[sid] = y
            samp_weights[sid] = 1.0
            continue

        heavy_hit = False
        expected_med_bleed = 0.0
        expected_lo_bleed = 0.0
        total_alpha = 0.0
        by_fam: dict[str, dict] = {}
        for e in tech_evs:
            sfam = organ_lineage_family(str(e["source_site"]))
            if sfam not in by_fam or float(e["alpha_hat"]) > float(by_fam[sfam]["alpha_hat"]):
                by_fam[sfam] = e

        for e in by_fam.values():
            ssrc = str(e["source_site"])
            a_hat = float(e["alpha_hat"])
            a_lo = float(e["alpha_lower_q25"])
            total_alpha += a_hat
            g_p95 = float(src_p95.get(ssrc, 0.0))
            g_med = float(src_med.get(ssrc, 0.0))
            bleed_p95 = a_hat * max(g_p95, g_med)
            # Additive contamination can only inflate y when y > 0 AND expected bleed is both
            # non-trivial (>= 0.15 TPM) and >= 20% of observed y (never drops zeros or ACTB/SOX9!).
            if (
                y > 0.0
                and bleed_p95 >= 0.20 * y
                and (
                    (a_hat >= heavy_splash_alpha and bleed_p95 >= 0.15)
                    or (mode == "host_exposure" and g_p95 >= 50.0 and bleed_p95 >= 0.15)
                )
            ):
                heavy_hit = True
            expected_med_bleed += a_hat * g_med
            expected_lo_bleed += a_lo * g_med

        if heavy_hit:
            dropped += 1
            continue

        if y > 0.0 and expected_med_bleed > 0.0 and expected_med_bleed >= 0.25 * y:
            y_clean = max(0.0, y - expected_lo_bleed)
        else:
            y_clean = y
        adj_vals[sid] = y_clean
        samp_weights[sid] = max(0.05, 1.0 - min(0.95, total_alpha))

    donor_groups: dict[str, list[tuple[float, float]]] = {}
    for sid, val in adj_vals.items():
        did = barcodes[sid].donor_id
        donor_groups.setdefault(did, []).append((val, samp_weights[sid]))

    if len(donor_groups) < min(5, len(raw_donors)):
        donor_groups.clear()
        dropped = 0
        for sid_raw, y_raw in s.items():
            sid = str(sid_raw)
            if any(e["regime"] == "tube_swap" for e in ev_by_samp.get(sid, ())):
                dropped += 1
                continue
            donor_groups.setdefault(barcodes[sid].donor_id, []).append((float(y_raw), 1.0))

    if not donor_groups:
        return {
            "value": raw_p95 if mode == "normal_p95" else raw_max,
            "raw_p95": raw_p95,
            "raw_max": raw_max,
            "clean_donors": 0,
            "dropped_samples": dropped,
        }

    d_vals = np.empty(len(donor_groups), dtype=np.float64)
    d_wts = np.empty(len(donor_groups), dtype=np.float64)
    for i, pairs in enumerate(donor_groups.values()):
        ws = np.array([p[1] for p in pairs], dtype=np.float64)
        vs = np.array([p[0] for p in pairs], dtype=np.float64)
        d_vals[i] = float(np.average(vs, weights=ws)) if mode == "normal_p95" else float(vs.max())
        d_wts[i] = float(ws.mean())

    n_donors = len(d_vals)
    if mode == "normal_p95":
        if n_donors < 5:
            est = float(d_vals.max())
        else:
            est = purity_weighted_harrell_davis_quantile(d_vals, d_wts, quantile=quantile)
    else:
        sorted_d = np.sort(d_vals)
        if n_donors >= 50:
            p99 = float(np.percentile(sorted_d, 99.0))
            second_hi = float(sorted_d[-2])
            est = min(float(sorted_d[-1]), max(p99, second_hi) if n_donors < 100 else p99)
        elif n_donors >= 30:
            est = float(sorted_d[-2])
        else:
            est = float(sorted_d[-1])

    return {
        "value": est,
        "raw_p95": raw_p95,
        "raw_max": raw_max,
        "clean_donors": n_donors,
        "dropped_samples": dropped,
    }


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


@functools.lru_cache(maxsize=1)
def _verify_deconvolution_receipt_cached() -> dict:
    data_root = files("oncoref").joinpath("data")
    receipt = json.loads(
        data_root.joinpath("gtex-deconvolution-receipt.json").read_text(encoding="utf-8")
    )
    for fname, expected_sha in receipt["files_sha256"].items():
        actual = _sha256_bytes(data_root.joinpath(fname).read_bytes())
        if actual != expected_sha:
            raise ValueError(f"{fname} does not match gtex-deconvolution-receipt.json SHA-256")
    return receipt


def gtex_deconvolution_receipt() -> dict:
    """Return SHA-256-verified provenance and summary counts for #778 GTEx deconvolution."""
    return copy.deepcopy(_verify_deconvolution_receipt_cached())


def gtex_deconvolution_anchors() -> pd.DataFrame:
    """Return the auto-selected lineage-disjoint GTEx v11 tissue anchor genes (#778)."""
    _verify_deconvolution_receipt_cached()
    return get_data("gtex-deconvolution-anchors").copy()


def gtex_deconvolution_tissue_summary() -> pd.DataFrame:
    """Return per-shard GTEx v11 swap, plate-splash, hopping, and boundary incidence (#778)."""
    _verify_deconvolution_receipt_cached()
    return get_data("gtex-deconvolution-tissue-summary").copy()


def gtex_quarantined_swaps() -> pd.DataFrame:
    """Return GTEx v11 biopsies quarantined globally as foreign-organ swaps (#778)."""
    _verify_deconvolution_receipt_cached()
    return get_data("gtex-quarantined-swaps").copy()


def gtex_donor_sex() -> pd.DataFrame:
    """Return splash-robust chromosomal sex (``XX`` vs ``XY``) across all 946 GTEx v11 donors."""
    _verify_deconvolution_receipt_cached()
    return get_data("gtex-donor-sex").copy()


def gtex_sample_splash_events(*, include_boundary: bool = True) -> pd.DataFrame:
    """Return per-sample GTEx v11 foreign splash, swap, carryover, and boundary events (#778)."""
    _verify_deconvolution_receipt_cached()
    with files("oncoref").joinpath("data", "gtex-sample-splash-events.csv.gz").open("rb") as fh:
        df = pd.read_csv(fh, compression="gzip", keep_default_na=False)
    if not include_boundary:
        df = df[df["regime"].ne("anatomical_boundary")].reset_index(drop=True)
    return df


__all__ = [
    "ANATOMICAL_BOUNDARY_PAIRS",
    "FEMALE_SPECIFIC_SITES",
    "GTEX_DECONVOLUTION_POLICY_ID",
    "MALE_SPECIFIC_SITES",
    "ORGAN_LINEAGE_FAMILIES",
    "Y_CHROMOSOME_SENTINEL_GENES",
    "GtexSampleBarcode",
    "deconvolve_gtex_samples",
    "gtex_deconvolution_anchors",
    "gtex_deconvolution_receipt",
    "gtex_deconvolution_tissue_summary",
    "gtex_donor_sex",
    "gtex_quarantined_swaps",
    "gtex_sample_splash_events",
    "infer_chromosomal_sex",
    "is_anatomical_boundary",
    "organ_lineage_family",
    "parse_gtex_sample_barcode",
    "purity_weighted_harrell_davis_quantile",
    "select_lineage_disjoint_anchors",
    "select_nnls_basis_genes",
    "solve_nonnegative_archetype_mixture",
    "summarize_deconvolved_donor_tpm",
]

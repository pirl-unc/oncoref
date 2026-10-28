"""Provisional cell-compartment policy for CTA research candidate prioritization.

Numeric limits express a review preference, not measured toxicity thresholds.
No gene identity or clinical product name participates in these rules. Source
units, anatomical labels and the prior assessment remain available for audit.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .cell_rna_policy import EXPECTED_UNITS

POLICY_ID = "oncoref.compartment-candidates.2026-10-27.1"
TIER_BASES = {"vital": 4.0, "protected": 8.0, "manageable": 16.0, "renewing_mucosa": 32.0}
IHC_RANKS = {
    "Missing": 0,
    "Missing/uncertain": 0,
    "Uncertain": 0,
    "Approved": 1,
    "Supported": 2,
    "Enhanced": 3,
}
LARGE_BOWEL_EPITHELIAL_CELLS = frozenset(
    {"distal enterocytes", "intestinal goblet cells", "paneth cells", "undifferentiated cells"}
)


def assess_contexts(contexts: pd.DataFrame) -> pd.DataFrame:
    """Apply linked ceilings to supported, normalized normal-expression contexts.

    Requires ``value`` to already contain the source-specific supported statistic;
    missing measurements are never changed to zero. ``category`` is the previous
    anatomical tier, and ``selection_tier`` records this proposal's explicit change.
    Only mapped adult human HPA tissue clusters can receive a cell-specific tier.
    Bulk, immune and corneal assays retain their existing anatomical tier. CP10K is
    rescaled numerically without claiming equivalence to length-normalized TPM.
    """
    required = {
        "source",
        "taxonomy_id",
        "stage",
        "tissue",
        "cell_type",
        "category",
        "mapping_status",
        "negative_ihc_grade",
        "value",
        "unit",
    }
    if not required <= set(contexts):
        raise ValueError(f"Missing context fields: {sorted(required - set(contexts))}")
    if not contexts.taxonomy_id.eq(9606).all() or not contexts.stage.eq("adult").all():
        raise ValueError("Require adult human contexts")
    if not contexts.unit.isin(["TPM", "nTPM", "CP10K"]).all():
        raise ValueError("Unsupported native unit")
    assays = contexts.source.where(~contexts.source.str.startswith("dice_"), "dice")
    expected = assays.map(EXPECTED_UNITS | {"hpa_rna_consensus": "nTPM"})
    if not contexts.unit.eq(expected).all():
        raise ValueError("Source and native unit do not match")
    if not contexts.negative_ihc_grade.isin(IHC_RANKS).all():
        raise ValueError("Unknown IHC grade")
    d = contexts.copy()
    values = pd.to_numeric(d.value, errors="raise")
    if values.lt(0).any() or np.isinf(values).any():
        raise ValueError("Require finite nonnegative supported RNA or missing")
    d["selection_tier"] = d.category
    d["compartment_basis"] = "anatomical_default"
    normal = d.category.isin(["vital", "protected", "manageable"]) & d.mapping_status.eq("mapped")
    eligible = normal & d.mapping_status.eq("mapped") & d.source.eq("hpa_single_cell_tissue")
    cells = d.cell_type.fillna("").str.strip().str.lower()
    tissues = d.tissue.fillna("").str.strip().str.lower()
    stroma = eligible & cells.eq("fibroblasts")
    mucosa = eligible & tissues.isin(["colon", "rectum"]) & cells.isin(LARGE_BOWEL_EPITHELIAL_CELLS)
    d.loc[stroma, ["selection_tier", "compartment_basis"]] = ["manageable", "identified_fibroblast"]
    d.loc[mucosa, ["selection_tier", "compartment_basis"]] = [
        "renewing_mucosa",
        "large_bowel_epithelial_population",
    ]
    d["stem_progenitor_caution"] = mucosa & cells.eq("undifferentiated cells")
    scale = d.unit.map({"TPM": 1.0, "nTPM": 1.0, "CP10K": 0.01})
    d["candidate_ceiling"] = (
        d.selection_tier.map(TIER_BASES) * 2.0 ** d.negative_ihc_grade.map(IHC_RANKS) * scale
    )
    d["strong_boundary"] = 2 * d.candidate_ceiling
    d["candidate_exclusion"] = normal & values.gt(d.candidate_ceiling)
    d["unsafe_evidence"] = normal & values.ge(d.strong_boundary)
    if "ihc_positive" in d:
        level = d.get("Level", pd.Series("Unquantified", index=d.index))
        d["low_normal_protein_flag"] = normal & d.ihc_positive.fillna(False) & level.eq("Low")
        positive = normal & d.ihc_positive.fillna(False) & ~level.eq("Low")
        d["candidate_exclusion"] |= positive
        d["unsafe_evidence"] |= positive
    d["candidate_exclusion"] |= ~d.mapping_status.eq("mapped")
    d["selection_policy"] = POLICY_ID
    return d


def assess_imported_protein_reviews(
    findings: pd.DataFrame, intensities: pd.DataFrame
) -> pd.DataFrame:
    """Permit explicitly Low protein as a flag; retain other accepted positives.

    Intensity must be explicitly extracted with source provenance for every reviewed
    finding. Unextracted intensity is not assumed Low. Full-text/abstract access is
    recorded independently and never serves as a reason to erase a known positive.
    Low protein cannot grant negative-IHC credit; the RNA fallback still applies.
    """
    d = findings.loc[findings.Finding.eq("Normal protein") & findings.Status.eq("Accepted")].copy()
    if not intensities["Review ID"].is_unique or set(intensities["Review ID"]) != set(
        d["Review ID"]
    ):
        raise ValueError(
            "Intensity audit must cover every accepted normal-protein finding exactly once"
        )
    if not intensities.protein_level.isin(["Low", "Medium", "High", "Detected_unquantified"]).all():
        raise ValueError("Unknown protein intensity")
    if intensities[["intensity_source", "intensity_locator"]].fillna("").eq("").any(axis=None):
        raise ValueError("Protein intensity provenance is required")
    d = d.merge(intensities, on="Review ID", validate="one_to_one")
    depth = d["Review depth"].fillna("").str.lower()
    d["primary_full_text_reviewed"] = depth.str.contains("primary full text", regex=False)
    d["low_normal_protein_flag"] = d.protein_level.eq("Low")
    d["retained_protein_exclusion"] = ~d.low_normal_protein_flag
    d["limited_review_flag"] = ~d.primary_full_text_reviewed
    d["selection_policy"] = POLICY_ID
    return d

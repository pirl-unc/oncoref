"""Provisional cell-compartment policy for CTA research candidate prioritization.

Numeric limits express a review preference, not measured toxicity thresholds.
No gene identity or clinical product name participates in these rules. Source
units, anatomical labels and the prior assessment remain available for audit.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from .cell_rna_policy import EXPECTED_UNITS, GRADES, NEUTRAL_REPRODUCTIVE_IMMUNE, POSITIVE

POLICY_ID = "oncoref.cta-selection.1.21.0"
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


def collapse_hpa_ihc(observations: pd.DataFrame, groups: pd.DataFrame) -> pd.DataFrame:
    """Preserve the strongest observed HPA level within each locus group/context.

    HPA categories are ordinal, never additive. The selected level retains only
    confidence from observations of that level. Every original level remains in
    ``observed_levels``. Negative evidence requires every expected locus and uses
    the weakest member's confidence. Missing members cannot grant a negative.
    Callers must additionally require a positive profile control for every member
    before using a group negative to increase its RNA threshold.
    """
    keys = ["Gene", "Tissue", "Cell type"]
    if observations.duplicated(keys).any():
        raise ValueError("Repeated locus/IHC context cannot be collapsed")
    members: dict[str, str] = {}
    for group in groups.itertuples():
        for gene_id in json.loads(group.member_gene_ids):
            if gene_id in members:
                raise ValueError("A locus belongs to multiple protein equivalence groups")
            members[gene_id] = group.protein_equivalence_group_id
    metadata = groups.set_index("protein_equivalence_group_id")
    data = observations[observations.Gene.isin(members)].copy()
    data["group_id"] = data.Gene.map(members)
    rows = []
    ranks = {"Low": 1, "Medium": 2, "High": 3, "Ascending": 4, "Descending": 4}
    for (group_id, tissue, cell), part in data.groupby(["group_id", "Tissue", "Cell type"]):
        complete = part.Gene.nunique() == metadata.loc[group_id, "member_count"]
        positive = part.Level.isin(POSITIVE)
        if positive.any():
            level = max(part.loc[positive, "Level"], key=lambda value: (ranks[value], value))
            confidence = part.loc[part.Level.eq(level), "Reliability"].map(GRADES).fillna(0).max()
        elif complete and part.Level.eq("Not detected").all():
            level = "Not detected"
            confidence = part.Reliability.map(GRADES).fillna(0).min()
        else:
            level, confidence = "Missing", 0
        rows.append(
            {
                "Gene": metadata.loc[group_id, "representative_gene_id"],
                "Tissue": tissue,
                "Cell type": cell,
                "Level": level,
                "Reliability": {v: k for k, v in GRADES.items()}.get(confidence, "Uncertain"),
                "observed_levels": json.dumps(sorted(part.Level.dropna().unique())),
                "complete_members": complete,
                "measured_members": part.Gene.nunique(),
                "expected_members": int(metadata.loc[group_id, "member_count"]),
            }
        )
    return pd.DataFrame(rows)


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
    categories = {
        "vital",
        "protected",
        "manageable",
        "reproductive",
        "developmental",
        "disposable",
        "unresolved",
        "permitted",
        NEUTRAL_REPRODUCTIVE_IMMUNE,
    }
    if not contexts.category.isin(categories).all():
        raise ValueError("Unknown anatomical category")
    d = contexts.copy()
    values = pd.to_numeric(d.value, errors="raise")
    if values.lt(0).any() or np.isinf(values).any():
        raise ValueError("Require finite nonnegative supported RNA or missing")
    d["value"] = values
    if (
        "ihc_positive" in d
        and (
            d.ihc_positive.fillna(False)
            & d.negative_ihc_grade.isin(["Approved", "Supported", "Enhanced"])
        ).any()
    ):
        raise ValueError("Positive IHC cannot receive negative-IHC credit")
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
    d["stem_progenitor_caution"] = mucosa & cells.eq("undifferentiated cells") & values.gt(0)
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
    A source's qualitative description remains on its own scale and is not a
    calibrated equivalent of an HPA staining level or a numerical H-score.
    """
    required = {
        "Review ID",
        "protein_level",
        "protein_scale",
        "intensity_source",
        "intensity_locator",
    }
    if not required <= set(intensities):
        raise ValueError(f"Missing protein review fields: {sorted(required - set(intensities))}")
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
    if not intensities.protein_scale.isin(
        ["author_description", "hpa_categorical_level", "not_extracted"]
    ).all():
        raise ValueError("Unknown protein reporting scale")
    if (
        intensities.protein_scale.eq("not_extracted")
        & ~intensities.protein_level.eq("Detected_unquantified")
    ).any():
        raise ValueError("Unextracted protein intensity cannot receive a Low category")
    d = d.merge(intensities, on="Review ID", validate="one_to_one")
    depth = d["Review depth"].fillna("").str.lower()
    d["primary_full_text_reviewed"] = depth.str.contains("primary full text", regex=False)
    d["low_normal_protein_flag"] = d.protein_level.eq("Low")
    d["retained_protein_exclusion"] = ~d.low_normal_protein_flag
    d["limited_review_flag"] = ~d.primary_full_text_reviewed
    d["selection_policy"] = POLICY_ID
    return d

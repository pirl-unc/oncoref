"""Experimental, explicit two-threshold assessment of normal-cell RNA.

This module has no default thresholds and does not change released CTA or
forbidden-protein membership. Cutoffs must be supplied for each assay, unit and
anatomical tier. An RNA exclusion is evidence for reference review, not proof of
protein translation, peptide presentation, or toxicity.
"""

from __future__ import annotations

import itertools
from collections.abc import Collection, Mapping
from functools import lru_cache

import numpy as np
import pandas as pd

from . import normal_tissues

NORMAL = {"vital", "protected", "manageable"}
PERMITTED = {"reproductive", "disposable"}
GRADES = {"Approved": 1, "Supported": 2, "Enhanced": 3}
POSITIVE = {"Low", "Medium", "High", "Ascending", "Descending"}
EXPECTED_UNITS = {
    "hpa_single_cell_tissue": "nTPM",
    "hpa_immune_cell": "nTPM",
    "hpa_immune_monaco": "TPM",
    "dice": "TPM",
    "cornea_collin_adult_counts": "CP10K",
}
NEUTRAL_REPRODUCTIVE_IMMUNE = "neutral_reproductive_immune"
LINKED_POLICY = "oncoref.linked-normal-rna.2026-10-27.1"
GEOMETRIC_POLICY = "oncoref.geometric-normal-rna.2026-10-27.1"
REFERENCE_REVIEW_MULTIPLIER = 2.0


def threshold_ladder(mode: str = "legacy") -> dict[str, dict[str, float]]:
    """One shared candidate ladder for bulk, clusters and sorted immune assays.

    Doubling modes are user-directed research parameterizations, not calibrated
    assay sensitivities. Missing and Uncertain share the first column. The
    historical ladder is retained for reproducible comparisons.
    """
    legacy = normal_tissues.normal_tissue_screening_policy("ihc_adaptive_cta")[
        "negative_ihc_rna_limits"
    ]
    if mode == "legacy":
        return legacy
    if mode not in {"ihc_doubling", "two_axis_doubling"}:
        raise ValueError("Unknown shared RNA ladder")
    tiers = ["vital", "protected", "manageable"]
    base = legacy["Missing"]["vital"]
    result = {}
    for j, grade in enumerate(["Missing", *GRADES]):
        result[grade] = {
            tier: (base * 2**i if mode == "two_axis_doubling" else legacy["Missing"][tier]) * 2**j
            for i, tier in enumerate(tiers)
        }
    result["Uncertain"] = result["Missing"].copy()
    return result


def apply_bulk_ladder(
    contexts: pd.DataFrame, *, mode: str = "legacy", max_ihc_grade: str = "Enhanced"
) -> pd.DataFrame:
    """Reapply a shared ladder to already matched, fully audited bulk evidence.

    This uses evidence grades and profile controls, never old numeric cutoffs
    or selection flags. Unresolved anatomy and missing observations remain holds.
    """
    d = contexts.copy()
    levels = {"Missing": 0, "Uncertain": 0, **GRADES}
    if max_ihc_grade not in levels:
        raise ValueError("Unknown IHC grade cap")
    grade = d.ihc_reliability_used.where(d.profile_positive_control, "Missing")
    if not grade.isin(levels).all():
        raise ValueError("Unknown bulk IHC evidence grade")
    rank = grade.map(levels).clip(upper=levels[max_ihc_grade])
    d["negative_ihc_grade"] = rank.map(dict(enumerate(["Missing", *GRADES])))
    ladder = threshold_ladder(mode)
    d["effective_ceiling_ntpm"] = [
        ladder[g].get(tier, np.nan)
        for g, tier in zip(d.negative_ihc_grade, d.screening_tier, strict=True)
    ]
    permitted = d.screening_tier.eq("permitted")
    d["local_exclusion"] = ~permitted & (
        d.effective_ceiling_ntpm.isna()
        | d.rna_ntpm.isna()
        | d.rna_ntpm.gt(d.effective_ceiling_ntpm)
    )
    d["reference_review_ceiling_ntpm"] = d.effective_ceiling_ntpm * REFERENCE_REVIEW_MULTIPLIER
    d["reference_review"] = ~permitted & d.rna_ntpm.ge(d.reference_review_ceiling_ntpm)
    d["threshold_ladder"] = mode
    return d


def protected_immune_contexts(contexts: pd.DataFrame) -> pd.DataFrame:
    """Use Protected for identified immune cells, retaining local neutrality.

    This is a selection-policy default, not a change to the anatomical registry
    or a claim of tolerable immune depletion. No subset has a reviewed Manageable
    exception in this policy. An organ's nonimmune cells keep their own tier.
    """
    d = contexts.copy()
    if "prior_immune_category" in d:
        raise ValueError("Immune tier policy has already been applied")
    concept = d.get("concept_id", pd.Series("", index=d.index)).fillna("")
    group = d.get("source_cell_group", pd.Series("", index=d.index)).fillna("")
    sorted_immune = d.source.isin(["hpa_immune_cell", "hpa_immune_monaco"]) | (
        d.source.str.startswith("dice_")
    )
    immune = concept.isin(["blood", "hofbauer_cell"]) | group.eq("Blood & immune cells")
    d["prior_immune_category"] = d.category
    d["immune_tier_default"] = (immune | sorted_immune) & d.category.isin(NORMAL)
    d.loc[d.immune_tier_default, "category"] = "protected"
    return d


def linked_thresholds(*, mode: str = "legacy") -> pd.DataFrame:
    """Generate every assay's limits from the existing bulk tissue/IHC ladder.

    The linkage is a provisional numerical policy, not assay equivalence. A
    per-million cutoff is divided by 100 to obtain its native CP10K cutoff; this
    only rescales counts and does not turn UMI CPM into length-adjusted TPM.
    The stronger review boundary is explicitly 2x the candidate ceiling, not a
    validated translation/toxicity threshold. Source-specific deviations require
    a separate evidence-backed policy version; none are configured here.
    """
    ladder = threshold_ladder(mode)
    rows = []
    for assay, unit in EXPECTED_UNITS.items():
        scale = 0.01 if unit == "CP10K" else 1.0
        for grade in ["Missing/uncertain", *GRADES]:
            for tier, ceiling in ladder[
                "Missing" if grade == "Missing/uncertain" else grade
            ].items():
                rows.append(
                    {
                        "assay": assay,
                        "unit": unit,
                        "category": tier,
                        "negative_ihc_grade": grade,
                        "cta_cutoff": ceiling * scale,
                        "reference_cutoff": ceiling * REFERENCE_REVIEW_MULTIPLIER * scale,
                        "nominal_per_million_ceiling": ceiling,
                        "native_units_per_nominal_unit": scale,
                        "policy": LINKED_POLICY if mode == "legacy" else GEOMETRIC_POLICY,
                    }
                )
    return pd.DataFrame(rows)


@lru_cache(maxsize=512)
def _reproductive_sampling_site(tissue: str) -> bool:
    annotation = normal_tissues.resolve_normal_tissue(tissue, species="human")
    return annotation.mapping_status == "mapped" and "reproductive" in annotation.scopes


def reproductive_immune_contexts(contexts: pd.DataFrame) -> pd.DataFrame:
    """Mark explicitly located reproductive/developmental immune cells neutral.

    This is a user-directed CTA selection convention, not evidence of uptake,
    transcription, protein absence or therapeutic tolerability. Source lineage,
    tissue and measurements are preserved. Blood or another sampled organ does
    not inherit reproductive permission from the same cell lineage elsewhere.
    Neutral observations cannot supply a positive reproductive admission route.

    Inputs are already anatomically mapped contexts. Immune identity requires
    the mapped blood/Hofbauer concept or an explicit source immune-cell group;
    neither an unknown cell label nor a reproductive donor label is sufficient.
    """
    required = {"taxonomy_id", "tissue", "stage", "category"}
    if not required <= set(contexts) or not contexts.taxonomy_id.eq(9606).all():
        raise ValueError("Human contexts with tissue, stage and category are required")
    if "prior_context_category" in contexts:
        raise ValueError("Context neutrality has already been applied")
    d = contexts.copy()
    concept = d.get("concept_id", pd.Series("", index=d.index)).fillna("")
    group = d.get("source_cell_group", pd.Series("", index=d.index)).fillna("")
    immune = concept.isin(["blood", "hofbauer_cell"]) | group.eq("Blood & immune cells")
    site = d.tissue.fillna("").map(_reproductive_sampling_site)
    developmental = d.stage.fillna("").str.lower().isin(["fetal", "foetal", "embryonic", "embryo"])
    d["prior_context_category"] = d.category
    d["neutral_reproductive_immune"] = immune & (site | developmental)
    d["neutral_basis"] = ""
    d.loc[d.neutral_reproductive_immune & site, "neutral_basis"] = (
        "explicit_reproductive_sampling_site"
    )
    d.loc[d.neutral_reproductive_immune & developmental, "neutral_basis"] = (
        "explicit_developmental_stage"
    )
    d.loc[d.neutral_reproductive_immune, "category"] = NEUTRAL_REPRODUCTIVE_IMMUNE
    return d


def match_negative_ihc(
    contexts: pd.DataFrame,
    ihc: pd.DataFrame,
    *,
    positive_controls: Collection[str],
    conflicting_genes: Collection[str] = (),
) -> pd.DataFrame:
    """Match exact adult tissue/cell labels, retaining all mismatches as missing.

    HPA tissue clusters are normal adult observations. Only these can use the
    normal-tissue IHC table here. Sorted/activated immune states and corneal
    subtypes have no reviewed equivalences and receive no IHC relaxation.
    Labels are case/space normalized only; no fuzzy or whole-organ substitution.
    A reliable positive profile control is necessary but does not establish
    assay sensitivity in the matched normal cell or pair RNA and IHC donors.
    """
    d = contexts.copy()
    h = ihc.rename(columns={"Gene": "gene_id", "Tissue": "tissue", "Cell type": "cell_type"})
    keys = ["gene_id", "tissue", "cell_type"]
    for frame in (d, h):
        for col in ["tissue", "cell_type"]:
            frame[col] = frame[col].fillna("").str.strip().str.lower()
    if h.duplicated(keys).any():
        raise ValueError("Curated IHC must contain one result per gene/tissue/cell")
    d = d.merge(h[[*keys, "Level", "Reliability"]], on=keys, how="left", validate="many_to_one")
    matched = (
        d.source.eq("hpa_single_cell_tissue")
        & d.stage.eq("adult")
        & d.state.isin(["normal", "unspecified", "unstimulated"])
        & d.Level.notna()
    )
    d["ihc_exact_match"] = matched
    d["ihc_positive"] = matched & d.Level.isin(POSITIVE) & d.Reliability.isin(GRADES)
    usable = (
        matched
        & d.Level.eq("Not detected")
        & d.Reliability.isin(GRADES)
        & d.gene_id.isin(positive_controls)
        & ~d.gene_id.isin(conflicting_genes)
    )
    d["negative_ihc_grade"] = d.Reliability.where(usable, "Missing/uncertain")
    return d


def assess(
    contexts: pd.DataFrame,
    thresholds: pd.DataFrame,
    *,
    ihc_factors: Mapping[str, float],
    inclusive_ceiling: bool = False,
) -> pd.DataFrame:
    """Apply independently specified CTA and reference-review cutoffs.

    `value` is an assay-specific supported statistic: a cluster maximum after
    the read guard, the third highest DICE value, or the second highest distinct
    corneal sample-group mean after cell/count guards. Values absent after these
    guards are NaN, never zero. No reproductive maximum can mask normal RNA.
    By default equality fails the respective cutoff. Linked bulk/cell policy
    uses inclusive_ceiling=True for candidate admission; reference review still
    starts at its upper boundary. Below-cutoff observations are not a
    claim of adequate coverage, protein absence or gene-level eligibility.
    """
    d = contexts.copy()
    if not d.taxonomy_id.eq(9606).all() or not d.stage.eq("adult").all():
        raise ValueError("This assessment requires human adult observations")
    d["assay"] = d.source.where(~d.source.str.startswith("dice_"), "dice")
    if not d.unit.eq(d.assay.map(EXPECTED_UNITS)).all():
        raise ValueError("Unsupported source or mismatched native RNA units")
    values = pd.to_numeric(d.value, errors="raise")
    if values.lt(0).any() or np.isinf(values).any():
        raise ValueError("RNA values must be nonnegative and finite, or missing")
    t = thresholds.copy()
    keys = ["assay", "unit", "category"]
    if "negative_ihc_grade" in t:
        keys.append("negative_ihc_grade")
        if any(v != 1 for v in ihc_factors.values()):
            raise ValueError("Grade-specific thresholds cannot also use IHC multipliers")
    if t.duplicated(keys).any():
        raise ValueError("Exactly one threshold pair per assay/unit/tissue tier")
    if (
        not np.isfinite(t[["cta_cutoff", "reference_cutoff"]].to_numpy()).all()
        or t.cta_cutoff.le(0).any()
        or t.reference_cutoff.le(t.cta_cutoff).any()
    ):
        raise ValueError("Require finite 0 < CTA cutoff < reference cutoff")
    if set(ihc_factors) != {"Missing/uncertain", *GRADES}:
        raise ValueError("Every IHC tier needs an explicit multiplier")
    factors = [ihc_factors[k] for k in ["Missing/uncertain", *GRADES]]
    if (
        factors[0] != 1
        or not np.isfinite(factors).all()
        or any(a > b for a, b in itertools.pairwise(factors))
    ):
        raise ValueError("IHC factors must start at 1 and increase with reliability")
    d = d.merge(t, on=keys, how="left", validate="many_to_one")
    if d.loc[d.category.isin(NORMAL), "cta_cutoff"].isna().any():
        raise ValueError("Missing threshold for a screened assay/tissue tier")
    multiplier = d.negative_ihc_grade.map(ihc_factors)
    if multiplier.isna().any():
        raise ValueError("Unknown negative IHC grade")
    d["cta_cutoff"] *= multiplier
    d["reference_cutoff"] *= multiplier
    d["assessment"] = "Unassessed"
    normal = d.category.isin(NORMAL) & d.mapping_status.eq("mapped")
    measured = normal & d.value.notna()
    d.loc[measured, "assessment"] = "Below CTA cutoff"
    exceeds = d.value.gt(d.cta_cutoff) if inclusive_ceiling else d.value.ge(d.cta_cutoff)
    d.loc[measured & exceeds, "assessment"] = "CTA exclusion"
    d.loc[measured & d.value.ge(d.reference_cutoff), "assessment"] = "RNA exclusion"
    d.loc[d.category.isin(PERMITTED) & d.mapping_status.eq("mapped"), "assessment"] = (
        "Permitted tissue"
    )
    d.loc[d.category.eq(NEUTRAL_REPRODUCTIVE_IMMUNE), "assessment"] = (
        "Neutral reproductive/developmental immune context"
    )
    d.loc[normal & d.ihc_positive, "assessment"] = "Protein exclusion"
    d["cta_excluded"] = d.assessment.isin(["CTA exclusion", "RNA exclusion", "Protein exclusion"])
    d["reference_review"] = d.assessment.isin(["RNA exclusion", "Protein exclusion"])
    return d

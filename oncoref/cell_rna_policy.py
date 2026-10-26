"""Experimental, explicit two-threshold assessment of normal-cell RNA.

This module has no default thresholds and does not change released CTA or
forbidden-protein membership. Cutoffs must be supplied for each assay, unit and
anatomical tier. An RNA exclusion is evidence for reference review, not proof of
protein translation, peptide presentation, or toxicity.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

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


def match_negative_ihc(contexts, ihc, *, positive_controls, conflicting_genes=()):
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


def assess(contexts, thresholds, *, ihc_factors):
    """Apply independently specified CTA and reference-review cutoffs.

    `value` is an assay-specific supported statistic: a cluster maximum after
    the read guard, the third highest DICE value, or the second highest distinct
    corneal sample-group mean after cell/count guards. Values absent after these
    guards are NaN, never zero. No reproductive maximum can mask normal RNA.
    Equality fails the respective cutoff. Below-cutoff observations are not a
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
    d.loc[measured & d.value.ge(d.cta_cutoff), "assessment"] = "CTA exclusion"
    d.loc[measured & d.value.ge(d.reference_cutoff), "assessment"] = "RNA exclusion"
    d.loc[d.category.isin(PERMITTED) & d.mapping_status.eq("mapped"), "assessment"] = (
        "Permitted tissue"
    )
    d.loc[normal & d.ihc_positive, "assessment"] = "Protein exclusion"
    d["cta_excluded"] = d.assessment.isin(["CTA exclusion", "RNA exclusion", "Protein exclusion"])
    d["reference_review"] = d.assessment.isin(["RNA exclusion", "Protein exclusion"])
    return d

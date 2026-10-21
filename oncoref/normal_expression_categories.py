# Licensed under the Apache License, Version 2.0.
"""Expression-based Reproductive, Normal and Other groups, with explicit gaps.

Publication nomination and therapeutic eligibility are independent annotations.
These human HPA v23 groups are research screens, not sequence-absence claims.
"""

from __future__ import annotations

import pandas as pd

from .cta_normal_screen import adaptive_normal_rna_limits

POLICY_ID = "oncoref.normal-expression-categories.2026-10-21.1"


def classify_normal_expression(
    observations, *, routine_rna_contexts, routine_ihc_contexts, reviewed_unreliable=False
):
    """Classify one gene using mapped observations and an independent coverage panel.

    Input uses the annotated report observation schema. Reproductive requires the
    approved admission, absolute, selectivity and IHC gates. Other requires no
    reproductive positive and passes absolute/IHC gates without a reproductive
    ratio. Its Disposable/Silent subgroups describe below-threshold expression.
    Missing routine RNA cannot enter either of these permissive groups. Missing
    IHC uses the approved baseline and stays explicitly unassessed for protein.
    """
    rows = observations
    audit = adaptive_normal_rna_limits(
        rows,
        routine_ihc_contexts=routine_ihc_contexts,
        reviewed_unreliable=reviewed_unreliable,
    )
    valid = rows.valid_measurement
    mapped = rows.mapping_status.eq("mapped")
    positive = valid & rows.strict_positive
    reproductive = rows.tissue_category.eq("reproductive") & mapped
    disposable = rows.tissue_category.eq("disposable") & mapped
    screened = rows.tissue_category.isin(["manageable", "protected", "vital"]) & mapped
    rna = rows.modality.eq("bulk_rna")
    ihc = rows.modality.eq("ihc")
    rpos = bool((positive & reproductive).any())
    dpos = bool((positive & disposable).any())
    outside_ihc = rows.loc[positive & ihc & screened]
    above = audit.loc[audit.rna_ntpm.gt(audit.rna_ceiling_ntpm)]
    below = audit.loc[audit.reproductive_selectivity.lt(audit.minimum_reproductive_selectivity)]
    measured_rna = set(zip(rows.loc[rna & valid, "tissue"], rows.loc[rna & valid, "cell_type"]))
    measured_ihc = set(zip(rows.loc[ihc & valid, "tissue"], rows.loc[ihc & valid, "cell_type"]))
    missing_rna = set(routine_rna_contexts) - measured_rna
    missing_ihc = set(routine_ihc_contexts) - measured_ihc
    quality = []
    if not routine_rna_contexts or missing_rna:
        quality.append("incomplete_routine_rna")
    if not valid.all():
        quality.append("unusable_measurement")
    if (positive & (~mapped | rows.screening_tier.eq("unresolved"))).any():
        quality.append("unresolved_positive_anatomy")
    if reviewed_unreliable:
        quality.append("reviewed_unreliable_ihc")
    reasons = []
    if len(above):
        reasons.append("rna_above_tissue_ceiling")
    if len(outside_ihc) and not reviewed_unreliable:
        reasons.append("positive_ihc_outside_permitted")
    if rpos and len(below):
        reasons.append("reproductive_selectivity_below_threshold")
    tested = audit.screening_tier.isin(["manageable", "protected", "vital"])
    if rpos and (not tested.any() or audit.loc[tested, "reproductive_selectivity"].isna().any()):
        quality.append("unresolved_reproductive_selectivity")
    category = (
        "Normal" if reasons else "Unassessed" if quality else "Reproductive" if rpos else "Other"
    )
    subtype = ("Disposable" if dpos else "Silent") if category == "Other" else ""
    observed_normal = bool((positive & screened).any())
    return {
        "category_policy_id": POLICY_ID,
        "category": category,
        "other_subgroup": subtype,
        "category_meaning": "below-threshold expression"
        if category == "Other"
        else "research expression screen",
        "reproductive_positive": rpos,
        "disposable_positive": dpos,
        "normal_expression_observed": observed_normal,
        "conservative_normal_reference": category in {"Normal", "Unassessed"} or observed_normal,
        "sequence_specificity_established": False,
        "reproductive_rna_max_ntpm": audit.reproductive_rna_max_ntpm.max(),
        "rna_routine_missing": len(missing_rna),
        "ihc_routine_missing": len(missing_ihc),
        "ihc_coverage": "unreliable"
        if reviewed_unreliable
        else "unavailable"
        if not measured_ihc
        else "partial"
        if missing_ihc
        else "routine_complete",
        "protein_absence_evidence": "not_established",
        "failure_reasons": ";".join(reasons),
        "quality_flags": ";".join(quality),
        "rna_exclusions": "; ".join(
            f"{r.tissue}: {r.rna_ntpm:g} nTPM (limit {r.rna_ceiling_ntpm:g})"
            for r in above.itertuples()
        ),
        "ihc_exclusions": "; ".join(
            f"{r.tissue} / {r.cell_type}: {r.ihc_level} ({r.ihc_reliability})"
            for r in outside_ihc.itertuples()
        ),
        "selectivity_exclusions": "; ".join(
            f"{r.tissue}: {r.reproductive_selectivity:.4%} (minimum {r.minimum_reproductive_selectivity:.0%})"
            for r in below.itertuples()
        )
        if rpos
        else "",
        "non_disposable_observed_rna_max": pd.to_numeric(
            rows.loc[rna & screened & valid, "rna_ntpm"]
        ).max(),
    }

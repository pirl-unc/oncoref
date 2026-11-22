# Licensed under the Apache License, Version 2.0.
"""Expression-based Reproductive, Normal and Other groups, with explicit gaps.

Publication nomination and therapeutic eligibility are independent annotations.
These human HPA v23 groups are research screens, not sequence-absence claims.
"""

from __future__ import annotations

import pandas as pd

from .cta_normal_screen import adaptive_normal_rna_limits

POLICY_ID = "oncoref.normal-expression-categories.2026-10-21.2"


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


RECONCILED_NORMAL_PROFILES: tuple[str, ...] = (
    "unsafe_normal",
    "review_or_unscreened",
    "antigen_candidate",
    "adult_disposable_only",
    "disposable_above_cap",
    "below_threshold_silent",
    "sub_ceiling_non_candidate",
)


def reconcile_with_target_selection(catalog_df: pd.DataFrame | None = None) -> pd.DataFrame:
    """Reconcile legacy ``Normal / Reproductive / Other(Disposable, Silent) / Unassessed``
    semantics onto active ``target_selection.results()`` without relabeling ``Other`` or
    weakening ``Unsafe`` precedence.

    Every row of ``target_selection.results()`` receives a deterministic
    ``reconciled_normal_profile`` and ``legacy_analog``:

    - ``"unsafe_normal"`` (legacy ``"Normal"``): ``category == "Unsafe"``. Always takes
      precedence regardless of any Disposable or reproductive expression.
    - ``"review_or_unscreened"`` (legacy ``"Unassessed"``): ``Other:Review``,
      ``Other:Sequence``, or ``Other:Unscreened tissue`` — missing/ambiguous evidence or
      unmeasured home tissue; never collapsed into ``Silent``.
    - ``"antigen_candidate"`` (legacy ``"Reproductive_or_Cancer_Antigen"``): qualified
      ``CTA``, ``Reproductive Antigen only``, or ``Cancer Antigen only``.
    - ``"adult_disposable_only"`` (legacy ``"Other:Disposable"``): ``Other:Expression``
      passing all shared checks and ``DISPOSABLE_CAP`` (32 TPM / 32 nTPM), with measured
      adult Disposable expression ``>= DETECTION_FLOOR`` (1.0) and non-reproductive
      normal reference ``< DETECTION_FLOOR`` on both native TPM and HPA nTPM scales.
    - ``"disposable_above_cap"`` (legacy ``"Other:Disposable_exceeds_cap"``):
      ``Other:Expression`` where ``disposable_max_tpm > DISPOSABLE_CAP`` or
      ``disposable_max_ntpm > DISPOSABLE_CAP``.
    - ``"below_threshold_silent"`` (legacy ``"Other:Silent"``): ``Other:Expression``
      passing all shared checks with both non-reproductive normal and Disposable
      expression ``< DETECTION_FLOOR`` (below-threshold screen, not proof of zero
      molecules or sequence absence).
    - ``"sub_ceiling_non_candidate"`` (legacy ``"Other:SubCeiling"``): remaining
      ``Other:Expression`` rows that pass tissue-tier ceilings (e.g. Manageable/Protected
      ceiling > 1) with sub-ceiling normal expression ``>= DETECTION_FLOOR`` or failed
      target-peak bounds.
    """
    from . import target_selection

    frame = target_selection.results() if catalog_df is None else catalog_df.copy()
    floor = float(target_selection.DETECTION_FLOOR)
    cap = float(target_selection.DISPOSABLE_CAP)

    disp_tpm = pd.to_numeric(frame["disposable_max_tpm"], errors="coerce").fillna(0.0)
    disp_ntpm = pd.to_numeric(frame["disposable_max_ntpm"], errors="coerce").fillna(0.0)
    norm_tpm = pd.to_numeric(frame["normal_reference_tpm"], errors="coerce").fillna(0.0)
    norm_ntpm = pd.to_numeric(frame["healthy_maximum_ntpm"], errors="coerce").fillna(0.0)

    is_unsafe = frame["category"].eq("Unsafe")
    is_candidate = frame["is_target_candidate"].astype(bool)
    is_other_expr = frame["classification"].eq("Other:Expression")
    shared_screen_ok = (
        frame["normal_screen_pass"].astype(bool)
        & frame["required_coverage_complete"].astype(bool)
        & frame["identity_translation_pass"].astype(bool)
        & frame["evidence_resolved"].astype(bool)
        & frame["sequence_pass"].astype(bool)
    )

    disp_positive = (disp_tpm >= floor) | (disp_ntpm >= floor)
    disp_over_cap = (disp_tpm > cap) | (disp_ntpm > cap)
    non_disp_quiet = (norm_tpm < floor) & (norm_ntpm < floor)

    profile = pd.Series("review_or_unscreened", index=frame.index, dtype="string")
    profile.loc[is_other_expr & shared_screen_ok] = "sub_ceiling_non_candidate"
    profile.loc[
        is_other_expr & shared_screen_ok & ~disp_over_cap & ~disp_positive & non_disp_quiet
    ] = "below_threshold_silent"
    profile.loc[
        is_other_expr & shared_screen_ok & ~disp_over_cap & disp_positive & non_disp_quiet
    ] = "adult_disposable_only"
    profile.loc[is_other_expr & shared_screen_ok & disp_over_cap] = "disposable_above_cap"
    profile.loc[is_candidate] = "antigen_candidate"
    # Unsafe precedence is absolute and applied last so no expression profile overrides it.
    profile.loc[is_unsafe] = "unsafe_normal"

    analog_map = {
        "unsafe_normal": "Normal",
        "review_or_unscreened": "Unassessed",
        "antigen_candidate": "Reproductive_or_Cancer_Antigen",
        "adult_disposable_only": "Other:Disposable",
        "disposable_above_cap": "Other:Disposable_exceeds_cap",
        "below_threshold_silent": "Other:Silent",
        "sub_ceiling_non_candidate": "Other:SubCeiling",
    }

    out = pd.DataFrame(
        {
            "entity_id": frame["entity_id"],
            "protein_group_label": frame["protein_group_label"],
            "member_gene_ids": frame["member_gene_ids"],
            "target_selection_category": frame["category"],
            "target_selection_classification": frame["classification"],
            "reconciled_normal_profile": profile,
            "legacy_analog": profile.map(analog_map).astype("string"),
            "disposable_max_tpm": disp_tpm,
            "disposable_max_ntpm": disp_ntpm,
            "normal_reference_tpm": norm_tpm,
            "healthy_maximum_ntpm": norm_ntpm,
            "expression_bounds_pass": frame["expression_bounds_pass"].astype(bool),
            "in_default_subtraction_background": frame["category"].isin(["Unsafe", "Other"]),
            "sequence_specificity_established": False,
        }
    )
    out.attrs.update(
        policy_id=target_selection.POLICY_ID,
        legacy_category_policy_id=POLICY_ID,
        detection_floor=floor,
        disposable_cap=cap,
        note=(
            "Reconciles legacy Disposable/Silent categories onto target_selection without "
            "relabeling Other or removing Other from default peptide-subtraction backgrounds."
        ),
    )
    return out.reset_index(drop=True)


def adult_disposable_only_proteins(catalog_df: pd.DataFrame | None = None) -> pd.DataFrame:
    """Return screened non-antigen proteins whose only >=1 TPM/nTPM normal expression is
    in adult Disposable tissues within ``target_selection.DISPOSABLE_CAP`` (32 TPM / 32 nTPM).
    """
    reconciled = reconcile_with_target_selection(catalog_df)
    subset = reconciled.loc[
        reconciled["reconciled_normal_profile"].eq("adult_disposable_only")
    ].reset_index(drop=True)
    subset.attrs.update(reconciled.attrs)
    return subset

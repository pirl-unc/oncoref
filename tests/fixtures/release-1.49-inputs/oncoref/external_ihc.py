"""Context-specific adjudication of external normal-tissue staining evidence.

External confidence is our explicit research policy, never an HPA reliability
label. Collection, review and eligibility are distinct. Reagent resellers and
reused images are not independent evidence, and WB validation is not IHC validation.
"""

from __future__ import annotations

import json

import pandas as pd

from .gene_matching import attach_gene_identities
from .load_dataset import get_data

POLICY = "oncoref.external-ihc.2026-10-27.1"
CONTEXT = ["gene_id", "tissue", "cell_type", "stage", "specimen_state"]
IHC_VALIDATION = {"genetic_ihc", "independent_epitopes_ihc", "orthogonal_same_specimens_ihc"}
REQUIRED = {
    "evidence_id",
    *CONTEXT,
    "taxonomy_id",
    "antibody_id",
    "reagent_lineage",
    "source_url",
    "source_locator",
    "source_kind",
    "assay",
    "result",
    "review_status",
    "positive_control",
    "specificity_validation",
    "validation_locator",
    "exact_cell_coverage",
    "normal_donors",
    "independent_of_hpa",
    "notes",
}


def evidence() -> pd.DataFrame:
    """Reviewed observations and explicit limitations; no implied missing negatives."""
    frame = get_data("external-ihc-evidence")
    return assess_evidence(frame)


def assess_evidence(frame: pd.DataFrame) -> pd.DataFrame:
    """Assess eligibility using the same rules for every gene and source.

    A context-matched negative can supersede an Uncertain HPA observation only
    with an identified independent reagent, successful positive control,
    specificity validation in tissue IHC, exact cell coverage and at least two
    normal adult human donors. Two donors is a stated review-policy minimum,
    not a calibrated guarantee of absence. Reviewed positives require the same
    specificity/context controls and at least one normal donor. Missing counts,
    culture/animal/tumor-only evidence and unverified supplier claims cannot qualify.
    """
    missing = REQUIRED - set(frame)
    if missing:
        raise ValueError(f"Missing external-IHC fields: {sorted(missing)}")
    d = frame.copy()
    if d.evidence_id.isna().any() or not d.evidence_id.is_unique:
        raise ValueError("Unique external evidence IDs required")
    if not d.result.isin(
        {"detected", "not_detected", "mixed", "reported_cross_reactivity", "not_assessed"}
    ).all():
        raise ValueError("Unknown staining result")
    if not d.review_status.isin({"reviewed", "pending"}).all():
        raise ValueError("Unknown review status")
    for column in ["positive_control", "exact_cell_coverage", "independent_of_hpa"]:
        if not d[column].isin([True, False]).all():
            raise ValueError(f"Explicit booleans required: {column}")
    for column in CONTEXT[1:]:
        d[column] = d[column].fillna("").str.strip().str.lower()
    donors = pd.to_numeric(d.normal_donors, errors="raise")
    if donors.lt(0).any() or (donors.notna() & donors.mod(1).ne(0)).any():
        raise ValueError("Nonnegative integer donor counts or missing required")
    d = attach_gene_identities(d)
    checks = {
        "no_interpretable_protein_result": ~d.result.isin(["detected", "not_detected", "mixed"]),
        "not_reviewed": d.review_status.ne("reviewed"),
        "not_human_adult_normal_tissue": ~(
            d.taxonomy_id.eq(9606) & d.stage.eq("adult") & d.specimen_state.eq("normal")
        ),
        "not_tissue_ihc": d.assay.ne("IHC"),
        "reagent_identity_or_independence_unresolved": (
            d.antibody_id.fillna("").eq("")
            | d.reagent_lineage.fillna("").eq("")
            | ~d.independent_of_hpa
        ),
        "no_positive_control": ~d.positive_control,
        "specificity_not_validated_in_ihc": ~d.specificity_validation.isin(IHC_VALIDATION),
        "validation_provenance_missing": d.validation_locator.fillna("").eq(""),
        "exact_cell_coverage_missing": ~d.exact_cell_coverage
        | d.cell_type.isin(["", "unspecified"]),
        "insufficient_normal_donors": donors.isna()
        | donors.lt(d.result.eq("not_detected").astype(int) + 1),
        "source_provenance_missing": d.source_url.fillna("").eq("")
        | d.source_locator.fillna("").eq(""),
    }
    failures = pd.DataFrame(checks)
    d["eligible"] = ~failures.any(axis=1)
    d["ineligible_reasons"] = failures.apply(lambda r: " | ".join(r.index[r]), axis=1)
    d["external_confidence"] = "documented"
    d.loc[d.review_status.eq("pending"), "external_confidence"] = "pending"
    d.loc[d.eligible, "external_confidence"] = "context_validated"
    d["external_policy"] = POLICY
    return d


def reconcile_hpa(contexts: pd.DataFrame, external: pd.DataFrame) -> pd.DataFrame:
    """Return an auditable context decision, retaining every original HPA field.

    Contexts require the five CONTEXT keys, taxonomy_id, Level and Reliability.
    Positive credible evidence is never erased by negative evidence. Whole-organ
    negatives do not match named cell subsets. A resolved context is not a gene-
    wide waiver: callers must account for every remaining hold and coverage gap.
    """
    required = {*CONTEXT, "taxonomy_id", "Level", "Reliability"}
    if not required <= set(contexts):
        raise ValueError("Explicit HPA context, stage, state and reliability required")
    if not contexts.taxonomy_id.eq(9606).all():
        raise ValueError("HPA reconciliation requires human observations")
    d = contexts.copy()
    for column in CONTEXT[1:]:
        d[column] = d[column].fillna("").str.strip().str.lower()
    e = assess_evidence(external)
    by_context = dict(tuple(e.groupby(CONTEXT, dropna=False)))
    decisions = []
    for _, row in d.iterrows():
        matched = by_context.get(tuple(row[c] for c in CONTEXT))
        positive, negative, identifiers = False, False, []
        if matched is not None:
            usable = matched[matched.eligible]
            if "Antibody" in row and pd.notna(row.Antibody):
                usable = usable[
                    ~usable.antibody_id.eq(row.Antibody) & ~usable.reagent_lineage.eq(row.Antibody)
                ]
            positive = usable.result.isin(["detected", "mixed"]).any()
            negative = usable.result.eq("not_detected").any()
            identifiers = sorted(usable.evidence_id)
        reliable_positive = row.Level in {
            "Low",
            "Medium",
            "High",
            "Ascending",
            "Descending",
        } and row.Reliability in {"Approved", "Supported", "Enhanced"}
        if (positive or reliable_positive) and negative:
            decision = "conflicting_validated_evidence"
        elif positive:
            decision = "normal_expression_supported"
        elif reliable_positive:
            decision = "retain_hpa_positive"
        elif negative and row.Reliability == "Uncertain":
            decision = "uncertain_superseded_by_negative"
        elif negative:
            decision = "external_negative_coverage"
        else:
            decision = "unresolved"
        decisions.append(
            {
                "external_decision": decision,
                "external_evidence_ids": json.dumps(identifiers),
                "external_policy": POLICY,
            }
        )
    return pd.concat([d.reset_index(drop=True), pd.DataFrame(decisions)], axis=1)

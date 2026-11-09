"""Cancer and Reproductive Antigen sets, with a shared eligibility screen.

This is the next-policy evaluator, not a replacement for the frozen catalog
returned by ``cta_selection.results()``. Inputs must be assessed under this
contract; legacy >0 tumor presence and 5-unit reproductive flags are not inputs.
"""

from __future__ import annotations

from copy import deepcopy

import numpy as np
import pandas as pd

POLICY_ID = "oncoref.target-selection.2026-11-09.1"
DETECTION_FLOOR = 1.0
SINGLE_TUMOR_FOLD = 5.0
REPEATED_TUMORS = 2
EVIDENCE_STATUSES = ("qualifies", "does_not_qualify", "unassessed")
CATEGORIES = ("CTA", "Cancer Antigen only", "Reproductive Antigen only", "Other", "Unsafe")
SHARED_CHECKS = (
    "normal_screen_pass",
    "required_coverage_complete",
    "identity_translation_pass",
    "evidence_resolved",
    "sequence_pass",
)


def definition() -> dict:
    """Return the versioned, machine-readable next-policy contract."""
    return {
        "policy_id": POLICY_ID,
        "status": "evaluator_available_catalog_rebuild_pending",
        "species": "human",
        "categories": list(CATEGORIES),
        "shared_checks": list(SHARED_CHECKS),
        "sets": {
            "Cancer Antigens": "eligible AND cancer evidence qualifies",
            "Reproductive Antigens": "eligible AND reproductive/developmental evidence qualifies",
            "CTA": "Cancer Antigens INTERSECTION Reproductive Antigens",
            "Target candidates": "Cancer Antigens UNION Reproductive Antigens",
        },
        "precedence": ["Unsafe", "Other:Review", "Other:Sequence", "set memberships"],
        "other_subtypes": {
            "Review": "Required evidence is missing, unresolved, or fails a shared admission check.",
            "Sequence": "Nested protein without independent support or >=95% Unsafe 9-mer containment.",
            "Expression": "Shared checks pass; neither positive expression route qualifies.",
        },
        "cancer_rna": {
            "detection_floor": DETECTION_FLOOR,
            "independent_participants_required": REPEATED_TUMORS,
            "repeated_rule": "T >= 1 AND T > N in at least two independent participants",
            "single_tumor_fold": SINGLE_TUMOR_FOLD,
            "single_rule": "T >= 5 * max(N, 1) in at least one qualified tumor sample",
            "normal_reference": (
                "Maximum QC-passing, identity-resolved measurement in the declared compatible "
                "normal panel, excluding specifically allowed Disposable and reproductive compartments."
            ),
            "comparability": (
                "One documented comparison scale and complete required normal contexts. Equal unit "
                "names alone do not establish comparability. No TPM/nTPM/pTPM/CP10K pooling."
            ),
            "missing_reference": "unassessed, never zero",
            "alternative": "Reviewed reliable tumor protein with exact identity and source scope.",
        },
        "reproductive_evidence": {
            "branches": ["adult reproductive", "placental", "extraembryonic", "embryonic", "fetal"],
            "rna_detection_floor": DETECTION_FLOOR,
            "bulk": "Mapped human reproductive compartment >=1 native TPM/nTPM.",
            "cell": "Mapped reproductive cells >=1 native TPM/nTPM and >=20 reads.",
            "developmental": (
                "Human stage/lineage mean >=1 TPM, >=20 cells, >=5 cells at >=1 TPM, "
                "and >=10% cells at >=1 TPM. Retain embryos, donors and stage."
            ),
            "protein": "Reliable mapped reproductive IHC or reviewed human reproductive protein.",
            "concentration": "Descriptive measurement; not an admission gate for either set.",
        },
        "shared_normal_policy": {
            "policy_id": "oncoref.cta-selection.2026-11-08.1",
            "scope": "Normal assessment only; legacy cancer, reproductive and concentration gates excluded.",
            "allowances": (
                "Same resolved reproductive and Disposable compartment allowances for both sets. "
                "Immune, endothelial and other protected cells do not inherit an organ's allowance."
            ),
            "threshold_interpretation": (
                "RNA/IHC cutoffs are screening conventions, not measured injury thresholds. "
                "Negative IHC does not establish a calibrated safety margin."
            ),
        },
        "translation": (
            "Qualifying coding reference or supported noncanonical ORF plus complete approved "
            "MANE-first/canonical representative; resolve conflicting coding evidence. "
            "Experimental protein detection is retained separately."
        ),
        "published_target_support": (
            "Positive evidence annotation for a curated, exact-scope target nomination. Preserve "
            "every protein-paper link and its attested claim; biology or background-table coverage "
            "alone is not a target nomination. Publication cannot override an exclusion or hold."
        ),
        "only_meaning": "Membership in one qualified set, not proven absence from the other context.",
    }


def _require_columns(frame: pd.DataFrame, columns: set[str]) -> None:
    if missing := columns - set(frame):
        raise ValueError(f"Missing fields: {sorted(missing)}")


def _require_booleans(frame: pd.DataFrame, columns) -> None:
    for column in columns:
        if not frame[column].map(lambda v: isinstance(v, (bool, np.bool_))).all():
            raise ValueError(f"{column} requires explicit booleans, without missing values")


def classify(facts: pd.DataFrame) -> pd.DataFrame:
    """Assign sets to assessed proteins, retaining every input evidence column.

    Required columns are ``unsafe``, the five ``SHARED_CHECKS``, and
    ``cancer_evidence``/``reproductive_evidence`` (``EVIDENCE_STATUSES``).
    Protein IDs, source links, findings and publication annotations pass through.
    An unassessed positive association does not erase the other set's membership;
    missing required *normal* coverage blocks both sets.
    """
    flags = ("unsafe", *SHARED_CHECKS)
    _require_columns(facts, {*flags, "cancer_evidence", "reproductive_evidence"})
    _require_booleans(facts, flags)
    for field in ("cancer_evidence", "reproductive_evidence"):
        if not facts[field].isin(EVIDENCE_STATUSES).all():
            raise ValueError(f"{field} must be one of {EVIDENCE_STATUSES}")
    d = facts.copy(deep=True)
    eligible = ~d.unsafe & d[list(SHARED_CHECKS)].all(axis=1)
    cancer = eligible & d.cancer_evidence.eq("qualifies")
    reproductive = eligible & d.reproductive_evidence.eq("qualifies")
    d["eligible"] = eligible
    d["is_cancer_antigen"] = cancer
    d["is_reproductive_antigen"] = reproductive
    d["is_cta"] = cancer & reproductive
    d["is_target_candidate"] = cancer | reproductive
    d["category"] = "Other"
    d.loc[cancer, "category"] = "Cancer Antigen only"
    d.loc[reproductive, "category"] = "Reproductive Antigen only"
    d.loc[cancer & reproductive, "category"] = "CTA"
    d.loc[d.unsafe, "category"] = "Unsafe"
    d["subtype"] = ""
    other = d.category.eq("Other")
    d.loc[other, "subtype"] = "Expression"
    d.loc[other & ~d.sequence_pass, "subtype"] = "Sequence"
    review = ~d[list(SHARED_CHECKS[:-1])].all(axis=1) | (
        eligible & d[["cancer_evidence", "reproductive_evidence"]].eq("unassessed").any(axis=1)
    )
    d.loc[other & review, "subtype"] = "Review"
    d["classification"] = d.category
    d.loc[other, "classification"] = "Other:" + d.loc[other, "subtype"]
    d["target_policy_id"] = POLICY_ID
    return d


def assess_cancer_rna(
    tumors: pd.DataFrame,
    normals: pd.DataFrame,
    *,
    required_normal_contexts: set[str],
    comparison_id: str,
    comparability_reference: str,
    unit: str,
) -> dict:
    """Evaluate one protein in one reviewed tumor/normal comparison panel.

    Both frames require ``entity_id, source_id, comparison_id, unit, value,
    qc_pass, identity_resolved``. Tumors also require ``sample_id, participant_id``;
    normals require ``context_id, scope``. Scope is ``normal``, ``disposable`` or
    ``reproductive`` and must come from reviewed compartment mappings, never from
    the organ name alone. Required contexts are the nonexempt normal contexts.

    The receipt retains observed maxima, supporting samples, and missing coverage.
    Missing/failed normal observations cannot lower the comparison to zero.
    Participant IDs must be resolved across source/cohort aliases. An unknown
    participant can support the single-sample route, but cannot count as replication.
    """
    common = {
        "entity_id",
        "source_id",
        "comparison_id",
        "unit",
        "value",
        "qc_pass",
        "identity_resolved",
    }
    _require_columns(tumors, common | {"sample_id", "participant_id"})
    _require_columns(normals, common | {"context_id", "scope"})
    if unit not in {"TPM", "nTPM", "pTPM"}:
        raise ValueError(
            "This rule requires a declared TPM/nTPM/pTPM comparison, not CP10K or counts"
        )
    if not comparison_id.strip() or not comparability_reference.strip():
        raise ValueError("A comparison ID and documented comparability reference are required")
    if not required_normal_contexts or any(
        not isinstance(c, str) or not c.strip() for c in required_normal_contexts
    ):
        raise ValueError("Declare the nonempty required nonexempt normal context panel")
    entity_ids = set()
    for frame in (tumors, normals):
        _require_booleans(frame, ["qc_pass", "identity_resolved"])
        if not frame.unit.eq(unit).all() or not frame.comparison_id.eq(comparison_id).all():
            raise ValueError("Do not combine different units or comparison panels")
        for field in ("entity_id", "source_id"):
            if not frame[field].map(lambda v: isinstance(v, str) and bool(v.strip())).all():
                raise ValueError(f"Every observation requires {field}")
        entity_ids.update(frame.entity_id)
        values = pd.to_numeric(frame.value, errors="raise")
        if ((values < 0) | (~np.isfinite(values) & values.notna())).any():
            raise ValueError("Expression must be finite and nonnegative, or missing")
    if len(entity_ids) != 1:
        raise ValueError("Assess one resolved protein identity at a time")
    if not normals.scope.isin(["normal", "disposable", "reproductive"]).all():
        raise ValueError("Unknown normal compartment scope")
    if not normals.context_id.map(lambda v: isinstance(v, str) and bool(v.strip())).all():
        raise ValueError("Every normal observation requires a context ID")
    if normals.groupby("context_id").scope.nunique().gt(1).any():
        raise ValueError("A normal context has conflicting compartment allowances")
    if not tumors.sample_id.map(lambda v: isinstance(v, str) and bool(v.strip())).all():
        raise ValueError("Every tumor observation requires a sample ID")
    for field in ("participant_id", "value"):
        if tumors.groupby(["source_id", "sample_id"])[field].nunique(dropna=False).gt(1).any():
            raise ValueError(f"Conflicting {field} for the same tumor observation")
    n = normals.copy()
    t = tumors.copy()
    n["value"] = pd.to_numeric(n.value)
    t["value"] = pd.to_numeric(t.value)
    valid_normal = n.scope.eq("normal") & n.qc_pass & n.identity_resolved & n.value.notna()
    reference = n[valid_normal]
    missing = sorted(required_normal_contexts - set(reference.context_id))
    normal_maximum = None if reference.empty else float(reference.value.max())
    receipt = {
        "policy_id": POLICY_ID,
        "entity_id": next(iter(entity_ids)),
        "comparison_id": comparison_id,
        "comparability_reference": comparability_reference,
        "unit": unit,
        "required_normal_contexts": sorted(required_normal_contexts),
        "missing_normal_contexts": missing,
        "normal_maximum": normal_maximum,
        "normal_maximum_observations": reference[reference.value.eq(normal_maximum)].to_dict(
            "records"
        ),
        "single_sample_threshold": None,
        "independent_participants_above_normal": 0,
        "supporting_tumors": [],
        "routes": [],
        "status": "unassessed",
    }
    if missing or normal_maximum is None:
        receipt["reason"] = "Required normal comparison is incomplete"
        return receipt
    threshold = SINGLE_TUMOR_FOLD * max(normal_maximum, DETECTION_FLOOR)
    receipt["single_sample_threshold"] = threshold
    valid_tumor = t.qc_pass & t.identity_resolved & t.value.notna()
    t = t[valid_tumor].drop_duplicates(["source_id", "sample_id"]).copy()
    if t.empty:
        receipt["reason"] = "No qualified tumor observations"
        return receipt
    t["above_normal"] = t.value.ge(DETECTION_FLOOR) & t.value.gt(normal_maximum)
    t["single_sample_pass"] = t.value.ge(threshold)
    known_participant = t.participant_id.map(lambda v: isinstance(v, str) and bool(v.strip()))
    participants = t.loc[t.above_normal & known_participant, "participant_id"].nunique()
    receipt["independent_participants_above_normal"] = int(participants)
    receipt["supporting_tumors"] = t[t.above_normal | t.single_sample_pass].to_dict("records")
    if participants >= REPEATED_TUMORS:
        receipt["routes"].append("two_independent_tumors_above_normal")
    if t.single_sample_pass.any():
        receipt["routes"].append("one_tumor_fivefold_above_normal")
    receipt["status"] = "qualifies" if receipt["routes"] else "does_not_qualify"
    receipt["reason"] = (
        "Cancer RNA route qualifies" if receipt["routes"] else "No cancer RNA route qualifies"
    )
    return deepcopy(receipt)

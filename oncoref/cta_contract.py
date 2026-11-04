"""Versioned meanings of CTA decisions, evidence states and consumer roles.

Classification and sequence roles have separate identities: an evidence-status
correction can change a role without changing the biological classification.
These roles apply to the selected representative protein, not every isoform or
every peptide/HLA product. The sequence bundle must additionally qualify records.
"""

from __future__ import annotations

from copy import deepcopy

import pandas as pd

CONTRACT_ID = "oncoref.cta-decision-contract.2026-11-04.1"
ROLE_POLICY_ID = "oncoref.sequence-roles.2026-11-04.1"

_CONTRACT = {
    "schema": "oncoref.cta-decision-contract",
    "schema_version": 1,
    "contract_id": CONTRACT_ID,
    "role_policy_id": ROLE_POLICY_ID,
    "title": "CTA and oncofetal candidate selection",
    "definition": (
        "An OncoRef CTA candidate is a reference-resolved protein with supported "
        "cancer expression and qualifying reproductive or developmental association "
        "that passes the named normal-tissue evidence, coverage and selection rules. "
        "Sequence realization requires separately qualified translation records."
    ),
    "categories": {
        "CTA": "Candidate under the named policy; not presentation or clinical clearance",
        "Unsafe": "Normal-tissue or applicable clinical evidence establishes policy exclusion",
        "Other": "Not admitted; inspect every reason and assessment status",
    },
    "assessment_states": {
        "screened": "Declared assessment completed; not exhaustive biological coverage",
        "incomplete_evidence": "Required coverage or evidence resolution is incomplete",
        "normal_evidence_hold": "Normal evidence requires review without establishing exclusion",
        "unassessed": "No applicable selection policy or no completed assessment",
    },
    "interpretation_states": {
        "accepted": "The specified observation is supported within its recorded scope",
        "unresolved": "A material identity, assay, context or interpretation question remains",
        "rejected": "The claimed observation is not supported; retain the rationale and source",
    },
    "interpretation_note": (
        "Review completion, interpretation and policy effect are separate. A completed "
        "review may retain unresolved evidence; a rejected claim is not proof of absence."
    ),
    "support_routes": ["testicular", "other_reproductive", "placental", "embryonic_fetal"],
    "support_route_note": "Routes may overlap; a developmental route is not adult testis restriction",
    "sequence_roles": {
        "exclude": "Retain applicable exclusion evidence; it takes precedence across occurrences",
        "review": "Evidence, identity, coverage or sequence qualification needs resolution",
        "allow": "No exclusion from this assessed record under the named policy",
    },
    "role_precedence": ["exclude", "review", "allow"],
    "role_scope": "selected_representative_protein",
    "translation_scope": "This decision interface does not qualify every translation or isoform",
    "consumer_scope": (
        "Use explicitly qualified records as a reference for exact subtraction or "
        "approximate near-self screening. Preserve every source occurrence and coverage "
        "gap. Patient-HLA prediction, similarity cutoffs and intervention rules are "
        "downstream policy; exact matching does not validate approximate matching."
    ),
    "prohibited_shortcuts": [
        "Blanket exemption of historical CTA nominations from normal-reference evidence",
        "A second legacy IHC filter that silently drops exclusions from other evidence",
        "Interpreting Other, absent records or unmeasured panels as negative normal evidence",
        "Transfer of human numerical rules to another species by symbol or orthology",
    ],
}

# Stable IDs describe independent questions. Display order never determines an
# entity's final class; all applicable evidence remains in the decision audit.
_GATES = (
    ("universe", "Reference and protein universe", "identity"),
    ("normal_bulk", "Normal bulk", "normal_evidence"),
    ("normal_clusters", "HPA tissue clusters", "normal_evidence"),
    ("normal_hpa_immune", "HPA immune", "normal_evidence"),
    ("normal_monaco", "Monaco immune", "normal_evidence"),
    ("normal_dice", "DICE immune", "normal_evidence"),
    ("normal_cornea", "Adult cornea", "normal_evidence"),
    ("normal_hpa_protein", "HPA normal protein", "normal_evidence"),
    ("normal_reviewed_protein", "Accepted normal-protein evidence", "normal_evidence"),
    ("normal_knowledge", "Reviewed normal-cell knowledge", "normal_evidence"),
    ("normal_external_protein", "External validated normal protein", "normal_evidence"),
    ("clinical_exclusion", "Clinical safety exclusion", "clinical_evidence"),
    ("rare_vital", "Rare vital-tissue cell expression", "normal_evidence"),
    ("rare_protected", "Rare Protected-cell expression review", "normal_evidence"),
    ("required_coverage", "Required RNA coverage", "coverage"),
    ("cancer_presence", "Cancer expression", "candidate_identity"),
    ("reproductive_support", "Reproductive or developmental expression", "candidate_identity"),
    ("reproductive_concentration", "Reproductive concentration", "candidate_identity"),
    ("nested_microprotein", "Nested microprotein", "sequence_identity"),
    ("unsafe_containment", "Unsafe sequence containment", "sequence_identity"),
    ("translation_qualification", "Translation qualification", "sequence_identity"),
    ("final_assignment", "CTA and Other assignments", "decision"),
)


def contract() -> dict:
    """Return detached definitions; numerical policy and source receipts are separate."""
    return deepcopy(_CONTRACT)


def gates() -> pd.DataFrame:
    """Return stable question IDs, not claims that every assay was measured."""
    return pd.DataFrame(_GATES, columns=["gate_id", "name", "axis"])


def sequence_roles(frame: pd.DataFrame) -> pd.DataFrame:
    """Assign representative-protein roles without erasing incomplete evidence.

    Unsafe exclusion has precedence. All other records need a recognized,
    completed assessment and resolved identity before receiving allowance.
    Missing or novel status values require review, regardless of Other subtype.
    This is not a complete translation or all-occurrence sequence qualification.
    """
    required = {"entity_id", "category", "assessment_status", "identity_status"}
    if not required <= set(frame):
        raise ValueError(f"Missing decision fields: {sorted(required - set(frame))}")
    if frame.entity_id.isna().any() or frame.entity_id.duplicated().any():
        raise ValueError("Unique nonmissing decision identities required")
    if not frame.category.dropna().isin(["CTA", "Unsafe", "Other"]).all():
        raise ValueError("Unknown decision category")
    out = frame.copy(deep=True)
    excluded = out.category.eq("Unsafe").fillna(False)
    resolved = out.identity_status.eq("resolved").fillna(False)
    complete = out.assessment_status.eq("screened").fillna(False)
    assessed = out.category.isin(["CTA", "Other"]) & resolved & complete
    out["direct_target_candidate"] = out.category.eq("CTA").fillna(False) & assessed
    out["shared_sequence"] = pd.Series("review", index=out.index, dtype="string")
    out.loc[assessed, "shared_sequence"] = "allow"
    out.loc[excluded, "shared_sequence"] = "exclude"

    def reasons(row):
        result = []
        if pd.notna(row.category) and row.category == "Unsafe":
            result.append("applicable_policy_exclusion")
        if pd.isna(row.category):
            result.append("selection_unassessed")
        if pd.isna(row.assessment_status) or row.assessment_status != "screened":
            result.append("assessment_requires_review")
        if pd.isna(row.identity_status) or row.identity_status != "resolved":
            result.append("identity_requires_review")
        return tuple(result) or ("assessed_without_policy_exclusion",)

    out["sequence_role_reasons"] = [reasons(row) for row in out.itertuples()]
    out["role_policy_id"] = ROLE_POLICY_ID
    out["role_scope"] = _CONTRACT["role_scope"]
    return out

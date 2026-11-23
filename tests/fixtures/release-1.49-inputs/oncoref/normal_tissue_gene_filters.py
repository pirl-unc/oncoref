# Licensed under the Apache License, Version 2.0.
"""Strict, context-aware gene restriction filters with explicit coverage gaps.

These filters consume assay-specific detection calls, not raw values or shared
cross-species thresholds. Passing means no counterevidence in a declared panel,
not proven absence, pregnancy safety, or preservation of fertility.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

from .normal_tissues import NormalTissueAnnotation, normal_tissue_taxonomy

POLICY_ID = "oncoref.normal-tissue-gene-filters.2026-10-20.1"
_PROFILES = {
    "male_reproductive_maternal_fetal_excluded": {
        "label": "Postnatal male reproductive only; maternal and developmental exclusion",
        "intent": "Pregnancy-oriented candidate restriction, not pregnancy safety",
        "allowed": "Postnatal male reproductive anatomy in male samples only",
        "excluded": "Every other normal context, including fetal male reproductive tissues",
        "required_negative_domains": (
            "female_somatic",
            "female_reproductive",
            "pregnant_maternal",
            "embryonic",
            "fetal",
            "fetal_male_gonad",
            "fetal_female_gonad",
            "placental",
            "other_extraembryonic",
        ),
    },
    "female_developmental_male_reproductive_excluded": {
        "label": "Female/developmental only; male reproductive exclusion",
        "intent": "Male reproductive preservation candidate restriction, not proven fertility preservation",
        "allowed": "Female reproductive samples plus embryonic, fetal and extraembryonic contexts",
        "excluded": "Male reproductive anatomy at every stage; other postnatal contexts",
        "required_negative_domains": (
            "male_reproductive",
            "male_somatic",
            "male_mammary",
        ),
    },
}


def normal_tissue_gene_filter_profiles() -> dict:
    """Detached, explicitly named gene-filter definitions."""
    return {key: {"policy_id": POLICY_ID, **value} for key, value in _PROFILES.items()}


@dataclass(frozen=True)
class NormalTissueEvidence:
    """One assay call with source context, gene and detection-rule provenance.

    ``context_id`` distinguishes source samples/stage bins within a dataset.
    ``not_detected`` must be an assessed negative under ``detection_policy_id``;
    missing values, absent rows, and unreliable IHC must be ``unassessed``.
    Gene IDs must be namespace-qualified; orthologs are evaluated separately.
    """

    gene_id: str
    annotation: NormalTissueAnnotation
    context_id: str
    assay: str
    detection_policy_id: str
    result: str
    evidence_id: str


@dataclass(frozen=True)
class NormalTissueCoverageRequirement:
    """A required source/sample/assay context, fixed independently of the gene."""

    annotation: NormalTissueAnnotation
    context_id: str
    assay: str
    detection_policy_id: str


def _key(row):
    a = row.annotation
    return (
        a.taxon_id,
        a.source_id,
        a.source_version,
        a.source_tissue,
        a.source_cell_type,
        a.source_stage,
        a.stage,
        a.sex,
        a.physiological_state,
        row.context_id,
        row.assay,
        row.detection_policy_id,
    )


def _domains(a):
    scopes = set(a.scopes)
    domains = set()
    if a.mapping_status not in {"mapped", "stage_only"}:
        return domains
    for name in ("embryonic", "fetal", "placental"):
        if name in scopes:
            domains.add(name)
    if "extraembryonic" in scopes and "placental" not in scopes:
        domains.add("other_extraembryonic")
    if a.stage == "fetal":
        if "testis" in scopes and a.sex == "male":
            domains.add("fetal_male_gonad")
        if "ovary" in scopes and a.sex == "female":
            domains.add("fetal_female_gonad")
    if a.stage == "postnatal" and not scopes & {"extraembryonic", "embryonic", "fetal"}:
        if a.sex == "female" and "male_reproductive" not in scopes:
            domains.add(
                "female_reproductive" if "female_reproductive" in scopes else "female_somatic"
            )
            if a.physiological_state == "pregnant":
                domains.add("pregnant_maternal")
        if a.sex == "male" and not ("female_reproductive" in scopes and "mammary" not in scopes):
            domains.add("male_reproductive" if "male_reproductive" in scopes else "male_somatic")
            if "mammary" in scopes:
                domains.add("male_mammary")
    return domains


def _allowed(a, profile):
    """True/False/None: admitted, excluded, or missing required metadata."""
    if a.mapping_status not in {"mapped", "stage_only"}:
        return None
    scopes = set(a.scopes)
    developmental = bool(scopes & {"extraembryonic", "embryonic", "fetal"})
    if profile == "male_reproductive_maternal_fetal_excluded":
        if developmental or "male_reproductive" not in scopes or a.sex == "female":
            return False
        return True if a.stage == "postnatal" and a.sex == "male" else None
    # An explicit male reproductive compartment vetoes the larger fetal union.
    if "male_reproductive" in scopes:
        return False
    if developmental:
        return True
    if "female_reproductive" not in scopes or a.sex == "male":
        return False
    return True if a.stage == "postnatal" and a.sex == "female" else None


def evaluate_normal_tissue_gene_filter(
    observations,
    *,
    profile,
    required_contexts=(),
    coverage_panel_id=None,
) -> dict:
    """Evaluate one gene against all observations and a declared negative panel.

    Every measured positive is considered, including sources outside the panel.
    One positive overrides any number of negatives. Negative calls match the
    exact source version, source labels, sex, stage, sample, assay and detection
    rule; they cannot clear another source or species. Required domains are a
    minimum sanity check, not a complete tissue/cell/gestational-age inventory.
    An empty or incomplete panel cannot pass. No default pregnancy panel is
    bundled: adult HPA data alone cannot establish maternal/fetal exclusion.
    """
    if profile not in _PROFILES:
        raise ValueError(f"Unknown normal-tissue gene filter {profile!r}")
    observations, required_contexts = tuple(observations), tuple(required_contexts)
    if any(not isinstance(x, NormalTissueEvidence) for x in observations):
        raise TypeError("Expected NormalTissueEvidence observations")
    if any(not isinstance(x, NormalTissueCoverageRequirement) for x in required_contexts):
        raise TypeError("Expected NormalTissueCoverageRequirement contexts")
    if len({x.gene_id for x in observations}) > 1:
        raise ValueError("Evaluate one gene at a time; do not pool genes or orthologs")
    rows = (*observations, *required_contexts)
    if any(not isinstance(x.annotation, NormalTissueAnnotation) for x in rows):
        raise TypeError("Expected NormalTissueAnnotation contexts")
    if len({x.annotation.taxon_id for x in rows}) > 1:
        raise ValueError("Cannot combine different species")
    current = normal_tissue_taxonomy()["taxonomy_id"]
    if any(x.annotation.taxonomy_id != current for x in rows):
        raise ValueError("Resolve observations against the current taxonomy")
    if any(not x.context_id or not x.assay or not x.detection_policy_id for x in rows):
        raise ValueError("Context, assay and detection policy provenance are required")
    if any(x.result not in {"detected", "not_detected", "unassessed"} for x in observations):
        raise ValueError("Unknown evidence result")
    if any(not x.gene_id or not x.evidence_id for x in observations):
        raise ValueError("Gene and evidence identifiers are required")
    if len({_key(x) for x in required_contexts}) != len(required_contexts):
        raise ValueError("Duplicate coverage requirements")
    excluded, unresolved, positive = [], [], []
    by_context = {}
    for i, row in enumerate(observations):
        by_context.setdefault(_key(row), []).append(row)
        if row.result == "unassessed":
            unresolved.append(i)
        elif row.result == "detected":
            admitted = _allowed(row.annotation, profile)
            (positive if admitted else unresolved if admitted is None else excluded).append(i)
    missing, covered_domains = [], set()
    for i, requirement in enumerate(required_contexts):
        calls = by_context.get(_key(requirement), ())
        valid = (
            requirement.annotation.mapping_status in {"mapped", "stage_only"}
            and _allowed(requirement.annotation, profile) is False
            and bool(calls)
            and all(x.result == "not_detected" for x in calls)
        )
        if valid:
            covered_domains.update(_domains(requirement.annotation))
        else:
            missing.append(i)
    missing_domains = sorted(set(_PROFILES[profile]["required_negative_domains"]) - covered_domains)
    panel_complete = (
        bool(required_contexts and coverage_panel_id) and not missing and not missing_domains
    )
    passed = bool(positive) and not excluded and not unresolved and panel_complete
    panel_hash = hashlib.sha256(
        json.dumps(
            {
                "taxonomy_id": current,
                "filter_policy_id": POLICY_ID,
                "profile": profile,
                "contexts": sorted(json.dumps(_key(row)) for row in required_contexts),
            },
            sort_keys=True,
        ).encode()
    ).hexdigest()
    return {
        "profile": profile,
        "filter_policy_id": POLICY_ID,
        "taxonomy_id": current,
        "gene_id": observations[0].gene_id if observations else None,
        "coverage_panel_id": coverage_panel_id,
        "coverage_panel_sha256": panel_hash,
        "status": "excluded"
        if excluded
        else "passes_declared_panel"
        if passed
        else "insufficient_evidence",
        "passes_declared_panel": passed,
        "positive_pattern_compatible": bool(positive) and not excluded and not unresolved,
        "allowed_positive_indices": tuple(positive),
        "excluded_positive_indices": tuple(excluded),
        "unassessed_indices": tuple(unresolved),
        "missing_required_context_indices": tuple(missing),
        "missing_negative_domains": tuple(missing_domains),
        "absence_evidence": "complete_declared_panel"
        if panel_complete
        else "incomplete_or_not_assessed",
        "clinical_safety": "not_assessed",
    }

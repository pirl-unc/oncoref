# Licensed under the Apache License, Version 2.0.
"""Frozen source semantics for CTA snapshot schema v2.

Ingestion and validation use these same rules. Replay checks archived records
against this versioned schema without replacing them or reading current curation.
Changing these semantics requires a new snapshot schema version.
"""

from __future__ import annotations

from copy import deepcopy

from .expression_contract import ExpressionMeasurement

# Freeze the v2 vocabulary independently of future changes to live CTA policy.
HUMAN_REPRODUCTIVE_TISSUES_V2 = frozenset(
    {
        "cervix",
        "endometrium",
        "endometrium 1",
        "endometrium 2",
        "epididymis",
        "fallopian tube",
        "ovary",
        "placenta",
        "prostate",
        "seminal vesicle",
        "testis",
        "vagina",
    }
)


def summary_records_v2(evidence, contract, source):
    """Expose native summary values separately from the saved human policy.

    Organ maxima can overlap the overall somatic maximum. They are summaries of
    one source, never extra samples or independent confirmation. Protein tissue
    flags are derived curation fields, not measured negative IHC observations.
    """
    if evidence["source_id"] != "oncoref-human-cta":
        return None
    if contract["taxon_id"] != 9606 or contract["reference_status"] != "legacy_human":
        raise ValueError("Historical human source requires its human summary contract")
    row = evidence["source_fields"]
    measurement = ExpressionMeasurement("gene", "nTPM", "source_native", False).as_dict()
    columns = [
        ("rna_testis_ntpm", "testis", "reproductive", "tissue_summary"),
        ("rna_ovary_ntpm", "ovary", "reproductive", "tissue_summary"),
        ("rna_placenta_ntpm", "placenta", "reproductive", "tissue_summary"),
        (
            "rna_max_somatic_ntpm",
            row.get("rna_max_somatic_tissue") or None,
            "somatic",
            "somatic_maximum",
        ),
        *[
            (f"rna_{organ}_max_ntpm", tissue, "somatic", "tissue_maximum")
            for organ, tissue in (
                ("brain", "brain safety-tissue group (including retina)"),
                ("heart", "heart muscle"),
                ("lung", "lung"),
                ("liver", "liver"),
                ("pancreas", "pancreas"),
            )
        ],
    ]
    facts = []
    for field, tissue, role, statistic in columns:
        if field not in row:
            raise ValueError(f"Historical CTA summary lacks source field {field}")
        value = row[field]
        if value is not None and (type(value) not in (int, float) or value < 0):
            raise ValueError(f"Historical CTA summary has invalid nTPM in {field}")
        facts.append(
            {
                "id": f"{evidence['id']}:{field}",
                "evidence_id": evidence["id"],
                "source_field": field,
                "tissue": tissue,
                "tissue_role": role,
                "statistic": statistic,
                "measurement": deepcopy(measurement),
                "observation": {
                    "status": "not_reported" if value is None else "reported",
                    "value": value,
                    "reported": "not reported" if value is None else str(value),
                    "positive_count": None,
                },
            }
        )
    if "rna_max_somatic_tissue" not in row or "protein_strict_expression" not in row:
        raise ValueError("Historical CTA summary lacks tissue source fields")
    protein = row["protein_strict_expression"]
    if protein is not None and not isinstance(protein, str):
        raise ValueError("Historical CTA protein summary must be text or null")
    tissues = (
        [t.strip() for t in protein.split(";") if t.strip()]
        if protein and protein.lower() not in {"no data", "not detected", "nan"}
        else []
    )
    for index, tissue in enumerate(tissues or [None]):
        facts.append(
            {
                "id": f"{evidence['id']}:protein_strict_expression:{index}",
                "evidence_id": evidence["id"],
                "source_field": "protein_strict_expression",
                "tissue": tissue,
                "tissue_role": "unknown"
                if tissue is None
                else ("reproductive" if tissue in HUMAN_REPRODUCTIVE_TISSUES_V2 else "somatic"),
                "statistic": "tissue_detection_summary",
                "measurement": ExpressionMeasurement(
                    "protein_group", None, "categorical IHC summary", None
                ).as_dict(),
                "observation": {
                    "status": "qualitative" if tissues else "not_reported",
                    "value": None,
                    "reported": protein or "not reported",
                    "positive_count": None,
                },
            }
        )
    policy_fields = (
        "specificity_status",
        "specificity_action",
        "specificity_source_anchor",
        "specificity_rationale",
        "passes_filters",
        "never_expressed",
        "restriction",
        "restriction_confidence",
        "protein_restriction",
        "protein_reliability",
        "protein_strict_expression",
        "protein_testis",
        "protein_ovary",
        "protein_placenta",
        "rna_restriction",
        "rna_restriction_level",
        "rna_reproductive_tissue_scope",
    )
    if any(field not in row for field in policy_fields):
        raise ValueError("Historical CTA curation lacks source fields")
    action = row["specificity_action"]
    curation = {
        "id": f"{evidence['id']}:historical-core-policy",
        "evidence_id": evidence["id"],
        "policy_id": "oncoref-human-cta-core",
        "policy_version": source["version"],
        "decision_field": "specificity_action",
        "included_values": ["include_default"],
        "decision": "unavailable"
        if not action
        else ("included" if action == "include_default" else "excluded"),
        "rule": "Historical HPA restriction, never-expressed rescue and source-anchored specificity reviews; frozen source decision, not a shared admission threshold.",
        "source_fields": {name: row.get(name) for name in policy_fields},
    }
    return facts, [curation]

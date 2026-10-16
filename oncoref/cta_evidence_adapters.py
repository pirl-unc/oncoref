# Licensed under the Apache License, Version 2.0.
"""Source adapters for the common CTA facts contract.

Only this ingestion layer knows the historical human column names. It reads
records already frozen in the snapshot: no current datasets or network access.
Native measurements are never converted to a common numerical scale.
"""

from __future__ import annotations

from copy import deepcopy
from functools import lru_cache

from .cta_tissues import EXTENDED_REPRODUCTIVE_TISSUES
from .expression_contract import ExpressionMeasurement


def native_observation_facts(evidence, contract, source):
    """Keep a study observation at its original assay and candidate scope."""
    return [
        {
            "id": f"{evidence['id']}:observation",
            "evidence_id": evidence["id"],
            "source_field": None,
            "tissue": evidence["tissue"],
            "tissue_role": evidence["tissue_role"],
            "statistic": "native_observation",
            "measurement": deepcopy(contract["measurement"]),
            "observation": deepcopy(evidence["observation"]),
        }
    ], []


def human_summary_facts(evidence, contract, source):
    """Expose native summary values separately from the saved human policy.

    Organ maxima can overlap the overall somatic maximum. They are summaries of
    one source, never extra samples or independent confirmation. Protein tissue
    flags are derived curation fields, not measured negative IHC observations.
    """
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
        value = row.get(field)
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
    protein = row.get("protein_strict_expression")
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
                else (
                    "reproductive"
                    if tissue in EXTENDED_REPRODUCTIVE_TISSUES | {"endometrium 1", "endometrium 2"}
                    else "somatic"
                ),
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
    action = row.get("specificity_action")
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


@lru_cache(maxsize=4)
def normalize_cta_snapshot(snapshot):
    """Freeze source-adapted facts in a v2 snapshot; existing v2 data is unchanged.

    Explicitly calling this on an archived v1 snapshot uses only that archive's
    records. Reading an archive alone never reinterprets it or reads current data.
    """
    from .species_antigens import CTAEvidenceSnapshot

    payload = snapshot.as_dict()
    if payload["schema_version"] == 2:
        return snapshot
    contracts = {row["contract_id"]: row for row in payload["contracts"]}
    sources = {row["id"]: row for row in payload["sources"]}
    payload["schema_version"] = 2
    payload["facts"] = []
    payload["curation"] = []
    for evidence in payload["evidence"]:
        contract = contracts[evidence["contract_id"]]
        adapter = native_observation_facts
        if evidence["source_id"] == "oncoref-human-cta":
            if payload["taxon_id"] != 9606 or contract["reference_status"] != "legacy_human":
                raise ValueError("Historical human source adapter requires its human contract")
            adapter = human_summary_facts
        facts, curation = adapter(evidence, contract, sources[evidence["source_id"]])
        payload["facts"].extend(facts)
        payload["curation"].extend(curation)
    payload.pop("snapshot_id")
    return CTAEvidenceSnapshot.create(payload)

# Licensed under the Apache License, Version 2.0.
"""Shared CTA fact validation and coverage reporting, independent of species.

The model keeps observations, candidate identity and historical policy separate.
It reports evidence gaps without imposing downstream admission thresholds.
"""

from __future__ import annotations

from collections import defaultdict

from .expression_contract import ExpressionMeasurement

FACT_FIELDS = {
    "id",
    "evidence_id",
    "source_field",
    "tissue",
    "tissue_role",
    "statistic",
    "measurement",
    "observation",
}
CURATION_FIELDS = {
    "id",
    "evidence_id",
    "policy_id",
    "policy_version",
    "decision",
    "rule",
    "source_fields",
    "decision_field",
    "included_values",
}
STATISTICS = {
    "native_observation",
    "tissue_summary",
    "somatic_maximum",
    "tissue_maximum",
    "tissue_detection_summary",
}
TISSUE_ROLES = {"somatic", "reproductive", "tumor", "mixed", "unknown"}


def _text(value, name):
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError(f"{name} must be nonempty text")


def _records(rows, fields, name):
    if not isinstance(rows, list):
        raise ValueError(f"{name} must be a list")
    seen = set()
    for row in rows:
        if not isinstance(row, dict) or set(row) != fields:
            raise ValueError(f"{name} requires its exact schema fields")
        _text(row["id"], f"{name} id")
        if row["id"] in seen:
            raise ValueError(f"Duplicate {name} id")
        seen.add(row["id"])
    return rows


def validate_cta_facts(record):
    """Validate frozen v2 facts against their parent evidence and source fields."""
    evidence = {r["id"]: r for r in record["evidence"]}
    contracts = {r["contract_id"]: r for r in record["contracts"]}
    for fact in _records(record["facts"], FACT_FIELDS, "CTA fact"):
        if fact["evidence_id"] not in evidence:
            raise ValueError("CTA fact has unknown evidence")
        parent = evidence[fact["evidence_id"]]
        if fact["tissue"] is not None:
            _text(fact["tissue"], "fact tissue")
        if fact["tissue_role"] not in TISSUE_ROLES or fact["statistic"] not in STATISTICS:
            raise ValueError("Unknown CTA fact tissue role/statistic")
        ExpressionMeasurement.from_dict(fact["measurement"])
        observation = fact["observation"]
        if not isinstance(observation, dict) or set(observation) != {
            "status",
            "value",
            "reported",
            "positive_count",
        }:
            raise ValueError("CTA fact requires a complete observation")
        if observation["status"] not in {
            "reported",
            "below_limit",
            "no_amplification",
            "not_reported",
            "qualitative",
        }:
            raise ValueError("Unknown CTA fact observation status")
        _text(observation["reported"], "reported fact")
        value = observation["value"]
        if value is not None and (
            type(value) not in (int, float)
            or observation["status"] not in {"reported", "below_limit"}
        ):
            raise ValueError("Missing/qualitative facts cannot have numeric values")
        if observation["status"] == "below_limit" and value is None:
            raise ValueError("Censored facts require their reported bound")
        count = observation["positive_count"]
        if count is not None and (
            type(count) is not int
            or count < 0
            or parent["sample_count"] is None
            or count > parent["sample_count"]
        ):
            raise ValueError("Fact positive count lacks sampled coverage")
        field = fact["source_field"]
        if field is None:
            contract = contracts[parent["contract_id"]]
            if (
                fact["statistic"] != "native_observation"
                or observation != parent["observation"]
                or fact["measurement"] != contract["measurement"]
                or fact["tissue"] != parent["tissue"]
                or fact["tissue_role"] != parent["tissue_role"]
            ):
                raise ValueError("Native CTA fact disagrees with its source observation")
        else:
            _text(field, "source field")
            source_value = parent["source_fields"].get(field)
            if fact["statistic"] == "tissue_detection_summary":
                if (
                    value is not None
                    or count is not None
                    or observation["status"] not in {"qualitative", "not_reported"}
                    or observation["reported"] != (source_value or "not reported")
                    or (
                        fact["tissue"] is not None
                        and (
                            not isinstance(source_value, str)
                            or fact["tissue"] not in [t.strip() for t in source_value.split(";")]
                        )
                    )
                ):
                    raise ValueError("Tissue detection fact disagrees with its source summary")
                continue
            if (
                fact["statistic"] == "native_observation"
                or value != source_value
                or observation["positive_count"] is not None
                or observation["status"] != ("not_reported" if source_value is None else "reported")
            ):
                raise ValueError("Summary CTA fact disagrees with its source field")
    for curation in _records(record["curation"], CURATION_FIELDS, "CTA curation"):
        if curation["evidence_id"] not in evidence:
            raise ValueError("CTA curation has unknown evidence")
        for key in ("policy_id", "policy_version", "rule"):
            _text(curation[key], key)
        if curation["decision"] not in {"included", "excluded", "unavailable"}:
            raise ValueError("Unknown historical curation decision")
        fields = curation["source_fields"]
        parent = evidence[curation["evidence_id"]]
        if (
            not isinstance(fields, dict)
            or not fields
            or any(value != parent["source_fields"].get(key) for key, value in fields.items())
        ):
            raise ValueError("CTA curation must preserve its source fields")
        field = curation["decision_field"]
        _text(field, "decision field")
        values = curation["included_values"]
        if field not in fields or not isinstance(values, list) or not values:
            raise ValueError("Historical curation requires its decision field and included values")
        for value in values:
            _text(value, "included value")
        action = fields[field]
        expected = (
            "unavailable"
            if action is None or action == ""
            else ("included" if action in values else "excluded")
        )
        if curation["decision"] != expected:
            raise ValueError("Historical curation decision disagrees with its source field")


def index_evidence_records(records):
    """Group fact/curation records by their source evidence without copying them."""
    result = defaultdict(list)
    for row in records:
        result[row["evidence_id"]].append(row)
    return result


def candidate_evidence_coverage(candidate, facts, evidence, contracts, reference, mapping):
    """Describe available observations and precise gaps using the same logic for every species.

    No atlas denominator or independent donor total is declared by the present
    snapshot schema. Keep those unknown; never sum technical replicates, overlapping
    assays, group observations or summary maxima as independent coverage.
    """
    groups = {role: [] for role in ("somatic", "reproductive", "tumor", "other", "not_reported")}
    tissues = set()
    for fact in facts:
        if fact["observation"]["status"] == "not_reported":
            groups["not_reported"].append(fact["id"])
            continue
        role = fact["tissue_role"]
        groups[role if role in groups else "other"].append(fact["id"])
        if fact["tissue"] is not None:
            tissues.add(fact["tissue"])
    missing = ["complete_tissue_panel", "independent_sample_total", "independent_donor_total"]
    if not facts:
        missing.append("normalized_observations")
    if not any(
        f["statistic"] == "native_observation"
        and f["observation"]["status"] in {"no_amplification", "below_limit"}
        for f in facts
    ):
        missing.append("measured_negative_or_censored_observations")
    for role in ("somatic", "reproductive", "tumor"):
        if not groups[role]:
            missing.append(f"{role}_observations")
    parents = [evidence[eid] for eid in candidate["evidence_ids"]]
    if not parents or any(
        contracts[p["contract_id"]]["reference_status"] != "verified" for p in parents
    ):
        missing.append("quantification_reference")
    elif any(
        contracts[p["contract_id"]]["reference"]["reference_id"] != reference["reference_id"]
        for p in parents
    ):
        missing.append("matching_quantification_reference")
    if mapping["status"] != "unique" or not mapping["complete"]:
        missing.append("unique_complete_gene_mapping")
    if not {"rna", "protein"} <= set(reference["content_provenance"]):
        missing.append("sequence_reference_assets")
    if not mapping["targets"] or any(
        not target["sequences"]
        or any(not pair["transcript_id"] or not pair["protein_id"] for pair in target["sequences"])
        for target in mapping["targets"]
    ):
        missing.append("transcript_protein_pairs")
    return {
        "status": "partial"
        if any(groups[r] for r in ("somatic", "reproductive", "tumor", "other"))
        else "unavailable",
        "fact_ids": groups,
        "reported_tissues": sorted(tissues),
        "unreported_tissues": None,
        "independent_sample_count": None,
        "independent_donor_count": None,
        "missing": missing,
    }

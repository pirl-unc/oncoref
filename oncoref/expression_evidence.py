# Licensed under the Apache License, Version 2.0.
"""Explicit, bounded assignment evidence through the common expression surface."""

from __future__ import annotations

from . import managed_expression as managed
from .expression_records import ExpressionRecords
from .expression_registry import expression_source_registry_entries


def _entry(expression):
    return next(
        (
            entry
            for entry in expression_source_registry_entries(species=expression.contract.taxon_id)
            if entry["id"] == expression.contract.source_id
        ),
        {},
    )


def evidence_sources(expression):
    result = []
    for profile in _entry(expression).get("evidence_sources", []):
        spec = managed.evidence_specification(profile, expression.contract.taxon_id)
        result.append(
            managed.ManagedExpressionAsset(
                profile["id"],
                expression.contract.taxon_id,
                spec["manifest"]["version"],
                spec,
            ).as_dict()
        )
    return tuple(result)


def with_evidence(
    expression,
    source_id,
    *,
    feature_ids,
    offline,
    cache_root,
    local_overrides,
    max_bytes,
    max_records,
):
    from .species_expression import ExpressionDataset

    if isinstance(feature_ids, (str, bytes)):
        raise TypeError("feature_ids must be an explicit collection of native identifiers")
    identifiers = tuple(feature_ids)
    if (
        not identifiers
        or any(not isinstance(x, str) or not x for x in identifiers)
        or len(set(identifiers)) != len(identifiers)
    ):
        raise ValueError("Unique nonempty feature_ids required")
    if type(max_records) is not int or max_records < 1:
        raise ValueError("max_records must be a positive integer")
    if len(identifiers) > max_records:
        raise ValueError("Feature selection exceeds max_records")
    if not set(identifiers) <= set(expression._features.index):
        raise ValueError("Evidence selection contains an unknown native feature")
    entry = _entry(expression)
    profile = next((p for p in entry.get("evidence_sources", []) if p["id"] == source_id), None)
    if profile is None:
        raise ValueError("No registered evidence source for this expression dataset")
    provenance = expression.provenance
    if provenance is None or provenance.get("producer") != "orfdata":
        raise ValueError("Evidence requires a verified managed expression snapshot")
    if provenance.get("enrichment"):
        raise ValueError(
            "Already enriched; select all required features from the original snapshot"
        )
    expected = managed.specification(entry)
    managed._validate_contract(expected, expression.contract)
    if any(provenance[key] != expected[key] for key in ("manifest", "download", "selection")):
        raise ValueError("Expression provenance differs from its pinned managed source")
    spec = managed.evidence_specification(profile, expression.contract.taxon_id)
    summary = provenance["manifest"]
    if (
        spec["manifest"]["references"] != summary["references"]
        or any(r["taxon_id"] != expression.contract.taxon_id for r in summary["references"])
        or spec["manifest"]["coverage"].get("accession") != provenance.get("experiment")
    ):
        raise ValueError("Evidence taxon/reference/experiment differs from expression")
    if profile["adapter"] != "pfizer-summary":
        raise ValueError("Unsupported expression evidence adapter")
    dataset = managed.open_managed(
        spec,
        taxon_id=expression.contract.taxon_id,
        offline=offline,
        cache_root=cache_root,
        local_overrides=local_overrides,
        max_bytes=max_bytes,
    )
    from .pfizer_canine_protein import assignment_contexts, assignment_records

    contexts, targets = assignment_contexts(expression, dataset, sorted(identifiers))
    retained = {"mappings": {}, "relations": {}}
    count = 0

    def retain(kind, row):
        nonlocal count
        value = row.to_dict()
        key = value["mapping_id" if kind == "mappings" else "relation_id"]
        previous = retained[kind].get(key)
        if previous is not None and previous != value:
            raise ValueError("Conflicting assignment record identities")
        if previous is None:
            count += 1
            if count > max_records:
                raise ValueError(
                    "Evidence exceeds max_records; select fewer features or raise the limit"
                )
            retained[kind][key] = value

    for kind, row in assignment_records(dataset, targets):
        retain(kind, row)
    for kind, added in retained.items():
        id_field = "mapping_id" if kind == "mappings" else "relation_id"
        if any(row[id_field] in added for row in expression.native_records(kind)):
            raise ValueError("Assignment identifiers collide with existing expression records")
    provenance["enrichment"] = {
        "source_id": source_id,
        **spec,
        "feature_ids": sorted(identifiers),
        "contexts": contexts,
        "record_ids": {kind: sorted(rows) for kind, rows in retained.items()},
        "interpretation": "Same-study representative context; no additional measurements, "
        "unique protein attribution or MHC presentation",
        "record_scope": "Assignments mentioning selected representatives, with all candidates "
        "and complete memberships of returned groups; not a full-source export",
    }

    def records(kind):
        yield from expression.native_records(kind)
        for key in sorted(retained.get(kind, {})):
            yield retained[kind][key]

    return ExpressionDataset(
        expression._values,
        expression._features,
        expression._samples,
        expression.contract,
        records=ExpressionRecords(provenance, records),
    )

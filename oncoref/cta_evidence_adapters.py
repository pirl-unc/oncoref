# Licensed under the Apache License, Version 2.0.
"""Source adapters for the common CTA facts contract.

Only this ingestion layer knows the historical human column names. It reads
records already frozen in the snapshot: no current datasets or network access.
Native measurements are never converted to a common numerical scale.
"""

from __future__ import annotations

from copy import deepcopy
from functools import lru_cache

from .cta_evidence_schema import summary_records_v2


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
    """Adapt a historical human summary using its frozen v2 source schema."""
    records = summary_records_v2(evidence, contract, source)
    if records is None:
        raise ValueError("Human summary adapter requires its historical source")
    return records


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

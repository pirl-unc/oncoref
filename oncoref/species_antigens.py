# Licensed under the Apache License, Version 2.0.
"""Offline CTA nominations and evidence, independent of patient admission policy.

The snapshot is immutable; query results are detached JSON records. A mapping to
an annotation is not a claim that an assay resolved that gene or its isoforms.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path

from .cta_evidence_model import (
    candidate_evidence_coverage,
    index_evidence_records,
    validate_cta_facts,
)
from .expression_contract import ExpressionContract, _fields
from .reference_identity import ReferenceIdentity, _canonical_json, _identity, _json_value, _label
from .species import resolve_species

NOMINATION_KINDS = {"published_nomination", "expression_derived", "assay_followup", "legacy_human"}
EVIDENCE_STATUSES = {"limited", "conflicting", "unresolved"}
_SCHEMA = "oncoref.cta-evidence"


def _records(rows, kind):
    if not isinstance(rows, list):
        raise ValueError(f"{kind} must be a list")
    result = {}
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError(f"{kind} must contain records")
        key = _label(row.get("id"), f"{kind} id")
        if key in result:
            raise ValueError(f"Duplicate {kind} id: {key}")
        result[key] = row
    return result


def _count(value, name):
    if value is not None and (type(value) is not int or value < 0):
        raise ValueError(f"{name} must be a nonnegative integer or null")


def _validate(record):
    if not isinstance(record, dict) or not _json_value(record):
        raise ValueError("Snapshot requires finite JSON-native records")
    version = record.get("schema_version")
    additions = {"facts", "curation"} if version == 2 else set()
    _fields(
        record,
        {
            "schema",
            "schema_version",
            "snapshot_id",
            "taxon_id",
            "source_version",
            "coverage",
            "references",
            "sources",
            "contracts",
            "candidates",
            "evidence",
        }
        | additions,
        "CTA snapshot",
    )
    if (
        record["schema"] != _SCHEMA
        or type(record["schema_version"]) is not int
        or record["schema_version"] not in (1, 2)
    ):
        raise ValueError("Unsupported CTA snapshot schema/version")
    taxon = resolve_species(record["taxon_id"]).taxon_id
    if type(record["taxon_id"]) is not int:
        raise ValueError("Snapshot taxon_id must be an integer")
    _label(record["source_version"], "source_version")
    coverage = record["coverage"]
    _fields(coverage, {"status", "limitations"}, "CTA coverage")
    if coverage["status"] not in {"partial", "unavailable"} or not isinstance(
        coverage["limitations"], list
    ):
        raise ValueError("CTA coverage must be explicitly partial or unavailable")
    for limitation in coverage["limitations"]:
        _label(limitation, "coverage limitation")
    refs = {}
    for row in _records(record["references"], "reference").values():
        _fields(row, {"id", "identity"}, "CTA reference")
        ref = ReferenceIdentity.from_dict(row["identity"])
        if ref.taxon_id != taxon or ref.reference_id in refs:
            raise ValueError("Duplicate reference identity or reference taxon mismatch")
        refs[ref.reference_id] = ref
    sources = _records(record["sources"], "source")
    for source in sources.values():
        _fields(
            source,
            {"id", "version", "url", "doi", "pmid", "sha256", "scope", "study_id", "lineage"},
            "CTA source",
        )
        for key in ("version", "url", "scope", "study_id"):
            _label(source[key], key)
        for key in ("doi", "pmid"):
            if source[key] is not None:
                _label(source[key], key)
        if not isinstance(source["lineage"], list) or not source["lineage"]:
            raise ValueError("Source study lineage is required")
        for ancestor in source["lineage"]:
            _label(ancestor, "lineage")
        digest = source["sha256"]
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or set(digest) - set("0123456789abcdef")
        ):
            raise ValueError("Source requires a SHA-256 digest")
    contracts = {}
    if not isinstance(record["contracts"], list):
        raise ValueError("contracts must be a list")
    for row in record["contracts"]:
        contract = ExpressionContract.from_dict(row)
        if contract.taxon_id != taxon or contract.source_id not in sources:
            raise ValueError("Expression contract species/source mismatch")
        if contract.contract_id in contracts:
            raise ValueError("Duplicate expression contract")
        if (
            contract.source_version != sources[contract.source_id]["version"]
            or contract.content_provenance.get("source") != sources[contract.source_id]["sha256"]
        ):
            raise ValueError("Expression contract/source provenance mismatch")
        contracts[contract.contract_id] = contract
    candidates = _records(record["candidates"], "candidate")
    evidence = _records(record["evidence"], "evidence")
    if candidates and not refs:
        raise ValueError("Candidate snapshots require an explicit target reference")
    if (coverage["status"] == "unavailable") != (not candidates and not evidence):
        raise ValueError("Coverage status disagrees with snapshot contents")
    for row in candidates.values():
        _fields(
            row,
            {
                "id",
                "label",
                "nomination_kind",
                "derivation",
                "evidence_status",
                "mappings",
                "evidence_ids",
                "notes",
            },
            "CTA candidate",
        )
        _label(row["label"], "candidate label")
        _label(row["notes"], "candidate notes")
        if (
            row["nomination_kind"] not in NOMINATION_KINDS
            or row["evidence_status"] not in EVIDENCE_STATUSES
        ):
            raise ValueError("Unknown CTA nomination/evidence status")
        if row["nomination_kind"] == "legacy_human" and taxon != 9606:
            raise ValueError("Legacy human curation cannot be applied to another species")
        if not isinstance(row["derivation"], list) or not row["derivation"]:
            raise ValueError("Candidate derivation is required")
        for derivation in row["derivation"]:
            _fields(derivation, {"kind", "source_id", "locator", "note"}, "Candidate derivation")
            if (
                derivation["kind"]
                not in {
                    "published_nomination",
                    "orthology_nomination",
                    "expression_derived",
                    "assay_followup",
                    "legacy_human",
                }
                or derivation["source_id"] not in sources
            ):
                raise ValueError("Unknown candidate derivation/source")
            for key in ("locator", "note"):
                _label(derivation[key], key)
        if row["nomination_kind"] == "expression_derived" and not any(
            d["kind"] == "expression_derived" for d in row["derivation"]
        ):
            raise ValueError("Expression-derived membership requires an explicit derivation")
        if not isinstance(row["mappings"], list):
            raise ValueError("mappings must be a list")
        seen = set()
        for mapping in row["mappings"]:
            _fields(
                mapping,
                {
                    "reference_id",
                    "gene_namespace",
                    "status",
                    "method",
                    "complete",
                    "targets",
                    "notes",
                },
                "Candidate mapping",
            )
            rid = mapping["reference_id"]
            if rid not in refs or rid in seen:
                raise ValueError("Unknown or duplicate mapping reference")
            seen.add(rid)
            if mapping["gene_namespace"] != refs[rid].identifier_namespaces["gene"]:
                raise ValueError("Mapping gene namespace disagrees with reference")
            if type(mapping["complete"]) is not bool:
                raise ValueError("Mapping completeness must be boolean")
            for key in ("method", "notes"):
                _label(mapping[key], key)
            targets = _records(mapping["targets"], "mapping target")
            status = mapping["status"]
            if (
                status not in {"unique", "ambiguous", "unresolved"}
                or (status == "unresolved") != (not targets)
                or (status == "unique" and len(targets) != 1)
                or (status == "unresolved" and mapping["complete"])
            ):
                raise ValueError("Mapping status disagrees with its targets")
            for target in targets.values():
                _fields(target, {"id", "symbol", "sequences"}, "Mapping target")
                _label(target["symbol"], "target symbol")
                if not isinstance(target["sequences"], list):
                    raise ValueError("sequences must be a list")
                for sequence in target["sequences"]:
                    _fields(sequence, {"transcript_id", "protein_id"}, "Sequence mapping")
                    if not any(sequence.values()):
                        raise ValueError("Sequence mapping requires an accession")
                    for key, value in sequence.items():
                        if value is not None:
                            _label(value, key)
            # A missing mapping is different from a resolved empty gene set.
        if seen != set(refs):
            raise ValueError("Every candidate needs a mapping status for every reference")
        ids = row["evidence_ids"]
        if (
            not isinstance(ids, list)
            or not ids
            or len(ids) != len(set(ids))
            or any(eid not in evidence for eid in ids)
        ):
            raise ValueError("Invalid candidate evidence links")
        if any(row["id"] not in evidence[eid]["candidate_ids"] for eid in ids):
            raise ValueError("Candidate/evidence links must be reciprocal")
    for row in evidence.values():
        _fields(
            row,
            {
                "id",
                "candidate_ids",
                "source_id",
                "study_id",
                "contract_id",
                "locator",
                "assay",
                "tissue",
                "tissue_role",
                "sample_ids",
                "sample_count",
                "donor_count",
                "technical_replicates",
                "observation",
                "source_fields",
                "notes",
            },
            "CTA evidence",
        )
        if not isinstance(row["source_fields"], dict):
            raise ValueError("source_fields must preserve a source record")
        contract = contracts.get(row["contract_id"])
        if contract is None or contract.source_id != row["source_id"]:
            raise ValueError("Evidence contract/source mismatch")
        if row["study_id"] != sources[row["source_id"]]["study_id"]:
            raise ValueError("Evidence study lineage mismatch")
        ids = row["candidate_ids"]
        if (
            not isinstance(ids, list)
            or not ids
            or len(ids) != len(set(ids))
            or any(
                cid not in candidates or row["id"] not in candidates[cid]["evidence_ids"]
                for cid in ids
            )
        ):
            raise ValueError("Invalid evidence/candidate links")
        if row["tissue_role"] not in {"somatic", "reproductive", "tumor", "mixed", "unknown"}:
            raise ValueError("Unknown tissue role")
        for key in ("locator", "assay", "notes"):
            _label(row[key], key)
        if row["tissue"] is not None:
            _label(row["tissue"], "tissue")
        if not isinstance(row["sample_ids"], list):
            raise ValueError("sample_ids must be a list")
        for sample in row["sample_ids"]:
            _label(sample, "sample id")
        for key in ("sample_count", "donor_count", "technical_replicates"):
            _count(row[key], key)
        observation = row["observation"]
        _fields(
            observation, {"status", "value", "reported", "positive_count"}, "Evidence observation"
        )
        if observation["status"] not in {
            "reported",
            "below_limit",
            "no_amplification",
            "not_reported",
            "qualitative",
        }:
            raise ValueError("Unknown observation status")
        value = observation["value"]
        if value is not None and (
            type(value) not in (int, float)
            or observation["status"] not in {"reported", "below_limit"}
        ):
            raise ValueError("Missing/qualitative observations cannot have numeric values")
        if observation["status"] == "below_limit" and value is None:
            raise ValueError("A censored observation requires its reported bound")
        _label(observation["reported"], "reported observation")
        _count(observation["positive_count"], "positive_count")
        if observation["positive_count"] is not None and (
            row["sample_count"] is None or observation["positive_count"] > row["sample_count"]
        ):
            raise ValueError("Positive count exceeds or lacks sampled coverage")
    if version == 2:
        validate_cta_facts(record)
    payload = {k: v for k, v in record.items() if k != "snapshot_id"}
    if record["snapshot_id"] != _identity("cta-evidence", version, payload):
        raise ValueError("CTA snapshot digest does not match its contents")


@dataclass(frozen=True)
class CTAEvidenceSnapshot:
    """Content-addressed evidence with strict, offline JSON replay."""

    _serialized: str

    def __post_init__(self):
        try:
            _validate(json.loads(self._serialized))
        except (KeyError, TypeError) as exc:
            raise ValueError("Malformed CTA snapshot record") from exc

    @classmethod
    def from_dict(cls, record):
        if not _json_value(record):
            raise ValueError("Snapshot requires finite JSON-native records")
        return cls(_canonical_json(record))

    @classmethod
    def create(cls, payload):
        """Validate and seal a newly built snapshot payload (without snapshot_id)."""
        if not isinstance(payload, dict) or "snapshot_id" in payload:
            raise ValueError("A new CTA snapshot payload must omit snapshot_id")
        return cls.from_dict(
            {
                **payload,
                "snapshot_id": _identity("cta-evidence", payload.get("schema_version"), payload),
            }
        )

    def as_dict(self):
        return json.loads(self._serialized)

    @property
    def snapshot_id(self):
        return self.as_dict()["snapshot_id"]

    @property
    def taxon_id(self):
        return self.as_dict()["taxon_id"]

    @property
    def coverage(self):
        return self.as_dict()["coverage"]

    def write(self, path):
        """Export the entire evidence/annotation snapshot without external assets."""
        Path(path).write_text(self._serialized + "\n", encoding="utf-8")

    @classmethod
    def read(cls, path, *, expected_snapshot_id=None):
        result = cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))
        if expected_snapshot_id is not None and result.snapshot_id != expected_snapshot_id:
            raise ValueError("Unexpected CTA snapshot identity")
        return result

    def _reference(self, reference, data):
        refs = data["references"]
        if reference is None:
            if not refs:
                return None
            if len(refs) == 1:
                return ReferenceIdentity.from_dict(refs[0]["identity"])
            raise ValueError("Select an explicit reference from this snapshot")
        if isinstance(reference, ReferenceIdentity):
            if reference.taxon_id != data["taxon_id"]:
                raise ValueError("CTA query/reference taxon mismatch")
            match = next(
                (r for r in refs if r["identity"]["reference_id"] == reference.reference_id), None
            )
        elif isinstance(reference, str):
            match = next((r for r in refs if r["id"] == reference), None)
        else:
            raise ValueError("An explicit ReferenceIdentity or snapshot reference name is required")
        if match is None:
            raise ValueError("CTA snapshot has no mapping for the requested annotation")
        return ReferenceIdentity.from_dict(match["identity"])

    def candidates(self, *, reference=None):
        """All nominations, including ambiguous/unresolved and counterevidence rows."""
        data = self.as_dict()
        ref = self._reference(reference, data)
        return self._candidates(data, ref)

    def _candidates(self, data, ref, candidate_ids=None):
        facts = index_evidence_records(data.get("facts", []))
        evidence = {r["id"]: r for r in data["evidence"]}
        contracts = {r["contract_id"]: r for r in data["contracts"]}
        result = []
        for row in data["candidates"]:
            if candidate_ids is not None and row["id"] not in candidate_ids:
                continue
            mapping = next(m for m in row.pop("mappings") if m["reference_id"] == ref.reference_id)
            result.append(
                {
                    **row,
                    "taxon_id": data["taxon_id"],
                    "snapshot_id": data["snapshot_id"],
                    "reference": ref.as_dict(),
                    "mapping": mapping,
                    "gene_ids": [t["id"] for t in mapping["targets"]],
                    "evidence_coverage": candidate_evidence_coverage(
                        row,
                        [fact for eid in row["evidence_ids"] for fact in facts[eid]],
                        evidence,
                        contracts,
                        ref.as_dict(),
                        mapping,
                    ),
                }
            )
        return tuple(result)

    def gene_ids(self, *, reference):
        """Unambiguously mapped nominee IDs; inspect candidates() for alternatives."""
        if reference is None:
            raise ValueError("Gene IDs require an explicit reference")
        return frozenset(
            gid
            for row in self.candidates(reference=reference)
            if (row["mapping"]["status"] == "unique" and row["mapping"]["complete"])
            for gid in row["gene_ids"]
        )

    def evidence(self, *, reference=None, candidate_id=None):
        """Native observations, common facts and separately named historical curation."""
        data = self.as_dict()
        ref = self._reference(reference, data)
        if candidate_id is not None and candidate_id not in {r["id"] for r in data["candidates"]}:
            raise ValueError(f"Unknown CTA candidate {candidate_id!r}")
        rows = [
            r
            for r in data["evidence"]
            if candidate_id is None or candidate_id in r["candidate_ids"]
        ]
        ids = {cid for row in rows for cid in row["candidate_ids"]}
        candidates = {r["id"]: r for r in self._candidates(data, ref, ids)}
        sources = {s["id"]: s for s in data["sources"]}
        contracts = {c["contract_id"]: c for c in data["contracts"]}
        facts = index_evidence_records(data.get("facts", []))
        curation = index_evidence_records(data.get("curation", []))
        return tuple(
            {
                **row,
                "taxon_id": data["taxon_id"],
                "snapshot_id": data["snapshot_id"],
                "source": sources[row["source_id"]],
                "contract": contracts[row["contract_id"]],
                "candidates": [candidates[cid] for cid in row["candidate_ids"]],
                "facts": facts[row["id"]],
                "curation": curation[row["id"]],
            }
            for row in rows
        )


def _snapshot_path(species):
    return files("oncoref").joinpath(f"data/species/{species.name}/cta-evidence.json")


def cta_snapshot(species, *, path=None, expected_snapshot_id=None):
    """Load a packaged or exported snapshot without downloading or human fallback."""
    selected = resolve_species(species)
    if path is not None:
        snapshot = CTAEvidenceSnapshot.read(path, expected_snapshot_id=expected_snapshot_id)
    elif selected.taxon_id == 9606:
        from .human_cta_snapshot import human_cta_snapshot

        snapshot = human_cta_snapshot()
    elif _snapshot_path(selected).is_file():
        snapshot = CTAEvidenceSnapshot.from_dict(
            json.loads(_snapshot_path(selected).read_text(encoding="utf-8"))
        )
    else:
        payload = {
            "schema": _SCHEMA,
            "schema_version": 1,
            "taxon_id": selected.taxon_id,
            "source_version": "unavailable-v1",
            "references": [],
            "sources": [],
            "contracts": [],
            "candidates": [],
            "evidence": [],
            "coverage": {
                "status": "unavailable",
                "limitations": [
                    "No versioned CTA evidence snapshot is curated for this species. Historical human APIs are unchanged."
                ],
            },
        }
        snapshot = CTAEvidenceSnapshot.from_dict(
            {**payload, "snapshot_id": _identity("cta-evidence", 1, payload)}
        )
    if snapshot.taxon_id != selected.taxon_id:
        raise ValueError("CTA snapshot/collection taxon mismatch")
    if path is None:
        from .cta_evidence_adapters import normalize_cta_snapshot

        snapshot = normalize_cta_snapshot(snapshot)
    if expected_snapshot_id is not None and snapshot.snapshot_id != expected_snapshot_id:
        raise ValueError("Unexpected CTA snapshot identity")
    return snapshot


def antigen_inventory(species):
    selected = resolve_species(species)
    if not _snapshot_path(selected).is_file():
        return []
    return [
        {
            "name": "cta-evidence",
            "held": "wheel",
            "category": "cta",
            "group": "antigens",
            "taxon_id": selected.taxon_id,
            "species": selected.name,
            "available": True,
            "cohorts": None,
            "description": "Versioned CTA nominations, native observations and counterevidence; not patient admission",
        }
    ]

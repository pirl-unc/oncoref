# Licensed under the Apache License, Version 2.0.
"""Exact protein identity and occurrence groups for any explicitly declared reference."""

from __future__ import annotations

import hashlib
import string
from collections import defaultdict
from collections.abc import Mapping
from copy import deepcopy
from dataclasses import fields

from .reference_identity import (
    ReferenceIdentity,
    _canonical_json,
    _identity,
    _json_value,
    _label,
    _string_map,
    _text,
)

PROTEIN_SEQUENCE_IDENTITY_VERSION = 1
PROTEIN_GROUP_SCHEMA_VERSION = 1
_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
_SEQUENCE_INPUT = _AMINO_ACIDS | frozenset("acdefghiklmnpqrstvwy*" + string.whitespace)


def normalize_protein_sequence(sequence: str) -> str:
    """v1: validate ASCII AA20 input, remove ASCII whitespace/one final *, uppercase.

    I/L remain distinct. Reject ambiguity codes, U/O, internal stops, empty input,
    non-ASCII whitespace and Unicode lookalikes BEFORE any case conversion. This
    identifies sequence only; callers establish translation completeness separately.
    """
    if not _text(sequence):
        raise ValueError("missing or non-string protein sequence")
    invalid = sorted(set(sequence) - _SEQUENCE_INPUT)
    if invalid:
        raise ValueError(f"invalid protein sequence characters={invalid!r}")
    normalized = (
        sequence.translate(str.maketrans("", "", string.whitespace)).removesuffix("*").upper()
    )
    if not normalized or set(normalized) - _AMINO_ACIDS:
        raise ValueError("empty protein sequence or internal/repeated stop")
    return normalized


def protein_sequence_id(sequence: str) -> str:
    """Typed v1 SHA-256 of normalized ASCII residues, independent of any reference."""
    digest = hashlib.sha256(normalize_protein_sequence(sequence).encode("ascii")).hexdigest()
    return f"openvax:protein-sequence:v{PROTEIN_SEQUENCE_IDENTITY_VERSION}:sha256:{digest}"


def _check_scope_fields(row, scope):
    for name in (field.name for field in fields(ReferenceIdentity)):
        if name not in row:
            continue
        value = row[name]
        if name == "content_provenance":
            value = _string_map(value, name, digests=True)
        if _canonical_json(value) != _canonical_json(scope[name]):
            raise ValueError(f"conflicting reference field {name}")


def protein_sequence_groups(occurrences, *, reference: ReferenceIdentity):
    """Group complete exact proteins; retain all source records in defensive copies.

    Every row must declare this full ``reference_id`` and nonempty string
    occurrence_id/gene_id/transcript_id/protein_id/source/sequence, plus complete=True.
    Gene/transcript/protein namespaces must be explicit in the reference. Optional
    redundant scope fields must agree. Original sequence and all JSON-native source
    fields are retained. Invalid rows raise ValueError naming the occurrence/reason;
    no partial grouping is returned. Exact duplicate records are deduplicated;
    contradictory reuse of an occurrence ID fails, including JSON type differences.

    Sequence IDs describe residues across references. Group IDs bind those IDs to
    a reference under a separate versioned schema. Output ordering is deterministic.
    Grouping establishes no RNA allocation, locus-specific translation, PTM-defined
    proteoform, reproductive restriction or target-panel membership.
    """
    if not isinstance(reference, ReferenceIdentity):
        raise ValueError("Explicit ReferenceIdentity required")
    if not {"gene", "transcript", "protein"} <= reference.identifier_namespaces.keys():
        raise ValueError("Gene, transcript and protein identifier namespaces required")
    scope = reference.as_dict()
    rid = reference.reference_id
    groups, sequences, seen = defaultdict(list), {}, {}
    for occurrence in occurrences:
        if not isinstance(occurrence, Mapping):
            raise ValueError("Occurrence must be a mapping")
        row = dict(occurrence)
        oid = row.get("occurrence_id")

        def reject(reason, oid=oid):
            raise ValueError(f"Occurrence {oid!r}: {reason}")

        if not _json_value(row):
            reject("source fields must contain finite JSON-native acyclic values")
        required = ("occurrence_id", "gene_id", "transcript_id", "protein_id", "source")
        if any(not _text(row.get(key)) for key in required):
            reject("incomplete occurrence identity or source provenance")
        if row.get("reference_id") != rid:
            reject("conflicting or missing reference_id")
        try:
            _check_scope_fields(row, scope)
        except ValueError as error:
            reject(str(error))
        if row.get("complete") is not True:
            reject("incomplete translation (complete=True required)")
        try:
            sequence = normalize_protein_sequence(row.get("sequence"))
        except ValueError as error:
            reject(str(error))
        if oid in seen:
            if _canonical_json(seen[oid]) != _canonical_json(row):
                reject("contradictory reuse of occurrence identity")
            continue
        row = deepcopy(row)
        seen[oid] = row
        sid = protein_sequence_id(sequence)
        sequences[sid] = sequence
        groups[sid].append(row)
    result = []
    for sid, rows in sorted(groups.items()):
        identity = {
            "schema": "openvax.protein-group",
            "schema_version": PROTEIN_GROUP_SCHEMA_VERSION,
            "reference_id": rid,
            "protein_sequence_id": sid,
        }
        result.append(
            {
                **identity,
                "group_id": _identity("protein-group", PROTEIN_GROUP_SCHEMA_VERSION, identity),
                "reference": reference.as_dict(),
                "sequence": sequences[sid],
                "occurrences": sorted(rows, key=lambda row: row["occurrence_id"]),
            }
        )
    return result


def adapt_reference_occurrences(
    occurrences,
    *,
    reference: ReferenceIdentity,
    source_reference,
    source_reference_id,
    reference_id_field,
    field_map,
    added_fields,
    evidence,
):
    """Explicit schema migration, preserving original records and mapping evidence.

    ``field_map`` maps target ReferenceIdentity fields to source-reference fields;
    ``added_fields`` explicitly supplies missing declarations. All target fields
    must be accounted for exactly once and must reproduce ``reference``. Taxon,
    assembly, annotation release and content hashes must come from source fields.
    ``evidence`` pins reviewed migration inputs by logical role -> SHA-256.

    Source IDs are opaque legacy IDs: callers verify their legacy derivation; this
    adapter checks each occurrence carries the declared source ID. It performs no
    taxon conversion, liftover, identifier crosswalk or reference-equivalence inference.
    Original records (including old IDs) remain in source_record; a migration ID
    links each output row to the returned audit envelope. Persist that envelope.
    """
    if not isinstance(reference, ReferenceIdentity):
        raise ValueError("Explicit ReferenceIdentity required")
    _label(source_reference_id, "source_reference_id")
    _label(reference_id_field, "reference_id_field")
    if not isinstance(source_reference, Mapping) or not _json_value(dict(source_reference)):
        raise ValueError("Source reference must be a finite JSON-native record")
    names = {field.name for field in fields(reference)}
    if not isinstance(field_map, Mapping) or not isinstance(added_fields, Mapping):
        raise ValueError("Explicit field_map and added_fields mappings required")
    required_mapped = {"taxon_id", "assembly_accession", "annotation_release", "content_provenance"}
    if (
        set(field_map) & set(added_fields)
        or set(field_map) | set(added_fields) != names
        or not required_mapped <= field_map.keys()
    ):
        raise ValueError("Map source scope/content and declare every reference field exactly once")
    if any(not _text(key) or key not in source_reference for key in field_map.values()):
        raise ValueError("Reference field map names an absent source field")
    mapped = {key: source_reference[value] for key, value in field_map.items()}
    mapped.update(added_fields)
    if ReferenceIdentity(**mapped) != reference:
        raise ValueError("Mapped source reference does not match the target reference")
    proof = _string_map(evidence, "migration evidence", digests=True)
    audit = {
        "schema": "openvax.reference-migration",
        "schema_version": 1,
        "source_reference_id": source_reference_id,
        "source_reference": deepcopy(dict(source_reference)),
        "reference_id_field": reference_id_field,
        "field_map": dict(field_map),
        "added_fields": deepcopy(dict(added_fields)),
        "reference": reference.as_dict(),
        "evidence": proof,
    }
    if not _json_value(audit):
        raise ValueError("Migration declarations must be finite JSON-native values")
    migration_id = _identity("reference-migration", 1, audit)
    output = []
    target_scope = reference.as_dict()
    target_id = reference.reference_id
    for occurrence in occurrences:
        if not isinstance(occurrence, Mapping) or not _json_value(dict(occurrence)):
            raise ValueError("Source occurrence must be a finite JSON-native record")
        row = deepcopy(dict(occurrence))
        if row.get(reference_id_field) != source_reference_id:
            raise ValueError(
                f"Occurrence {row.get('occurrence_id')!r}: conflicting source reference ID"
            )
        if {"source_record", "reference_migration_id"} & row.keys():
            raise ValueError("Occurrence already contains reserved migration fields")
        if reference_id_field != "reference_id" and "reference_id" in row:
            raise ValueError("Occurrence already declares a target reference_id")
        try:
            _check_scope_fields(row, target_scope)
            # Reject conflicting per-row source declarations as well as canonical
            # fields: a schema migration must never relabel another reference.
            _check_scope_fields(
                {name: row[source] for name, source in field_map.items() if source in row},
                target_scope,
            )
        except ValueError as error:
            raise ValueError(f"Occurrence {row.get('occurrence_id')!r}: {error}") from error
        output.append(
            {
                **row,
                "source_record": deepcopy(row),
                "reference_id": target_id,
                "reference_migration_id": migration_id,
            }
        )
    return {**audit, "migration_id": migration_id, "occurrences": output}

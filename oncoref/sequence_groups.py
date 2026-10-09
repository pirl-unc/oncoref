# Licensed under the Apache License, Version 2.0.
"""Exact complete protein sequences within an explicit species and reference.

These are annotation facts, independent of the human CTA registries. Grouping
infers neither RNA/protein abundance, locus-specific translation, PTM-defined
proteoforms, nor canine target-panel membership.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from collections.abc import Mapping
from copy import deepcopy

_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _text(value):
    return isinstance(value, str) and bool(value.strip())


def _json_value(value):
    """Require JSON-native provenance so serialization loses no fields or types."""
    if value is None or type(value) in (str, int, bool):
        return True
    if type(value) is float:
        return value == value and abs(value) != float("inf")
    if type(value) is list:
        return all(_json_value(item) for item in value)
    if type(value) is dict:
        return all(type(key) is str and _json_value(item) for key, item in value.items())
    return False


def protein_sequence_groups(
    occurrences, *, taxon_id, assembly_accession, annotation_release, content_provenance
):
    """Group normalized, exact proteins and retain every distinct source occurrence.

    ``taxon_id`` must be a positive integer. Assembly accession/version and
    annotation release must be explicit nonempty strings. ``content_provenance``
    maps source asset names to their SHA-256 digests, binding reference identity to
    actual content; filenames alone are insufficient. Inputs assert a common
    reference, and each row must carry matching taxon/assembly/annotation fields.
    A row's optional ``content_provenance`` must match the supplied asset map.

    Each occurrence requires nonempty string ``occurrence_id``, ``gene_id``,
    ``transcript_id``, ``protein_id``, ``source``, and ``sequence`` fields, plus
    ``complete=True``. IDs retain their versions. All source fields (including
    coordinates, hashes, and nested provenance) must be JSON-native and are copied.
    Whitespace and case are normalized and one terminal stop is removed; the
    original source sequence is retained in the occurrence. Only the 20 standard
    amino acids are accepted: I/L remain distinct, while ambiguity codes, U/O,
    internal stops and incomplete translations are rejected explicitly.

    Returns a deterministically ordered list of dictionaries with sequence digest
    ``protein_sequence_id``, reference digest ``reference_id``, scoped ``group_id``,
    normalized ``sequence``, reference fields and all source ``occurrences``.
    Equal sequences across references share a sequence digest, but not a group ID.
    JSON round trips preserve the result; use ``sort_keys=True`` for canonical bytes.

    Identical duplicate source records are harmless. Reusing an occurrence ID with
    contradictory source fields raises ``ValueError``, as does any invalid row,
    with its occurrence ID and reason. No partial grouping is returned. Callers
    can record these errors in their own quarantine workflow. Independently
    allocated RNA contributions may be summed downstream; one ambiguous RNA
    observation must never be multiplied by the occurrence count.
    """
    if type(taxon_id) is not int or taxon_id <= 0:
        raise ValueError("Positive integer taxonomy identifier required")
    if not _text(assembly_accession) or not _text(annotation_release):
        raise ValueError("Explicit assembly accession/version and annotation release required")
    if not isinstance(content_provenance, Mapping) or not content_provenance:
        raise ValueError("Content provenance must map source assets to SHA-256 digests")
    provenance = {}
    for asset, digest in content_provenance.items():
        if (
            not _text(asset)
            or not isinstance(digest, str)
            or len(digest) != 64
            or set(digest.lower()) - set("0123456789abcdef")
        ):
            raise ValueError("Content provenance requires asset names and 64-digit SHA-256 digests")
        provenance[asset] = digest.lower()
    scope = {
        "taxon_id": taxon_id,
        "assembly_accession": assembly_accession,
        "annotation_release": annotation_release,
        "content_provenance": dict(sorted(provenance.items())),
    }
    reference_id = hashlib.sha256(
        json.dumps(scope, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    ).hexdigest()
    groups = defaultdict(list)
    sequences = {}
    seen = {}
    for occurrence in occurrences:
        if not isinstance(occurrence, Mapping):
            raise ValueError("Occurrence must be a mapping")
        row = dict(occurrence)
        oid = row.get("occurrence_id")

        def reject(reason, oid=oid):
            raise ValueError(f"Occurrence {oid!r}: {reason}")

        if not _json_value(row):
            reject("source fields must contain finite JSON-native values")
        required = ("occurrence_id", "gene_id", "transcript_id", "protein_id", "source")
        if any(not _text(row.get(key)) for key in required):
            reject("incomplete occurrence identity or source provenance")
        if type(row.get("taxon_id")) is not int or any(
            row.get(key) != scope[key]
            for key in ("taxon_id", "assembly_accession", "annotation_release")
        ):
            reject("conflicting species/reference namespace")
        if "content_provenance" in row:
            row_content = row["content_provenance"]
            if (
                not isinstance(row_content, dict)
                or any(not isinstance(digest, str) for digest in row_content.values())
                or {asset: digest.lower() for asset, digest in row_content.items()} != provenance
            ):
                reject("conflicting reference content provenance")
        if row.get("complete") is not True:
            reject("incomplete translation (complete=True required)")
        if not _text(row.get("sequence")):
            reject("missing or non-string protein sequence")
        sequence = "".join(row["sequence"].split()).upper().removesuffix("*")
        invalid = sorted(set(sequence) - _AMINO_ACIDS)
        if not sequence or invalid:
            reject(f"empty or unresolved amino-acid sequence; invalid residues={invalid!r}")
        if oid in seen:
            # Python equality conflates True, 1 and 1.0; provenance must retain
            # JSON types and must not depend on which duplicate arrived first.
            if json.dumps(seen[oid], sort_keys=True) != json.dumps(row, sort_keys=True):
                reject("contradictory reuse of occurrence identity")
            continue
        row = deepcopy(row)
        seen[oid] = row
        sid = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        sequences[sid] = sequence
        groups[sid].append(row)
    return [
        {
            "group_id": f"{reference_id}:{sid}",
            "reference_id": reference_id,
            "protein_sequence_id": sid,
            "sequence": sequences[sid],
            **deepcopy(scope),
            "occurrences": sorted(rows, key=lambda row: row["occurrence_id"]),
        }
        for sid, rows in sorted(groups.items())
    ]

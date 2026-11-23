# Licensed under the Apache License, Version 2.0.
"""Versioned reference identities shared by OpenVax consumers (standard library only)."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from dataclasses import dataclass, fields
from types import MappingProxyType

REFERENCE_IDENTITY_SCHEMA = "openvax.reference"
REFERENCE_IDENTITY_VERSION = 1


def _text(value):
    return type(value) is str and bool(value.strip())


def _json_value(value):
    """Require finite JSON-native, acyclic source records without coercing types."""

    def visit(item, ancestors):
        if item is None or type(item) in (str, int, bool):
            return True
        if type(item) is float:
            return item == item and abs(item) != float("inf")
        if type(item) not in (list, dict) or id(item) in ancestors:
            return False
        parents = ancestors | {id(item)}
        if type(item) is list:
            return all(visit(child, parents) for child in item)
        return all(type(key) is str and visit(child, parents) for key, child in item.items())

    return visit(value, set())


def _canonical_json(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
    )


def _identity(kind, version, payload):
    digest = hashlib.sha256(_canonical_json(payload).encode("ascii")).hexdigest()
    return f"openvax:{kind}:v{version}:sha256:{digest}"


def _label(value, name):
    # Preserve case and spelling; do not invent aliases or normalize namespaces.
    if (
        not _text(value)
        or value != value.strip()
        or any(ord(c) < 32 or ord(c) == 127 for c in value)
    ):
        raise ValueError(
            f"{name} must be nonempty text without edge whitespace or control characters"
        )
    return value


def _string_map(value, name, *, digests=False):
    if not isinstance(value, Mapping) or not value:
        raise ValueError(f"{name} must be a nonempty mapping")
    result = {}
    for key, item in value.items():
        _label(key, f"{name} key")
        _label(item, f"{name} value")
        if digests:
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", key):
                raise ValueError(f"{name} keys must be logical asset roles, not paths or URLs")
            if len(item) != 64 or set(item.lower()) - set("0123456789abcdef"):
                raise ValueError(f"{name} requires 64-digit SHA-256 digests")
            item = item.lower()
        result[key] = item
    return dict(sorted(result.items()))


@dataclass(frozen=True)
class ReferenceIdentity:
    """One explicitly declared species/assembly/annotation/content namespace.

    All fields affect identity. Asset keys are stable logical roles, never local
    paths or retrieval URLs; digest values identify exact source bytes. The source
    version names the upstream snapshot, not the consuming package or pipeline.
    Identifier namespaces map entity kinds (e.g. gene/transcript/protein) to their
    naming authorities. Extra provenance belongs alongside this closed schema.

    Checksums are validated structurally; callers verify the actual source bytes.
    Equal identity requires identical declared scope and assets. Different IDs do
    not prove biological difference, and matching symbols/orthology never imply
    equality. There is no default taxon, annotation provider or reference.
    """

    taxon_id: int
    assembly_accession: str
    annotation_source: str
    annotation_release: str
    source_version: str
    identifier_namespaces: Mapping[str, str]
    content_provenance: Mapping[str, str]

    def __post_init__(self):
        if type(self.taxon_id) is not int or self.taxon_id <= 0:
            raise ValueError("Positive integer taxonomy identifier required")
        for name in (
            "assembly_accession",
            "annotation_source",
            "annotation_release",
            "source_version",
        ):
            _label(getattr(self, name), name)
        for name in ("identifier_namespaces", "content_provenance"):
            value = _string_map(getattr(self, name), name, digests=name == "content_provenance")
            object.__setattr__(self, name, MappingProxyType(value))

    def _payload(self):
        payload = {
            "schema": REFERENCE_IDENTITY_SCHEMA,
            "schema_version": REFERENCE_IDENTITY_VERSION,
        }
        for field in fields(self):
            value = getattr(self, field.name)
            payload[field.name] = dict(value) if isinstance(value, Mapping) else value
        return payload

    def __hash__(self):
        return hash(self.reference_id)

    @property
    def reference_id(self):
        return _identity("reference", REFERENCE_IDENTITY_VERSION, self._payload())

    def as_dict(self):
        """Return an independent JSON-native record, including its verifiable ID."""
        return {**self._payload(), "reference_id": self.reference_id}

    @classmethod
    def from_dict(cls, record):
        """Read v1 only; reject unknown fields/versions and mismatched stored IDs."""
        expected = {field.name for field in fields(cls)} | {
            "schema",
            "schema_version",
            "reference_id",
        }
        if not isinstance(record, Mapping) or set(record) != expected:
            raise ValueError("Reference identity requires the exact versioned schema fields")
        if (
            record["schema"] != REFERENCE_IDENTITY_SCHEMA
            or type(record["schema_version"]) is not int
            or record["schema_version"] != REFERENCE_IDENTITY_VERSION
        ):
            raise ValueError("Unsupported reference identity schema/version")
        result = cls(**{field.name: record[field.name] for field in fields(cls)})
        if record["reference_id"] != result.reference_id:
            raise ValueError("Reference identity digest does not match its contents")
        return result

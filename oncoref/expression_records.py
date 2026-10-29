# Licensed under the Apache License, Version 2.0.
"""Streaming source records accompanying a native expression matrix.

Records retain the producer's schema and identifiers. They are provenance, not
an instruction to allocate group/gene measurements to candidate proteins.
"""

from __future__ import annotations

import copy
import gzip
import json

RECORD_KINDS = (
    "samples",
    "assays",
    "measurements",
    "reference_mappings",
    "mappings",
    "relations",
    "pairings",
)


def json_copy(value):
    return json.loads(json.dumps(value, allow_nan=False))


class ExpressionRecords:
    """Repeatable iterators over verified, immutable source records."""

    def __init__(self, provenance, reader):
        self._provenance = json_copy(provenance)
        self._reader = reader

    @property
    def provenance(self):
        return copy.deepcopy(self._provenance)

    def records(self, kind):
        if kind not in RECORD_KINDS:
            raise ValueError(f"Unknown expression record kind {kind!r}; choose {RECORD_KINDS}")
        for record in self._reader(kind):
            if not isinstance(record, dict):
                raise ValueError("Native expression records must be JSON objects")
            yield copy.deepcopy(record)

    def write(self, directory):
        names = []
        for kind in RECORD_KINDS:
            name = f"{kind}.jsonl.gz"
            names.append(name)
            # No filename or clock in the gzip header: repeated exports are identical.
            with (
                (directory / name).open("wb") as raw,
                gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as stream,
            ):
                for record in self.records(kind):
                    stream.write(
                        (
                            json.dumps(
                                record, sort_keys=True, separators=(",", ":"), allow_nan=False
                            )
                            + "\n"
                        ).encode("utf-8")
                    )
        return names

    @classmethod
    def from_directory(cls, directory, provenance):
        def read(kind):
            with gzip.open(directory / f"{kind}.jsonl.gz", "rt", encoding="utf-8") as stream:
                for line in stream:
                    yield json.loads(line)

        return cls(provenance, read)

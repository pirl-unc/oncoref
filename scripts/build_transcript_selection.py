#!/usr/bin/env python3
"""Build the shared transcript-selection registry from pinned offline inputs."""

import gzip
import hashlib
import json
from pathlib import Path

from oncoref import transcript_selection, translation_qualification

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "oncoref/data"


def build():
    table = translation_qualification.records()
    selected = transcript_selection.select(table, transcript_selection.evidence())
    columns = [c for c in selected if selected[c].map(lambda v: isinstance(v, list)).any()]
    encoded = selected.copy()
    for col in columns:
        encoded[col] = encoded[col].map(lambda v: json.dumps(v, separators=(",", ":")))
    out = DATA / "transcript-selection.csv.gz"
    out.write_bytes(
        gzip.compress(encoded.to_csv(index=False, lineterminator="\n").encode(), mtime=0)
    )
    meta = {
        "schema": "oncoref.transcript-selection",
        "schema_version": 1,
        "policy_id": transcript_selection.POLICY_ID,
        "records_sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
        "json_columns": columns,
        "dependencies": {
            p: hashlib.sha256((DATA / p).read_bytes()).hexdigest()
            for p in [
                "translation-qualification.csv.gz",
                "translation-qualification.json",
                "transcript-evidence.json",
            ]
        },
        "status_counts": selected.selection_status.value_counts().to_dict(),
        "basis_counts": selected.selection_basis.value_counts().to_dict(),
    }
    (DATA / "transcript-selection.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2), flush=True)


if __name__ == "__main__":
    build()

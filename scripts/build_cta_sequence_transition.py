#!/usr/bin/env python3
"""Compare pinned historical forbidden sequences with current full-record roles."""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from build_cta_contract_report import write_table  # noqa: E402

from oncoref.sequence_bundle import ROLE_ORDER, digest, load  # noqa: E402


def sequences(path):
    with gzip.open(path, "rt") as handle:
        pieces = []
        for line in handle:
            if line.startswith(">"):
                if pieces:
                    yield "".join(pieces)
                pieces = []
            else:
                pieces.append(line.strip())
        if pieces:
            yield "".join(pieces)


def build(historical, current):
    policy = json.loads((ROOT / "oncoref/data/cta-policy.json").read_text())
    pins = policy["sequence_artifacts"]
    old = set()
    for name in ["unsafe.fasta.gz", "additional-forbidden.fasta.gz", "reference.json"]:
        p = historical / name
        if digest(p.read_bytes()) != pins["sha256"][name]:
            raise ValueError("Historical sequence asset mismatch: " + name)
        if name.endswith(".fasta.gz"):
            old.update(digest(seq.encode()) for seq in sequences(p))
    bundle = load(current)
    rows = {}
    for record in bundle.records():
        key = record["sequence_sha256"]
        row = rows.setdefault(
            key,
            {
                "sequence_sha256": key,
                "historical_forbidden": key in old,
                "current_shared_sequence": "allow",
                "current_records": 0,
                "current_gene_ids": set(),
            },
        )
        row["current_shared_sequence"] = max(
            row["current_shared_sequence"], record["shared_sequence"], key=ROLE_ORDER.get
        )
        row["current_records"] += 1
        row["current_gene_ids"].update(o["gene_id"] for o in record["occurrences"])
    for key in old - set(rows):
        rows[key] = {
            "sequence_sha256": key,
            "historical_forbidden": True,
            "current_shared_sequence": "absent_requires_review",
            "current_records": 0,
            "current_gene_ids": set(),
        }
    for row in rows.values():
        row["current_gene_ids"] = sorted(row["current_gene_ids"])
    frame = pd.DataFrame(rows.values()).sort_values("sequence_sha256")
    output = ROOT / "docs/audits/sequence-bundle-20261104"
    write_table(frame, output / "historical-current-sequence-roles.csv.gz")
    counts = Counter(
        (r["historical_forbidden"], r["current_shared_sequence"]) for r in rows.values()
    )
    receipt = {
        "historical_policy_id": policy.get("policy_id", policy.get("policy")),
        "historical_assets": pins,
        "historical_reference": json.loads((historical / "reference.json").read_text()),
        "current_bundle_id": bundle.identity,
        "comparison_unit": "SHA-256 of full original amino-acid strings; no partial-segment or approximate equivalence",
        "counts": [
            {"historical_forbidden": old_flag, "current_role": role, "sequences": count}
            for (old_flag, role), count in sorted(counts.items())
        ],
        "notes": [
            "Historical Unsafe plus additional-forbidden remains a distinct named product.",
            "Absent sequences require coverage review, never implicit allowance.",
            "Current record roles are combined with exclude > review > allow across all occurrences.",
            "This comparison does not validate an approximate peptide-matching rule.",
        ],
        "output_sha256": digest((output / "historical-current-sequence-roles.csv.gz").read_bytes()),
        "generator_sha256": digest(Path(__file__).read_bytes()),
    }
    (output / "historical-transition.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt["counts"]))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("historical", type=Path)
    p.add_argument("current", type=Path)
    a = p.parse_args()
    build(a.historical, a.current)

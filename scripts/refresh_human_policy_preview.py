#!/usr/bin/env python3
"""Refresh human-integrated-policy-preview.{csv,json} and curation-build-receipt.json."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import curation, human_screening_inputs  # noqa: E402


def main() -> None:
    data_dir = ROOT / "oncoref/data"
    changed, summary = human_screening_inputs.compute_next_policy_transitions()
    csv_path = data_dir / "human-integrated-policy-preview.csv"
    changed.to_csv(csv_path, index=False)
    summary["sha256"] = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    (data_dir / "human-integrated-policy-preview.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    curation.write_receipt(ROOT)
    verified = curation.check(ROOT)
    print(
        json.dumps(
            {
                "changed_proteins": summary["changed_proteins"],
                "transition_counts": summary["transition_counts"],
                "curation_check": verified,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

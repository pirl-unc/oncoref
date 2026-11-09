"""Current breadth exports bind the ontology and preserve missing coverage."""

import hashlib
import json
import zipfile
from pathlib import Path

import pandas as pd

from oncoref import normal_expression_breadth

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/audits/normal-breadth-20261104"


def test_current_breadth_receipt_and_policy():
    receipt = json.loads((REPORT / "manifest.json").read_text())
    assert receipt["eligibility_changes"] == 0
    assert receipt["policy"] == normal_expression_breadth.policy()
    for kind, base in [("inputs", ROOT), ("outputs", REPORT)]:
        for name, expected in receipt[kind].items():
            if kind == "inputs" and name in receipt.get("input_snapshots", []):
                saved = receipt["input_snapshot_archive"]
                path = REPORT / saved["path"]
                assert hashlib.sha256(path.read_bytes()).hexdigest() == saved["sha256"]
                with zipfile.ZipFile(path) as archive:
                    content = archive.read(name)
            else:
                content = (base / name).read_bytes()
            assert hashlib.sha256(content).hexdigest() == expected, name
    panels = pd.read_csv(REPORT / "panels.csv")
    assert (
        panels.groupby(["scenario", "brain_weight", "category"])
        .weight.sum()
        .between(1 - 1e-12, 1 + 1e-12)
        .all()
    )
    donors = pd.read_csv(REPORT / "donor-breadth.csv")
    unknown = donors[donors.summary.eq("q90")]
    assert unknown.lower.eq(0).all() and unknown.upper.eq(1).all()

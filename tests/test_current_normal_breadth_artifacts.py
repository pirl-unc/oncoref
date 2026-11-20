"""Current breadth exports bind the ontology and preserve missing coverage."""

import hashlib
import json
import zipfile
from pathlib import Path

import pandas as pd

from oncoref import data_file_renames, normal_expression_breadth, tissue_tiers
from tests import released_inputs

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/audits/normal-breadth-20261104"


def test_current_breadth_receipt_and_policy():
    receipt = json.loads((REPORT / "manifest.json").read_text())
    assert receipt["eligibility_changes"] == 0
    # The receipt names the registry as released; the live file is its renamed copy.
    recorded, current = dict(receipt["policy"]), normal_expression_breadth.policy()
    released = hashlib.sha256(
        data_file_renames.as_released("normal-tissue-hazards.yaml")
    ).hexdigest()
    assert recorded.pop("hazard_registry_id").endswith(released)
    assert current.pop("tier_registry_id") == tissue_tiers.policy_id()
    assert recorded.pop("policy_id") == "oncoref.normal-expression-breadth.2026-11-04.1"
    assert current.pop("policy_id") == normal_expression_breadth.POLICY_ID
    assert recorded == current
    for kind, base in [("inputs", ROOT), ("outputs", REPORT)]:
        for name, expected in receipt[kind].items():
            if kind == "inputs" and name in receipt.get("input_snapshots", []):
                saved = receipt["input_snapshot_archive"]
                path = REPORT / saved["path"]
                assert hashlib.sha256(path.read_bytes()).hexdigest() == saved["sha256"]
                with zipfile.ZipFile(path) as archive:
                    content = archive.read(name)
            else:
                content = (
                    released_inputs.read(name, expected)
                    if base == ROOT
                    else (base / name).read_bytes()
                )
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

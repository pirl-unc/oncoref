"""The published hazard review must match its registry and portable receipts."""

import hashlib
import json
from pathlib import Path

import pandas as pd

from oncoref import cta_rare_cell, tissue_hazards

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/audits/tissue-hazards-1.29.0"


def test_report_receipts_and_complete_assessment_exports():
    receipt = json.loads((REPORT / "manifest.json").read_text())
    for role, base in [("inputs", ROOT), ("outputs", REPORT)]:
        for name, digest in receipt[role].items():
            assert not Path(name).is_absolute()
            assert hashlib.sha256((base / name).read_bytes()).hexdigest() == digest, name
    doc = tissue_hazards.registry()
    assert receipt["registry_id"] == doc["registry_id"]
    exported = pd.read_csv(REPORT / "assessments.csv")
    assert set(exported.id) == set(tissue_hazards.assessments().id)
    assert exported.evidence_scope.notna().all()
    assert set(pd.read_csv(REPORT / "references.csv").reference_id) == set(doc["references"])
    assert len(pd.read_csv(REPORT / "canonical-anatomy.csv")) == len(
        tissue_hazards.anatomy_inventory()
    )


def test_rare_cell_provenance_names_the_same_policy():
    assert (
        cta_rare_cell.provenance()["rules"]["hazard_registry_id"]
        == tissue_hazards.registry()["registry_id"]
    )
    assert cta_rare_cell.provenance()["rules"]["protected_action"] == "candidate_review_hold"


def test_recorded_release_baseline_and_change_table_reconcile():
    baseline = pd.read_csv(REPORT / "baseline-classifications.csv.gz")
    assert baseline.category.value_counts().to_dict() == {"Unsafe": 18278, "Other": 689, "CTA": 637}
    current = pd.read_csv(
        ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz",
        low_memory=False,
    )
    joined = baseline.merge(
        current,
        on="protein_equivalence_group_id",
        suffixes=("_before", "_after"),
        validate="one_to_one",
    )
    changed = joined[joined.classification_before.ne(joined.classification_after)]
    exported = pd.read_csv(REPORT / "classification-changes.csv")
    assert set(changed.protein_equivalence_group_id) == set(exported.protein_equivalence_group_id)

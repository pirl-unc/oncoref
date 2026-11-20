"""The published tier review must match its registry and portable receipts."""

import hashlib
import io
import json
import zipfile
from pathlib import Path

import pandas as pd

from oncoref import cta_rare_cell, data_file_renames, tissue_tiers

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/audits/tissue-hazards-1.29.0"


def test_report_receipts_and_complete_assessment_exports():
    receipt = json.loads((REPORT / "manifest.json").read_text())
    saved = receipt["input_snapshot_archive"]
    archive_path = REPORT / saved["path"]
    assert hashlib.sha256(archive_path.read_bytes()).hexdigest() == saved["sha256"]
    for role, base in [("inputs", ROOT), ("outputs", REPORT)]:
        for name, digest in receipt[role].items():
            assert not Path(name).is_absolute()
            if role == "inputs":
                with zipfile.ZipFile(archive_path) as archive:
                    content = archive.read(name)
            else:
                content = (base / name).read_bytes()
            assert hashlib.sha256(content).hexdigest() == digest, name
    doc = tissue_tiers.registry()
    # The 1.29.0 receipt names the registry as released; the live file is its renamed copy.
    released = data_file_renames.as_released("normal-tissue-hazards.yaml")
    assert receipt["registry_id"].endswith(hashlib.sha256(released).hexdigest())
    exported = pd.read_csv(REPORT / "assessments.csv")
    assert set(exported.id) == set(tissue_tiers.registry_assessments().id)
    assert exported.evidence_scope.notna().all()
    assert set(pd.read_csv(REPORT / "references.csv").reference_id) == set(doc["references"])
    historical = set(pd.read_csv(REPORT / "canonical-anatomy.csv").concept_id)
    assert set(tissue_tiers.anatomy_inventory().concept_id) == historical | set(
        tissue_tiers.updates().id
    )


def test_rare_cell_provenance_names_the_same_policy():
    assert (
        cta_rare_cell.provenance()["rules"]["tier_registry_id"]
        == tissue_tiers.registry()["registry_id"]
    )
    assert cta_rare_cell.provenance()["rules"]["protected_action"] == "candidate_review_hold"


def test_recorded_release_baseline_and_change_table_reconcile():
    baseline = pd.read_csv(REPORT / "baseline-classifications.csv.gz")
    assert baseline.category.value_counts().to_dict() == {"Unsafe": 18278, "Other": 689, "CTA": 637}
    receipt = json.loads((REPORT / "manifest.json").read_text())
    with zipfile.ZipFile(REPORT / receipt["input_snapshot_archive"]["path"]) as archive:
        current = pd.read_csv(
            io.BytesIO(
                archive.read(
                    "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz"
                )
            ),
            compression="gzip",
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

"""The public explanation must remain bound to the exact released decisions."""

import hashlib
import json
from collections import Counter
from pathlib import Path

from oncoref import cta_selection, sequence_bundle
from oncoref.version import __version__

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/audits/cta-definition-20261104"


def test_public_counts_statuses_and_worked_cases_match_the_current_api():
    report = json.loads((REPORT / "report-data.json").read_text())
    roles = cta_selection.sequence_dispositions()
    assert report["release"] == __version__
    assert report["coverage"] == cta_selection.coverage()
    assert report["counts"] == roles.category.value_counts().to_dict()
    assert len(report["rows"]) == len(roles)
    assert Counter(row["shared_sequence"] for row in report["rows"]) == Counter(
        roles.shared_sequence
    )
    assert report["bundle"] == sequence_bundle.descriptor()
    current = roles.set_index("label")
    for case in report["worked"]:
        row = current.loc[case["label"]]
        assert case["classification"] == row.classification
        assert case["assessment_status"] == row.assessment_status
        assert case["sequence_role"] == row.shared_sequence
    assert {"MAGEC2", "LGALS13", "PLCZ1", "MAGEA12", "SERPINA2", "MROH5"} <= {
        row["label"] for row in report["worked"]
    }


def test_final_receipt_binds_source_audits_figures_and_public_documentation():
    receipt = json.loads((REPORT / "release-receipt.json").read_text())
    assert receipt["release"] == __version__
    assert receipt["bundle_id"] == sequence_bundle.descriptor()["bundle_id"]
    assert receipt["selection_policy_id"] == cta_selection.coverage()["policy_id"]
    for key, base in [("inputs", ROOT), ("outputs", REPORT), ("generated_documentation", ROOT)]:
        for name, expected in receipt[key].items():
            assert hashlib.sha256((base / name).read_bytes()).hexdigest() == expected, name
    assert {"report.html", "report-data.json", "cta-definition.pdf"} <= set(receipt["outputs"])
    assert len([name for name in receipt["outputs"] if name.endswith(".svg")]) == 4


def test_walkthrough_slide_status_counts_cover_every_classified_group():
    presentation = json.loads(
        (ROOT / "docs/audits/cta-source-walkthrough-20261028/presentation-data.json").read_text()
    )
    assert (
        presentation["assessment_statuses"]
        == cta_selection.results().assessment_status.value_counts().to_dict()
    )
    example = presentation["concentration_example"]
    current = cta_selection.results().set_index("label").loc[example["label"]]
    assert example["classification"] == current.classification
    profiles = cta_selection.profiles()
    matches = profiles[profiles.preliminary_tumor_restriction_profile]
    assert presentation["preliminary_tumor_restriction_profile_matches"] == len(matches)
    assert (
        presentation["tumor_restriction_statuses"]
        == matches.tumor_restriction_status.value_counts().to_dict()
    )


def test_inventory_audit_retains_baseline_and_discloses_restriction_conflicts():
    audit = ROOT / "docs/audits/cta-inventory-20261105"
    receipt = json.loads((audit / "manifest.json").read_text())
    for key, base in [("inputs", ROOT), ("outputs", audit)]:
        for name, expected in receipt[key].items():
            assert hashlib.sha256((base / name).read_bytes()).hexdigest() == expected, name
    summary = json.loads((audit / "summary.json").read_text())
    assert summary["input_records"] == len(cta_selection.results())
    assert summary["existing_classification_changes"] == 0
    assert sum(summary["restored_classifications"].values()) == summary["restored_records"]
    preliminary = summary["preliminary_expression_profile"]
    assert sum(preliminary["disposition_counts"].values()) == preliminary["matches"]
    assert preliminary["disposition_counts"]["excluded_by_normal_evidence"] == 6
    zero = summary["stored_zero_rna_screen"]
    assert zero["adult_and_reproductive_bulk_zero"] == (
        zero["also_developmental_measured"] + zero["developmental_measurement_missing"]
    )
    assert zero["developmental_measurement_missing"] > 0

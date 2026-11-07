"""The public explanation must remain bound to the exact released decisions."""

import hashlib
import json
from collections import Counter
from pathlib import Path

from oncoref import cta_selection, sequence_bundle

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/audits/cta-definition-20261104"


def test_public_counts_statuses_and_worked_cases_match_the_current_api():
    report = json.loads((REPORT / "report-data.json").read_text())
    roles = cta_selection.sequence_dispositions()
    assert report["release"] == cta_selection._metadata()["release"]
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
    assert receipt["release"] == cta_selection._metadata()["release"]
    assert receipt["bundle_id"] == sequence_bundle.descriptor()["bundle_id"]
    assert receipt["selection_policy_id"] == cta_selection.coverage()["policy_id"]
    for key, base in [("inputs", ROOT), ("outputs", REPORT), ("generated_documentation", ROOT)]:
        for name, expected in receipt[key].items():
            assert hashlib.sha256((base / name).read_bytes()).hexdigest() == expected, name
    assert {"report.html", "report-data.json", "cta-definition.pdf"} <= set(receipt["outputs"])
    assert len([name for name in receipt["outputs"] if name.endswith(".svg")]) == 5


def test_publication_funnel_preserves_all_sources_overlap_and_unannotated_candidates():
    import pandas as pd

    audit = cta_selection.audit()
    report = json.loads((REPORT / "report-data.json").read_text())
    table = pd.DataFrame(report["publication_funnel"]).set_index("source_id")
    papers = table[table.source_kind.eq("paper")]
    assert set(papers.index) == set(audit.sources.source_id)
    current = audit.sources.set_index("source_id")
    for source, row in papers.iterrows():
        assert row.Inventory == current.loc[source, "mapped_entities"]
        assert row["Unsafe sequence containment"] == current.loc[source, "CTA"]
        assert row["Final Unsafe"] == current.loc[source, "Unsafe"]
        assert row["Final Other"] == current.loc[source, "Other"]
    count_columns = table.select_dtypes("number").columns
    assert (
        table.loc["registered_union", count_columns] + table.loc["unannotated", count_columns]
    ).equals(table.loc["inventory", count_columns])
    assert table.loc["registered_union", "Inventory"] == audit.memberships.entity_id.nunique()
    assert (
        table.loc["inventory", "Unsafe sequence containment"]
        == (cta_selection.results().category == "CTA").sum()
    )
    # A source is allowed to overlap another source; the union is deduplicated.
    assert papers.Inventory.sum() > table.loc["registered_union", "Inventory"]
    assert table.loc["unannotated", "Unsafe sequence containment"] > 0


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

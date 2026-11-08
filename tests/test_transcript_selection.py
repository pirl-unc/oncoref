"""Reference preferences cannot bypass sequence quality or unresolved evidence."""

import json

import pandas as pd
import pytest

from oncoref import transcript_selection as selection
from oncoref.translation_qualification import annotate


def tx(name, seq, tags=(), **changes):
    result = {
        "gene_id": "G",
        "reference_id": "reference:1",
        "transcript_accession": name,
        "protein_id": name + "p",
        "sequence_sha256": seq,
        "length": 100,
        "contig": "1",
        "tags": list(tags),
        "target_qualified": True,
        "qualification_reasons": [],
        "transcript_support_level": "1",
    }
    result.update(changes)
    return result


def ledger(*findings):
    return {
        "schema": "oncoref.transcript-evidence",
        "schema_version": 1,
        "findings": list(findings),
    }


def finding(name="finding", **changes):
    result = {
        "id": name,
        "gene_id": "G",
        "reference_id": "reference:1",
        "scope": "gene",
        "disposition": "hold",
        "rationale": "Conflicting coding evidence",
        "evidence_kind": "translation_challenge",
        "reviewed_on": "2026-11-08",
        "limitations": "Tested cells only",
        "sources": [{"title": "Study", "url": "https://example.org/study"}],
    }
    result.update(changes)
    return result


def test_priority_is_independent_of_length_and_row_order():
    table = pd.DataFrame(
        [
            tx("mane.1", "aaa", ["MANE_Select"]),
            tx("canonical.1", "bbb", ["Ensembl_canonical"], length=200),
            tx("long.1", "ccc", ["gencode_basic"], length=300),
        ]
    )
    assert selection.select(table).selected_sequence_sha256.tolist() == ["aaa"]
    pd.testing.assert_frame_equal(selection.select(table), selection.select(table.iloc[::-1]))


def test_incomplete_mane_does_not_fall_back_to_canonical_or_patch():
    table = pd.DataFrame(
        [
            tx("mane.1", "aaa", ["MANE_Select"]),
            tx("canonical.1", "bbb", ["Ensembl_canonical"]),
            tx("patch.1", "ccc", ["MANE_Select"], contig="HSCHR1_FIX"),
        ]
    )
    table.loc[0, "target_qualified"] = False
    table.at[0, "qualification_reasons"] = ["incomplete_stop"]
    result = selection.select(table).iloc[0]
    assert result.selection_status == "review"
    assert result.selected_sequence_sha256 is None
    assert result.candidate_sequence_sha256s == ["aaa"]
    assert "incomplete_stop" in result.selection_reasons


def test_conflicting_designations_and_no_designation_require_review():
    table = pd.DataFrame(
        [tx("a.1", "aaa", ["Ensembl_canonical"]), tx("b.1", "bbb", ["Ensembl_canonical"])]
    )
    assert selection.select(table).selection_reasons.iloc[0] == ["conflicting_designated_sequences"]
    table["sequence_sha256"] = "aaa"
    assert selection.select(table).selected_transcripts.iloc[0] == ["a.1", "b.1"]
    table["tags"] = [[], []]
    assert selection.select(table).selection_reasons.iloc[0] == ["no_designated_representative"]


def test_support_does_not_clear_hold_and_resolution_is_explicit():
    table = pd.DataFrame([tx("a.1", "aaa", ["MANE_Select"])])
    hold = finding()
    support = finding("positive", disposition="support", evidence_kind="gene_rna")
    assert selection.select(table, ledger(hold, support)).selection_status.iloc[0] == "review"
    resolution = finding("resolution", disposition="resolved", resolves=["finding"])
    assert (
        selection.select(table, ledger(hold, support, resolution)).selection_status.iloc[0]
        == "selected"
    )
    resolution.update(scope="exact_sequence", sequence_sha256s=["aaa"])
    with pytest.raises(ValueError, match="gene-wide hold"):
        selection.select(table, ledger(hold, resolution))


def test_reviewed_preference_needs_exact_identity_and_cannot_bypass_quality():
    table = pd.DataFrame([tx("mane.1", "aaa", ["MANE_Select"]), tx("other.1", "bbb")])
    preference = finding(
        disposition="prefer",
        scope="exact_sequence",
        sequence_sha256s=["bbb"],
        evidence_kind="isoform_resolved_protein",
    )
    assert selection.select(table, ledger(preference)).selected_sequence_sha256.iloc[0] == "bbb"
    preference["evidence_kind"] = "gene_rna"
    with pytest.raises(ValueError, match="cannot select an isoform"):
        selection.select(table, ledger(preference))
    preference["evidence_kind"] = "reviewed_annotation_reconciliation"
    table.loc[1, "target_qualified"] = False
    assert selection.select(table, ledger(preference)).selection_status.iloc[0] == "review"
    preference["sequence_sha256s"] = ["absent"]
    with pytest.raises(ValueError, match="absent from locus"):
        selection.select(table, ledger(preference))
    preference.update(sequence_sha256s=["bbb"], reference_id="old-reference")
    with pytest.raises(ValueError, match="reference mismatch"):
        selection.select(table, ledger(preference))


@pytest.mark.parametrize(
    "value,flag",
    [
        ("", "transcript_support_not_supplied"),
        ("NA", "transcript_support_not_assessed"),
        ("5", "limited_full_transcript_alignment_support"),
        ("1 (assigned to previous version 1)", "support_assigned_to_previous_version"),
    ],
)
def test_support_details_preserved(value, flag):
    row = tx("a.2", "aaa", ["MANE_Select"])
    row["transcript_support_level"] = value
    selected = selection.select(pd.DataFrame([row])).iloc[0]
    assert flag in selected.support_flags
    assert selected.candidate_support[0]["transcript_support_level"] == value
    assert selected.protein_production == "not_established_by_reference_selection"


def test_every_member_must_share_an_approved_default():
    a, b = tx("a.1", "aaa", ["MANE_Select"]), tx("b.1", "bbb", ["MANE_Select"])
    b["gene_id"] = "H"
    facts = pd.DataFrame({"member_gene_ids": [json.dumps(["G", "H"]), '["absent"]', '["G"]']})
    result = annotate(facts, pd.DataFrame([a, b]), ledger=ledger(finding()))
    assert result.translation_qualified.tolist() == [False, False, True]
    assert not result.representative_selection_approved.any()
    assert result.representative_sequence_sha256s.tolist() == ["[]", "[]", "[]"]
    assert json.loads(result.representative_candidate_sha256s.iloc[2]) == ["aaa"]


def test_pinned_ledger_holds_are_applied_without_claiming_sequence_failure():
    rows = selection.records().set_index("gene_id")
    for claim in selection.evidence()["findings"]:
        row = rows.loc[claim["gene_id"]]
        assert row.selection_status == "review"
        assert row.qualification_passed
        assert claim["id"] in row.unresolved_evidence_ids
        assert row.selected_sequence_sha256 == ""

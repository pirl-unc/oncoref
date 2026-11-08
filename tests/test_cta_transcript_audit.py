"""Transcript support cannot be replaced by length or an alternate-locus rescue."""

import gzip
import json
from pathlib import Path

import pandas as pd
import pytest
from audit_cta_transcripts import proposal, support_summary

from oncoref.translation_qualification import records


def transcript(tid, seq, tags=(), **changes):
    return {
        "transcript_accession": tid,
        "transcript_id": tid.split(".")[0],
        "sequence_sha256": seq,
        "contig": "1",
        "tags": list(tags),
        "target_qualified": True,
        **changes,
    }


def test_mane_beats_longer_complete_canonical_and_basic_alternatives():
    table = pd.DataFrame(
        [
            transcript("short.1", "aaa", ["MANE_Select"], length=100),
            transcript("canonical.1", "bbb", ["Ensembl_canonical"], length=120),
            transcript("long.1", "ccc", ["gencode_basic"], length=300),
        ]
    )
    assert proposal(table)["sequence_sha256"] == "aaa"
    assert proposal(table.iloc[::-1])["sequence_sha256"] == "aaa"


def test_no_implicit_rescue_when_designated_transcript_fails():
    table = pd.DataFrame(
        [
            transcript("mane.1", "aaa", ["MANE_Select"], target_qualified=False),
            transcript("canonical.1", "bbb", ["Ensembl_canonical"]),
            transcript("patch.1", "ccc", ["MANE_Select"], contig="HSCHR1_FIX"),
        ]
    )
    result = proposal(table)
    assert result["sequence_sha256"] is None
    assert result["status"] == "review: representative fails qualification"


def test_different_proteins_require_review_identical_proteins_keep_transcripts():
    table = pd.DataFrame(
        [
            transcript("one.1", "aaa", ["Ensembl_canonical"]),
            transcript("two.1", "bbb", ["Ensembl_canonical"]),
        ]
    )
    assert proposal(table)["status"] == "review: conflicting representative proteins"
    table["sequence_sha256"] = "aaa"
    result = proposal(table)
    assert result["sequence_sha256"] == "aaa"
    assert result["qualified_transcripts"] == ["one.1", "two.1"]


def test_unsupported_longest_is_not_an_automatic_default():
    table = pd.DataFrame([transcript("long.1", "aaa", ["gencode_basic"])])
    assert proposal(table)["sequence_sha256"] is None


def test_missing_na_and_previous_version_support_remain_distinct():
    table = pd.DataFrame([transcript("one.2", "aaa")])
    assert support_summary(table, {"one": {"transcript_support_level": ""}}) == "missing"
    assert support_summary(table, {"one": {"transcript_support_level": "NA"}}) == "NA"
    prior = {"one": {"transcript_support_level": "1 (assigned to previous version 1)"}}
    assert support_summary(table, prior) == "1"
    assert prior["one"]["transcript_support_level"].endswith("version 1)")


@pytest.fixture(scope="module")
def reference():
    return records()


def test_real_galp_and_psg11_defaults_are_unique(reference):
    for symbol, transcript_id in [("GALP", "ENST00000357330.7"), ("PSG11", "ENST00000320078.12")]:
        result = proposal(reference.loc[reference.symbol.eq(symbol)])
        assert result["qualified_transcripts"] == [transcript_id]


def test_cdy_identity_split_is_detectable_but_rbmy_stays_identical(reference):
    def selected(symbol):
        return proposal(reference.loc[reference.symbol.eq(symbol)])["sequence_sha256"]

    assert selected("CDY1") != selected("CDY1B")
    assert selected("RBMY1F") == selected("RBMY1J")


@pytest.mark.parametrize("symbol", ["SERPINA2", "MROH5"])
def test_existing_qualification_failures_are_not_promoted(reference, symbol):
    assert proposal(reference.loc[reference.symbol.eq(symbol)])["sequence_sha256"] is None


def test_historical_report_covers_134_ctas_and_retains_alternative_translations():
    root = Path(__file__).resolve().parents[1]
    audit = json.loads(
        gzip.decompress((root / "docs/audits/cta-transcripts-20261108/audit.json.gz").read_bytes())
    )
    # This is the frozen pre-migration audit, not the current selection registry.
    assert len({r["entity_id"] for r in audit["proteins"]}) == 602
    assert len(audit["transcripts"]) == 1452
    proteins = {r["protein"]: r for r in audit["proteins"]}
    assert sum(r["default_change"] for r in proteins.values()) == 53
    assert proteins["CDY1/B"]["shared_identity_splits"]
    assert len(proteins["PSG11"]["current_sequence_sha256s"]) == 2
    assert len(proteins["PSG11"]["proposed_sequence_sha256s"]) == 1
    assert proteins["SPAM1"]["proposed_matches_surface_uniprot"] is True
    assert proteins["GDPD4"]["proposed_matches_surface_uniprot"] is False
    assert proteins["GDPD4"]["surface_uniprot_accession"] == "Q6W3E5"
    assert proteins["SLCO1B3-SLCO1B7"]["proposed_matches_surface_uniprot"] is False
    assert proteins["GALP"]["proposed_matches_surface_uniprot"] is None
    assert "Proposed transcript has weaker alignment support" in proteins["CIMIP3"]["review_flags"]
    alternate = [r for r in audit["transcripts"] if not r["proposed_default_transcript"]]
    assert any(not r["target_qualified"] for r in alternate)
    assert any("assigned to previous version" in r["transcript_support_level"] for r in alternate)

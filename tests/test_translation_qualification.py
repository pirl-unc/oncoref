"""Qualification cannot erase reference sequence or promote an alternate locus."""

import json

import pandas as pd

from oncoref.translation_qualification import annotate, qualify, records, valid_segments


def annotation(**changes):
    return {
        "contig": "1",
        "gene_biotype": "protein_coding",
        "transcript_biotype": "protein_coding",
        "cds_bases": 12,
        "start_codon_bases": 3,
        "stop_codon_bases": 3,
        "tags": [],
        **changes,
    }


def test_complete_cds_and_explicit_holds():
    assert qualify("MART*", annotation())["target_qualified"]
    for changes in [
        {"contig": "HSCHR14_7_CTG1"},
        {"stop_codon_bases": 0},
        {"start_codon_bases": 0},
        {"transcript_biotype": "nonsense_mediated_decay"},
        {"transcript_biotype": "protein_coding_LoF"},
        {"tags": ["cds_end_NF"]},
    ]:
        assert not qualify("MART", annotation(**changes))["target_qualified"]
    assert not qualify("MAXRT", annotation())["target_qualified"]


def test_segments_preserve_positions_without_crossing_gaps():
    segments = valid_segments("MARTXACD*EFGUHIK")
    assert [(s["start"], s["end"], s["sequence"]) for s in segments] == [
        (0, 4, "MART"),
        (5, 8, "ACD"),
        (9, 12, "EFG"),
        (13, 16, "HIK"),
    ]
    assert valid_segments("XXX") == ()


def test_primary_record_cannot_be_rescued_by_longer_alternate():
    table = records()
    ser = table[table.symbol.eq("SERPINA2")]
    assert set(ser.length) == {309, 421, 64}
    primary = ser[ser.contig.eq("14")].iloc[0]
    assert primary.transcript_accession.startswith("ENST00000553483.")
    assert not primary.target_qualified
    f = pd.DataFrame({"member_gene_ids": [json.dumps([primary.gene_id])]})
    result = annotate(f, ser)
    assert not result.translation_qualified.any()
    assert json.loads(result.representative_sequence_sha256s.iloc[0]) == [primary.sequence_sha256]
    mroh = table[table.symbol.eq("MROH5") & table.contig.eq("8")]
    assert not mroh.target_qualified.any()
    assert len(mroh) == 3


def test_missing_and_disagreeing_group_sequences_require_review():
    base = {"contig": "1", "length": 4, "target_qualified": True, "qualification_reasons": []}
    table = pd.DataFrame(
        [
            {**base, "gene_id": "A", "sequence_sha256": "aaa"},
            {**base, "gene_id": "B", "sequence_sha256": "bbb"},
        ]
    )
    f = pd.DataFrame({"member_gene_ids": [json.dumps(["A", "B"]), json.dumps(["unknown"])]})
    result = annotate(f, table)
    assert not result.translation_qualified.any()
    assert result.translation_qualification_reasons.tolist() == [
        "representative_sequence_mismatch",
        "missing_primary_translation",
    ]

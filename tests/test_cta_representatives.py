"""Consumers must receive approved transcripts, never recompute longest CDS."""

import pandas as pd
import pytest

from oncoref import cta_representatives


@pytest.fixture(scope="module")
def representatives():
    return cta_representatives.records()


def test_all_current_cta_loci_have_exact_verified_accessions(representatives):
    verified = pd.read_csv("docs/audits/cta-current-review-20261108/verified-cta-transcripts.csv")
    observed = representatives.set_index("transcript_id")
    assert set(observed.index) == set(verified.transcript)
    for row in verified.to_dict("records"):
        choice = observed.loc[row["transcript"]]
        assert choice.gene_id == row["gene_id"]
        assert choice.protein_id == row["protein_id"]
        assert choice.sequence_sha256 == row["sequence_sha256"]
        assert choice.length == row["length"]
    assert set(representatives.annotation_release) == {"115"}


def test_reported_longest_transcript_failures_have_approved_choices(representatives):
    observed = representatives.set_index("symbol")
    for symbol, length in {
        "IZUMO2": 221,
        "C16orf90": 182,
        "ARL13A": 256,
        "CATSPER4": 472,
        "GAGE1": 117,
        "BOLL": 283,
    }.items():
        assert observed.loc[symbol, "length"] == length
        assert observed.loc[symbol, "selection_basis"] == "MANE Select"
    assert "LLCFC1" not in observed.index  # Current normal-bulk exclusion, not a CTA.


def test_gene_accessor_accepts_exact_symbols_and_rejects_other_releases():
    symbol = cta_representatives.records("BOLL", annotation_release=115)
    by_id = cta_representatives.records("ENSG00000152430", annotation_release="115")
    pd.testing.assert_frame_equal(symbol, by_id)
    assert len(symbol) == 1
    assert cta_representatives.records("GDPD4").empty  # Published hold is still in force.
    assert cta_representatives.records("LLCFC1").empty  # No historical-CTA fallback.
    assert cta_representatives.records("no-such-gene").empty
    with pytest.raises(ValueError, match="requested release"):
        cta_representatives.records("BOLL", annotation_release=114)
    with pytest.raises(ValueError, match="nonempty"):
        cta_representatives.records("")

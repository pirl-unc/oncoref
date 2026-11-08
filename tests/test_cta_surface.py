"""Guard sequence attribution, inclusive sizes, and measured versus absent tumors."""

import pytest

from oncoref import cta_literature
from oncoref._cta_surface import dimensions, interval


def test_inclusive_coordinates_and_unknown_endpoints():
    assert dimensions([[5, 5], [20, 30]]) == (12, 11)
    assert dimensions([]) == (None, None)
    assert (
        interval(
            {
                "location": {
                    "start": {"value": 5, "modifier": "EXACT"},
                    "end": {"value": 30, "modifier": "UNCERTAIN"},
                }
            }
        )
        is None
    )


def test_all_surface_candidates_have_sequence_scoped_sizes():
    rows = cta_literature.surface_domains().set_index("protein")
    assert len(rows) == 37
    assert rows.sequence_match.eq("exact_sequence").sum() == 35
    assert rows.extracellular_total_aa.notna().sum() == 34
    assert rows.loc["BTN1A1", "extracellular_total_aa"] == 216
    assert rows.loc["PKDREJ", "extracellular_total_aa"] == 1487
    assert rows.loc["PKDREJ", "extracellular_longest_aa"] == 1165
    assert rows.loc["RTP2", "extracellular_longest_aa"] == 6
    assert rows.loc["OR10G3", "extracellular_total_aa"] == 97
    assert rows.loc["OR10G3", "extracellular_longest_aa"] == 37
    for protein in ["TEX38", "XKR3", "XKR5"]:
        assert rows.loc[protein, "extracellular_intervals"] == []
    assert cta_literature.surface_domains("not-a-gene").empty


def test_gpi_processing_and_readthrough_coordinates_are_not_silently_transferred():
    rows = cta_literature.surface_domains().set_index("protein")
    assert rows.loc["IZUMO1R", "extracellular_intervals"] == [[20, 228]]
    assert rows.loc["LYPD4", "extracellular_total_aa"] == 199
    assert rows.loc["SPAM1", "extracellular_total_aa"] == 455
    assert rows.loc["SPAM1", "sequence_match"] == "mature_chain_identical_anchor_tail_differs"
    assert rows.loc["ZP4", "extracellular_intervals"] == [[19, 505]]
    readthrough = rows.loc["SLCO1B3-SLCO1B7"]
    assert readthrough.accession == "F5H094"
    assert readthrough.selected_lengths == [748]
    assert readthrough.sequence_match == "different_sequence"
    assert readthrough.selected_extracellular_intervals == []


def test_prevalence_counts_are_native_and_cohorts_remain_separate():
    rows = cta_literature.surface_cancer_rna()
    assert set(rows.cohort) == {"TCGA", "validation"}
    assert rows.protein.nunique() == 37
    assert len(rows) == 1112
    assert rows[rows.cohort.eq("TCGA")].cancer.nunique() == 21
    assert (rows.expressed_samples_ptpm_ge_5 <= rows.expressed_samples_ptpm_ge_1).all()
    assert (rows.expressed_samples_ptpm_ge_1 <= rows.samples).all()
    assert (rows.samples + rows.missing_samples == rows.source_samples).all()
    stomach = (
        cta_literature.surface_cancer_rna("CT83", cohort="TCGA")
        .set_index("cancer_code")
        .loc["STAD"]
    )
    assert stomach.expressed_samples_ptpm_ge_5 == 174
    assert stomach.samples == 346
    zeros = cta_literature.surface_cancer_rna("OR4M2", cohort="TCGA")
    assert zeros.expressed_samples_ptpm_ge_1.sum() == 0
    assert zeros.samples.min() > 0
    assert cta_literature.surface_cancer_ihc("TEX38").empty
    ihc = cta_literature.surface_cancer_ihc()
    assert ihc.protein.nunique() == 14
    assert (ihc.high + ihc.medium + ihc.low + ihc.not_detected == ihc.total).all()
    with pytest.raises(ValueError, match="cohort"):
        cta_literature.surface_cancer_rna(cohort="pooled")


def test_concern_details_distinguish_coding_gaps_from_experimental_counterevidence():
    concerns = cta_literature.concerns().set_index("protein")
    assert concerns.evidence_details.str.len().min() > 100
    assert "FLAG" in concerns.loc["FRG2", "evidence_details"]
    assert "9.4-pg/mL" in concerns.loc["RLN1", "interpretation"]
    assert "unnamed locus" in concerns.loc["ENSG00000237378", "interpretation"]
    assert "founding SOX" in concerns.loc["SRY", "interpretation"]

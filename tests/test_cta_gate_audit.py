"""Explanations must reproduce the owner gate, including exact boundaries."""

import pytest

from oncoref.cta_curation_plots import RELIABILITY_THRESHOLD
from oncoref.cta_gate_audit import bradley_candidate_audit, cta_gate_audit
from oncoref.cta_tissues import HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS
from oncoref.load_dataset import get_data


def test_every_candidate_gate_is_reproduced():
    result = cta_gate_audit()
    assert result.hpa_gate_recorded.equals(result.hpa_gate_recomputed)
    assert not (result.default_panel & ~result.hpa_gate_recomputed).any()


def test_bradley_outcomes_and_independent_failure_reasons():
    rows = bradley_candidate_audit().set_index("Symbol")
    assert len(rows) == 10 and set(rows.loc[rows.default_panel].index) == {"PLAC1"}
    assert rows.loc["VGLL1", "protein_gate_pass"]
    assert rows.loc["VGLL1", "rna_deflated_reproductive_frac"] == 0.828
    assert rows.loc["VGLL1", "required_rna_fraction"] == 0.9
    assert set(rows.loc[~rows.protein_gate_pass].index) == {
        "IGF2BP3",
        "DEPDC1B",
        "SLC38A9",
        "MMP11",
    }
    assert not rows.loc["CAPN6", "has_protein_evidence"]
    assert rows.loc["CAPN6", "required_rna_fraction"] == 0.97
    assert (rows.loc[~rows.default_panel, "gate_reasons"].str.contains("RNA_fraction_below")).all()


@pytest.mark.parametrize("fraction,passed", [(0.97, True), (0.9699, False), (None, False)])
def test_missing_protein_boundary_and_missing_rna(fraction, passed):
    row = get_data("cancer-testis-antigens").iloc[:1].copy()
    row["biotype"] = "protein_coding"
    row["protein_reliability"] = "no data"
    row["protein_reproductive"] = "no data"
    row["rna_deflated_reproductive_frac"] = fraction
    row["passes_filters"] = passed
    result = cta_gate_audit(row).iloc[0]
    assert result.hpa_gate_recomputed == passed
    if fraction is None:
        assert "missing_RNA" in result.gate_reasons


def test_plot_thresholds_cannot_drift_from_the_filter():
    assert RELIABILITY_THRESHOLD["no data"] == HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS["Missing"]
    for tier in ["Enhanced", "Supported", "Approved", "Uncertain"]:
        assert RELIABILITY_THRESHOLD[tier] == HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS[tier]


def test_corrupt_stored_gate_is_not_silently_explained():
    row = get_data("cancer-testis-antigens").iloc[:1].copy()
    row["passes_filters"] = not row.iloc[0].passes_filters
    with pytest.raises(ValueError, match="does not reproduce"):
        cta_gate_audit(row)

"""Evidence matching and tumor recurrence must not silently change meaning."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scripts.cta_quantitative_threshold_report import bulk_limits, tumor_comparison


def test_sample_maximum_can_admit_rare_tumors_without_prevalence_or_unit_conversion():
    facts = pd.DataFrame(
        {
            "gene_id": ["rare", "common", "absent", "normal_positive"],
            "positive_context": True,
            "rna_coverage_complete": True,
            "ihc_detected": [False, False, False, True],
            "adult_bulk_exclusion": False,
            "adult_cell_exclusion": False,
            "immune_exclusion": False,
            "cornea_exclusion": False,
            "tumor_positive": [False, True, False, True],
            "selected": [False, True, False, False],
        }
    )
    tumor = pd.DataFrame(
        {
            "gene_id": ["rare", "common", "absent", "normal_positive"],
            "cancer_code": "LUAD",
            "cohort": "TCGA",
            "samples": [1000, 1000, 0, 1000],
            "max_ptpm": [30, 29.9, np.nan, 100],
        }
    )
    result = tumor_comparison(facts, tumor)
    assert result.maximum_tumor_selected.tolist() == [True, False, False, False]
    assert result.selected.tolist() == facts.selected.tolist()
    assert result.tumor_maximum_max_ptpm.iloc[0] == 30
    assert "maximum_tumor_selected" not in facts
    with pytest.raises(ValueError, match="Duplicate"):
        tumor_comparison(facts, pd.concat([tumor, tumor.iloc[:1]]))
    tumor["max_ptpm"] = [10, 9.99, np.nan, 100]
    assert tumor_comparison(facts, tumor, floor=10).maximum_tumor_selected.tolist() == [
        True,
        False,
        False,
        False,
    ]


def ihc(level="Not detected", reliability="Enhanced", tissue="heart muscle"):
    return pd.DataFrame(
        {
            "Gene": ["g"],
            "Tissue": [tissue],
            "Cell type": ["myocytes"],
            "Level": [level],
            "Reliability": [reliability],
        }
    )


def bulk_result(protein, *, positive=True, conflict=False):
    return bulk_limits(
        pd.DataFrame({"Gene": ["g"], "Tissue": ["heart muscle"], "nTPM": [5.0]}),
        protein,
        positive_controls={"g"} if positive else set(),
        conflicts={"g"} if conflict else set(),
    ).iloc[0]


def test_local_negative_grade_selects_bulk_ceiling_with_original_cap():
    approved = bulk_result(ihc(reliability="Approved"))
    enhanced = bulk_result(ihc())
    assert approved.effective_ceiling_ntpm == 4 and approved.local_exclusion
    assert enhanced.rna_ceiling_ntpm == 12
    assert enhanced.effective_ceiling_ntpm == 10 and not enhanced.local_exclusion


def test_linked_bulk_uses_full_shared_ceiling_and_retains_unresolved_hold():
    bulk = pd.DataFrame(
        {"Gene": ["g", "g"], "Tissue": ["heart muscle", "unknown"], "nTPM": [12.0, 0]}
    )
    r = bulk_limits(bulk, ihc(), positive_controls={"g"}, conflicts=set(), legacy_cap=None)
    heart = r[r.tissue.eq("heart muscle")].iloc[0]
    assert heart.effective_ceiling_ntpm == 12 and not heart.local_exclusion
    assert heart.reference_review_ceiling_ntpm == 24
    assert r[r.tissue.eq("unknown")].iloc[0].local_exclusion


@pytest.mark.parametrize(
    "kind", ["other_tissue", "partial", "uncertain", "positive", "no_control", "conflict"]
)
def test_unusable_local_ihc_never_earns_relaxation(kind):
    protein = ihc()
    if kind == "other_tissue":
        protein = ihc(tissue="liver")
    elif kind == "partial":
        # Another gene establishes that this routine tissue has a second
        # expected context. Its absence for g must not be called complete.
        other = ihc().assign(Gene="other", **{"Cell type": "endothelial cells"})
        protein = pd.concat([protein, other], ignore_index=True)
    elif kind == "uncertain":
        protein = ihc(reliability="Uncertain")
    elif kind == "positive":
        protein = ihc(level="Low")
    row = bulk_result(protein, positive=kind != "no_control", conflict=kind == "conflict")
    assert row.effective_ceiling_ntpm == 2 and row.local_exclusion


def test_genome_context_revision_preserves_coverage_and_independent_critical_signals():
    root = Path(__file__).resolve().parents[1] / "docs/audits/cta-tissue-context-20261027"
    genes = pd.read_csv(root / "gene-comparison.csv.gz")
    assert len(genes) == 20097 and genes.gene_id.is_unique
    complete = genes[genes.cell_records.eq(426)]
    assert complete.cell_neutral_records.eq(27).all()
    assert complete.cell_screened_records.eq(399).all()
    assert (
        (complete.cell_neutral_records + complete.cell_screened_records)
        .eq(complete.cell_records)
        .all()
    )
    assert genes.selected.sum() == 457 and genes.historical_uniform_selected.sum() == 447
    assert genes.ihc_conditioned_selected.sum() == 201
    assert genes.maximum_tumor_ihc_conditioned_selected.sum() == 167
    named = genes.set_index("symbol")
    for symbol in ["PRAME", "CTAG2"]:
        assert not named.loc[symbol, "protected30_selected"]
    for symbol in ["CTAG1A", "CTAG1B", "MAGEA4", "CT83"]:
        assert named.loc[symbol, "ihc_conditioned_selected"]
    cells = pd.read_csv(root / "named-cell-measurements.csv")
    ct83 = cells[cells.gene_id.eq(named.loc["CT83", "gene_id"]) & cells.tissue.eq("testis")]
    assert ct83.neutral_reproductive_immune.all() and ct83.value.max() == 8.1
    heart = cells[cells.gene_id.eq(named.loc["PRAME", "gene_id"]) & cells.tissue.eq("heart muscle")]
    assert not heart.neutral_reproductive_immune.any()
    assert heart.loc[heart.read_count.ge(20), "value"].max() == 12


def test_linked_audit_recovers_immune_candidates_without_hidden_uniform_caps():
    root = Path(__file__).resolve().parents[1] / "docs/audits/cta-linked-thresholds-20261027"
    genes = pd.read_csv(root / "gene-comparison.csv.gz")
    assert len(genes) == 20097 and genes.gene_id.is_unique
    assert genes.linked_recurrence_selected.sum() == 244
    assert genes.linked_maximum_tumor_selected.sum() == 337
    assert genes.linked_tumor30_selected.sum() == 197
    assert not (genes.linked_tumor30_selected & ~genes.linked_maximum_tumor_selected).any()
    named = genes.set_index("symbol")
    for symbol in ["SMC1B", "PAGE1", "PSMA8", "ENTHD1", "RDM1"]:
        assert named.loc[symbol, "linked_maximum_tumor_selected"]
    # This real-data regression catches accidental retention of the old 5 TPM cap.
    assert named.loc["PSMA8", "dice_third_max_tpm"] > 5
    assert named.loc["PSMA8", "immune_exclusion"]
    for symbol in ["PRAME", "CTAG2"]:
        assert not named.loc[symbol, "linked_maximum_tumor_selected"]
    for symbol in ["DKKL1", "INSL4", "DCAF4L1", "RBAK-RBAKDN"]:
        assert named.loc[symbol, "linked_maximum_tumor_legacy_cornea_selected"]
        assert not named.loc[symbol, "linked_maximum_tumor_selected"]

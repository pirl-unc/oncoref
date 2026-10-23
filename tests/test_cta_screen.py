"""Protect donor exclusion, sample provenance and the reproductive IHC gate."""

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from oncoref import cta, cta_screen

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "docs/audits/cta-curation-20261021-r2"
REPORT = ROOT / "docs/audits/cta-screen-20261023"
GID = "ENSG00000126890"


@pytest.fixture(scope="module")
def sources():
    candidates = pd.read_csv(BASE / "candidate-decisions.csv.gz", keep_default_na=False)
    audit = pd.read_csv(BASE / "tissue-evidence.csv.gz", keep_default_na=False)
    return candidates[candidates.Symbol.eq("CTAG2")], audit[audit.Symbol.eq("CTAG2")]


def test_variable_is_excluded_despite_negative_ihc_and_passing_mean(sources):
    candidates, audit = sources
    result, tissues, donors = cta_screen.screen_candidates(candidates, audit)
    assert result.iloc[0].Category == "Variable"
    assert not result.iloc[0].Selected
    assert tissues.loc[tissues.screening_tier.eq("vital"), "point_result"].eq("Pass").all()
    row = donors.iloc[0]
    assert row["Reported zero"] == row["Failing donors"] == 2
    assert row["RNA/IHC pairing"] == "Different donors"
    assert "paired RNA/protein" in result.iloc[0].Reconsideration


def test_selectivity_ceiling_applies_to_donors_below_absolute_ceiling(sources):
    candidates, audit = sources
    donors = cta_screen.reviewed_donor_rna().iloc[:2].copy()
    donors["rna_ntpm"] = [7.1, 0]
    result, _, evidence = cta_screen.screen_candidates(candidates, audit, donors=donors)
    assert evidence.iloc[0]["Failing donors"] == 1  # 7.1 < 12, but above 5.538...
    assert result.iloc[0].Category == "Variable"
    donors["rna_ntpm"] = [3.0, 0]
    result, _, evidence = cta_screen.screen_candidates(candidates, audit, donors=donors)
    assert evidence.iloc[0]["Failing donors"] == 0
    assert result.iloc[0].Category == "Retained"


def test_zero_and_high_cannot_be_pooled_across_cohorts(sources):
    candidates, audit = sources
    donors = cta_screen.reviewed_donor_rna().iloc[:2].copy()
    donors.loc[donors.rna_ntpm.eq(0), "cohort"] = "different cohort"
    result, _, evidence = cta_screen.screen_candidates(candidates, audit, donors=donors)
    assert not evidence["Variable pattern"].any()
    assert result.iloc[0].Category == "Excluded"


def test_donor_limits_follow_the_same_local_ihc_confidence(sources):
    candidates, audit = sources
    donors = cta_screen.reviewed_donor_rna().iloc[:2].copy()
    donors["rna_ntpm"] = [4.5, 0]
    result, _, evidence = cta_screen.screen_candidates(candidates, audit, donors=donors)
    assert result.iloc[0].Category == "Retained"
    assert evidence.iloc[0]["Absolute ceiling (nTPM)"] == 12
    approved = audit.copy()
    heart = approved.tissue.eq("heart muscle")
    approved.loc[heart, "ihc_evidence_status"] = "negative_approved"
    approved.loc[heart, "rna_ceiling_ntpm"] = "4"
    approved.loc[heart, "minimum_reproductive_selectivity"] = "0.70"
    result, _, evidence = cta_screen.screen_candidates(candidates, approved, donors=donors)
    assert result.iloc[0].Category == "Variable"
    assert evidence.iloc[0]["Absolute ceiling (nTPM)"] == 4
    assert evidence.iloc[0]["Failing donors"] == 1


def test_revised_audit_cannot_be_tightened_twice(sources):
    _, audit = sources
    with pytest.raises(ValueError, match="original HPA audit"):
        cta_screen.screen_tissues(cta_screen.screen_tissues(audit))


def test_missing_tissue_cannot_pass_a_partial_panel(sources):
    candidates, audit = sources
    partial = audit.loc[~audit.tissue.eq("liver")]
    result, _, _ = cta_screen.screen_candidates(candidates, partial, donors=pd.DataFrame())
    assert result.iloc[0].Category == "Unassessed"
    assert not result.iloc[0].Selected
    assert "missing normal RNA: liver" in result.iloc[0].Reason


def test_default_donor_registry_can_screen_a_different_gene():
    candidates = pd.read_csv(BASE / "candidate-decisions.csv.gz", keep_default_na=False)
    audit = pd.read_csv(BASE / "tissue-evidence.csv.gz", keep_default_na=False)
    result, _, evidence = cta_screen.screen_candidates(
        candidates.loc[candidates.Symbol.eq("XAGE1B")], audit.loc[audit.Symbol.eq("XAGE1B")]
    )
    assert result.iloc[0].Category == "Retained"
    assert result.iloc[0]["Donor review"] == "Not assessed"
    assert evidence.empty


def test_independent_tissue_inventory_matches_pinned_manifest():
    from cta_continuous_report import panels

    assert set(panels()[0]) == cta_screen.NORMAL_PANEL


@pytest.mark.parametrize("bad", [float("nan"), -1, float("inf")])
def test_invalid_donor_measurement_is_not_a_zero(sources, bad):
    _, audit = sources
    donors = cta_screen.reviewed_donor_rna()
    donors.loc[0, "rna_ntpm"] = bad
    with pytest.raises(ValueError, match="Invalid donor RNA"):
        cta_screen.screen_donors(audit, donors)


def test_donor_duplicates_and_incompatible_units_fail_closed(sources):
    _, audit = sources
    donors = cta_screen.reviewed_donor_rna()
    duplicate = pd.concat([donors, donors.iloc[:1]], ignore_index=True)
    with pytest.raises(ValueError, match="Duplicate donors"):
        cta_screen.screen_donors(audit, duplicate)
    donors["unit"] = "TPM"
    with pytest.raises(ValueError, match="do not substitute raw TPM"):
        cta_screen.screen_donors(audit, donors)


def test_no_donor_review_does_not_mean_donor_negative(sources):
    candidates, audit = sources
    result, _, _ = cta_screen.screen_candidates(candidates, audit, donors=pd.DataFrame())
    assert result.iloc[0]["Donor review"] == "Not assessed"
    assert result.iloc[0].Category == "Retained"  # counterfactual mean-only screen


def test_negative_ihc_requires_a_reproductive_positive_control(sources):
    _, audit = sources
    audit = audit.copy()
    audit.loc[audit.tissue.eq("testis"), "supporting_ihc"] = "[]"
    row = cta_screen.screen_tissues(audit).set_index("tissue").loc["heart muscle"]
    assert row.rna_ceiling_ntpm == 2
    assert row.minimum_reproductive_selectivity == 0.8
    assert not row.negative_ihc_relaxation_used
    assert row.point_result == "Fail"


def test_missing_ihc_cannot_create_variable_category(sources):
    candidates, audit = sources
    audit = audit.copy()
    audit.loc[audit.tissue.eq("heart muscle"), "ihc_evidence_status"] = "unavailable"
    result, _, evidence = cta_screen.screen_candidates(candidates, audit)
    assert result.iloc[0].Category == "Excluded"
    assert not evidence.iloc[0]["Variable pattern"]


def test_positive_ihc_and_external_holds_override_variable(sources):
    candidates, audit = sources
    held = candidates.copy()
    held["external_safety_hold"] = True
    result, _, _ = cta_screen.screen_candidates(held, audit)
    assert result.iloc[0].Category == "Excluded"
    assert result.iloc[0].Reconsideration == ""
    positive = audit.copy()
    cells = json.loads(positive.loc[positive.tissue.eq("heart muscle"), "supporting_ihc"].iloc[0])
    cells[0]["ihc_level"] = "Low"
    positive.loc[positive.tissue.eq("heart muscle"), "supporting_ihc"] = json.dumps(cells)
    positive.loc[positive.tissue.eq("heart muscle"), "ihc_evidence_status"] = "positive_ihc"
    result, _, _ = cta_screen.screen_candidates(candidates, positive)
    assert result.iloc[0].Category == "Excluded"
    assert "positive normal IHC" in result.iloc[0].Reason


def test_published_selection_and_benchmark():
    result = cta_screen.results().set_index("Gene")
    assert result.Selected.sum() == 740
    assert len(cta_screen.retained_gene_ids()) == 740
    assert cta_screen.variable_gene_ids() == {GID}
    assert not (cta_screen.retained_gene_ids() & cta_screen.variable_gene_ids())
    assert all(isinstance(name, str) and name for name in cta_screen.retained_gene_names())
    assert cta_screen.variable_gene_names() == {"CTAG2"}
    assert "CTAG2" not in cta_screen.retained_gene_names()
    assert "CTAG2" not in cta.cta_gene_names()
    assert result.loc["ZNF730", "Category"] == "Excluded"
    benchmark = {
        "AKAP4",
        "CT83",
        "CTAG1A",
        "CTAG1B",
        "CTCFL",
        "GAGE1",
        "MAGEA1",
        "MAGEA10",
        "MAGEA3",
        "MAGEA4",
        "MAGEC1",
        "MAGEC2",
        "PRAME",
        "ROPN1",
        "SAGE1",
        "SPANXB1",
        "SSX2",
        "SSX4",
        "SSX5",
        "XAGE1A",
        "XAGE1B",
    }
    assert benchmark <= cta_screen.retained_gene_names()
    assert result.loc[list(benchmark), "Donor review"].eq("Not assessed").all()


def test_report_hashes_and_published_dataset_agree():
    manifest = json.loads((REPORT / "run-manifest.json").read_text())
    assert manifest["donor_reviewed_genes"] == ["CTAG2"]
    for section, root in [("inputs", ROOT), ("implementation_sha256", ROOT), ("outputs", REPORT)]:
        for name, expected in manifest[section].items():
            assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
    dataset = ROOT / "oncoref/data/cta-screen-results.csv"
    assert hashlib.sha256(dataset.read_bytes()).hexdigest() == manifest["published_dataset_sha256"]
    assert dataset.read_bytes() == (REPORT / "gene-results.csv").read_bytes()

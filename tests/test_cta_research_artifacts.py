"""The shipped shortlist and audit must preserve evidence and coverage boundaries."""

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from oncoref import cta_policy, cta_research_policy

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "docs/audits/cta-research-20261026"
BREADTH = ROOT / "docs/audits/normal-breadth-20261026"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize("report", [RESEARCH, BREADTH])
def test_pinned_implementations_inputs_and_outputs(report):
    manifest = json.loads((report / "run-manifest.json").read_text())
    for name, expected in manifest["inputs"].items():
        assert digest(ROOT / name) == expected, name
    for name, expected in manifest["outputs"].items():
        assert digest(report / name) == expected, name
    for name, expected in manifest.get("retrieval_inputs", {}).items():
        assert digest(report / name) == expected, name


def test_shipped_selection_rebuilds_and_retains_every_expression_decision():
    facts = pd.read_csv(ROOT / "docs/audits/cta-policy-20261025/decision-facts.csv.gz")
    findings = cta_policy.findings()
    rebuilt = cta_research_policy.classify(facts, findings, cta_research_policy.clinical_evidence())
    shipped = cta_research_policy.results()
    pd.testing.assert_frame_equal(shipped, rebuilt)
    expression = cta_policy.classify(facts, findings)
    assert shipped["Expression category"].equals(expression.Category)
    assert shipped["Expression forbidden"].equals(expression.Forbidden)
    assert shipped["Expression reasons"].equals(expression["All reasons"])
    assert shipped["Expression reason"].equals(expression.Reason)
    added = shipped.loc[shipped.Target & ~expression.Target]
    assert set(added.Gene) == {"PRAME", "MAGEA4"}
    assert len(cta_research_policy.gene_ids()) == 215
    assert {"PRAME", "MAGEA4", "CTAG1A", "CTAG1B"} <= set(shipped.loc[shipped.Target, "Gene"])
    assert not (expression.Target & ~shipped.Target).any()
    assert digest(ROOT / "oncoref/data/cta-research-gene-sets.csv.gz") == digest(
        RESEARCH / "gene-results.csv.gz"
    )


def test_literature_completion_is_separate_from_retrieval_and_antibody_adjudication():
    summary = json.loads((RESEARCH / "summary.json").read_text())
    ledger = pd.read_csv(RESEARCH / "literature-review-status.csv", keep_default_na=False)
    assert len(ledger) == 2540 and ledger["Gene ID"].is_unique
    assert ledger["Remaining work"].str.len().gt(0).all()
    assert not summary["issue_638_complete"]
    assert summary["review_status_counts"]["Unassessed primary papers"] == 2492
    retrieval = pd.read_csv(RESEARCH / "expanded-searches.csv")
    assert retrieval["Complete pagination"].all()
    assert not retrieval["Paper review complete"].any()
    assert set(pd.read_csv(RESEARCH / "unmapped-identity-reviews.csv").iloc[:, 0]) == {
        "BAGE",
        "GAGE2B",
    }


def test_breadth_keeps_nomination_gaps_controls_and_donor_uncertainty():
    scores = pd.read_csv(BREADTH / "breadth.csv.gz")
    fixed = scores.loc[scores.scenario.eq("fixed_baseline") & scores.brain_weight.eq(1)]
    assert fixed.gene_id.nunique() == 2543
    assert set(fixed.loc[fixed.Category.eq("Normal control"), "Gene"]) == {"ACTB", "GAPDH", "PPIA"}
    applicable = scores.loc[scores.scenario.eq("applicable_ihc") & scores.brain_weight.eq(1)]
    assert applicable.gene_id.nunique() == 2540
    absent = applicable.groupby("gene_id").observed_compartments.sum().eq(0)
    assert absent.sum() == 6
    entirely_missing = applicable.loc[applicable.gene_id.isin(absent.index[absent])]
    assert entirely_missing.lower.eq(0).all() and entirely_missing.upper.eq(1).all()
    donor = pd.read_csv(BREADTH / "donor-summaries.csv").set_index("summary")
    assert donor.donors.eq(4).all()
    assert pd.isna(donor.loc["q90", "value"])
    assert donor.loc["median", "bootstrap_lower"] <= donor.loc["median", "value"]
    assert donor.loc["median", "bootstrap_upper"] >= donor.loc["median", "value"]
    donor_scores = pd.read_csv(BREADTH / "donor-breadth.csv")
    unavailable = donor_scores.loc[donor_scores.summary.eq("q90")]
    assert unavailable.lower.eq(0).all() and unavailable.upper.eq(1).all()

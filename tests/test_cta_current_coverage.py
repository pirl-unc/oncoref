"""Current paper coverage must retain locus scope and adverse source roles."""

import json
from pathlib import Path

import pytest

from oncoref import cta_literature


def test_remapping_keeps_partial_coverage_and_does_not_borrow_labels():
    profiles = [
        {"entity_id": "new", "label": "A/B", "gene_ids": ["a", "b"], "category": "CTA"},
        {"entity_id": "split", "label": "C", "gene_ids": ["c"], "category": "CTA"},
        {"entity_id": "excluded", "label": "D", "gene_ids": ["d"], "category": "Other"},
    ]
    links = [
        {"source_id": "human", "gene_ids": ["a"], "eligible_for_coverage": True},
        {"source_id": "mouse", "gene_ids": ["b"], "eligible_for_coverage": False},
    ]
    result = cta_literature._reconcile_coverage(profiles, [], links, []).set_index("protein")
    assert set(result.index) == {"A/B", "C"}
    assert result.loc["A/B", "paper_coverage"] == "some_members"
    assert result.loc["A/B", "uncovered_gene_ids"] == ["b"]
    assert result.loc["A/B", "related_sources"] == ["mouse"]
    assert result.loc["C", "paper_coverage"] == "no_paper"


def test_normal_paper_closes_coverage_but_does_not_become_nomination():
    additions = cta_literature.followup_findings()
    gene = "ENSG00000183128"
    profiles = [{"entity_id": "calhm3", "label": "CALHM3", "gene_ids": [gene], "category": "CTA"}]
    row = cta_literature._reconcile_coverage(profiles, [], [], additions["memberships"]).iloc[0]
    assert row.paper_coverage == "all_members"
    assert row.nomination_sources == []
    assert row.exact_paper_sources == ["10.1371/journal.pone.0007682"]
    assert additions["sources"][0]["roles"] == ["normal_tissue"]
    assert "unchanged" in additions["classification_effect"]


def test_current_nomination_adds_new_locus_without_changing_frozen_memberships():
    profiles = [{"entity_id": "new", "label": "new", "gene_ids": ["new"], "category": "CTA"}]
    nominations = [{"source_id": "catalog", "gene_ids": ["new"]}]
    row = cta_literature._reconcile_coverage(profiles, nominations, [], []).iloc[0]
    assert row.paper_coverage == "all_members"
    assert row.nomination_sources == row.exact_paper_sources == ["catalog"]


def test_gap_report_retains_hpa_signal_when_other_tumor_tpm_is_low():
    path = (
        Path(__file__).resolve().parents[1]
        / "docs/audits/cta-current-review-20261108/paper-gap-judgments.json"
    )
    rows = {r["protein"]: r for r in json.loads(path.read_text())}
    for gene, ptpm in {"CIMIP3": 55.3853, "DEFB113": 63.7402, "MSANTD5": 12.2279}.items():
        assert rows[gene]["tumor_max_tpm"] < 4
        assert rows[gene]["tumor_max_ptpm"] == pytest.approx(ptpm)
    assert rows["ENSG00000289809"]["tumor_max_ptpm"] is None

"""Curated observations keep identity, species and decision strength separate."""

import copy
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from oncoref import cta_selection, gene_knowledge
from oncoref import cta_selection_policy as policy
from oncoref.cta_funnel import EvidenceLayer


def test_other_review_covers_every_baseline_protein_and_preserves_evidence():
    root = Path(__file__).resolve().parents[1]
    report = root / "docs/audits/cta-other-literature-20261030"
    baseline = pd.read_csv(report / "baseline-proteins.csv.gz")
    decisions = pd.read_csv(report / "protein-decisions.csv")
    assert decisions.protein_id.is_unique
    assert set(decisions.protein_id) == set(baseline.protein_equivalence_group_id)
    assert len(decisions) == 600
    assert decisions.decision.notna().all() and decisions.review_scope.notna().all()
    assert decisions.shared_sequence.isin(["allow", "exclude", "review"]).all()
    high = pd.read_csv(report / "high-raw-normal-review.csv")
    assert set(high.protein_label) == {"OR2B6", "SCGB1C2", "OR6C2"}
    assert high.classification.eq("Unassessed").all()
    receipt = json.loads((report / "manifest.json").read_text())
    for section, folder in [("inputs", root), ("outputs", report)]:
        for name, expected in receipt[section].items():
            assert hashlib.sha256((folder / name).read_bytes()).hexdigest() == expected, name


def test_findings_preserve_native_units_and_do_not_transfer_species():
    row = gene_knowledge.findings("IFNA4").iloc[0]
    assert row.value == 549.7
    assert row.unit == "activated/unstimulated protein ratio"
    assert row.effect == "exclude" and row.state == "activated"
    assert gene_knowledge.findings("IFNA4", taxonomy_id=10090).empty
    assert gene_knowledge.findings("unknown").empty
    first = gene_knowledge.findings("BPY2")
    first.iloc[0].gene_ids.append("mutation")
    assert "mutation" not in gene_knowledge.findings("BPY2").iloc[0].gene_ids


def test_review_is_not_protein_exclusion_and_requires_provenance():
    row = gene_knowledge.findings("OR10A6").iloc[0].to_dict()
    assert row["effect"] == "review" and row["unit"] == "positive neurons"
    wrong = copy.deepcopy(row)
    wrong["effect"] = "exclude"
    with pytest.raises(ValueError, match="established adult human"):
        gene_knowledge.validate([wrong])
    wrong = copy.deepcopy(row)
    wrong["limitations"] = ""
    with pytest.raises(ValueError, match=r"limitations|context"):
        gene_knowledge.validate([wrong])
    wrong = copy.deepcopy(row)
    wrong["unit"] = None
    with pytest.raises(ValueError, match="native value and unit"):
        gene_knowledge.validate([wrong])


def test_shared_sequence_any_normal_member_but_all_reproductive_members():
    normal = gene_knowledge.findings("GFRAL").iloc[0].to_dict()
    repro = gene_knowledge.findings("BPY2").iloc[0].to_dict()
    genes = repro["gene_ids"]
    facts = pd.DataFrame(
        {
            "member_gene_ids": [
                json.dumps(genes),
                json.dumps([*genes, "ENSG00000999999"]),
                json.dumps([normal["gene_ids"][0], genes[0]]),
            ]
        }
    )
    out = gene_knowledge.annotate(facts, pd.DataFrame([normal, repro]))
    assert out.knowledge_reproductive_support.tolist() == [True, False, False]
    assert out.knowledge_normal_exclusion.tolist() == [False, False, True]
    assert not out.knowledge_normal_hold.any()
    assert json.loads(out.iloc[2].knowledge_finding_ids) == sorted([normal["id"], repro["id"]])
    objects = gene_knowledge.findings().astype(object)
    objects = objects.where(pd.notna(objects), None)
    gene_knowledge.validate(objects.to_dict("records"))
    with pytest.raises(ValueError, match="explicit member identity"):
        gene_knowledge.annotate(pd.DataFrame({"member_gene_ids": ["[]"]}))


def test_reproductive_protein_below_rna_floor_does_not_override_safety_hold():
    facts = pd.DataFrame(
        {
            "entity": ["protein", "uncertain", "normal"],
            "reproductive_ntpm": [0.1] * 3,
            "developmental_mean_tpm": [0.0] * 3,
            "knowledge_reproductive_support": [True] * 3,
            "required_coverage": [True] * 3,
            "unresolved_evidence": [False] * 3,
            "other_disposable": [False] * 3,
            "other_low_expression": [True] * 3,
        }
    )
    facts = policy.association_evidence(facts)
    false = pd.Series(False, index=facts.index)
    layers = [
        EvidenceLayer("Reproductive / embryonic association", ~facts.positive_context, false),
        EvidenceLayer("Cancer expression", false, false),
        EvidenceLayer(
            "Reviewed normal-cell knowledge",
            pd.Series([False, True, True]),
            pd.Series([False, False, True]),
        ),
    ]
    result = policy.classify(facts, layers, entity_column="entity").classifications
    assert result.category.iloc[0] == "CTA"
    assert pd.isna(result.category.iloc[1])
    assert result.assessment_status.iloc[1] == "normal_evidence_hold"
    assert result.category.iloc[2] == "Unsafe"


def test_sequence_dispositions_keep_unresolved_background_separate(monkeypatch):
    monkeypatch.setattr(
        cta_selection,
        "results",
        lambda **kwargs: pd.DataFrame(
            {
                "category": pd.Series(["CTA", "Unsafe", "Other", None], dtype="string"),
                "entity_id": ["a", "b", "c", "d"],
            }
        ),
    )
    result = cta_selection.sequence_dispositions()
    assert result.direct_target_candidate.tolist() == [True, False, False, False]
    assert result.shared_sequence.tolist() == ["allow", "exclude", "allow", "review"]

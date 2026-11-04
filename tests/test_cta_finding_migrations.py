"""Historical findings and clinical counterevidence cannot vanish during migration."""

import copy
import hashlib
import json
from importlib.resources import files

import pandas as pd
import pytest

from oncoref import cta_selection, cta_selection_audit, gene_knowledge, tcr_trials


def test_every_legacy_finding_preserves_original_and_explicit_disposition():
    root = files("oncoref").joinpath("data")
    receipt = json.loads(root.joinpath("cta-finding-migrations.json").read_text())
    assert (
        hashlib.sha256(root.joinpath("cta-policy-findings.csv").read_bytes()).hexdigest()
        == receipt["historical_sha256"]
    )
    old = pd.read_csv(root.joinpath("cta-policy-findings.csv"), keep_default_na=False)
    migrated = gene_knowledge.historical_dispositions()
    assert migrated.historical_finding_id.is_unique
    assert set(migrated.historical_finding_id) == set(old["Review ID"])
    by_id = old.set_index("Review ID", drop=False).to_dict("index")
    for row in migrated.itertuples():
        assert row.original == by_id[row.historical_finding_id]
        assert row.disposition in {"retained", "revised", "unresolved", "rejected", "superseded"}
        assert row.rationale
        assert row.replacement_finding_id in set(gene_knowledge.findings(row.gene_id).id)
        assert pd.isna(row.current_gate_id) or row.current_gate_id in set(
            cta_selection.gates().gate_id
        )


def test_named_conflicts_are_review_holds_and_rejected_protein_is_not_reinstated():
    for gene in ["MAGEC2", "LGALS13", "PLCZ1"]:
        rows = gene_knowledge.findings(gene)
        held = rows[rows.effect.eq("review")]
        assert not held.empty and held.interpretation_status.eq("unresolved").all()
        assert held.review_status.eq("accepted").all()
        assert not rows.effect.eq("exclude").any()
    rejected = gene_knowledge.findings("PLCZ1").query("interpretation_status == 'rejected'")
    assert len(rejected) == 1 and rejected.effect.eq("context").all()
    assert "CFTR" in rejected.iloc[0].decision_reason
    bad = copy.deepcopy(rejected.iloc[0].to_dict())
    bad["effect"] = "exclude"
    with pytest.raises(ValueError, match="Unresolved or rejected"):
        gene_knowledge.validate([bad])


def test_clinical_injury_keeps_competing_causes_and_product_scope():
    event = tcr_trials.events().set_index("event_id").loc["magea3_neurological_deaths"]
    assert event.patients == 2 and event.grade == 5
    assert pd.isna(event.recognized_protein)
    assert {r["protein"] for r in event.causal_findings} == {"MAGEA12", "EPS8L2"}
    assert all(r["interpretation_status"] == "unresolved" for r in event.causal_findings)
    assert not tcr_trials.products().set_index("product_id").loc["magea3_a2_neuro", "rna_credit"]
    assert tcr_trials.gene_exclusions().effect.eq("context").all()


def test_current_audit_retains_sources_even_when_they_do_not_determine_category():
    audit = cta_selection_audit.audit()
    rejected = audit.findings("PLCZ1").query("interpretation_status == 'rejected'")
    assert len(rejected) == 1
    result = cta_selection.sequence_dispositions().set_index("label")
    for gene in ["MAGEC2", "LGALS13", "PLCZ1"]:
        assert result.loc[gene, "category"] == "Other"
        assert result.loc[gene, "shared_sequence"] == "review"
    mage = audit.explain("MAGEA12")
    assert result.loc["MAGEA12", "category"] == "Unsafe"
    assert "Rare vital-tissue cell expression" in set(mage.layer)

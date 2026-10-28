"""Cross-species source accounting and faithful released-decision explanations."""

import gzip
import json
from pathlib import Path

import pandas as pd
import pytest

from oncoref import canine, cta_selection, human, mouse
from oncoref import cta_selection_audit as module

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "docs/audits/cta-selection-20261028"


@pytest.fixture(scope="module")
def audit():
    return human.cta_selection_audit()


def test_paper_unions_and_overlaps(audit):
    members = audit.memberships
    assert members.entity_id.nunique() == 2482
    assert len({gene for ids in members.gene_ids for gene in ids}) == 2538
    assert (
        members.loc[members.category.eq("CTA"), "entity_id"].nunique()
        <= human.cta_selection_coverage()["counts"]["CTA"]
    )
    sets = members.groupby("source_id").entity_id.agg(set)
    for row in audit.overlaps.itertuples():
        assert row.mapped_overlap == len(
            sets.get(row.source_a, set()) & sets.get(row.source_b, set())
        )
        assert row.selected_overlap <= row.mapped_overlap
    assert audit.sources.citation.str.match(r".+ et al\. \(\d{4}\)").all()


def test_prior_cta_context_expands_paper_annotations_without_admission(audit):
    """Prior-antigen recognition is provenance, never a positive-control bypass."""
    sources = audit.sources.set_index("source_id")
    garg = "10.1158/1078-0432.CCR-06-2340"
    boel = "10.1016/s1074-7613(95)80053-0"
    assert (sources.loc[garg, "mapped_entities"], sources.loc[garg, "CTA"]) == (5, 2)
    assert (sources.loc[boel, "mapped_entities"], sources.loc[boel, "CTA"]) == (2, 0)
    members = audit.memberships[audit.memberships.source_id.eq(garg)]
    assert set(members.label) == {"SPAG9", "SPA17", "SYCP1", "NY-ESO-1", "SSX4/B"}
    assert audit.memberships[audit.memberships.source_id.eq(boel)].category.isna().all()
    context = pd.read_csv(ROOT / "oncoref/data/cta-publication-context.csv")
    assert context.evidence_role.eq("previously_recognized_CTA").all()
    assert context.source_location.notna().all()


def test_every_released_hold_and_unsafe_layer_is_explained(audit):
    facts = pd.read_csv(RELEASE / "protein-group-classifications.csv.gz", low_memory=False)
    for column in facts.filter(regex="^hold:"):
        layer = column[5:]
        label = layer.replace(
            "Reproductive / embryonic association", "Reproductive or developmental expression"
        )
        for decision, mask in [
            ("unsafe", facts["unsafe:" + layer]),
            ("hold", facts[column] & ~facts["unsafe:" + layer]),
        ]:
            expected = set(facts.loc[mask, "protein_equivalence_group_id"])
            explained = audit.explanations[
                audit.explanations.layer.eq(label) & audit.explanations.decision.eq(decision)
            ]
            assert set(explained.entity_id) == expected, (layer, decision)
    end = audit.funnel.iloc[-1]
    expected = human.cta_selection_coverage()["counts"]
    assert (end.eligible, end.Unsafe, end.Neither, end.Unassigned) == tuple(
        expected[c] for c in ["CTA", "Unsafe", "Neither", "Unassigned"]
    )
    assert audit.funnel[["eligible", "Unsafe", "Neither", "Unassigned"]].sum(axis=1).eq(19863).all()
    assert not audit.explanations.reason.str.contains(r"\bnan\b", regex=True).any()


def test_specific_contexts_and_same_tissue_negative_ihc(audit):
    row = (
        audit.explain("CTAG2")
        .query("layer == 'HPA tissue clusters' and decision == 'observation'")
        .iloc[0]
    )
    assert (row.tissue, row.cell_type, row.value, row.unit) == (
        "bronchus",
        "basal respiratory cells",
        24.9,
        "nTPM",
    )
    assert row.candidate_ceiling == 32
    assert "reviewed_compartment" in row.ihc_coverage
    colon = audit.explain("PRAME").query("tissue == 'colon'").iloc[0]
    assert (colon.value, colon.candidate_ceiling, colon.selection_tier) == (20.3, 128, "manageable")
    assert "glandular cells" in colon.ihc_coverage
    assert audit.explain("MAGEA1").query("decision == 'hold'").reason.str.contains("4.4 nTPM").any()
    assert audit.explain("MAGEA1").category.isna().all()
    assert audit.explain("TBXT").iloc[0].category == "CTA"
    for name in ("TBXT", "CT83"):
        assert "external_review_flag" in set(audit.explain(name).query("decision == 'flag'").value)
    assert audit.explain("ENSG00000185686").entity_id.nunique() == 1
    with pytest.raises(KeyError):
        audit.explain("not-a-gene")


def test_unassessed_species_and_detached_tables(audit):
    dog = canine.cta_selection_audit()
    assert dog.status == "unassessed" and dog.policy_id is None
    assert sorted(dog.sources.mapped_entities) == [5, 10]
    assert dog.memberships.entity_id.nunique() == 14
    assert dog.funnel.eligible.isna().all()
    assert dog.overlaps.selected_overlap.isna().all()
    assert dog.overlaps.query("source_a != source_b").mapped_overlap.eq(1).all()
    unavailable = mouse.cta_selection_audit()
    assert unavailable.explanations.empty and unavailable.memberships.empty
    assert not unavailable.funnel.note.str.contains("Dog").any()
    for name in ("sources", "memberships", "overlaps", "criteria", "explanations", "funnel"):
        assert (
            list(getattr(dog, name))
            == list(getattr(audit, name))
            == list(getattr(unavailable, name))
        )
    first = cta_selection.audit(species="dog")
    first.explanations.loc[0, "reason"] = "modified"
    assert canine.cta_selection_audit().explanations.loc[0, "reason"] != "modified"


def test_checksum_and_species_validation(tmp_path, monkeypatch):
    source = ROOT / "oncoref/data"
    dest = tmp_path / "data"
    dest.mkdir()
    for name in ("cta-selection-audit.json", "cta-selection-audit-dog.json.gz"):
        (dest / name).write_bytes((source / name).read_bytes())
    monkeypatch.setattr(module, "files", lambda _: tmp_path)
    package = dest / "cta-selection-audit-dog.json.gz"
    package.write_bytes(package.read_bytes() + b"x")
    with pytest.raises(ValueError, match="checksum"):
        module.audit(species="dog")
    data = json.loads(gzip.decompress((source / package.name).read_bytes()))
    data["taxon_id"] = 9606
    package.write_bytes(gzip.compress(json.dumps(data).encode()))
    import hashlib

    receipt = json.loads((dest / "cta-selection-audit.json").read_text())
    receipt["files"][package.name] = hashlib.sha256(package.read_bytes()).hexdigest()
    (dest / "cta-selection-audit.json").write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="species mismatch"):
        module.audit(species="dog")


def test_neither_profiles_and_unassigned_counts_reconcile(audit):
    facts = pd.read_csv(RELEASE / "protein-group-classifications.csv.gz", low_memory=False)
    neither = facts[facts.category.eq("Neither")]
    assert (neither.neither_low_expression | neither.neither_diffuse_manageable).all()
    assert neither.required_coverage.all()
    assert not neither.unresolved_evidence.any()
    assert not neither.vital_protected_exclusion.any()
    assert facts.loc[facts.category.isna(), "neither_profile"].isna().all()
    assert canine.cta_selection().category.isna().all()
    assert audit.funnel.iloc[:-1].Neither.eq(0).all()


def test_current_replay_and_walkthrough_receipts_match_artifacts():
    import hashlib

    for directory in [RELEASE, ROOT / "docs/audits/cta-source-walkthrough-20261028"]:
        receipt = json.loads((directory / "manifest.json").read_text())
        for role, base in [("inputs", ROOT), ("outputs", directory)]:
            for name, expected in receipt[role].items():
                assert hashlib.sha256((base / name).read_bytes()).hexdigest() == expected, name
    summary = json.loads((RELEASE / "summary.json").read_text())
    presentation = json.loads(
        (ROOT / "docs/audits/cta-source-walkthrough-20261028/presentation-data.json").read_text()
    )
    assert presentation["summary"] == summary
    assert human.cta_selection_coverage()["counts"] == summary["categories"]

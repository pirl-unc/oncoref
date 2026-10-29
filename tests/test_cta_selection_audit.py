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
    assert members.entity_id.nunique() >= 2482
    assert len({gene for ids in members.gene_ids for gene in ids}) >= 2538
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


def test_primary_source_consolidation_preserves_the_complete_union(audit):
    sources = audit.sources
    sets = audit.memberships.groupby("source_id").entity_id.agg(set)
    primary = sources[sources.primary_source]
    assert set().union(*(sets[sid] for sid in primary.source_id)) == set(
        audit.memberships.entity_id
    )
    for row in sources[~sources.primary_source].itertuples():
        assert row.mapped_entities < 10
        assert set(row.covered_by) <= set(primary.source_id)
        assert sets[row.source_id] <= set().union(*(sets[sid] for sid in row.covered_by))
    egfl6 = audit.memberships[audit.memberships.label.eq("EGFL6")]
    assert {"10.4161/cbt.12.3.16260", "10.1371/journal.pone.0144155"} <= set(egfl6.source_id)
    intake = pd.read_csv(
        ROOT
        / "docs/audits/cta-source-walkthrough-20261028/publication-intake/mapped-source-rows.csv"
    )
    hirbo = intake[intake.source_tag.eq("Hirbo2015_placenta")]
    assert len(hirbo) == 86 and hirbo.identity_status.eq("matched").sum() == 82
    assert (
        "placenta"
        in sources.set_index("source_id").loc["10.1371/journal.pone.0144155", "evidence_scope"]
    )


def test_prior_cta_context_expands_paper_annotations_without_admission(audit):
    """Prior-antigen recognition is provenance, never a positive-control bypass."""
    sources = audit.sources.set_index("source_id")
    garg = "10.1158/1078-0432.CCR-06-2340"
    boel = "10.1016/s1074-7613(95)80053-0"
    assert (sources.loc[garg, "mapped_entities"], sources.loc[garg, "CTA"]) == (5, 3)
    assert (sources.loc[boel, "mapped_entities"], sources.loc[boel, "CTA"]) == (2, 1)
    members = audit.memberships[audit.memberships.source_id.eq(garg)]
    assert set(members.label) == {"SPAG9", "SPA17", "SYCP1", "NY-ESO-1", "SSX4/B"}
    boel_members = audit.memberships[audit.memberships.source_id.eq(boel)].set_index("label")
    assert boel_members.loc["MAGEA1", "category"] == "CTA"
    assert pd.isna(boel_members.loc["MAGEA3", "category"])
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


def test_independent_assessments_and_only_normal_safety_can_establish_unsafe(audit):
    axes = audit.explanations[audit.explanations.layer.str.startswith("Assessment: ")]
    assert axes.groupby("entity_id").size().eq(4).all()
    assert not audit.explanations.layer.str.contains("Adult restriction", regex=False).any()
    assert not axes.loc[~axes.layer.eq("Assessment: Safety"), "decision"].eq("unsafe").any()
    assert (
        axes.loc[
            axes.category.eq("CTA") & ~axes.layer.eq("Assessment: Reproductive concentration"),
            "decision",
        ]
        .eq("pass")
        .all()
    )
    for label in ["MAGEA4", "NY-ESO-1", "PRAME", "CTAG2", "XAGE1A/B"]:
        assert audit.explain(label).category.eq("CTA").all()
    for label in ["TGIF2-RAB5IF", "ATF7-NPFF"]:
        rows = audit.explain(label)
        assert rows.category.isna().all()
        assert (
            rows.loc[rows.layer.eq("Required RNA coverage"), "reason"]
            .str.contains("GTEx v11")
            .all()
        )


def test_required_gtex_core_cannot_be_replaced_with_other_bulk_sources(audit):
    """Coverage failures with other passing evidence remain unknown, not Unsafe."""
    coverage = pd.read_csv(RELEASE / "combined-rna-coverage.csv.gz").set_index("entity_id")
    selected = human.cta_selection(category="CTA")
    assert coverage.loc[selected.entity_id, "gtex_core_complete"].all()
    for label in ["TSSK2", "NEDD8-MDP1", "NPHP3-ACAD11"]:
        explanation = audit.explain(label)
        assert explanation.category.isna().all()
        reason = explanation[explanation.layer.eq("Required RNA coverage")]
        assert reason.reason.str.contains("GTEx v11").all()
        assert reason.decision.eq("hold").all()


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
    assert audit.explain("MAGEA1").category.eq("CTA").all()
    mage = (
        pd.read_csv(RELEASE / "protein-group-classifications.csv.gz", low_memory=False)
        .set_index("protein_group_label")
        .loc["MAGEA1"]
    )
    assert mage.developmental_route and mage.developmental_mean_tpm == pytest.approx(
        6.315480240244203
    )
    assert not mage.reproductive_route  # admission comes from measured developmental RNA
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
    assert (
        neither.neither_low_expression
        | neither.neither_diffuse_manageable
        | neither.neither_profile.eq("diffuse_reproductive")
    ).all()
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
    facts = pd.read_csv(RELEASE / "protein-group-classifications.csv.gz", low_memory=False)
    assert sum(row["entrants"] for row in presentation["routes"]) == facts.positive_context.sum()
    assert sum(row["CTA"] for row in presentation["routes"]) == facts.category.eq("CTA").sum()


def test_missing_tumor_measurements_are_not_incomplete_measured_groups():
    facts = pd.read_csv(RELEASE / "protein-group-classifications.csv.gz", low_memory=False)
    measured = (
        facts[["oncoref_maximum_tpm", "hpa_tumor_group_lower_bound_ptpm"]].notna().any(axis=1)
    )
    pd.testing.assert_series_equal(facts.tumor_measured, measured, check_names=False)
    assert not facts.loc[facts.member_count.gt(1), "tumor_profile_complete"].any()
    gaps = pd.read_csv(
        ROOT / "docs/audits/cta-source-walkthrough-20261028/unassigned-evidence-gaps.csv"
    )
    assert not (gaps.tumor_measurement_missing & gaps.tumor_group_profile_incomplete).any()


def test_developmental_alias_recovery_retains_conflicting_vectors():
    recovery = pd.read_csv(RELEASE / "developmental-identity-recoveries.csv")
    lineages = pd.read_csv(RELEASE / "developmental-lineages.csv.gz")
    measured = set(lineages.loc[lineages.source.eq("tyser_E-MTAB-9388"), "Gene ID"])
    recovered = set(recovery.loc[recovery.status.eq("registered_alias"), "gene_id"])
    conflicting = set(recovery.loc[recovery.status.eq("conflicting_vectors"), "gene_id"])
    assert recovered <= measured
    assert not conflicting & measured
    assert {
        "ENSG00000164871",
        "ENSG00000179172",
        "ENSG00000204694",
        "ENSG00000212725",
        "ENSG00000293560",
    } <= recovered


def test_alternative_sources_restore_measurements_without_filling_missing_zeros():
    coverage = pd.read_csv(RELEASE / "combined-rna-coverage.csv.gz").set_index("label")
    assert coverage.loc["POLR2A", "bulk_complete"]
    assert coverage.loc["POLR2A", "required_coverage"]
    assert coverage.loc["CIMIP3", "bulk_complete"]
    assert json.loads(coverage.loc["CIMIP3", "missing_bulk_tissues"]) == []
    assert coverage.required_bulk_tissues.eq(36).all()
    bulk = pd.read_csv(RELEASE / "alternative-evidence/normal-bulk.csv.gz")
    pituitary = bulk[
        bulk.gene_id.eq("ENSG00000181222")
        & bulk.source.eq("gtex_tissue_rna")
        & bulk.tissue.eq("pituitary gland")
    ].iloc[0]
    assert pituitary.value == pytest.approx(109.319) and pituitary.unit == "TPM"
    assert pituitary.source_context == "Pituitary"
    cimip3 = bulk[
        bulk.gene_id.eq("ENSG00000287363")
        & bulk.source.eq("gtex_tissue_rna")
        & bulk.tissue.eq("pituitary gland")
    ].iloc[0]
    assert cimip3.value == 0 and cimip3.version == "v11"
    assert bulk.loc[bulk.source.eq("hpa_brain"), "unit"].eq("nTPM").all()
    lineage = pd.read_csv(RELEASE / "developmental-lineages.csv.gz")
    petro = set(lineage.loc[lineage.source.eq("petropoulos_E-MTAB-3929"), "Gene ID"])
    for symbol in ["PRAME", "MAGEA4", "NY-ESO-1", "XAGE1A/B", "CT83"]:
        row = human.cta_selection().set_index("label").loc[symbol]
        assert set(row.gene_ids) <= petro
    assert not lineage.duplicated(["Gene ID", "source", "Lineage"]).any()
    assert (
        lineage.loc[lineage.source.eq("petropoulos_E-MTAB-3929"), "stage"]
        .eq("preimplantation, 3–7 days")
        .all()
    )


def test_normal_evidence_keeps_native_source_and_concrete_release(audit):
    gtex = audit.explanations[audit.explanations.source_id.eq("gtex_tissue_rna")]
    assert len(gtex) > 0
    assert gtex.source_version.eq("v11").all()
    assert gtex.source_context.notna().all()
    assert gtex.source_url.str.contains("adult-gtex/bulk-gex/v11/", regex=False).all()
    assert audit.explanations.label.notna().all()


def test_diffuseness_changes_identity_without_changing_safety(audit):
    before = pd.read_csv(
        ROOT / "docs/audits/cta-unassigned-20261029/before-diffuseness.csv.gz", low_memory=False
    ).set_index("protein_equivalence_group_id")
    after = pd.read_csv(
        RELEASE / "protein-group-classifications.csv.gz", low_memory=False
    ).set_index("protein_equivalence_group_id")
    assert set(before[before.category.eq("Unsafe")].index) == set(
        after[after.category.eq("Unsafe")].index
    )
    changed = after.category.fillna("Unassigned").ne(before.category.fillna("Unassigned"))
    assert after.loc[changed, "protein_group_label"].tolist() == ["MTCP1"]
    assert after.loc[changed, "category"].tolist() == ["Neither"]
    row = after.loc[changed].iloc[0]
    assert row.hpa_diffuseness == pytest.approx(0.24533644300531762)
    assert row.neither_profile == "diffuse_reproductive"
    explanation = audit.explain("MTCP1")
    assert explanation[explanation.layer.eq("Assessment: Safety")].decision.eq("pass").all()
    assert (
        explanation[explanation.layer.eq("Reproductive concentration")]
        .reason.str.contains("24.53%")
        .all()
    )
    selected = human.cta_selection(category="CTA").set_index("label")
    assert "Developmental concentration unassessed" in selected.loc["FAM72D", "evidence_flags"]
    assert "Adult concentration unassessed" in selected.loc["ZC3H11C", "evidence_flags"]

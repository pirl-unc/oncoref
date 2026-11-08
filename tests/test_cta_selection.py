"""One current-selection API preserves species, identity, units and unknowns."""

import json
from pathlib import Path

import pandas as pd
import pytest

from oncoref import canine, cta_selection, human, mouse


@pytest.mark.parametrize("module,taxon", [(human, 9606), (canine, 9615), (mouse, 10090)])
def test_common_selection_schema_species_and_counts(module, taxon):
    frame = module.cta_selection()
    assert tuple(frame.columns) == cta_selection.RESULT_COLUMNS
    assert frame.taxon_id.eq(taxon).all()
    assert frame.entity_id.is_unique
    assert str(frame.taxon_id.dtype) == "int64"
    assert str(frame.category.dtype) == "string"
    summary = module.cta_selection_coverage()
    assert sum(summary["counts"].values()) == len(frame)
    assert set(summary["evidence_categories"]) == set(cta_selection.EVIDENCE_CATEGORIES)
    pd.testing.assert_frame_equal(frame, cta_selection.results(species=taxon))


def test_human_current_selection_matches_the_complete_released_funnel():
    audit = Path(__file__).resolve().parents[1] / "docs/audits/cta-selection-20261028"
    expected = json.loads((audit / "summary.json").read_text())
    assert human.cta_selection_coverage()["counts"] == expected["categories"]
    selected = human.cta_selection(category="CTA")
    assert selected.gene_ids.map(len).sum() == expected["selected_member_loci"]
    assert {"PRAME", "MAGEA4", "NY-ESO-1", "TBXT", "XAGE1A/B", "CT83"} <= set(selected.label)
    assert "CTAG2" in set(selected.label)
    assert selected.identity_status.eq("resolved").all()


def test_restored_reproductive_protein_route_preserves_old_ctas_and_cancer_gate():
    root = Path(__file__).resolve().parents[1]
    old = pd.read_csv(root / "docs/audits/cta-neither-20261030/baseline-1.24.0.csv.gz")
    current = human.cta_selection().set_index("label")
    previous_genes = {
        g for value in old.loc[old.category.eq("CTA"), "member_gene_ids"] for g in json.loads(value)
    }
    kept = current[current.gene_ids.map(lambda ids: bool(set(ids) & previous_genes))]
    # Earlier CTAs can leave under a revised tissue policy. Each removal must
    # retain a concrete exclusion or hold rather than silently losing evidence.
    profiles = cta_selection.profiles().set_index("label")
    contained = kept.classification.eq("Other:Contained")
    facts = pd.read_csv(
        root / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz",
        low_memory=False,
    ).set_index("protein_group_label")
    removed = kept.index[~kept.category.eq("CTA") & ~contained]
    safety = [c for c in facts if c.startswith("unsafe:")]
    unsafe = removed[kept.loc[removed, "category"].eq("Unsafe")]
    assert facts.loc[unsafe, safety].any(axis=1).all()
    held = removed[kept.loc[removed, "classification"].eq("Other:Held")]
    assert (facts.loc[held, "rare_cell_hold"] | facts.loc[held, "knowledge_normal_hold"]).all()
    translation = removed[kept.loc[removed, "classification"].eq("Other:Translation")]
    assert (
        ~facts.loc[translation, "translation_qualified"]
        | ~facts.loc[translation, "representative_selection_approved"]
    ).all()
    assert set(removed) == set(unsafe) | set(held) | set(translation)
    overlap = profiles.unsafe_sequence_overlap_pct
    assert overlap.loc[kept.index[contained]].ge(95).all()
    assert contained.sum() <= 40
    assert current.loc[["SAGE1", "GAGE10", "CDY2A/B"], "category"].eq("CTA").all()
    assert current.loc["CT45A2/8/9", "category"] == "CTA"  # below the 95% containment floor
    assert current.loc[["CDY1", "CDY1B"], "category"].eq("CTA").all()
    assert current.loc[["USP17L3", "FAM90A9"], "category"].eq("Other").all()
    candidates = cta_selection.candidates()
    assert set(candidates.category) <= set(cta_selection.CANDIDATE_CATEGORIES)
    assert not candidates.entity_id.duplicated().any()


def test_native_nominations_never_inherit_human_admission_or_human_gene_ids():
    frame = canine.cta_selection()
    assert len(frame) == 14
    assert frame.category.isna().all()
    assert frame.assessment_status.eq("unassessed").all()
    assert frame.policy_id.isna().all()
    assert frame.entity_kind.eq("native_nomination").all()
    assert not any(gid.startswith("ENSG") for ids in frame.gene_ids for gid in ids)
    assert any(len(ids) > 1 for ids in frame.gene_ids)  # Ambiguous family assays remain visible.
    assert canine.cta_selection(category="CTA").empty
    assert len(mouse.cta_selection()) == 1
    assert mouse.cta_selection().category.isna().all()
    assert mouse.cta_selection_coverage()["selection_status"] == "unavailable"


def test_sources_distinguish_native_units_and_unintegrated_resources():
    dog = canine.cta_sources().drop_duplicates("source_id").set_index("source_id")
    assert dog.loc["barkbase", "native_unit"] == "TMM-CPM"
    assert dog.loc["pfizer-dog-protein", "native_unit"] == "iBAQ"
    assert dog.loc["doga", "status"] == "metadata_only"
    assert dog.loc["gse304066-osteosarcoma", "status"] == "blocked_provenance"
    assert set(mouse.cta_sources().status) == {"not_integrated", "native_observations"}
    assert not canine.cta_sources(category="normal_protein").source_id.eq("nemec2019").any()


def test_nested_identities_are_detached_and_missing_history_remains_unknown():
    first = human.cta_selection(category="CTA")
    first.iloc[0].gene_identities[0]["synonyms"].append("test-mutation")
    second = human.cta_selection(category="CTA")
    assert "test-mutation" not in second.iloc[0].gene_identities[0]["synonyms"]
    assert canine.cta_selection().iloc[0].gene_identities[0]["historical_gene_ids"] is None


def test_invalid_queries_fail_without_falling_back_to_human():
    with pytest.raises(ValueError, match="Unknown species"):
        cta_selection.results(species="rat")
    with pytest.raises(ValueError, match="Unknown category"):
        human.cta_selection(category="safe")
    with pytest.raises(ValueError, match="Unknown evidence"):
        mouse.cta_sources(category="typo")
    with pytest.raises(ValueError, match="antigens"):
        canine.expression.cta_selection()


def test_exhaustive_categorical_labels_and_legacy_query_alias():
    all_rows = cta_selection.results()
    other = all_rows[all_rows.category.eq("Other")]
    assert set(other.subtype) <= set(cta_selection.OTHER_SUBTYPES)
    assert other.classification.eq("Other:" + other.subtype).all()
    assert all_rows.category.notna().all(), "every human protein resolves to CTA, Unsafe or Other"
    assert all_rows.loc[~all_rows.category.eq("Other").fillna(False), "subtype"].isna().all()
    with pytest.warns(DeprecationWarning, match="Other"):
        legacy = cta_selection.results(category="Neither")
    pd.testing.assert_frame_equal(legacy, other.reset_index(drop=True))
    assert set(other.subtype) <= {
        "Nonreproductive",
        "Reproductive",
        "Diffuse",
        "Held",
        "Contained",
        "Nested",
        "Translation",
    }
    with pytest.raises(ValueError, match="subtype"):
        cta_selection.results(category="Unsafe", subtype="Low")
    with pytest.raises(ValueError, match="subtype"):
        cta_selection.results(subtype="Other")


def test_evidence_attributes_are_independent_of_the_categorical_decision():
    profiles = cta_selection.profiles().set_index("label")
    assert profiles.loc["DBX1", "category"] == "Other"
    assert profiles.loc["DBX1", "replicated_tumor_rna"]
    assert not profiles.loc["DBX1", "reproductive_evidence"]
    assert profiles.loc["DBX1", "tumor_max_tpm"] > 290
    assert profiles.loc["DBX1", "normal_max_ntpm"] < 1
    assert profiles.loc[["LCE3B", "LCE3C"], "adult_disposable_only"].all()
    assert not profiles.loc[["LCE3B", "LCE3C"], "developmental_measured"].any()
    assert cta_selection.profiles(species="dog").adult_disposable_only.isna().all()
    assert cta_selection.profiles(species="mouse").normal_max_ntpm.isna().all()


def test_full_reproductive_and_cell_support_preserve_independent_gates():
    rows = cta_selection.results().set_index("label")
    assert rows.loc[["MAGEA11", "FOXR2", "CAGE1", "PRDM9", "XKR5"], "category"].eq("CTA").all()
    assert rows.loc[["IFNA2", "IFNA14"], "classification"].eq("Unsafe").all()
    assert rows.loc[["CDY1", "CDY1B"], "classification"].eq("CTA").all()
    assert rows.loc["OR2T7", "classification"] == "Other:Translation"
    profile = cta_selection.profiles().set_index("label").loc["OR2T7"]
    assert profile.reproductive_evidence and not profile.translation_qualified

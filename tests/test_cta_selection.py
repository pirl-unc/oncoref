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
    audit = Path(__file__).resolve().parents[1] / "docs/audits/cta-release-1.21.0"
    expected = json.loads((audit / "summary.json").read_text())
    assert human.cta_selection_coverage()["counts"] == expected["categories"]
    selected = human.cta_selection(category="CTA")
    assert selected.gene_ids.map(len).sum() == expected["selected_member_loci"]
    assert {"PRAME", "MAGEA4", "NY-ESO-1", "TBXT", "XAGE1A/B", "CT83"} <= set(selected.label)
    assert "CTAG2" in set(human.cta_selection(category="Unsafe").label)
    assert selected.identity_status.eq("resolved").all()


def test_native_nominations_never_inherit_human_admission_or_human_gene_ids():
    frame = canine.cta_selection()
    assert len(frame) == 14
    assert frame.category.eq("Neither").all()
    assert frame.assessment_status.eq("unassessed").all()
    assert frame.policy_id.isna().all()
    assert frame.entity_kind.eq("native_nomination").all()
    assert not any(gid.startswith("ENSG") for ids in frame.gene_ids for gid in ids)
    assert any(len(ids) > 1 for ids in frame.gene_ids)  # Ambiguous family assays remain visible.
    assert canine.cta_selection(category="CTA").empty
    assert mouse.cta_selection().empty
    assert mouse.cta_selection_coverage()["selection_status"] == "unavailable"


def test_sources_distinguish_native_units_and_unintegrated_resources():
    dog = canine.cta_sources().drop_duplicates("source_id").set_index("source_id")
    assert dog.loc["barkbase", "native_unit"] == "TMM-CPM"
    assert dog.loc["pfizer-dog-protein", "native_unit"] == "iBAQ"
    assert dog.loc["doga", "status"] == "metadata_only"
    assert dog.loc["gse304066-osteosarcoma", "status"] == "blocked_provenance"
    assert mouse.cta_sources().status.eq("not_integrated").all()
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

"""The explicit module API preserves species boundaries and native source contracts."""

import importlib
from types import ModuleType

import pytest

import oncoref
from oncoref import hpa, hpa_cancer, species_expression
from oncoref.species import DATA_GROUPS


@pytest.mark.parametrize("name,taxon", [("human", 9606), ("canine", 9615), ("mouse", 10090)])
def test_species_modules_are_importable_and_groups_stay_scoped(name, taxon):
    module = importlib.import_module(f"oncoref.{name}")
    assert isinstance(module, ModuleType)
    assert getattr(oncoref, name) is module
    assert name in oncoref.__all__
    assert module.species.taxon_id == taxon
    assert tuple(group.name for group in module.groups) == DATA_GROUPS
    for group in module.groups:
        assert getattr(module, group.name) is group
        assert group.species.taxon_id == taxon
        assert all(row["taxon_id"] == taxon for row in group.datasets())
    assert module.inventory() == oncoref.for_species(taxon).inventory()
    for dataset in module.datasets:
        name = dataset.__name__.rsplit(".", 1)[1]
        assert isinstance(dataset, ModuleType)
        assert importlib.import_module(dataset.__name__) is dataset
        assert getattr(oncoref, name) is getattr(module, name) is dataset
        assert name in oncoref.__all__


def test_species_modules_do_not_expose_other_species_datasets():
    assert not hasattr(oncoref.canine, "human_protein_atlas")
    assert not hasattr(oncoref.human, "barkbase")
    assert oncoref.canine.antigens.coverage()["coverage"] == "partial"
    assert {m.source().id for m in oncoref.canine.datasets} == {
        s.id for s in oncoref.canine.expression.sources()
    }
    assert all(m.source().taxon_id == 9615 for m in oncoref.canine.datasets)


def test_human_atlas_reuses_existing_loaders_and_their_caches():
    atlas = oncoref.human.human_protein_atlas
    assert atlas.rna_consensus is hpa.hpa_rna_consensus
    assert atlas.normal_tissue is hpa.hpa_normal_tissue
    assert atlas.single_cell is hpa.hpa_single_cell
    assert atlas.resolve_safety_tissue_group is hpa.resolve_safety_tissue_group
    assert atlas.cancer_rna_prevalence is hpa_cancer.hpa_cancer_rna_prevalence
    assert atlas.cancer_ihc_prevalence is hpa_cancer.hpa_cancer_ihc_prevalence
    assert atlas.cancer_sources() == hpa_cancer.hpa_cancer_sources()


@pytest.mark.parametrize("dataset", oncoref.canine.datasets)
def test_named_datasets_cannot_be_redirected_to_other_species_or_sources(dataset):
    with pytest.raises(TypeError, match="species"):
        dataset.read(species="human")
    with pytest.raises(TypeError, match="source_id"):
        dataset.read(source_id="other-source")


@pytest.mark.parametrize(
    "dataset,status",
    [
        (oncoref.doga, "metadata_only"),
        (oncoref.pfizer_canine_rna, "metadata_only"),
        (oncoref.gse304066, "blocked_provenance"),
    ],
)
def test_pending_dataset_modules_report_real_gaps_before_fetching(dataset, status, monkeypatch):
    def no_fetch(*args, **kwargs):
        pytest.fail("Pending dataset tried to fetch an asset")

    monkeypatch.setattr(species_expression, "ensure_asset", no_fetch)
    assert dataset.source().acquisition_status == status
    with pytest.raises(ValueError, match=status):
        dataset.read()

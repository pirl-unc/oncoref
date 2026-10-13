# Licensed under the Apache License, Version 2.0.
import hashlib
import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest
import yaml

import oncoref
from oncoref import expression_registry as registry
from oncoref import species_assets as sa
from oncoref.expression_contract import ExpressionContract, ExpressionMeasurement
from oncoref.reference_identity import ReferenceIdentity
from oncoref.species_expression import (
    ExpressionDataset,
    combine_expression,
    read_species_expression,
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def dataset():
    reference = ReferenceIdentity(
        9615, "fixture.1", "fixture", "1", "v1", {"gene": "fixture"}, {"gtf": "a" * 64}
    )
    contract = ExpressionContract(
        9615,
        "fixture",
        "v1",
        "normal",
        reference,
        "verified",
        "fixture",
        ExpressionMeasurement("gene", "TPM", "linear", False),
        {"matrix": "b" * 64},
        "Synthetic",
    )
    values = pd.DataFrame(
        [[0, np.nan, 2], [np.nan, np.nan, 3]],
        index=pd.Index(["g1", "g2"], name="feature_id"),
        columns=pd.Index(["s1", "s2", "pool"], name="sample_id"),
    )
    samples = pd.DataFrame(
        [
            {
                "taxon_id": 9615,
                "library_id": "l1",
                "run_ids": ["r1"],
                "donor_id": "d1",
                "donor_scope": "fixture",
                "pooling": "single",
                "tissue": "heart",
            },
            {
                "taxon_id": 9615,
                "library_id": "l2",
                "run_ids": ["r2"],
                "donor_id": "d1",
                "donor_scope": "fixture",
                "pooling": "single",
                "tissue": "heart",
            },
            {
                "taxon_id": 9615,
                "library_id": "l3",
                "run_ids": ["r3"],
                "donor_id": "d2;d3",
                "donor_scope": "fixture",
                "pooling": "pooled",
                "tissue": "heart",
            },
        ],
        index=values.columns,
    )
    features = pd.DataFrame({"source_feature_id": values.index}, index=values.index)
    return ExpressionDataset(values, features, samples, contract)


def test_missing_measurements_zero_and_pooled_donors_are_distinct(dataset):
    summary = dataset.tissue_summary().set_index("feature_id")
    assert summary.loc["g1", "measured_libraries"] == 2
    assert summary.loc["g1", "known_individual_donors_measured"] == 1
    assert summary.loc["g2", "measured_libraries"] == 1
    assert summary.loc["g2", "known_individual_donors_measured"] == 0
    assert summary.loc["g2", "libraries_with_unresolved_or_pooled_donors"] == 1
    assert np.isnan(dataset.values.loc["g2", "s1"])
    assert dataset.values.loc["g1", "s1"] == 0


def test_artifact_roundtrip_checksum_and_no_overwrite(dataset, tmp_path):
    directory = dataset.save(tmp_path / "artifact")
    restored = ExpressionDataset.load(directory)
    pd.testing.assert_frame_equal(restored.values, dataset.values)
    pd.testing.assert_frame_equal(restored.samples, dataset.samples)
    assert restored.contract.as_dict() == dataset.contract.as_dict()
    with pytest.raises(FileExistsError):
        dataset.save(directory)
    (directory / "values.parquet").write_bytes(b"changed")
    with pytest.raises(ValueError, match="checksum mismatch"):
        ExpressionDataset.load(directory)


def test_nested_metadata_is_defensive(dataset):
    samples = dataset.samples
    samples.iloc[0]["run_ids"].append("bad")
    assert dataset.samples.iloc[0]["run_ids"] == ["r1"]
    other = ExpressionDataset(dataset.values, dataset.features, samples, dataset.contract)
    samples.iloc[0]["run_ids"].append("later")
    assert other.samples.iloc[0]["run_ids"] == ["r1", "bad"]


def test_combining_rejects_reused_libraries_and_unverified_references(dataset):
    with pytest.raises(ValueError, match="Shared runs"):
        combine_expression([dataset, dataset])
    unresolved = ExpressionDataset(
        dataset.values,
        dataset.features,
        dataset.samples,
        replace(dataset.contract, reference=None, reference_status="unresolved"),
    )
    with pytest.raises(ValueError, match="unresolved"):
        combine_expression([unresolved])
    samples = dataset.samples
    samples["run_ids"] = [[], [], []]
    unknown = ExpressionDataset(dataset.values, dataset.features, samples, dataset.contract)
    with pytest.raises(ValueError, match="run lineage"):
        combine_expression([unknown])


@pytest.fixture
def source_fixture(tmp_path, monkeypatch):
    root = tmp_path / "species" / "dog"
    root.mkdir(parents=True)
    monkeypatch.setattr(registry, "_DATA_DIR", tmp_path)
    monkeypatch.setattr(sa, "_MANIFEST_ROOT", tmp_path / "species")
    registry.load_registry.cache_clear()

    def make(adapter, text):
        path = tmp_path / "matrix.tsv"
        path.write_text(text)
        metadata = {
            "schema": "oncoref.expression-samples",
            "schema_version": 1,
            "taxon_id": 9615,
            "source_id": "fixture",
            "samples": [
                {
                    "sample_id": s,
                    "taxon_id": 9615,
                    "library_id": s,
                    "run_ids": [],
                    "donor_id": None,
                    "donor_scope": "fixture",
                    "pooling": "unknown",
                    "tissue": "testis",
                }
                for s in ["s1", "s2"]
            ],
        }
        sample_file = root / "samples.json"
        sample_file.write_text(json.dumps(metadata))
        unit = "TPM" if adapter == "bgee" else "TMM-CPM"
        contract = ExpressionContract(
            9615,
            "fixture",
            "v1",
            "normal",
            None,
            "unresolved",
            "fixture",
            ExpressionMeasurement("gene", unit, "linear", False),
            {"matrix": sha(path), "samples": sha(sample_file)},
            "Synthetic unresolved reference",
        )
        source = {
            "id": "fixture",
            "taxon_id": 9615,
            "category": "expression",
            "context": "normal",
            "cancer_codes": [],
            "source_type": adapter,
            "accession": "STUDY",
            "source_version": "v1",
            "unit": unit,
            "source_scale_class": "linear",
            "tpm_proxy": False,
            "contract": contract.as_dict(),
            "matrix_asset": "fixture",
            "sample_metadata": "samples.json",
            "adapter": adapter,
        }
        (root / "expression_sources.yaml").write_text(yaml.safe_dump({"sources": [source]}))
        asset = sa.SpeciesAsset(
            "fixture",
            9615,
            "expression",
            "matrix.tsv",
            "https://example.org/matrix.tsv",
            sha(path),
            path.stat().st_size,
            "v1",
            "Synthetic",
            "Synthetic",
        )
        (root / "assets.json").write_text(
            json.dumps(
                {
                    "schema": "oncoref.species-assets",
                    "schema_version": 1,
                    "taxon_id": 9615,
                    "assets": [asset.as_dict()],
                    "references": [],
                }
            )
        )
        registry.load_registry.cache_clear()
        return {"9615:fixture": path}

    yield make
    registry.load_registry.cache_clear()


def test_bgee_native_reader_preserves_absence(source_fixture):
    overrides = source_fixture(
        "bgee",
        "Experiment ID\tLibrary ID\tGene ID\tTPM\nSTUDY\ts1\tg1\t0\nSTUDY\ts2\tg1\t2\nSTUDY\ts1\tg2\t\n",
    )
    data = read_species_expression(
        "fixture", species="dog", offline=True, local_overrides=overrides
    )
    assert data.values.loc["g1", "s1"] == 0
    assert np.isnan(data.values.loc["g2", "s2"])
    assert data.contract.reference_status == "unresolved"
    with pytest.raises(KeyError):
        read_species_expression("fixture", species="mouse", offline=True)


def test_bgee_rejects_duplicate_measurements(source_fixture):
    overrides = source_fixture(
        "bgee",
        "Experiment ID\tLibrary ID\tGene ID\tTPM\nSTUDY\ts1\tg1\t0\nSTUDY\ts1\tg1\t2\nSTUDY\ts2\tg1\t3\n",
    )
    with pytest.raises(ValueError, match="Duplicate Bgee"):
        read_species_expression("fixture", species="dog", offline=True, local_overrides=overrides)


def test_barkbase_keeps_unmapped_features_and_aliases(source_fixture):
    overrides = source_fixture(
        "barkbase",
        "StringTieID\tEnsemblGeneID\tGeneSymbol\ts1\ts2\nMSTRG.1\t\t\t0\t\nMSTRG.2\tg1;g2\talias\t1\t2\n",
    )
    data = read_species_expression(
        "fixture", species="dog", offline=True, local_overrides=overrides
    )
    assert len(data.values) == 2
    assert data.features.loc["MSTRG.2", "EnsemblGeneID"] == "g1;g2"
    assert data.contract.measurement.unit == "TMM-CPM"
    assert np.isnan(data.values.loc["MSTRG.1", "s2"])
    with pytest.raises(ValueError, match="cannot establish"):
        data.contract.require_feature_level("transcript")


def test_registered_corrections_and_lineage_evidence():
    dog = oncoref.for_species("dog")
    assert len(dog.expression.sources(context="normal")) == 6
    assert [s.id for s in dog.expression.sources(context="tumor")] == ["gse304066-osteosarcoma"]
    source = registry.expression_source("bgee-srp181218", species="dog")
    file = (
        registry.expression_source_registry_path(species="dog").parent
        / "bgee-srp181218-samples.json"
    )
    assert sha(file) == source.expression_contract.content_provenance["samples"]
    samples = json.loads(file.read_text())["samples"]
    corrections = [s for s in samples if s["tissue_correction"]]
    assert len(corrections) == 5
    assert all(
        s["original_tissue"] == "colon" and s["tissue"] == "heart" and s["evidence"]
        for s in corrections
    )
    overlap = dog.expression.lineage_overlaps()
    assert len(overlap["shared_run_ids"]) == 46
    assert overlap["independent_donors"] is None


def test_protein_groups_and_missing_ibaq_survive_without_allocation(tmp_path):
    from oncoref.species_expression import _read_pfizer_protein

    file = tmp_path / "protein.tsv"
    file.write_text(
        "Run\tProtein\tChannel\tBioReplicate\tCondition\tTechRepMixture\tMixture\tsample.type\tProteinGroup\tIBAQ\n"
        "Pool1_1\tP1\t126\tpool-a\theart\t1\tPool1\tsample\tP1;P2\t0\n"
        "Pool1_1\tP1\t127\tpool-b\tliver\t1\tPool1\tsample\tP1;P2\tNA\n"
    )
    values, features = _read_pfizer_protein(file)
    assert values.shape == (1, 2)
    assert values.isna().sum().sum() == 1
    assert values.iloc[0].dropna().tolist() == [0]
    assert features.index.tolist() == ["P1;P2"]
    assert features.Protein.tolist() == ["P1"]
    file.write_text(file.read_text() + file.read_text().splitlines()[1] + "\n")
    with pytest.raises(ValueError, match="silently aggregated"):
        _read_pfizer_protein(file)


def test_registered_protein_scope_and_pooled_samples():
    source = registry.expression_source("pfizer-dog-protein", species="dog")
    assert source.expression_contract.measurement.feature_level == "protein_group"
    assert source.unit == "iBAQ"
    assert not source.linear_tpm_comparable
    assert source.expression_contract.reference is None
    assert (
        source.expression_contract.content_provenance["searched_fasta"]
        == sa.species_asset("pfizer-dog-search-fasta", species="dog").sha256
    )
    file = (
        registry.expression_source_registry_path(species="dog").parent
        / "pfizer-dog-protein-samples.json"
    )
    samples = json.loads(file.read_text())["samples"]
    assert len(samples) == 38
    assert all(row["pooling"] == "pooled" and row["donor_id"] is None for row in samples)
    assert all(row["source_condition"]["presentation_evidence"] is False for row in samples)


@pytest.mark.parametrize("name", ["doga", "pfizer-dog-rna", "gse304066-osteosarcoma"])
def test_pending_sources_are_inspectable_but_cannot_be_read(name, monkeypatch):
    import oncoref.species_expression as module

    monkeypatch.setattr(
        module,
        "ensure_asset",
        lambda *a, **kw: pytest.fail("unavailable source attempted download"),
    )
    source = registry.expression_source(name, species="dog")
    assert source.acquisition_status in {"metadata_only", "blocked_provenance"}
    assert source.expression_contract.measurement.unit is None
    with pytest.raises(ValueError, match="unavailable"):
        read_species_expression(name, species="dog")


def test_reused_run_inside_one_dataset_cannot_inflate_coverage(dataset):
    samples = dataset.samples
    samples.at["s2", "run_ids"] = ["r1"]
    with pytest.raises(ValueError, match="Run reused"):
        ExpressionDataset(dataset.values, dataset.features, samples, dataset.contract)

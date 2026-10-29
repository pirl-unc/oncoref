# Licensed under the Apache License, Version 2.0.
"""Exercise the public producer contract and its failure boundaries using tiny releases."""

import copy
import hashlib
from dataclasses import replace

import pandas as pd
import pytest

from oncoref import expression_registry
from oncoref import managed_expression as managed
from oncoref.expression_contract import ExpressionContract
from oncoref.species_expression import combine_expression


def test_optional_dataset_discovery_and_cli_failure_without_dependency(monkeypatch, capsys):
    from oncoref import catalog, cli
    from oncoref.species_assets import asset_inventory, inspect_asset, species_asset

    def missing(cache_root):
        raise ImportError("Install oncoref[quantitative]")

    monkeypatch.setattr(managed, "_registry", missing)
    asset = species_asset("pfizer-dog-protein", species="dog")
    assert asset.size_bytes == 84379381
    assert asset.as_dict()["installed_size_bytes"] == 1404428288
    assert inspect_asset(asset)["status"] == "dependency_missing"
    assert any(row["name"] == asset.id for row in asset_inventory("dog"))
    assert catalog.versions(asset.id, species="dog") == ("2026-10-12.1",)
    assert catalog.path(asset.id, species="dog") is None
    assert catalog.provenance(asset.id, species="dog")["verification_state"] == "dependency_missing"
    assert catalog.status(asset.id, species="dog")[0]["present"] is False
    assert cli.main(["data", "metadata", asset.id, "--species", "dog"]) == 0
    capsys.readouterr()
    assert cli.main(["data", "fetch", asset.id, "--species", "dog"]) == 1
    output = capsys.readouterr()
    assert "oncoref[quantitative]" in output.err and "Traceback" not in output.err


@pytest.fixture
def release(tmp_path, monkeypatch):
    pytest.importorskip("orfdata")
    from orfdata import DatasetRegistry, QuantitativeDataset
    from orfdata.measurements import Assay, DatasetManifest, Measurement, Sample, TargetMapping
    from orfdata.reference import ScopedTarget

    entry = next(
        row
        for row in expression_registry.expression_source_registry_entries(species="dog")
        if row["id"] == "pfizer-dog-protein"
    )
    spec = managed.specification(entry)
    manifest = copy.deepcopy(spec["manifest"])
    manifest.update(dataset_id="fixture-protein", version="fixture-1", title="Synthetic release")
    manifest["coverage"] = {"fixture": True}
    manifest = DatasetManifest.from_dict(manifest)
    reference = manifest.references[0]
    target = ScopedTarget(reference, "protein_group", "group-1", spec["selection"]["namespace"])
    attrs = {"Run": "pool1", "Channel": "126"}
    samples = [
        Sample(name, reference, tissue="heart", attributes={"source_attributes": attrs})
        for name in ("s1", "s2", "s3")
    ]
    assays = [
        Assay(s.sample_id, s.sample_id, "digest_ms", reference, run_id="pool1", channel=str(i))
        for i, s in enumerate(samples)
    ]
    observations = [
        Measurement(
            f"row:{i}:IBAQ",
            target,
            sid,
            sid,
            "protein_abundance",
            value,
            "iBAQ",
            manifest.dataset_id,
            manifest.version,
            status=status,
            provenance={
                "source_asset": "protein.tsv",
                "line": i,
                "native_value": native,
                "representative": "P1",
            },
        )
        for i, sid, value, status, native in [
            (2, "s1", 0.0, "measured", "0"),
            (3, "s2", None, "missing", "NA"),
        ]
    ]
    mapping = TargetMapping(
        "fixture:group:1",
        target.namespace,
        target.identifier,
        reference,
        (
            ScopedTarget(reference, "protein", "P1", "uniprot"),
            ScopedTarget(reference, "protein", "P2", "uniprot"),
        ),
        "ambiguous",
        "reported membership",
        evidence={"unique_protein_detection": False},
    )
    path = tmp_path / "dataset.sqlite"
    QuantitativeDataset.create(
        path,
        manifest,
        samples=samples,
        assays=assays,
        measurements=observations,
        mappings=[mapping],
    )
    pin = {
        "url": "https://example.org/fixture.sqlite",
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "size": path.stat().st_size,
    }
    spec["manifest"] = manifest.to_dict()
    spec["download"] = {
        "assets": {"dataset.sqlite": pin},
        "distribution": {"access_terms": "Synthetic fixture", "redistributable": True},
    }
    definitions = {
        manifest.dataset_id: {
            "default_version": manifest.version,
            "versions": {manifest.version: {"manifest": manifest.to_dict(), **spec["download"]}},
        }
    }
    monkeypatch.setattr(
        managed,
        "_registry",
        lambda cache_root: DatasetRegistry(definitions, cache_root=tmp_path / "cache"),
    )
    monkeypatch.setattr(managed, "specification", lambda unused: copy.deepcopy(spec))
    contract = ExpressionContract.from_dict(entry["contract"])
    contract = replace(
        contract,
        content_provenance={
            **contract.content_provenance,
            "quantitative_manifest": managed._digest(manifest.to_dict()),
            "quantitative_artifact": pin["sha256"],
        },
    )
    curated = pd.DataFrame(
        [
            {
                "sample_id": sample.sample_id,
                "taxon_id": 9615,
                "library_id": sample.sample_id,
                "run_ids": [],
                "donor_id": None,
                "donor_scope": "fixture",
                "pooling": "pooled",
                "tissue": "heart",
                "source_tissue": "heart",
                "source_attributes": attrs,
            }
            for sample in samples
        ]
    ).set_index("sample_id")
    options = {
        "offline": True,
        "cache_root": None,
        "local_overrides": {f"9615:{manifest.dataset_id}": path},
        "max_bytes": 0,
    }
    return entry, spec, contract, curated, options, path


def test_native_summary_missing_absent_shared_and_portable(release, tmp_path):
    from oncoref.species_expression import ExpressionDataset

    entry, spec, contract, samples, options, _ = release
    dataset = managed.read_managed_expression(entry, samples, contract, **options)
    assert dataset.values.loc["group-1", "s1"] == 0
    assert dataset.values.loc["group-1", ["s2", "s3"]].isna().all()
    observations = list(dataset.native_records("measurements"))
    assert len(observations) == 2  # Absent s3 does not acquire a fabricated measurement.
    assert observations[1]["status"] == "missing"
    assert observations[1]["provenance"]["native_value"] == "NA"
    (mapping,) = dataset.native_records("mappings")
    assert mapping["status"] == "ambiguous" and len(mapping["targets"]) == 2
    assert dataset.features.index.tolist() == ["group-1"]
    assert dataset.features.Protein.tolist() == ["P1"]
    assert dataset.tissue_summary().known_individual_donors_measured.tolist() == [0]
    restored = ExpressionDataset.load(dataset.save(tmp_path / "export"))
    assert list(restored.native_records("measurements")) == observations
    assert list(restored.native_records("mappings")) == [mapping]
    assert restored.provenance["manifest"] == spec["manifest"]
    with pytest.raises(ValueError, match="unresolved"):
        combine_expression([dataset, restored])


@pytest.mark.parametrize(
    "field", ["taxon", "artifact", "manifest", "source", "measurement", "sample"]
)
def test_incompatible_identity_rejected_before_measurements(release, field):
    entry, _spec, contract, samples, options, _ = release
    if field == "taxon":
        contract = replace(contract, taxon_id=10090)
    elif field in {"artifact", "manifest", "source"}:
        role = {
            "artifact": "quantitative_artifact",
            "manifest": "quantitative_manifest",
            "source": "matrix",
        }[field]
        contract = replace(
            contract, content_provenance={**contract.content_provenance, role: "0" * 64}
        )
    elif field == "measurement":
        contract = replace(contract, measurement=replace(contract.measurement, unit="TPM"))
    else:
        samples.at["s1", "source_tissue"] = "liver"
    with pytest.raises(ValueError, match=r"disagree|matching reference"):
        managed.read_managed_expression(entry, samples, contract, **options)


def test_override_checksum_failure_never_repairs_or_downloads(release):
    entry, _, contract, samples, options, path = release
    path.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="checksum mismatch"):
        managed.read_managed_expression(entry, samples, contract, **{**options, "offline": False})
    assert path.read_bytes() == b"corrupt"


def test_registry_drift_is_not_an_implicit_upgrade(release):
    _, spec, _, _, options, _ = release
    changed = copy.deepcopy(spec)
    changed["manifest"]["schema"] = "orfdata.measurements.v999"
    with pytest.raises(ValueError, match="pins disagree"):
        managed.open_managed(changed, taxon_id=9615, **options)


def test_missing_managed_release_is_offline_and_budget_bounded(release):
    _, spec, _, _, options, _ = release
    options = {**options, "local_overrides": None}
    with pytest.raises(FileNotFoundError):
        managed.open_managed(spec, taxon_id=9615, **options)
    with pytest.raises(ValueError, match="budget"):
        managed.open_managed(spec, taxon_id=9615, **{**options, "offline": False})

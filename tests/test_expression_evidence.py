# Licensed under the Apache License, Version 2.0.
"""Bounded common APIs against the producer's audited native/synthetic fixture."""

import copy
import gzip
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from oncoref import expression_registry
from oncoref import managed_expression as managed
from oncoref.expression_contract import ExpressionContract
from oncoref.expression_records import ExpressionRecords
from oncoref.species_expression import ExpressionDataset, protein_sample_id

DIRECTORY = Path(__file__).parent / "data/pfizer-assignments"
SOURCE = "pfizer-dog-peptide-evidence"
FEATURE = "A0A5F4CET6"


@pytest.fixture
def evidence(tmp_path, monkeypatch):
    pytest.importorskip("orfdata")
    from orfdata import DatasetRegistry, QuantitativeDataset

    receipt = json.loads((DIRECTORY / "receipt.json").read_text())
    compressed = (DIRECTORY / "dataset.sqlite.gz").read_bytes()
    assert hashlib.sha256(compressed).hexdigest() == receipt["compressed_sha256"]
    raw = gzip.decompress(compressed)
    assert hashlib.sha256(raw).hexdigest() == receipt["dataset_sha256"]
    path = tmp_path / "dataset.sqlite"
    path.write_bytes(raw)
    native = QuantitativeDataset(path)
    entry = next(
        e
        for e in expression_registry.expression_source_registry_entries(species="dog")
        if e["id"] == "pfizer-dog-protein"
    )
    summary = managed.specification(entry)
    baseline = dict(native.manifest.processing["protein_summary"])
    manifest = summary["manifest"]
    manifest.update(
        dataset_id=baseline["dataset_id"],
        version=baseline["version"],
        references=[r.to_dict() for r in native.manifest.references],
    )
    manifest["coverage"]["source_counts"] = dict(baseline["source_counts"])
    for source in manifest["sources"]:
        if source["asset_id"] == "protein.tsv":
            source.update(sha256=baseline["source_sha256"])
        else:
            source.update(sha256=native.manifest.references[0].sequence_sha256)
    spec = {
        "manifest": native.manifest.to_dict(),
        "download": {
            "assets": {
                "dataset.sqlite": {
                    "url": "https://example.org/fixture.sqlite",
                    "sha256": receipt["dataset_sha256"],
                    "size": len(raw),
                }
            },
            "distribution": {"redistributable": True, "access_terms": "Explicit test fixture"},
        },
    }
    # The fixture importer explicitly uses regression coverage; it is not the full release.
    spec["manifest"]["coverage"]["accession"] = entry["accession"]
    # Regenerate the fixture's manifest so the real registry/open gate remains exercised.
    import sqlite3

    with sqlite3.connect(path) as connection:
        connection.execute(
            "UPDATE metadata SET value=? WHERE name='manifest'", (json.dumps(spec["manifest"]),)
        )
    pin = spec["download"]["assets"]["dataset.sqlite"]
    pin.update(sha256=hashlib.sha256(path.read_bytes()).hexdigest(), size=path.stat().st_size)
    native = QuantitativeDataset(path)
    name, version = spec["manifest"]["dataset_id"], spec["manifest"]["version"]
    registry = DatasetRegistry(
        {
            name: {
                "default_version": version,
                "versions": {version: {"manifest": spec["manifest"], **spec["download"]}},
            }
        },
        cache_root=tmp_path / "cache",
    )
    monkeypatch.setattr(managed, "_registry", lambda unused: registry)
    monkeypatch.setattr(managed, "specification", lambda unused: copy.deepcopy(summary))
    monkeypatch.setattr(managed, "evidence_specification", lambda *unused: copy.deepcopy(spec))
    contract = ExpressionContract.from_dict(entry["contract"])
    contract = replace(
        contract,
        content_provenance={
            **contract.content_provenance,
            "quantitative_manifest": managed._digest(manifest),
            "quantitative_artifact": summary["download"]["assets"]["dataset.sqlite"]["sha256"],
            "matrix": baseline["source_sha256"],
            "searched_fasta": native.manifest.references[0].sequence_sha256,
        },
    )
    rows = pd.read_csv(DIRECTORY / "protein.tsv", sep="\t", dtype=str, keep_default_na=False)
    rows["sample_id"] = rows.apply(lambda r: protein_sample_id(r.to_dict()), axis=1)
    values = pd.DataFrame(
        [[None if v == "NA" else float(v) for v in rows.IBAQ]],
        index=[FEATURE],
        columns=rows.sample_id,
    )
    values.index.name, values.columns.name = "feature_id", "sample_id"
    samples = (
        pd.DataFrame(
            json.loads(
                (
                    expression_registry.expression_source_registry_path(species="dog").parent
                    / entry["sample_metadata"]
                ).read_text()
            )["samples"]
        )
        .set_index("sample_id")
        .loc[values.columns]
    )
    features = pd.DataFrame({"Protein": [FEATURE]}, index=values.index)
    provenance = {
        "producer": "orfdata",
        "manifest": manifest,
        "download": summary["download"],
        "selection": summary["selection"],
        "experiment": entry["accession"],
    }
    original = ExpressionDataset(
        values,
        features,
        samples,
        contract,
        records=ExpressionRecords(provenance, lambda kind: iter(())),
    )
    return original, spec, native, {f"9615:{name}": path}


def enrich(evidence, **options):
    expression, _, _, overrides = evidence
    return expression.with_evidence(
        SOURCE, feature_ids=[FEATURE], local_overrides=overrides, **options
    )


def test_complete_candidates_native_ids_and_offline_replay(evidence, tmp_path):
    original, _, native, _ = evidence
    result = enrich(evidence)
    pd.testing.assert_frame_equal(result.values, original.values)
    pd.testing.assert_frame_equal(result.samples, original.samples)
    pd.testing.assert_frame_equal(result.features, original.features)
    assert result.contract == original.contract
    assert result.values.isna().sum().sum() == 13
    mappings = list(result.native_records("mappings"))
    reported = [m for m in mappings if m["evidence"]["basis"] == "reported"]
    assert any(m["status"] == "unresolved" and len(m["targets"]) == 2 for m in reported)
    assert all(m["evidence"]["unique_protein_detection"] is False for m in mappings)
    original_mappings = [r.to_dict() for r in native.mappings()]
    for mapping in mappings:
        restored = copy.deepcopy(mapping)
        context = restored["provenance"].pop("oncoref_native_context", None)
        assert restored in original_mappings
        if mapping["evidence"]["basis"] == "reported":
            assert context and context["record"]
            assert len(context["reported_links"]) == mapping["evidence"]["membership_records"]
    relations = list(result.native_records("relations"))
    assert relations and relations == [
        r.to_dict() for r in native.relations() if r.left.identifier == "1"
    ]
    (context,) = result.provenance["enrichment"]["contexts"]
    assert context["group_id"] == FEATURE
    assert context["native_proteins"][0]["group_links"][0]["left"]["record"]["ProteinGroupID"] == 1
    assert list(result.native_records("measurements")) == []  # No extra abundance observation.
    assert original.provenance.get("enrichment") is None
    replay = ExpressionDataset.load(result.save(tmp_path / "portable"))
    assert replay.provenance == result.provenance
    for kind in ("mappings", "relations", "measurements"):
        assert list(replay.native_records(kind)) == list(result.native_records(kind))
    pd.testing.assert_frame_equal(replay.values, original.values)
    with pytest.raises(ValueError, match="Already enriched"):
        replay.with_evidence(SOURCE, feature_ids=[FEATURE])


def test_source_discovery_is_optional_and_transfer_is_explicit(evidence, monkeypatch):
    expression, spec, _, _ = evidence
    assert expression.evidence_sources()[0]["id"] == SOURCE

    def absent(*args, **kwargs):
        raise FileNotFoundError("missing fixture")

    registry = managed._registry(None)
    monkeypatch.setattr(registry, "open", absent)
    with pytest.raises(FileNotFoundError):
        expression.with_evidence(SOURCE, feature_ids=[FEATURE])
    with pytest.raises(ValueError, match="transfer bytes"):
        expression.with_evidence(SOURCE, feature_ids=[FEATURE], offline=False)
    with pytest.raises(ValueError, match="max_records"):
        enrich(evidence, max_records=1)
    spec["download"]["assets"]["dataset.sqlite"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match=r"pins disagree|checksum"):
        enrich(evidence)


@pytest.mark.parametrize("field", ["reference", "taxon", "experiment", "summary", "representative"])
def test_incompatible_identity_rejected(evidence, field):
    expression, spec, _, _ = evidence
    if field == "reference":
        spec["manifest"]["references"][0]["sequence_sha256"] = "f" * 64
    elif field == "taxon":
        spec["manifest"]["references"][0]["taxon_id"] = 10090
    elif field == "experiment":
        spec["manifest"]["coverage"]["accession"] = "different study"
    elif field == "summary":
        expression._records._provenance["manifest"]["version"] = "different generation"
    else:
        expression._features.at[FEATURE, "Protein"] = "different representative"
    with pytest.raises(ValueError, match=r"differs|disagrees"):
        enrich(evidence)


@pytest.mark.parametrize("selection", [[], [FEATURE, FEATURE], ["absent"], [None], FEATURE])
def test_explicit_native_selection_required(evidence, selection):
    expression, _, _, _ = evidence
    with pytest.raises((ValueError, TypeError)):
        expression.with_evidence(SOURCE, feature_ids=selection)


def test_native_scalar_encoding_is_lossless():
    import base64
    import struct

    from oncoref.pfizer_canine_protein import _portable_native

    original = {
        "blob": b"\x00\xff",
        "identifier": -9222281927842617047,
        "missing": None,
        "nonfinite": float("inf"),
        "zero": 0,
    }
    result = _portable_native(original)
    json.dumps(result, allow_nan=False)
    assert base64.b64decode(result["blob"]["base64"]) == original["blob"]
    assert int(result["identifier"]["decimal"]) == original["identifier"]
    assert struct.unpack("<d", bytes.fromhex(result["nonfinite"]["hex"]))[0] == float("inf")
    assert result["missing"] is None and result["zero"] == 0

# Licensed under the Apache License, Version 2.0.
import hashlib
import json
from dataclasses import replace

import pytest

import oncoref
from oncoref import cli
from oncoref import species_assets as sa
from oncoref.reference_identity import ReferenceIdentity


@pytest.fixture
def asset():
    return sa.SpeciesAsset(
        "fixture",
        9615,
        "references",
        "fixture.gtf.gz",
        "https://example.org/fixture.gtf.gz",
        hashlib.sha256(b"fixture").hexdigest(),
        7,
        "v1",
        "Synthetic fixture",
        "Synthetic test data",
    )


def test_real_reference_manifests_are_distinct_and_complete():
    rows = oncoref.for_species("dog").references.registered_references()
    assert {r["id"] for r in rows} == {"refseq106", "ensembl98"}
    a, b = [ReferenceIdentity.from_dict(r["identity"]) for r in rows]
    assert a.reference_id != b.reference_id
    assert {a.assembly_accession, b.assembly_accession} == {"GCF_011100685.1", "GCA_000002285.2"}
    assert len(sa.species_assets("dog")) == 8
    assert sa.species_assets("mouse") == ()
    for row in rows:
        identity = ReferenceIdentity.from_dict(row["identity"])
        for role, name in row["assets"].items():
            assert sa.species_asset(name, species="dog").sha256 == identity.content_provenance[role]
    with pytest.raises(KeyError):
        sa.species_asset("refseq106-gtf", species="mouse")
    with pytest.raises(ValueError, match="belongs to"):
        oncoref.for_species("dog").expression.ensure("refseq106-gtf")


def test_override_reuse_is_read_only_and_rejects_corruption(asset, tmp_path, monkeypatch):
    import datacache

    monkeypatch.setattr(datacache, "fetch_file", lambda *a, **kw: pytest.fail("network"))
    original = tmp_path / "existing.gz"
    original.write_bytes(b"fixture")
    options = {"local_overrides": {"9615:fixture": original}, "cache_root": tmp_path / "unused"}
    assert sa.ensure_asset(asset, offline=True, max_bytes=0, **options) == original
    assert not (tmp_path / "unused").exists()
    original.write_bytes(b"corrupt")
    for offline in (True, False):
        with pytest.raises(ValueError, match="fixture"):
            sa.ensure_asset(asset, offline=offline, **options)
    assert original.read_bytes() == b"corrupt"


def test_cached_offline_replay_verifies_bytes(asset, tmp_path):
    path = sa.asset_path(asset, cache_root=tmp_path)
    assert not tmp_path.joinpath("species").exists()
    with pytest.raises(ValueError, match="missing"):
        sa.ensure_asset(asset, offline=True, cache_root=tmp_path)
    path.parent.mkdir(parents=True)
    path.write_bytes(b"fixture")
    assert sa.ensure_asset(asset, offline=True, cache_root=tmp_path) == path
    path.write_bytes(b"corrupt")
    with pytest.raises(ValueError):
        sa.ensure_asset(asset, offline=True, cache_root=tmp_path)


def test_transfer_budget_and_datacache_ownership(asset, tmp_path, monkeypatch):
    import datacache

    calls = []

    def fetch(url, **options):
        calls.append((url, options))
        assert options["raw"] and options["resume"]
        assert options["expected_sha256"] == asset.sha256
        assert options["expected_size"] == asset.size_bytes
        return options["destination"]

    monkeypatch.setattr(datacache, "fetch_file", fetch)
    with pytest.raises(ValueError, match="transfer budget"):
        sa.ensure_asset(asset, cache_root=tmp_path, max_bytes=6)
    assert not calls
    sa.ensure_asset(asset, cache_root=tmp_path, max_bytes=7)
    assert len(calls) == 1


def test_reference_replay_rejects_missing_roles_and_changed_hashes(asset, tmp_path, monkeypatch):
    root = tmp_path / "manifests" / "dog"
    root.mkdir(parents=True)
    identity = ReferenceIdentity(
        9615, "fixture.1", "fixture", "1", "v1", {"gene": "fixture"}, {"gtf": asset.sha256}
    )
    manifest = {
        "schema": "oncoref.species-assets",
        "schema_version": 1,
        "taxon_id": 9615,
        "assets": [asset.as_dict()],
        "references": [
            {
                "id": "fixture",
                "identity": identity.as_dict(),
                "assets": {"gtf": asset.id},
                "notes": "Synthetic",
            }
        ],
    }
    file = root / "assets.json"
    file.write_text(json.dumps(manifest))
    monkeypatch.setattr(sa, "_MANIFEST_ROOT", root.parent)
    local = tmp_path / "existing"
    local.write_bytes(b"fixture")
    reference, paths = sa.ensure_reference(
        "fixture", species="dog", offline=True, local_overrides={"9615:fixture": local}, max_bytes=0
    )
    assert reference == identity and paths["gtf"] == local
    manifest["references"][0]["assets"] = {}
    file.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="roles disagree"):
        sa.registered_references("dog")
    manifest["references"][0]["assets"] = {"gtf": asset.id}
    manifest["assets"][0]["sha256"] = "b" * 64
    file.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="digest disagree"):
        sa.registered_references("dog")


@pytest.mark.parametrize(
    "changes",
    [
        {"id": "../escape"},
        {"filename": "/tmp/escape"},
        {"size_bytes": True},
        {"sha256": "bad"},
        {"url": "file:///tmp/x"},
        {"taxon_id": True},
    ],
)
def test_bad_asset_contracts_fail(asset, changes):
    with pytest.raises(ValueError):
        replace(asset, **changes)


def test_cli_offline_does_not_fall_back_to_human(tmp_path, capsys):
    assert (
        cli.main(
            [
                "data",
                "fetch",
                "refseq106-gtf",
                "--species",
                "dog",
                "--offline",
                "--cache-root",
                str(tmp_path),
            ]
        )
        == 1
    )
    assert "missing" in capsys.readouterr().err
    assert cli.main(["data", "metadata", "ensembl98-gtf", "--species", "dog"]) == 0
    assert json.loads(capsys.readouterr().out)[0]["taxon_id"] == 9615

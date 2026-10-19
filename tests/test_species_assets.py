# Licensed under the Apache License, Version 2.0.
import hashlib
import io
import json
import os
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

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


@pytest.fixture(params=["native", "without_resume"])
def acquisition_mode(request, monkeypatch):
    if request.param == "without_resume":
        from datacache import resume

        # Simulate the platform decision without changing os.name globally:
        # pathlib and the real datacache file I/O must still use the host OS.
        monkeypatch.setattr(sa, "os", SimpleNamespace(name="nt", environ=os.environ))

        def unsupported_resume(*args, **kwargs):
            raise NotImplementedError("Resumable downloads require a POSIX local filesystem")

        monkeypatch.setattr(resume, "validate_resume", unsupported_resume)


@pytest.fixture
def download_response(monkeypatch):
    import requests

    state = SimpleNamespace(body=b"fixture", calls=[])

    def get(url, **kwargs):
        state.calls.append((url, kwargs))
        response = requests.Response()
        response.status_code = 200
        response.url = url
        response.headers["Content-Length"] = str(len(state.body))
        response.raw = io.BytesIO(state.body)
        return response

    # Replace only HTTP transport; datacache still stages, checks, publishes,
    # and records provenance using its production implementation.
    monkeypatch.setattr(requests, "get", get)
    return state


@pytest.mark.parametrize("cached_bytes", [None, b"corrupt", b"short", b"fixture"])
def test_online_acquisition_verifies_and_repairs_managed_files(
    asset, tmp_path, acquisition_mode, download_response, cached_bytes
):
    path = sa.asset_path(asset, cache_root=tmp_path)
    if cached_bytes is not None:
        path.parent.mkdir(parents=True)
        path.write_bytes(cached_bytes)
    if cached_bytes != b"fixture":
        with pytest.raises(ValueError, match="transfer budget"):
            sa.ensure_asset(asset, cache_root=tmp_path, max_bytes=asset.size_bytes - 1)
        assert not download_response.calls
        assert path.read_bytes() == cached_bytes if cached_bytes is not None else not path.exists()

    assert sa.ensure_asset(asset, cache_root=tmp_path, max_bytes=asset.size_bytes) == path
    assert path.read_bytes() == b"fixture"
    assert sa.inspect_asset(asset, cache_root=tmp_path)["status"] == "available"
    assert len(download_response.calls) == (0 if cached_bytes == b"fixture" else 1)
    if download_response.calls:
        from datacache import provenance

        record = json.loads(Path(provenance.sidecar_path(path)).read_text())
        assert record["sha256"] == asset.sha256
    count = len(download_response.calls)
    assert sa.ensure_asset(asset, cache_root=tmp_path, offline=True, max_bytes=0) == path
    assert sa.ensure_asset(asset, cache_root=tmp_path, max_bytes=0) == path
    assert len(download_response.calls) == count


@pytest.mark.parametrize("cached_bytes", [None, b"corrupt"])
@pytest.mark.parametrize("response_bytes", [b"changed", b"too many bytes"])
def test_failed_acquisition_preserves_destination(
    asset, tmp_path, acquisition_mode, download_response, cached_bytes, response_bytes
):
    from datacache import FileValidationError

    path = sa.asset_path(asset, cache_root=tmp_path)
    if cached_bytes is not None:
        path.parent.mkdir(parents=True)
        path.write_bytes(cached_bytes)
    download_response.body = response_bytes
    with pytest.raises(FileValidationError):
        sa.ensure_asset(asset, cache_root=tmp_path)
    assert download_response.calls
    assert path.read_bytes() == cached_bytes if cached_bytes is not None else not path.exists()


def test_cli_fetch_repairs_managed_corruption_without_force(
    asset, tmp_path, monkeypatch, capsys, acquisition_mode, download_response
):
    monkeypatch.setattr(sa, "species_assets", lambda species: (asset,))
    path = sa.asset_path(asset, cache_root=tmp_path)
    path.parent.mkdir(parents=True)
    path.write_bytes(b"corrupt")
    assert (
        cli.main(
            [
                "data",
                "fetch",
                "fixture",
                "--species",
                "dog",
                "--cache-root",
                str(tmp_path),
                "--max-bytes",
                str(asset.size_bytes),
            ]
        )
        == 0
    )
    output = capsys.readouterr()
    assert output.out.strip() == str(path)
    assert not output.err
    assert path.read_bytes() == b"fixture"
    assert len(download_response.calls) == 1


def test_real_reference_manifests_are_distinct_and_complete():
    rows = oncoref.for_species("dog").references.registered_references()
    assert {r["id"] for r in rows} == {"refseq106", "ensembl98"}
    a, b = [ReferenceIdentity.from_dict(r["identity"]) for r in rows]
    assert a.reference_id != b.reference_id
    assert {a.assembly_accession, b.assembly_accession} == {"GCF_011100685.1", "GCA_000002285.2"}
    assert len(sa.species_assets("dog", group="references")) == 9
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
    original.unlink()
    for offline in (True, False):
        with pytest.raises(ValueError, match="missing"):
            sa.ensure_asset(asset, offline=offline, **options)


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
        assert options["raw"]
        assert options["resume"] == (os.name == "posix")
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

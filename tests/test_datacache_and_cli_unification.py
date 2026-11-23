"""Tests for unified datacache roots, OPENVAX_DATA_CACHE, install_archive store, and CLI (#747)."""

from __future__ import annotations

import hashlib
import json
import tarfile

import datacache
import pytest

from oncoref import (
    _downloads,
    cli,
    data_bundle,
    reference_data,
    sequence_bundle,
    source_matrices,
    species_assets,
)
from oncoref.version import DATA_VERSION, SOURCE_MATRIX_VERSION


def _clear_cache_env(monkeypatch):
    for key in (
        "OPENVAX_DATA_CACHE",
        "ONCOREF_CACHE_DIR",
        "CANCERDATA_DATA_DIR",
        "CANCERDATA_BUNDLED_DATA",
        "PIRLYGENES_BUNDLED_DATA",
        "CANCERDATA_SOURCE_MATRICES",
        "ONCOREF_SEQUENCE_CACHE_DIR",
        "ONCOREF_ASSET_CACHE",
        "ONCOREF_SHARD_CACHE_DIR",
        "OPENVAX_OFFLINE",
        "ONCOREF_OFFLINE",
    ):
        monkeypatch.delenv(key, raising=False)


def _write_bundle_fixture(root):
    for p in data_bundle.DOWNLOADABLE_PATHS:
        target = root / p
        if target.suffix:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(f"{p}\n")
        else:
            target.mkdir(parents=True, exist_ok=True)
            (target / "shard.parquet").write_text(f"{p}\n")


def _bundle_tarball_and_manifest(tmp_path):
    src = tmp_path / "src"
    _write_bundle_fixture(src)
    tar_path = tmp_path / "bundle.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tf:
        for p in data_bundle.DOWNLOADABLE_PATHS:
            tf.add(src / p, arcname=p)
    manifest = {
        "manifest_version": data_bundle.BUNDLE_MANIFEST_VERSION,
        "data_version": DATA_VERSION,
        "source": "oncoref",
        "repo": data_bundle.GITHUB_REPO,
        "manifest_url": data_bundle.RELEASE_MANIFEST_URL,
        "source_matrix_version": SOURCE_MATRIX_VERSION,
        "tarball": {
            "filename": data_bundle.TARBALL_FILENAME,
            "url": data_bundle.RELEASE_URL,
            "bytes": tar_path.stat().st_size,
            "sha256": hashlib.sha256(tar_path.read_bytes()).hexdigest(),
            "downloadable_paths": list(data_bundle.DOWNLOADABLE_PATHS),
        },
    }
    return tar_path, manifest


def test_openvax_data_cache_unifies_every_oncoref_cache_family(monkeypatch, tmp_path):
    _clear_cache_env(monkeypatch)
    shared = tmp_path / "openvax_shared"
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(shared))

    expected_base = shared / "oncoref"
    assert _downloads.oncoref_cache_root() == expected_base
    assert data_bundle.cache_root() == expected_base / "bundled_data"
    assert data_bundle.store_dir() == expected_base / "bundled_data" / f"v{DATA_VERSION}"
    assert reference_data.cache_base_dir() == expected_base
    assert source_matrices.cache_root() == expected_base / "source-matrices"
    assert sequence_bundle.cache_root() == expected_base / "sequence-bundles"
    assert species_assets.asset_cache_dir() == expected_base / "species"


def test_legacy_dot_cache_reused_when_platform_cache_empty(monkeypatch, tmp_path):
    _clear_cache_env(monkeypatch)
    fake_home = tmp_path / "home"
    legacy_root = fake_home / ".cache" / "oncoref"
    (legacy_root / "sources" / "hpa_rna_consensus" / "v23").mkdir(parents=True)
    (legacy_root / "sources" / "hpa_rna_consensus" / "v23" / "rna_tissue_consensus.tsv").write_text(
        "Gene\tTissue\tnTPM\n"
    )
    monkeypatch.setenv("HOME", str(fake_home))
    monkeypatch.setattr(
        datacache.common,
        "_platform_cache_dir",
        lambda name: str(fake_home / "Library" / "Caches" / name),
    )

    # Populated ~/.cache/oncoref is reused rather than forcing a re-download into ~/Library/Caches
    assert _downloads.oncoref_cache_root() == legacy_root

    # Once ~/Library/Caches/oncoref holds species data, split-brain subdirs in ~/.cache/oncoref
    # STILL resolve without forced re-downloads!
    plat_species = fake_home / "Library" / "Caches" / "oncoref" / "species" / "9615" / "marker.txt"
    plat_species.parent.mkdir(parents=True)
    plat_species.write_text("active")
    assert _downloads.oncoref_cache_root() == fake_home / "Library" / "Caches" / "oncoref"
    assert reference_data.cache_base_dir() == legacy_root
    (legacy_root / "bundled_data" / f"v{DATA_VERSION}").mkdir(parents=True)
    assert data_bundle.cache_root() == legacy_root / "bundled_data"


def test_data_bundle_uses_datacache_archive_store_and_prunes_generations(
    monkeypatch, tmp_path, http_downloads
):
    _clear_cache_env(monkeypatch)
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(tmp_path / "ovx"))
    tar_path, manifest = _bundle_tarball_and_manifest(tmp_path)
    manifest["tarball"]["url"] = "https://example.test/bundle.tar.gz"
    http_downloads(lambda _url: tar_path.read_bytes())

    sdir = data_bundle.store_dir()
    data_bundle._download_and_extract(
        "https://example.test/bundle.tar.gz",
        sdir,
        verbose=False,
        release_manifest=manifest,
    )

    # Installed as a real datacache archive store with atomic generation directory
    inspected = datacache.inspect_archive(sdir, verify_files=True)
    assert inspected.status == "available"
    assert inspected.recorded_sha256 == manifest["tarball"]["sha256"]
    reg = data_bundle.archive_registry(manifest)
    reg_inspected = reg.inspect("expression-bundle", verify_files=True)
    assert reg_inspected.status == "available"
    assert reg_inspected.verified is True
    assert len(datacache.list_bundles(sdir)) == 1
    assert (
        data_bundle.cache_dir()
        == tmp_path
        / "ovx"
        / "oncoref"
        / "bundled_data"
        / f"v{DATA_VERSION}"
        / "bundles"
        / datacache.list_bundles(sdir)[0].name
    )
    assert data_bundle.verify_local()["completion_marker"]["valid"] is True

    # Re-installing publishes a new generation and prunes the superseded one via prune_bundles
    data_bundle._download_and_extract(
        "https://example.test/bundle.tar.gz",
        sdir,
        verbose=False,
        release_manifest=manifest,
    )
    assert len(datacache.list_bundles(sdir)) == 1

    # Tampering with an installed file is caught by datacache.inspect_archive inside verify_local()
    target_csv = data_bundle.cache_dir() / "pan-cancer-expression.csv"
    target_csv.write_text("tampered\n")
    with pytest.raises(RuntimeError, match="datacache archive verification"):
        data_bundle.verify_local()


def test_streamlined_cli_actions_universal_json_and_global_offline(monkeypatch, tmp_path, capsys):
    _clear_cache_env(monkeypatch)
    cache_root = tmp_path / "cli_cache"

    # dirs --json --cache-root
    assert cli.main(["data", "dirs", "--cache-root", str(cache_root), "--json"]) == 0
    dirs = json.loads(capsys.readouterr().out)
    assert dirs["root"] == str(cache_root)
    assert dirs["bundle"].startswith(str(cache_root))
    assert dirs["species"] == str(cache_root / "species")

    # status --json
    assert cli.main(["data", "status", "hpa", "--cache-root", str(cache_root), "--json"]) == 0
    status_rows = json.loads(capsys.readouterr().out)
    assert isinstance(status_rows, list) and status_rows
    assert all(r["kind"] == "hpa" and r["present"] is False for r in status_rows)

    # info --json (bundle and wheel/reference datasets)
    assert cli.main(["data", "info", "bundle", "--cache-root", str(cache_root), "--json"]) == 0
    binfo = json.loads(capsys.readouterr().out)
    assert binfo["data_version"] == DATA_VERSION
    assert binfo["release_manifest"] is None

    assert cli.main(["data", "info", "housekeeping-genes", "--json"]) == 0
    hk_info = json.loads(capsys.readouterr().out)
    assert hk_info["name"] == "housekeeping-genes"
    assert hk_info["exists"] is True
    assert hk_info["available_versions"]

    # verify --json on wheel dataset succeeds; on missing bundle returns 1 with structured JSON
    assert cli.main(["data", "verify", "housekeeping-genes", "--json"]) == 0
    v_ok = json.loads(capsys.readouterr().out)
    assert v_ok["verified"] is True

    assert cli.main(["data", "verify", "bundle", "--cache-root", str(cache_root), "--json"]) == 1
    v_fail = json.loads(capsys.readouterr().out)
    assert v_fail["verified"] is False

    # Global --offline and OPENVAX_OFFLINE=1 block network fetches across human datasets too
    assert (
        cli.main(
            ["data", "fetch", "hpa_rna_consensus", "--cache-root", str(cache_root), "--offline"]
        )
        == 1
    )
    assert "Offline mode" in capsys.readouterr().err

    monkeypatch.setenv("OPENVAX_OFFLINE", "1")
    assert cli.main(["data", "path", "hpa_rna_consensus", "--cache-root", str(cache_root)]) == 1
    assert "Offline mode" in capsys.readouterr().err

    # Local wheel dataset path still resolves cleanly while OPENVAX_OFFLINE=1 and supports --json
    assert cli.main(["data", "path", "housekeeping-genes", "--json"]) == 0
    p_payload = json.loads(capsys.readouterr().out)
    assert p_payload["name"] == "housekeeping-genes"
    assert p_payload["path"].endswith(".csv")

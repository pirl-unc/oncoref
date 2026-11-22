# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

from __future__ import annotations

import json
import os

from filelock import FileLock

from oncoref import (
    cache_retention,
    cli,
    data_bundle,
    reference_data,
    sequence_bundle,
    source_matrices,
)
from oncoref.version import DATA_VERSION


def test_hardlink_inode_accounting_and_overlay_base_retention(monkeypatch, tmp_path):
    bundle_root = tmp_path / "bundled_data"
    base_ver = "5.23.0"
    obsolete_overlay_ver = "5.23.10"
    current_ver = DATA_VERSION

    base_dir = bundle_root / f"v{base_ver}"
    obs_dir = bundle_root / f"v{obsolete_overlay_ver}"
    cur_dir = bundle_root / f"v{current_ver}"
    for d in (base_dir, obs_dir, cur_dir):
        d.mkdir(parents=True)

    # Shared 10 KB base file hard-linked across all three bundle versions!
    shared_payload = b"S" * 10_000
    (base_dir / "pan-cancer-expression.csv").write_bytes(shared_payload)
    os.link(
        base_dir / "pan-cancer-expression.csv",
        obs_dir / "pan-cancer-expression.csv",
    )
    os.link(
        base_dir / "pan-cancer-expression.csv",
        cur_dir / "pan-cancer-expression.csv",
    )
    # Unique 2 KB delta in obsolete overlay and 3 KB delta in current overlay.
    (obs_dir / "delta.csv").write_bytes(b"O" * 2_000)
    (cur_dir / "delta.csv").write_bytes(b"C" * 3_000)

    # Current bundle is an overlay depending on v5.23.0 -> v5.23.0 MUST be retained!
    (cur_dir / data_bundle.CACHE_COMPLETE_FILENAME).write_text(
        json.dumps(
            {
                "manifest_version": 1,
                "data_version": current_ver,
                "release_manifest": {
                    "bundle_layout": "overlay",
                    "base_bundle": {"data_version": base_ver, "tarball": {"sha256": "a" * 64}},
                    "tarball": {"sha256": "b" * 64},
                },
            }
        )
    )
    monkeypatch.setenv("CANCERDATA_BUNDLED_DATA", str(cur_dir))

    inv = cache_retention.cache_inventory(targets="bundle")
    by_ver = {e["version"]: e for e in inv["entries"]}
    assert by_ver[f"v{current_ver}"]["status"] == "retained"
    assert by_ver[f"v{base_ver}"]["status"] == "retained"
    assert by_ver[f"v{base_ver}"]["reason"] == "overlay_base_dependency"
    assert by_ver[f"v{obsolete_overlay_ver}"]["prunable"] is True

    # Logical size counts the 10 KB shared file three times (~35 KB), while unique physical
    # bytes counts it once (~15 KB), and reclaimable_bytes is ONLY the 2 KB unique delta!
    assert inv["logical_bytes"] >= 35_000
    assert inv["unique_bytes"] < inv["logical_bytes"]
    assert by_ver[f"v{obsolete_overlay_ver}"]["reclaimable_bytes"] == 2_000
    assert inv["reclaimable_bytes"] == 2_000

    # Compatibility data_bundle.prune_cache also honors overlay base retention!
    compat_plan = data_bundle.prune_cache(dry_run=True)
    assert [e["version"] for e in compat_plan] == [f"v{obsolete_overlay_ver}"]

    res = cache_retention.prune_caches(targets="bundle", dry_run=False)
    assert res["deleted_count"] == 1
    assert not obs_dir.exists()
    assert base_dir.exists() and cur_dir.exists()
    assert (cur_dir / "pan-cancer-expression.csv").read_bytes() == shared_payload


def test_per_cohort_source_matrix_pins_and_newest_retention(monkeypatch, tmp_path):
    sm_root = tmp_path / "source-matrices"
    monkeypatch.setenv(source_matrices.CACHE_DIR_ENV_VAR, str(sm_root))

    luad_pin = source_matrices.source_matrix_version("LUAD")
    old_dir = sm_root / "v1.0.0"
    pin_dir = sm_root / f"v{luad_pin}"
    newer_dir = sm_root / "v99.0.0"
    for d in (old_dir, pin_dir, newer_dir):
        d.mkdir(parents=True)

    (old_dir / "LUAD.parquet").write_bytes(b"old")
    (pin_dir / "LUAD.parquet").write_bytes(b"pinned")
    (newer_dir / "LUAD.parquet").write_bytes(b"newer")
    # Unmanaged custom file must be protected, never deleted.
    (old_dir / "CUSTOM_COHORT.parquet").write_bytes(b"custom")

    inv = cache_retention.cache_inventory(targets="source")
    by_id = {e["id"]: e for e in inv["entries"]}
    assert by_id["LUAD@v1.0.0"]["prunable"] is True
    assert by_id[f"LUAD@v{luad_pin}"]["status"] == "retained"
    assert by_id["LUAD@v99.0.0"]["status"] == "retained"
    assert by_id["LUAD@v99.0.0"]["reason"] == "newest_cached_cohort_version"
    assert by_id["CUSTOM_COHORT@v1.0.0"]["status"] == "protected"

    pruned = cache_retention.prune_caches(targets="source", dry_run=False)
    assert pruned["deleted_count"] == 1
    assert not (old_dir / "LUAD.parquet").exists()
    assert (old_dir / "CUSTOM_COHORT.parquet").exists()
    assert (pin_dir / "LUAD.parquet").exists()
    assert (newer_dir / "LUAD.parquet").exists()


def test_reference_sequence_symlink_and_lock_protection_plus_cli(monkeypatch, tmp_path, capsys):
    data_root = tmp_path / "oncoref_data"
    seq_root = tmp_path / "seq_bundles"
    monkeypatch.setenv(reference_data.CACHE_DIR_ENV_KEY, str(data_root))
    monkeypatch.setenv("ONCOREF_SEQUENCE_CACHE_DIR", str(seq_root))
    monkeypatch.setenv("CANCERDATA_BUNDLED_DATA", str(tmp_path / "bundled" / f"v{DATA_VERSION}"))
    monkeypatch.setenv(source_matrices.CACHE_DIR_ENV_VAR, str(tmp_path / "empty_sm"))
    monkeypatch.setenv("ONCOREF_ASSET_CACHE", str(tmp_path / "empty_species"))

    # Reference: v23 is active pin; v19 is obsolete; symlinked v18 is protected.
    ref_active = data_root / "sources" / "hpa_rna_consensus" / "v23"
    ref_old = data_root / "sources" / "hpa_rna_consensus" / "v19"
    ref_active.mkdir(parents=True)
    ref_old.mkdir(parents=True)
    (ref_active / "rna_tissue_consensus.tsv").write_bytes(b"active")
    (ref_old / "rna_tissue_consensus.tsv").write_bytes(b"old_ref")

    outside_target = tmp_path / "outside_dir"
    outside_target.mkdir()
    (outside_target / "precious.txt").write_bytes(b"keep_me")
    (data_root / "sources" / "hpa_rna_consensus" / "v18").symlink_to(outside_target)

    # Sequence bundles: pinned digest kept; locked obsolete digest protected; unlocked obsolete pruned.
    pin_hex = sequence_bundle.descriptor()["bundle_id"].rsplit(":", 1)[-1]
    locked_hex = "1" * 64
    stale_hex = "2" * 64
    for h in (pin_hex, locked_hex, stale_hex):
        d = seq_root / h
        d.mkdir(parents=True)
        (d / "manifest.json").write_bytes(b"{}")

    held_lock = FileLock(str(seq_root / f"{locked_hex}.lock"))
    held_lock.acquire(timeout=1)
    try:
        # CLI dry-run JSON output
        rc = cli.main(["data", "prune", "all", "--json"])
        assert rc == 0
        payload = json.loads(capsys.readouterr().out)
        assert payload["dry_run"] is True
        statuses = {e["id"]: (e["status"], e["reason"]) for e in payload["entries"]}
        assert statuses["hpa_rna_consensus@v23"] == ("retained", "active_reference_pin")
        assert statuses["hpa_rna_consensus@v19"][0] == "prunable"
        assert statuses["hpa_rna_consensus@v18"] == ("protected", "symlink_protected")
        assert statuses[pin_hex] == ("retained", "active_sequence_bundle_pin")
        assert statuses[locked_hex] == ("protected", "active_lock")
        assert statuses[stale_hex][0] == "prunable"

        # CLI --yes deletes only unlocked, non-symlinked prunable items
        rc_del = cli.main(["data", "prune", "all", "--yes", "--json"])
        assert rc_del == 0
        del_payload = json.loads(capsys.readouterr().out)
        assert del_payload["deleted_count"] == 2
        assert not ref_old.exists()
        assert not (seq_root / stale_hex).exists()
        assert ref_active.exists()
        assert (outside_target / "precious.txt").exists()
        assert (seq_root / locked_hex).exists()
    finally:
        held_lock.release()

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Tests for unified reference-data management of externalized wheel tables (#738)."""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import pytest

from oncoref import (
    cache_retention,
    catalog,
    cta_rare_cell,
    cta_selection_audit,
    data_manifest,
    load_dataset,
    reference_data,
    translation_qualification,
)

_DATA_DIR = Path(__file__).resolve().parents[1] / "oncoref" / "data"


def test_externalized_originated_tables_are_managed_in_reference_data_and_catalog():
    inv = {row["name"]: row for row in catalog.inventory()}
    ds_names = {d.name for d in catalog.datasets()}

    assert len(data_manifest.EXTERNALIZED_ORIGINATED) >= 10
    for stem, (expected_cat, _) in data_manifest.EXTERNALIZED_ORIGINATED.items():
        assert stem not in data_manifest.CANCERDATA_ORIGINATED
        assert stem in reference_data.REFERENCE_SOURCES
        assert stem in ds_names
        assert stem in inv
        assert inv[stem]["held"] == "reference"
        assert inv[stem]["category"] == expected_cat

        spec = reference_data.REFERENCE_SOURCES[stem]
        ver = reference_data.resolve_version(stem)
        pin = spec["pins"][ver]
        local_file = _DATA_DIR / spec["filename"]
        assert local_file.is_file()
        assert local_file.stat().st_size == pin["bytes"]
        assert hashlib.sha256(local_file.read_bytes()).hexdigest() == pin["sha256"]
        assert pin["url"].startswith(
            "https://github.com/pirl-unc/oncoref/releases/download/reference-"
        )


def test_oncoref_cache_dir_env_fallback_and_precedence(monkeypatch, tmp_path):
    oncoref_dir = tmp_path / "oncoref_cache"
    cancerdata_dir = tmp_path / "cancerdata_cache"

    monkeypatch.delenv("CANCERDATA_DATA_DIR", raising=False)
    monkeypatch.setenv("ONCOREF_CACHE_DIR", str(oncoref_dir))
    assert reference_data.cache_base_dir() == oncoref_dir
    assert reference_data.cache_dir() == oncoref_dir / "sources"

    # Populate a cached active entry under ONCOREF_CACHE_DIR and verify cache_retention sees it
    ver = reference_data.resolve_version("cta-screen-results")
    dest = reference_data.local_path("cta-screen-results", ver)
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(_DATA_DIR / "cta-screen-results.csv", dest)

    report = cache_retention.cache_inventory(targets="reference")
    entries = {e["id"]: e for e in report["entries"]}
    assert f"cta-screen-results@{ver}" in entries
    assert entries[f"cta-screen-results@{ver}"]["status"] == "retained"

    # CANCERDATA_DATA_DIR takes precedence when both are set
    monkeypatch.setenv("CANCERDATA_DATA_DIR", str(cancerdata_dir))
    assert reference_data.cache_base_dir() == cancerdata_dir


def test_wheel_install_fallback_resolves_from_reference_cache_and_rejects_tampering(
    monkeypatch, tmp_path
):
    cache_root = tmp_path / "cache"
    empty_bundled = tmp_path / "empty_wheel_data"
    empty_bundled.mkdir()
    monkeypatch.delenv("CANCERDATA_DATA_DIR", raising=False)
    monkeypatch.setenv("ONCOREF_CACHE_DIR", str(cache_root))

    # Seed reference cache for target-gene-sets and cta-rare-cell-screen without network
    for name in ("target-gene-sets", "cta-rare-cell-screen", "cta-selection-audit-human"):
        spec = reference_data.REFERENCE_SOURCES[name]
        ver = reference_data.resolve_version(name)
        dest = reference_data.local_path(name, ver)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(_DATA_DIR / spec["filename"], dest)

    # Simulate an installed wheel whose bundled data dir omits externalized files
    resolved = reference_data.bundled_or_managed_path("target-gene-sets", bundled_dir=empty_bundled)
    assert resolved == reference_data.local_path("target-gene-sets")

    monkeypatch.setattr(load_dataset, "_BUNDLED_DATA_DIR", empty_bundled)
    load_dataset._clear_cache()
    try:
        df = load_dataset.get_data("target-gene-sets")
        assert not df.empty
        assert {"Gene ID", "Category"} <= set(df.columns)
    finally:
        monkeypatch.setattr(load_dataset, "_BUNDLED_DATA_DIR", _DATA_DIR)
        load_dataset._clear_cache()

    # Tampered local file in bundled_dir must fail closed on SHA-256 mismatch
    bad_bundled = tmp_path / "tampered_data"
    bad_bundled.mkdir()
    (bad_bundled / "target-gene-sets.csv.gz").write_bytes(b"corrupted-bytes")
    with pytest.raises(reference_data.ReferenceDataError, match="checksum/size"):
        reference_data.bundled_or_managed_path("target-gene-sets", bundled_dir=bad_bundled)

    # Direct readers work when bundled_or_managed_path resolves to cached reference files
    assert not cta_rare_cell.screen_table().empty
    assert cta_selection_audit.audit(species="human").species == "human"
    assert not translation_qualification.records().empty

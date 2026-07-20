# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

import csv
import gzip
import importlib.util
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).parents[1] / "scripts" / "canonicalize_reference_summary_sources.py"
_SPEC = importlib.util.spec_from_file_location("canonicalize_reference_summary_sources", _SCRIPT)
assert _SPEC and _SPEC.loader
migration = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(migration)


def _write_shard(path: Path, code: str, source: str) -> bytes:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Ensembl_Gene_ID",
                "cancer_code",
                "source_cohort",
                "source_version",
                "processing_pipeline",
                "notes",
                "TPM_median",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "Ensembl_Gene_ID": "ENSG000001",
                "cancer_code": code,
                "source_cohort": source,
                "source_version": (
                    "Treehouse Tumor Compendium 25.01 PolyA, TCGA subset; release 2025"
                ),
                "processing_pipeline": "treehouse_polya_25_01_tcga_subset_clean_tpm",
                "notes": "TCGA subset only: selected from the Treehouse compendium.",
                "TPM_median": "1.5",
            }
        )
    return path.read_bytes()


def _read_rows(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_canonicalize_sarc_histology_sources_rewrites_only_two_shards(tmp_path):
    summary = tmp_path / migration.SUMMARY_DIR
    for code in migration.SARC_HISTOLOGY_CODES:
        _write_shard(
            summary / f"{migration.LEGACY_SOURCE}__{code}.csv.gz",
            code,
            migration.LEGACY_SOURCE,
        )
    pleolps = summary / f"{migration.LEGACY_SOURCE}__SARC_PLEOLPS.csv.gz"
    pleolps_bytes = _write_shard(pleolps, "SARC_PLEOLPS", migration.LEGACY_SOURCE)

    result = migration.canonicalize_sarc_histology_sources(tmp_path)

    assert [row["cancer_code"] for row in result] == ["SARC_DDLPS", "SARC_WDLPS"]
    assert all(row["changed"] is True and row["rows"] == 1 for row in result)
    for code in migration.SARC_HISTOLOGY_CODES:
        assert not (summary / f"{migration.LEGACY_SOURCE}__{code}.csv.gz").exists()
        canonical = summary / f"{migration.CANONICAL_SOURCE}__{code}.csv.gz"
        assert _read_rows(canonical)[0]["source_cohort"] == migration.CANONICAL_SOURCE
    assert pleolps.read_bytes() == pleolps_bytes

    repeated = migration.canonicalize_sarc_histology_sources(tmp_path)
    assert all(row["changed"] is False for row in repeated)


def test_canonicalize_sarc_histology_sources_preserves_invalid_input(tmp_path):
    summary = tmp_path / migration.SUMMARY_DIR
    for code in migration.SARC_HISTOLOGY_CODES:
        _write_shard(
            summary / f"{migration.LEGACY_SOURCE}__{code}.csv.gz",
            code,
            migration.LEGACY_SOURCE,
        )
    invalid = summary / f"{migration.LEGACY_SOURCE}__SARC_DDLPS.csv.gz"
    original = _write_shard(invalid, "SARC_PLEOLPS", migration.LEGACY_SOURCE)

    with pytest.raises(ValueError, match="cancer_code='SARC_PLEOLPS'"):
        migration.canonicalize_sarc_histology_sources(tmp_path)

    assert invalid.read_bytes() == original
    assert not (summary / f"{migration.CANONICAL_SOURCE}__SARC_DDLPS.csv.gz").exists()


def test_canonicalize_reference_summary_sources_migrates_generic_and_sarc_destinations(
    tmp_path,
):
    summary = tmp_path / migration.SUMMARY_DIR
    for code in ("LUAD", "SARC_PLEOLPS", "SARC_DDLPS"):
        _write_shard(
            summary / f"{migration.LEGACY_SOURCE}__{code}.csv.gz",
            code,
            migration.LEGACY_SOURCE,
        )

    result = migration.canonicalize_reference_summary_sources(tmp_path)

    assert [row["cancer_code"] for row in result] == ["LUAD", "SARC_DDLPS", "SARC_PLEOLPS"]
    for code in ("LUAD", "SARC_PLEOLPS"):
        path = summary / f"{migration.TCGA_SAMPLES_SOURCE}__{code}.csv.gz"
        row = _read_rows(path)[0]
        assert row["source_cohort"] == migration.TCGA_SAMPLES_SOURCE
        assert row["source_version"] == (
            "Treehouse Tumor Compendium 25.01 PolyA samples selected by TCGA provenance; "
            "release 2025"
        )
        assert row["processing_pipeline"] == "treehouse_polya_25_01_tcga_samples_clean_tpm"
        assert row["notes"] == (
            "TCGA-provenance samples only: selected from the Treehouse compendium."
        )
    sarc = summary / f"{migration.SARC_HISTOLOGY_SOURCE}__SARC_DDLPS.csv.gz"
    assert _read_rows(sarc)[0]["source_cohort"] == migration.SARC_HISTOLOGY_SOURCE
    assert migration.canonicalize_reference_summary_sources(tmp_path) == []

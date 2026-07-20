# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Rewrite legacy Treehouse TCGA summary shards to canonical source identities.

This migration streams one gzip row at a time. It is intended for a staged data bundle,
not an installed cache that callers may be reading concurrently.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
from pathlib import Path

SUMMARY_DIR = "cancer-reference-expression"
LEGACY_SOURCE = "TREEHOUSE_POLYA_25_01_TCGA_SUBSET"
TCGA_SAMPLES_SOURCE = "TREEHOUSE_POLYA_25_01_TCGA_SAMPLES"
SARC_HISTOLOGY_SOURCE = "TREEHOUSE_POLYA_25_01_TCGA_SARC_HISTOLOGY"
# Backward-compatible name used by the original two-shard migration.
CANONICAL_SOURCE = SARC_HISTOLOGY_SOURCE
SARC_HISTOLOGY_CODES = ("SARC_DDLPS", "SARC_WDLPS")


def _shard_path(root: Path, source: str, code: str) -> Path:
    return root / SUMMARY_DIR / f"{source}__{code}.csv.gz"


def _validate_rows(path: Path, *, code: str, source: str) -> int:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"cancer_code", "source_cohort"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} lacks required columns: {', '.join(sorted(missing))}")
        n_rows = 0
        for row_number, row in enumerate(reader, start=2):
            if row["cancer_code"] != code:
                raise ValueError(
                    f"{path}:{row_number} has cancer_code={row['cancer_code']!r}; expected {code!r}"
                )
            if row["source_cohort"] != source:
                raise ValueError(
                    f"{path}:{row_number} has source_cohort={row['source_cohort']!r}; "
                    f"expected {source!r}"
                )
            n_rows += 1
    if n_rows == 0:
        raise ValueError(f"{path} contains no data rows")
    return n_rows


def _canonical_source(code: str) -> str:
    if code in SARC_HISTOLOGY_CODES:
        return SARC_HISTOLOGY_SOURCE
    return TCGA_SAMPLES_SOURCE


def _canonicalize_metadata(row: dict[str, str]) -> None:
    replacements = {
        "treehouse_polya_25_01_tcga_subset": "treehouse_polya_25_01_tcga_samples",
        "Treehouse Tumor Compendium 25.01 PolyA, TCGA subset": (
            "Treehouse Tumor Compendium 25.01 PolyA samples selected by TCGA provenance"
        ),
        "TCGA subset only": "TCGA-provenance samples only",
        "TCGA subset": "samples selected by TCGA provenance",
    }
    for column in ("source_version", "processing_pipeline", "notes"):
        value = row.get(column)
        if not value:
            continue
        for legacy, canonical in replacements.items():
            value = value.replace(legacy, canonical)
        row[column] = value


def _rewrite_shard(
    root: Path, code: str, *, canonical_source: str | None = None
) -> dict[str, object]:
    canonical_source = canonical_source or _canonical_source(code)
    legacy_path = _shard_path(root, LEGACY_SOURCE, code)
    canonical_path = _shard_path(root, canonical_source, code)
    if legacy_path.exists() and canonical_path.exists():
        raise FileExistsError(f"both legacy and canonical shards exist for {code}")
    if canonical_path.exists():
        return {
            "cancer_code": code,
            "rows": _validate_rows(canonical_path, code=code, source=canonical_source),
            "changed": False,
            "path": str(canonical_path),
        }
    if not legacy_path.exists():
        raise FileNotFoundError(f"missing legacy summary shard for {code}: {legacy_path}")

    temporary_path = canonical_path.with_name(f".{canonical_path.name}.tmp")
    try:
        with gzip.open(legacy_path, "rt", encoding="utf-8", newline="") as source_handle:
            reader = csv.DictReader(source_handle)
            required = {"cancer_code", "source_cohort"}
            missing = required - set(reader.fieldnames or ())
            if missing:
                raise ValueError(
                    f"{legacy_path} lacks required columns: {', '.join(sorted(missing))}"
                )
            with (
                temporary_path.open("wb") as raw_handle,
                gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as gzip_handle,
                io.TextIOWrapper(gzip_handle, encoding="utf-8", newline="") as target_handle,
            ):
                writer = csv.DictWriter(
                    target_handle,
                    fieldnames=reader.fieldnames,
                    lineterminator="\n",
                )
                writer.writeheader()
                n_rows = 0
                for row_number, row in enumerate(reader, start=2):
                    if row["cancer_code"] != code:
                        raise ValueError(
                            f"{legacy_path}:{row_number} has "
                            f"cancer_code={row['cancer_code']!r}; expected {code!r}"
                        )
                    if row["source_cohort"] != LEGACY_SOURCE:
                        raise ValueError(
                            f"{legacy_path}:{row_number} has "
                            f"source_cohort={row['source_cohort']!r}; "
                            f"expected {LEGACY_SOURCE!r}"
                        )
                    row["source_cohort"] = canonical_source
                    _canonicalize_metadata(row)
                    writer.writerow(row)
                    n_rows += 1
        if n_rows == 0:
            raise ValueError(f"{legacy_path} contains no data rows")
        temporary_path.replace(canonical_path)
        legacy_path.unlink()
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise

    return {
        "cancer_code": code,
        "rows": n_rows,
        "changed": True,
        "path": str(canonical_path),
    }


def canonicalize_sarc_histology_sources(root: Path) -> list[dict[str, object]]:
    """Canonicalize the DDLPS/WDLPS physical shards under a staged bundle root."""
    return [
        _rewrite_shard(root, code, canonical_source=SARC_HISTOLOGY_SOURCE)
        for code in SARC_HISTOLOGY_CODES
    ]


def canonicalize_reference_summary_sources(root: Path) -> list[dict[str, object]]:
    """Canonicalize every physical shard that still uses the legacy generic ID."""
    summary_dir = root / SUMMARY_DIR
    prefix = f"{LEGACY_SOURCE}__"
    suffix = ".csv.gz"
    legacy_paths = sorted(summary_dir.glob(f"{LEGACY_SOURCE}__*{suffix}"))
    codes = [path.name[len(prefix) : -len(suffix)] for path in legacy_paths]
    return [_rewrite_shard(root, code) for code in codes]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle_root", type=Path, help="staged data-bundle directory")
    args = parser.parse_args()
    print(json.dumps(canonicalize_reference_summary_sources(args.bundle_root), indent=2))


if __name__ == "__main__":
    main()

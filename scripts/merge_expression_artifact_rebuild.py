# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Merge a focused expression-artifact rebuild into a complete bundle directory.

The focused rebuild must come from ``rebuild_expression_artifacts.py``. Rows and
shards for its cancer codes replace any existing versions in the bundle; all
other cohorts remain unchanged. Bundle-level sample totals are recomputed from
the merged per-cohort metadata.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from oncoref.data_bundle import DOWNLOADABLE_PATHS
from oncoref.expression import (
    EXPRESSION_ARTIFACT_BUILD_METADATA_JSON_PATH,
    EXPRESSION_ARTIFACT_BUILD_METADATA_PATH,
    EXPRESSION_ARTIFACT_BUILD_METADATA_SCHEMA_VERSION,
    SOURCE_MATRIX_SAMPLE_QC_MANIFEST_PATH,
)

_REPRESENTATIVE_PROVENANCE = "cancer-reference-expression-representatives/_provenance.csv"
_REFERENCE_SUMMARY_DIR = "cancer-reference-expression"
_FOCUSED_SHARD_DIRS = tuple(
    path for path in DOWNLOADABLE_PATHS if path.startswith("cancer-reference-expression-")
)


def _replace_rows(
    bundle_path: Path,
    rebuild_path: Path,
    *,
    cancer_codes: set[str],
) -> pd.DataFrame:
    existing = pd.read_csv(bundle_path)
    rebuilt = pd.read_csv(rebuild_path)
    kept = existing[~existing["cancer_code"].astype(str).isin(cancer_codes)]
    merged = pd.concat([kept, rebuilt], ignore_index=True, sort=False)
    merged.to_csv(bundle_path, index=False)
    return merged


def _merge_representative_provenance(
    bundle_dir: Path,
    rebuild_dir: Path,
    *,
    cancer_codes: set[str],
) -> None:
    bundle_path = bundle_dir / _REPRESENTATIVE_PROVENANCE
    rebuild_path = rebuild_dir / _REPRESENTATIVE_PROVENANCE
    existing = pd.read_csv(bundle_path)
    rebuilt = pd.read_csv(rebuild_path)
    existing_codes = existing["representative_id"].astype(str).str.split("__").str[0]
    kept = existing[~existing_codes.isin(cancer_codes)]
    pd.concat([kept, rebuilt], ignore_index=True, sort=False).to_csv(bundle_path, index=False)


def _copy_rebuilt_shards(
    bundle_dir: Path,
    rebuild_dir: Path,
    *,
    cancer_codes: set[str],
) -> None:
    copies: list[tuple[Path, Path]] = []
    for relative in _FOCUSED_SHARD_DIRS:
        rebuilt_subdir = rebuild_dir / relative
        if not rebuilt_subdir.is_dir():
            raise FileNotFoundError(f"focused rebuild lacks artifact directory: {relative}")
        bundle_subdir = bundle_dir / relative
        for code in sorted(cancer_codes):
            source = rebuilt_subdir / f"{code}.parquet"
            if not source.exists():
                raise FileNotFoundError(f"focused rebuild lacks {relative}/{code}.parquet")
            copies.append((source, bundle_subdir / source.name))

    for source, destination in copies:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def _code_specific_summary_path(summary_dir: Path, source_stem: str, code: str) -> Path | None:
    candidates = (
        summary_dir / f"{source_stem}__{code}.csv.gz",
        summary_dir / f"{source_stem}__{code}.csv",
    )
    return next((path for path in candidates if path.exists()), None)


def _without_summary_codes(frame: pd.DataFrame, *, codes: set[str]) -> pd.DataFrame:
    if frame.empty:
        return frame
    return frame[~frame["cancer_code"].astype(str).isin(codes)]


def _replace_summary_codes(
    existing: pd.DataFrame,
    rebuilt: pd.DataFrame,
    *,
    codes: set[str],
) -> pd.DataFrame:
    """Replace ``codes`` in one summary frame, preserving all other rows."""
    kept = _without_summary_codes(existing, codes=codes)
    replacements = rebuilt[rebuilt["cancer_code"].astype(str).isin(codes)]
    return pd.concat([kept, replacements], ignore_index=True, sort=False)


def _merge_reference_summaries(
    bundle_dir: Path,
    rebuild_dir: Path,
    *,
    cancer_codes: set[str],
) -> None:
    """Replace focused cancer-code rows while preserving other rows in each source shard."""
    rebuilt_dir = rebuild_dir / _REFERENCE_SUMMARY_DIR
    rebuilt_paths = sorted(rebuilt_dir.glob("*.csv")) + sorted(rebuilt_dir.glob("*.csv.gz"))
    if not rebuilt_paths:
        raise FileNotFoundError(
            f"focused rebuild lacks reference summary shards under {_REFERENCE_SUMMARY_DIR}"
        )

    seen_codes: set[str] = set()
    writes: list[tuple[Path, pd.DataFrame]] = []
    deletes: list[Path] = []
    for source in rebuilt_paths:
        rebuilt = pd.read_csv(source)
        shard_codes = set(rebuilt["cancer_code"].astype(str))
        unexpected = sorted(shard_codes - cancer_codes)
        if unexpected:
            raise ValueError(
                f"focused summary shard {source.name} contains unexpected codes: {unexpected}"
            )
        duplicates = sorted(shard_codes & seen_codes)
        if duplicates:
            raise ValueError(f"focused summary rows occur in multiple shards: {duplicates}")
        seen_codes.update(shard_codes)

        bundle_summary_dir = bundle_dir / _REFERENCE_SUMMARY_DIR
        source_stem = source.name.removesuffix(".gz").removesuffix(".csv")
        code_specific: dict[str, Path] = {}
        for code in sorted(shard_codes):
            destination = _code_specific_summary_path(bundle_summary_dir, source_stem, code)
            if destination is not None:
                code_specific[code] = destination

        for code, destination in code_specific.items():
            existing = pd.read_csv(destination)
            writes.append((destination, _replace_summary_codes(existing, rebuilt, codes={code})))

        consolidated = bundle_summary_dir / source.name
        consolidated_codes = shard_codes - set(code_specific)
        if consolidated.exists() or consolidated_codes:
            existing = pd.read_csv(consolidated) if consolidated.exists() else pd.DataFrame()
            kept = _without_summary_codes(existing, codes=shard_codes)
            replacements = rebuilt[rebuilt["cancer_code"].astype(str).isin(consolidated_codes)]
            merged = pd.concat(
                [kept, replacements],
                ignore_index=True,
                sort=False,
            )
            if merged.empty:
                deletes.append(consolidated)
            else:
                writes.append((consolidated, merged))

    missing = sorted(cancer_codes - seen_codes)
    if missing:
        raise ValueError(f"focused rebuild lacks reference summary rows for: {missing}")
    for destination in deletes:
        destination.unlink(missing_ok=True)
    for destination, rows in writes:
        destination.parent.mkdir(parents=True, exist_ok=True)
        rows.to_csv(destination, index=False, compression="infer")


def _metadata_sum(metadata: pd.DataFrame, column: str) -> int:
    if column not in metadata:
        return 0
    return int(pd.to_numeric(metadata[column], errors="coerce").fillna(0).sum())


def merge(bundle_dir: Path, rebuild_dir: Path) -> set[str]:
    """Merge ``rebuild_dir`` into ``bundle_dir`` and return replaced cancer codes."""
    rebuilt_metadata_path = rebuild_dir / EXPRESSION_ARTIFACT_BUILD_METADATA_PATH
    rebuilt_metadata = pd.read_csv(rebuilt_metadata_path)
    cancer_codes = set(rebuilt_metadata["cancer_code"].astype(str))
    if not cancer_codes:
        raise ValueError("focused rebuild metadata contains no cancer codes")

    _copy_rebuilt_shards(bundle_dir, rebuild_dir, cancer_codes=cancer_codes)
    _merge_reference_summaries(bundle_dir, rebuild_dir, cancer_codes=cancer_codes)
    _merge_representative_provenance(
        bundle_dir,
        rebuild_dir,
        cancer_codes=cancer_codes,
    )
    _replace_rows(
        bundle_dir / SOURCE_MATRIX_SAMPLE_QC_MANIFEST_PATH,
        rebuild_dir / SOURCE_MATRIX_SAMPLE_QC_MANIFEST_PATH,
        cancer_codes=cancer_codes,
    )
    metadata = _replace_rows(
        bundle_dir / EXPRESSION_ARTIFACT_BUILD_METADATA_PATH,
        rebuilt_metadata_path,
        cancer_codes=cancer_codes,
    )
    if "build_source_cohort" not in metadata:
        metadata["build_source_cohort"] = metadata["source_cohort"]
    else:
        build_source = metadata["build_source_cohort"].astype("string")
        missing_build_source = build_source.isna() | build_source.str.strip().eq("")
        metadata.loc[missing_build_source, "build_source_cohort"] = metadata.loc[
            missing_build_source, "source_cohort"
        ]
    metadata.to_csv(bundle_dir / EXPRESSION_ARTIFACT_BUILD_METADATA_PATH, index=False)

    summary_path = bundle_dir / EXPRESSION_ARTIFACT_BUILD_METADATA_JSON_PATH
    summary = json.loads(summary_path.read_text())
    fallback = metadata.get("sample_qc_fallback_reason", pd.Series(dtype="string"))
    derived_artifacts = list(summary.get("derived_artifacts") or [])
    if _REFERENCE_SUMMARY_DIR not in derived_artifacts:
        derived_artifacts.append(_REFERENCE_SUMMARY_DIR)
    summary.update(
        {
            "schema_version": EXPRESSION_ARTIFACT_BUILD_METADATA_SCHEMA_VERSION,
            "n_cohorts": len(metadata),
            "n_source_samples": _metadata_sum(metadata, "n_source_samples"),
            "n_cohort_samples": _metadata_sum(metadata, "n_cohort_samples"),
            "n_negative_values_clipped": _metadata_sum(metadata, "n_negative_values_clipped"),
            "sample_qc_fallbacks": int(fallback.fillna("").astype(str).str.strip().ne("").sum()),
            "derived_artifacts": derived_artifacts,
        }
    )
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return cancer_codes


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle_dir", type=Path, help="Complete bundle staging directory")
    parser.add_argument("rebuild_dir", type=Path, help="Focused rebuild output directory")
    args = parser.parse_args(argv)
    codes = merge(args.bundle_dir.expanduser(), args.rebuild_dir.expanduser())
    print(f"merged {len(codes)} cohort(s): {', '.join(sorted(codes))}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

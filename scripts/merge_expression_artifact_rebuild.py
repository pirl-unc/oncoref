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
    for relative in DOWNLOADABLE_PATHS:
        rebuilt_subdir = rebuild_dir / relative
        if not rebuilt_subdir.is_dir():
            continue
        bundle_subdir = bundle_dir / relative
        bundle_subdir.mkdir(parents=True, exist_ok=True)
        for code in sorted(cancer_codes):
            source = rebuilt_subdir / f"{code}.parquet"
            if not source.exists():
                raise FileNotFoundError(f"focused rebuild lacks {relative}/{code}.parquet")
            shutil.copy2(source, bundle_subdir / source.name)


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
    summary.update(
        {
            "schema_version": EXPRESSION_ARTIFACT_BUILD_METADATA_SCHEMA_VERSION,
            "n_cohorts": len(metadata),
            "n_source_samples": _metadata_sum(metadata, "n_source_samples"),
            "n_cohort_samples": _metadata_sum(metadata, "n_cohort_samples"),
            "n_negative_values_clipped": _metadata_sum(metadata, "n_negative_values_clipped"),
            "sample_qc_fallbacks": int(fallback.fillna("").astype(str).str.strip().ne("").sum()),
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

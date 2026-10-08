#!/usr/bin/env python
"""Refresh CTA protein summaries while preserving a published bundle's gene basis.

Read the selected raw source matrices and the base bundle's sample-QC manifest.
Before changing protein groups, require the rebuilt gene-level percentiles and
prevalence to reproduce the base artifacts. Output contains only the two refreshed
proteoform directories and updated bundle metadata; copy these into a staged full
bundle before running build_data_overlay.py. Source matrices must be cached (use
CANCERDATA_SOURCE_MATRICES to select a separate staging cache).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.rebuild_expression_artifacts import (
    _clip_negative_expression,
    _drop_technical,
    _sample_selection_for_qc,
    build_clean,
    read_raw,
)

from oncoref import source_matrices
from oncoref.expression import SHARD_DATASETS
from oncoref.expression_builders import (
    cohort_percentile_vectors,
    sample_columns,
    within_sample_top_fractions,
)
from oncoref.proteoforms import _proteoform_registry_metadata, collapse_to_proteoforms


def refresh(base: Path, out: Path) -> None:
    """Rebuild every base cohort; reject source, normalization or sample-QC drift."""
    if base.resolve() == out.resolve():
        raise ValueError("output must differ from the published base bundle")
    metadata = pd.read_csv(base / "expression-artifact-build-metadata.csv")
    qc = pd.read_csv(base / "source-matrix-sample-qc.csv")
    expected_codes = set(metadata["cancer_code"].astype(str))
    builds = {
        "percentiles": cohort_percentile_vectors,
        "within_sample": within_sample_top_fractions,
    }
    for name in builds:
        dataset = SHARD_DATASETS[name]
        for proteoform in (False, True):
            codes = {
                p.stem for p in (base / dataset.subdir(proteoform=proteoform)).glob("*.parquet")
            }
            if codes != expected_codes:
                raise ValueError(f"{name}: base shard inventory differs from cohort metadata")
        (out / dataset.subdir(proteoform=True)).mkdir(parents=True, exist_ok=True)

    audit = []
    for row in metadata.itertuples(index=False):
        code = str(row.cancer_code)
        source = source_matrices.local_path(code)
        if not source.is_file():
            raise FileNotFoundError(f"{code}: cache the published source matrix at {source}")
        raw = read_raw(source)
        raw_samples = sample_columns(raw)
        raw, n_clipped = _clip_negative_expression(raw, raw_samples)
        cohort_qc = qc.loc[qc["cancer_code"].astype(str) == code]
        samples, effective, fallback = _sample_selection_for_qc(
            raw_samples, cohort_qc, row.sample_qc
        )
        if len(raw_samples) != row.n_source_samples or len(samples) != row.n_cohort_samples:
            raise ValueError(f"{code}: source or QC-selected sample count drift")
        if effective != row.sample_qc_effective or fallback != (
            row.sample_qc_fallback_reason if pd.notna(row.sample_qc_fallback_reason) else ""
        ):
            raise ValueError(f"{code}: sample-QC policy drift")
        clean = build_clean(raw)[["Ensembl_Gene_ID", "Symbol", *samples]]
        bio = _drop_technical(clean)
        grouped = collapse_to_proteoforms(bio, scope="cta", sample_cols=samples)
        for name, build in builds.items():
            dataset = SHARD_DATASETS[name]
            expected = pd.read_parquet(base / dataset.gene_dir / f"{code}.parquet")
            actual = build(bio, samples)
            pd.testing.assert_frame_equal(actual, expected, check_dtype=False, check_exact=True)
            result = build(grouped, samples)
            destination = out / dataset.subdir(proteoform=True) / f"{code}.parquet"
            result.to_parquet(destination, index=False, compression="zstd")
            restored = pd.read_parquet(destination)
            pd.testing.assert_frame_equal(result, restored)
            if result.attrs != restored.attrs:
                raise ValueError(f"{code}: proteoform registry metadata did not round-trip")
        audit.append(
            {
                "cancer_code": code,
                "source_matrix": str(source),
                "n_samples": len(samples),
                "n_genes": len(bio),
                "n_proteoforms": len(grouped),
                "n_negative_values_clipped": n_clipped,
            }
        )
        print(
            f"{code}: gene artifacts unchanged; {len(samples)} samples, {len(grouped)} protein rows",
            flush=True,
        )

    bundle_metadata = json.loads((base / "expression-artifact-build-metadata.json").read_text())
    bundle_metadata["proteoform_registry"] = _proteoform_registry_metadata("cta")
    (out / "expression-artifact-build-metadata.json").write_text(
        json.dumps(bundle_metadata, indent=2, sort_keys=True) + "\n"
    )
    (out / "proteoform-refresh-audit.json").write_text(json.dumps(audit, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    refresh(args.base, args.out)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Re-extract the named-candidate audit from its pinned local source matrices.

Set CANCERDATA_SOURCE_MATRICES to the cache containing v<version>/<CODE>.parquet.
This command never downloads or edits the matrices. Invalid matrices are logged,
not repaired. Build the report separately after reviewing changes to its inputs.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cta_report_common import sha256, write_csv  # noqa: E402

from oncoref import expression, source_matrices  # noqa: E402

DEFAULT = ROOT / "docs/audits/cta-candidate-evidence-20261108/inputs"
QC_COLUMNS = [
    "cancer_code",
    "source_cohort",
    "sample_id",
    "source_scale_class",
    "linear_tpm_comparable",
    "tpm_proxy",
    "sample_qc_status",
    "sample_qc_reasons",
    "recommended_for_absolute_tpm_floor",
    "total_tpm",
    "n_measured_genes",
    "n_detected_genes",
    "top_gene_fraction",
    "top10_fraction",
    "housekeeping_genes_above_floor",
]


def collect(inputs):
    records = json.loads((inputs / "candidates.json").read_text())
    ids = [r["gene_id"] for r in records]
    provenance = json.loads((inputs / "matrix-provenance.json").read_text())
    pinned = {r["cohort"]: r for r in provenance}
    registry = source_matrices.registry()
    if set(pinned) != set(registry.cancer_code):
        raise ValueError("Source registry changed; review the audit scope before recollecting")
    parts, coverage, quality, errors = [], [], [], []
    for r in registry.itertuples(index=False):
        code = r.cancer_code
        asset = f"source-v{r.source_matrix_version}/{code}.parquet"
        path = source_matrices.local_path(code)
        if pinned[code]["asset"] != asset or sha256(path) != pinned[code]["sha256"]:
            raise ValueError(f"Pinned source identity mismatch: {code}")
        data = pd.read_parquet(path, filters=[("Ensembl_Gene_ID", "in", ids)])
        for gid in ids:
            coverage.append(
                {
                    "cohort": code,
                    "gene_id": gid,
                    "status": "measured" if data.Ensembl_Gene_ID.eq(gid).any() else "gene_absent",
                    "matrix_version": r.source_matrix_version,
                    "matrix_asset": asset,
                }
            )
        rows = data.melt(id_vars=["Ensembl_Gene_ID", "Symbol"], var_name="sample", value_name="tpm")
        rows["cohort"] = code
        rows["source_cohort"] = r.source_cohort
        rows["namespace"] = source_matrices.source_sample_namespace(r.source_cohort)
        rows["matrix_version"] = r.source_matrix_version
        parts.append(rows)
        try:
            quality.append(source_matrices.sample_qc(code, auto_fetch=False)[QC_COLUMNS])
        except ValueError as error:
            # The raw values stay in the ledger; none is eligible without QC.
            errors.append({"cohort": code, "reason": str(error).replace(str(path), asset)})
        finally:
            expression._load_per_sample_matrix.cache_clear()
        print(code, flush=True)
    write_csv(pd.concat(parts, ignore_index=True), inputs / "tumor-samples.csv.gz")
    write_csv(pd.DataFrame(coverage), inputs / "tumor-coverage.csv")
    write_csv(pd.concat(quality, ignore_index=True), inputs / "sample-qc.csv.gz")
    (inputs / "sample-qc-errors.json").write_text(json.dumps(errors, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=DEFAULT)
    collect(parser.parse_args().inputs)

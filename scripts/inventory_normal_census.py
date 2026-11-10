#!/usr/bin/env python3
"""Cache primary, normal-labelled human tissue metadata from a pinned Census.

This is an inventory, not a declaration that every retained cell passes QC.
It preserves age, assay, donor, cell type, tissue and study identity for review.
Requires the optional cellxgene-census analysis environment.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import cellxgene_census
import pandas as pd
import pyarrow.parquet as pq

FILTER = "is_primary_data == True and disease == 'normal' and tissue_type == 'tissue'"
COLUMNS = [
    "soma_joinid",
    "dataset_id",
    "donor_id",
    "assay",
    "cell_type",
    "cell_type_ontology_term_id",
    "development_stage",
    "development_stage_ontology_term_id",
    "sex",
    "suspension_type",
    "tissue",
    "tissue_general",
    "raw_sum",
    "nnz",
    "n_measured_vars",
]
GROUPS = [
    "dataset_id",
    "donor_id",
    "assay",
    "cell_type",
    "cell_type_ontology_term_id",
    "development_stage",
    "development_stage_ontology_term_id",
    "sex",
    "suspension_type",
    "tissue",
    "tissue_general",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", required=True, help="Explicit build date, never latest/stable")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.release in {"stable", "latest"}:
        parser.error("Pin an explicit release date")
    args.output.mkdir(parents=True, exist_ok=True)
    counts = []
    files = []
    with cellxgene_census.open_soma(census_version=args.release) as census:
        datasets = census["census_info"]["datasets"].read().concat().to_pandas()
        datasets.to_csv(args.output / "datasets.csv", index=False)
        census["census_data"]["homo_sapiens"].ms["RNA"].var.read().concat().to_pandas().to_parquet(
            args.output / "features.parquet", index=False
        )
        read = census["census_data"]["homo_sapiens"].obs.read(
            value_filter=FILTER, column_names=COLUMNS
        )
        for i, table in enumerate(read):
            path = args.output / f"observations-{i:04d}.parquet"
            pq.write_table(table, path, compression="zstd")
            frame = table.to_pandas()
            count = (
                frame.groupby(GROUPS, observed=True, dropna=False)
                .size()
                .rename("n_cells")
                .reset_index()
            )
            counts.append(count)
            files.append(
                {
                    "path": path.name,
                    "n_rows": len(frame),
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                }
            )
            print(f"metadata chunk {i}: {len(frame):,} cells", flush=True)
    summary = (
        pd.concat(counts).groupby(GROUPS, observed=True, dropna=False).n_cells.sum().reset_index()
    )
    summary.to_csv(
        args.output / "populations.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )
    for field in ["development_stage", "tissue_general", "assay"]:
        summary.groupby(field, observed=True).n_cells.sum().sort_values(ascending=False).to_csv(
            args.output / f"{field}-counts.csv"
        )
    manifest = {
        "census_release": args.release,
        "filter": FILTER,
        "normal_label_is_not_QC": True,
        "n_cells": int(summary.n_cells.sum()),
        "n_datasets": int(summary.dataset_id.nunique()),
        "files": files,
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: v for k, v in manifest.items() if k != "files"}, indent=2), flush=True)


if __name__ == "__main__":
    main()

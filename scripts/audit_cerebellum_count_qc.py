#!/usr/bin/env python3
"""Inspect count support and sequencing-depth dependence of sparse cerebellum hits.

This diagnoses sensitivity to counting rules; it does not classify RNA as ambient
or establish protein production. A Poisson comparison is descriptive only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from audit_cta_cerebellum import extract, selections, sha256
from audit_plac1_brain import read_column
from scipy import sparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    proteins = pd.read_csv(args.comparison / "adult-granule-all-proteins.csv")
    requested = {g for ids in proteins.member_gene_ids for g in json.loads(ids)}
    cells, genes, matrix = extract(args.h5ad, requested)
    with h5py.File(args.h5ad) as handle:
        for name in ["batch", "assay", "size_factor"]:
            cells[name] = read_column(handle, f"obs/{name}")
        pointers = handle["X/indptr"][:]
        total, detected = [], []
        for start in range(0, len(cells), 10000):
            end = min(start + 10000, len(cells))
            left, right = int(pointers[start]), int(pointers[end])
            block = sparse.csr_matrix(
                (
                    handle["X/data"][left:right],
                    handle["X/indices"][left:right],
                    pointers[start : end + 1] - left,
                ),
                shape=(end - start, len(handle["var/_index"])),
            )
            block.sum_duplicates()
            block.eliminate_zeros()
            total.extend(np.asarray(block.sum(axis=1)).ravel())
            detected.extend(np.diff(block.indptr))
    cells["total_counts"] = total
    cells["detected_genes"] = detected
    mask = selections(cells)["adult_granule"].to_numpy()
    obs = cells.loc[mask].reset_index(drop=True)
    mat = matrix[mask].tocsc()
    obs["depth_quartile"] = (
        obs.groupby(["donor_id", "batch"], observed=True)
        .total_counts.transform(lambda x: pd.qcut(x.rank(method="first"), 4, labels=[1, 2, 3, 4]))
        .astype(int)
    )
    summaries, by_depth = [], []
    for j, gene in enumerate(genes.itertuples()):
        values = mat[:, j].toarray().ravel()
        for scope, indices in [
            ("all", np.arange(len(obs))),
            *list(obs.groupby("donor_id").indices.items()),
        ]:
            v, depth = values[indices], obs.iloc[indices].total_counts.to_numpy()
            mean = float(v.mean())
            summaries.append(
                {
                    "gene_id": gene.gene_id,
                    "symbol": gene.dataset_symbol,
                    "scope": scope,
                    "n_cells": len(v),
                    "n_positive": int((v > 0).sum()),
                    "n_single": int((v == 1).sum()),
                    "n_ge2": int((v >= 2).sum()),
                    "n_ge3": int((v >= 3).sum()),
                    "max_count": int(v.max()),
                    "total_counts": int(v.sum()),
                    "mean_count": mean,
                    "pct_positive": float(100 * (v > 0).mean()),
                    "pct_ge2": float(100 * (v >= 2).mean()),
                    "constant_poisson_expected_ge2": float(
                        len(v) * (-np.expm1(-mean) - mean * np.exp(-mean))
                    ),
                    "median_depth_positive": float(np.median(depth[v > 0]))
                    if (v > 0).any()
                    else None,
                    "median_depth_negative": float(np.median(depth[v == 0]))
                    if (v == 0).any()
                    else None,
                }
            )
        for (donor, batch, quartile), idx in obs.groupby(
            ["donor_id", "batch", "depth_quartile"]
        ).indices.items():
            v = values[idx]
            by_depth.append(
                {
                    "gene_id": gene.gene_id,
                    "symbol": gene.dataset_symbol,
                    "donor_id": donor,
                    "batch": batch,
                    "depth_quartile": quartile,
                    "n_cells": len(idx),
                    "n_positive": int((v > 0).sum()),
                    "pct_positive": float(100 * (v > 0).mean()),
                    "median_total_counts": float(obs.iloc[idx].total_counts.median()),
                }
            )
    pd.DataFrame(summaries).to_csv(args.output / "count-support.csv", index=False)
    pd.DataFrame(by_depth).to_csv(
        args.output / "depth-sensitivity.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    obs.groupby(["donor_id", "batch", "assay"]).agg(
        nuclei=("total_counts", "size"),
        median_counts=("total_counts", "median"),
        median_genes=("detected_genes", "median"),
    ).to_csv(args.output / "libraries.csv")
    (args.output / "manifest.json").write_text(
        json.dumps(
            {
                "source_sha256": sha256(args.h5ad),
                "n_adult_normal_granule": len(obs),
                "gene_mapping": "exact Ensembl IDs",
                "classification_changed": False,
                "interpretation": "Depth dependence and Poisson compatibility do not distinguish genuine weak RNA from background. No alignment or empty-droplet data were examined.",
            },
            indent=2,
        )
        + "\n"
    )
    print(
        pd.DataFrame(summaries)
        .query(
            "scope == 'all' and symbol in ['PLAC1','CBY3','CAGE1','DMRT1','NUTM1','TMEM269','SYCP1','PRAME']"
        )
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()

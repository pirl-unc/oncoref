#!/usr/bin/env python3
"""Stream raw counts for the rare-cell panel only, by cell, from a pinned Census.

``audit_normal_census.py`` scans every normal observation once per gene block.
This reads only the adult, donor-resolved panel cells selected by
``build_rare_cell_evidence.panel`` and all requested genes in one pass. Output
rows have the same columns as the block files. Resumable by cell chunk.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import cellxgene_census as cc
import numpy as np
import pandas as pd


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def summarize(tables, mapping, var_map, n_populations, n_genes):
    shape = (n_populations, n_genes)
    n_positive = np.zeros(shape, np.int32)
    n_ge2 = np.zeros(shape, np.int32)
    n_fractional = np.zeros(shape, np.int32)
    count_sum = np.zeros(shape, np.float64)
    count_max = np.zeros(shape, np.float32)
    for table in tables:
        rows = mapping[table["soma_dim_0"].to_numpy()]
        cols = var_map[table["soma_dim_1"].to_numpy()]
        values = table["soma_data"].to_numpy()
        if (rows < 0).any() or (cols < 0).any():
            raise ValueError("Expression outside the requested panel")
        if not np.all(np.isfinite(values) & (values > 0)):
            raise ValueError("Invalid sparse counts")
        idx = (rows, cols)
        np.add.at(n_positive, idx, 1)
        np.add.at(n_ge2, idx, values >= 2)
        np.add.at(n_fractional, idx, values != np.floor(values))
        np.add.at(count_sum, idx, values)
        np.maximum.at(count_max, idx, values)
    return n_positive, n_ge2, n_fractional, count_sum, count_max


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--chunk", type=int, default=100_000)
    parser.add_argument("--shard", default="0/1", help="K/N: take every Nth chunk from K")
    parser.add_argument("--max-chunks", type=int)
    args = parser.parse_args()
    k, n = map(int, args.shard.split("/"))
    out = args.counts / "panel-counts"
    out.mkdir(exist_ok=True)
    base = json.loads((args.counts / "manifest.json").read_text())
    features = pd.read_csv(args.counts / "selected-features.csv").sort_values("soma_joinid")
    populations = pd.read_parquet(args.counts / "panel-populations.parquet")
    obs = np.load(args.counts / "panel-obs-joinids.npy")
    mapping = np.load(args.counts / "observation-population-map.npy", mmap_mode="r")
    if not np.isin(mapping[obs], populations.population_id.to_numpy()).all():
        raise ValueError("Panel cells and populations disagree")
    if len(obs) != int(populations.n_cells.sum()):
        raise ValueError("Panel denominator mismatch")
    # Dense local population index keeps the accumulators small.
    local = np.full(int(mapping.max()) + 1, -1, np.int32)
    local[populations.population_id.to_numpy()] = np.arange(len(populations))
    var_ids = features.soma_joinid.to_numpy()
    var_map = np.full(int(var_ids.max()) + 1, -1, np.int32)
    var_map[var_ids] = np.arange(len(var_ids))
    manifest = {
        "census_release": base["census_release"],
        "base_manifest_sha256": digest(args.counts / "manifest.json"),
        "panel_populations_sha256": digest(args.counts / "panel-populations.parquet"),
        "panel_obs_sha256": digest(args.counts / "panel-obs-joinids.npy"),
        "n_cells": len(obs),
        "n_genes": len(var_ids),
        "chunk": args.chunk,
        "chunks": -(-len(obs) // args.chunk),
        "layer": "raw",
        "method_version": 1,
    }
    path = out / "manifest.json"
    if path.exists() and json.loads(path.read_text()) != manifest:
        raise ValueError("Resume identity mismatch; use a new output directory")
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    ctx = cc.get_default_soma_context(
        {"soma.init_buffer_bytes": 1024**3, "py.init_buffer_bytes": 1024**3}
    )
    done = 0
    with cc.open_soma(census_version=manifest["census_release"], context=ctx) as census:
        x = census["census_data"]["homo_sapiens"].ms["RNA"].X["raw"]
        for chunk_id in range(k, manifest["chunks"], n):
            part = out / f"part-{chunk_id:05d}.parquet"
            if part.exists():
                continue
            if args.max_chunks is not None and done >= args.max_chunks:
                break
            cells = obs[chunk_id * args.chunk : (chunk_id + 1) * args.chunk]
            row_map = np.full(int(cells.max()) + 1, -1, np.int32)
            row_map[cells] = local[mapping[cells]]
            tables = x.read(coords=(cells, var_ids)).tables()
            n_positive, n_ge2, n_fractional, count_sum, count_max = summarize(
                tables, row_map, var_map, len(populations), len(var_ids)
            )
            row, col = np.nonzero(n_positive)
            result = pd.DataFrame(
                {
                    "population_id": populations.population_id.to_numpy()[row],
                    "gene_id": features.feature_id.to_numpy()[col],
                    "n_positive": n_positive[row, col],
                    "n_ge2": n_ge2[row, col],
                    "count_sum": count_sum[row, col],
                    "count_max": count_max[row, col],
                    "n_fractional": n_fractional[row, col],
                }
            )
            temp = part.with_suffix(".partial")
            result.to_parquet(temp, index=False)
            temp.replace(part)
            done += 1
            print(f"chunk {chunk_id}: {len(cells):,} cells, {len(result):,} rows", flush=True)
    print("Requested chunks complete", flush=True)


if __name__ == "__main__":
    main()

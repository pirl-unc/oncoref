#!/usr/bin/env python3
"""Independent adult-cerebellum audit of every protein record, by donor.

Use the three disjoint Siletti cerebellum dissections, not overlapping cell-type
exports. Cell IDs are checked for duplication; missing features remain missing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from audit_plac1_brain import read_column
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    catalog = pd.read_csv(
        ROOT / "docs/audits/target-selection-20261109/category-profiles/all-protein-profiles.csv.gz"
    )
    requested = {g for ids in catalog.member_gene_ids for g in json.loads(ids)}
    collection = json.loads((args.inputs / "collection.json").read_text())
    datasets = [
        d for d in collection["datasets"] if d["title"].startswith("Dissection: Cerebellum (CB)")
    ]
    seen = set()
    summaries, sources, cell_summaries = [], [], []
    for dataset in datasets:
        path = args.inputs / (dataset["dataset_id"] + ".h5ad")
        with path.open("rb") as stream:
            checksum = hashlib.file_digest(stream, "sha256").hexdigest()
        sources.append(
            {
                "dataset_id": dataset["dataset_id"],
                "title": dataset["title"],
                "url": next(a["url"] for a in dataset["assets"] if a["filetype"] == "H5AD"),
                "sha256": checksum,
                "bytes": path.stat().st_size,
            }
        )
        with h5py.File(path) as handle:
            obs = pd.DataFrame(
                {
                    name: read_column(handle, "obs/" + name)
                    for name in [
                        "CellID",
                        "donor_id",
                        "cell_type",
                        "cell_type_ontology_term_id",
                        "development_stage",
                        "disease_ontology_term_id",
                        "is_primary_data",
                        "sample_id",
                        "fraction_mitochondrial",
                        "fraction_unspliced",
                        "total_UMIs",
                        "total_genes",
                    ]
                }
            )
            if not obs.CellID.is_unique or seen.intersection(obs.CellID):
                raise ValueError("Overlapping dissection exports")
            seen.update(obs.CellID)
            if not (obs.is_primary_data & obs.disease_ontology_term_id.eq("PATO:0000461")).all():
                raise ValueError("Unexpected nonprimary or non-normal observations")
            if not obs.development_stage.isin(
                ["29-year-old stage", "42-year-old stage", "50-year-old stage"]
            ).all():
                raise ValueError("Unexpected age")
            ids = read_column(handle, "var/_index")
            names = read_column(handle, "var/feature_name")
            available = np.isin(ids, list(requested)) & ~read_column(
                handle, "var/feature_is_filtered"
            )
            gene_ids, symbols = ids[available], names[available]
            if len(set(gene_ids)) != len(gene_ids):
                raise ValueError("Duplicated feature IDs")
            mapping = np.full(len(ids), -1, np.int32)
            mapping[available] = np.arange(available.sum())
            labels = (
                obs[["donor_id", "cell_type", "cell_type_ontology_term_id"]]
                .drop_duplicates()
                .reset_index(drop=True)
            )
            labels["population"] = np.arange(len(labels))
            obs = obs.merge(
                labels,
                on=["donor_id", "cell_type", "cell_type_ontology_term_id"],
                how="left",
                sort=False,
                validate="many_to_one",
            )
            shape = (len(labels), len(gene_ids))
            pos = np.zeros(shape, np.int32)
            ge2 = np.zeros(shape, np.int32)
            total = np.zeros(shape, np.float64)
            maximum = np.zeros(shape, np.float32)
            if handle["X"].attrs["encoding-type"] != "csr_matrix":
                raise ValueError("Expected raw CSR counts")
            ptr = handle["X/indptr"][:]
            computed_totals = np.zeros(len(obs))
            for start in range(0, len(obs), 5000):
                end = min(start + 5000, len(obs))
                left, right = int(ptr[start]), int(ptr[end])
                vals = handle["X/data"][left:right]
                indices = handle["X/indices"][left:right]
                if not np.all(np.isfinite(vals) & (vals >= 0) & (vals == np.floor(vals))):
                    raise ValueError("Not integer raw counts")
                block = sparse.csr_matrix(
                    (vals, indices, ptr[start : end + 1] - left), shape=(end - start, len(ids))
                )
                computed_totals[start:end] = np.asarray(block.sum(axis=1)).ravel()
                block = block[:, available]
                for pop in obs.iloc[start:end].population.unique():
                    selected = block[obs.iloc[start:end].population.to_numpy() == pop]
                    pos[pop] += np.asarray((selected > 0).sum(axis=0)).ravel()
                    ge2[pop] += np.asarray((selected >= 2).sum(axis=0)).ravel()
                    total[pop] += np.asarray(selected.sum(axis=0)).ravel()
                    maximum[pop] = np.maximum(maximum[pop], selected.max(axis=0).toarray().ravel())
            obs["computed_total_counts"] = computed_totals
            obs["dataset_id"] = dataset["dataset_id"]
            # Metadata totals may refer to the original reference; retain both.
            cell_summaries.append(
                obs.groupby(["dataset_id", "donor_id", "cell_type"])
                .agg(
                    n_cells=("CellID", "size"),
                    median_count=("computed_total_counts", "median"),
                    median_reported_UMIs=("total_UMIs", "median"),
                    median_genes=("total_genes", "median"),
                    median_mito=("fraction_mitochondrial", "median"),
                    median_unspliced=("fraction_unspliced", "median"),
                )
                .reset_index()
            )
            denominators = obs.groupby("population").size()
            for pop, label in labels.iterrows():
                result = pd.DataFrame(
                    {
                        "gene_id": gene_ids,
                        "symbol": symbols,
                        "n_positive": pos[pop],
                        "n_ge2": ge2[pop],
                        "count_sum": total[pop],
                        "max_count": maximum[pop],
                        "n_cells": denominators[pop],
                        "dataset_id": dataset["dataset_id"],
                        "donor_id": label.donor_id,
                        "cell_type": label.cell_type,
                        "cell_type_ontology_term_id": label.cell_type_ontology_term_id,
                    }
                )
                summaries.append(result)
        print(dataset["title"], len(obs), "cells", len(gene_ids), "catalog genes", flush=True)
    all_rows = pd.concat(summaries, ignore_index=True)
    all_rows.to_parquet(args.output / "gene-populations.parquet", index=False)
    pd.concat(cell_summaries, ignore_index=True).to_csv(
        args.output / "population-qc.csv", index=False
    )
    granule = all_rows[all_rows.cell_type_ontology_term_id.eq("CL:0001031")]
    combined = []
    for scope, frame in [("all", granule), *list(granule.groupby("donor_id"))]:
        x = (
            frame.groupby(["gene_id", "symbol"])
            .agg(
                n_cells=("n_cells", "sum"),
                n_positive=("n_positive", "sum"),
                n_ge2=("n_ge2", "sum"),
                count_sum=("count_sum", "sum"),
                max_count=("max_count", "max"),
            )
            .reset_index()
        )
        x["scope"] = scope
        x["pct_positive"] = 100 * x.n_positive / x.n_cells
        x["pct_ge2"] = 100 * x.n_ge2 / x.n_cells
        combined.append(x)
    gene_result = pd.concat(combined, ignore_index=True)
    gene_result.to_csv(
        args.output / "granule-gene-detection.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    protein_rows = []
    for scope, frame in gene_result.groupby("scope"):
        frame = frame.set_index("gene_id")
        for protein in catalog.itertuples():
            members = json.loads(protein.member_gene_ids)
            found = [g for g in members if g in frame.index]
            missing = [g for g in members if g not in frame.index]
            row = {
                "entity_id": protein.entity_id,
                "protein": protein.protein,
                "category": protein.category,
                "scope": scope,
                "available_gene_ids": json.dumps(found),
                "missing_gene_ids": json.dumps(missing),
                "coverage": "unmeasured" if not found else "partial" if missing else "complete",
            }
            if found:
                chosen = frame.loc[found].sort_values("pct_positive", ascending=False).iloc[0]
                row.update(chosen.drop("scope").to_dict())
                row["gene_id"] = chosen.name
            protein_rows.append(row)
    proteins = pd.DataFrame(protein_rows)
    proteins.to_csv(
        args.output / "granule-protein-detection.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    (args.output / "manifest.json").write_text(
        json.dumps(
            {
                "sources": sources,
                "unique_nuclei": len(seen),
                "catalog_proteins": len(catalog),
                "classification_changed": False,
                "raw_count_unit": "integer counts from 10x UMI assay; X is source raw layer",
                "replication": "one independent study, three donors, disjoint anatomical partitions; not three studies",
                "remaining_QC": "alignment, ambient correction, upstream gene filtering and donor pathology review",
            },
            indent=2,
        )
        + "\n"
    )
    print(
        gene_result.query(
            "scope == 'all' and symbol in ['PLAC1','CBY3','DMRT1','CAGE1','TMEM269','NUTM1','SYCP1','MAGEB10','PRAME']"
        ).to_string(index=False),
        flush=True,
    )


if __name__ == "__main__":
    main()

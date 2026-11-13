#!/usr/bin/env python3
"""Compare current antigen sets and familiar targets in the Sepp brain dataset.

All gene mapping is by exact Ensembl ID. Missing features remain unmeasured.
RNA counts are not protein detections. Requires h5py, scipy, numpy and pandas.
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
BASE = ROOT / "docs/audits/target-selection-20261109/category-profiles"
TARGET_CATEGORIES = {"CTA", "Cancer Antigen only", "Reproductive Antigen only"}


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def extract(path: Path, requested: set[str]):
    """Read CSR chunks once, retaining just the requested, unfiltered features."""
    with h5py.File(path, "r") as handle:
        if handle["X"].attrs["encoding-type"] != "csr_matrix":
            raise ValueError("Expected CSR raw counts in X")
        var_ids = read_column(handle, "var/_index")
        names = read_column(handle, "var/feature_name")
        filtered = read_column(handle, "var/feature_is_filtered")
        keep = np.array([gene in requested for gene in var_ids]) & ~filtered
        genes = pd.DataFrame({"gene_id": var_ids[keep], "dataset_symbol": names[keep]})
        if not genes.gene_id.is_unique:
            raise ValueError("Duplicate feature IDs")
        mapping = np.full(len(var_ids), -1, dtype=np.int32)
        mapping[keep] = np.arange(keep.sum(), dtype=np.int32)
        columns = [
            "donor_id",
            "author_stage",
            "development_stage",
            "cell_type",
            "cell_type_ontology_term_id",
            "disease_ontology_term_id",
            "is_primary_data",
        ]
        cells = pd.DataFrame({key: read_column(handle, f"obs/{key}") for key in columns})
        rows, cols, values = [], [], []
        pointers = handle["X/indptr"][:]
        for start in range(0, len(cells), 10_000):
            end = min(start + 10_000, len(cells))
            left, right = int(pointers[start]), int(pointers[end])
            mapped = mapping[handle["X/indices"][left:right]]
            selected = np.flatnonzero(mapped >= 0)
            data = handle["X/data"][left:right][selected]
            if not np.all(np.isfinite(data) & (data >= 0) & (data == np.floor(data))):
                raise ValueError("Input is not nonnegative integer raw counts")
            rows.append(np.searchsorted(pointers, selected + left, side="right") - 1)
            cols.append(mapped[selected])
            values.append(data.astype(np.int32))
    matrix = sparse.csr_matrix(
        (np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))),
        shape=(len(cells), len(genes)),
    )
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    return cells, genes, matrix


def selections(cells: pd.DataFrame):
    normal = cells.is_primary_data & cells.disease_ontology_term_id.eq("PATO:0000461")
    adult = normal & cells.author_stage.str.startswith("adult")
    granule = cells.cell_type_ontology_term_id.eq("CL:0001031")
    out = {"adult_granule": adult & granule, "adult_all_cells": adult}
    for donor in sorted(cells.loc[adult, "donor_id"].unique()):
        out[f"adult_granule_donor_{donor}"] = adult & granule & cells.donor_id.eq(donor)
    for cell_type in sorted(cells.loc[adult, "cell_type"].unique()):
        out[f"adult_type:{cell_type}"] = adult & cells.cell_type.eq(cell_type)
    prenatal = cells.author_stage.str.startswith("CS") | cells.author_stage.str.contains("wpc")
    child = cells.author_stage.str.startswith(("infant", "newborn", "toddler"))
    out["prenatal_granule"] = normal & granule & prenatal
    out["child_granule"] = normal & granule & child
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, default=BASE / "all-protein-profiles.csv.gz")
    parser.add_argument("--familiar", type=Path, default=BASE / "familiar-cta-profiles.csv")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    catalog = pd.read_csv(args.profiles, low_memory=False)
    familiar = set(pd.read_csv(args.familiar).entity_id)
    proteins = catalog[
        catalog.category.isin(TARGET_CATEGORIES) | catalog.entity_id.isin(familiar)
    ].copy()
    proteins["familiar_target"] = proteins.entity_id.isin(familiar)
    requested = {gene for value in proteins.member_gene_ids for gene in json.loads(value)}
    cells, genes, matrix = extract(args.h5ad, requested)
    genes["column"] = np.arange(len(genes))
    by_id = genes.set_index("gene_id")
    summaries = []
    for scope, mask in selections(cells).items():
        selected = matrix[mask.to_numpy()]
        if selected.shape[0] == 0:
            continue
        summary = genes.copy()
        summary["scope"] = scope
        summary["n_cells"] = selected.shape[0]
        summary["n_positive"] = np.asarray((selected > 0).sum(axis=0)).ravel()
        summary["n_ge2"] = np.asarray((selected >= 2).sum(axis=0)).ravel()
        summary["total_counts"] = np.asarray(selected.sum(axis=0)).ravel()
        summary["max_raw_count"] = selected.max(axis=0).toarray().ravel()
        summary["pct_positive"] = 100 * summary.n_positive / summary.n_cells
        summary["pct_ge2"] = 100 * summary.n_ge2 / summary.n_cells
        summary["detection_screen_hit"] = (
            summary.n_cells.ge(3000) & summary.n_positive.ge(10) & summary.pct_positive.ge(1)
        )
        summaries.append(summary.drop(columns="column"))
    gene_results = pd.concat(summaries, ignore_index=True)
    gene_results.to_csv(
        args.output / "gene-detection.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    records = []
    for scope, summary in gene_results.groupby("scope", sort=False):
        lookup = summary.set_index("gene_id")
        for protein in proteins.itertuples():
            members = json.loads(protein.member_gene_ids)
            available = [gene for gene in members if gene in by_id.index]
            missing = [gene for gene in members if gene not in by_id.index]
            row = {
                "entity_id": protein.entity_id,
                "protein": protein.protein,
                "category": protein.category,
                "classification": protein.classification,
                "familiar_target": protein.familiar_target,
                "published_nomination": protein.published_nomination,
                "member_gene_ids": protein.member_gene_ids,
                "available_gene_ids": json.dumps(available),
                "missing_gene_ids": json.dumps(missing),
                "coverage": "complete"
                if not missing
                else ("partial" if available else "unmeasured"),
                "scope": scope,
            }
            if available:
                # Match OncoRef's any-member rule; never add gene detection frequencies.
                best = (
                    lookup.loc[available]
                    .sort_values(["n_positive", "total_counts"], ascending=False)
                    .iloc[0]
                )
                row["max_member_gene_id"] = best.name
                for key in [
                    "dataset_symbol",
                    "n_cells",
                    "n_positive",
                    "n_ge2",
                    "total_counts",
                    "max_raw_count",
                    "pct_positive",
                    "pct_ge2",
                    "detection_screen_hit",
                ]:
                    row[key] = best[key]
            records.append(row)
    result = pd.DataFrame(records)
    result.to_csv(
        args.output / "protein-detection.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    adult = result[result.scope.eq("adult_granule")].copy()
    adult.to_csv(args.output / "adult-granule-all-proteins.csv", index=False)
    adult[adult.category.eq("CTA")].sort_values("pct_positive", ascending=False).to_csv(
        args.output / "adult-granule-current-ctas.csv", index=False
    )
    adult[adult.familiar_target].sort_values("pct_positive", ascending=False).to_csv(
        args.output / "adult-granule-familiar-targets.csv", index=False
    )
    counts = []
    for category, frame in adult.groupby("category"):
        measured = frame[frame.coverage.ne("unmeasured")]
        counts.append(
            {
                "category": category,
                "n_proteins": len(frame),
                "complete": int(frame.coverage.eq("complete").sum()),
                "partial": int(frame.coverage.eq("partial").sum()),
                "unmeasured": int(frame.coverage.eq("unmeasured").sum()),
                "zero_detected_measured": int(measured.n_positive.eq(0).sum()),
                "below_1pct_positive": int(
                    (measured.n_positive.gt(0) & measured.pct_positive.lt(1)).sum()
                ),
                "at_least_1pct": int(measured.pct_positive.ge(1).sum()),
                "at_least_1pct_ge2_counts": int(measured.pct_ge2.ge(1).sum()),
            }
        )
    pd.DataFrame(counts).to_csv(args.output / "category-summary.csv", index=False)
    manifest = {
        "gene_mapping": "exact Ensembl IDs; unavailable or filtered features are unmeasured",
        "protein_summary": "highest-detection available member; never sum frequencies across paralogs",
        "classification_changed": False,
        "screen_scope": "1% detection, >=3000 nuclei, >=10 positive; local population comparison only, not the full global rare-cell policy",
        "raw_count_two_sensitivity": "descriptive only; not an adopted exclusion threshold",
        "h5ad_sha256": sha256(args.h5ad),
        "profiles_sha256": sha256(args.profiles),
        "familiar_sha256": sha256(args.familiar),
        "selected_proteins": len(proteins),
        "requested_genes": len(requested),
        "measured_genes": len(genes),
        "scope_sizes": {key: int(mask.sum()) for key, mask in selections(cells).items()},
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(counts).to_string(index=False))
    print(
        adult[adult.familiar_target][
            ["protein", "category", "coverage", "n_positive", "pct_positive", "n_ge2"]
        ]
        .sort_values("pct_positive", ascending=False)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()

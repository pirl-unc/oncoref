#!/usr/bin/env python3
"""Extract PLAC1 counts and donor metadata from the Sepp human cerebellum H5AD.

Requires h5py, numpy, pandas and matplotlib. The input is the public raw-count
asset named in the audit manifest; no environment packages are changed here.
This is an evidence audit, not a catalog classification update.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

GENE = "ENSG00000170965"
GRANULE = "CL:0001031"


def read_column(handle: h5py.File, path: str) -> np.ndarray:
    column = handle[path]
    if isinstance(column, h5py.Group):
        codes = column["codes"][:]
        categories = column["categories"].asstr()[:]
        return np.where(codes >= 0, categories[np.maximum(codes, 0)], None)
    return column.asstr()[:] if column.dtype.kind in "OSU" else column[:]


def summarize(cells: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    result = cells.groupby(by, observed=True, dropna=False).agg(
        n_cells=("positive", "size"),
        n_positive=("positive", "sum"),
        n_ge2=("two_or_more", "sum"),
        max_raw_count=("PLAC1_raw_count", "max"),
        total_counts=("PLAC1_raw_count", "sum"),
    )
    result["pct_positive"] = 100 * result.n_positive / result.n_cells
    return result.reset_index()


def extract(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as handle:
        gene_indices = np.flatnonzero(read_column(handle, "var/_index") == GENE)
        if len(gene_indices) != 1:
            raise ValueError(f"Expected one {GENE} feature, found {len(gene_indices)}")
        gene_index = int(gene_indices[0])
        if read_column(handle, "var/feature_name")[gene_index] != "PLAC1":
            raise ValueError("Gene ID/symbol mismatch")
        if handle["X"].attrs["encoding-type"] != "csr_matrix":
            raise ValueError("Expected CSR raw counts in X")
        columns = [
            "donor_id",
            "development_stage",
            "development_stage_ontology_term_id",
            "author_stage",
            "cell_type",
            "cell_type_ontology_term_id",
            "author_cell_type",
            "dev_state",
            "tissue",
            "tissue_fragment",
            "batch",
            "sex",
            "disease",
            "disease_ontology_term_id",
            "is_primary_data",
        ]
        cells = pd.DataFrame({key: read_column(handle, f"obs/{key}") for key in columns})
        cells.insert(0, "cell_id", read_column(handle, "obs/_index"))
        if not cells.cell_id.is_unique:
            raise ValueError("Duplicate cell IDs")
        counts = np.zeros(len(cells), dtype=np.float64)
        pointers = handle["X/indptr"][:]
        indices, values = handle["X/indices"], handle["X/data"]
        for start in range(0, len(cells), 10_000):
            end = min(start + 10_000, len(cells))
            left, right = int(pointers[start]), int(pointers[end])
            positions = np.flatnonzero(indices[left:right] == gene_index)
            rows = np.searchsorted(pointers, positions + left, side="right") - 1
            np.add.at(counts, rows, values[left:right][positions])
        if not np.all(np.isfinite(counts) & (counts >= 0) & (counts == np.floor(counts))):
            raise ValueError("Expected nonnegative integer raw counts")
        cells["PLAC1_raw_count"] = counts.astype(np.int64)
    stage = cells.author_stage
    cells["age_group"] = np.select(
        [
            stage.str.startswith("adult"),
            stage.str.startswith("CS") | stage.str.contains("wpc"),
            stage.str.startswith(("infant", "newborn", "toddler")),
        ],
        ["Adult", "Prenatal", "Postnatal child"],
        default="Unknown",
    )
    cells["positive"] = cells.PLAC1_raw_count > 0
    cells["two_or_more"] = cells.PLAC1_raw_count >= 2
    return cells


def plot_donors(donors: pd.DataFrame, output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"Prenatal": "#6386a8", "Postnatal child": "#c89242", "Adult": "#a63d40"}
    donors = donors.sort_values(["age_group", "development_stage", "donor_id"])
    labels = [f"{row.donor_id} · {row.development_stage}" for row in donors.itertuples()]
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.barh(labels, donors.pct_positive, color=[colors.get(x, "#777777") for x in donors.age_group])
    ax.axvline(1, color="#444444", linestyle="--", linewidth=1, label="OncoRef 1% screen")
    for y, row in enumerate(donors.itertuples()):
        ax.text(
            row.pct_positive + 0.12,
            y,
            f"{row.n_positive:,}/{row.n_cells:,}",
            va="center",
            fontsize=9,
        )
    ax.set_xlim(0, 14)
    ax.invert_yaxis()
    ax.set_xlabel("Granule nuclei with at least one PLAC1 count (%)")
    ax.set_title(
        "PLAC1 RNA is detected in both normal-labelled adult donors",
        loc="left",
        weight="bold",
        pad=14,
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower right", frameon=False)
    fig.text(
        0.01,
        0.015,
        "Sepp et al. · CELLxGENE asset 30626c44… · primary, normal-labelled tissue only.\nCounts are RNA detections, not protein abundance. Labels show positive / total nuclei.",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(output / "granule-by-donor.png", dpi=160)
    fig.savefig(output / "granule-by-donor.svg")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    cells = extract(args.h5ad)
    cells.to_csv(
        args.output / "sepp-plac1-cells.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    normal = cells[cells.is_primary_data & cells.disease_ontology_term_id.eq("PATO:0000461")]
    granule = normal[normal.cell_type_ontology_term_id.eq(GRANULE)]
    donors = summarize(granule, ["age_group", "donor_id", "development_stage", "sex"])
    donors.to_csv(args.output / "granule-by-donor.csv", index=False)
    summarize(granule, ["age_group", "author_stage", "development_stage"]).to_csv(
        args.output / "granule-by-age.csv", index=False
    )
    summarize(normal, ["age_group", "cell_type"]).to_csv(
        args.output / "all-cells-by-age-type.csv", index=False
    )
    adult = granule[granule.age_group.eq("Adult")]
    with args.h5ad.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    summary = {
        "gene_id": GENE,
        "source_sha256": digest,
        "source_bytes": args.h5ad.stat().st_size,
        "input_cells": len(cells),
        "normal_granule_cells": len(granule),
        "normal_granule_positive": int(granule.positive.sum()),
        "adult_normal_granule_cells": len(adult),
        "adult_normal_granule_positive": int(adult.positive.sum()),
        "adult_normal_granule_positive_donors": int(
            adult.loc[adult.positive, "donor_id"].nunique()
        ),
        "adult_positive_count_distribution": adult.loc[adult.positive, "PLAC1_raw_count"]
        .value_counts()
        .sort_index()
        .to_dict(),
    }
    (args.output / "sepp-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    plot_donors(donors, args.output)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

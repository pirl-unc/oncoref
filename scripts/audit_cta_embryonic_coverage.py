#!/usr/bin/env python3
"""Audit missing atlas rows against the authors' independent count matrix.

This does not impute TPM, alter selection, or infer fetal-wide absence. Rscript
reads the source RDS without installing packages; all output retains native counts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-source-walkthrough-20261028"
AUTHOR_COMMIT = "3e6ab8546e57decc5ec09770720709b5fa6428ef"
TARGETS = ["PRAME", "MAGEA4", "CTAG1A", "CTAG1B", "XAGE1A", "XAGE1B", "CT83", "TBXT"]
R_CODE = r"""
args <- commandArgs(trailingOnly=TRUE)
d <- readRDS(args[1])
stopifnot(is.data.frame(d), nrow(d) == 1195)
targets <- c("PRAME", "MAGEA4", "CTAG1A", "CTAG1B", "XAGE1A", "XAGE1B", "CT83", "TBXT")
rows <- lapply(targets, function(g) {
  present <- g %in% names(d)
  x <- if(present) d[[g]] else NA_real_
  if(present) stopifnot(all(is.finite(x)), all(x >= 0))
  data.frame(source_symbol=g, gene_column_present=present, cells=nrow(d),
    total_reported_reads=if(present) sum(x) else NA_real_,
    positive_cells=if(present) sum(x>0) else NA_integer_,
    maximum_cell_reads=if(present) max(x) else NA_real_)
})
write.csv(do.call(rbind, rows), args[2], row.names=FALSE, na="")
"""


def digest(path: Path) -> str:
    """Hash exact downloaded source bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(sources: Path, atlas: Path) -> None:
    """Record source membership and direct count observations for comparators."""
    with tempfile.TemporaryDirectory(prefix="oncoref-embryo-") as directory:
        script, table = Path(directory) / "audit.R", Path(directory) / "counts.csv"
        script.write_text(R_CODE)
        subprocess.run(
            ["Rscript", str(script), str(sources / "author-raw.rds"), str(table)], check=True
        )
        counts = pd.read_csv(table)
    from oncoref import cta_selection

    identity = {}
    for row in cta_selection.results().itertuples():
        for member in row.gene_identities:
            if member["symbol"] in TARGETS:
                identity[member["symbol"]] = (member["gene_id"], row.entity_id)
    counts["gene_id"] = counts.source_symbol.map(lambda symbol: identity[symbol][0])
    counts["protein_equivalence_group_id"] = counts.source_symbol.map(
        lambda symbol: identity[symbol][1]
    )
    for name, path in [("atlas_tpm_row", atlas), ("atlas_count_row", sources / "atlas-raw.zip")]:
        with zipfile.ZipFile(path) as archive:
            filename = next(n for n in archive.namelist() if n.endswith("_rows"))
            genes = set(pd.read_csv(archive.open(filename), sep="\t", header=None)[0])
        counts[name] = counts.gene_id.isin(genes)
    counts["unit"] = "author-reported unnormalized reads"
    counts["scope"] = "One Carnegie-stage-7 embryo; not fetal-wide absence"
    counts.to_csv(OUT / "embryonic-comparator-coverage.csv", index=False)
    receipt = {
        "author_commit": AUTHOR_COMMIT,
        "sources": {
            "author_raw_counts": {
                "url": f"https://github.com/ScialdoneLab/human-gastrula-shiny/blob/{AUTHOR_COMMIT}/raw_reads.rds",
                "sha256": digest(sources / "author-raw.rds"),
                "matrix_shape": {"cells": 1195, "gene_columns": 57490},
            },
            "atlas_counts": {
                "url": "https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-9388/download/zip?fileType=quantification-raw",
                "sha256": digest(sources / "atlas-raw.zip"),
                "matrix_shape": {"gene_rows": 25330, "cells": 1195},
                "note": "The Raw counts download itself contains aggregated_filtered_counts; it does not restore excluded gene rows.",
            },
            "atlas_tpm": {
                "sha256": digest(atlas),
                "receipt": "../cta-quantitative-20261026/developmental-source.json",
            },
        },
        "script_sha256": digest(Path(__file__)),
        "output_sha256": digest(OUT / "embryonic-comparator-coverage.csv"),
        "interpretation": "Zero reported reads is a measured nondetection in this source; absent columns remain missing. Counts are not converted to TPM or inserted in the selection funnel.",
    }
    (OUT / "embryonic-coverage-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, default=ROOT / "tmp/cta-source-review")
    parser.add_argument("--atlas", type=Path, default=ROOT / "tmp/cta-quantitative/atlas_tpm")
    args = parser.parse_args()
    build(args.sources, args.atlas)

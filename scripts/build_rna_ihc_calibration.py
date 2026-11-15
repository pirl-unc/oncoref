#!/usr/bin/env python3
"""Calibrate bulk RNA against antibody staining, to read an RNA value as a protein grade.

Pairs HPA tissue RNA (consensus nTPM) with the strongest HPA staining level of the
same gene in the same tissue, using only antibodies of Enhanced reliability. Writes
the share of pairs that stain at each RNA band and the grade cutoffs drawn from it.

    python scripts/build_rna_ihc_calibration.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from oncoref import hpa

OUT = Path(__file__).resolve().parents[1] / "oncoref/data/rna-ihc-calibration.json"
LEVELS = {"Not detected": 0, "Low": 1, "Medium": 2, "High": 3}
BANDS = [0, 0.5, 1, 2, 5, 10, 20, 50, 100, np.inf]
# Grade floors, in nTPM. Chosen from the bands below: under 1, fewer than one pair in
# six stains; from 5, most pairs stain; from 20, three in four stain and most of those
# at medium or high level.
GRADES = {"low": 1.0, "medium": 5.0, "high": 20.0}


def main():
    stain = hpa.hpa_normal_tissue()
    stain = stain[
        stain.Reliability.eq("Enhanced") & stain.Level.isin(LEVELS) & stain.Tissue.notna()
    ]
    stain = stain.assign(
        tissue=stain.Tissue.str.replace(r" \d+$", "", regex=True), level=stain.Level.map(LEVELS)
    )
    strongest = stain.groupby(["Gene", "tissue"])["level"].max().reset_index()
    rna = hpa.hpa_rna_consensus().rename(columns={"Tissue": "tissue"})
    pairs = strongest.merge(rna[["Gene", "tissue", "nTPM"]], on=["Gene", "tissue"])
    pairs["band"] = pd.cut(pairs.nTPM, BANDS, right=False)
    bands = [
        {
            "ntpm_from": float(band.left),
            "ntpm_below": None if np.isinf(band.right) else float(band.right),
            "pairs": len(rows),
            "stained": round(float(rows.level.ge(1).mean()), 3),
            "medium_or_high": round(float(rows.level.ge(2).mean()), 3),
            "high": round(float(rows.level.ge(3).mean()), 3),
        }
        for band, rows in pairs.groupby("band", observed=True)
    ]
    medians = {
        name: round(float(pairs[pairs.level.eq(value)].nTPM.median()), 1)
        for name, value in LEVELS.items()
    }
    document = {
        "schema": "oncoref.rna-ihc-calibration",
        "schema_version": 1,
        "source": "Human Protein Atlas normal-tissue staining and tissue RNA consensus, as "
        "loaded by oncoref.hpa",
        "reliability": "Enhanced",
        "pairs": len(pairs),
        "genes": int(pairs.Gene.nunique()),
        "tissues": int(pairs.tissue.nunique()),
        "median_ntpm_by_staining_level": medians,
        "bands": bands,
        "grade_floors_ntpm": GRADES,
        "limits": [
            "Staining levels overlap widely in RNA. RNA separates stained from unstained "
            "tissue far better than low from high staining.",
            "The cutoffs are in HPA nTPM. FPKM from another study and counts per million "
            "from single cells are read on the same scale only as an approximation.",
        ],
    }
    OUT.write_text(json.dumps(document, indent=2) + "\n")
    print(json.dumps({"pairs": document["pairs"], "medians": medians, "grades": GRADES}))


if __name__ == "__main__":
    main()

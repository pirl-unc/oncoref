#!/usr/bin/env python3
"""Build ``oncoref/data/tissue-marker-panels.csv`` from GTEx v11 tissue medians.

Each panel holds up to six protein-coding genes that a normal tissue expresses at
100 TPM or more and that no tissue outside its group reaches 1% of. Library QC uses
the panels to recognise normal tissue inside a tumor RNA library: tissue brings all of
a panel's transcripts at the same fraction of their tissue level, whereas a tumor that
switches on one tissue-specific gene does not.

The testis panel is the late-spermatid program (seminiferous tubules with complete
spermatogenesis). Panels with fewer than three qualifying genes are left out.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import gtex  # noqa: E402

OUT = ROOT / "oncoref/data/tissue-marker-panels.csv"
SOURCE = "GTEx v11 gene median TPM (GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3)"
ENSEMBL_RELEASE = 115
MIN_TPM = 100.0
MAX_OUTSIDE_FRACTION = 0.01
PANEL_SIZE = 6
MIN_PANEL = 3


def groups(columns: list[str]) -> dict[str, list[str]]:
    """GTEx tissue columns per panel; tissues sharing a program are one group."""

    def starting(prefix):
        return [c for c in columns if c.startswith(prefix)]

    return {
        "testis": ["Testis"],
        "pancreas": starting("Pancreas"),
        "liver": starting("Liver"),
        "stomach": [*starting("Stomach"), "Esophagus_Gastroesophageal_Junction"],
        "striated_muscle": ["Muscle_Skeletal", *starting("Heart")],
        "thyroid": ["Thyroid"],
        "salivary_gland": ["Minor_Salivary_Gland"],
        "pituitary": ["Pituitary"],
        "squamous_epithelium": [
            "Esophagus_Mucosa",
            "Vagina",
            "Cervix_Ectocervix",
            *starting("Skin"),
        ],
        "prostate": ["Prostate"],
        "brain": starting("Brain"),
        "lung": ["Lung"],
        "kidney": starting("Kidney"),
        "adrenal_gland": ["Adrenal_Gland"],
    }


def main() -> None:
    from pyensembl import EnsemblRelease

    genome = EnsemblRelease(ENSEMBL_RELEASE)
    table = gtex.tissue_rna()
    table["gene_id"] = table.Name.str.split(".").str[0]
    table = table.drop_duplicates("gene_id").set_index("gene_id")
    symbols = table.Description
    tpm = table.drop(columns=["Name", "Description"]).astype(float)

    def protein_coding(gene_id: str) -> bool:
        try:
            return genome.gene_by_id(gene_id).biotype == "protein_coding"
        except ValueError:
            return False

    rows = []
    for panel, cols in groups(list(tpm.columns)).items():
        inside = tpm[cols].max(axis=1)
        outside = tpm[[c for c in tpm.columns if c not in cols]].max(axis=1)
        ok = inside.ge(MIN_TPM) & outside.le(MAX_OUTSIDE_FRACTION * inside)
        ranked = inside[ok].sort_values(ascending=False)
        chosen = [g for g in ranked.index if protein_coding(g)][:PANEL_SIZE]
        if len(chosen) < MIN_PANEL:
            print(f"{panel}: {len(chosen)} genes, left out")
            continue
        for gene_id in chosen:
            rows.append(
                {
                    "panel": panel,
                    "gene_id": gene_id,
                    "symbol": symbols[gene_id],
                    "reference_tpm": round(float(inside[gene_id]), 1),
                    "reference_tissue": tpm.loc[gene_id, cols].idxmax(),
                    "max_outside_tpm": round(float(outside[gene_id]), 2),
                    "source": SOURCE,
                }
            )
    frame = pd.DataFrame(rows)
    if frame.gene_id.duplicated().any():
        raise ValueError("A gene marks two panels")
    frame.to_csv(OUT, index=False)
    print(frame.groupby("panel").symbol.apply(", ".join).to_string())


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Regenerate six unpaired RNA/IHC plots from the pinned oncoref reference.

The checked-in gene list preserves the 253-gene legacy plotting scope; the
underlying oncoref references and APIs are genome-wide. No tsarina dependency.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from oncoref import hpa_cancer_rna_ihc_comparison, hpa_cancer_sources


def plot(gene_list: Path, output: Path):
    output.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(gene_list)
    frame = hpa_cancer_rna_ihc_comparison(genes.gene_id, cohort="TCGA", threshold=1)
    frame = frame.merge(genes, on="gene_id", how="left", validate="many_to_one")
    frame.to_csv(output / "comparison.csv", index=False)
    valid = frame.loc[frame.comparison_status == "comparable"].copy()
    groups = sorted(valid.cancer.unique())
    if len(groups) != 16:
        raise ValueError(f"expected 16 comparable groups, found {len(groups)}")
    pairs = [
        (
            "mean-rna-any-ihc",
            "rna_mean_ptpm",
            "prevalence_detected",
            "Mean RNA (native pTPM; log1p axis)",
            "Any IHC-positive patients (%)",
        ),
        (
            "rna-positive-any-ihc",
            "rna_prevalence",
            "prevalence_detected",
            "RNA samples with ≥1 pTPM (%)",
            "Any IHC-positive patients (%)",
        ),
        (
            "rna-positive-medium-high-ihc",
            "rna_prevalence",
            "prevalence_medium_high",
            "RNA samples with ≥1 pTPM (%)",
            "Medium/high IHC-positive patients (%)",
        ),
    ]
    for name, xcol, ycol, xlabel, ylabel in pairs:
        for page in range(2):
            fig, axes = plt.subplots(2, 4, figsize=(15, 8.2))
            for ax, cancer in zip(axes.flat, groups[page * 8 : (page + 1) * 8]):
                data = valid.loc[valid.cancer == cancer]
                x = data[xcol].astype(float)
                x = np.log1p(x) if xcol == "rna_mean_ptpm" else x * 100
                ax.scatter(
                    x,
                    data[ycol].astype(float) * 100,
                    s=9 + data.total.astype(float) * 1.5,
                    color="#166b86",
                    alpha=0.42,
                    linewidths=0.2,
                    edgecolors="white",
                )
                if xcol == "rna_mean_ptpm":
                    ticks = [0, 1, 10, 100, 1000, 10000]
                    ticks = [t for t in ticks if np.log1p(t) < max(x.max() * 1.1, 1.2)]
                    ax.set_xticks(np.log1p(ticks), [str(t) for t in ticks])
                    ax.set_xlim(left=-0.15)
                else:
                    ax.set_xlim(-3, 103)
                    ax.set_xticks([0, 25, 50, 75, 100])
                ax.set_ylim(-3, 103)
                ax.set_yticks([0, 25, 50, 75, 100])
                ax.set_title(f"{cancer.capitalize()}\nn = {len(data)} genes", fontsize=11)
                ax.grid(alpha=0.15)
                ax.set_axisbelow(True)
                ax.tick_params(labelsize=9)
                for spine in ["top", "right"]:
                    ax.spines[spine].set_visible(False)
            fig.suptitle(f"HPA 25.1 • RNA and cancer IHC • {page + 1}/2", fontsize=17, y=0.98)
            fig.supxlabel(xlabel, y=0.075, fontsize=12)
            fig.supylabel(ylabel, x=0.016, fontsize=12)
            fig.text(
                0.5,
                0.024,
                "Unpaired cohorts • IHC fractions count scored patients, not stained cells • Marker area increases with IHC patient count\n"
                "TCGA RNA only • Colorectal/lung/renal pooled by measured sample counts • Legacy 253-gene plotting scope; legacy HPA release unknown",
                ha="center",
                fontsize=9,
                color="#444444",
            )
            fig.subplots_adjust(
                left=0.065, right=0.99, top=0.86, bottom=0.16, hspace=0.42, wspace=0.27
            )
            fig.savefig(output / f"{name}-{page + 1}.png", dpi=150)
            plt.close(fig)
    audit = {
        "hpa_release": hpa_cancer_sources()["hpa_release"],
        "plot_scope": "gene IDs with IHC in the legacy tsarina cache (frozen plot-genes.csv)",
        "requested_genes": len(genes),
        "comparable_genes": valid.gene_id.nunique(),
        "comparable_groups": len(groups),
        "comparable_gene_group_rows": len(valid),
        "status_counts": {k: int(v) for k, v in frame.comparison_status.value_counts().items()},
        "genes_by_group": {k: int(v) for k, v in valid.groupby("cancer").size().items()},
        "cohort": "TCGA",
        "threshold_ptpm": 1,
        "cohorts_paired": False,
    }
    (output / "plot-audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gene-list", type=Path, default=Path("docs/audits/hpa-cancer/plot-genes.csv")
    )
    parser.add_argument("--output", type=Path, default=Path("docs/audits/hpa-cancer/figures"))
    args = parser.parse_args()
    plot(args.gene_list, args.output)


if __name__ == "__main__":
    main()

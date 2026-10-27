#!/usr/bin/env python3
"""Export scientific figures for the uniform quantitative candidate audit."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.hashsalt": "cta-quantitative-v1",
    }
)


def save(fig, directory, name: str) -> None:
    for ext in ["png", "svg"]:
        fig.savefig(
            directory / f"{name}.{ext}",
            dpi=180,
            bbox_inches="tight",
            metadata={"Date": None} if ext == "svg" else {},
        )
    plt.close(fig)


def plot(audit: Path):
    d = pd.read_csv(audit / "gene-results.csv.gz", low_memory=False)
    figures = audit / "figures"
    figures.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.3))
    for ax, ycol, positive, units, title in [
        (axes[0], "reproductive_ntpm", "reproductive_route", "nTPM", "Reproductive route"),
        (axes[1], "developmental_mean_tpm", "developmental_route", "TPM", "Developmental route"),
    ]:
        ok = d.adult_bulk_max_ntpm.notna() & d[ycol].notna()
        ax.scatter(
            d.loc[ok, "adult_bulk_max_ntpm"] + 0.05,
            d.loc[ok, ycol] + 0.05,
            s=6,
            c="#cbd5e1",
            alpha=0.4,
            rasterized=True,
            label="Other genes",
        )
        m = ok & d.selected & d[positive]
        ax.scatter(
            d.loc[m, "adult_bulk_max_ntpm"] + 0.05,
            d.loc[m, ycol] + 0.05,
            s=15,
            c="#c45228",
            alpha=0.8,
            label="Selected through this route",
        )
        ax.axvline(10.05, color="#334155", ls="--", lw=1)
        ax.axhline(5.05, color="#334155", ls=":", lw=1)
        ax.set(
            xscale="log",
            yscale="log",
            xlabel="Maximum adult bulk RNA (nTPM; +0.05 for plotting)",
            ylabel=f"Maximum qualifying context ({units}; +0.05)",
            title=title,
        )
        names = (
            ["PRAME", "MAGEA4", "CTAG1A", "CTAG1B", "CTAG2", "CT83", "SSX2", "XAGE1A", "XAGE1B"]
            if ycol == "reproductive_ntpm"
            else ["TBXT", "SALL4", "AFP", "CLDN6", "GPC3"]
        )
        texts = []
        for row in d[d.symbol.isin(names) & ok].itertuples():
            x = row.adult_bulk_max_ntpm + 0.05
            y = getattr(row, ycol) + 0.05
            ax.scatter([x], [y], s=35, facecolors="none", edgecolors="#172033", lw=1)
            texts.append(ax.text(x, y, row.symbol, fontsize=9, color="#172033"))
        adjust_text(
            texts,
            ax=ax,
            arrowprops={"arrowstyle": "-", "color": "#64748b", "lw": 0.5},
            expand=(1.2, 1.4),
        )
        ax.legend(loc="lower right", fontsize=8)
    axes[0].plot([0.05, 10000], [0.05, 10000], color="#334155", ls="-.", lw=0.7)
    fig.suptitle(
        "447 quantitative research candidates — no gene-specific admission", fontsize=15, y=1.02
    )
    fig.text(
        0.5,
        -0.035,
        "Panels show positive-context and bulk gates only. Cell, immune, IHC and coverage gates also apply. Native TPM and nTPM are not ratioed.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout()
    save(fig, figures, "expression-distribution")
    labels = [
        "Canonical universe",
        "Positive reproductive / developmental context",
        "Measured tumor expression",
        "Adult bulk at most 10 nTPM",
        "No detected adult IHC",
        "Supported adult cells below 30 nTPM",
        "Immune RNA below 5 native units",
        "Recurrent corneal RNA below 1 CP10K",
        "Complete required RNA panels",
    ]
    masks = [
        pd.Series(True, index=d.index),
        d.positive_context,
        d.tumor_positive,
        ~d.adult_bulk_exclusion,
        ~d.ihc_detected,
        ~d.adult_cell_exclusion,
        ~d.immune_exclusion,
        ~d.cornea_exclusion,
        d.rna_coverage_complete,
    ]
    counts = []
    m = pd.Series(True, index=d.index)
    for mask in masks:
        m &= mask
        counts.append(int(m.sum()))
    fig, ax = plt.subplots(figsize=(11, 5.3))
    y = np.arange(len(labels))
    ax.barh(y, counts, color=["#94a3b8"] * 8 + ["#c45228"])
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 22000)
    ax.set_xlabel("Genes remaining (sequential, non-overlapping removals)")
    for yi, n in zip(y, counts, strict=True):
        ax.text(n + 150, yi, f"{n:,}", va="center", fontsize=10)
    ax.set_title("Filtering all 20,097 genes with the same rules", pad=16)
    fig.tight_layout()
    save(fig, figures, "filtering-funnel")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("audit", type=Path)
    plot(p.parse_args().audit)

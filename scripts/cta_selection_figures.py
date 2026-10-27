#!/usr/bin/env python3
"""Export source overlaps and measured CTA distributions from the frozen audit."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LABELS = {
    "10.18632/oncotarget.21715": "da Silva 2017",
    "10.1038/ncomms10499": "Wang 2016",
    "10.1002/cam4.2223": "Chang 2019",
    "10.1038/s41467-021-22695-y": "Gong 2021",
    "10.1038/s41388-018-0357-2": "Bruggeman 2018",
    "10.1002/1878-0261.12900": "Jamin 2021",
    "10.1371/journal.pgen.1011734": "Loriot 2025",
    "10.1136/jitc-2023-007935": "Carter 2023",
    "10.1186/s12967-024-04918-0": "Seager 2024",
    "10.1038/s41467-020-19141-w": "Bradley 2020",
    "10.1158/0008-5472.CAN-16-0225": "Bai 2016",
    "10.1158/1078-0432.CCR-06-2340": "Garg 2007",
    "10.1016/s1074-7613(95)80053-0": "Boel 1995",
}
HIGHLIGHTS = {"PRAME": "#B94430", "MAGEA4": "#B17312", "CTAG1A": "#5058AC", "CTAG1B": "#2783AA"}


def build(audit: Path, out: Path) -> None:
    """Build the versioned audit tables, report and source provenance."""
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 13,
            "font.family": "DejaVu Sans",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
        }
    )

    def save(fig, name):
        fig.savefig(out / f"{name}.png", dpi=150, bbox_inches="tight", facecolor="white")
        fig.savefig(out / f"{name}.svg", bbox_inches="tight", facecolor="white")
        plt.close(fig)

    sources = pd.read_csv(audit / "nomination-sources.csv")
    overlaps = pd.read_csv(audit / "source-overlaps.csv").pivot(
        index="source_a", columns="source_b", values="overlap"
    )
    order = sources.doi.tolist()
    matrix = overlaps.loc[order, order].to_numpy()
    fig, ax = plt.subplots(figsize=(9.4, 7.3))
    ax.imshow(np.log1p(matrix), cmap="Blues", vmin=0)
    ax.set_xticks(
        range(len(order)), [LABELS[d].rsplit(" ", 1)[0] for d in order], rotation=60, ha="right"
    )
    ax.set_yticks(range(len(order)), [LABELS[d] for d in order])
    for i in range(len(order)):
        for j in range(len(order)):
            ax.text(
                j,
                i,
                str(matrix[i, j]),
                ha="center",
                va="center",
                fontsize=9,
                color="white" if matrix[i, j] > 80 else "#243A44",
            )
    ax.set_title("Shared mapped coding genes", loc="left", pad=16, weight="bold")
    fig.tight_layout()
    save(fig, "source-overlaps")

    d = pd.read_csv(audit / "plot-data.csv")
    fig, ax = plt.subplots(figsize=(10.2, 6.2))
    colors = {"Unsafe": "#D3D6D7", "Neither": "#9CAAB0", "CTA": "#158177"}
    for category in ["Unsafe", "Neither", "CTA"]:
        g = d[d.Category.eq(category)].dropna(
            subset=["Reproductive RNA (nTPM)", "adult_bulk_max_ntpm"]
        )
        ax.scatter(
            np.log10(1 + g["Reproductive RNA (nTPM)"]),
            np.log10(1 + g.adult_bulk_max_ntpm),
            s=17 if category != "CTA" else 28,
            color=colors[category],
            alpha=0.65,
            label=f"{category} (n={len(g):,} plotted)",
            linewidths=0,
        )
    offsets = {"PRAME": (15, 22), "MAGEA4": (15, 25), "CTAG1A": (6, 50), "CTAG1B": (-85, 22)}
    for name, color in HIGHLIGHTS.items():
        r = d[d.Gene.eq(name)].iloc[0]
        xy = (np.log10(1 + r["Reproductive RNA (nTPM)"]), np.log10(1 + r.adult_bulk_max_ntpm))
        ax.scatter(*xy, color=color, s=90, edgecolors="white", zorder=5)
        ax.annotate(
            name,
            xy,
            xytext=offsets[name],
            textcoords="offset points",
            fontsize=12,
            color=color,
            weight="bold",
            arrowprops={"arrowstyle": "-", "color": color},
        )
    ticks = [0, 1, 10, 100, 1000, 10000, 100000]
    ax.set_xticks(np.log10(1 + np.array(ticks[:-1])), [str(x) for x in ticks[:-1]])
    ax.set_yticks(np.log10(1 + np.array(ticks)), [str(x) for x in ticks])
    ax.set_xlim(-0.08, 4.3)
    ax.set_ylim(-0.08, 5.3)
    ax.set_xlabel("Reproductive bulk RNA maximum (nTPM, log[1+x] scale)")
    ax.set_ylabel("Adult normal bulk RNA maximum (nTPM)")
    ax.legend(loc="upper left", frameon=False)
    ax.grid(alpha=0.15)
    fig.tight_layout()
    save(fig, "bulk-rna-scatter")

    selected = d[d.Target & d.maximum_cell_cutoff_ratio.notna()].sort_values(
        "maximum_cell_cutoff_ratio"
    )
    x = np.log10(1 + selected.maximum_cell_cutoff_ratio.to_numpy())
    y = np.arange(1, len(x) + 1) / len(x)
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.step(x, y, where="post", color="#158177", lw=2.5)
    ax.axvline(np.log10(2), color="#B94430", linestyle="--", lw=1.5)
    ax.text(
        np.log10(2) + 0.02,
        1.015,
        "CTA cutoff",
        color="#B94430",
        fontsize=11,
    )
    for name, color in HIGHLIGHTS.items():
        r = selected[selected.Gene.eq(name)].iloc[0]
        ratio = r.maximum_cell_cutoff_ratio
        pct = float(selected.maximum_cell_cutoff_ratio.le(ratio).mean())
        xpos = np.log10(1 + ratio)
        ax.scatter(xpos, pct, s=70, color=color, zorder=5)
        offset = {"PRAME": (10, -20), "MAGEA4": (20, 25), "CTAG1A": (20, 5), "CTAG1B": (20, -24)}[
            name
        ]
        ax.annotate(
            f"{name}: {ratio:.3g}",
            (xpos, pct),
            xytext=offset,
            textcoords="offset points",
            color=color,
            fontsize=12,
            weight="bold",
        )
    ticks = [0, 1, 2, 5, 10, 20]
    ax.set_xticks(np.log10(1 + np.array(ticks)), [str(v) for v in ticks])
    ax.set_xlim(-0.03, max(x) * 1.08)
    ax.set_xlabel("Maximum supported cell RNA / its source-and-tissue CTA cutoff")
    ax.set_ylabel("Cumulative fraction of current CTA genes")
    ax.set_ylim(0, 1.06)
    ax.grid(alpha=0.15)
    fig.tight_layout()
    save(fig, "cell-cutoff-distribution")

    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    s = pd.read_csv(audit / "threshold-sensitivity.csv")
    for uf, group in s.groupby("upper_factor"):
        ax.plot(
            group.lower_factor,
            group.reference_review,
            marker="o",
            label=f"Reference cutoff ×{uf:g}",
        )
    ax.set_xlabel("CTA cutoff multiplier")
    ax.set_ylabel("Genes in the RNA/protein reference-review layer")
    ax.set_xticks([0.5, 1, 1.5])
    ax.legend(frameon=False)
    ax.grid(alpha=0.15)
    fig.tight_layout()
    save(fig, "independent-threshold-sensitivity")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.audit, args.output)

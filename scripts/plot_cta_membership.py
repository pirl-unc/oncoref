#!/usr/bin/env python3
"""Scientific figures for the biological membership audit; native units stay separate."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def build(audit):
    out = audit / "figures"
    out.mkdir(exist_ok=True)
    d = pd.read_csv(audit / "gene-results.csv.gz", low_memory=False)
    selected = d.loc[d.Member]
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
            "svg.hashsalt": "cta-membership-v1",
        }
    )

    def save(fig, name):
        fig.savefig(out / (name + ".png"), dpi=160, bbox_inches="tight", facecolor="white")
        fig.savefig(
            out / (name + ".svg"), bbox_inches="tight", facecolor="white", metadata={"Date": None}
        )
        plt.close(fig)

    tiers = ["Established", "Publication consensus", "Expression-supported"]
    counts = selected["Evidence tier"].value_counts().reindex(tiers)
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={"width_ratios": [1.3, 1]})
    a.barh(tiers[::-1], counts.iloc[::-1], color=["#8CB8B3", "#317D80", "#214C5C"])
    for i, n in enumerate(counts.iloc[::-1]):
        a.text(n + 12, i, str(n), va="center")
    a.set_xlim(0, counts.max() * 1.2)
    a.set_xlabel("Genes (tiers are disjoint)")
    a.set_title(
        f"{len(selected):,} biological candidates\n{len(d):,} canonical genes reprocessed",
        loc="left",
        weight="bold",
    )
    transition = [
        int((d.Member & d["Previous member"]).sum()),
        int((d.Member & ~d["Previous member"]).sum()),
        int((~d.Member & d["Previous member"]).sum()),
    ]
    b.barh(
        ["Retained from 215", "Added", "Removed"],
        transition,
        color=["#317D80", "#8CB8B3", "#C66B51"],
    )
    for i, n in enumerate(transition):
        b.text(n + 15, i, str(n), va="center")
    b.invert_yaxis()
    b.set_xlim(0, max(transition) * 1.2)
    b.set_xlabel("Genes")
    b.set_title("Change from the earlier shortlist", loc="left", weight="bold")
    fig.tight_layout()
    save(fig, "membership-and-changes")

    fig, ax = plt.subplots(figsize=(10.4, 6.2))
    for tier, color in zip(tiers[::-1], ["#8CB8B3", "#317D80", "#C66B51"]):
        rows = selected.loc[selected["Evidence tier"].eq(tier)]
        ax.scatter(
            rows["Normal breadth vital lower"],
            np.log10(rows["Core RNA (nTPM)"] + 1),
            s=16,
            alpha=0.5,
            color=color,
            label=tier,
        )
    highlights = {
        "PRAME": (0.20, 2.3),
        "MAGEA4": (0.07, 1.55),
        "CTAG1A": (0.08, 0.9),
        "CTAG1B": (0.08, 0.62),
        "MAGEA3": (0.08, 1.2),
        "TBXT": (0.075, 0.1),
        "GPC3": (0.61, 2.9),
    }
    for gene, offset in highlights.items():
        row = selected.loc[selected.Gene.eq(gene)].iloc[0]
        point = (row["Normal breadth vital lower"], np.log10(row["Core RNA (nTPM)"] + 1))
        color = {
            "Established": "#C66B51",
            "Publication consensus": "#317D80",
            "Expression-supported": "#8CB8B3",
        }[row["Evidence tier"]]
        ax.scatter(*point, s=45, color=color, zorder=4, edgecolors="black", linewidths=0.5)
        ax.annotate(
            gene,
            point,
            xytext=offset,
            textcoords="data",
            fontsize=10,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1},
            arrowprops={"arrowstyle": "-", "color": color, "lw": 0.6},
        )
    ax.set_xlabel("Vital-tissue expression breadth: observed lower bound (0–1)")
    ax.set_ylabel("log₁₀(1 + core reproductive RNA in nTPM)")
    ax.set_title("Membership does not erase normal expression", loc="left", weight="bold")
    ax.legend(frameon=False, loc="upper right")
    fig.text(
        0.1,
        0.01,
        "HPA v23 consensus. Missing organs increase the upper bound; this is not a safety score.\nReproductive axis = testis/ovary/placenta; developmental RNA is not included in this measurement.",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    save(fig, "expression-distribution")

    sources = pd.read_csv(audit / "source-counts.csv").sort_values("Mapped nominees")
    names = {
        "10.18632/oncotarget.21715": "da Silva 2017",
        "10.1038/ncomms10499": "Wang 2016",
        "10.1002/cam4.2223": "Chang 2019 (TGCT)",
        "10.1038/s41467-021-22695-y": "Gong 2021 (normal)",
        "10.1038/s41388-018-0357-2": "Bruggeman 2018",
        "10.1002/1878-0261.12900": "Jamin 2021",
        "10.1371/journal.pgen.1011734": "Loriot 2025",
        "10.1136/jitc-2023-007935": "Carter 2023",
        "10.1186/s12967-024-04918-0": "Seager 2024 (mixed)",
        "10.1038/s41467-020-19141-w": "Bradley 2020 (placenta)",
        "10.1158/0008-5472.CAN-16-0225": "Bai 2016",
        "10.1158/1078-0432.CCR-06-2340": "Garg 2007",
    }
    fig, ax = plt.subplots(figsize=(10.4, 6.2))
    labels = [names.get(doi, doi) for doi in sources.DOI]
    ax.barh(labels, sources["Mapped nominees"], color="#D9E4E7", label="All mapped nominees")
    ax.barh(labels, sources.Selected, color="#317D80", label="Selected in holistic union")
    for i, row in enumerate(sources.itertuples()):
        ax.text(row._2 + 10, i, f"{row.Selected}/{row._2}", va="center", fontsize=9)
    ax.set_xlim(0, sources["Mapped nominees"].max() * 1.22)
    ax.set_title(
        "Source coverage overlaps; counts do not sum to the union", loc="left", weight="bold"
    )
    ax.set_xlabel("Canonical genes")
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "source-coverage")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit", type=Path)
    build(parser.parse_args().audit)

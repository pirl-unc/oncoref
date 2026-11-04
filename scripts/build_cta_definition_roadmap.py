#!/usr/bin/env python3
"""Render the proposed CTA PR dependencies; this figure contains no result counts."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/plans/cta-definition-roadmap"


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "font.size": 12})
    fig, ax = plt.subplots(figsize=(14, 10))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    fig.patch.set_facecolor("#F7F8FA")

    def box(x, y, width, height, title, detail, *, support=False):
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                width,
                height,
                boxstyle="round,pad=0.009,rounding_size=0.012",
                linewidth=1.3,
                edgecolor="#A6AFBD" if support else "#4A6685",
                facecolor="#EDF0F4" if support else "#FFFFFF",
                zorder=2,
            )
        )
        ax.text(
            x + 0.017,
            y + height - 0.025,
            title,
            fontsize=13.5 if support else 15,
            fontweight="bold",
            color="#203A57",
            va="top",
            zorder=3,
        )
        ax.text(
            x + 0.017,
            y + height - 0.062,
            detail,
            fontsize=11.5,
            linespacing=1.65,
            color="#354555",
            va="top",
            zorder=3,
        )

    def arrow(start, end):
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=16,
                linewidth=1.5,
                color="#69829C",
                zorder=1,
            )
        )

    ax.text(
        0.055,
        0.96,
        "A consistent, explainable CTA definition",
        fontsize=24,
        fontweight="bold",
        color="#203A57",
        va="top",
    )
    ax.text(
        0.055,
        0.91,
        "Proposed PR dependencies  •  One decision contract, then evidence and sequence qualification",
        fontsize=12,
        color="#526477",
        va="top",
    )
    box(
        0.055,
        0.705,
        0.61,
        0.145,
        "1  Decision contract",
        "Definition, evidence status and sequence roles\n"
        "Held or incomplete does not imply allowed  •  #640 / #647 / #607",
    )
    box(
        0.71,
        0.705,
        0.235,
        0.145,
        "Independent  •  #626",
        "Isolate shard caches\nLand before large replays",
        support=True,
    )
    box(
        0.04,
        0.465,
        0.285,
        0.16,
        "2  Evidence reconciliation",
        "MAGEC2, LGALS13, PLCZ1\nMAGE-A12 / EPS8L2\n#633–635 / #654",
    )
    box(
        0.3575,
        0.465,
        0.285,
        0.16,
        "3  One tissue policy",
        "Current breadth uses ontology\nPreserve historical replay\n#637",
    )
    box(
        0.675,
        0.465,
        0.285,
        0.16,
        "4  Sequence qualification",
        "Locus, CDS and translation QC\nSERPINA2 / MROH5\n#657",
    )
    box(
        0.2,
        0.245,
        0.6,
        0.135,
        "5  Verified sequence bundle",
        "Policy, reference, hashes, occurrences and coverage\n"
        "Consumer migration and offline replay  •  #641 / #647",
    )
    box(
        0.2,
        0.045,
        0.6,
        0.135,
        "6  Public definition and report",
        "One receipt supplies narrative, flow, counts and evidence views\n"
        "HTML / PDF / SVG / CSV / JSON  •  #607",
    )
    for start, end in [
        ((0.2, 0.695), (0.1825, 0.637)),
        ((0.4, 0.695), (0.5, 0.637)),
        ((0.61, 0.695), (0.8175, 0.637)),
        ((0.1825, 0.454), (0.32, 0.392)),
        ((0.5, 0.454), (0.5, 0.392)),
        ((0.8175, 0.454), (0.68, 0.392)),
        ((0.5, 0.233), (0.5, 0.192)),
    ]:
        arrow(start, end)
    ax.text(
        0.055,
        0.015,
        "4 November 2026  •  Proposed work, not released policy or scientific results",
        fontsize=9.5,
        color="#64748B",
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT.with_suffix(".svg"), metadata={"Date": None})
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()

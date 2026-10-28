#!/usr/bin/env python3
"""Make exportable figures for the linked two-axis CTA threshold comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap

INK = "#172B3A"
TEAL = "#087E83"
RED = "#BF583C"
MUTED = "#617481"
PALE = "#DCE7EA"
BG = "#F7F9F8"


def setup() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "text.color": INK,
            "axes.labelcolor": INK,
            "xtick.color": MUTED,
            "ytick.color": INK,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.spines.left": False,
            "axes.spines.bottom": False,
            "figure.facecolor": BG,
            "axes.facecolor": BG,
            "svg.fonttype": "none",
            "svg.hashsalt": "oncoref-cta-geometric-20261027",
        }
    )


def title(fig, headline: str, subtitle: str) -> None:
    fig.text(0.045, 0.955, headline, fontsize=21, weight="bold", va="top")
    fig.text(0.045, 0.900, subtitle, fontsize=11, color=MUTED, va="top")


def footer(fig, note: str) -> None:
    fig.text(0.045, 0.025, note, fontsize=9, color=MUTED, va="bottom")


def save(fig, out: Path, name: str) -> None:
    for extension in ["png", "svg"]:
        fig.savefig(
            out / f"{name}.{extension}",
            dpi=200,
            facecolor=BG,
            metadata={"Date": None} if extension == "svg" else {},
        )
    plt.close(fig)


def threshold_figure(audit: Path, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.1))
    fig.subplots_adjust(left=0.13, right=0.97, top=0.73, bottom=0.22, wspace=0.48)
    title(
        fig,
        "One rule. Two doubling axes.",
        "Candidate ceiling = 2 × 2ᵗⁱˢˢᵘᵉ ʳᵃⁿᵏ × 2ᴵᴴᶜ ʳᵃⁿᵏ   ·   every step right or down doubles the allowance",
    )
    thresholds = pd.read_csv(audit / "cell-thresholds.csv")
    thresholds = thresholds[thresholds.assay.eq("hpa_single_cell_tissue")]
    matrices = [
        thresholds.pivot(index="category", columns="negative_ihc_grade", values=column)
        .loc[
            ["vital", "protected", "manageable"],
            ["Missing/uncertain", "Approved", "Supported", "Enhanced"],
        ]
        .to_numpy(dtype=int)
        for column in ["cta_cutoff", "reference_cutoff"]
    ]
    for ax, matrix, heading in zip(
        axes,
        matrices,
        ["Candidate ceiling  ·  pass at or below", "Stronger RNA review  ·  starts at"],
        strict=True,
    ):
        ax.imshow(np.log2(matrix), cmap="YlGnBu", vmin=0, vmax=8, aspect="auto")
        ax.set_xticks(
            range(4),
            [
                "Missing /\nunusable",
                "Approved\nnegative",
                "Supported\nnegative",
                "Enhanced\nnegative",
            ],
            fontsize=10,
        )
        ax.set_yticks(range(3), ["Vital", "Protected", "Manageable"])
        ax.tick_params(length=0, pad=12)
        ax.set_title(heading, loc="left", pad=19, fontsize=12, weight="bold")
        for i, j in np.ndindex(matrix.shape):
            ax.text(
                j,
                i,
                str(matrix[i, j]),
                ha="center",
                va="center",
                fontsize=23,
                weight="bold",
                color="white" if matrix[i, j] >= 16 else INK,
            )
        ax.set_xticks(np.arange(-0.5, 4, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, 3, 1), minor=True)
        ax.grid(which="minor", color=BG, linewidth=3)
        ax.tick_params(which="minor", length=0)
    footer(
        fig,
        "Shared provisional settings in native TPM/nTPM. Corneal CP10K cutoffs = table value ÷ 100.\nIHC must match the measured context; no independent single-cell adjustment. Stronger review does not assert protein or toxicity.",
    )
    save(fig, out, "threshold-grid")


def cascade_figure(audit: Path, out: Path) -> None:
    d = pd.read_csv(audit / "filter-cascade.csv")
    d = d[d["mode"].eq("two_axis_doubling")]
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 7.4), gridspec_kw={"width_ratios": [1, 1.08]})
    fig.subplots_adjust(left=0.22, right=0.96, top=0.76, bottom=0.15, wspace=0.90)
    title(
        fig,
        "20,097 genes → 319 candidates",
        "Every bar shows genes remaining. Red annotations show genes removed at that step only.",
    )
    labels = {
        "Reference genes": "All reference genes",
        "Reproductive or embryonic association": "Reproductive / embryonic",
        "Tumor sample ≥10 pTPM": "Tumor sample ≥10 pTPM",
        "No detected nonpermitted IHC": "Normal IHC gate",
        "Complete required RNA panels": "RNA coverage",
        "HPA tissue clusters": "HPA tissue clusters",
    }
    for ax, part, heading, maximum in [
        (axes[0], d.iloc[:5], "01  Association, tumor, protein & coverage", 20097),
        (axes[1], d.iloc[4:], "02  Tissue and immune RNA", 2948),
    ]:
        y = np.arange(len(part))
        ax.barh(
            y,
            part.remaining,
            color=[TEAL if n == 319 else PALE for n in part.remaining],
            height=0.56,
        )
        ax.set_yticks(y, [labels.get(s, s) for s in part["filter"]], fontsize=10.5)
        ax.invert_yaxis()
        ax.set_xlim(0, maximum * 1.30)
        ax.set_xticks([])
        ax.tick_params(axis="y", length=0, pad=9)
        ax.set_title(heading, loc="left", fontsize=11, weight="bold", pad=24)
        for k, r in enumerate(part.itertuples()):
            ax.text(
                r.remaining + maximum * 0.025,
                k,
                f"{r.remaining:,}",
                va="center",
                fontsize=11,
                weight="bold",
            )
            removed = 0 if k == 0 and ax is axes[1] else r.removed
            if removed:
                ax.text(
                    maximum * 1.29,
                    k,
                    f"−{removed:,}",
                    va="center",
                    ha="right",
                    color=RED,
                    fontsize=10,
                )
    footer(
        fig,
        "Tumor admission: ≥10 pTPM in one HPA TCGA/validation sample. Sequential counts depend on order; each gene is removed once.\nThe 2,948-gene RNA entry pool differs from the 2,960-gene evidence union, which also covers the alternative tumor recurrence rule.",
    )
    save(fig, out, "filter-cascade")


def choice_figure(audit: Path, out: Path) -> None:
    s = pd.read_csv(audit / "scenario-counts.csv")
    h = pd.read_csv(audit / "ihc-choice-counts.csv")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.8))
    fig.subplots_adjust(left=0.08, right=0.96, top=0.75, bottom=0.24, wspace=0.35)
    title(
        fig,
        "Simpler thresholds, 18 fewer candidates",
        "The lower Protected baseline drives the loss; greater IHC allowances do not add qualifying genes here.",
    )
    counts = s[s.tumor.eq("maximum10")].selected.to_numpy()
    for ax, values, names, heading in [
        (
            axes[0],
            counts,
            ["Previous\nlinked grid", "Double across\nIHC only", "Double across\nboth axes"],
            "Change one axis at a time",
        ),
        (
            axes[1],
            h.selected.to_numpy(),
            ["No IHC\nrelaxation", "Up to\nApproved", "Up to\nSupported", "Up to\nEnhanced"],
            "What does matching negative IHC add?",
        ),
    ]:
        x = np.arange(len(values))
        ax.bar(x, values, width=0.58, color=[PALE] * (len(values) - 1) + [TEAL])
        ax.set_xticks(x, names, fontsize=10)
        ax.set_ylim(0, 380)
        ax.set_ylabel("Selected genes")
        ax.set_title(heading, loc="left", weight="bold", fontsize=12, pad=15)
        ax.grid(axis="y", alpha=0.15)
        ax.set_axisbelow(True)
        for i, value in enumerate(values):
            ax.text(i, value + 8, str(value), ha="center", weight="bold", fontsize=17)
    footer(
        fig,
        "All panels use tumor ≥10 pTPM and all other gates. IHC panel caps the allowance earned by usable local negative evidence;\nit does not assume missing IHC is negative. 319 retained · 18 lost · 0 gained relative to the previous 337-gene result.",
    )
    save(fig, out, "policy-choices")


def overlap_figure(audit: Path, out: Path) -> None:
    d = pd.read_csv(audit / "filter-independent-effects.csv")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.5), sharey=True)
    fig.subplots_adjust(left=0.18, right=0.96, top=0.74, bottom=0.18, wspace=0.30)
    title(
        fig,
        "Which RNA filters do the work?",
        "Each filter is evaluated against the same 2,948 genes entering normal-RNA screening.",
    )
    for ax, col, heading, color in [
        (axes[0], "independent_failures", "Fail this filter, including overlap", MUTED),
        (axes[1], "uniquely_excluded", "Added if only this filter were removed", TEAL),
    ]:
        ax.barh(range(len(d)), d[col], height=0.54, color=color)
        ax.set_yticks(range(len(d)), d["filter"])
        ax.set_xlim(0, d[col].max() * 1.22)
        ax.set_title(heading, loc="left", fontsize=11, weight="bold", pad=18)
        ax.set_xlabel("Genes")
        ax.tick_params(axis="y", length=0)
        for i, value in enumerate(d[col]):
            ax.text(
                value + d[col].max() * 0.025,
                i,
                f"{value:,}",
                va="center",
                fontsize=11,
                weight="bold",
            )
    axes[0].invert_yaxis()
    footer(
        fig,
        "Independent failures overlap and must not be summed. The right panel changes one filter at a time, keeping all other gates.\nFor example, removing the HPA tissue-cluster filter alone would add 118 genes to the 319 selected.",
    )
    save(fig, out, "filter-overlap")


def named_figure(audit: Path, out: Path) -> None:
    m = pd.read_csv(audit / "named-normal-rna-margins.csv")
    g = pd.read_csv(audit / "named-gene-comparison.csv")
    names = [
        "PRAME",
        "CTAG2",
        "MAGEA4",
        "CTAG1A",
        "CTAG1B",
        "XAGE1A",
        "XAGE1B",
        "SSX2",
        "CT83",
        "TBXT",
        "SMC1B",
        "PAGE1",
    ]
    cols = [
        "Normal bulk",
        "HPA tissue clusters",
        "HPA immune",
        "Monaco immune",
        "DICE immune",
        "Adult cornea",
    ]
    ratios = m.pivot(index="symbol", columns="filter", values="ratio").reindex(
        index=names, columns=cols
    )
    arr = ratios.to_numpy()
    cmap = ListedColormap(["#B9DDD5", "#DEEEE9", "#F0F4D4", "#E5AE91", "#B94F3A"])
    cmap.set_bad("#E9EDEF")
    fig, ax = plt.subplots(figsize=(13, 8.7))
    fig.subplots_adjust(left=0.20, right=0.90, top=0.78, bottom=0.20)
    title(
        fig,
        "Where the named genes sit relative to the limits",
        "Each cell shows the largest supported RNA value ÷ its own applicable ceiling. Ratios >1 fail that RNA gate.",
    )
    ax.imshow(
        np.ma.masked_invalid(np.clip(arr, 0.0625, 16)),
        norm=BoundaryNorm([0, 0.25, 0.5, 1.000000000001, 4, 16], 5, clip=True),
        cmap=cmap,
        aspect="auto",
    )
    status = g.set_index("symbol").two_axis_doubling_maximum10_selected
    ax.set_yticks(range(len(names)), [n + ("  ✓" if status[n] else "  ×") for n in names])
    ax.set_xticks(
        range(len(cols)),
        ["Bulk", "HPA tissue\nclusters", "HPA\nimmune", "Monaco\nimmune", "DICE\nimmune", "Cornea"],
        fontsize=11,
    )
    ax.xaxis.tick_top()
    ax.tick_params(length=0, pad=10)
    ax.set_xticks(np.arange(-0.5, 6, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 12, 1), minor=True)
    ax.grid(which="minor", color=BG, linewidth=2)
    ax.tick_params(which="minor", length=0)
    for i, j in np.ndindex(arr.shape):
        value = arr[i, j]
        label = (
            "—"
            if np.isnan(value)
            else "0"
            if value == 0
            else "<0.01×"
            if value < 0.005
            else f"{value:.2f}×"
        )
        ax.text(
            j,
            i,
            label,
            ha="center",
            va="center",
            fontsize=11,
            color="white" if value > 4 else INK,
            weight="bold" if value > 1 else "normal",
        )
    footer(
        fig,
        "✓ selected · × excluded. Gray dash = no supported statistic, not a measured zero. Zero = supported source measurement of zero.\nOnly ratios to each source's policy ceiling are compared; native RNA units are retained. No missing panel is interpreted as absence.\nPRAME: heart fibroblasts 12/2 = 6×. CTAG2: bronchial basal cells 24.9/2 = 12.45×. SMC1B: DICE 3.50/4 = 0.87×.",
    )
    save(fig, out, "named-gene-margins")


def build(audit: Path) -> None:
    """Build the versioned audit tables, report and source provenance."""
    setup()
    out = audit / "figures"
    out.mkdir(exist_ok=True)
    threshold_figure(audit, out)
    cascade_figure(audit, out)
    choice_figure(audit, out)
    overlap_figure(audit, out)
    named_figure(audit, out)

    def h(p):
        return hashlib.sha256(p.read_bytes()).hexdigest()

    (out / "manifest.json").write_text(
        json.dumps(
            {
                "script_sha256": h(Path(__file__)),
                "audit_manifest_sha256": h(audit / "run-manifest.json"),
                "outputs": {
                    p.name: h(p) for p in sorted(out.iterdir()) if p.suffix in [".png", ".svg"]
                },
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit", type=Path)
    build(parser.parse_args().audit)

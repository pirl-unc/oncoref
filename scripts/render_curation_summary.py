#!/usr/bin/env python3
"""Draw the curation summary page: set overlap, funnels, thresholds and distributions.

Every number is read from the rebuilt catalog and the policy constants; nothing
here classifies a protein.
"""

from __future__ import annotations

import json
import sys
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyBboxPatch, Patch
from matplotlib.ticker import FuncFormatter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import (  # noqa: E402
    cta_literature,
    cta_selection_policy,
    target_filters,
    target_selection,
)

OUT = ROOT / "docs/audits/target-selection-20261110/rebuild"
FIG = OUT / "summary"
PEAKS = ROOT / "docs/audits/cta-selection-20261028/normal-peak-evidence.csv.gz"

UNSAFE, REVIEW, SEQUENCE, EXPRESSION = "#ac5552", "#b88a3c", "#8d7aa8", "#9aabb5"
CANDIDATE, CANCER, REPRODUCTIVE, PUBLISHED, GREY = (
    "#3f8f6b",
    "#3e8bba",
    "#438a69",
    "#c99a46",
    "#6f7f89",
)

#: Filtering ideas in the order drawn; each groups the evidence layers that implement it.
IDEAS = [
    ("Bulk tissue RNA above the tissue ceiling", {"Normal bulk"}),
    ("Cell-type cluster RNA (HPA single cell)", {"HPA tissue clusters"}),
    ("Sorted immune-cell RNA", {"HPA immune", "Monaco immune", "DICE immune"}),
    ("Rare-cell RNA established by donor", {"Established rare normal-cell RNA"}),
    ("Normal-tissue protein (IHC)", {"HPA normal protein", "Accepted normal-protein evidence"}),
    (
        "Curated or published normal expression",
        {"Reviewed normal-cell knowledge", "Published normal-tissue expression"},
    ),
    ("Adult cornea RNA", {"Adult cornea"}),
]
OTHER = [
    ("Coverage, identity, translation or unresolved evidence", "Other:Review", REVIEW),
    ("Nested, Unsafe-contained or read-through sequence", "Other:Sequence", SEQUENCE),
    ("No qualifying expression or bounds fail", "Other:Expression", EXPRESSION),
    ("Home tissue no panel resolves", "Other:Unscreened tissue", EXPRESSION),
]


def reasons(data):
    return data.Unsafe_reasons.fillna("").map(lambda v: set(v.split(" | ")) if v else set())


def funnel(data, start_label):
    """Sequential attribution: a protein leaves at the first idea that removes it."""
    flagged = reasons(data)
    left = pd.Series(True, index=data.index)
    rows = [{"step": start_label, "removed": 0, "remaining": len(data), "color": GREY}]
    for label, layers in IDEAS:
        any_layer = flagged.map(lambda s, layers=layers: bool(s & layers))
        hit = left & any_layer
        total = int(any_layer.sum())
        left &= ~hit
        rows.append(
            {
                "step": label,
                "removed": int(hit.sum()),
                "remaining": int(left.sum()),
                "flagged_in_total": total,
                "color": UNSAFE,
            }
        )
    for label, classification, color in OTHER:
        hit = left & data.classification.eq(classification)
        left &= ~hit
        rows.append(
            {"step": label, "removed": int(hit.sum()), "remaining": int(left.sum()), "color": color}
        )
    if not data[left].is_target_candidate.all() or left.sum() != data.is_target_candidate.sum():
        raise ValueError("Funnel does not end at the Candidate antigen proteins")
    return pd.DataFrame(rows)


def draw_funnel(table, title, path, zoom_from=1):
    """Single-panel attrition waterfall; uses a broken x-axis only when step 0 dwarfs step 1."""
    table = table[(table.removed > 0) | (table.index == 0)].reset_index(drop=True)
    # Append a visual-only final bar highlighting the surviving candidate union.
    display = pd.concat(
        [
            table,
            pd.DataFrame(
                [
                    {
                        "step": "Retained candidates (CA ∪ RA)",
                        "removed": 0,
                        "remaining": int(table.remaining.iloc[-1]),
                        "color": CANDIDATE,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    y = np.arange(len(display))[::-1]
    comma = FuncFormatter(lambda val, _: f"{int(val):,}")
    first_val = int(display.remaining.iloc[0])
    second_val = int(display.remaining.iloc[zoom_from]) if len(display) > zoom_from else first_val
    use_break = first_val > 0 and second_val / first_val < 0.25

    legend_items = [
        (GREY, "Initial pool"),
        (UNSAFE, "Unsafe (normal RNA / protein)"),
        (REVIEW, "Other:Review"),
        (SEQUENCE, "Other:Sequence"),
        (EXPRESSION, "Other:Expression / Unscreened"),
        (CANDIDATE, "Retained candidates (CA ∪ RA)"),
    ]

    if not use_break:
        fig, ax = plt.subplots(figsize=(10.5, 0.45 * len(display) + 2.0))
        ax.barh(y, display.remaining, color=display.color)
        for yi, row in zip(y, display.itertuples()):
            is_final = row.step == "Retained candidates (CA ∪ RA)"
            note = f"{row.remaining:,}" + (f"  (−{row.removed:,})" if row.removed else "")
            ax.text(
                row.remaining + first_val * 0.018,
                yi,
                note,
                va="center",
                fontsize=9.2,
                fontweight="bold" if is_final else "normal",
                color="#1f5e41" if is_final else "#1e2e38",
            )
        ax.set_xlim(0, first_val * 1.30)
        ax.xaxis.set_major_formatter(comma)
        ax.set_yticks(y, display.step, fontsize=9.2)
        ax.set_xlabel("Proteins remaining")
        ax.set_title(title, loc="left", fontsize=11.5, fontweight="bold", pad=10)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color="#e3e8eb", linewidth=0.6)
        ax.set_axisbelow(True)
    else:
        fig, (ax_l, ax_r) = plt.subplots(
            1,
            2,
            figsize=(11.2, 0.45 * len(display) + 2.0),
            sharey=True,
            gridspec_kw={"width_ratios": [5.0, 1.35], "wspace": 0.045},
        )
        left_max = max(second_val * 1.46, 3000)
        right_min = max(first_val - 1800, left_max + 1000)
        right_max = first_val + 2600

        ax_l.barh(y, display.remaining, color=display.color)
        ax_r.barh([y[0]], [first_val], color=[display.color.iloc[0]])

        for yi, row in zip(y[1:], display.iloc[1:].itertuples()):
            is_final = row.step == "Retained candidates (CA ∪ RA)"
            note = f"{row.remaining:,}" + (f"  (−{row.removed:,})" if row.removed else "")
            ax_l.text(
                row.remaining + left_max * 0.014,
                yi,
                note,
                va="center",
                fontsize=9.2,
                fontweight="bold" if is_final else "normal",
                color="#1f5e41" if is_final else "#1e2e38",
            )
        ax_r.text(
            first_val + 220,
            y[0],
            f"{first_val:,}",
            va="center",
            fontsize=9.2,
            fontweight="bold",
            color="#1e2e38",
        )

        ax_l.set_xlim(0, left_max)
        ax_r.set_xlim(right_min, right_max)
        ax_l.xaxis.set_major_formatter(comma)
        ax_r.set_xticks([20000])
        ax_r.xaxis.set_major_formatter(comma)
        ax_l.set_yticks(y, display.step, fontsize=9.2)
        ax_l.spines[["top", "right"]].set_visible(False)
        ax_r.spines[["top", "left", "right"]].set_visible(False)
        ax_r.tick_params(axis="y", length=0)
        for ax in (ax_l, ax_r):
            ax.grid(axis="x", color="#e3e8eb", linewidth=0.6)
            ax.set_axisbelow(True)

        # Draw diagonal break marks on the bottom axis spine.
        d = 0.015
        kw = {"color": "#334854", "lw": 1.2, "clip_on": False}
        ax_l.plot((1 - d, 1 + d), (-d, +d), transform=ax_l.transAxes, **kw)
        ax_r.plot((-d * 3.7, +d * 3.7), (-d, +d), transform=ax_r.transAxes, **kw)
        # Subtle break marks across the top bar itself.
        bar_y = (y[0] - ax_l.get_ylim()[0]) / (ax_l.get_ylim()[1] - ax_l.get_ylim()[0])
        ax_l.plot((1 - d, 1 + d), (bar_y - 0.03, bar_y + 0.03), transform=ax_l.transAxes, **kw)
        ax_r.plot(
            (-d * 3.7, +d * 3.7), (bar_y - 0.03, bar_y + 0.03), transform=ax_r.transAxes, **kw
        )

        ax_l.set_xlabel("Proteins remaining (broken x-axis above 3,000)", labelpad=8, fontsize=9.8)
        ax_l.set_title(title, loc="left", fontsize=11.5, fontweight="bold", pad=10)

    fig.legend(
        handles=[Patch(color=c, label=lbl) for c, lbl in legend_items],
        loc="lower center",
        bbox_to_anchor=(0.56, -0.02),
        ncol=3,
        frameon=False,
        fontsize=8.6,
    )
    if use_break:
        fig.subplots_adjust(bottom=0.17, top=0.92)
    else:
        fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return table


def venn(data, path):
    a, b, c = (
        data.is_cancer_antigen.astype(bool),
        data.is_reproductive_antigen.astype(bool),
        data.published_nomination.fillna(False).astype(bool),
    )
    n = lambda mask: f"{int(mask.sum()):,}"  # noqa: E731
    cta_total = int((a & b).sum())
    cta_pub = int((a & b & c).sum())
    cta_new = int((a & b & ~c).sum())
    ca_only = int((a & ~b).sum())
    ra_only = int((~a & b).sum())
    union_total = int((a | b).sum())
    pub_outside = int((~a & ~b & c).sum())

    fig, ax = plt.subplots(figsize=(9.2, 8.0))
    for (x, y), color in zip(
        [(-0.75, 0.35), (0.75, 0.35), (0, -0.9)], [CANCER, REPRODUCTIVE, PUBLISHED]
    ):
        ax.add_patch(Circle((x, y), 1.35, color=color, alpha=0.33, lw=1.6, ec=color))
    for (x, y), text in {
        (-1.45, 0.75): n(a & ~b & ~c),
        (1.45, 0.75): n(~a & b & ~c),
        (0, -1.65): n(~a & ~b & c),
        (0, 0.86): n(a & b & ~c),
        (-0.95, -0.5): n(a & ~b & c),
        (0.95, -0.5): n(~a & b & c),
        (0, -0.05): n(a & b & c),
    }.items():
        ax.text(
            x, y, text, ha="center", va="center", fontsize=15, fontweight="bold", color="#182730"
        )
    ax.annotate(
        f"CTA = CA ∩ RA = {cta_total:,}\n({cta_pub:,} published + {cta_new:,} newly identified)",
        xy=(0, 1.25),
        xytext=(0, 2.48),
        ha="center",
        va="center",
        fontsize=10.2,
        fontweight="bold",
        color="#1d5b42",
        bbox={
            "boxstyle": "round,pad=0.38",
            "facecolor": "#eef6f1",
            "edgecolor": "#438a69",
            "lw": 1.1,
        },
        arrowprops={"arrowstyle": "->", "color": "#438a69", "lw": 1.3},
    )
    ax.text(
        -2.25,
        1.78,
        f"Cancer Antigens\n{n(a)}",
        ha="center",
        va="center",
        color=CANCER,
        fontsize=11.5,
        fontweight="bold",
    )
    ax.text(
        2.25,
        1.78,
        f"Reproductive Antigens\n{n(b)}",
        ha="center",
        va="center",
        color=REPRODUCTIVE,
        fontsize=11.5,
        fontweight="bold",
    )
    ax.text(
        0,
        -2.58,
        f"Published CTA / oncofetal nominations\n{n(c)}",
        ha="center",
        color="#8a6a25",
        fontsize=11.5,
        fontweight="bold",
    )
    ax.text(
        0,
        -3.08,
        f"Candidate union CA ∪ RA = {union_total:,} ({cta_total:,} CTA · {ca_only:,} CA only · {ra_only:,} RA only)  |  "
        f"{pub_outside:,} published nominations excluded or held",
        ha="center",
        va="center",
        fontsize=9.2,
        color="#475b66",
    )
    ax.set_xlim(-3.05, 3.05)
    ax.set_ylim(-3.3, 2.95)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(
        f"Rule-derived sets against published nominations ({len(data):,} proteins; areas not to scale)",
        fontsize=11.5,
        fontweight="bold",
        pad=8,
    )
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return {
        "cta_published": cta_pub,
        "cta_new": cta_new,
        "published_outside": pub_outside,
    }


def layers(data, path):
    flagged = reasons(data)
    names = sorted({r for s in flagged for r in s})
    rows = pd.DataFrame(
        {
            "layer": names,
            "flagged": [int(flagged.map(lambda s, r=r: r in s).sum()) for r in names],
            "sole": [int(flagged.map(lambda s, r=r: s == {r}).sum()) for r in names],
        }
    ).sort_values("flagged")
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), sharey=True)
    for ax, column, title in zip(
        axes,
        ["flagged", "sole"],
        ["Unsafe proteins flagged by each layer", "…for which it is the only reason"],
    ):
        ax.barh(rows.layer, rows[column], color=UNSAFE)
        for i, value in enumerate(rows[column]):
            ax.text(value + rows[column].max() * 0.01, i, f"{value:,}", va="center", fontsize=9)
        ax.set_xlim(0, rows[column].max() * 1.15)
        ax.set_title(title, loc="left", fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return rows.sort_values("flagged", ascending=False)


def heat(ax, values, rows, columns, title, fmt, cmap="YlOrRd_r", log=True):
    shade = np.log2(values) if log else values
    ax.imshow(shade, cmap=cmap, aspect="auto", vmin=shade.min() - 1, vmax=shade.max() + 1)
    ax.set_xticks(range(len(columns)), columns)
    ax.set_yticks(range(len(rows)), rows)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            ax.text(j, i, fmt(values[i, j]), ha="center", va="center", fontsize=12)
    ax.set_title(title, loc="left", fontsize=10.5)
    ax.tick_params(length=0)


def thresholds(path):
    """Colour tables of the ceilings actually applied, read from the policy constants."""
    tiers = list(cta_selection_policy.TIER_BASES)
    grades = ["Missing/uncertain", "Approved", "Supported", "Enhanced"]
    bulk = np.array(
        [
            [
                cta_selection_policy.TIER_BASES[t] * 2.0 ** cta_selection_policy.IHC_RANKS[g]
                for g in grades
            ]
            for t in tiers
        ]
    )
    rare = np.array(
        [
            [f * 1.0, f * target_selection.RARE_CELL_SINGLE_STUDY_PCT]
            for f in (1.0, target_selection.RARE_CELL_PROTECTED_FACTOR)
        ]
    )
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.4), gridspec_kw={"width_ratios": [4, 4, 2.4]})
    heat(
        axes[0],
        bulk,
        [t.capitalize() for t in tiers],
        [f"{g}\nnegative IHC" for g in grades],
        "RNA ceiling, TPM or nTPM (bulk, cluster, immune)",
        lambda v: f"{v:g}",
    )
    heat(
        axes[1],
        bulk * cta_selection_policy.UNIT_SCALE["CP10K"],
        [t.capitalize() for t in tiers],
        [f"{g}\nnegative IHC" for g in grades],
        "Same ceiling for adult cornea, CP10K",
        lambda v: f"{v:g}",
    )
    heat(
        axes[2],
        rare,
        ["Vital", "Protected"],
        ["≥ 2 studies", "1 study"],
        "Rare-cell signal: % of a donor's cells\nwith ≥ 2 molecules",
        lambda v: f"{v:g}%",
        cmap="PuBu_r",
    )
    fig.suptitle(
        "Darker = stricter. Negative IHC in every routine context of the tissue raises the RNA "
        "ceiling; positive IHC is assessed in its own layer and earns no credit.",
        fontsize=9.5,
        y=0.02,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(bulk, index=tiers, columns=grades)


def applied(path):
    """How often each ceiling cell decided an exclusion, from the recorded peak observations."""
    peaks = pd.read_csv(
        PEAKS,
        usecols=[
            "assay",
            "unit",
            "selection_tier",
            "source",
            "tissue",
            "cell_type",
            "taxonomy_id",
            "stage",
            "category",
            "mapping_status",
            "negative_ihc_grade",
            "unsafe_evidence",
            "value",
        ],
        low_memory=False,
    )
    peaks = peaks.loc[peaks.assay.eq("Normal bulk")]
    current_tier = cta_selection_policy.compartment_tiers(peaks)
    peaks = peaks.loc[
        peaks.selection_tier.eq(current_tier) & current_tier.isin(cta_selection_policy.TIER_BASES)
    ].copy()
    peaks["grade"] = peaks.negative_ihc_grade.replace({"Missing": "Missing/uncertain"})
    tiers, grades = (
        ["vital", "protected", "manageable"],
        ["Missing/uncertain", "Approved", "Supported", "Enhanced"],
    )
    n = (
        peaks.groupby(["selection_tier", "grade"])
        .size()
        .unstack()
        .reindex(index=tiers, columns=grades)
        .fillna(0)
    )
    over = (
        peaks[peaks.unsafe_evidence]
        .groupby(["selection_tier", "grade"])
        .size()
        .unstack()
        .reindex(index=tiers, columns=grades)
        .fillna(0)
    )
    share = (over / n.where(n > 0)).to_numpy()
    fig, ax = plt.subplots(figsize=(8.2, 3.2))
    ax.imshow(np.nan_to_num(share), cmap="Reds", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(4), [f"{g}\nnegative IHC" for g in grades])
    ax.set_yticks(range(3), [t.capitalize() for t in tiers])
    for i in range(3):
        for j in range(4):
            text = "–" if not n.iloc[i, j] else f"{int(over.iloc[i, j]):,} / {int(n.iloc[i, j]):,}"
            ax.text(j, i, text, ha="center", va="center", fontsize=10)
    ax.tick_params(length=0)
    ax.set_title(
        "Recorded bulk RNA peaks with current tiers: Unsafe evidence / observations",
        loc="left",
        fontsize=10.5,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def groups(data):
    return [
        ("Unsafe", data.category.eq("Unsafe"), UNSAFE),
        ("Other", data.category.eq("Other"), EXPRESSION),
        ("Target candidate", data.is_target_candidate.astype(bool), CANDIDATE),
    ]


def distributions(data, path):
    """Per-protein values on each source scale; vertical lines are the rule thresholds."""
    base = cta_selection_policy.TIER_BASES
    tier_lines = [
        (base["vital"], "Vital"),
        (base["protected"], "Protected"),
        (base["manageable"], "Manageable"),
    ]
    scale = cta_selection_policy.UNIT_SCALE["CP10K"]
    panels = [
        (
            "healthy_maximum_ntpm",
            "HPA v23 nTPM: highest non-reproductive normal context",
            tier_lines,
            0.01,
        ),
        (
            "healthy_maximum_tpm",
            "Native TPM (GTEx v11, Monaco, DICE): highest normal context",
            tier_lines,
            0.01,
        ),
        (
            "healthy_maximum_cp10k",
            "Adult cornea CP10K: highest cell population",
            [(base["vital"] * scale, "Vital"), (base["protected"] * scale, "Protected")],
            0.001,
        ),
        (
            "reproductive_fold",
            "Reproductive observation / shared normal reference",
            [(target_selection.REPRODUCTIVE_FOLD, "RA gate")],
            0.01,
        ),
        (
            "developmental_all_lineage_max_tpm",
            "Embryo TPM: highest developmental lineage mean",
            [(1.0, "floor")],
            0.01,
        ),
        (
            "tumor_max_tpm",
            "Tumor TPM: highest qualified sample (proteins passing shared checks)",
            [(1.0, "floor")],
            0.01,
        ),
        (
            "disposable_max_tpm",
            "Disposable-tissue TPM (proteins passing shared checks)",
            [(target_selection.DISPOSABLE_CAP, "cap")],
            0.01,
        ),
        (
            "max_multi_molecule_pct",
            "Rare-cell panel: highest % of a donor population with ≥ 2 molecules",
            [(1.0, "replicated"), (target_selection.RARE_CELL_SINGLE_STUDY_PCT, "one study")],
            0.001,
        ),
    ]
    fig, axes = plt.subplots(4, 2, figsize=(14, 13))
    for ax, (column, title, lines, floor) in zip(axes.flat, panels):
        values = pd.to_numeric(data[column], errors="coerce")
        top = max(values.max(), max(v for v, _ in lines) * 3)
        bins = np.logspace(np.log10(floor), np.log10(top), 60)
        sets = groups(data)
        ax.hist(
            [values[m & values.notna()].clip(lower=floor) for _, m, _ in sets],
            bins=bins,
            stacked=True,
            color=[c for _, _, c in sets],
            label=[f"{name} ({int((m & values.notna()).sum()):,})" for name, m, _ in sets],
        )
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=10)
        for value, label in lines:
            ax.axvline(value, color="black", ls="--", lw=1)
            ax.text(
                value,
                0.55,
                f"{label} {value:g}",
                transform=ax.get_xaxis_transform(),
                va="center",
                ha="center",
                fontsize=8,
                rotation=90,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1},
            )
        ax.set_title(title, loc="left", fontsize=10)
        ax.set_ylabel("Proteins")
        ax.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Values at or below the left edge are drawn in the first bin. Ceilings shown without negative-IHC credit.",
        fontsize=9.5,
        y=0.005,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def _draw_enrichment_panel(ax, data, column, membership, title, color, callouts, floor=0.05):
    """Draw a single enrichment scatter plot against the shared non-thymic baseline."""
    d = data[data[column].notna()].copy()
    x = d[["normal_reference_tpm", "disposable_max_tpm"]].max(axis=1).clip(lower=floor)
    y = d[column].clip(lower=floor)
    held = ~d[membership] & d.classification.isin(["Other:Review", "Other:Unscreened tissue"])
    below = ~d[membership] & ~held
    for label, mask, shade, z in [
        ("Below threshold / non-member", below, "#b4c0c7", 1),
        ("Held (Review / Unscreened)", held, REVIEW, 2),
        ("Qualified membership", d[membership], color, 3),
    ]:
        ax.scatter(
            x[mask],
            y[mask],
            s=18 if z > 1 else 12,
            alpha=0.75 if z > 1 else 0.55,
            color=shade,
            zorder=z,
            label=f"{label} ({int(mask.sum())})",
        )
    grid = np.logspace(np.log10(floor), np.log10(max(100, x.max())), 100)
    ax.plot(
        grid,
        target_selection.ENRICHMENT_FOLD * np.maximum(grid, 1),
        "k:",
        lw=1.5,
        zorder=4,
        label="10 × max(N, D, 1 TPM)",
    )
    indexed = d.assign(_x=x, _y=y).set_index("protein_group_label")
    for gene, (dx, dy) in callouts.items():
        if gene not in indexed.index:
            continue
        row = indexed.loc[gene]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        is_held = row.classification in ("Other:Review", "Other:Unscreened tissue")
        ax.scatter(
            [row._x],
            [row._y],
            s=36,
            facecolors="none",
            edgecolors="#1d2f3a",
            linewidths=1.1,
            zorder=5,
        )
        ax.annotate(
            gene,
            xy=(row._x, row._y),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=8.8,
            fontweight="bold",
            color="#7d5918" if is_held else "#1a2d38",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.0},
            arrowprops={"arrowstyle": "-", "color": "#5f707a", "lw": 0.85},
            zorder=6,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="major", color="#dce3e7", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    ax.set_xlabel("Highest non-thymic, non-reproductive reference max(N, D) (TPM)", fontsize=10)
    ax.set_ylabel("Highest supported expression (TPM)", fontsize=10)
    ax.set_title(
        f"{title}: {len(d):,} measured proteins",
        loc="left",
        fontsize=11.5,
        fontweight="bold",
        pad=8,
    )
    ax.legend(fontsize=8.8, frameon=True, facecolor="white", edgecolor="#e0e5e8", loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)


def tumor_normal(data, path):
    """Both positive compartments against their identical nonreproductive reference.

    Writes both standalone full-width figures (`cancer-enrichment.png`,
    `reproductive-enrichment.png`) and the legacy combined figure (`path`).
    """
    floor = 0.05
    cancer_callouts = {
        "LALBA": (16, -2),
        "TGM4": (16, 4),
        "PRAME": (18, 12),
        "MAGEA4": (-48, 18),
        "NY-ESO-1": (-66, 4),
        "AMELX": (-54, 12),
        "GFY": (26, -14),
        "CTAG2": (24, -14),
    }
    repro_callouts = {
        "LALBA": (-48, 2),
        "SEMG2": (18, -4),
        "TGM4": (-46, 10),
        "PRAME": (24, 10),
        "CT83": (20, 22),
        "NY-ESO-1": (-64, -16),
        "CTAG2": (24, -12),
    }
    specs = [
        (
            "tumor_max_tpm",
            "is_cancer_antigen",
            "Cancer enrichment (CA)",
            CANCER,
            cancer_callouts,
            path.parent / "cancer-enrichment.png",
        ),
        (
            "reproductive_supported_tpm",
            "is_reproductive_antigen",
            "Reproductive and developmental enrichment (RA)",
            REPRODUCTIVE,
            repro_callouts,
            path.parent / "reproductive-enrichment.png",
        ),
    ]
    # Individual full-width figures for crisp readability in the whitepaper.
    for column, membership, title, color, callouts, out_path in specs:
        fig, ax = plt.subplots(figsize=(9.2, 6.2))
        _draw_enrichment_panel(ax, data, column, membership, title, color, callouts, floor=floor)
        fig.text(
            0.01,
            -0.02,
            "Values below 0.05 TPM are clipped to the axis floor. Gold points above the dotted boundary clear the 10× rule but are held in Other:Review or Other:Unscreened tissue.",
            ha="left",
            fontsize=8.5,
            color="#4e616b",
        )
        fig.tight_layout()
        fig.savefig(out_path, dpi=160, bbox_inches="tight")
        plt.close(fig)

    # Legacy combined side-by-side figure retained for summary.md.
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 6.5), sharex=True, sharey=True)
    for ax, (column, membership, title, color, callouts, _) in zip(axes, specs, strict=True):
        _draw_enrichment_panel(ax, data, column, membership, title, color, callouts, floor=floor)
    fig.suptitle(
        "Shared 10× enrichment rule across cancer and reproductive compartments",
        fontsize=12.5,
        fontweight="bold",
        y=0.99,
    )
    fig.text(
        0.5,
        0.005,
        "Values below 0.05 TPM are clipped to the left axis. Gold points above the dotted boundary clear the fold rule but are held for review or unscreened normal-home coverage.",
        ha="center",
        fontsize=8.8,
        color="#4e616b",
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.955))
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def rare(data, path):
    d = data[data.rare_cell_status.isin(["established", "unresolved", "dismissed", "below_screen"])]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.8), gridspec_kw={"width_ratios": [3, 2]})
    for status, color in [
        ("below_screen", "#c4ccd1"),
        ("unresolved", REVIEW),
        ("dismissed", GREY),
        ("established", UNSAFE),
    ]:
        part = d[d.rare_cell_status.eq(status)]
        axes[0].scatter(
            part.max_detection_pct.clip(lower=0.01),
            part.max_multi_molecule_pct.clip(lower=0.001),
            s=10,
            alpha=0.7,
            color=color,
            label=f"{status.replace('_', ' ')} ({len(part):,})",
        )
    for value, label in [
        (1.0, "1%: replicated"),
        (target_selection.RARE_CELL_SINGLE_STUDY_PCT, "5%: one study"),
    ]:
        axes[0].axhline(value, color="black", ls="--", lw=1)
        axes[0].text(0.011, value * 1.1, label, fontsize=8)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Highest % of a donor population with any molecule")
    axes[0].set_ylabel("Highest % with ≥ 2 molecules")
    axes[0].set_title("Rare-cell panel: frequency by donor", loc="left", fontsize=11)
    axes[0].legend(fontsize=8, frameon=False)
    basis = (
        d[d.rare_cell_status.ne("below_screen")]
        .groupby(["rare_cell_status", "rare_cell_basis"])
        .size()
    )
    basis = basis.sort_values()
    axes[1].barh(
        [textwrap.fill(b, 42) for _, b in basis.index],
        basis.values,
        color=[
            {"established": UNSAFE, "unresolved": REVIEW, "dismissed": GREY}[s]
            for s, _ in basis.index
        ],
    )
    for i, value in enumerate(basis.values):
        axes[1].text(value + 2, i, str(value), va="center", fontsize=9)
    axes[1].set_title(
        "Reason for established, unresolved\nor dismissed signals", loc="left", fontsize=11
    )
    axes[1].set_xlim(0, basis.max() * 1.15)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def stages(data, path):
    review = data[data.classification.eq("Other:Review")]
    coverage = ~review.required_coverage_complete.astype(bool)
    identity = ~coverage & ~review.identity_translation_pass.astype(bool)
    normal = ~coverage & ~identity & ~review.normal_screen_pass.astype(bool)
    rest = ~coverage & ~identity & ~normal
    why = pd.Series(
        {
            "Required normal coverage incomplete": coverage.sum(),
            "Identity or translation unresolved": identity.sum(),
            "Normal signal blocks admission (not Unsafe)": normal.sum(),
            "Other unresolved evidence or comparison": rest.sum(),
        }
    )
    candidates = data[data.is_target_candidate.astype(bool)]
    cancer = (
        candidates.cancer_routes.dropna().str.split(" | ", regex=False).explode().value_counts()
    )
    reproductive = (
        candidates.reproductive_routes.dropna()
        .str.split(" | ", regex=False)
        .explode()
        .value_counts()
    )
    fig, axes = plt.subplots(1, 3, figsize=(17, 3.8))
    for ax, series, color, title in [
        (axes[0], why, REVIEW, f"Other:Review ({len(review)}): first failing check"),
        (
            axes[1],
            cancer,
            CANCER,
            f"Cancer routes among {int(candidates.is_cancer_antigen.sum())} Cancer Antigens",
        ),
        (
            axes[2],
            reproductive,
            REPRODUCTIVE,
            f"Reproductive routes among {int(candidates.is_reproductive_antigen.sum())} Reproductive Antigens",
        ),
    ]:
        series = series.sort_values()
        ax.barh([s.replace("_", " ") for s in series.index], series.values, color=color)
        for i, value in enumerate(series.values):
            ax.text(value + series.max() * 0.01, i, f"{int(value):,}", va="center", fontsize=9)
        ax.set_xlim(0, series.max() * 1.15)
        ax.set_title(title, loc="left", fontsize=10.5)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return why


def markdown(table, columns):
    head = "| " + " | ".join(columns.values()) + " |\n|" + "---|" * len(columns) + "\n"
    return head + "\n".join(
        "| "
        + " | ".join(f"{v:,}" if isinstance(v, (int, np.integer)) else str(v) for v in row)
        + " |"
        for row in table[list(columns)].itertuples(index=False)
    )


def antigen_venn(data, path):
    """The two actual membership sets, separate from nomination provenance."""
    ca, ra = data.is_cancer_antigen, data.is_reproductive_antigen
    fig, ax = plt.subplots(figsize=(9, 5.4))
    ax.add_patch(Circle((-0.7, 0), 1.35, color=CANCER, alpha=0.28))
    ax.add_patch(Circle((0.7, 0), 1.35, color=REPRODUCTIVE, alpha=0.28))
    for x, text in [
        (-1.2, f"CA only\n{int((ca & ~ra).sum()):,}"),
        (0, f"CTA\n{int((ca & ra).sum()):,}"),
        (1.2, f"RA only\n{int((ra & ~ca).sum()):,}"),
    ]:
        ax.text(x, 0, text, ha="center", va="center", fontsize=14, weight="bold")
    ax.text(
        -1.2, 1.55, f"Cancer Antigens\n{int(ca.sum()):,}", ha="center", color=CANCER, fontsize=12
    )
    ax.text(
        1.2,
        1.55,
        f"Reproductive Antigens\n{int(ra.sum()):,}",
        ha="center",
        color=REPRODUCTIVE,
        fontsize=12,
    )
    ax.text(
        0,
        -1.65,
        f"Union: {int((ca | ra).sum()):,} candidate antigen proteins",
        ha="center",
        fontsize=14,
    )
    ax.text(
        0,
        -2,
        f"Outside both: {int(data.category.eq('Unsafe').sum()):,} Unsafe; "
        f"{int(data.category.eq('Other').sum()):,} Other",
        ha="center",
        fontsize=11,
    )
    ax.set(xlim=(-2.3, 2.3), ylim=(-2.2, 2.05), aspect="equal")
    ax.axis("off")
    ax.set_title("Current antigen sets (schematic; areas not to scale)", fontsize=14, pad=16)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def _save_both(fig, stem_path: Path, dpi: int = 170):
    with plt.rc_context({"svg.hashsalt": f"oncoref-{stem_path.name}-v1"}):
        svg_path = stem_path.with_suffix(".svg")
        fig.savefig(svg_path, bbox_inches="tight", metadata={"Date": None})
        svg_path.write_text(
            "\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n"
        )
        fig.savefig(stem_path.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def source_to_tiers_architecture(data, everything, stem_path):
    """Minimal visual architecture diagram: proteome -> safety screen -> CA/RA/CTA -> patient filter."""
    ca = data.is_cancer_antigen.astype(bool)
    ra = data.is_reproductive_antigen.astype(bool)
    n_total = len(data)
    n_unsafe = int(data.category.eq("Unsafe").sum())
    n_other = int(data.category.eq("Other").sum())
    n_cta = int((ca & ra).sum())
    n_ca_only = int((ca & ~ra).sum())
    n_ra_only = int((~ca & ra).sum())
    n_union = int((ca | ra).sum())
    n_eligible = int(
        (data[list(target_selection.SHARED_CHECKS)].all(axis=1) & ~data.category.eq("Unsafe")).sum()
    )

    fig, ax = plt.subplots(figsize=(10.2, 5.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.0)
    ax.axis("off")

    def node(x, y, w, h, title, count, fc, ec, title_c="#102a43", count_c="#102a43", note=""):
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.05,rounding_size=0.16",
                facecolor=fc,
                edgecolor=ec,
                linewidth=1.5,
            )
        )
        if count and note:
            ax.text(
                x + w / 2,
                y + h - 0.17,
                title,
                ha="center",
                va="top",
                fontsize=10.2,
                fontweight="bold",
                color=title_c,
            )
            ax.text(
                x + w / 2,
                y + h * 0.48,
                count,
                ha="center",
                va="center",
                fontsize=14.5,
                fontweight="bold",
                color=count_c,
            )
            ax.text(
                x + w / 2, y + 0.14, note, ha="center", va="bottom", fontsize=8.5, color="#486581"
            )
        elif count:
            ax.text(
                x + w / 2,
                y + h * 0.68,
                title,
                ha="center",
                va="center",
                fontsize=10.5,
                fontweight="bold",
                color=title_c,
            )
            ax.text(
                x + w / 2,
                y + h * 0.28,
                count,
                ha="center",
                va="center",
                fontsize=14.5,
                fontweight="bold",
                color=count_c,
            )
        else:
            ax.text(
                x + w / 2,
                y + h / 2,
                title,
                ha="center",
                va="center",
                fontsize=10.2,
                fontweight="bold",
                color=title_c,
            )

    def arrow(x0, y0, x1, y1, color="#627d98"):
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops={
                "arrowstyle": "-|>,head_width=0.28,head_length=0.40",
                "color": color,
                "lw": 1.6,
            },
        )

    # Top: Proteome
    node(3.15, 4.95, 3.70, 0.88, "Coding proteome", f"{n_total:,}", "#f4f7fa", "#9fb3c8")

    # Middle: Safety & eligibility screen branching to Unsafe (left) and Other (right)
    arrow(5.0, 4.90, 5.0, 4.35)
    node(
        3.15,
        3.22,
        3.70,
        1.06,
        "Normal safety & eligibility",
        f"{n_eligible:,}",
        "#eef4f8",
        "#5b7c92",
        note="Bulk · Single-cell · IHC",
    )
    arrow(3.10, 3.75, 2.52, 3.75, color=UNSAFE)
    node(
        0.20,
        3.22,
        2.25,
        1.06,
        "Unsafe",
        f"{n_unsafe:,}",
        "#fbeeed",
        UNSAFE,
        title_c="#8c2d2a",
        count_c="#8c2d2a",
        note="Healthy tissue",
    )

    arrow(6.90, 3.75, 7.48, 3.75, color=REVIEW)
    node(
        7.55,
        3.22,
        2.25,
        1.06,
        "Other",
        f"{n_other:,}",
        "#fbf5e7",
        REVIEW,
        title_c="#785416",
        count_c="#785416",
        note="Held / low fold",
    )

    # Enrichment step: symmetric 10x -> CA only, CTA, RA only
    arrow(4.10, 3.17, 1.75, 2.40, color=CANCER)
    arrow(5.00, 3.17, 5.00, 2.40, color=REPRODUCTIVE)
    arrow(5.90, 3.17, 8.25, 2.40, color="#6aa48b")
    ax.text(
        5.0,
        2.79,
        "≥ 10× healthy baseline",
        ha="center",
        va="center",
        fontsize=9.4,
        fontweight="bold",
        color="#243b53",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 2.5},
    )

    node(
        0.40,
        1.25,
        2.70,
        1.08,
        "Cancer Antigen only",
        f"{n_ca_only:,}",
        "#ebf4fa",
        CANCER,
        title_c="#1a567b",
        count_c="#1a567b",
        note="Tumor ≥ 10×",
    )
    node(
        3.55,
        1.25,
        2.90,
        1.08,
        "CTA  (CA ∩ RA)",
        f"{n_cta:,}",
        "#ebf5ef",
        REPRODUCTIVE,
        title_c="#1d573c",
        count_c="#1d573c",
        note="Tumor & reproductive",
    )
    node(
        6.90,
        1.25,
        2.70,
        1.08,
        "Reproductive only",
        f"{n_ra_only:,}",
        "#f0f7f4",
        "#7db19a",
        title_c="#266049",
        count_c="#266049",
        note="Reproductive ≥ 10×",
    )

    # Bottom: Optional patient-specific filtering with 3 clean choice pills
    arrow(5.0, 1.20, 5.0, 0.78, color="#627d98")
    ax.add_patch(
        FancyBboxPatch(
            (0.90, 0.08),
            8.20,
            0.64,
            boxstyle="round,pad=0.04,rounding_size=0.14",
            facecolor="#f7f9fb",
            edgecolor="#bcccdc",
            linewidth=1.2,
        )
    )
    ax.text(
        1.18,
        0.40,
        f"Patient-specific filter  ({n_union:,} candidates)",
        ha="left",
        va="center",
        fontsize=9.6,
        fontweight="bold",
        color="#243b53",
    )
    for px, lbl, fc, ec, tc in [
        (5.45, "Protect", "#fdf2f2", "#e4b8b6", "#8c2d2a"),
        (6.62, "Accept loss", "#fbf6e9", "#e4d1a3", "#785416"),
        (8.02, "Absent", "#eef5fa", "#b4d2e6", "#1a567b"),
    ]:
        w_pill = 1.02 if lbl != "Accept loss" else 1.26
        ax.add_patch(
            FancyBboxPatch(
                (px, 0.21),
                w_pill,
                0.38,
                boxstyle="round,pad=0.02,rounding_size=0.10",
                facecolor=fc,
                edgecolor=ec,
                linewidth=1.1,
            )
        )
        ax.text(
            px + w_pill / 2,
            0.40,
            lbl,
            ha="center",
            va="center",
            fontsize=8.8,
            fontweight="bold",
            color=tc,
        )

    fig.tight_layout()
    _save_both(fig, stem_path, dpi=165)


def _draw_host_matrix(ax, matrix, compartments, labels, colors, legend_labels):
    ax.imshow(
        matrix.T, cmap=matplotlib.colors.ListedColormap(colors), vmin=0, vmax=3, aspect="auto"
    )
    ax.set_yticks(
        range(len(compartments)), [c.replace("_", " ") for c in compartments], fontsize=9.5
    )
    ax.set_xticks(range(len(labels)), labels, rotation=32, ha="right", fontsize=9.0)
    ax.set_xticks(np.arange(-0.5, len(labels)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(compartments)), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)
    ax.tick_params(which="minor", length=0)
    for i in range(len(labels)):
        for j in range(len(compartments)):
            ax.text(
                i,
                j,
                ["—", "Protect", "Loss", "Absent"][matrix[i, j]],
                ha="center",
                va="center",
                fontsize=8.4,
                color="#2c3d47",
            )
    ax.set_title("Patient anatomy and organ-loss choices", loc="left", weight="bold", pad=12)
    ax.tick_params(length=0)
    used_states = sorted(set(matrix.ravel().tolist()))
    ax.legend(
        handles=[Patch(color=colors[idx], label=legend_labels[idx]) for idx in used_states],
        loc="upper left",
        bbox_to_anchor=(0, -0.24),
        ncol=len(used_states),
        frameon=False,
        fontsize=9.2,
    )


def _draw_host_outcomes(ax, profiles, labels):
    y = np.arange(len(profiles))
    left = np.zeros(len(profiles))
    total_cta = int(profiles.CTA.max())
    small_notes = [[] for _ in range(len(profiles))]
    for column, label, short, color in [
        ("CTA", "Retained for research screen", "retained", "#438a69"),
        ("CTA_unresolved", "Unresolved coverage / evidence", "unresolved", "#b88a3c"),
        ("CTA_excluded", "Excluded by patient-anatomy rule", "excluded", "#ac5552"),
    ]:
        values = profiles[column].to_numpy()
        ax.barh(y, values, left=left, label=label, color=color, height=0.68)
        for yi, offset, v in zip(y, left, values):
            if v >= 16:
                ax.text(
                    offset + v / 2,
                    yi,
                    str(v),
                    ha="center",
                    va="center",
                    fontsize=9.2,
                    fontweight="bold",
                    color="white",
                )
            elif v > 0:
                small_notes[yi].append(f"{v} {short}")
        left += values
    for yi, notes in enumerate(small_notes):
        if notes:
            ax.text(
                total_cta + 5,
                yi,
                " · ".join(notes),
                va="center",
                ha="left",
                fontsize=8.6,
                color="#4d5f69",
            )
    ax.set_xlim(0, total_cta * 1.24)
    ax.set_yticks(y, [s.replace("\n", " ") for s in labels], fontsize=9.4)
    ax.invert_yaxis()
    ax.set_xlabel("Biological CTAs (identity unchanged)", fontsize=10)
    ax.set_title(
        "Patient-profile screening outcome for each CTA", loc="left", weight="bold", pad=12
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper left", bbox_to_anchor=(0, -0.15), frameon=False, fontsize=9.2)


def cta_tiers_and_profiles(data, stem_path):
    """Host choices and outcomes; biological origins are not toxicity tiers.

    Writes both standalone full-width figures (`host-profiles-matrix`,
    `host-profiles-outcomes`) and the legacy combined figure (`stem_path`).
    """
    profiles = target_filters.profile_summary(data)
    states = ["Not protected", "Protect", "Accept loss", "Absent"]
    legend_labels = ["Not protected (—)", "Protect", "Accept loss", "Absent"]
    colors = ["#e4e8eb", "#efc7c4", "#e9d69c", "#c7ddea"]
    compartments = list(target_selection.REPRODUCTIVE_COMPARTMENTS)
    matrix = np.array(
        [
            [
                states.index(json.loads(row.host_states).get(c, "Not protected"))
                for c in compartments
            ]
            for row in profiles.itertuples()
        ]
    )
    labels = [r.replace("-patient", " anatomy").replace(" (", "\n(") for r in profiles.profile]
    labels[0] = "Catalog only"

    # Individual Figure A: Host-profile compartment choice matrix.
    fig_m, ax_m = plt.subplots(figsize=(10.2, 6.4))
    _draw_host_matrix(ax_m, matrix, compartments, labels, colors, legend_labels)
    fig_m.tight_layout()
    _save_both(fig_m, stem_path.parent / "host-profiles-matrix", dpi=165)

    # Individual Figure B: Host-screen outcomes for biological CTAs.
    fig_b, ax_b = plt.subplots(figsize=(9.8, 5.4))
    _draw_host_outcomes(ax_b, profiles, labels)
    fig_b.tight_layout()
    _save_both(fig_b, stem_path.parent / "host-profiles-outcomes", dpi=165)

    # Legacy combined figure retained for index.md / summary.md / report.html.
    fig, (grid, bars) = plt.subplots(
        1, 2, figsize=(17.2, 7.6), gridspec_kw={"width_ratios": [1.22, 1.05]}
    )
    _draw_host_matrix(grid, matrix, compartments, labels, colors, legend_labels)
    _draw_host_outcomes(bars, profiles, labels)
    fig.suptitle(
        "Protein origin and tolerated tissue loss are separate questions",
        fontsize=15.5,
        weight="bold",
    )
    fig.text(
        0.02,
        0.015,
        "No reproductive protection by default; select tissues to protect or accept loss explicitly via API/CLI. Presets record assumed anatomy; pregnancy coverage is incomplete.",
        fontsize=9.5,
        color="#475b66",
    )
    fig.tight_layout(rect=(0, 0.055, 1, 0.94))
    _save_both(fig, stem_path, dpi=165)


LITERATURE_SOURCE_META = {
    "10.18632/oncotarget.21715": (
        "CTA / cancer-germline",
        "da Silva et al. (2017) · pan-cancer CT",
    ),
    "10.1038/ncomms10499": (
        "CTA / cancer-germline",
        "Wang et al. (2016) · pan-cancer CT (19 types)",
    ),
    "10.1002/cam4.2223": ("CTA / cancer-germline", "Chang et al. (2019) · TGCT testis-specific"),
    "10.1038/s41388-018-0357-2": (
        "CTA / cancer-germline",
        "Bruggeman et al. (2018) · germ-cell in cancer",
    ),
    "10.1002/1878-0261.12900": ("CTA / cancer-germline", "Jamin et al. (2021) · microarray CT"),
    "10.1371/journal.pgen.1011734": (
        "CTA / cancer-germline",
        "Loriot et al. (2025) · cancer-germline",
    ),
    "10.1136/jitc-2023-007935": (
        "CTA / cancer-germline",
        "Carter et al. (2023) · thymus-absent pan-CTA",
    ),
    "10.1186/s12967-024-04918-0": (
        "CTA / cancer-germline",
        "Seager et al. (2024) · 17-target lung CTA",
    ),
    "10.1158/1078-0432.CCR-06-2340": (
        "CTA / cancer-germline",
        "Garg et al. (2007) · SPAG9 ovarian CTA",
    ),
    "10.1016/j.celrep.2023.113114": (
        "Oncofetal / DUX4 / embryo",
        "Smith et al. (2023) · DUX4 / ZGA / 8-cell",
    ),
    "10.1016/j.devcel.2019.06.011": (
        "Oncofetal / DUX4 / embryo",
        "Chew et al. (2019) · DUX4 cancer program",
    ),
    "10.1038/s41467-021-22695-y": (
        "Placental / reproductive",
        "Gong et al. (2021) · placenta & reproductive",
    ),
    "10.1371/journal.pone.0144155": (
        "Placental / reproductive",
        "Hirbo et al. (2015) · HPA placenta-enriched",
    ),
    "10.1038/s41467-020-19141-w": (
        "Placental / reproductive",
        "Bradley et al. (2020) · cancer-placenta (VGLL1)",
    ),
    "10.4161/cbt.12.3.16260": (
        "Tumor-vasculature / focused",
        "Sasaroli et al. (2011) · ovarian vasculature",
    ),
}
LITERATURE_CLASS_ORDER = [
    "CTA / cancer-germline",
    "Oncofetal / DUX4 / embryo",
    "Placental / reproductive",
    "Tumor-vasculature / focused",
]


def _classify_literature_source(source_id, roles):
    if source_id in LITERATURE_SOURCE_META:
        return LITERATURE_SOURCE_META[source_id][0]
    role_set = set(roles)
    if "cancer_embryo_program" in role_set:
        return "Oncofetal / DUX4 / embryo"
    if "catalog_nomination" in role_set:
        return "CTA / cancer-germline"
    if "reproductive_developmental" in role_set:
        return "Placental / reproductive"
    return "Tumor-vasculature / focused"


def literature_sources(data, fig_path, csv_path):
    """Map curated CTA/oncofetal/placental literature sources to current OncoRef outcomes."""
    att = cta_literature.attestations()
    cat = data.set_index("entity_id")
    per_src = att.drop_duplicates(["source_id", "entity_id"]).copy()
    mapped = per_src[per_src["entity_id"].isin(cat.index)].copy()
    mapped["classification"] = mapped["entity_id"].map(cat["classification"])
    mapped["category"] = mapped["entity_id"].map(cat["category"])
    mapped["is_nom"] = mapped["published_nomination"].fillna(False).astype(bool)

    rows = []
    for source_id, grp in mapped.groupby("source_id", sort=False):
        roles = sorted(
            {
                r
                for x in grp["source_roles"].dropna()
                for r in (x if isinstance(x, (list, tuple)) else [x])
            }
        )
        n = int(grp["entity_id"].nunique())
        cta = int(grp["classification"].eq("CTA").sum())
        ca_only = int(grp["classification"].eq("Cancer Antigen only").sum())
        ra_only = int(grp["classification"].eq("Reproductive Antigen only").sum())
        other = int(grp["category"].eq("Other").sum())
        unsafe = int(grp["category"].eq("Unsafe").sum())
        antigen_class = _classify_literature_source(source_id, roles)
        rows.append(
            {
                "source_id": source_id,
                "citation": grp["citation"].iloc[0],
                "antigen_class": antigen_class,
                "display_label": LITERATURE_SOURCE_META.get(
                    source_id, (antigen_class, grp["citation"].iloc[0])
                )[1],
                "source_roles": " | ".join(roles),
                "published_nomination_proteins": int(grp["is_nom"].sum()),
                "mapped_proteins": n,
                "cta": cta,
                "ca_only": ca_only,
                "ra_only": ra_only,
                "other": other,
                "unsafe": unsafe,
                "cta_pct": round(100.0 * cta / n, 1) if n else 0.0,
                "candidate_pct": round(100.0 * (cta + ca_only + ra_only) / n, 1) if n else 0.0,
                "source_url": grp["source_url"].iloc[0],
            }
        )
    full = pd.DataFrame(rows).sort_values(
        ["mapped_proteins", "cta", "citation"], ascending=[False, False, True]
    )
    full.to_csv(csv_path, index=False)

    major = full[full["mapped_proteins"].ge(5)].copy()
    small_ids = set(full[full["mapped_proteins"].lt(5)]["source_id"])
    small_unique = mapped[mapped["source_id"].isin(small_ids)].drop_duplicates("entity_id")
    small_n = len(small_unique)
    small_row = pd.DataFrame(
        [
            {
                "source_id": "focused_reviews",
                "citation": f"Focused reviews ({len(small_ids)} papers)",
                "antigen_class": "Tumor-vasculature / focused",
                "display_label": f"Focused reviews ({len(small_ids)} papers · unique proteins)",
                "mapped_proteins": small_n,
                "cta": int(small_unique["classification"].eq("CTA").sum()),
                "ca_only": int(small_unique["classification"].eq("Cancer Antigen only").sum()),
                "ra_only": int(
                    small_unique["classification"].eq("Reproductive Antigen only").sum()
                ),
                "other": int(small_unique["category"].eq("Other").sum()),
                "unsafe": int(small_unique["category"].eq("Unsafe").sum()),
            }
        ]
    )
    plot_df = pd.concat([major, small_row], ignore_index=True)
    plot_df["class_rank"] = plot_df["antigen_class"].map(
        {c: i for i, c in enumerate(LITERATURE_CLASS_ORDER)}
    )
    plot_df = plot_df.sort_values(
        ["class_rank", "mapped_proteins", "cta"], ascending=[True, False, False]
    ).reset_index(drop=True)

    categories = [
        ("cta", "CTA", REPRODUCTIVE),
        ("ca_only", "CA only", CANCER),
        ("ra_only", "RA only", "#94c4ad"),
        ("other", "Other (Review / Seq / Expr / Unscreened)", REVIEW),
        ("unsafe", "Unsafe (normal-tissue excluded)", UNSAFE),
    ]
    y = np.arange(len(plot_df))[::-1]
    totals = plot_df["mapped_proteins"].to_numpy(dtype=float)
    max_n = totals.max()
    comma = FuncFormatter(lambda val, _: f"{int(val):,}")

    def _draw_class_dividers(ax, x_badge):
        prev_class = None
        for idx, cls in enumerate(plot_df["antigen_class"]):
            if prev_class is not None and cls != prev_class:
                ax.axhline(y[idx] + 0.5, color="#495f6b", lw=1.0, ls="--", alpha=0.65)
            prev_class = cls
        for cls in LITERATURE_CLASS_ORDER:
            idxs = plot_df.index[plot_df["antigen_class"].eq(cls)].tolist()
            if not idxs:
                continue
            mid_y = (y[idxs[0]] + y[idxs[-1]]) / 2.0
            ax.text(
                x_badge,
                mid_y,
                cls,
                va="center",
                ha="left",
                fontsize=8.6,
                fontweight="bold",
                color="#2c414d",
                clip_on=False,
            )

    # Individual Figure A: Mapped protein counts per curated literature source.
    fig_c, ax_c = plt.subplots(figsize=(10.8, 0.45 * len(plot_df) + 2.3))
    left_c = np.zeros(len(plot_df))
    for col, label, color in categories:
        vals = plot_df[col].to_numpy(dtype=float)
        ax_c.barh(y, vals, left=left_c, color=color, label=label, height=0.70)
        left_c += vals
    for yi, row in zip(y, plot_df.itertuples()):
        ax_c.text(
            row.mapped_proteins + max_n * 0.015,
            yi,
            f"{int(row.mapped_proteins):,}",
            va="center",
            fontsize=8.8,
            color="#20323c",
        )
    _draw_class_dividers(ax_c, max_n * 1.15)
    ax_c.set_xlim(0, max_n * 1.13)
    ax_c.xaxis.set_major_formatter(comma)
    ax_c.set_xlabel("Mapped protein groups in source")
    ax_c.set_yticks(y, plot_df["display_label"], fontsize=9.2)
    ax_c.set_title(
        f"Curated literature sources: mapped protein count ({len(full)} sources)",
        loc="left",
        fontsize=11.5,
        fontweight="bold",
        pad=10,
    )
    ax_c.grid(axis="x", color="#e3e8eb", linewidth=0.6)
    ax_c.set_axisbelow(True)
    ax_c.spines[["top", "right"]].set_visible(False)
    handles, labels = ax_c.get_legend_handles_labels()
    fig_c.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.48, -0.02),
        ncol=3,
        frameon=False,
        fontsize=8.8,
    )
    fig_c.subplots_adjust(left=0.32, right=0.78, top=0.90, bottom=0.15)
    fig_c.savefig(fig_path.parent / "literature-sources-counts.png", dpi=160, bbox_inches="tight")
    plt.close(fig_c)

    # Individual Figure B: Share of mapped proteins retained vs excluded.
    fig_r, ax_r = plt.subplots(figsize=(10.8, 0.45 * len(plot_df) + 2.3))
    left_r = np.zeros(len(plot_df))
    for col, label, color in categories:
        vals = plot_df[col].to_numpy(dtype=float)
        pcts = np.where(totals > 0, 100.0 * vals / totals, 0.0)
        ax_r.barh(y, pcts, left=left_r, color=color, label=label, height=0.70)
        left_r += pcts
    for yi, row in zip(y, plot_df.itertuples()):
        cta_pct = 100.0 * row.cta / row.mapped_proteins if row.mapped_proteins else 0.0
        cand_pct = (
            100.0 * (row.cta + row.ca_only + row.ra_only) / row.mapped_proteins
            if row.mapped_proteins
            else 0.0
        )
        ax_r.text(
            101.8,
            yi,
            f"{cta_pct:.0f}% CTA ({cand_pct:.0f}% CA∪RA)",
            va="center",
            fontsize=8.7,
            fontweight="bold" if cta_pct >= 30 else "normal",
            color="#1f5e41" if cta_pct >= 30 else "#3e515c",
        )
    _draw_class_dividers(ax_r, 146.0)
    ax_r.set_xlim(0, 142)
    ax_r.set_xticks([0, 25, 50, 75, 100])
    ax_r.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax_r.set_xlabel("Share of mapped proteins by OncoRef disposition")
    ax_r.set_yticks(y, plot_df["display_label"], fontsize=9.2)
    ax_r.set_title(
        f"Curated literature sources: retention vs. normal-tissue exclusion ({len(full)} sources)",
        loc="left",
        fontsize=11.5,
        fontweight="bold",
        pad=10,
    )
    ax_r.grid(axis="x", color="#e3e8eb", linewidth=0.6)
    ax_r.set_axisbelow(True)
    ax_r.spines[["top", "right"]].set_visible(False)
    fig_r.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.48, -0.02),
        ncol=3,
        frameon=False,
        fontsize=8.8,
    )
    fig_r.subplots_adjust(left=0.32, right=0.78, top=0.90, bottom=0.15)
    fig_r.savefig(
        fig_path.parent / "literature-sources-retention.png", dpi=160, bbox_inches="tight"
    )
    plt.close(fig_r)

    # Legacy combined side-by-side figure retained for summary.md.
    fig, (ax_n, ax_p) = plt.subplots(
        1,
        2,
        figsize=(16.2, 0.44 * len(plot_df) + 2.4),
        sharey=True,
        gridspec_kw={"width_ratios": [1.28, 1.12], "wspace": 0.08},
    )
    left_n = np.zeros(len(plot_df))
    left_p = np.zeros(len(plot_df))
    for col, label, color in categories:
        vals = plot_df[col].to_numpy(dtype=float)
        pcts = np.where(totals > 0, 100.0 * vals / totals, 0.0)
        ax_n.barh(y, vals, left=left_n, color=color, label=label, height=0.72)
        ax_p.barh(y, pcts, left=left_p, color=color, height=0.72)
        left_n += vals
        left_p += pcts

    for yi, row in zip(y, plot_df.itertuples()):
        ax_n.text(
            row.mapped_proteins + max_n * 0.016,
            yi,
            f"{int(row.mapped_proteins):,}",
            va="center",
            fontsize=8.6,
            color="#20323c",
        )
        cta_pct = 100.0 * row.cta / row.mapped_proteins if row.mapped_proteins else 0.0
        cand_pct = (
            100.0 * (row.cta + row.ca_only + row.ra_only) / row.mapped_proteins
            if row.mapped_proteins
            else 0.0
        )
        summary_txt = f"{cta_pct:.0f}% CTA ({cand_pct:.0f}% CA∪RA)"
        ax_p.text(
            101.8,
            yi,
            summary_txt,
            va="center",
            fontsize=8.4,
            fontweight="bold" if cta_pct >= 30 else "normal",
            color="#1f5e41" if cta_pct >= 30 else "#3e515c",
        )

    prev_class = None
    for idx, cls in enumerate(plot_df["antigen_class"]):
        if prev_class is not None and cls != prev_class:
            sep_y = y[idx] + 0.5
            ax_n.axhline(sep_y, color="#495f6b", lw=1.0, ls="--", alpha=0.65)
            ax_p.axhline(sep_y, color="#495f6b", lw=1.0, ls="--", alpha=0.65)
        prev_class = cls

    for cls in LITERATURE_CLASS_ORDER:
        idxs = plot_df.index[plot_df["antigen_class"].eq(cls)].tolist()
        if not idxs:
            continue
        mid_y = (y[idxs[0]] + y[idxs[-1]]) / 2.0
        ax_p.text(
            156.0,
            mid_y,
            cls,
            va="center",
            ha="left",
            fontsize=8.5,
            fontweight="bold",
            color="#2c414d",
            bbox={"boxstyle": "round,pad=0.28", "facecolor": "#f1f4f6", "edgecolor": "#cfd8dd"},
            clip_on=False,
        )

    ax_n.set_xlim(0, max_n * 1.16)
    ax_n.xaxis.set_major_formatter(comma)
    ax_n.set_xlabel("Mapped protein groups in source")
    ax_n.set_yticks(y, plot_df["display_label"], fontsize=8.8)
    ax_n.set_title(
        "Mapped protein count by OncoRef disposition", loc="left", fontsize=11, fontweight="bold"
    )
    ax_n.grid(axis="x", color="#e3e8eb", linewidth=0.6)
    ax_n.set_axisbelow(True)

    ax_p.set_xlim(0, 152)
    ax_p.set_xticks([0, 25, 50, 75, 100])
    ax_p.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax_p.set_xlabel("Share of mapped proteins by disposition")
    ax_p.set_title(
        "Retention vs. normal-tissue exclusion", loc="left", fontsize=11, fontweight="bold"
    )
    ax_p.tick_params(axis="y", length=0)
    ax_p.grid(axis="x", color="#e3e8eb", linewidth=0.6)
    ax_p.set_axisbelow(True)

    for ax in (ax_n, ax_p):
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        f"Curated CTA, oncofetal, and placental literature sources in OncoRef ({len(full)} sources)",
        fontsize=12.5,
        fontweight="bold",
        y=0.985,
    )
    handles, labels = ax_n.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.48, -0.01),
        ncol=3,
        frameon=False,
        fontsize=8.8,
    )
    fig.subplots_adjust(left=0.26, right=0.81, top=0.90, bottom=0.14)
    fig.savefig(fig_path, dpi=155, bbox_inches="tight")
    plt.close(fig)
    return full


def classification_overview(data, path):
    """Write the minimal visual proteome-partition diagram to `path`."""
    source_to_tiers_architecture(data, None, path.with_suffix(""))


def normal_evidence_layers(data, path):
    """Single-panel visual showing what each non-bulk safety lens catches beyond bulk RNA."""
    flagged = reasons(data)
    bulk_mask = flagged.map(lambda s: "Normal bulk" in s)
    modalities = [
        ("Bulk tissue RNA\n(GTEx v11 + HPA v23)", {"Normal bulk"}),
        ("Single-cell clusters\n(HPA scRNA)", {"HPA tissue clusters"}),
        (
            "Sorted immune subsets\n(HPA / DICE / Monaco)",
            {"HPA immune", "Monaco immune", "DICE immune"},
        ),
        (
            "Normal protein staining\n(HPA IHC + curated)",
            {"HPA normal protein", "Accepted normal-protein evidence"},
        ),
        ("Donor rare-cell UMI\n(CELLxGENE + Siletti)", {"Established rare normal-cell RNA"}),
        (
            "Curated normal literature\n(912-protein rubric)",
            {"Reviewed normal-cell knowledge", "Published normal-tissue expression"},
        ),
        ("Adult cornea scRNA\n(Collin CP10K)", {"Adult cornea"}),
    ]
    rows = []
    # Sequential marginal additions beyond bulk RNA
    seen = pd.Series(False, index=data.index)
    for idx, (label, layer_set) in enumerate(modalities):
        hit = flagged.map(lambda s, ls=layer_set: bool(s & ls))
        added = hit & ~seen
        seen |= hit
        rows.append(
            {
                "label": label,
                "total": int(hit.sum()),
                "added": int(added.sum()),
                "beyond_bulk": int((hit & ~bulk_mask).sum()) if idx > 0 else 0,
            }
        )
    df = pd.DataFrame(rows)
    # Show the non-bulk modalities and how many Unsafe proteins each catches that bulk RNA missed!
    non_bulk = df.iloc[1:].copy().sort_values("beyond_bulk", ascending=True)
    y = np.arange(len(non_bulk))

    fig, ax = plt.subplots(figsize=(9.4, 4.8))
    bars = ax.barh(y, non_bulk["beyond_bulk"], color=UNSAFE, height=0.62, alpha=0.88)
    max_v = max(non_bulk["beyond_bulk"].max(), 1)
    for bar, row in zip(bars, non_bulk.itertuples(), strict=True):
        w = bar.get_width()
        ax.text(
            w + max_v * 0.02,
            bar.get_y() + bar.get_height() / 2,
            f"+{row.beyond_bulk:,} missed by bulk RNA   ({row.total:,} flagged overall)",
            va="center",
            ha="left",
            fontsize=8.8,
            color="#1f2933",
        )
    ax.set_yticks(y, non_bulk["label"], fontsize=9.2)
    ax.set_xlim(0, max_v * 1.92)
    ax.set_xlabel(
        "Unsafe proteins flagged despite passing bulk normal-tissue RNA ceilings", fontsize=9.8
    )
    ax.set_title(
        "Why bulk RNA alone is insufficient: Unsafe proteins caught by complementary assays",
        loc="left",
        fontsize=11.2,
        fontweight="bold",
        color="#102a43",
        pad=10,
    )
    ax.text(
        0.0,
        -0.19,
        f"Bulk tissue RNA alone flags {df.iloc[0]['total']:,} Unsafe proteins; the six complementary modalities above catch an additional "
        f"{int((seen & ~bulk_mask).sum()):,} Unsafe proteins (bars overlap where multiple non-bulk assays agree).",
        transform=ax.transAxes,
        fontsize=8.4,
        color="#486581",
    )
    ax.grid(axis="x", color="#e3e8eb", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def sex_restricted_antigens(data, fig_counts_path, fig_scatter_path, csv_path):
    """Partition Reproductive Antigens and CTAs by male vs. female adult organ restriction."""
    female_sex = (
        "ovary",
        "fallopian tube",
        "uterus",
        "endometrium",
        "cervix",
        "vagina",
        "lactating_breast",
    )
    female_all = (*female_sex, "resting_breast")
    male_sex = ("testis", "epididymis", "prostate", "seminal vesicle")
    male_all = (*male_sex, "resting_breast")

    ra = data[data.is_reproductive_antigen.astype(bool)].copy()
    f_all = target_filters.apply(
        target_filters.TargetFilter(
            host_profile="female-patient", protect_tissues=female_all, lactating=True
        ),
        ra,
    ).set_index("entity_id")
    f_sex = target_filters.apply(
        target_filters.TargetFilter(
            host_profile="female-patient", protect_tissues=female_sex, lactating=True
        ),
        ra,
    ).set_index("entity_id")
    m_all = target_filters.apply(
        target_filters.TargetFilter(host_profile="male-patient", protect_tissues=male_all),
        ra,
    ).set_index("entity_id")
    m_sex = target_filters.apply(
        target_filters.TargetFilter(host_profile="male-patient", protect_tissues=male_sex),
        ra,
    ).set_index("entity_id")

    exposures_df = pd.read_csv(ROOT / "oncoref/data/host-exposure.csv.gz").set_index("entity_id")

    def max_comp_tpm(eid, compartments):
        if eid not in exposures_df.index:
            return np.nan
        rec = json.loads(exposures_df.at[eid, "exposures"])
        vals = [
            rec[c]["max_tpm"]
            for c in compartments
            if c in rec and rec[c].get("max_tpm") is not None
        ]
        return float(max(vals)) if vals else np.nan

    rows = []
    for r in ra.itertuples():
        eid = r.entity_id
        fs_all = f_all.at[eid, "host_status"] == "eligible"
        ms_all = m_all.at[eid, "host_status"] == "eligible"
        fs_sex = f_sex.at[eid, "host_status"] == "eligible"
        ms_sex = m_sex.at[eid, "host_status"] == "eligible"
        if fs_all and not ms_all:
            sex_class = "Male-restricted"
            code = "CTA-M" if r.is_cta else "MRA-only"
        elif ms_all and not fs_all:
            sex_class = "Female-restricted"
            code = "CTA-F" if r.is_cta else "FRA-only"
        elif fs_all and ms_all:
            sex_class = "Embryo / placenta-only"
            code = "CTA-E" if r.is_cta else "ERA-only"
        else:
            sex_class = "Shared male & female"
            code = "CTA-MF" if r.is_cta else "SRA-only"
        rows.append(
            {
                "entity_id": eid,
                "protein_group_label": r.protein_group_label,
                "classification": r.classification,
                "is_cta": bool(r.is_cta),
                "sex_restriction_class": sex_class,
                "sex_restriction_code": code,
                "female_all_adult_eligible": fs_all,
                "male_all_adult_eligible": ms_all,
                "female_sex_organ_only_eligible": fs_sex,
                "male_sex_organ_only_eligible": ms_sex,
                "male_sex_organ_max_tpm": max_comp_tpm(eid, male_sex),
                "female_sex_organ_max_tpm": max_comp_tpm(eid, female_sex),
                "resting_breast_max_tpm": max_comp_tpm(eid, ("resting_breast",)),
                "tumor_max_tpm": r.tumor_max_tpm,
                "tumor_top_cancer_code": r.tumor_top_cancer_code,
                "reproductive_qualifying_compartments": r.reproductive_qualifying_compartments,
            }
        )
    table = pd.DataFrame(rows).sort_values(
        ["is_cta", "sex_restriction_class", "tumor_max_tpm", "protein_group_label"],
        ascending=[False, True, False, True],
    )
    table.to_csv(csv_path, index=False)

    # Figure 1: Stacked horizontal breakdown across All RA, CTA (CA ∩ RA), and RA only.
    class_specs = [
        ("Male-restricted", "Male-restricted (MRA / CTA-M · female-safe)", "#2a729e"),
        (
            "Embryo / placenta-only",
            "Embryo / placenta-only (ERA / CTA-E · adult-absent)",
            "#3f8f6b",
        ),
        ("Female-restricted", "Female-restricted (FRA / CTA-F · male-safe)", "#9a4c7a"),
        ("Shared male & female", "Shared male & female (SRA / CTA-MF)", "#c4a362"),
    ]
    subsets = [
        ("All Reproductive Antigens\n(RA)", table),
        ("Cancer–Testis Antigens\n(CTA = CA ∩ RA)", table[table["is_cta"]]),
        ("Reproductive only\n(RA only)", table[~table["is_cta"]]),
    ]
    fig, ax = plt.subplots(figsize=(10.2, 4.4))
    y = np.arange(len(subsets))[::-1]
    totals = np.array([len(df) for _, df in subsets], dtype=float)
    left = np.zeros(len(subsets))
    for key, label, color in class_specs:
        vals = np.array(
            [int(df["sex_restriction_class"].eq(key).sum()) for _, df in subsets], dtype=float
        )
        pcts = 100.0 * vals / totals
        ax.barh(y, pcts, left=left, color=color, label=label, height=0.56)
        for yi, off, v, p in zip(y, left, vals, pcts, strict=True):
            if p >= 7.5:
                ax.text(
                    off + p / 2,
                    yi,
                    f"{int(v):,}",
                    ha="center",
                    va="center",
                    fontsize=9.5,
                    fontweight="bold",
                    color="white",
                )
        left += pcts
    for yi, (_, df) in zip(y, subsets, strict=True):
        f_safe = int(df["female_all_adult_eligible"].sum())
        m_safe = int(df["male_all_adult_eligible"].sum())
        f_only = int(df["sex_restriction_class"].eq("Female-restricted").sum())
        ax.text(
            102.0,
            yi,
            f"n={len(df):,}  |  {f_safe} female-safe · {m_safe} male-safe ({f_only} FRA/CTA-F)",
            va="center",
            ha="left",
            fontsize=8.6,
            color="#243b53",
        )
    ax.set_xlim(0, 158)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax.set_yticks(y, [lbl for lbl, _ in subsets], fontsize=9.4)
    ax.set_xlabel(
        "Share of proteins under strict protection of all opposite-sex adult organs (including breast)"
    )
    ax.set_title(
        "Sex-restricted Reproductive Antigens (MRA / FRA / ERA) and cross-sex CTA subsets",
        loc="left",
        fontsize=11.2,
        fontweight="bold",
        color="#102a43",
        pad=10,
    )
    ax.grid(axis="x", color="#e3e8eb", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0, -0.22),
        ncol=2,
        frameon=False,
        fontsize=8.7,
    )
    fig.tight_layout()
    fig.savefig(fig_counts_path, dpi=160, bbox_inches="tight")
    plt.close(fig)

    # Figure 2: 2D male vs. female sex-organ exposure scatter across all 327 CTAs.
    ctas = table[table["is_cta"]].copy()
    floor = 0.05
    ctas["_mx"] = ctas["male_sex_organ_max_tpm"].fillna(0).clip(lower=floor)
    ctas["_fy"] = ctas["female_sex_organ_max_tpm"].fillna(0).clip(lower=floor)

    fig2, ax2 = plt.subplots(figsize=(9.4, 6.4))
    z_order = {
        "Shared male & female": 1,
        "Male-restricted": 2,
        "Embryo / placenta-only": 3,
        "Female-restricted": 4,
    }
    for key, label, color in class_specs:
        part = ctas[ctas["sex_restriction_class"].eq(key)]
        size = 42 if key == "Female-restricted" else 22
        ax2.scatter(
            part["_mx"],
            part["_fy"],
            s=size,
            alpha=0.82 if key != "Shared male & female" else 0.58,
            color=color,
            edgecolors="white" if key == "Female-restricted" else "none",
            linewidths=0.6,
            zorder=z_order[key],
            label=f"{label.split(' · ')[0]} ({len(part)})",
        )
    ax2.axvline(4.0, color="#627d98", ls=":", lw=1.2, zorder=0)
    ax2.axhline(4.0, color="#627d98", ls=":", lw=1.2, zorder=0)
    ax2.text(
        2200.0,
        0.07,
        "CTA-M: male-restricted\n(absent from female organs)",
        fontsize=8.4,
        color="#1a567b",
        fontweight="bold",
        va="bottom",
    )
    ax2.text(
        0.055,
        3.2,
        "CTA-E: embryo /\nplacenta-only",
        fontsize=8.4,
        color="#1d573c",
        fontweight="bold",
        va="top",
    )
    ax2.text(
        0.055,
        6500,
        "CTA-F: female-restricted\n(absent from male organs)",
        fontsize=8.4,
        color="#7a325d",
        fontweight="bold",
        va="top",
    )
    callouts = {
        "XAGE1A/B": (14, 14),
        "MAGEA4": (14, -22),
        "NY-ESO-1": (-28, 22),
        "TGM4": (-48, 14),
        "SEMG2": (14, -16),
        "TBXT": (14, 14),
        "MBD3L2": (16, -10),
        "ZP4": (16, 10),
        "ERVV-1": (-56, 16),
        "LALBA": (14, -8),
        "PRAME": (14, 10),
        "SCGB1D4": (-62, 10),
    }
    indexed = ctas.set_index("protein_group_label")
    for gene, (dx, dy) in callouts.items():
        if gene not in indexed.index:
            continue
        row = indexed.loc[gene]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        ax2.scatter(
            [row._mx],
            [row._fy],
            s=40,
            facecolors="none",
            edgecolors="#102a43",
            linewidths=1.0,
            zorder=5,
        )
        ax2.annotate(
            gene,
            xy=(row._mx, row._fy),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=8.6,
            fontweight="bold",
            color="#102a43",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.8},
            arrowprops={"arrowstyle": "-", "color": "#627d98", "lw": 0.8},
            zorder=6,
        )
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(0.035, 1.5e5)
    ax2.set_ylim(0.035, 1.5e5)
    ax2.grid(True, which="major", color="#dce3e7", linewidth=0.6, alpha=0.7)
    ax2.set_axisbelow(True)
    ax2.set_xlabel(
        "Max male sex-organ RNA (TPM: testis, epididymis, prostate, seminal vesicle)", fontsize=9.8
    )
    ax2.set_ylabel("Max female sex-organ RNA (TPM: ovary, tract, lactating breast)", fontsize=9.8)
    ax2.set_title(
        f"Male vs. female reproductive-organ exposure across all {len(ctas):,} CTAs",
        loc="left",
        fontsize=11.3,
        fontweight="bold",
        color="#102a43",
        pad=10,
    )
    ax2.legend(
        fontsize=8.3,
        frameon=True,
        facecolor="white",
        edgecolor="#dce3e7",
        loc="upper right",
    )
    ax2.spines[["top", "right"]].set_visible(False)
    fig2.text(
        0.01,
        -0.02,
        "Dotted lines mark the 4 TPM protected-organ screening threshold. Points in the bottom-right colored ochre (e.g., TGM4, SEMG2) have <4 TPM in female sex organs but ≥4 TPM in resting breast.",
        ha="left",
        fontsize=8.3,
        color="#486581",
    )
    fig2.tight_layout()
    fig2.savefig(fig_scatter_path, dpi=160, bbox_inches="tight")
    plt.close(fig2)
    return table


def build():
    data = pd.read_csv(OUT / "protein-decisions.csv.gz", low_memory=False)
    FIG.mkdir(exist_ok=True)
    classification_overview(data, FIG / "classification-overview.png")
    normal_evidence_layers(data, FIG / "normal-evidence-layers.png")
    sex_restricted_antigens(
        data,
        FIG / "sex-restricted-ra-cta.png",
        FIG / "sex-restricted-landscape.png",
        FIG / "sex-restricted-antigens.csv",
    )
    antigen_venn(data, FIG / "antigen-sets.png")
    counts = data.classification.value_counts()
    sets = venn(data, FIG / "venn.png")
    everything = draw_funnel(
        funnel(data, "All assessed proteins"),
        "Every assessed protein, by filtering idea",
        FIG / "funnel-all.png",
    )
    source_to_tiers_architecture(data, everything, FIG / "source-to-tiers-architecture")
    cta_tiers_and_profiles(data, FIG / "cta-tiers-and-profiles")
    profile_table = target_filters.profile_summary(data)
    profile_table.to_csv(FIG / "cta-profiles.csv", index=False)
    published = data[data.published_nomination.fillna(False).astype(bool)]
    nominated = draw_funnel(
        funnel(published, "Published CTA / oncofetal nominations"),
        "Published nominations through the same filters",
        FIG / "funnel-published.png",
    )
    if everything.removed.iloc[1:].eq(0).any():
        raise ValueError("A filtering idea removes nothing; reorder or merge it")
    everything.drop(columns="color").to_csv(FIG / "funnel-all.csv", index=False)
    nominated.drop(columns="color").to_csv(FIG / "funnel-published.csv", index=False)
    layer_table = layers(data, FIG / "unsafe-layers.png")
    layer_table.to_csv(FIG / "unsafe-layers.csv", index=False)
    ceilings = thresholds(FIG / "thresholds.png")
    applied(FIG / "thresholds-applied.png")
    distributions(data, FIG / "distributions.png")
    tumor_normal(data, FIG / "tumor-normal.png")
    rare(data, FIG / "rare-cell.png")
    why = stages(data, FIG / "stages.png")
    literature_sources(data, FIG / "literature-sources.png", FIG / "literature-sources.csv")
    fates = published.classification.value_counts()
    steps = {"step": "Filtering idea", "removed": "Removed here", "remaining": "Remaining"}
    (OUT / "summary.md").write_text(
        "# Curation summary\n\n"
        f"Policy `{target_selection.POLICY_ID}`, {len(data):,} human protein records. "
        f"**{counts['CTA']} CTA, {counts['Cancer Antigen only']} Cancer Antigen only, "
        f"{counts['Reproductive Antigen only']} Reproductive Antigen only.** These are outputs of "
        "screening rules, not validated tumor antigens. Every figure is drawn from the "
        "[rebuilt catalog](index.md); thresholds are read from the policy constants.\n\n"
        "## Catalog evidence, biological sets and host selection\n\n"
        "![Curation architecture](summary/source-to-tiers-architecture.png)\n\n"
        "Milk-program labels are inferred from expression, not provider-confirmed lactation. "
        "Negative signature evidence does not establish universal resting-state absence.\n\n"
        "## Host choices and screening outcomes\n\n"
        "![Host choices and outcomes](summary/cta-tiers-and-profiles.png)\n\n"
        "Reproductive protection is opt-in for every anatomy preset. Not protected, protect, "
        "accept loss and absent retain distinct meanings. Unresolved protected-tissue coverage fails the host screen, while biological "
        "CA/RA/CTA membership remains unchanged. [Rules and limitations](../../../host-selection.md).\n\n"
        + markdown(
            profile_table,
            {
                "profile": "Host profile",
                "protected_compartments": "Protected compartments",
                "CTA": "CTA",
                "CA_only": "CA only",
                "RA_only": "RA only",
                "Union": "CA ∪ RA",
            },
        )
        + "\n\n[Profile table CSV](summary/cta-profiles.csv)\n\n"
        "## Antigen membership\n\n![Antigen sets](summary/antigen-sets.png)\n\n"
        "CTA is the intersection. Reproductive-only proteins have no qualifying cancer evidence; "
        "the union name does not certify tumor targeting. Unsafe is the accepted-exclusion background; Unsafe|Other also subtracts unresolved and nonqualifying proteins.\n\n"
        "## Rule-derived sets and published nominations\n\n![Set overlap](summary/venn.png)\n\n"
        f"Of {counts['CTA']} CTA, {sets['cta_published']} are published nominations and "
        f"{sets['cta_new']} are not. {sets['published_outside']:,} published nominations are in "
        f"neither set: {fates.get('Unsafe', 0):,} of all {len(published):,} nominations are Unsafe, "
        f"{fates.get('Other:Review', 0)} Other:Review and {fates.get('Other:Sequence', 0)} "
        "Other:Sequence. A published nomination is an annotation; it never overrides an exclusion.\n\n"
        "### Curated literature sources by antigen class\n\n"
        "![Curated literature sources and dispositions](summary/literature-sources.png)\n\n"
        "[All curated literature sources and disposition counts (CSV)](summary/literature-sources.csv)\n\n"
        "## Funnel by filtering idea\n\n"
        "A protein leaves at the first idea that removes it, so later steps show only what each "
        "idea adds. The order is chosen so that every step removes something; the layers overlap "
        "heavily, as the next section shows.\n\n![All proteins](summary/funnel-all.png)\n\n"
        + markdown(everything, steps)
        + "\n\n![Published nominations](summary/funnel-published.png)\n\n"
        + markdown(nominated, steps)
        + "\n\nSteps that remove no published nomination are omitted from the second funnel. "
        "[All proteins](summary/funnel-all.csv) · [Published](summary/funnel-published.csv)\n\n"
        "## Normal-expression layers\n\n![Unsafe layers](summary/unsafe-layers.png)\n\n"
        + markdown(
            layer_table,
            {"layer": "Layer", "flagged": "Unsafe proteins flagged", "sole": "Only reason"},
        )
        + "\n\n## Thresholds\n\n![Thresholds](summary/thresholds.png)\n\n"
        "Ceilings are tier base × 2 per negative-IHC reliability step "
        f"(Vital {ceilings.iloc[0, 0]:g}, Protected {ceilings.iloc[1, 0]:g}, Manageable "
        f"{ceilings.iloc[2, 0]:g} without IHC credit). Exceeding a Vital or Protected ceiling is "
        f"Unsafe; Manageable is Unsafe only above {cta_selection_policy.MANAGEABLE_UNSAFE_FLOOR:g}. "
        "These are screening conventions, not measured injury thresholds.\n\n"
        "![Ceilings applied](summary/thresholds-applied.png)\n\n"
        "The second table counts recorded bulk-RNA peak observations, not proteins: one protein "
        "can contribute a peak per source.\n\n"
        "## Expression distributions and where the thresholds fall\n\n"
        "![Distributions](summary/distributions.png)\n\n"
        "One value per protein on each source's own scale; scales are never pooled. Counts use a "
        "compressed axis so the small candidate group stays visible.\n\n"
        "## Shared enrichment rule\n\nBoth CA and RA independently require 10× max(N, D, 1 TPM). Both normal compartments exclude thymus. The cancer panel below shows the tumor side.\n\n![Tumor against normal](summary/tumor-normal.png)\n\n"
        "## Rare-cell RNA\n\n![Rare-cell grading](summary/rare-cell.png)\n\n"
        "[Donor-level evidence](../../normal-cell-atlas-20261110/rare-cell/index.md)\n\n"
        "## Inside the later stages\n\n![Stage detail](summary/stages.png)\n\n"
        + "\n".join(f"- {label}: {int(value)}" for label, value in why.items())
        + "\n\nA protein can use several routes, so route counts exceed set sizes.\n"
    )


if __name__ == "__main__":
    build()

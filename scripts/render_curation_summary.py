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


def tumor_normal(data, path):
    """Both positive compartments against their identical nonreproductive reference."""
    floor = 0.05
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 6.5), sharex=True, sharey=True)
    cancer_callouts = {
        "LALBA": (16, -2),
        "TGM4": (16, 4),
        "PRAME": (18, 12),
        "MAGEA4": (-44, 18),
        "NY-ESO-1": (-62, 4),
        "AMELX": (-52, 12),
        "GFY": (26, -14),
        "CTAG2": (24, -14),
    }
    repro_callouts = {
        "LALBA": (-48, 2),
        "SEMG2": (18, -4),
        "TGM4": (-46, 10),
        "PRAME": (24, 10),
        "CT83": (20, 22),
        "NY-ESO-1": (-62, -16),
        "CTAG2": (24, -12),
    }
    for ax, column, membership, title, color, callouts in [
        (axes[0], "tumor_max_tpm", "is_cancer_antigen", "Cancer (CA)", CANCER, cancer_callouts),
        (
            axes[1],
            "reproductive_supported_tpm",
            "is_reproductive_antigen",
            "Reproductive / developmental (RA)",
            REPRODUCTIVE,
            repro_callouts,
        ),
    ]:
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
                s=13 if z > 1 else 9,
                alpha=0.72 if z > 1 else 0.55,
                color=shade,
                zorder=z,
                label=f"{label} ({int(mask.sum())})",
            )
        grid = np.logspace(np.log10(floor), np.log10(max(100, x.max())), 100)
        ax.plot(
            grid,
            target_selection.ENRICHMENT_FOLD * np.maximum(grid, 1),
            "k:",
            lw=1.4,
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
                s=28,
                facecolors="none",
                edgecolors="#1d2f3a",
                linewidths=1.0,
                zorder=5,
            )
            ax.annotate(
                gene,
                xy=(row._x, row._y),
                xytext=(dx, dy),
                textcoords="offset points",
                fontsize=8.2,
                fontweight="bold",
                color="#7d5918" if is_held else "#1a2d38",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.8},
                arrowprops={"arrowstyle": "-", "color": "#5f707a", "lw": 0.8},
                zorder=6,
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(True, which="major", color="#dce3e7", linewidth=0.6, alpha=0.7)
        ax.set_axisbelow(True)
        ax.set_xlabel("Highest non-thymic, non-reproductive reference max(N, D) (TPM)")
        ax.set_ylabel("Highest supported expression (TPM)")
        ax.set_title(
            f"{title}: {len(d):,} measured proteins", loc="left", fontsize=11, fontweight="bold"
        )
        ax.legend(
            fontsize=8.2, frameon=True, facecolor="white", edgecolor="#e0e5e8", loc="lower right"
        )
        ax.spines[["top", "right"]].set_visible(False)
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
    """Three separate decisions, with measured counts and explicitly scoped evidence."""
    ca, ra = data.is_cancer_antigen, data.is_reproductive_antigen
    eligible = data[list(target_selection.SHARED_CHECKS)].all(axis=1) & ~data.category.eq("Unsafe")
    boxes = [
        (
            "1  SHARED CATALOG SCREEN",
            f"{len(data):,} protein groups → {int(eligible.sum()):,} shared-eligible",
            "HPA / GTEx bulk RNA · cell clusters · donor-resolved rare cells\nImmune RNA · normal IHC · reviewed findings · identity / sequence\nAccepted normal exclusions → Unsafe; gaps and technical holds → Other",
        ),
        (
            "2  BIOLOGICAL MEMBERSHIP",
            f"CA {int(ca.sum())} · RA {int(ra.sum())} · CTA {int((ca & ra).sum())} · union {int((ca | ra).sum())}",
            "QC-passing tumor samples and reproductive / developmental RNA\nBoth require ≥10 × max(N, D, 1 TPM). N and D exclude thymus.\nCTA = CA ∩ RA; origin compartments describe evidence, not tolerable injury",
        ),
        (
            "3  OPTIONAL HOST SELECTION",
            "Protect  /  Accept loss  /  Absent",
            "Separate exposure evidence, including zeros and missing coverage\nProtected TPM and nTPM screened separately; unresolved protein findings held\nHost rejection preserves biological identity and expands peptide background",
        ),
    ]
    fig, ax = plt.subplots(figsize=(13, 8.2))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    for i, (title, count, body) in enumerate(boxes):
        y = 0.96 - i * 0.31
        ax.add_patch(
            FancyBboxPatch(
                (0.02, y - 0.265),
                0.96,
                0.25,
                boxstyle="round,pad=.01",
                facecolor=["#edf1f4", "#e8f1ec", "#f6efe0"][i],
                edgecolor="#b8c3c8",
            )
        )
        ax.text(0.045, y - 0.045, title, fontsize=12, weight="bold", va="top")
        ax.text(0.045, y - 0.09, count, fontsize=16, weight="bold", va="top")
        ax.text(0.045, y - 0.145, body, fontsize=11, va="top", linespacing=1.5)
        if i < 2:
            ax.annotate(
                "",
                xy=(0.5, y - 0.31),
                xytext=(0.5, y - 0.275),
                arrowprops={"arrowstyle": "->", "color": "#536671", "lw": 2},
            )
    ax.text(
        0.02,
        0.005,
        "Research screening only. RNA enrichment does not establish protein presentation, recognition or clinical safety.",
        fontsize=10,
        color="#536671",
    )
    fig.tight_layout()
    _save_both(fig, stem_path, dpi=165)


def cta_tiers_and_profiles(data, stem_path):
    """Host choices and outcomes; biological origins are not toxicity tiers."""
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
    fig, (grid, bars) = plt.subplots(
        1, 2, figsize=(17.2, 7.6), gridspec_kw={"width_ratios": [1.22, 1.05]}
    )
    grid.imshow(
        matrix.T, cmap=matplotlib.colors.ListedColormap(colors), vmin=0, vmax=3, aspect="auto"
    )
    grid.set_yticks(
        range(len(compartments)), [c.replace("_", " ") for c in compartments], fontsize=9
    )
    grid.set_xticks(range(len(labels)), labels, rotation=32, ha="right", fontsize=8.5)
    grid.set_xticks(np.arange(-0.5, len(labels)), minor=True)
    grid.set_yticks(np.arange(-0.5, len(compartments)), minor=True)
    grid.grid(which="minor", color="white", linewidth=1.0)
    grid.tick_params(which="minor", length=0)
    for i in range(len(labels)):
        for j in range(len(compartments)):
            grid.text(
                i,
                j,
                ["—", "Protect", "Loss", "Absent"][matrix[i, j]],
                ha="center",
                va="center",
                fontsize=8,
                color="#2c3d47",
            )
    grid.set_title("Declared anatomy and tissue-loss choices", loc="left", weight="bold", pad=12)
    grid.tick_params(length=0)
    used_states = sorted(set(matrix.ravel().tolist()))
    grid.legend(
        handles=[Patch(color=colors[idx], label=legend_labels[idx]) for idx in used_states],
        loc="upper left",
        bbox_to_anchor=(0, -0.24),
        ncol=len(used_states),
        frameon=False,
        fontsize=9,
    )
    y = np.arange(len(profiles))
    left = np.zeros(len(profiles))
    total_cta = int(profiles.CTA.max())
    small_notes = [[] for _ in range(len(profiles))]
    for column, label, short, color in [
        ("CTA", "Retained for research screen", "retained", "#438a69"),
        ("CTA_unresolved", "Unresolved coverage / evidence", "unresolved", "#b88a3c"),
        ("CTA_excluded", "Fails a host requirement", "excluded", "#ac5552"),
    ]:
        values = profiles[column].to_numpy()
        bars.barh(y, values, left=left, label=label, color=color)
        for yi, offset, v in zip(y, left, values):
            if v >= 16:
                bars.text(
                    offset + v / 2,
                    yi,
                    str(v),
                    ha="center",
                    va="center",
                    fontsize=9,
                    fontweight="bold",
                    color="white",
                )
            elif v > 0:
                small_notes[yi].append(f"{v} {short}")
        left += values
    for yi, notes in enumerate(small_notes):
        if notes:
            bars.text(
                total_cta + 5,
                yi,
                " · ".join(notes),
                va="center",
                ha="left",
                fontsize=8.2,
                color="#4d5f69",
            )
    bars.set_xlim(0, total_cta * 1.24)
    bars.set_yticks(y, [s.replace("\n", " ") for s in labels], fontsize=9)
    bars.invert_yaxis()
    bars.set_xlabel("Biological CTAs (identity unchanged)")
    bars.set_title("Host-screen outcome for each CTA", loc="left", weight="bold", pad=12)
    bars.spines[["top", "right"]].set_visible(False)
    bars.legend(loc="upper left", bbox_to_anchor=(0, -0.14), frameon=False, fontsize=9)
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
    fig, (ax_n, ax_p) = plt.subplots(
        1,
        2,
        figsize=(16.2, 0.44 * len(plot_df) + 2.4),
        sharey=True,
        gridspec_kw={"width_ratios": [1.28, 1.12], "wspace": 0.08},
    )
    left_n = np.zeros(len(plot_df))
    left_p = np.zeros(len(plot_df))
    totals = plot_df["mapped_proteins"].to_numpy(dtype=float)
    for col, label, color in categories:
        vals = plot_df[col].to_numpy(dtype=float)
        pcts = np.where(totals > 0, 100.0 * vals / totals, 0.0)
        ax_n.barh(y, vals, left=left_n, color=color, label=label, height=0.72)
        ax_p.barh(y, pcts, left=left_p, color=color, height=0.72)
        left_n += vals
        left_p += pcts

    max_n = totals.max()
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

    # Draw subtle horizontal dividers and class headers between antigen classes.
    prev_class = None
    for idx, cls in enumerate(plot_df["antigen_class"]):
        if prev_class is not None and cls != prev_class:
            sep_y = y[idx] + 0.5
            ax_n.axhline(sep_y, color="#495f6b", lw=1.0, ls="--", alpha=0.65)
            ax_p.axhline(sep_y, color="#495f6b", lw=1.0, ls="--", alpha=0.65)
        prev_class = cls

    # Annotate antigen-class group labels cleanly in the right margin of ax_p.
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

    comma = FuncFormatter(lambda val, _: f"{int(val):,}")
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


def build():
    data = pd.read_csv(OUT / "protein-decisions.csv.gz", low_memory=False)
    FIG.mkdir(exist_ok=True)
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

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
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import cta_selection_policy, target_filters, target_selection  # noqa: E402

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
    ("Home tissue no panel resolves", "Other:Tissue", EXPRESSION),
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
    """Full scale on the left; the right panel rescales after the first, dominant step."""
    table = table[(table.removed > 0) | (table.index == 0)].reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 0.46 * len(table) + 1.6), sharey=True)
    y = np.arange(len(table))[::-1]
    for ax, start in zip(axes, [0, zoom_from]):
        part = table.iloc[start:]
        ax.barh(y[start:], part.remaining, color=part.color)
        top = part.remaining.max()
        for yi, row in zip(y[start:], part.itertuples()):
            note = f"{row.remaining:,}" + (f"  (−{row.removed:,})" if row.removed else "")
            ax.text(row.remaining + top * 0.01, yi, note, va="center", fontsize=9)
        ax.set_xlim(0, top * 1.25)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlabel("Proteins remaining")
    axes[0].set_yticks(y, table.step)
    axes[0].set_title(title, loc="left", fontsize=11)
    axes[1].set_title("Same steps, rescaled after the first filter", loc="left", fontsize=10)
    fig.tight_layout()
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
    fig, ax = plt.subplots(figsize=(8.4, 7))
    for (x, y), color in zip(
        [(-0.75, 0.45), (0.75, 0.45), (0, -0.8)], [CANCER, REPRODUCTIVE, PUBLISHED]
    ):
        ax.add_patch(Circle((x, y), 1.35, color=color, alpha=0.33, lw=1.5, ec=color))
    for (x, y), text in {
        (-1.45, 0.85): n(a & ~b & ~c),
        (1.45, 0.85): n(~a & b & ~c),
        (0, -1.55): n(~a & ~b & c),
        (0, 1.0): n(a & b & ~c),
        (-0.95, -0.4): n(a & ~b & c),
        (0.95, -0.4): n(~a & b & c),
        (0, 0.05): n(a & b & c),
    }.items():
        ax.text(x, y, text, ha="center", va="center", fontsize=15, fontweight="bold")
    ax.text(-1.5, 2.0, f"Cancer Antigens\n{n(a)}", ha="center", color=CANCER, fontsize=12)
    ax.text(
        1.5, 2.0, f"Reproductive Antigens\n{n(b)}", ha="center", color=REPRODUCTIVE, fontsize=12
    )
    ax.text(
        0,
        -2.55,
        f"Published CTA / oncofetal nominations\n{n(c)}",
        ha="center",
        color="#8a6a25",
        fontsize=12,
    )
    ax.set_xlim(-2.6, 2.6)
    ax.set_ylim(-2.9, 2.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(
        f"Rule-derived sets against published nominations ({len(data):,} proteins; areas not to scale)",
        fontsize=11,
    )
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return {
        "cta_published": int((a & b & c).sum()),
        "cta_new": int((a & b & ~c).sum()),
        "published_outside": int((~a & ~b & c).sum()),
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
            "negative_ihc_grade",
            "unsafe_evidence",
            "value",
        ],
        low_memory=False,
    )
    peaks = peaks[
        peaks.assay.eq("Normal bulk")
        & peaks.selection_tier.isin(["vital", "protected", "manageable"])
    ]
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
        "Bulk-RNA peak observations counted as Unsafe evidence / recorded, by ceiling cell",
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
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2), sharex=True, sharey=True)
    for ax, column, membership, title, color in [
        (axes[0], "tumor_max_tpm", "is_cancer_antigen", "Cancer", CANCER),
        (
            axes[1],
            "reproductive_supported_tpm",
            "is_reproductive_antigen",
            "Reproductive / developmental",
            REPRODUCTIVE,
        ),
    ]:
        d = data[data[column].notna()]
        x = d[["normal_reference_tpm", "disposable_max_tpm"]].max(axis=1).clip(lower=floor)
        y = d[column].clip(lower=floor)
        for label, mask, shade in [
            ("Qualified membership", d[membership], color),
            ("No membership / held", ~d[membership], EXPRESSION),
        ]:
            ax.scatter(
                x[mask], y[mask], s=9, alpha=0.6, color=shade, label=f"{label} ({int(mask.sum())})"
            )
        grid = np.logspace(np.log10(floor), np.log10(max(100, x.max())), 100)
        ax.plot(
            grid,
            target_selection.ENRICHMENT_FOLD * np.maximum(grid, 1),
            "k:",
            lw=1.3,
            label="10 × max(N, Disposable, 1)",
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Highest nonreproductive reference, including Disposable (TPM)")
        ax.set_ylabel("Highest supported expression (TPM)")
        ax.set_title(f"{title}: {len(d)} measured proteins", loc="left", fontsize=11)
        ax.legend(fontsize=8, frameon=False, loc="lower right")
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
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


def source_to_tiers_architecture(data: pd.DataFrame, funnel_all: pd.DataFrame, stem_path: Path):
    """End-to-end schematic: upstream data sources -> induced sets -> gradual funnel -> CTA tiers."""
    n_proteins = len(data)
    n_genes = sum(len(json.loads(m)) for m in data.member_gene_ids)
    n_nominated = int(data.published_nomination.fillna(False).astype(bool).sum())
    n_unsafe = int(data.category.eq("Unsafe").sum())
    after_unsafe = n_proteins - n_unsafe
    n_review = int(data.classification.eq("Other:Review").sum())
    n_sequence = int(data.classification.eq("Other:Sequence").sum())
    n_shared_elig = int(data.eligible.astype(bool).sum())
    n_bounds_fail = after_unsafe - n_review - n_sequence - n_shared_elig
    n_tissue = int(data.classification.eq("Other:Tissue").sum())
    n_expr_total = int(data.classification.eq("Other:Expression").sum())
    n_neither_10x = n_expr_total - n_bounds_fail
    n_ca = int(data.is_cancer_antigen.astype(bool).sum())
    n_ra = int(data.is_reproductive_antigen.astype(bool).sum())
    n_cta = int(data.is_cta.astype(bool).sum())
    n_ca_only = int(data.category.eq("Cancer Antigen only").sum())
    n_ra_only = int(data.category.eq("Reproductive Antigen only").sum())
    n_cand = int(data.is_target_candidate.astype(bool).sum())

    coverage = pd.read_csv(OUT / "tumor-source-coverage.csv")
    assessed_cohorts = coverage[coverage.status.eq("assessed")]
    n_cohorts = len(assessed_cohorts)
    n_tumor_samples = int(assessed_cohorts.samples.sum())
    repro_obs = pd.read_csv(
        OUT / "reproductive-comparison.csv.gz",
        usecols=["scope", "sample_id", "compartment"],
        low_memory=False,
    )
    n_repro_biopsies = int(repro_obs.loc[repro_obs.scope.eq("reproductive"), "sample_id"].nunique())

    tiers = data.loc[data.is_cta.astype(bool), "cta_tier"].value_counts()
    n_core_tier = int(tiers.get("core", 0))
    n_core_ext = int(tiers.get("core_with_extended_exposure", 0))
    n_ext_only = int(tiers.get("extended_only", 0))
    profiles = target_filters.profile_summary(data).set_index("profile")

    source_specs = [
        ("HPA v23 (186 biopsies + bulk) & GTEx v11 (54 tissues, brain)", "TPM / nTPM"),
        ("HPA v23 Single-Cell Tissue Atlas (81 cell clusters)", "nTPM"),
        ("Sorted Immune Panels (DICE 15 + Monaco 29 + HPA subsets)", "TPM / nTPM"),
        ("Donor scRNA/snRNA UMI (CELLxGENE 93 studies + Siletti 13.6M)", ">=2 UMI/donor"),
        ("Normal-Tissue Protein Staining (HPA v23 IHC, 45 organs)", "graded IHC"),
        ("Reviewed Normal-Expression Literature (912-protein rubric)", "quoted assays"),
        ("Adult Human Cornea scRNA (Collin et al.)", "CP10K"),
    ]

    fig, ax = plt.subplots(figsize=(16.5, 15.2))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    def box(x, y, w, h, fc="#ffffff", ec="#b8c4cc", lw=1.2, pad=0.35):
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle=f"round,pad={pad},rounding_size=0.7",
                facecolor=fc,
                edgecolor=ec,
                lw=lw,
            )
        )

    def arrow(x0, y0, x1, y1, color="#5c6c75", lw=1.4):
        ax.add_patch(
            FancyArrowPatch(
                (x0, y0),
                (x1, y1),
                arrowstyle="-|>,head_length=5,head_width=3",
                color=color,
                lw=lw,
                mutation_scale=1.0,
            )
        )

    # Header
    ax.text(
        50,
        98.4,
        "OncoRef Antigen Curation Architecture: Data Sources → Induced Protein Sets → CTA Tiers",
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
        color="#1d2d35",
    )
    ax.text(
        50,
        96.5,
        f"Policy {target_selection.POLICY_ID} · {n_proteins:,} human protein groups from {n_genes:,} Ensembl 115 coding loci · Every count computed from the packaged catalog",
        ha="center",
        va="center",
        fontsize=9.6,
        color="#4e6069",
    )

    # Stage 0: Reference Universe & Literature Nominations (y = 89.2 .. 94.8)
    box(2.0, 89.4, 30.0, 5.1, fc="#f2f6f8", ec="#8fa2ad")
    ax.text(
        3.2,
        93.3,
        "Genome & Transcript Identity Sources",
        fontsize=9.5,
        fontweight="bold",
        color="#243742",
    )
    ax.text(
        3.2,
        90.4,
        f"Ensembl 115 ({n_genes:,} coding genes)\n+ MANE Select / Canonical (transcript-selection)",
        fontsize=8.4,
        color="#3d505b",
        linespacing=1.35,
    )
    arrow(32.5, 92.0, 35.8, 92.0)

    box(36.3, 89.4, 27.4, 5.1, fc="#e6eef2", ec="#5c7684", lw=1.5)
    ax.text(
        50.0,
        92.9,
        f"{n_proteins:,} Human Protein Groups",
        ha="center",
        fontsize=11.5,
        fontweight="bold",
        color="#1d2d35",
    )
    ax.text(
        50.0,
        90.5,
        "Identical reference sequences grouped;\nall member Ensembl gene IDs retained",
        ha="center",
        fontsize=8.2,
        color="#3d505b",
        linespacing=1.3,
    )

    box(68.0, 89.4, 30.0, 5.1, fc="#fcf7eb", ec="#c99a46")
    ax.text(
        69.2,
        93.3,
        f"Literature Nominations ({n_nominated:,} proteins)",
        fontsize=9.5,
        fontweight="bold",
        color="#7a5b1b",
    )
    ax.text(
        69.2,
        90.4,
        "16 curated CTA & oncofetal catalogs (CTpedia,\nBruggeman, Da Silva, ...); annotation only",
        fontsize=8.3,
        color="#5c4922",
        linespacing=1.35,
    )
    arrow(67.5, 92.0, 64.2, 92.0, color="#b38732")
    arrow(50.0, 89.0, 50.0, 86.7)

    # Stage 1: Normal Somatic Exclusion Screen (y = 57.2 .. 86.2)
    box(2.0, 57.2, 96.0, 28.8, fc="#fdf8f7", ec="#cfa3a0", lw=1.3)
    ax.text(
        3.2,
        84.6,
        f"1. Shared Somatic Safety Screen — 7 Normal Data Panels  ({n_proteins:,} → {after_unsafe:,} survive;  −{n_unsafe:,} Unsafe)",
        fontsize=10.5,
        fontweight="bold",
        color="#7d2b29",
    )
    ax.text(
        3.4,
        82.6,
        "Upstream Normal Data Source & Scale",
        fontsize=8.5,
        fontweight="bold",
        color="#526068",
    )
    ax.text(
        47.5,
        82.6,
        "Sequential Funnel (surviving proteins)",
        fontsize=8.5,
        fontweight="bold",
        color="#526068",
    )
    ax.text(79.0, 82.6, "Total Flagged by Layer", fontsize=8.5, fontweight="bold", color="#526068")

    unsafe_steps = funnel_all.iloc[1 : 1 + len(IDEAS)].reset_index(drop=True)
    y_rows = np.linspace(80.1, 61.2, len(unsafe_steps))
    for idx, (y_i, row, (src_label, scale_label)) in enumerate(
        zip(y_rows, unsafe_steps.itertuples(), source_specs)
    ):
        ax.text(3.4, y_i, f"{idx + 1}. {src_label}", va="center", fontsize=8.2, color="#243742")
        ax.text(
            45.8, y_i, f"[{scale_label}]", va="center", ha="right", fontsize=7.6, color="#687780"
        )
        # Draw rescaled bar for surviving proteins (max 2,500 after step 1; step 1 starts from 19,893)
        bar_w = 18.0 * (row.remaining / 2600.0)
        ax.add_patch(
            FancyBboxPatch(
                (47.5, y_i - 0.85),
                max(bar_w, 0.4),
                1.7,
                boxstyle="round,pad=0.02,rounding_size=0.25",
                facecolor="#5f7d8c",
                edgecolor="none",
            )
        )
        ax.text(
            47.5 + max(bar_w, 0.4) + 0.6,
            y_i,
            f"{int(row.remaining):,}  (−{int(row.removed):,})",
            va="center",
            fontsize=8.2,
            fontweight="bold" if idx == len(unsafe_steps) - 1 else "normal",
            color="#243742",
        )
        flagged_total = int(row.flagged_in_total)
        ax.text(
            79.0,
            y_i,
            f"{flagged_total:,} Unsafe",
            va="center",
            fontsize=8.1,
            color="#8f3532",
        )

    box(78.5, 58.0, 18.3, 2.3, fc="#ac5552", ec="#8c3b38", pad=0.2)
    ax.text(
        87.65,
        59.15,
        f"Unsafe Set: {n_unsafe:,} proteins",
        ha="center",
        va="center",
        fontsize=9.2,
        fontweight="bold",
        color="#ffffff",
    )
    ax.text(
        3.4,
        58.6,
        "Note: Many genes fail multiple normal panels; sequential (−removed) attributes each protein to the first panel it fails.",
        fontsize=7.9,
        color="#5e6d75",
    )
    arrow(50.0, 56.8, 50.0, 54.5)

    # Stage 2: Technical, Sequence & Bounds Eligibility (y = 43.2 .. 54.0)
    box(2.0, 43.2, 30.0, 10.6, fc="#f4f6f8", ec="#9aa8b2")
    ax.text(
        3.2,
        52.3,
        "2. Eligibility & Sequence Sources",
        fontsize=9.5,
        fontweight="bold",
        color="#243742",
    )
    ax.text(
        3.2,
        44.4,
        "• Complete 36-tissue native-TPM normal panel\n"
        "• MANE/canonical ORF & translation qualification\n"
        "• Ensembl 115 9-mers (<95% in Unsafe)\n"
        "• Curated readthrough-loci.csv\n"
        "• Disposable tissues <= 32 TPM & <= 32 nTPM",
        fontsize=8.0,
        color="#3d505b",
        linespacing=1.35,
    )
    arrow(32.5, 48.5, 35.8, 48.5)

    box(36.3, 43.2, 27.4, 10.6, fc="#eaf2f0", ec="#4c8473", lw=1.4)
    ax.text(
        50.0,
        52.2,
        "Shared-Eligible Gate",
        ha="center",
        fontsize=10.2,
        fontweight="bold",
        color="#1f4e40",
    )
    ax.text(
        50.0,
        47.9,
        f"{after_unsafe:,} surviving somatic screen\n"
        f"  − {n_review:,} Other:Review\n"
        f"  − {n_sequence:,} Other:Sequence\n"
        f"  − {n_bounds_fail:,} Bounds fail",
        ha="center",
        va="center",
        fontsize=8.3,
        color="#2c443d",
        linespacing=1.32,
    )
    ax.text(
        50.0,
        44.2,
        f"→ {n_shared_elig:,} Shared-Eligible (H₀ = max(N, D, 1))",
        ha="center",
        fontsize=9.1,
        fontweight="bold",
        color="#1f4e40",
    )

    arrow(64.2, 50.8, 67.5, 50.8, color="#9e7430")
    box(68.0, 48.8, 30.0, 5.0, fc="#fbf5ea", ec="#b88a3c")
    ax.text(
        69.2,
        52.3,
        f"Other:Review ({n_review:,} proteins)",
        fontsize=9.2,
        fontweight="bold",
        color="#855f1f",
    )
    ax.text(
        69.2,
        49.5,
        "Incomplete 36-tissue coverage, uncertain ORF,\nor unresolved donor rare-cell / mapping signal",
        fontsize=7.9,
        color="#5c4926",
        linespacing=1.28,
    )

    arrow(64.2, 45.6, 67.5, 45.6, color="#73608f")
    box(68.0, 43.2, 30.0, 4.9, fc="#f4f1f8", ec="#8d7aa8")
    ax.text(
        69.2,
        46.6,
        f"Other:Sequence ({n_sequence:,} proteins)",
        fontsize=9.2,
        fontweight="bold",
        color="#574473",
    )
    ax.text(
        69.2,
        43.9,
        ">=95% 9-mer overlap with Unsafe proteins,\nnested microprotein, or read-through locus",
        fontsize=7.9,
        color="#473b59",
        linespacing=1.28,
    )

    arrow(50.0, 42.8, 50.0, 40.4)

    # Stage 3: Symmetric 10x Compartment Enrichment (y = 24.0 .. 40.0)
    box(2.0, 24.0, 29.0, 15.8, fc="#edf5fa", ec="#3e8bba", lw=1.3)
    ax.text(
        3.2,
        38.2,
        "3A. Tumor RNA Sources (CA)",
        fontsize=9.6,
        fontweight="bold",
        color="#1d5b82",
    )
    ax.text(
        3.2,
        31.8,
        f"• {n_cohorts} linear-TPM cancer cohorts\n"
        f"  ({n_tumor_samples:,} QC-passing tumor libraries;\n"
        "  library_qc screens spermatid/OR/ectopic)\n"
        "• nose-skin-reference.csv & family homes\n"
        "• Rule: >=1 tumor >= 10 x H₀",
        fontsize=8.0,
        color="#274659",
        linespacing=1.33,
    )
    ax.text(
        3.2,
        25.3,
        f"Induces {n_ca:,} Cancer Antigens (CA)\n({n_tissue:,} unmeasured home -> Other:Tissue)",
        fontsize=8.5,
        fontweight="bold",
        color="#1d5b82",
        linespacing=1.3,
    )
    arrow(31.5, 32.0, 34.3, 32.0, color="#3e8bba")

    box(69.0, 24.0, 29.0, 15.8, fc="#edf7f1", ec="#438a69", lw=1.3)
    ax.text(
        70.2,
        38.2,
        "3B. Reproductive & Embryo Sources (RA)",
        fontsize=9.6,
        fontweight="bold",
        color="#225e42",
    )
    ax.text(
        70.2,
        31.8,
        f"• HPA v23 ({n_repro_biopsies} biopsies across 11 adult\n"
        "  compartments: lactating breast.410 split\n"
        "  from 3 resting breast biopsies)\n"
        "• Petropoulos E3–E7 & Tyser CS7 embryos\n"
        "• Rule: sample or lineage >= 10 x H₀",
        fontsize=8.0,
        color="#264a39",
        linespacing=1.33,
    )
    ax.text(
        70.2,
        25.3,
        f"Induces {n_ra:,} Reproductive Antigens (RA)\n(embryo replication enforced per stage)",
        fontsize=8.5,
        fontweight="bold",
        color="#225e42",
        linespacing=1.3,
    )
    arrow(68.5, 32.0, 65.7, 32.0, color="#438a69")

    box(34.8, 24.0, 30.4, 15.8, fc="#f7fbf9", ec="#3f8f6b", lw=1.6)
    ax.text(
        50.0,
        38.1,
        f"{n_cand:,} Candidate Antigens (CA ∪ RA)",
        ha="center",
        fontsize=10.6,
        fontweight="bold",
        color="#1d3537",
    )
    # Three partition cards
    box(36.0, 29.3, 8.2, 6.8, fc="#dcebf5", ec="#3e8bba", pad=0.2)
    ax.text(40.1, 33.8, "CA only", ha="center", fontsize=8.6, fontweight="bold", color="#1d5b82")
    ax.text(
        40.1, 31.1, f"{n_ca_only:,}", ha="center", fontsize=12.0, fontweight="bold", color="#1d5b82"
    )

    box(45.2, 29.0, 9.6, 7.4, fc="#2f7f6f", ec="#1f5c50", lw=1.5, pad=0.2)
    ax.text(
        50.0, 34.4, "CTA (CA ∩ RA)", ha="center", fontsize=8.6, fontweight="bold", color="#ffffff"
    )
    ax.text(
        50.0, 31.7, f"{n_cta:,}", ha="center", fontsize=13.5, fontweight="bold", color="#ffffff"
    )
    ax.text(50.0, 29.8, 'default "all"', ha="center", fontsize=7.5, color="#e2f5f0")

    box(55.8, 29.3, 8.2, 6.8, fc="#dff0e6", ec="#438a69", pad=0.2)
    ax.text(59.9, 33.8, "RA only", ha="center", fontsize=8.6, fontweight="bold", color="#225e42")
    ax.text(
        59.9, 31.1, f"{n_ra_only:,}", ha="center", fontsize=12.0, fontweight="bold", color="#225e42"
    )

    ax.text(
        50.0,
        25.4,
        f"Excluded here: {n_neither_10x:,} below 10x (+{n_bounds_fail:,} bounds = {n_expr_total:,} Other:Expression)\nand {n_tissue:,} Other:Tissue",
        ha="center",
        fontsize=7.8,
        color="#55656e",
        linespacing=1.28,
    )

    arrow(50.0, 23.6, 50.0, 21.3, color="#2f7f6f", lw=1.6)

    # Stage 4: CTA Biological Tiers & Host Toxicity Profiles (y = 1.4 .. 21.0)
    box(2.0, 1.4, 96.0, 19.6, fc="#f6faf8", ec="#2f7f6f", lw=1.4)
    ax.text(
        3.2,
        19.4,
        f"4. Slicing the {n_cta:,} Default CTAs — Biological Compartment Tiers & Patient Toxicity Profiles",
        fontsize=10.5,
        fontweight="bold",
        color="#1f5c50",
    )

    # Left subpanel: 3 biological CTA tiers
    ax.text(
        3.4,
        17.3,
        f"Biological Origin & Exposure Tiers (sum = {n_cta} CTAs)",
        fontsize=9.0,
        fontweight="bold",
        color="#243742",
    )
    tier_cards = [
        (
            12.4,
            f"Core CTA: {n_core_tier}",
            "Qualifies in embryo, placenta, testis or ovary; all 8 extended adult organs <=32 TPM & >=10x below tumor (MAGEA4, NY-ESO-1, PRAME, SSX2/B, MAGEC2)",
            "#2f7f6f",
        ),
        (
            7.2,
            f"Core + Extended Exposure: {n_core_ext}",
            "Core origin, plus >32 TPM or <10x tumor margin in epididymis, tract or resting breast (MAGEA11, MAGEB4, RNASE11, KIAA1210, ADAM2, TCP11X2)",
            "#5a8f7d",
        ),
        (
            2.0,
            f"Extended-Only CTA: {n_ext_only}",
            "Qualifies only via accessory/mammary/tract: 21 male accessory (TGM4, RLN1, SEMG2, SPAG11A/B), 2 lactating breast (LALBA, BTN1A1), 1 tract (SCGB1D4)",
            "#9c7228",
        ),
    ]
    for y_card, title_txt, body_txt, badge_col in tier_cards:
        box(3.4, y_card, 44.5, 4.3, fc="#ffffff", ec=badge_col, lw=1.1, pad=0.2)
        ax.text(
            4.3,
            y_card + 3.45,
            title_txt,
            va="center",
            fontsize=8.5,
            fontweight="bold",
            color=badge_col,
        )
        ax.text(
            4.3,
            y_card + 2.75,
            textwrap.fill(body_txt, 69),
            va="top",
            fontsize=7.35,
            color="#374952",
            linespacing=1.18,
        )

    # Right subpanel: Host / Patient Toxicity Profiles
    ax.text(
        50.2,
        17.3,
        "Configurable Host Profiles (select protected vs tolerated compartments)",
        fontsize=9.0,
        fontweight="bold",
        color="#243742",
    )
    prof_rows = [
        ("all", "all (default: no compartments protected)", "#2f7f6f"),
        ("female-patient", "female-patient (protects tract + resting breast)", "#3d8c78"),
        (
            "female-patient (lactating)",
            "female-patient --lactating (drops LALBA, BTN1A1)",
            "#4c9985",
        ),
        (
            "male-patient (post-prostatectomy)",
            "male-patient --allow-tissue prostate,sem. vesicle",
            "#4a7da8",
        ),
        ("male-patient", "male-patient (protects male acc. + resting breast)", "#3e8bba"),
        ("core", "core (protects all 8 extended adult compartments)", "#5c6e7a"),
    ]
    for idx, (pkey, plabel, pcol) in enumerate(prof_rows):
        y_p = 15.0 - idx * 2.38
        cta_val = int(profiles.loc[pkey, "CTA"])
        bw = 11.8 * (cta_val / max(n_cta, 1))
        ax.text(50.2, y_p, plabel, va="center", fontsize=7.8, color="#243742")
        ax.add_patch(
            FancyBboxPatch(
                (79.4, y_p - 0.7),
                bw,
                1.4,
                boxstyle="round,pad=0.02,rounding_size=0.2",
                facecolor=pcol,
                edgecolor="none",
            )
        )
        ax.text(
            79.4 + bw + 0.5,
            y_p,
            f"{cta_val} CTA",
            va="center",
            fontsize=8.1,
            fontweight="bold",
            color="#1d2d35",
        )

    _save_both(fig, stem_path, dpi=170)


def cta_tiers_and_profiles(data: pd.DataFrame, stem_path: Path):
    """Two-panel visual of the 12 compartments x host profiles and resulting CTA sets."""
    scenarios = [
        ("all (default)", target_filters.TargetFilter(host_profile="all")),
        ("female-patient", target_filters.TargetFilter(host_profile="female-patient")),
        (
            "female-patient\n(--lactating)",
            target_filters.TargetFilter(host_profile="female-patient", lactating=True),
        ),
        (
            "male-patient\n(post-prostatectomy)",
            target_filters.TargetFilter(
                host_profile="male-patient", allow_tissues=("prostate", "seminal vesicle")
            ),
        ),
        ("male-patient", target_filters.TargetFilter(host_profile="male-patient")),
        ("core", target_filters.TargetFilter(host_profile="core")),
        (
            "female-patient\n(--pregnant)",
            target_filters.TargetFilter(host_profile="female-patient", pregnant=True),
        ),
    ]
    compartments = list(target_selection.REPRODUCTIVE_COMPARTMENTS)
    comp_labels = [
        "embryo (Core)",
        "placenta (Core)",
        "testis (Core)",
        "ovary (Core)",
        "epididymis (Male)",
        "prostate (Male)",
        "seminal vesicle (Male)",
        "fallopian tube (Female)",
        "endometrium (Female)",
        "cervix (Female)",
        "lactating_breast [breast.410]",
        "resting_breast [3 biopsies]",
    ]
    matrix = np.zeros((len(scenarios), len(compartments)), dtype=int)
    tier_counts = []
    measured = target_filters.measures(data)
    for i, (_, filt) in enumerate(scenarios):
        prot = set(target_filters.protected_compartments(filt))
        for j, comp in enumerate(compartments):
            matrix[i, j] = 1 if comp in prot else 0
        ev = target_filters.apply(
            target_filters.TargetFilter(
                membership="cta",
                host_profile=filt.host_profile,
                lactating=filt.lactating,
                pregnant=filt.pregnant,
                allow_tissues=filt.allow_tissues,
                protect_tissues=filt.protect_tissues,
            ),
            measured,
        )
        kept = ev[ev.keep]
        vc = kept.cta_tier.value_counts()
        tier_counts.append(
            (
                int(vc.get("core", 0)),
                int(vc.get("core_with_extended_exposure", 0)),
                int(vc.get("extended_only", 0)),
            )
        )

    fig, (ax_grid, ax_bars) = plt.subplots(
        1, 2, figsize=(15.5, 6.4), gridspec_kw={"width_ratios": [1.35, 1.0]}
    )

    # Left: Compartment allowance vs protection grid (compartments on y, profiles on x)
    grid_t = matrix.T  # shape: (12 compartments, 7 scenarios)
    cmap = matplotlib.colors.ListedColormap(["#d9efe6", "#f4d4d2"])
    ax_grid.imshow(grid_t, cmap=cmap, aspect="auto", vmin=0, vmax=1)
    ax_grid.set_yticks(range(len(comp_labels)), comp_labels, fontsize=9)
    ax_grid.set_xticks(
        range(len(scenarios)),
        [name for name, _ in scenarios],
        fontsize=8.3,
        rotation=25,
        ha="right",
    )
    for r_i in range(grid_t.shape[0]):
        for c_i in range(grid_t.shape[1]):
            protected = bool(grid_t[r_i, c_i])
            ax_grid.text(
                c_i,
                r_i,
                "Protected" if protected else "Allowed",
                ha="center",
                va="center",
                fontsize=7.8,
                fontweight="bold" if protected else "normal",
                color="#8a2c28" if protected else "#1f5c48",
            )
    for sep in (3.5, 6.5, 9.5):
        ax_grid.axhline(sep, color="#5c6c75", lw=1.1, ls="--")
    ax_grid.set_title(
        "Which compartments are protected from toxicity vs allowed as CTA origins?",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
        pad=10,
    )
    ax_grid.tick_params(length=0)

    # Right: Stacked horizontal bars of CTA counts per scenario
    y_pos = np.arange(len(scenarios))[::-1]
    core_vals = np.array([t[0] for t in tier_counts])
    cext_vals = np.array([t[1] for t in tier_counts])
    ext_vals = np.array([t[2] for t in tier_counts])
    totals = core_vals + cext_vals + ext_vals

    ax_bars.barh(
        y_pos, core_vals, color="#2f7f6f", label="Core CTA (germline/fetoplacental, clean extended)"
    )
    ax_bars.barh(
        y_pos,
        cext_vals,
        left=core_vals,
        color="#6ea899",
        label="Core + extended organ exposure",
    )
    ax_bars.barh(
        y_pos,
        ext_vals,
        left=core_vals + cext_vals,
        color="#c99a46",
        label="Extended-only CTA (prostate, epididymis, lactating breast, tract)",
    )
    for yi, total, c_v in zip(y_pos, totals, core_vals):
        detail = f"{int(total)} CTA" + (
            f"  ({int(c_v)} core + {int(total - c_v)} ext)" if total != c_v else "  (core)"
        )
        ax_bars.text(
            total + 4, yi, detail, va="center", fontsize=8.6, fontweight="bold", color="#1d2d35"
        )
    ax_bars.set_yticks(y_pos, [name.replace("\n", " ") for name, _ in scenarios], fontsize=8.8)
    ax_bars.set_xlim(0, max(totals) * 1.48)
    ax_bars.set_xlabel("Retained CTA proteins")
    ax_bars.set_title(
        "Induced CTA set size by host / patient profile",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
        pad=10,
    )
    ax_bars.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, -0.14),
        fontsize=8.0,
        framealpha=0.95,
    )
    ax_bars.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    _save_both(fig, stem_path, dpi=165)


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
    fates = published.classification.value_counts()
    steps = {"step": "Filtering idea", "removed": "Removed here", "remaining": "Remaining"}
    (OUT / "summary.md").write_text(
        "# Curation summary\n\n"
        f"Policy `{target_selection.POLICY_ID}`, {len(data):,} human protein records. "
        f"**{counts['CTA']} CTA, {counts['Cancer Antigen only']} Cancer Antigen only, "
        f"{counts['Reproductive Antigen only']} Reproductive Antigen only.** These are outputs of "
        "screening rules, not validated tumor antigens. Every figure is drawn from the "
        "[rebuilt catalog](index.md); thresholds are read from the policy constants.\n\n"
        "## End-to-end architecture: data sources to CTA tiers\n\n"
        "![End-to-end curation architecture](summary/source-to-tiers-architecture.png)\n\n"
        "Every upstream panel induces a defined exclusion, review gate, or enrichment set. "
        "Lactating mammary (`breast.410`, expressing coordinated milk proteins `LALBA`, `CSN1S1`, "
        "`CSN2`, `CSN3`) is separated from resting mammary (`breast.373`, `breast.390`, `breast.405`), "
        "keeping lactation-specific antigens available for non-lactating women and men while allowing "
        "explicit exclusion via `--lactating`.\n\n"
        "## CTA compartment tiers and host toxicity profiles\n\n"
        "![CTA tiers and host profiles](summary/cta-tiers-and-profiles.png)\n\n"
        + markdown(
            profile_table,
            {
                "profile": "Host profile",
                "protected_compartments": "Protected compartments (<=32 TPM & >=10x below tumor)",
                "CTA": "CTA",
                "CA_only": "CA only",
                "RA_only": "RA only",
                "Union": "CA | RA",
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
        "## Shared enrichment rule\n\nBoth CA and RA independently require 10× max(N, Disposable, 1). The cancer panel below shows the tumor side.\n\n![Tumor against normal](summary/tumor-normal.png)\n\n"
        "## Rare-cell RNA\n\n![Rare-cell grading](summary/rare-cell.png)\n\n"
        "[Donor-level evidence](../../normal-cell-atlas-20261110/rare-cell/index.md)\n\n"
        "## Inside the later stages\n\n![Stage detail](summary/stages.png)\n\n"
        + "\n".join(f"- {label}: {int(value)}" for label, value in why.items())
        + "\n\nA protein can use several routes, so route counts exceed set sizes.\n"
    )


if __name__ == "__main__":
    build()

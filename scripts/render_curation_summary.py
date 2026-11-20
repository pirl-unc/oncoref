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
from matplotlib.patches import Circle, FancyBboxPatch

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
    from matplotlib.patches import Patch

    profiles = target_filters.profile_summary(data)
    states = ["Not protected", "Protect", "Accept loss", "Absent"]
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
    fig, (grid, bars) = plt.subplots(1, 2, figsize=(17, 8), gridspec_kw={"width_ratios": [1.25, 1]})
    grid.imshow(
        matrix.T, cmap=matplotlib.colors.ListedColormap(colors), vmin=0, vmax=3, aspect="auto"
    )
    grid.set_yticks(
        range(len(compartments)), [c.replace("_", " ") for c in compartments], fontsize=9
    )
    grid.set_xticks(range(len(labels)), labels, rotation=32, ha="right", fontsize=8)
    for i in range(len(labels)):
        for j in range(len(compartments)):
            grid.text(
                i,
                j,
                ["—", "Protect", "Loss", "Absent"][matrix[i, j]],
                ha="center",
                va="center",
                fontsize=8,
            )
    grid.set_title("Declared anatomy and tissue-loss choices", loc="left", weight="bold", pad=14)
    grid.tick_params(length=0)
    grid.legend(
        handles=[Patch(color=c, label=s) for c, s in zip(colors, states)],
        loc="upper left",
        bbox_to_anchor=(0, -0.29),
        ncol=2,
        fontsize=9,
    )
    y = np.arange(len(profiles))
    left = np.zeros(len(profiles))
    for column, label, color in [
        ("CTA", "Retained for research screen", "#438a69"),
        ("CTA_unresolved", "Unresolved coverage / evidence", "#b88a3c"),
        ("CTA_excluded", "Fails a host requirement", "#ac5552"),
    ]:
        values = profiles[column].to_numpy()
        bars.barh(y, values, left=left, label=label, color=color)
        for yi, offset, v in zip(y, left, values):
            if v:
                bars.text(
                    offset + v / 2, yi, str(v), ha="center", va="center", fontsize=9, color="white"
                )
        left += values
    bars.set_yticks(y, [s.replace("\n", " ") for s in labels], fontsize=9)
    bars.invert_yaxis()
    bars.set_xlabel("Biological CTAs (identity unchanged)")
    bars.set_title("Host-screen outcome for each CTA", loc="left", weight="bold", pad=14)
    bars.spines[["top", "right"]].set_visible(False)
    bars.legend(loc="upper left", bbox_to_anchor=(0, -0.12), fontsize=9)
    fig.suptitle(
        "Protein origin and tolerated tissue loss are separate questions",
        fontsize=16,
        weight="bold",
    )
    fig.text(
        0.02,
        0.015,
        "No reproductive protection by default. Select tissues to protect explicitly. Presets record assumed anatomy; pregnancy coverage is incomplete.",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
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

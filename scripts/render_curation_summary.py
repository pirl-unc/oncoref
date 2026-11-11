#!/usr/bin/env python3
"""Draw the curation summary page: set overlap, funnels, thresholds and distributions.

Every number is read from the rebuilt catalog and the policy constants; nothing
here classifies a protein.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import cta_selection_policy, target_selection  # noqa: E402

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
    ("Incomplete coverage, identity or translation", "Other:Review", REVIEW),
    ("Nested or Unsafe-contained sequence", "Other:Sequence", SEQUENCE),
    ("No qualifying expression or bounds fail", "Other:Expression", EXPRESSION),
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
        raise ValueError("Funnel does not end at the Target candidates")
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
        ("reproductive_ntpm", "HPA nTPM: highest reproductive tissue", [(1.0, "floor")], 0.01),
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
    d = data[pd.to_numeric(data.tumor_max_tpm, errors="coerce").notna()]
    floor = 0.05
    x = d.normal_reference_tpm.clip(lower=floor)
    y = d.tumor_max_tpm.clip(lower=floor)
    fig, ax = plt.subplots(figsize=(7.6, 6.6))
    for label, mask, color in [
        ("Cancer Antigen", d.is_cancer_antigen.astype(bool), CANCER),
        (
            "Reproductive Antigen only",
            d.classification.eq("Reproductive Antigen only"),
            REPRODUCTIVE,
        ),
        ("Other", d.category.eq("Other"), EXPRESSION),
    ]:
        ax.scatter(
            x[mask], y[mask], s=9, alpha=0.6, color=color, label=f"{label} ({int(mask.sum())})"
        )
    grid = np.logspace(np.log10(floor), 2, 100)
    fold = target_selection.SINGLE_TUMOR_FOLD
    ax.plot(grid, np.maximum(grid, 1), "k--", lw=1, label="T = max(N, 1): repeated-patient rule")
    ax.plot(
        grid,
        fold * np.maximum(grid, 1),
        "k:",
        lw=1.3,
        label=f"T = {fold:g} × max(N, 1): single-tumor rule",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Normal reference N: highest compatible normal context (TPM)")
    ax.set_ylabel("Highest qualified tumor sample T (TPM)")
    ax.set_title(
        f"Cancer rule on the {len(d)} proteins assessed for cancer evidence",
        loc="left",
        fontsize=11,
    )
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def rare(data, path):
    d = data[data.rare_cell_status.isin(["established", "dismissed", "below_screen"])]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4), gridspec_kw={"width_ratios": [3, 2]})
    for status, color in [
        ("below_screen", "#c4ccd1"),
        ("dismissed", REVIEW),
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
    axes[0].set_title("Rare-cell panel: single molecules do not decide", loc="left", fontsize=11)
    axes[0].legend(fontsize=8, frameon=False)
    basis = (
        d[d.rare_cell_status.ne("below_screen")]
        .groupby(["rare_cell_status", "rare_cell_basis"])
        .size()
    )
    basis = basis.sort_values()
    axes[1].barh(
        [b for _, b in basis.index],
        basis.values,
        color=[UNSAFE if s == "established" else REVIEW for s, _ in basis.index],
    )
    for i, value in enumerate(basis.values):
        axes[1].text(value + 2, i, str(value), va="center", fontsize=9)
    axes[1].set_title("Why a signal was established (red) or dismissed", loc="left", fontsize=11)
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
            "Patient identity unresolved in tumor samples": rest.sum(),
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


def build():
    data = pd.read_csv(OUT / "protein-decisions.csv.gz", low_memory=False)
    FIG.mkdir(exist_ok=True)
    counts = data.classification.value_counts()
    sets = venn(data, FIG / "venn.png")
    everything = draw_funnel(
        funnel(data, "All assessed proteins"),
        "Every assessed protein, by filtering idea",
        FIG / "funnel-all.png",
    )
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
        "## Cancer rule\n\n![Tumor against normal](summary/tumor-normal.png)\n\n"
        "## Rare-cell RNA\n\n![Rare-cell grading](summary/rare-cell.png)\n\n"
        "[Donor-level evidence](../../normal-cell-atlas-20261110/rare-cell/index.md)\n\n"
        "## Inside the later stages\n\n![Stage detail](summary/stages.png)\n\n"
        + "\n".join(f"- {label}: {int(value)}" for label, value in why.items())
        + "\n\nA protein can use several routes, so route counts exceed set sizes.\n"
    )


if __name__ == "__main__":
    build()

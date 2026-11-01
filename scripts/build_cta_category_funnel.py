"""Render a coarse view of the released policy, with normal assessment first.

Run from the OncoRef checkout using its local Python environment. This is a
presentation of released decisions; it does not alter classification rules.
"""

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

from oncoref.cta_selection_policy import CANDIDATE_ONLY_LAYERS
from oncoref.version import __version__


def build(root: Path, output: Path) -> None:
    source = root / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz"
    facts = pd.read_csv(source, low_memory=False)
    unsafe = facts[[c for c in facts if c.startswith("unsafe:")]].any(axis=1)
    normal_holds = facts[
        [
            c
            for c in facts
            if c.startswith("hold:") and c.removeprefix("hold:") not in CANDIDATE_ONLY_LAYERS
        ]
    ].any(axis=1)
    normal_pass = ~unsafe & ~normal_holds & facts.required_coverage & ~facts.unresolved_evidence
    unassessed = ~unsafe & ~normal_pass
    eligible = normal_pass.copy()
    other = pd.Series(False, index=facts.index)
    rows = []
    for label, failure in [
        ("Normal safety + coverage", pd.Series(False, index=facts.index)),
        ("Reproductive presence", ~facts.positive_context),
        ("Cancer presence", ~facts.cancer_association),
        ("Reproductive concentration", facts.reproductive_concentration_hold),
    ]:
        new_other = eligible & failure
        other |= new_other
        eligible &= ~failure
        rows.append(
            {
                "stage": label,
                "eligible": int(eligible.sum()),
                "Other": int(other.sum()),
                "new_Other": int(new_other.sum()),
                "Unsafe": int(unsafe.sum()),
                "Unassessed": int(unassessed.sum()),
            }
        )
    assert eligible.equals(facts.category.eq("CTA"))
    assert other.equals(facts.category.eq("Other"))
    assert unsafe.equals(facts.category.eq("Unsafe"))
    assert unassessed.equals(facts.category.isna())
    table = pd.DataFrame(rows)
    assert table[["eligible", "Other", "Unsafe", "Unassessed"]].sum(axis=1).eq(len(facts)).all()
    output.mkdir(parents=True, exist_ok=True)
    table.to_csv(output / "category-funnel.csv", index=False)
    colors = {
        "eligible": "#087E8B",
        "Unsafe": "#AC493B",
        "Other": "#B18335",
        "Unassessed": "#8995A4",
    }
    ink, muted = "#173B55", "#536675"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(15, 8.8), facecolor="white")
    fig.text(
        0.055, 0.955, "Where proteins become CTA or Other", fontsize=26, weight="bold", color=ink
    )
    fig.text(
        0.055,
        0.918,
        f"OncoRef {__version__}  |  {len(facts):,} representative proteins  |  Normal assessment first",
        fontsize=12,
        color=muted,
    )
    counts = {
        "CTA": int(eligible.sum()),
        "Unsafe": int(unsafe.sum()),
        "Other": int(other.sum()),
        "Unassessed": int(unassessed.sum()),
    }
    for i, (name, count) in enumerate(counts.items()):
        x = 0.055 + i * 0.232
        color = colors["eligible" if name == "CTA" else name]
        fig.add_artist(
            FancyBboxPatch(
                (x, 0.762),
                0.218,
                0.121,
                transform=fig.transFigure,
                boxstyle="round,pad=0.009,rounding_size=0.008",
                facecolor="#F3F6F8",
                edgecolor="none",
            )
        )
        fig.text(x + 0.014, 0.84, f"{count:,}", fontsize=25, color=color, weight="bold")
        fig.text(x + 0.014, 0.795, name, fontsize=13, color=ink)
    fig.text(
        0.055,
        0.709,
        f"Follow the {normal_pass.sum():,} proteins that pass normal safety and coverage",
        fontsize=17,
        color=ink,
        weight="bold",
    )
    fig.text(
        0.055,
        0.676,
        f"{unsafe.sum():,} Unsafe and {unassessed.sum():,} Unassessed remain outside this candidate track at every stage below.",
        fontsize=12,
        color=muted,
    )
    ax = fig.add_axes([0.285, 0.265, 0.59, 0.35])
    yy = np.arange(len(table))
    ax.barh(yy, table.eligible, color=colors["eligible"], height=0.61)
    ax.barh(yy, table.Other, left=table.eligible, color=colors["Other"], height=0.61)
    ax.set_yticks(yy, table.stage, fontsize=12, color=ink)
    ax.tick_params(axis="y", length=0, pad=15)
    ax.set_xlim(0, int(normal_pass.sum()))
    ax.set_ylim(len(table) - 0.4, -0.6)
    ax.set_xticks([0, 250, 500, 750, 1000, int(normal_pass.sum())])
    ax.tick_params(axis="x", colors=muted, labelsize=10, length=0, pad=8)
    ax.set_xlabel("Proteins in the normal-pass subset", color=muted, fontsize=11, labelpad=10)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i, row in table.iterrows():
        ax.text(
            row.eligible / 2,
            i,
            f"{row.eligible:,}",
            va="center",
            ha="center",
            color="white",
            weight="bold",
            fontsize=15,
        )
        if row.Other:
            ax.text(
                row.eligible + row.Other / 2,
                i,
                f"{row.Other:,}",
                va="center",
                ha="center",
                color="white",
                weight="bold",
                fontsize=15,
            )
        ax.text(
            1.035,
            i,
            "" if i == 0 else f"+{row.new_Other:,} Other",
            transform=ax.get_yaxis_transform(),
            va="center",
            color=ink,
            fontsize=11,
        )
    fig.text(
        0.29, 0.63, "CTA-eligible so far", fontsize=11, color=colors["eligible"], weight="bold"
    )
    fig.text(0.69, 0.63, "Other (cumulative)", fontsize=11, color=colors["Other"], weight="bold")
    fig.text(
        0.055,
        0.152,
        f"{table.iloc[1].new_Other:,} lack qualifying reproductive evidence; {table.iloc[-1].new_Other:,} additional proteins fail concentration.",
        fontsize=14,
        color=ink,
        weight="bold",
    )
    fig.text(
        0.055,
        0.112,
        f"Cancer presence moves {table.iloc[2].new_Other:,} additional proteins into Other; {table.iloc[2].eligible:,} remain eligible.",
        fontsize=11,
        color=muted,
    )
    fig.text(
        0.055,
        0.066,
        "Reproductive includes placental and embryonic/fetal contexts. Other requires a complete, passing normal screen.\nThe figure changes display order only; all final categories exactly match the released classification.",
        fontsize=10,
        color=muted,
        linespacing=1.6,
    )
    for extension in ["png", "svg"]:
        fig.savefig(output / ("category-funnel." + extension), dpi=180, facecolor="white")
    plt.close(fig)
    (output / "manifest.json").write_text(
        json.dumps(
            {
                "source": str(source.relative_to(root)),
                "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "policy": facts.policy.unique().tolist(),
                "final_counts": counts,
                "normal_first_display_only": True,
            },
            indent=2,
        )
        + "\n"
    )
    print(table.to_string(index=False))


if __name__ == "__main__":
    build(Path.cwd(), Path("docs/audits/cta-source-walkthrough-20261028/figures"))

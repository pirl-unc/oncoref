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

from oncoref.cta_selection_policy import (
    CANDIDATE_ONLY_LAYERS,
    CONTAINED_FLOOR,
    translation_admission,
)
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
    # Proteins with a held, incomplete or unresolved normal screen are Other:Held
    # or Other by gate; they stay outside the candidate track but are not rejected.
    held = ~unsafe & ~normal_pass
    layer = "unsafe:Rare vital-tissue cell expression"
    other_unsafe = facts[[c for c in facts if c.startswith("unsafe:") and c != layer]].any(axis=1)
    rare_only = int((facts[layer] & ~other_unsafe).sum())
    eligible = normal_pass.copy()
    other = pd.Series(False, index=facts.index)
    rows = []
    for label, failure in [
        ("Normal safety + coverage", pd.Series(False, index=facts.index)),
        (
            "Reproductive/developmental presence",
            ~facts.positive_context | facts.nested_microprotein,
        ),
        ("Cancer presence", ~facts.cancer_association),
        ("Reproductive concentration", facts.reproductive_concentration_hold),
        ("Translation qualification", ~translation_admission(facts)),
        (
            "Unsafe sequence containment",
            facts.unsafe_sequence_overlap_pct.ge(100 * CONTAINED_FLOOR),
        ),
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
                "Unresolved": int(held.sum()),
            }
        )
    assert eligible.equals(facts.category.eq("CTA"))
    assert (other | held).equals(facts.category.eq("Other"))
    assert unsafe.equals(facts.category.eq("Unsafe"))
    assert not facts.category.isna().any()
    table = pd.DataFrame(rows)
    assert table[["eligible", "Other", "Unsafe", "Unresolved"]].sum(axis=1).eq(len(facts)).all()
    output.mkdir(parents=True, exist_ok=True)
    table.to_csv(output / "category-funnel.csv", index=False)
    colors = {
        "eligible": "#087E8B",
        "Unsafe": "#AC493B",
        "Other": "#B18335",
        "Unresolved": "#8995A4",
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
        "Unresolved": int(held.sum()),
    }
    sublabels = {
        "Other": "passed the normal screen, failed a target gate",
        "Unresolved": "normal screen unresolved; also Other, not rejected",
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
        if name in sublabels:
            fig.text(x + 0.014, 0.768, sublabels[name], fontsize=7.5, color=muted)
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
        f"{unsafe.sum():,} Unsafe and {held.sum():,} Other proteins with an unresolved normal screen stay outside this track; none is unassessed.",
        fontsize=12,
        color=muted,
    )
    ax = fig.add_axes([0.285, 0.265, 0.59, 0.35])
    yy = np.arange(len(table))
    ax.barh(yy, table.eligible, color=colors["eligible"], height=0.61)
    ax.barh(yy, table.Other, left=table.eligible, color=colors["Other"], height=0.61)
    labels = table.stage.str.replace(
        "Reproductive/developmental presence", "Reproductive/developmental\npresence", regex=False
    )
    ax.set_yticks(yy, labels, fontsize=12, color=ink)
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
        0.178,
        f"{table.iloc[1].new_Other:,} lack reproductive evidence / are nested; {table.iloc[3].new_Other:,} fail concentration; {table.iloc[4].new_Other:,} need translation review; {table.iloc[5].new_Other:,} have Unsafe sequence containment.",
        fontsize=12,
        color=ink,
        weight="bold",
    )
    fig.text(
        0.055,
        0.142,
        f"Cancer presence moves {table.iloc[2].new_Other:,} additional proteins into Other; {table.iloc[2].eligible:,} remain eligible.",
        fontsize=11,
        color=muted,
    )
    fig.text(
        0.055,
        0.108,
        f"Rare vital-tissue cell expression (at least 1% in a named Vital population in the existing brain/spinal-cord/heart/lung panel) is the only Unsafe evidence for {rare_only:,} proteins.",
        fontsize=10,
        color=muted,
    )
    fig.text(
        0.055,
        0.05,
        "Reproductive includes placental and embryonic/fetal contexts. Unresolved proteins pass no CTA gate but are not Unsafe; Other:Held is the subset passing every target gate.\nThe figure changes display order only; all final categories exactly match the released classification.",
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
                "rare_vital_only_unsafe": rare_only,
            },
            indent=2,
        )
        + "\n"
    )
    print(table.to_string(index=False))


if __name__ == "__main__":
    build(Path.cwd(), Path("docs/audits/cta-source-walkthrough-20261028/figures"))

#!/usr/bin/env python3
"""Render the human/dog source walkthrough from the shared selection audit API."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import plot_cta_geometric as style  # noqa: E402
from build_cta_selection_audit import COMPARATORS, OUTPUT, RELEASE, records  # noqa: E402

from oncoref import cta_selection  # noqa: E402


def markdown(frame: pd.DataFrame) -> str:
    """Write a compact Markdown table without an optional tabulate dependency."""

    def cell(value):
        if isinstance(value, float) and np.isnan(value):
            return "—"
        return str(value).replace("|", " / ").replace("\n", " ")

    return "\n".join(
        [
            "| " + " | ".join(map(cell, frame.columns)) + " |",
            "| " + " | ".join(["---"] * len(frame.columns)) + " |",
            *[
                "| " + " | ".join(map(cell, row)) + " |"
                for row in frame.itertuples(index=False, name=None)
            ],
        ]
    )


def figures(facts: pd.DataFrame, audit) -> None:
    """Use citation labels, explicit decision thresholds and unshifted RNA axes."""
    style.setup()
    out = OUTPUT / "figures"
    out.mkdir(exist_ok=True)
    colors = {"CTA": style.TEAL, "Unsafe": style.RED, "Neither": "#A28C68"}
    fig, ax = plt.subplots(figsize=(15, 9))
    fig.subplots_adjust(left=0.29, right=0.77, top=0.82, bottom=0.16)
    style.title(
        fig,
        "Human selection funnel",
        "19,863 protein-equivalence groups; each layer can change the category",
    )
    left = np.zeros(len(audit.funnel))
    for key, label in [
        ("eligible", "Eligible so far"),
        ("Unsafe", "Unsafe"),
        ("Neither", "Neither"),
    ]:
        ax.barh(
            audit.funnel.stage,
            audit.funnel[key],
            left=left,
            color=colors.get(key, style.TEAL),
            label=label,
        )
        left += audit.funnel[key].to_numpy()
    ax.invert_yaxis()
    ax.set_xlabel("Protein-equivalence groups")
    ax.legend(frameon=False, ncol=3, loc="upper left", bbox_to_anchor=(0, -0.08), fontsize=9)
    ax.text(1.02, 1.025, "Eligible   Unsafe   Neither", transform=ax.transAxes, fontsize=9)
    for i, row in enumerate(audit.funnel.itertuples()):
        ax.text(
            1.02,
            i,
            f"{row.eligible:>6,}  {row.Unsafe:>6,}  {row.Neither:>6,}",
            transform=ax.get_yaxis_transform(),
            va="center",
            fontfamily="monospace",
            fontsize=9,
        )
    style.footer(
        fig,
        "Current policy: 463 CTA / 17,388 Unsafe / 2,012 Neither. Intermediate eligibility is not a final CTA assignment.\nSequential removals depend on display order; final classifications do not.",
    )
    style.save(fig, out, "human-funnel")

    fig, ax = plt.subplots(figsize=(12, 7))
    fig.subplots_adjust(left=0.24, right=0.96, top=0.79, bottom=0.24)
    style.title(fig, "Active RNA cutoffs", "Candidate ceilings under the released 1.21.0 policy")
    grid = np.array([[base * 2**i for i in range(4)] for base in (4, 8, 16, 32)])
    ax.imshow(np.log2(grid), cmap="YlGnBu", aspect="auto")
    ax.set_xticks(
        range(4),
        ["Unmatched / unusable", "Approved negative", "Supported negative", "Enhanced negative"],
    )
    ax.set_yticks(range(4), ["Vital", "Protected", "Manageable", "Specified bowel epithelium"])
    for i, j in np.ndindex(grid.shape):
        ax.text(
            j,
            i,
            str(grid[i, j]),
            ha="center",
            va="center",
            fontsize=20,
            color="white" if grid[i, j] >= 32 else style.INK,
        )
    style.footer(
        fig,
        "At/below the ceiling passes; above it holds. The released Unsafe RNA threshold is 2× each value.\nTPM/nTPM retain native assay units. Corneal CP10K uses nominal limits /100, not a TPM conversion.\nSame-tissue negative IHC is shown separately when cell-type coverage is unresolved.",
    )
    style.save(fig, out, "active-rna-cutoffs")

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    fig.subplots_adjust(left=0.075, right=0.97, top=0.80, bottom=0.20, wspace=0.21)
    style.title(
        fig,
        "Widely studied antigens",
        "Reproductive and developmental routes remain separate measurements",
    )
    from adjustText import adjust_text

    comparators = facts[facts.protein_group_label.isin(COMPARATORS)]
    for ax, column, label in zip(
        axes,
        ["reproductive_ntpm", "developmental_mean_tpm"],
        ["Reproductive RNA (nTPM)", "Embryonic lineage RNA (TPM)"],
        strict=True,
    ):
        valid = facts.adult_bulk_max_ntpm.notna() & facts[column].notna()
        for category, color in colors.items():
            part = facts[valid & facts.category.eq(category)]
            ax.scatter(
                part.adult_bulk_max_ntpm,
                part[column],
                s=8,
                color=color,
                alpha=0.22,
                rasterized=True,
                label=f"{category}: {len(part):,} plotted",
            )
        ax.set_xscale("symlog", linthresh=0.1)
        ax.set_yscale("symlog", linthresh=0.1)
        labels = []
        for row in comparators[comparators[column].notna()].to_dict("records"):
            x, y = row["adult_bulk_max_ntpm"], row[column]
            ax.scatter(x, y, color=colors[row["category"]], s=38, edgecolor=style.INK, zorder=3)
            display = row["protein_group_label"].replace("MAGEA", "MAGE-A")
            labels.append(ax.text(x, y, display, fontsize=9))
        adjust_text(
            labels,
            ax=ax,
            iter_lim=200,
            arrowprops={"arrowstyle": "-", "color": style.MUTED, "lw": 0.6},
        )
        ax.set_xlabel("Adult bulk RNA (nTPM)")
        ax.set_ylabel(label)
        ax.legend(frameon=False, fontsize=8, loc="upper left")
    style.footer(
        fig,
        "Axes include zero on a symmetric-log scale. Only groups with measurements are plotted; colors use the complete funnel.\nTBXT enters through the developmental route. MAGE-A1, MAGE-A3 and MAGE-A12 are comparators, not forced positive controls.",
    )
    style.save(fig, out, "widely-studied-antigens")

    ranked = audit.sources.sort_values("mapped_entities", ascending=False)
    fig, ax = plt.subplots(figsize=(13, 8.5))
    fig.subplots_adjust(left=0.25, right=0.94, top=0.82, bottom=0.18)
    style.title(
        fig, "Published candidate lists", "Tables from the same paper are unioned before counting"
    )
    ax.barh(ranked.citation, ranked.mapped_entities, color=style.PALE, label="Mapped groups")
    ax.barh(ranked.citation, ranked.CTA, color=style.TEAL, label="Selected CTA groups")
    ax.set_xlim(0, 1220)
    for i, row in enumerate(ranked.itertuples()):
        ax.text(
            row.mapped_entities + 15,
            i,
            f"{row.CTA} / {row.mapped_entities:,}",
            va="center",
            fontsize=9,
        )
    ax.invert_yaxis()
    ax.set_xlabel("Protein-equivalence groups")
    ax.legend(frameon=False, ncol=2, loc="upper left", bbox_to_anchor=(0, -0.07), fontsize=10)
    style.footer(
        fig,
        "Union: 2,482 mapped groups / 2,538 mapped loci; 418 selected groups have a registered list annotation.\n45 selected groups have no annotation in these lists. Counts overlap; list membership is not an eligibility requirement.",
    )
    style.save(fig, out, "published-list-contributions")

    top = ranked.head(6)
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    fig.subplots_adjust(left=0.18, right=0.96, top=0.77, bottom=0.28, wspace=0.28)
    style.title(
        fig,
        "Overlap between published lists",
        "The six largest mapped paper sets; full pairwise tables are downloadable",
    )
    short = top.citation.str.replace(" et al.", "", regex=False).tolist()
    for ax, field, title in zip(
        axes,
        ["mapped_overlap", "selected_overlap"],
        ["All mapped groups", "Selected CTA groups"],
        strict=True,
    ):
        matrix = (
            audit.overlaps.pivot(index="source_a", columns="source_b", values=field)
            .reindex(index=top.source_id, columns=top.source_id)
            .to_numpy(dtype=int)
        )
        ax.imshow(matrix, cmap="YlGnBu", aspect="equal")
        ax.set_xticks(range(6), short, rotation=40, ha="right", fontsize=10)
        ax.set_yticks(range(6), short, fontsize=10)
        ax.set_title(title, fontsize=13)
        for i, j in np.ndindex(matrix.shape):
            ax.text(
                j,
                i,
                str(matrix[i, j]),
                ha="center",
                va="center",
                fontsize=11,
                color="white" if matrix[i, j] > matrix.max() * 0.55 else style.INK,
            )
    style.footer(
        fig,
        "Diagonal cells give each paper's set size. Off-diagonal cells count shared groups, not independent validation studies.\nCitations and original source-table IDs remain linked in the source inventory.",
    )
    style.save(fig, out, "published-list-overlaps")

    fig, ax = plt.subplots(figsize=(12.5, 6.5))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.77, bottom=0.18)
    style.title(
        fig,
        "Dog evidence funnel",
        "14 native nominations; quantitative selection has not been evaluated",
    )
    ax.axis("off")
    values = [
        ["Nominations", "5 + 10 − 1 = 14", "NT5C1B is shared between the two studies"],
        ["Identity", "13 unique; 1 family", "MAGE assay remains ambiguous"],
        [
            "Normal / developmental RNA",
            "Native matrices available",
            "Reference reconciliation and scale policy needed",
        ],
        ["Tumor RNA", "Cohort blocked", "Count/FPKM and annotation conflict"],
        ["Normal IHC / MS", "Limited IHC; pooled iBAQ", "Not a complete tissue-cell panel"],
        ["Selection", "14 Neither / unassessed", "No passing/failing count is inferred"],
    ]
    table = ax.table(
        cellText=values,
        colLabels=["Stage", "Available evidence", "Interpretation"],
        colWidths=[0.25, 0.25, 0.50],
        cellLoc="left",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2.1)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor(style.PALE)
        cell.set_facecolor(style.INK if row == 0 else style.BG)
        cell.set_text_props(color="white" if row == 0 else style.INK)
    style.footer(
        fig,
        "Hindriks et al. (2025): eight original assay targets plus two MAGE follow-ups. Nemec et al. (2019): five curated nominees.\nMissing canine evidence does not inherit human HPA grades, thresholds, or normal-tissue absence.",
    )
    style.save(fig, out, "dog-funnel")


def bowel_comparison(facts: pd.DataFrame) -> pd.DataFrame:
    """Evaluate only the monotonic change from the bowel 32 row to Manageable 16."""
    contexts = pd.read_csv(RELEASE / "changed-compartment-contexts.csv.gz")
    bowel = contexts[contexts.selection_tier.eq("renewing_mucosa")].copy()
    bowel["hold"] = bowel.value.gt(bowel.candidate_ceiling / 2)
    bowel["unsafe"] = bowel.value.ge(bowel.candidate_ceiling)
    changed = bowel.groupby("gene_id")[["hold", "unsafe"]].any()
    result = facts.copy()
    result["hold:HPA tissue clusters"] |= result.representative_gene_id.map(changed.hold).fillna(
        False
    )
    result["unsafe:HPA tissue clusters"] |= result.representative_gene_id.map(
        changed.unsafe
    ).fillna(False)
    held = result.filter(regex="^hold:").any(axis=1)
    unsafe = result.filter(regex="^unsafe:").any(axis=1)
    result["comparison_category"] = np.select([unsafe, held], ["Unsafe", "Neither"], default="CTA")
    delta = result[result.category.ne(result.comparison_category)]
    delta[
        ["protein_equivalence_group_id", "protein_group_label", "category", "comparison_category"]
    ].to_csv(OUTPUT / "bowel-manageable16-comparison.csv", index=False)
    return result


def build() -> None:
    """Render documents and plots against the current, unchanged selection snapshot."""
    audit = cta_selection.audit()
    dog = cta_selection.audit(species="dog")
    facts = pd.read_csv(RELEASE / "protein-group-classifications.csv.gz", low_memory=False)
    comparison = bowel_comparison(facts)
    figures(facts, audit)
    known = facts[facts.protein_group_label.isin(COMPARATORS)].copy()
    known["detailed_reason"] = known.protein_equivalence_group_id.map(
        audit.explanations[audit.explanations.decision.isin(["hold", "unsafe"])]
        .groupby("entity_id")
        .reason.agg(" | ".join)
    ).fillna("Passes the released rules; retained observations remain in the detailed evidence")
    known[
        [
            "protein_equivalence_group_id",
            "protein_group_label",
            "category",
            "detailed_reason",
            "reproductive_ntpm",
            "developmental_mean_tpm",
        ]
    ].to_csv(OUTPUT / "widely-studied-antigens.csv", index=False)
    paper = audit.sources.copy()
    paper["Paper"] = paper.apply(lambda r: f"[{r.citation}]({r.source_url})", axis=1)
    paper = paper.rename(
        columns={
            "mapped_entities": "Mapped groups",
            "mapped_gene_loci": "Gene loci",
            "unique_entities": "Unique to this paper",
        }
    )
    association = []
    for name, mask in [
        ("Reproductive only", facts.reproductive_route & ~facts.developmental_route),
        ("Developmental only", ~facts.reproductive_route & facts.developmental_route),
        ("Both", facts.reproductive_route & facts.developmental_route),
    ]:
        association.append(
            {
                "Route": name,
                "Entrants": int(mask.sum()),
                "Final CTA": int((mask & facts.category.eq("CTA")).sum()),
            }
        )
    sources = cta_selection.sources()
    used = sources[sources.status.eq("used_in_selection")].copy()
    used["Source"] = used.apply(lambda r: f"[{r['name']}]({r.source_url})", axis=1)
    used = used.rename(
        columns={
            "evidence_category": "Role",
            "native_unit": "Unit",
            "limitations": "Scope / limitations",
        }
    )
    overview = f"""# Human and dog CTA evidence walkthrough

This walkthrough explains the **released 1.21.0 decisions** with more specific reasons and consistent labels. It does not change the thresholds or claim a completed canine quantitative selection. [Dog funnel](dog.md) · [Policy questions and developmental sources](policy-review.md) · [Shared API](api.md).

[Download the 17-slide evidence walkthrough](slides/cta-evidence-walkthrough.pptx).

## Published candidate sources

Counts below are OncoRef mappings, not the papers' original row counts. Tables from one paper are unioned; distinct papers can share genes. Protein-equivalence groups and gene loci remain separate denominators. The union is **2,482 groups / 2,538 loci**. Of the 463 selected groups, **418** have a registered paper annotation and **45** do not. A paper annotation does not bypass any quantitative gate. Boël's unresolved BAGE record has no mapped eligible group in this registry.

{markdown(paper[["Paper", "Gene loci", "Mapped groups", "CTA", "Unsafe", "Neither", "Unique to this paper"]])}

![Published source contributions](figures/published-list-contributions.png)
![Paper overlaps](figures/published-list-overlaps.png)

[All source memberships](human-memberships.csv.gz) · [All pairwise overlaps](human-overlaps.csv.gz). Unique means present in only that paper within this registry, not uniquely discovered or independently validated.

## Why 12,139 groups enter

The starting universe is **19,863 protein-equivalence groups / 20,097 loci**. The first gate uses measured RNA across this universe, not the 2,482-group literature union. It accepts reproductive maximum ≥5 nTPM and >adult bulk maximum **or** a supported Tyser embryonic lineage mean ≥5 TPM.

{markdown(pd.DataFrame(association))}

The developmental criterion tests expression, not restriction. Consequently it admits many broadly expressed genes before normal-tissue filtering. Tyser measures one Carnegie-stage-7 embryo; later fetal atlases are not integrated. Reproductive and developmental routes share a permissive-context role but retain distinct stage, tissue and units. TBXT has reproductive 0.2 nTPM, adult bulk maximum 2.5 nTPM, and supported embryonic lineage mean 245.088 TPM: it enters through the developmental route.

## Normal and tumor evidence

{markdown(used[["Role", "Source", "Unit", "Scope / limitations"]])}

The OncoRef tumor union uses **132 compatible cohort views of 143 registered**, with **17,763 sample occurrences**, **14,452 namespace-deduplicated sample IDs**, and **34 source namespaces**. These are not 14,452 proven independent patients. HPA cancer can overlap the underlying cohorts. [Per-cohort availability and provenance](../cta-expanded-funnel-20261027/tumor-union/cohort-availability.csv).

## Active criteria and funnel

{markdown(audit.criteria[["stage", "rule"]].rename(columns={"stage": "Stage", "rule": "Criterion"}))}

![Active RNA cutoffs](figures/active-rna-cutoffs.png)

“Unsafe threshold” is the released second RNA boundary, **2× the candidate ceiling**. It is not an IHC intensity category. Above the candidate ceiling but below that boundary creates a hold; another source can independently make the final decision Unsafe. These are the active policy cutoffs, not clinically measured toxicity limits.

{markdown(audit.funnel[["stage", "candidate_removals", "eligible", "Unsafe", "Neither"]].rename(columns={"stage": "Stage", "candidate_removals": "New removals", "eligible": "Eligible so far"}))}

![Human filtering funnel](figures/human-funnel.png)

## Widely studied antigens

MAGE-A1, MAGE-A3 and MAGE-A12 are comparison genes, not forced positive controls. Detailed explanations retain the source, tissue, cell type, value, units, threshold and evidence uncertainty.

{markdown(known[["protein_group_label", "category", "detailed_reason"]].rename(columns={"protein_group_label": "Antigen / group", "category": "Decision", "detailed_reason": "Reason"}))}

![Widely studied antigen distributions](figures/widely-studied-antigens.png)

[All detailed human decisions](human-explanations.csv.gz) · [Comparison table](widely-studied-antigens.csv). Each RNA layer links a representative supported context and the complete source export; it does not claim that this is the only normal expression.

## Folding bowel into Manageable

Moving specified bowel epithelium from its current base ceiling 32 to the existing Manageable base 16, preserving all other released rules, gives **{int(comparison.comparison_category.eq("CTA").sum())} CTA / {int(comparison.comparison_category.eq("Unsafe").sum()):,} Unsafe / {int(comparison.comparison_category.eq("Neither").sum()):,} Neither**. PRAME moves to Neither; MAGE-A3, HOXA13 and FAM72B move from Neither to Unsafe. This is a comparison, not a replacement release. [Every changed decision](bowel-manageable16-comparison.csv).

PRAME has Enhanced “Not detected” HPA colon glandular, endothelial and peripheral nerve/ganglion staining. Its 20.3 nTPM signal is in undifferentiated colon cells. The current match requires cell-type coverage, so this is **same-tissue negative IHC with unresolved cell-type coverage**, not absent IHC. Assigning Enhanced credit to the undifferentiated cells requires evidence that the stained compartment covers them. [HPA PRAME tissue evidence](https://v23.proteinatlas.org/ENSG00000185686-PRAME/tissue/colon).
"""
    (OUTPUT / "index.md").write_text(overview)
    dog_sources = cta_selection.sources(species="dog").copy()
    dog_sources["Source"] = dog_sources.apply(lambda r: f"[{r['name']}]({r.source_url})", axis=1)
    dog_sources = (
        dog_sources.groupby(["Source", "status", "native_unit"], sort=False)
        .agg(
            Roles=("evidence_category", lambda x: ", ".join(x)),
            Limitations=("limitations", "first"),
        )
        .reset_index()
    )
    dog_doc = f"""# Dog CTA sources and evidence funnel

Dogs use the same source categories, identity contracts and audit schema as humans. **Quantitative selection is not yet evaluated.** The current output is **14 Neither / unassessed nominations**, not a genome-wide statement that dogs have zero CTAs or zero unsafe genes.

## Candidate studies and overlap

{markdown(dog.sources[["citation", "mapped_entities", "mapped_gene_loci", "unique_entities"]].rename(columns={"citation": "Study", "mapped_entities": "Nominations", "mapped_gene_loci": "Mapped loci / alternatives", "unique_entities": "Not in the other study"}))}

[Nemec et al. (2019)](https://doi.org/10.1111/vco.12475) contributes five curated nominees: AKR1E2, IGF2BP3, NT5C1B, SPECC1 and TPX2. [Hindriks et al. (2025)](https://doi.org/10.1016/j.vetimm.2025.110946) contributes eight assay targets and two MAGE follow-ups: CAGE1, CTCFL, DDX53, MAGE, MAGEA13P, MAGEB10, NT5C1B, PAGE3-like, PRAME and SSX6. **NT5C1B is the one shared nomination: 5 + 10 − 1 = 14.** The MAGE family assay overlaps its gene follow-ups; it is not an additional resolved gene. Thirteen nomination identities are unique; one is an ambiguous family assay.

The snapshot retains 209 observations, with 10 nominations described as conflicting and four as limited. These source-evidence statuses do not replace a quantitative CTA/Unsafe decision. Nemec is abstract-level curation; its grouped normal-tissue statement cannot allocate each tissue to an individual gene. Hindriks retains numeric PCR percentages and reported categorical signs separately because some disagree.

## Normal, developmental, tumor and protein resources

{markdown(dog_sources.rename(columns={"status": "Integration status", "native_unit": "Native unit"}))}

BarkBase's local matrix has **37,149 features × 150 TMM-CPM columns**. Bgee/BarkBase reuse 46 libraries, which must not count as independent dogs. Pfizer protein has **9,929 reported groups × 38 tissue pools**, including **64,047 missing matrix cells**. These are native features/protein groups, not reconciled CTA candidate counts. A pooled protein signal does not establish detection in each donor or identify a cell type.

## Funnel status

{markdown(dog.funnel[["stage", "status", "note"]].rename(columns={"stage": "Stage", "status": "Status", "note": "Evidence / gap"}))}

![Dog evidence funnel](figures/dog-funnel.png)

No human RNA cutoff, HPA antibody grade or clinical finding transfers by orthology. The eventual dog filter can use the same logical gates, with native source adapters, explicit developmental stages, canine tissue panels and unit-aware thresholds. Until those inputs are reconciled, per-stage pass/fail numbers are **unknown**, not zero.

## Tracked work

- [OncoRef #574](https://github.com/pirl-unc/oncoref/issues/574): Bgee/BarkBase, complete DoGA and Pfizer RNA, reference reconciliation, tumor cohorts and source/sample lineage.
- [OncoRef #585](https://github.com/pirl-unc/oncoref/issues/585): consume the shared Pfizer canine bulk-MS artifact; richer protein/group assignments and missingness.
- [orfdata #7](https://github.com/openvax/orfdata/issues/7), [#23](https://github.com/openvax/orfdata/issues/23), [#24](https://github.com/openvax/orfdata/issues/24): generic atlas extraction and compact versioned source data.
- [Hitlist #666](https://github.com/pirl-unc/hitlist/issues/666), [#668](https://github.com/pirl-unc/hitlist/issues/668): peptide/DLA-facing coverage and consumer integration. Bulk MS is not evidence of normal DLA presentation.

No dedicated OncoRef issue for systematic expansion of canine IHC was found in searches of issue titles and canine/IHC issue text. The missing work is a source/reagent/tissue-cell inventory with measured staining, validation evidence, cross-reactivity and explicit negative/missing states; it should not inherit HPA grades.

[Shared API](api.md) · [All native source memberships](dog-memberships.csv.gz) · [Human walkthrough](index.md).
"""
    (OUTPUT / "dog.md").write_text(dog_doc)
    presentation = json.loads((RELEASE / "presentation-data.json").read_text())
    presentation["named"] = records(known)
    presentation["papers"] = [
        {
            "paper": r.citation,
            "doi": r.source_id,
            "mapped_groups": r.mapped_entities,
            "CTA": r.CTA,
            "Unsafe": r.Unsafe,
            "Neither": r.Neither,
        }
        for r in audit.sources.itertuples()
    ]
    presentation["overlap_order"] = audit.sources.head(6).citation.tolist()
    presentation["funnel"][0]["layer"] = "Reproductive or developmental expression"
    presentation["dog_funnel"] = records(dog.funnel)
    (OUTPUT / "presentation-data.json").write_text(json.dumps(presentation, indent=2) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    build()

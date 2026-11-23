#!/usr/bin/env python3
"""Reclassify gene loci and protein equivalence groups and rebuild the CTA report."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import matplotlib.pyplot as plt  # noqa: E402
import plot_cta_geometric as style  # noqa: E402
from cta_expanded_funnel import ASSAYS, OPTIONAL, markdown_table  # noqa: E402
from cta_report_common import (  # noqa: E402
    input_hashes,
    sha256,
    write_csv,
)

from oncoref import cell_rna_policy, cta_tissues  # noqa: E402
from oncoref.cta_funnel import EvidenceLayer, classify_layers  # noqa: E402
from oncoref.external_ihc import evidence as external_evidence  # noqa: E402
from oncoref.gene_matching import attach_gene_identities  # noqa: E402
from oncoref.protein_equivalence import equivalence_groups  # noqa: E402

NAMED = [
    "PRAME",
    "MAGEA4",
    "CTAG1A",
    "CTAG1B",
    "XAGE1A",
    "XAGE1B",
    "CTAG2",
    "SSX2",
    "CT83",
    "TBXT",
    "SMC1B",
    "PAGE1",
]
POLICY_ID = "oncoref.harmonized-cta.2026-10-27.2"


def rebuild_groups(facts: pd.DataFrame, baseline: Path, output: Path) -> pd.DataFrame:
    """Use measured group sums for RNA; preserve member protein and safety evidence."""
    groups = equivalence_groups(facts.gene_id)
    members = groups[["protein_equivalence_group_id", "member_gene_ids"]].copy()
    members["gene_id"] = members.member_gene_ids.map(json.loads)
    members = members.explode("gene_id").drop(columns="member_gene_ids")
    joined = members.merge(facts, on="gene_id", how="left", validate="many_to_one")
    base_columns = [
        c
        for c in facts
        if not c.startswith(("hold:", "unsafe:"))
        and c not in ["category", "exclusion_reasons", "Unsafe_reasons", "identity_catalog_sha256"]
    ]
    # Representative fields initialize singleton rows only. Every multi-locus
    # measurement or decision used below is replaced by its explicit group rule.
    result = groups.merge(
        facts[base_columns],
        left_on="representative_gene_id",
        right_on="gene_id",
        how="left",
        validate="one_to_one",
    )
    result["gene_id"] = result.representative_gene_id
    result["symbol"] = result.protein_group_label
    indexed = joined.groupby("protein_equivalence_group_id")
    for name in [
        "ihc_detected",
        "ihc_reliable_detected",
        "Accepted normal protein",
        "Safety veto",
        "Antibody concern",
        "external_review_flag",
        "external_eligible_positive",
    ]:
        result[name] = result.protein_equivalence_group_id.map(indexed[name].any()).astype(bool)
    required = indexed.required_coverage.agg(lambda s: s.notna().all() and s.all())
    result["required_coverage"] = result.protein_equivalence_group_id.map(required).astype(bool)
    for name in ["optional_DICE_complete", "optional_cornea_complete"]:
        complete = indexed[name].agg(lambda s: s.notna().all() and s.all())
        result[name] = result.protein_equivalence_group_id.map(complete).astype(bool)
    result["member_gene_categories"] = result.protein_equivalence_group_id.map(
        indexed.apply(
            lambda d: json.dumps(
                dict(zip(d.gene_id, d.category.fillna("outside_reference"))), sort_keys=True
            )
        )
    )
    result["member_named_controls"] = result.protein_equivalence_group_id.map(
        indexed.symbol.agg(lambda s: ";".join(sorted(set(s.dropna()) & set(NAMED))))
    )
    multi = result.member_count.gt(1)
    normal = pd.read_csv(output / "group-normal/normal-assessments.csv.gz")
    old_normal = pd.read_csv(baseline / "normal/normal-assessments.csv.gz")
    for assay in ASSAYS:
        lookup = old_normal[old_normal.assay.eq(assay)].set_index("gene_id")
        grouped = normal[normal.assay.eq(assay)].set_index("gene_id")
        for target, value, missing in [
            ("hold:", "candidate_exclusion", assay not in OPTIONAL),
            ("unsafe:", "unsafe_evidence", False),
        ]:
            column = target + assay
            result[column] = result.gene_id.map(lookup[value]).fillna(missing).astype(bool)
            result.loc[multi, column] = (
                result.loc[multi, "gene_id"].map(grouped[value]).fillna(missing).astype(bool)
            )
    bulk = pd.read_csv(output / "group-normal/bulk-contexts.csv.gz")
    bulkid = "Ensembl_Gene_ID"
    reproductive = (
        bulk[bulk.tissue.isin(cta_tissues.CORE_REPRODUCTIVE_TISSUES)].groupby(bulkid).rna_ntpm.max()
    )
    adult = (
        bulk[~bulk.tissue.isin(cta_tissues.EXTENDED_REPRODUCTIVE_TISSUES)]
        .groupby(bulkid)
        .rna_ntpm.max()
    )
    result.loc[multi, "reproductive_ntpm"] = result.loc[multi, "gene_id"].map(reproductive)
    result.loc[multi, "adult_bulk_max_ntpm"] = result.loc[multi, "gene_id"].map(adult)
    complete_bulk = bulk.groupby(bulkid).agg(
        contexts=("tissue", "nunique"), complete=("complete_members", "all")
    )
    result.loc[multi, "required_coverage"] &= (
        result.loc[multi, "gene_id"]
        .map(complete_bulk.contexts.eq(50) & complete_bulk.complete)
        .fillna(False)
    )
    result["reproductive_route"] = result.reproductive_ntpm.ge(5) & result.reproductive_ntpm.gt(
        result.adult_bulk_max_ntpm
    )
    dev = pd.read_csv(output / "group-normal/developmental-lineage-bounds.csv.gz")
    devbest = (
        dev[dev.supported]
        .sort_values("mean_tpm", ascending=False)
        .drop_duplicates("gene_id")
        .set_index("gene_id")
    )
    result.loc[multi, "developmental_route"] = (
        result.loc[multi, "gene_id"].map(dev.groupby("gene_id").admitted.any()).fillna(False)
    )
    result.loc[multi, "developmental_mean_tpm"] = result.loc[multi, "gene_id"].map(devbest.mean_tpm)
    result["positive_context"] = result.reproductive_route | result.developmental_route
    tumor = pd.read_csv(baseline / "tumor-union/protein-group-tumor-maxima.csv.gz").set_index(
        "protein_equivalence_group_id"
    )
    for target, source in [
        ("oncoref_maximum_tpm", "maximum_tpm"),
        ("oncoref_maximum_sample", "maximum_sample"),
        ("oncoref_cancer_code", "cancer_code"),
    ]:
        result.loc[multi, target] = result.loc[multi, "protein_equivalence_group_id"].map(
            tumor[source]
        )
    # One measured member >=2 establishes a group lower bound in that real
    # sample. We never call a sum of unpaired gene maxima a group maximum.
    hpa_lower = indexed.tumor_maximum_max_ptpm.max()
    result["hpa_tumor_group_lower_bound_ptpm"] = result.protein_equivalence_group_id.map(hpa_lower)
    result.loc[multi, "tumor_maximum_max_ptpm"] = np.nan
    result["tumor_union_ge2"] = result.oncoref_maximum_tpm.ge(
        2
    ) | result.hpa_tumor_group_lower_bound_ptpm.ge(2)
    # Do not expose the representative locus's historical measurements or
    # selection flags as if they described the summed group.
    columns = [
        "protein_equivalence_group_id",
        "protein_group_label",
        "proteoform_key",
        "representative_gene_id",
        "member_gene_ids",
        "member_count",
        "member_identities",
        "group_registry_sha256",
        "identity_catalog_sha256",
        "member_gene_categories",
        "member_named_controls",
        "reproductive_ntpm",
        "adult_bulk_max_ntpm",
        "developmental_mean_tpm",
        "reproductive_route",
        "developmental_route",
        "positive_context",
        "required_coverage",
        "optional_DICE_complete",
        "optional_cornea_complete",
        "ihc_detected",
        "ihc_reliable_detected",
        "Accepted normal protein",
        "Safety veto",
        "Antibody concern",
        "external_review_flag",
        "external_eligible_positive",
        "oncoref_maximum_tpm",
        "oncoref_maximum_sample",
        "oncoref_cancer_code",
        "hpa_tumor_group_lower_bound_ptpm",
        "tumor_union_ge2",
        *[c for c in result if c.startswith(("hold:", "unsafe:"))],
    ]
    return result[columns].copy()


def decision_layers(facts: pd.DataFrame) -> list[EvidenceLayer]:
    """Common order-independent rule set for loci and protein equivalence groups."""
    empty = pd.Series(False, index=facts.index)
    return [
        EvidenceLayer("Reproductive / embryonic association", ~facts.positive_context, empty),
        EvidenceLayer("Tumor sample ≥2 TPM/pTPM", ~facts.tumor_union_ge2, empty),
        EvidenceLayer(
            "HPA normal protein",
            facts.ihc_detected,
            facts.ihc_reliable_detected & ~facts["Antibody concern"],
        ),
        EvidenceLayer("Required RNA coverage", ~facts.required_coverage, empty),
        *[
            EvidenceLayer(assay, facts["hold:" + assay], facts["unsafe:" + assay])
            for assay in ASSAYS
        ],
        EvidenceLayer(
            "Accepted normal-protein evidence",
            facts["Accepted normal protein"],
            facts["Accepted normal protein"],
        ),
        EvidenceLayer("Clinical safety exclusion", facts["Safety veto"], facts["Safety veto"]),
        EvidenceLayer(
            "External validated normal protein",
            facts.external_eligible_positive,
            facts.external_eligible_positive,
        ),
    ]


def figures(
    output: Path,
    result: pd.DataFrame,
    funnel: pd.DataFrame,
    independent: pd.DataFrame,
    named: pd.DataFrame,
    sources: pd.DataFrame,
    overlap: pd.DataFrame,
) -> None:
    """Generate figures with current counts and an explicit entity denominator."""
    style.setup()
    figdir = output / "figures"
    figdir.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(16, 8.5))
    fig.subplots_adjust(left=0.24, right=0.77, top=0.82, bottom=0.17)
    style.title(
        fig,
        "CTA selection at every evidence layer",
        f"{len(result):,} protein equivalence groups; evidence remains active after earlier holds",
    )
    left = np.zeros(len(funnel))
    for col, label, color in [
        ("CTA_eligible_so_far", "CTA / eligible so far", style.TEAL),
        ("Unsafe", "Unsafe", style.RED),
        ("Neither", "Neither", "#BEA477"),
    ]:
        ax.barh(funnel.layer, funnel[col], left=left, label=label, color=color, height=0.7)
        left += funnel[col]
    ax.invert_yaxis()
    ax.text(
        1.03,
        1.025,
        "CTA*       Unsafe      Neither",
        transform=ax.transAxes,
        fontsize=9,
        fontfamily="monospace",
        weight="bold",
    )
    for i, row in enumerate(funnel.itertuples()):
        ax.text(
            1.03,
            i,
            f"{row.CTA_eligible_so_far:>6,}     {row.Unsafe:>6,}     {row.Neither:>6,}",
            transform=ax.get_yaxis_transform(),
            va="center",
            fontsize=9,
            fontfamily="monospace",
            weight="bold" if i == len(funnel) - 1 else "normal",
        )
    ax.set_xlabel("Protein equivalence groups · disjoint categories")
    ax.legend(frameon=False, ncol=3, loc="upper left", bbox_to_anchor=(0, -0.08), fontsize=9)
    style.footer(
        fig,
        "Unsafe is a research-policy exclusion, not demonstrated toxicity. Neither includes missing or conflicting evidence.\n*CTA means eligible so far until the final row. Every layer evaluates every group; later evidence can change Neither to Unsafe.",
    )
    style.save(fig, figdir, "classification-funnel")

    fig, ax = plt.subplots(figsize=(12.5, 7))
    fig.subplots_adjust(left=0.28, right=0.92, top=0.8, bottom=0.18)
    style.title(
        fig,
        "How many candidates each choice removes",
        "Bars show sequential removals; the source table also reports overlapping independent effects",
    )
    ax.barh(funnel.layer, funnel.new_candidate_removals, color=style.TEAL)
    for i, n in enumerate(funnel.new_candidate_removals):
        ax.text(n + max(funnel.new_candidate_removals) * 0.01, i, f"{n:,}", va="center", fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, max(funnel.new_candidate_removals) * 1.18)
    ax.set_xlabel("Groups newly removed from CTA eligibility")
    style.footer(
        fig,
        "The order changes which layer gets first-removal credit; it does not change the final CTA / Unsafe / Neither assignments.",
    )
    style.save(fig, figdir, "filter-removals")

    ladder = cell_rna_policy.threshold_ladder("two_axis_doubling")
    table = np.array(
        [
            [ladder[g][tier] for g in ["Missing", "Approved", "Supported", "Enhanced"]]
            for tier in ["vital", "protected", "manageable"]
        ]
    )
    fig, ax = plt.subplots(figsize=(11, 6.5))
    fig.subplots_adjust(left=0.17, right=0.95, top=0.76, bottom=0.25)
    style.title(
        fig,
        "One linked normal-expression threshold grid",
        "Candidate ceilings double across IHC confidence and down anatomical tiers",
    )
    ax.imshow(np.log2(table), cmap="YlGnBu", aspect="auto")
    ax.set_xticks(
        range(4),
        ["Missing / unusable", "Approved negative", "Supported negative", "Enhanced negative"],
    )
    ax.set_yticks(range(3), ["Vital", "Protected", "Manageable"])
    for i, j in np.ndindex(table.shape):
        ax.text(
            j,
            i,
            f"{table[i, j]:g}\nUnsafe ≥{2 * table[i, j]:g}",
            ha="center",
            va="center",
            color="white" if table[i, j] >= 16 else style.INK,
            fontsize=15,
        )
    style.footer(
        fig,
        "Candidate equality passes. RNA above the ceiling and below 2× is Neither; RNA ≥2× is Unsafe.\nNative TPM/nTPM remain separate assays. CP10K uses nominal values /100, not an empirical TPM conversion.\nNegative IHC requires matching tissue/cell coverage, a positive profile control and no unresolved antibody conflict.",
    )
    style.save(fig, figdir, "threshold-grid")

    colors = {"CTA": style.TEAL, "Unsafe": style.RED, "Neither": "#BEA477"}
    fig, ax = plt.subplots(figsize=(12, 7.5))
    fig.subplots_adjust(left=0.12, right=0.96, top=0.8, bottom=0.19)
    style.title(
        fig,
        "Where the named antigens land",
        "Maximum nonpermitted adult bulk versus reproductive RNA · HPA native nTPM",
    )
    valid = result.adult_bulk_max_ntpm.notna() & result.reproductive_ntpm.notna()
    for category, color in colors.items():
        d = result[valid & result.category.eq(category)]
        ax.scatter(
            d.adult_bulk_max_ntpm + 0.05,
            d.reproductive_ntpm + 0.05,
            s=9,
            color=color,
            alpha=0.23,
            label=f"{category} ({int(result.category.eq(category).sum()):,})",
            rasterized=True,
        )
    # Label placement operates in display coordinates, so finalize the axis
    # transforms first; changing to log scale afterward can clip the labels.
    ax.set_xscale("log")
    ax.set_yscale("log")
    labels = []
    for row in named.drop_duplicates("protein_equivalence_group_id").itertuples():
        x, y = row.adult_bulk_max_ntpm + 0.05, row.reproductive_ntpm + 0.05
        ax.scatter([x], [y], s=45, color=colors[row.category], edgecolor=style.INK, zorder=4)
        labels.append(ax.text(x, y, row.protein_group_label, fontsize=9))
    from adjustText import adjust_text

    adjust_text(
        labels,
        ax=ax,
        iter_lim=200,
        arrowprops={"arrowstyle": "-", "color": style.MUTED, "lw": 0.6},
    )
    ax.set_xlabel("Highest nonpermitted adult bulk nTPM (+0.05 for zero display)")
    ax.set_ylabel("Highest reproductive nTPM (+0.05 for zero display)")
    ax.legend(frameon=False, fontsize=9)
    style.footer(
        fig,
        "Color uses the complete evidence cascade, not these two RNA values alone. TBXT can enter through embryonic expression.\nCTAG1A/B and XAGE1A/B are summed within matched source contexts before filtering; named genes have no automatic admission.",
    )
    style.save(fig, figdir, "named-antigen-distribution")

    fig, axes = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={"width_ratios": [1.3, 1]})
    fig.subplots_adjust(left=0.25, right=0.95, top=0.79, bottom=0.29, wspace=0.22)
    style.title(
        fig,
        "Published gene-list provenance and overlap",
        "Membership annotates the quantitative selection; it never bypasses the expression filters",
    )
    ranked = sources.sort_values("mapped_groups", ascending=False).head(15)
    axes[0].barh(
        ranked.source_tag, ranked.mapped_groups, color=style.PALE, label="All mapped groups"
    )
    axes[0].barh(ranked.source_tag, ranked.CTA, color=style.TEAL, label="Final CTA")
    maximum = ranked.mapped_groups.max()
    axes[0].set_xlim(0, maximum * 1.27)
    for i, row in enumerate(ranked.itertuples()):
        axes[0].text(
            row.mapped_groups + maximum * 0.015,
            i,
            f"{row.CTA} / {row.mapped_groups:,}",
            va="center",
            fontsize=8,
        )
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Protein equivalence groups")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0, -0.10), ncol=2)
    order = ranked.source_tag.tolist()
    matrix = (
        overlap.pivot(index="source_a", columns="source_b", values="selected_overlap")
        .reindex(index=order, columns=order)
        .fillna(0)
    )
    im = axes[1].imshow(matrix, cmap="YlGnBu", aspect="auto")
    for i, j in np.ndindex(matrix.shape):
        value = int(matrix.iloc[i, j])
        if value:
            axes[1].text(
                j,
                i,
                str(value),
                ha="center",
                va="center",
                fontsize=6,
                color="white" if value > matrix.to_numpy().max() * 0.55 else style.INK,
            )
    axes[1].set_xticks(range(len(order)), range(1, len(order) + 1))
    axes[1].set_yticks(range(len(order)), range(1, len(order) + 1))
    axes[1].set_title("Shared selected groups\nRows/columns follow list order", fontsize=11)
    fig.colorbar(im, ax=axes[1], shrink=0.65, label="Final CTA groups")
    style.footer(
        fig,
        "Top 15 lists by mapped group count are shown; the downloads include every registered source and pairwise overlap.\nEmbryonic and placental RNA atlases nominate biological association; they are not themselves validated oncofetal-antigen lists.",
    )
    style.save(fig, figdir, "published-list-overlaps")


def build(output: Path, baseline: Path) -> None:
    """Rebuild classifications, known-CTA concordance, figures and provenance."""
    output.mkdir(parents=True, exist_ok=True)
    facts = attach_gene_identities(
        pd.read_csv(baseline / "gene-classifications.csv.gz", low_memory=False)
    )
    ext = external_evidence()
    positive = set(
        ext.loc[
            ext.eligible
            & ext.result.isin(["detected", "mixed"])
            & ext.stage.eq("adult")
            & ext.specimen_state.eq("normal"),
            "gene_id",
        ]
    )
    reviewed = set(
        ext.loc[
            ext.review_status.eq("reviewed")
            & ext.result.isin(["detected", "mixed", "reported_cross_reactivity"]),
            "gene_id",
        ]
    )
    facts["external_eligible_positive"] = facts.gene_id.isin(positive)
    facts["external_review_flag"] = facts.gene_id.isin(reviewed)
    # There are no eligible external negatives in this release; do not silently
    # claim the reconciliation engine has cleared a gene-wide antibody hold.
    if (ext.eligible & ext.result.eq("not_detected")).any():
        raise ValueError(
            "New eligible negatives require context reconciliation before this report is rebuilt"
        )
    gene_result, gene_funnel, _ = classify_layers(facts, decision_layers(facts))
    group_facts = rebuild_groups(gene_result, baseline, output)
    result, funnel, independent = classify_layers(
        group_facts, decision_layers(group_facts), entity_column="protein_equivalence_group_id"
    )
    reverse = classify_layers(
        group_facts,
        list(reversed(decision_layers(group_facts))),
        entity_column="protein_equivalence_group_id",
    )
    if not reverse.classifications.category.equals(result.category):
        raise ValueError("Order-dependent classifications")
    result["policy"] = POLICY_ID
    gene_result["policy"] = POLICY_ID
    gene_columns = [
        "gene_id",
        "symbol",
        "canonical_symbol",
        "historical_ensembl_gene_ids",
        "alternate_ensembl_gene_ids",
        "gene_name_synonyms",
        "identity_catalog_sha256",
        "taxonomy_id",
        "reference",
        "policy",
        "category",
        "exclusion_reasons",
        "Unsafe_reasons",
        "Expression category",
        "Expression detail",
        "Expression forbidden",
        "reproductive_ntpm",
        "adult_bulk_max_ntpm",
        "developmental_mean_tpm",
        "developmental_lineage",
        "reproductive_route",
        "developmental_route",
        "positive_context",
        "oncoref_maximum_tpm",
        "oncoref_maximum_sample",
        "oncoref_cancer_code",
        "tumor_maximum_max_ptpm",
        "tumor_union_ge2",
        "required_coverage",
        "optional_DICE_complete",
        "optional_cornea_complete",
        "Antibody concern",
        "external_review_flag",
        "external_eligible_positive",
        *[c for c in gene_result if c.startswith(("hold:", "unsafe:"))],
    ]
    gene_result = gene_result[gene_columns].copy()
    named = result[result.member_named_controls.ne("")].copy()
    write_csv(result, output / "protein-group-classifications.csv.gz")
    write_csv(gene_result, output / "gene-classifications.csv.gz")
    write_csv(funnel, output / "filtering-funnel.csv")
    write_csv(gene_funnel, output / "gene-filtering-funnel.csv")
    write_csv(independent, output / "independent-layer-effects.csv")
    write_csv(named, output / "named-antigen-decisions.csv")
    for category in ["CTA", "Unsafe", "Neither"]:
        write_csv(
            result[result.category.eq(category)],
            output / f"{category.lower()}-protein-groups.csv.gz",
        )
    members = result[
        ["protein_equivalence_group_id", "member_gene_ids", "category", "exclusion_reasons"]
    ].copy()
    members["gene_id"] = members.member_gene_ids.map(json.loads)
    members = members.explode("gene_id").drop(columns="member_gene_ids")
    member_decisions = gene_result.rename(
        columns={"category": "gene_category", "exclusion_reasons": "gene_exclusion_reasons"}
    ).merge(members, on="gene_id", validate="one_to_one")
    member_decisions = member_decisions.rename(
        columns={
            "category": "protein_group_category",
            "exclusion_reasons": "protein_group_exclusion_reasons",
        }
    )
    write_csv(member_decisions, output / "gene-to-protein-group-decisions.csv.gz")
    write_csv(
        member_decisions[
            member_decisions.gene_category.ne(member_decisions.protein_group_category)
        ],
        output / "aggregation-category-changes.csv",
    )
    nominations = pd.read_csv(
        ROOT / "docs/audits/cta-membership-20261026/nomination-membership.csv.gz", low_memory=False
    )
    nominations = nominations[nominations.protein_candidate_eligible.eq(True)].copy()
    nominations = nominations.merge(
        members,
        on="gene_id" if "gene_id" in nominations else None,
        left_on=None if "gene_id" in nominations else "Ensembl_Gene_ID",
        right_on=None if "gene_id" in nominations else "gene_id",
        how="left",
        validate="many_to_one",
    )
    source_rows = []
    source_sets = {}
    for source, part in nominations.groupby("source_tag"):
        d = part.dropna(subset=["protein_equivalence_group_id"]).drop_duplicates(
            "protein_equivalence_group_id"
        )
        source_sets[source] = set(d.loc[d.category.eq("CTA"), "protein_equivalence_group_id"])
        source_rows.append(
            {
                "source_tag": source,
                "mapped_groups": len(d),
                **{c: int(d.category.eq(c).sum()) for c in ["CTA", "Unsafe", "Neither"]},
            }
        )
    sources = pd.DataFrame(source_rows)
    overlap = pd.DataFrame(
        [
            {"source_a": a, "source_b": b, "selected_overlap": len(av & bv)}
            for a, av in source_sets.items()
            for b, bv in source_sets.items()
        ]
    )
    known = set(nominations.gene_id.dropna())
    known_decisions = member_decisions[
        member_decisions.gene_id.isin(known) | member_decisions.symbol.isin(NAMED)
    ].copy()
    known_decisions["membership_sources"] = known_decisions.gene_id.map(
        nominations.groupby("gene_id").source_tag.agg(lambda s: ";".join(sorted(set(s))))
    )
    known_decisions["named_control"] = known_decisions.symbol.isin(NAMED)
    known_decisions["historical_CTA"] = known_decisions["Expression category"].eq("CTA")
    write_csv(known_decisions, output / "known-cta-decisions.csv.gz")
    write_csv(sources, output / "source-category-counts.csv")
    write_csv(overlap, output / "source-overlaps.csv")
    historical = member_decisions[member_decisions["Expression category"].eq("CTA")]
    write_csv(historical, output / "historical-cta-decisions.csv")
    figures(output, result, funnel, independent, named, sources, overlap)
    counts = result.category.value_counts().to_dict()
    gene_counts = gene_result.category.value_counts().to_dict()
    summary = {
        "policy": POLICY_ID,
        "protein_group_categories": counts,
        "gene_categories_before_grouping": gene_counts,
        "reference_gene_loci": len(gene_result),
        "protein_equivalence_groups": len(result),
        "selected_member_loci": int(member_decisions.protein_group_category.eq("CTA").sum()),
        "historical_CTA_categories": historical.protein_group_category.value_counts().to_dict(),
        "named_control_categories": named.set_index("protein_group_label").category.to_dict(),
        "eligible_external_observations": int(ext.eligible.sum()),
        "order_independent": True,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    simple = named[["protein_group_label", "category", "exclusion_reasons"]].copy()
    simple.exclusion_reasons = (
        simple.exclusion_reasons.str.replace(" | ", "; ", regex=False)
        .fillna("")
        .replace("", "All required gates pass")
    )
    (output / "index.md").write_text(f"""# Current CTA and oncofetal selection funnel

The complete protein-group funnel selects **{counts.get("CTA", 0):,} CTA**, with
**{counts.get("Unsafe", 0):,} Unsafe** and **{counts.get("Neither", 0):,} Neither** among
{len(result):,} protein equivalence groups. These represent {len(gene_result):,}
reference gene loci; {summary["selected_member_loci"]} reference loci belong to selected groups.
Before grouping, the same evidence policy assigns **{gene_counts.get("CTA", 0):,} CTA /
{gene_counts.get("Unsafe", 0):,} Unsafe / {gene_counts.get("Neither", 0):,} Neither** at the gene-locus level.
The two denominators are intentionally reported separately.

## Uniform rules

1. Reproductive HPA RNA ≥5 nTPM and greater than nonpermitted adult bulk, or a
   supported Tyser embryonic lineage mean ≥5 TPM (≥20 cells, ≥5 positive cells,
   ≥10% at TPM≥1). No gene-specific literature admission route.
2. At least one qualifying tumor sample ≥2 native TPM/pTPM. The complete compatible
   OncoRef matrix union is used; same-sample group sums precede maxima.
3. Require bulk/HPA-cluster/HPA-immune/Monaco coverage. DICE and cornea are optional
   for coverage; observed expression still contributes.
4. Apply the linked 2×2 normal-RNA/IHC grid. Above the candidate ceiling is a hold;
   ≥2× the ceiling supplies Unsafe evidence. Preserve exact anatomy and assay units.
5. Retain credible normal-protein and clinical-safety exclusions and externally
   reviewed findings. An antibody concern disables negative-IHC threshold boosts;
   it does not impose an additional gene-wide veto after the stricter RNA screen.
   Context-specific positive findings remain active. The current external collection supplies
   {int(ext.eligible.sum())} fully context-validated observations and no automatic negative overrides.

CTA means the research selection rules pass, Unsafe means a research-policy
exclusion, and Neither means an admission, coverage or unresolved-evidence hold.
Every layer evaluates every entity; later Unsafe evidence can supersede an earlier
Neither. Released annotations and forbidden sequence references are not overwritten.

![Linked thresholds](figures/threshold-grid.png)

## Complete cascade

{markdown_table(funnel)}

![Three-way funnel](figures/classification-funnel.png)
![Candidates removed by each choice](figures/filter-removals.png)

[Sequential counts](filtering-funnel.csv), [independent overlapping effects](independent-layer-effects.csv),
[gene-locus counts](gene-filtering-funnel.csv).

## Named antigens and known CTA tracking

{markdown_table(simple)}

[Current decisions for mapped CTA/source-list genes](known-cta-decisions.csv.gz)
retain each gene's source lists, locus decision, group decision and all exclusion reasons.
The current funnel evaluates the complete reference universe; previous CTA labels
are comparison annotations and do not define its candidate set or admission rules.
Imported lists include broad nomination lists; membership is not a claim that every
entry is an established CTA. Other MAGEs receive ordinary evaluation, not control status.

![Named-antigen distribution](figures/named-antigen-distribution.png)
![Source counts and overlaps](figures/published-list-overlaps.png)

[Source category counts](source-category-counts.csv), [source overlap matrix](source-overlaps.csv).
HPA reproductive and Tyser embryonic measurements are biological nomination routes;
Tyser is one early embryo, not a comprehensive fetal atlas or validated oncofetal target list.
Adult normal-reference coverage does not establish absence in pediatric tissues.

## Identity, aggregation and evidence limits

- Every gene result carries its canonical ENSG, registered historical/alternate IDs,
  synonyms and identity-catalog hash. Ambiguous names and contradictory IDs are
  quarantined rather than assigned to the first match.
- Protein equivalence groups describe identical chosen reference proteins. Their
  IDs hash member identity; they are not sequence hashes or all-isoform equivalence.
  Complete member identity bundles and the pinned group-registry hash are retained.
- CTAG1A/B and XAGE1A/B are summed only within matched samples or tissue contexts.
  Bulk HPA and GTEx components are summed separately before taking a source maximum.
  Aliases of one locus are never added as extra expression.
- DICE lacks donor pairing across genes. Ranked-sum bounds retain this uncertainty;
  a possible upper exceedance causes review, not invented Unsafe evidence. Missing
  optional DICE coverage alone does not hold a group.
- Embryonic lineage means and corneal context means are additive. Detection-cell
  union bounds remain explicit where raw cell unions are unavailable. Embryonic
  admission and strong corneal evidence use the conservative supported bound.
- New external findings retain categorical staining and numeric scores, assay,
  specimen state, reagent lineage, missing donor counts and exact-context limitations.
  CT83's positive tumor-adjacent gastric evidence is flagged; adult disc findings
  remain attached to TBXT. Neither is silently recoded as healthy normal absence.
- TBXT now passes the strict RNA fallback, with the disputed HPA reagent and adult
  disc evidence retained as review flags. This does not claim an experimentally
  validated negative in every adult compartment or clinical safety.
- PRAME has substantial clinical target-specific evidence from IMA203. Its
  expression-screen exclusion must not be read as demonstrated clinical unsafety.
  See the [PRAME, CTAG2 and TBXT evidence assessment](../cta-target-evidence-20261027/index.md).

## Downloads and reproduction

[Selected CTA groups](cta-protein-groups.csv.gz), [Unsafe](unsafe-protein-groups.csv.gz),
[Neither](neither-protein-groups.csv.gz), [all group decisions](protein-group-classifications.csv.gz),
[all locus decisions](gene-classifications.csv.gz), [locus/group crosswalk](gene-to-protein-group-decisions.csv.gz),
[changes caused by aggregation](aggregation-category-changes.csv), [named controls](named-antigen-decisions.csv).

[Validation record](validation.json) records the checks against this report's manifest.

Run the expanded tumor and normal builders, then the expanded funnel, followed by
`scripts/cta_protein_group_normal.py` and `scripts/cta_harmonized_report.py` with this
audit directory. External evidence is rebuilt by `scripts/collect_external_ihc.py`
and `scripts/curate_external_ihc.py`. Source receipts, input hashes and output hashes
are in the [manifest](manifest.json) and component manifests. Every figure also has
an editable SVG alongside its PNG.
""")
    inputs = [
        Path(__file__),
        ROOT / "oncoref/cta_funnel.py",
        ROOT / "oncoref/gene_matching.py",
        ROOT / "oncoref/protein_equivalence.py",
        ROOT / "oncoref/external_ihc.py",
        ROOT / "oncoref/data/external-ihc-evidence.csv.gz",
        ROOT / "scripts/plot_cta_geometric.py",
        ROOT / "docs/audits/cta-membership-20261026/nomination-membership.csv.gz",
        baseline / "manifest.json",
        output / "group-normal/manifest.json",
    ]
    (output / "manifest.json").write_text(
        json.dumps(
            {
                "inputs": input_hashes(inputs),
                "outputs": {
                    str(p.relative_to(output)): sha256(p)
                    for p in [
                        *output.glob("*.csv*"),
                        output / "summary.json",
                        output / "index.md",
                        *(output / "figures").glob("*"),
                    ]
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--baseline", type=Path, default=ROOT / "docs/audits/cta-expanded-funnel-20261027"
    )
    args = parser.parse_args()
    build(args.output, args.baseline)

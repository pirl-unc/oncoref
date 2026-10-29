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
from build_cta_selection_audit import COMPARATORS, OUTPUT, RELEASE, citation, records  # noqa: E402

from oncoref import cta_selection  # noqa: E402


def markdown(frame: pd.DataFrame) -> str:
    """Write a compact Markdown table without an optional tabulate dependency."""

    def cell(value):
        if isinstance(value, float):
            return "—" if np.isnan(value) else f"{value:.4g}"
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


COLORS = {"CTA": "#087E8B", "Unsafe": "#AC493B", "Neither": "#A28C68", "Unassigned": "#B9C3CD"}


def antigen_label(label: str) -> str:
    """Display familiar names together with the grouped gene identities."""
    return {
        "NY-ESO-1": "NY-ESO-1 (CTAG1A/B)",
        "CTAG2": "LAGE-1 (CTAG2)",
        "TBXT": "TBXT (brachyury)",
    }.get(label, label.replace("MAGEA", "MAGE-A"))


def expression_axes(ax, *, x: str, y: str) -> None:
    """Keep zeros visible on native expression axes without display offsets."""
    ax.set_xscale("symlog", linthresh=0.1)
    ax.set_yscale("symlog", linthresh=0.1)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(x)
    ax.set_ylabel(y)


def highlight_antigens(ax, frame: pd.DataFrame, x: str, y: str, *, high_both=False):
    """Label comparators and, optionally, measured groups high on both axes."""
    from adjustText import adjust_text

    mask = frame.protein_group_label.isin(COMPARATORS)
    if high_both:
        mask |= frame[x].ge(5) & frame[y].ge(5)
    labels = []
    for _, row in frame[mask].iterrows():
        ax.scatter(row[x], row[y], s=85, facecolors="none", edgecolors=style.INK, lw=1)
        labels.append(
            ax.text(
                row[x],
                row[y],
                antigen_label(
                    row.protein_group_label
                    if pd.notna(row.protein_group_label)
                    else row.representative_gene_id
                ),
                fontsize=11,
                bbox={"facecolor": style.BG, "edgecolor": "none", "alpha": 0.85, "pad": 0.8},
            )
        )
    adjust_text(
        labels,
        ax=ax,
        iter_lim=300,
        expand=(1.12, 1.25),
        arrowprops={"arrowstyle": "-", "color": style.MUTED, "lw": 0.6},
    )
    return labels


def missing_expression_strip(ax, frame: pd.DataFrame, x: str, *, color=None, norm=None):
    """Show missing y measurements in a separate categorical band, never at zero.

    Vertical staggering is for legibility only and has no expression meaning.
    Every missing-y group with measured x appears, including missing color values.
    """
    band = frame[frame[x].notna()].sort_values([x, "protein_group_label"]).copy()
    band["band_position"] = 0.2 + (np.arange(len(band)) % 5) * 0.11
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#A8B1B8")
    colors = band[color] if color else band.category.fillna("Unassigned").map(COLORS)
    ax.scatter(
        band[x],
        band.band_position,
        c=colors,
        cmap=cmap if color else None,
        norm=norm,
        plotnonfinite=True,
        s=27,
        alpha=0.8,
        edgecolors="white",
        linewidths=0.3,
    )
    ax.set_facecolor("#EAF0F2")
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_title(
        f"Embryonic measurement unavailable · {len(band):,} groups", loc="left", fontsize=11, pad=8
    )
    labels = highlight_antigens(ax, band, x, "band_position")
    return labels


def save_expression_figure(fig, ax, band, out: Path, name: str, note: str, *, colorbar=None):
    """Export a full-width report figure and a landscape slide version."""
    style.footer(fig, note)
    style.save(fig, out, name)
    fig.set_size_inches(16, 6.5)
    for annotation in fig.texts:
        annotation.set_visible(False)
    if band is not None:
        ax.set_position([0.085, 0.43, 0.75 if colorbar else 0.86, 0.53])
        band.set_position([0.085, 0.14, 0.75 if colorbar else 0.86, 0.15])
    else:
        ax.set_position([0.085, 0.16, 0.75 if colorbar else 0.86, 0.80])
    if colorbar is not None:
        colorbar.ax.set_position([0.89, 0.15, 0.018, 0.81])
    ax.tick_params(labelsize=12)
    ax.yaxis.label.set_size(13)
    if band is None:
        ax.xaxis.label.set_size(13)
    else:
        band.xaxis.label.set_size(13)
    # The landscape canvas needs its own label placement after resizing.
    from adjustText import adjust_text
    from matplotlib.patches import FancyArrowPatch

    for patch in list(ax.patches):
        if isinstance(patch, FancyArrowPatch):
            patch.remove()
    adjust_text(
        list(ax.texts),
        ax=ax,
        iter_lim=300,
        expand=(1.15, 1.3),
        arrowprops={"arrowstyle": "-", "color": style.MUTED, "lw": 0.6},
    )
    fig.savefig(out / (name + "-slide.png"), dpi=200, facecolor=style.BG)


def expression_figures(facts: pd.DataFrame, out: Path) -> None:
    """Render full-width expression plots with explicit missing-measurement bands."""
    from matplotlib.colors import SymLogNorm

    embryo = "developmental_all_lineage_max_tpm"
    reproductive = "reproductive_ntpm"
    bulk = "adult_bulk_max_ntpm"
    selected = facts[facts.category.eq("CTA")].copy()
    for column, label, suffix in [
        (reproductive, "Reproductive RNA (nTPM)", "reproductive"),
        (embryo, "Embryonic lineage mean (TPM)", "embryonic"),
    ]:
        fig = plt.figure(figsize=(16, 10))
        has_band = column == embryo
        ax = fig.add_axes([0.095, 0.39 if has_band else 0.20, 0.85, 0.43 if has_band else 0.62])
        style.title(
            fig,
            f"Widely studied antigens: {suffix} expression",
            "Current categories; highlighted antigens follow the same selection rules",
        )
        measured = facts[facts[column].notna() & facts[bulk].notna()]
        for category, color in COLORS.items():
            part = measured[measured.category.fillna("Unassigned").eq(category)]
            ax.scatter(
                part[bulk],
                part[column],
                s=10,
                alpha=0.24,
                color=color,
                label=category,
                rasterized=True,
            )
        expression_axes(ax, x="Adult bulk RNA (nTPM)", y=label)
        highlight_antigens(ax, measured, bulk, column)
        ax.legend(frameon=False, fontsize=11, markerscale=2)
        band = None
        if has_band:
            band = fig.add_axes([0.095, 0.16, 0.85, 0.13], sharex=ax)
            missing_expression_strip(band, facts[facts[column].isna()], bulk)
            ax.set_xlabel("")
            ax.tick_params(labelbottom=False)
            band.set_xlabel("Adult bulk RNA (nTPM)")
        save_expression_figure(
            fig,
            ax,
            band,
            out,
            "widely-studied-antigens-" + suffix,
            "Each dot is a protein-equivalence group. Missing embryonic values are shown separately, not imputed as zero.\nLabels include PRAME, CTAG1A/B, XAGE1A/B, MAGE-A4, CT83, SSX2/B and TBXT where the corresponding measurement exists.",
        )
    for field, unit, suffix in [
        ("healthy_maximum_ntpm", "nTPM", ""),
        ("healthy_maximum_tpm", "TPM", "-immune"),
        ("healthy_maximum_cp10k", "CP10K", "-cornea"),
    ]:
        fig = plt.figure(figsize=(16, 10))
        ax = fig.add_axes([0.095, 0.39, 0.74, 0.43])
        band = fig.add_axes([0.095, 0.16, 0.74, 0.13], sharex=ax)
        style.title(
            fig,
            "Selected CTA: reproductive versus embryonic RNA",
            f"All {len(selected)} selected groups · color: non-disposable healthy maximum ({unit}); gray: color unavailable",
        )
        norm = SymLogNorm(
            linthresh=0.1 if unit != "CP10K" else 0.001,
            vmin=0,
            vmax=max(selected[field].max(), 0.1),
        )
        cmap = plt.get_cmap("viridis").copy()
        cmap.set_bad("#A8B1B8")
        part = selected[selected[reproductive].notna() & selected[embryo].notna()]
        points = ax.scatter(
            part[reproductive],
            part[embryo],
            c=part[field],
            norm=norm,
            cmap=cmap,
            plotnonfinite=True,
            s=40,
            alpha=0.85,
            edgecolors="white",
            linewidths=0.4,
        )
        colorbar = fig.colorbar(
            points, cax=fig.add_axes([0.89, 0.16, 0.018, 0.66]), label=f"Healthy maximum ({unit})"
        )
        ax.axvline(5, ls=":", color=style.MUTED, lw=0.8)
        ax.axhline(5, ls=":", color=style.MUTED, lw=0.8)
        expression_axes(ax, x="", y="Embryonic lineage mean (TPM)")
        ax.set_xlim(0, selected[reproductive].max() * 1.3)
        ax.set_ylim(0, part[embryo].max() * 1.3)
        missing = selected[selected[embryo].isna()]
        missing_expression_strip(band, missing, reproductive, color=field, norm=norm)
        ax.tick_params(labelbottom=False)
        band.set_xlabel("Reproductive RNA (nTPM)")
        highlight_antigens(ax, part, reproductive, embryo)
        save_expression_figure(
            fig,
            ax,
            band,
            out,
            "selected-reproductive-embryonic" + suffix,
            f"{len(part)} measured on both axes; {len(missing)} shown in the missing-embryonic band. Band height is categorical, not expression.\nEmbryonic TPM combines source-separated preimplantation and gastrulation contexts. Missing rows are not zero.\nLabels: widely studied antigens. Companion panels label every selected group ≥5 on both axes. No mixed-unit color maximum.",
            colorbar=colorbar,
        )
    # Label all high-on-both groups in three readable panels; the overview keeps
    # the widely studied labels visible instead of overlaying 82 gene names.
    high = selected[selected[reproductive].ge(5) & selected[embryo].ge(5)].sort_values(reproductive)
    norm = SymLogNorm(linthresh=0.1, vmin=0, vmax=max(selected.healthy_maximum_ntpm.max(), 0.1))
    for number, indices in enumerate(np.array_split(np.arange(len(high)), 3), 1):
        part = high.iloc[indices]
        fig = plt.figure(figsize=(16, 10))
        ax = fig.add_axes([0.095, 0.18, 0.74, 0.63])
        style.title(
            fig,
            "Selected CTA with reproductive and embryonic RNA ≥5",
            f"Labels {number} of 3 · {len(part)} of {len(high)} groups · color: healthy maximum nTPM",
        )
        ax.scatter(high[reproductive], high[embryo], c="#C8D0D6", alpha=0.35, s=28)
        points = ax.scatter(
            part[reproductive],
            part[embryo],
            c=part.healthy_maximum_ntpm,
            cmap="viridis",
            norm=norm,
            s=48,
            edgecolors="white",
            linewidths=0.5,
        )
        colorbar = fig.colorbar(
            points, cax=fig.add_axes([0.89, 0.18, 0.018, 0.63]), label="Healthy maximum (nTPM)"
        )
        expression_axes(ax, x="Reproductive RNA (nTPM)", y="Embryonic lineage mean (TPM)")
        ax.set_xlim(3.5, high[reproductive].max() * 1.5)
        ax.set_ylim(3.5, high[embryo].max() * 1.7)
        highlight_antigens(ax, part, reproductive, embryo, high_both=True)
        save_expression_figure(
            fig,
            ax,
            None,
            out,
            f"selected-high-both-{number}",
            "Colored points and labels cover one third of the selected groups ≥5 on both axes; gray points show the other groups.\nPanels share the same axes. This descriptive threshold does not replace developmental cell-support requirements.",
            colorbar=colorbar,
        )
    selected["embryonic_plot_status"] = np.where(
        selected[embryo].notna(), "measured", "unavailable"
    )
    selected.to_csv(
        OUTPUT / "selected-expression-scatter.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )


def figures(facts: pd.DataFrame, audit) -> None:
    """Render counts and native-unit expression without imputing missing values."""
    style.setup()
    np.random.seed(0)  # deterministic label placement, independent of replay mode
    out = OUTPUT / "figures"
    out.mkdir(exist_ok=True)
    facts.category.fillna("Unassigned").value_counts()
    selected = facts[facts.category.eq("CTA")]
    annotated = audit.memberships.loc[audit.memberships.category.eq("CTA"), "entity_id"].nunique()
    unannotated = len(selected) - annotated
    fig, ax = plt.subplots(figsize=(16, 9))
    fig.subplots_adjust(left=0.29, right=0.77, top=0.83, bottom=0.16)
    style.title(
        fig,
        "Human selection funnel",
        f"{len(facts):,} protein-equivalence groups; every layer evaluates every group",
    )
    left = np.zeros(len(audit.funnel))
    for key, label in [
        ("eligible", "Eligible so far"),
        ("Unsafe", "Unsafe"),
        ("Neither", "Neither"),
        ("Unassigned", "Unassigned"),
    ]:
        ax.barh(
            audit.funnel.stage,
            audit.funnel[key],
            left=left,
            color=COLORS.get(key, COLORS["CTA"]),
            label=label,
        )
        left += audit.funnel[key].to_numpy()
    ax.invert_yaxis()
    ax.set_xlabel("Protein-equivalence groups")
    ax.legend(frameon=False, ncol=4, loc="upper left", bbox_to_anchor=(0, -0.08), fontsize=9)
    ax.text(
        1.02, 1.025, "Eligible / Unsafe / Neither / Unassigned", transform=ax.transAxes, fontsize=9
    )
    for i, row in enumerate(audit.funnel.itertuples()):
        ax.text(
            1.02,
            i,
            f"{row.eligible:,} / {row.Unsafe:,} / {row.Neither:,} / {row.Unassigned:,}",
            transform=ax.get_yaxis_transform(),
            va="center",
            fontsize=9,
        )
    style.footer(
        fig,
        "Neither is assigned by a measured expression profile at the final stage. Earlier holds remain unassigned.\nSequential removals depend on display order; final assignments do not. Unsafe denotes policy exclusion, not demonstrated toxicity.",
    )
    style.save(fig, out, "human-funnel")
    fig, ax = plt.subplots(figsize=(12, 7))
    fig.subplots_adjust(left=0.22, right=0.96, top=0.78, bottom=0.24)
    style.title(fig, "RNA cutoffs", "One ceiling per tissue tier and usable negative-IHC grade")
    grid = np.array([[b * 2**i for i in range(4)] for b in (4, 8, 16)])
    ax.imshow(np.log2(grid), cmap="YlGnBu", aspect="auto")
    ax.set_xticks(
        range(4),
        ["Missing / unusable", "Approved negative", "Supported negative", "Enhanced negative"],
    )
    ax.set_yticks(range(3), ["Vital", "Protected", "Manageable"])
    for i, j in np.ndindex(grid.shape):
        ax.text(
            j,
            i,
            str(grid[i, j]),
            ha="center",
            va="center",
            fontsize=24,
            color="white" if grid[i, j] >= 32 else style.INK,
        )
    style.footer(
        fig,
        "Above ceiling: Vital/Protected supplies Unsafe evidence; Manageable blocks CTA admission. Equality passes.\nBowel epithelium and identified fibroblasts use Manageable; bulk retains its organ tier.\nNative TPM/nTPM numbers are linked. Corneal CP10K uses nominal /100, not a TPM conversion.",
    )
    style.save(fig, out, "active-rna-cutoffs")
    expression_figures(facts, out)
    ranked = audit.sources[audit.sources.primary_source].sort_values(
        "mapped_entities", ascending=False
    )
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
        f"Union: {audit.memberships.entity_id.nunique():,} mapped groups / {len({g for ids in audit.memberships.gene_ids for g in ids}):,} mapped loci; {annotated} selected groups have a registered list annotation.\n{unannotated} selected groups have no annotation in these lists. Counts overlap; list membership is not an eligibility requirement.",
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
        ["Selection", "14 unassigned / unassessed", "No passing/failing count is inferred"],
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


def build() -> None:
    """Publish consistent reports, figure data and slide inputs from one replay."""
    import hashlib

    facts = pd.read_csv(RELEASE / "protein-group-classifications.csv.gz", low_memory=False)
    summary = json.loads((RELEASE / "summary.json").read_text())
    audit, dog = cta_selection.audit(), cta_selection.audit(species="dog")
    figures(facts, audit)
    counts = summary["categories"]
    selected = facts[facts.category.eq("CTA")]
    embryo_measured_count = int(selected.developmental_all_lineage_max_tpm.notna().sum())
    embryo_missing_count = len(selected) - embryo_measured_count
    comparators = facts[facts.protein_group_label.isin(COMPARATORS)][
        [
            "protein_group_label",
            "category",
            "assessment_status",
            "reproductive_ntpm",
            "developmental_all_lineage_max_tpm",
            "adult_bulk_max_ntpm",
            "exclusion_reasons",
        ]
    ].copy()
    comparators["category"] = comparators.category.fillna("Unassigned")
    detailed_reasons = {}
    for label in comparators.protein_group_label:
        evidence = audit.explain(label)
        reasons = []
        for row in evidence[evidence.decision.isin(["unsafe", "hold"])].itertuples():
            if pd.notna(row.value) and pd.notna(row.tissue):
                context = f"{row.tissue}: {row.cell_type}; {row.value} {row.unit}"
                if pd.notna(row.candidate_ceiling):
                    context += f" (ceiling {row.candidate_ceiling:g})"
                reasons.append(context)
            else:
                reasons.append(row.reason)
        detailed_reasons[label] = " | ".join(reasons)
    comparators["exclusion_reasons"] = comparators.protein_group_label.map(detailed_reasons)
    comparators.to_csv(OUTPUT / "widely-studied-antigens.csv", index=False)
    unassigned = facts[facts.category.isna()].copy()
    gap_columns = {
        "normal_panel_missing": ~unassigned.required_coverage,
        "tumor_measurement_missing": ~unassigned.tumor_measured,
        "tumor_group_profile_incomplete": unassigned.tumor_measured
        & ~unassigned.tumor_profile_complete,
        "embryonic_measurement_missing": ~unassigned.developmental_measured,
        "unresolved_evidence": unassigned.unresolved_evidence,
    }
    gaps = unassigned[
        [
            "protein_equivalence_group_id",
            "protein_group_label",
            "member_gene_ids",
            "assessment_status",
        ]
    ].assign(**gap_columns)
    gaps.to_csv(OUTPUT / "unassigned-evidence-gaps.csv", index=False)
    incomplete = gaps[gaps.assessment_status.eq("incomplete_evidence")]
    gap_counts = pd.DataFrame(
        [
            {"Evidence gap": name.replace("_", " "), "Groups": int(incomplete[name].sum())}
            for name in gap_columns
        ]
    )
    normal_gaps = pd.read_csv(OUTPUT / "normal-panel-gap-counts.csv")
    coverage = pd.read_csv(OUTPUT / "tumor-cohort-coverage.csv.gz")
    incomplete_ids = set(incomplete.protein_equivalence_group_id)
    coverage = coverage[coverage.entity_id.isin(incomplete_ids)]
    tumor_source_counts = []
    for provider, rows in coverage.groupby("provider"):
        present = rows.groupby("entity_id").measured_members.max().gt(0)
        tumor_source_counts.append(
            {
                "Source": provider,
                "Cohort views": rows[["cohort", "source_cohort"]].drop_duplicates().shape[0],
                "Groups without any mapped measurement": int((~present).sum()),
            }
        )
    recovery = pd.read_csv(RELEASE / "developmental-identity-recoveries.csv")
    recovery["gene"] = recovery.gene_id.map(
        facts.set_index("representative_gene_id").protein_group_label
    ).fillna(recovery.gene_id)
    recovery.to_csv(OUTPUT / "developmental-identity-recoveries.csv", index=False)
    coverage_text = f"""# Which sources have missing evidence?

These counts describe the **{len(incomplete):,} unassigned groups with incomplete evidence**, not every gene in the reference. Sources overlap: their counts must not be added.

{markdown(gap_counts)}

The earlier **71 tumor-measurement gaps** combined 41 groups without a mapped measurement and 30 measured groups with an incomplete low-expression profile. After the identity repair, **{int(incomplete.tumor_measurement_missing.sum())}** groups lack a mapped measurement, while **{int(incomplete.tumor_group_profile_incomplete.sum())}** have measured tumor values but incomplete group coverage for the low-expression profile. Multi-locus group lower bounds cannot establish low expression. This restriction does not erase observed values and does not block tumor admission when a sample exceeds the floor.

The repair also fixes a source-name sentinel bug: the text `NaN` was interpreted as a historical SCN11A synonym instead of a missing name, contradicting valid ENSG IDs. Seven affected cohort views were reprocessed from their pinned sample matrices and original passing-sample QC: ACINIC, ADCC, CRANIO, DIPG, MENINGIOMA, SARC_MMNST and VSCC. Their corrected measurements replace the older cohort summaries; source samples and units are unchanged. [Repaired cohort evidence and hashes](../cta-selection-20261028/tumor-identity-repair/manifest.json).

## Required normal panels

{markdown(normal_gaps.rename(columns={"source": "Source", "groups": "Groups", "example_genes": "Examples"})[["Source", "Groups", "Examples"]])}

Bulk completeness now requires every nonreproductive tissue concept in the HPA consensus panel, using compatible measurements from HPA consensus, the separate HPA brain atlas, or direct GTEx v11. Reproductive tissues do not block normal-organ coverage. The separate requirements remain 426 HPA tissue/cell contexts, 19 HPA immune populations and 30 Monaco immune populations for every group member. The table distinguishes single-export gaps from gaps remaining after the source union. A zero-valued row is measured coverage. A missing tissue row, unresolved identity or absent group member is not. DICE and adult cornea remain optional and do not contribute to this mandatory-panel gap.

[Every missing group/source/context](normal-panel-coverage-gaps.csv.gz) · [Every member, current/historical IDs, synonyms and observed/missing contexts](normal-panel-member-coverage.csv.gz) · [Source-row identity matches and unresolved alternatives](coverage-identity-matches.csv.gz).

## Cancer and embryonic sources

{markdown(pd.DataFrame(tumor_source_counts))}

[Per-group, per-cohort coverage](tumor-cohort-coverage.csv.gz) identifies each OncoRef cohort view and HPA cancer/cohort, the measured-member count and missing member IDs. Cohort views can share biological samples. [Tumor identity alternatives](tumor-missing-identity-candidates.csv.gz) retains source rows matching the missing groups through any registered current/historical ID or known name candidate, including conflicts that cannot be substituted.

Embryonic coverage now combines **Tyser et al. (2021), E-MTAB-9388** (gastrulation) and **Petropoulos et al. (2016), E-MTAB-3929** (preimplantation days 3–7). Both supply filtered gene-level TPM, with stages and source-specific lineages kept separate. Petropoulos contributes 1,519 measured cells from 88 embryos; the metadata-only ten additional cells are not inserted into the matrix. It recovers **683 of the previous 1,006 Tyser-missing groups**. The current **{int(incomplete.embryonic_measurement_missing.sum()):,}** combined gaps are after reclassification and therefore have a changed denominator. [Coverage and missing member IDs for each atlas and their union](embryonic-coverage-gaps.csv.gz), [Petropoulos source](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-3929/Results), [download hashes and scope](../cta-selection-20261028/alternative-evidence/manifest.json). No later fetal atlas contributes yet. The author raw-count table remains a separate diagnostic, not a TPM substitute.

The new identity pass recovers {int(recovery.status.ne("conflicting_vectors").sum())} gene rows omitted by the earlier exact-ID filter. It keeps {int(recovery.status.eq("conflicting_vectors").sum())} conflicting alternate-row cases unresolved. [Every attempted recovery and source ID](developmental-identity-recoveries.csv). Canonical rows take precedence; multiple alternative rows are accepted only when their entire cell-expression vectors agree. Different values are neither added nor arbitrarily selected.

## Unresolved evidence

The **{int(incomplete.unresolved_evidence.sum())}** current unresolved groups have contributing **HPA v23 normal-tissue IHC rated Uncertain**, with Medium or High staining. These are measured but unresolved observations. [Exact gene, tissue, cell type, staining and reliability](unresolved-protein-contexts.csv.gz).

## Identity matching

The shared matcher checks canonical ENSG IDs, registered historical and alternate IDs, and unique official/synonym names. It preserves the source identifier and matching method. A conflicting explicit ENSG ID is not overwritten with a symbol guess, and ambiguous aliases do not get assigned to the first matching paralog. An empty history list means no replacement is registered in this snapshot, not that a gene has never had another identifier.

For example, **POLR2A / ENSG00000181222** and **CIMIP3 / ENSG00000287363** each have 41 of 50 HPA bulk tissue rows. CIMIP3 already matches the source name **GUCA1ANB** through its stable ID. Both lack amygdala, basal ganglia, cerebellum, hippocampal formation, hypothalamus, midbrain, pituitary gland, spinal cord and vagina rows. Eight missing tissues are Vital (the seven brain/spinal regions plus pituitary); vagina is reproductive. The separate HPA brain file supplies all seven missing brain/spinal regions for both genes. Direct GTEx v11 also supplies POLR2A pituitary (**109.319 TPM**) and the other missing regions, completing its normal bulk coverage. GTEx v11 also contains CIMIP3 as **ENSG00000287363.4 / GUCA1ANB**, with a measured pituitary median of **0 TPM**, completing its normal bulk coverage. The zero is an observation; no GUCA1A measurement is borrowed. [Combined anatomical coverage](combined-rna-coverage.csv.gz). Exact whole-organ GTEx measurements receive the same reviewed bulk IHC grade as their HPA anatomical equivalent; finer brain subregions do not inherit parent-organ negative staining. HPA nTPM and direct GTEx median TPM remain separate measurements with the same nominal cutoffs; they are not pooled or converted. Conversely, a tumor row named NOX5 but carrying the SPESP1-NOX5 readthrough ID remains an explicit ID/name conflict; it is not assigned to either locus by name alone.

This audit distinguishes absent accepted rows, partial tissue coverage and unresolved identity alternatives. An omitted source row can be supplemented by another compatible source; no missing source row is filled with zero. The registry also contains an explicit unresolved protein-group member, ENSG00000275869 in the TP53TG3E/F group; its absence from the canonical catalog remains visible rather than being silently deleted from the group.

[Source hashes and audit code receipt](coverage-audit-receipt.json).
"""
    (OUTPUT / "coverage.md").write_text(coverage_text)
    intake = pd.read_csv(ROOT / "oncoref/data/cta-publication-membership.csv").fillna("")
    registry = pd.read_csv(ROOT / "oncoref/data/cta-publication-sources.csv")
    context = pd.read_csv(ROOT / "oncoref/data/cta-publication-context.csv")
    context.to_csv(OUTPUT / "publication-context.csv", index=False)
    registry["citation"] = registry.source_tag.map(citation)
    registry[
        ["citation", "source_tag", "published_rows", "source_table", "evidence_scope", "source_url"]
    ].to_csv(OUTPUT / "paper-source-scope.csv", index=False)
    targeted, expanded_members = [], []
    source_labels = audit.sources.set_index("source_id").citation
    for tag in ["Bai2016_EGFL6", "Garg2007_SPAG9", "Boel1995_BAGE"]:
        row = intake[intake.source_tag.eq(tag)].iloc[0]
        source = registry[registry.source_tag.eq(tag)].iloc[0]
        match = audit.memberships[audit.memberships.source_id.eq(source.doi)]
        extra = context[context.source_tag.eq(tag)]
        targeted.append(
            {
                "Paper": citation(tag),
                "Principal target": row.source_symbol,
                "Previously recognized CTAs": "; ".join(
                    extra.label.map(antigen_label).str.replace("\n", " ")
                )
                or "None curated",
                "Mapped groups": len(match),
                "Selected CTA": int(match.category.eq("CTA").sum()),
                "Source": "https://doi.org/" + source.doi,
            }
        )
        for target in match.itertuples():
            others = audit.memberships.loc[
                audit.memberships.entity_id.eq(target.entity_id)
                & ~audit.memberships.source_id.eq(source.doi),
                "source_id",
            ]
            expanded_members.append(
                {
                    "Paper": citation(tag),
                    "Target": target.label,
                    "Category": target.category if pd.notna(target.category) else "Unassigned",
                    "Other registered papers": "; ".join(source_labels.loc[others]) or "None",
                    "entity_id": target.entity_id,
                    "member_gene_ids": json.dumps(list(target.gene_ids)),
                }
            )
    targeted = pd.DataFrame(targeted)
    targeted.to_csv(OUTPUT / "targeted-paper-contributions.csv", index=False)
    expanded_members = pd.DataFrame(expanded_members)
    expanded_members.to_csv(OUTPUT / "targeted-paper-overlaps.csv", index=False)
    embryonic_coverage = pd.read_csv(OUTPUT / "embryonic-comparator-coverage.csv")
    routes = pd.DataFrame(
        [
            {
                "route": label,
                "entrants": int(mask.sum()),
                "CTA": int((mask & facts.category.eq("CTA")).sum()),
            }
            for label, mask in [
                (
                    "Reproductive only",
                    facts.reproductive_association & ~facts.developmental_association,
                ),
                (
                    "Developmental only",
                    ~facts.reproductive_association & facts.developmental_association,
                ),
                ("Both", facts.reproductive_association & facts.developmental_association),
            ]
        ]
    )
    inventory = cta_selection.sources()
    inventory.to_csv(OUTPUT / "human-measurement-sources.csv", index=False)
    dog_sources = cta_selection.sources(species="dog")
    dog_sources.to_csv(OUTPUT / "dog-measurement-sources.csv", index=False)
    known = pd.read_csv(RELEASE / "known-cta-decisions.csv.gz", low_memory=False)
    # This is a comparison cohort, never a second set of current CTA assignments.
    [c for c in known if "histor" in c or "Expression category" in c]
    historical = known[
        known.get("historical_expression_category", pd.Series("", index=known.index)).eq("CTA")
    ]
    if historical.empty:
        historical = known[
            known.get("historical_category", pd.Series("", index=known.index)).eq("CTA")
        ]
    historical_counts = (
        historical.protein_group_category.fillna("Unassigned").value_counts().to_dict()
    )
    annotated = audit.memberships.loc[audit.memberships.category.eq("CTA"), "entity_id"].nunique()
    missing = facts[facts.category.isna()].assessment_status.value_counts().to_dict()
    funnel = pd.read_csv(RELEASE / "filtering-funnel.csv")
    overview = f"""# CTA and oncofetal selection

**{counts.get("CTA", 0):,} CTA groups ({summary["selected_member_loci"]:,} loci), {counts.get("Unsafe", 0):,} Unsafe, {counts.get("Neither", 0):,} Neither and {counts.get("Unassigned", 0):,} unassigned**, across 19,863 protein-equivalence groups. Policy: `{summary["policy"]}`. This is the current policy; the [1.21.0 release](../cta-release-1.21.0/index.md) remains a frozen comparison.

[Slides](slides/cta-evidence-walkthrough-expanded.pptx) · [Dog evidence](dog.md) · [Shared API](api.md) · [Policy and IHC matching](policy-review.md) · [All decisions](../cta-selection-20261028/protein-group-classifications.csv.gz) · [Changed decisions](../cta-selection-20261028/classification-changes.csv.gz)

## Categories

**Safety, cancer expression, adult restriction and reproductive/developmental
association are separate requirements. Only safety evidence can establish Unsafe.**
Cancer needs one qualifying measured sample, without a tumor-versus-normal ratio.
Restriction compares a qualifying association signal against supported adult
normal expression in its own native units. [Definitions and logic](../../cta-evidence-axes.md).

![Four separate CTA requirements](../../images/cta-evidence-axes.png)

- **CTA:** association, adult restriction, one qualifying cancer sample, required coverage and acceptable normal expression all pass.
- **Unsafe:** above the single RNA ceiling in Vital/Protected, accepted non-Low protein in those compartments, or retained clinical exclusion evidence. This is a research-policy exclusion, not demonstrated toxicity.
- **Neither:** an affirmative low-expression profile or diffuse-Manageable profile, with complete required RNA coverage and no unresolved evidence. Failed candidate gates alone do not establish Neither.
- **Unassigned:** missing/uncertain evidence or an observed profile that fits none of those definitions. In the API, category is null and assessment_status explains why.

The separate **2× Strong boundary is removed**. Every layer evaluates every group; Neither is assigned only after its complete-profile check. DICE and corneal coverage remain optional, but observed counterevidence contributes.

## Neither profiles

The shared appreciable-expression floor is **2** in native TPM/nTPM, matching the tumor gate. Corneal CP10K uses the existing nominal /100 policy. This is an operational threshold, not a source-validated detection limit.

| Profile | Required evidence | Groups |
|---|---|---:|
| Low expression | All measured normal, reproductive, embryonic and tumor RNA below 2; required normal panels and mapped tumor/developmental measurements present; no detected protein | {summary["neither_profiles"].get("low_expression", 0)} |
| Diffuse Manageable | Fails both reproductive/developmental association routes; RNA ≥2 in at least two distinct Manageable tissue concepts; required normal coverage; no above-cutoff or accepted non-Low Vital/Protected protein evidence | {summary["neither_profiles"].get("diffuse_manageable", 0)} |

Colon and rectum count as one tissue for breadth. A high reproductive measurement cannot conceal broad normal expression. Missing optional measurements are disclosed, not set to zero. {int(summary.get("neither_profiles", {}).get("low_expression", 0))} groups meet the complete low-expression profile. Unassigned statuses: {missing}.

## What the assessment statuses mean

| API status | Meaning | Human groups |
|---|---|---:|
| `screened` | A CTA, Unsafe or Neither category has been assigned by this policy. It does not mean safe, clinically validated or complete coverage of every source. | {int(facts.category.notna().sum()):,} |
| `incomplete_evidence` | No category was assigned, and required normal coverage, mapped tumor/developmental measurements or evidence resolution is missing. | {missing.get("incomplete_evidence", 0):,} |
| `outside_category_definitions` | No category was assigned; the tracked coverage/resolution checks pass but the observed profile meets none of the three definitions. | {missing.get("outside_category_definitions", 0):,} |
| `unassessed` | No species-specific quantitative policy has been applied; used for native dog nominations. | 0 |

An available exclusion can establish Unsafe despite another missing assay. A reproductive-route CTA need not have embryonic measurements. The low-expression Neither profile needs both tumor and embryonic coverage; the diffuse-Manageable profile has different requirements.

Among the {len(incomplete):,} incomplete assessments, gaps overlap:

{markdown(gap_counts)}

[Exact data sources, missing tissue/cell contexts, member identities and cohort coverage](coverage.md) · [Per-group gap flags](unassigned-evidence-gaps.csv). The previous 71 tumor gaps combined 41 missing measurements and 30 measured groups with incomplete low-expression-profile coverage. The table above reflects the subsequent identity repair. The comparator table below reflects the new developmental and normal-bulk sources; MAGE family members are not admission controls.

## Sources of candidate annotations

These published lists annotate the selection; they do not provide a gene-specific admission bypass. The union is **{audit.memberships.entity_id.nunique():,} mapped groups / {len({g for ids in audit.memberships.gene_ids for g in ids}):,} loci**. **{annotated}** selected groups have a listed-paper annotation; **{counts.get("CTA", 0) - annotated}** do not. Source contribution counts overlap. **Mapped groups means curated nominations linked to the current reference, not newly discovered genes or selected CTAs.** The CTA column reports the subset passing the current funnel.

{markdown(audit.sources[audit.sources.primary_source][["citation", "mapped_entities", "CTA", "Unsafe", "Neither", "Unassigned", "source_url"]])}

### Broader nominations and a compact primary-source list

[Sasaroli et al. (2011)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3230481/) contributes its 13 detailed tumor-vascular candidates, all mapped. [Hirbo et al. (2015)](https://doi.org/10.1371/journal.pone.0144155) contributes 86 placenta-enriched nominations: 82 mapped identities and four unresolved names. Its 30 fast-evolving genes are a nested subset, not 30 additional genes. Noncoding mapped rows remain in the intake audit but do not enter protein-group counts. These are vascular and placental nomination scopes, not proof of CTA restriction.

EGFL6 now occurs in Sasaroli, Hirbo and Bai. The main source table and plots omit fully redundant single-digit sets: Bai, Boël and Bradley. Their memberships and citations remain in the downloadable provenance. Garg stays because it uniquely contributes SPAG9 in the current registry. A source is retired from the main view only if the retained sources cover its entire mapped set; the union is unchanged by this presentation choice.

{markdown(audit.sources[~audit.sources.primary_source][["citation", "mapped_entities", "covered_by"]])}

Bai's 538 EGFL6-correlated genes are not a CTA list. Garg recognizes SPA17, SYCP1, NY-ESO-1 and SSX4 alongside SPAG9; Boël recognizes MAGE-A1/A3 alongside the unresolved BAGE label. BAGE's textual overlap with Seager remains outside mapped overlap counts. [Original targeted sets and overlaps](targeted-paper-overlaps.csv) · [Evidence roles](publication-context.csv) · [Every paper's extraction scope](paper-source-scope.csv) · [New source rows, aliases and unresolved names](publication-intake/mapped-source-rows.csv).

The larger papers already contribute their registered gene sets, including Wang's 1,019 CT-pattern genes, Bruggeman's 756 germ-cell genes, da Silva's 1,103 testis-biased and nested 745 CT genes, Chang's 1,036 TGCT genes, Carter's 103 final CT genes, and Loriot's 146 strict plus 134 preferential annotations. Input rows, gene loci and collapsed protein groups differ; these numbers must not be summed. Gong contributes reproductive/placental expression nominations rather than a claimed CTA list. This audit adds the declared vascular and placental sets; it does not claim a fresh exhaustive full-text review of every large source.

EGFL6 is excluded by lung bulk RNA (24.8 nTPM), lung smooth-muscle RNA (363.3 nTPM) and corneal evidence. SPAG9 has accepted bronchial protein staining and broad RNA counterevidence, including spinal cord (66.3 nTPM) and brain oligodendrocytes (287.1 nTPM). These are current quantitative-policy exclusions, irrespective of the original paper's terminology.

![Published source contributions](figures/published-list-contributions.png)
![Published source overlaps](figures/published-list-overlaps.png)

[All pairwise overlaps](human-overlaps.csv.gz) · [Memberships and gene IDs](human-memberships.csv.gz)

## Why roughly 12,000 groups have reproductive or developmental association

The starting association set comes from whole-reference **RNA measurements**, not a 12,000-gene literature list. Reproductive association requires HPA ≥5 nTPM and greater than adult nonreproductive bulk maximum. The developmental route requires a supported source-specific stage/lineage mean ≥5 TPM in Tyser or Petropoulos, ≥20 cells, ≥5 positive cells and ≥10% positive at ≥1 TPM. It establishes embryonic presence, not embryonic restriction. Later fetal and child/adult developmental atlases are not yet integrated.

{markdown(routes)}

[Tyser et al. (2021)](https://doi.org/10.1038/s41586-021-04158-y), [HPA v23 downloads](https://v23.proteinatlas.org/about/download).

## Measurement sources and filtering rules

{markdown(inventory[["name", "evidence_category", "native_unit", "status", "source_url"]])}

{markdown(audit.criteria[["stage", "rule"]])}

![RNA cutoffs](figures/active-rna-cutoffs.png)

The negative-IHC columns require a usable profile control, complete group-member coverage and a supported tissue/cell match. Same-tissue IHC for an unrelated cell cannot increase a ceiling. Values retain their native units; pTPM denotes protein-coding **RNA** TPM. Tumor eligibility uses ≥2 in at least one compatible OncoRef TPM or HPA pTPM sample, without summing maxima from different samples.

## Filtering funnel

{markdown(funnel)}

![Human funnel](figures/human-funnel.png)

[Independent layer effects](../cta-selection-20261028/independent-layer-effects.csv) distinguish overlapping evidence from order-dependent first removals. All counts use protein-equivalence groups; gene-locus counts are explicitly labelled.

## Selected CTA expression

These plots compare reproductive nTPM against the highest measured embryonic lineage mean TPM. **{embryo_measured_count} of {len(selected)} selected groups ({embryo_measured_count / len(selected):.0%}) have both coordinates; {embryo_missing_count} ({embryo_missing_count / len(selected):.0%}) lack an embryonic coordinate.** The plotted embryonic coordinate is the maximum across separate stage/lineage means from Tyser (gastrulation) and Petropoulos (days 3–7). It is not an average across studies. Petropoulos restores measured rows for PRAME, MAGE-A4, CT83, CTAG1A/B and XAGE1A/B. Missing rows remain unavailable; neither atlas establishes expression throughout later fetal development.

The follow-up checked the study authors' **57,490-column raw-count matrix**, separately from the atlas's 25,330-row filtered matrix. The atlas download labelled “Raw counts” is also filtered and retains the same gene rows. Direct author observations are:

{markdown(embryonic_coverage[["source_symbol", "gene_column_present", "cells", "total_reported_reads", "positive_cells"]])}

PRAME, CT83, CTAG1A/B and XAGE1B have zero reported reads across these 1,195 cells; MAGE-A4 has 14 reads in one cell. TBXT has 31,149 reported reads in 226 cells. XAGE1A has no separate column, so complete XAGE1A/B group coverage is not established in this raw-count source; the separate Petropoulos TPM matrix does measure both loci. These are observed nondetections or sparse signals in one early embryo, **not absence throughout fetal development**. The exact per-gene reason for exclusion from the atlas remains unavailable; the authors' counts support low detection but come from a separate quantification. They are not converted to TPM or inserted into the selection funnel. [Author data](https://github.com/ScialdoneLab/human-gastrula-shiny), [count evidence and gene IDs](embryonic-comparator-coverage.csv), [source hashes](embryonic-coverage-receipt.json).

Each figure now shows the {embryo_measured_count} measured groups above and the {embryo_missing_count} missing-embryonic groups in a separate categorical band below. All selected groups are represented; the band has no embryonic-expression value. Gray points mean the healthy-expression color is unavailable. Overview labels highlight the widely studied antigens, including their gene-group names. Three companion panels label every selected group ≥5 on both axes, with the same axes and one third of the names per panel. The latter is a descriptive expression criterion; developmental admission also requires cell support.

![Selected CTA reproductive versus embryonic expression](figures/selected-reproductive-embryonic.png)
![Selected CTA colored by healthy TPM](figures/selected-reproductive-embryonic-immune.png)
![Selected CTA colored by adult corneal CP10K](figures/selected-reproductive-embryonic-cornea.png)

Color uses the maximum supported value across non-disposable healthy compartments within each unit: HPA bulk/brain/cell/immune nTPM, direct GTEx/Monaco/DICE TPM, or corneal CP10K. No maximum mixes units. Reproductive/developmental locations are excluded from the color scale. Missing embryonic measurements are not plotted at zero. [Every selected group's plot data](selected-expression-scatter.csv.gz).

### Groups high in both expression compartments

![Selected groups ≥5 on both axes, labels 1 of 3](figures/selected-high-both-1.png)
![Selected groups ≥5 on both axes, labels 2 of 3](figures/selected-high-both-2.png)
![Selected groups ≥5 on both axes, labels 3 of 3](figures/selected-high-both-3.png)

## Widely studied antigens

{markdown(comparators)}

Each comparison has its own full-width figure. The embryonic comparison retains missing-measurement antigens in a separate band.

![Widely studied antigens: reproductive expression](figures/widely-studied-antigens-reproductive.png)
![Widely studied antigens: embryonic expression](figures/widely-studied-antigens-embryonic.png)

[Detailed tissue/cell/value explanations](human-explanations.csv.gz), [all comparator contexts](../cta-selection-20261028/widely-studied-contexts.csv.gz), and [every previously annotated CTA and its current category](../cta-selection-20261028/known-cta-decisions.csv.gz). Previous CTA labels are comparison annotations only. Historical-cohort counts, where available: {historical_counts}.

## Evidence limits

Adult normal panels do not establish pediatric safety. Tyser contains one gastrulating embryo; Petropoulos adds preimplantation embryos. Later fetal and pediatric coverage remains incomplete. RNA expression does not measure peptide-HLA presentation; HPA negative staining has finite sensitivity. PRAME informed earlier parameter choices and is not independent validation. All gene IDs, equivalence-group members, source versions and native units remain in the exports.
"""
    (OUTPUT / "index.md").write_text(overview)
    (OUTPUT / "policy-review.md").write_text("""# Current policy and IHC matching

The active RNA grid is Vital 4/8/16/32, Protected 8/16/32/64, Manageable 16/32/64/128 across missing, Approved-negative, Supported-negative and Enhanced-negative IHC. There is one boundary; equality passes. Above-ceiling Vital/Protected evidence is Unsafe. Manageable excess blocks CTA but does not automatically establish Neither.

## RNA/IHC correspondence

The matcher records **exact_cell**, **reviewed_compartment**, **same_tissue_unresolved_cell**, or **unmatched**. Its source/version-specific crosswalk is gene-independent. Exact adult tissue/cell matches and declared histological compartment matches can receive negative-IHC credit; unrelated same-tissue cells, stimulated immune cells, fetal samples and other species cannot.

The new compartment mapping relates intestinal epithelial populations to HPA glandular cells, and bronchial basal/club/ciliated cells to respiratory epithelial cells. This is an anatomical interpretation of the histological compartment, not direct proof that every rare/stem cell was observed. The source grade is retained, and the coverage limitation is explicit. Undifferentiated intestinal populations retain a progenitor caution. This assumption affects selection and requires empirical validation; a crosswalk does not improve antibody sensitivity.

[HPA colon histology](https://www.proteinatlas.org/learn/dictionary/normal/colon), [HPA bronchus histology](https://www.proteinatlas.org/learn/dictionary/normal/bronchus), [all matched measurements](../cta-selection-20261028/compartment-ihc-matches.csv.gz).

PRAME colon RNA is 20.3 nTPM with 96 reads. It receives the declared epithelial-compartment match to Enhanced-negative glandular IHC. Its Manageable ceiling is therefore 128. PRAME's adrenal bulk 3.6 nTPM retains the missing-IHC ceiling 4: glandular-cell staining does not supply complete coverage of the specified adrenal medullary and cortical subcompartments. Low reported adrenal protein remains a positive flag, not negative evidence.

CTAG2/LAGE-1 has Enhanced-negative heart cardiomyocyte IHC and bronchial respiratory-epithelial IHC. The revised matcher applies the latter to the basal respiratory compartment; the report retains the parent-compartment coverage assumption. Negative IHC is counterevidence, not proof of clinical safety.

## Anatomy and protein

Identified fibroblasts and large-bowel epithelial cells use Manageable in both RNA and HPA protein assessment. They are not called disposable. Other cell types retain their anatomical tiers; organ-level bulk is not assigned the tier of one constituent cell type. Accepted non-Low Vital/Protected protein can independently supply Unsafe evidence. Manageable protein blocks CTA only. Uncertain staining is unresolved; explicitly Low protein remains flagged under the existing policy.

Selective fibroblast depletion in mouse contexts is not evidence that all human tissue-resident fibroblasts are dispensable. [Kuwabara et al. (2022)](https://elifesciences.org/articles/69854). Epithelial targeting can cause severe colitis: [Parkhurst et al. (2011)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3048186/). The tier choices describe this screening policy's tolerance assumptions, not validated toxicity predictions.

## Curated developmental nominations

The current developmental inputs are Tyser gastrulation and Petropoulos preimplantation expression. Newly identified literature sources remain outside current counts: [Shaw et al. (2024)](https://doi.org/10.1038/s41467-024-47649-y) reports 157 cancer-specific targets with six oncofetal annotations, including exon-specific FN1/TNC evidence; [the TRIM71 liver study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11388079/) contains a liver oncofetal signature. [Cardoso-Moreira et al. (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6658352/) supplies developmental organ trajectories with child/adult comparisons. Their identities, feature scope and developmental/postnatal measurements require ingestion before contributing nominations or exclusion evidence.
""")
    (OUTPUT / "dog.md").write_text(f"""# Dog CTA evidence funnel

**14 native nominations, all unassigned/unassessed.** There is no completed canine quantitative selection or genome-wide denominator. None is automatically called Neither. Human RNA ceilings, HPA grades and clinical evidence do not transfer by orthology.

{markdown(dog.sources[["citation", "mapped_entities", "mapped_gene_loci", "source_url"]])}

Nemec et al. contributes 5 nominations; Hindriks et al. contributes 10; NT5C1B is shared, giving 14. Thirteen have unique gene mappings; one MAGE family assay remains ambiguous. A unique gene mapping does not establish identical protein isoforms.

![Dog funnel](figures/dog-funnel.png)

{markdown(dog.funnel[["stage", "status", "note"]])}

## Sources and gaps

{markdown(dog_sources[["name", "evidence_category", "native_unit", "status", "source_url"]])}

Bgee and BarkBase include 46 overlapping libraries, not independent replication. BarkBase has 37,149 features × 150 TMM-CPM columns; Pfizer MS has 9,929 reported groups × 38 tissue pools and missing matrix cells. Pooled iBAQ is not cell-specific IHC or antigen presentation. Tumor GSE304066 count/FPKM and reference provenance remain unresolved. These quantities must retain their native scales.

[RNA/reference/lineage integration issue #574](https://github.com/pirl-unc/oncoref/issues/574) and [Pfizer MS integration issue #585](https://github.com/pirl-unc/oncoref/issues/585) track the main work. Systematic canine IHC and normal-cell atlases remain coverage gaps. [Pairwise nomination overlaps](dog-overlaps.csv.gz).

The same API exposes sources, memberships, overlaps, criteria, funnel and explanations for human, dog and mouse; unavailable evaluations have null counts and null biological categories.
""")
    (OUTPUT / "api.md").write_text("""# Shared CTA selection API

```python
from oncoref import cta_selection

results = cta_selection.results(species="human")
selected = cta_selection.results(species="human", category="CTA")
unassigned = results[results.category.isna()]
audit = cta_selection.audit(species="human")
audit.sources
audit.overlaps
audit.funnel
audit.criteria
audit.explain("PRAME")
```

The same calls accept `dog` and `mouse`. `category` is CTA, Unsafe, Neither, or null. Null is an assessment status, not a fourth biological class. `assessment_status="screened"` means a category was assigned, not that safety was established or every source was complete. `incomplete_evidence` means an unassigned group lacks required coverage, a tumor/developmental measurement or evidence resolution; `outside_category_definitions` means those tracked checks pass but no category definition fits. `unassessed` denotes native nominations without an applied species policy. [Status definitions and per-group gaps](index.md#what-the-assessment-statuses-mean) are in the report. Coverage counts include an `Unassigned` count so the denominator reconciles.

Human entities are protein-equivalence groups with all member ENSG IDs and curated synonyms/history. Dog entities are native nominations with their actual identity status. Species-specific units, references and coverage remain explicit. `sources()` inventories measurement roles; `audit.sources` counts all published list annotations; `primary_source` selects the compact source view and `covered_by` records retained coverage of retired small sets. Every source table retains its citation/URL.

`audit.explain()` accepts an exact entity ID, gene ID or case-insensitive display label. RNA explanations contain the actual tissue, cell type, value, native unit, tier, single ceiling and IHC-match scope. HPA protein and imported categorical observations remain on their own scales. Unknown queries fail rather than fuzzy-matching a different gene.

Raw data share a version selector across the catalog:

```python
from oncoref import catalog, gtex
catalog.versions("gtex_tissue_rna")             # v11, v10, v8
catalog.ensure("gtex_tissue_rna")               # newest registered: v11
catalog.ensure("gtex_tissue_rna", version="v10")
catalog.provenance("gtex_tissue_rna", "v11")
gtex.tissue_rna("v11")                          # native matrix + provenance
```

[Source versions, historical caches and species-specific data](../../data-versions.md) explains the common contract. This saved analysis pins GTEx v11, HPA normal panels v23, Tyser and Petropoulos snapshots. Public HPA RNA summary loaders default to v25.1; that does not silently migrate this saved analysis. The raw GTEx reader preserves all 68 native contexts; this whole-organ adapter uses 52 and separately lists the 14 microdissected and two cultured/transformed contexts outside its scope.

All results are packaged and checksum-verified. Frozen release audits remain separately accessible; changing the current policy does not rewrite their categories or figures.
""")
    presentation = {
        "summary": summary,
        "report_sha256": hashlib.sha256(
            (RELEASE / "protein-group-classifications.csv.gz").read_bytes()
        ).hexdigest(),
        "papers": records(audit.sources[audit.sources.primary_source]),
        "paper_union_groups": audit.memberships.entity_id.nunique(),
        "paper_union_loci": len({g for ids in audit.memberships.gene_ids for g in ids}),
        "retired_papers": records(audit.sources[~audit.sources.primary_source]),
        "funnel": records(funnel),
        "comparators": records(comparators),
        "routes": records(routes),
        "annotated": annotated,
        "unannotated": counts.get("CTA", 0) - annotated,
        "neither_profiles": summary["neither_profiles"],
        "assessment_statuses": missing,
        "normal_gap_counts": records(normal_gaps),
        "evidence_gap_counts": {name: int(incomplete[name].sum()) for name in gap_columns},
    }
    (OUTPUT / "presentation-data.json").write_text(json.dumps(presentation, indent=2) + "\n")
    print(f"Rendered current report and figures: {counts}")


if __name__ == "__main__":
    build()

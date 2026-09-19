#!/usr/bin/env python3
"""Scientific figures for the focused proteoform CTA analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

KEY = "proteoform_key"
ROOT = Path(__file__).resolve().parents[1]
FOCUSED = [(50, 70), (70, 70), (50, 90), (70, 90)]
TEAL, ORANGE, GRAY, PURPLE = "#127C80", "#C8772B", "#BAC6CF", "#725C9E"
COLLAPSE_TITLE = "Identical-protein gene groups and their effect on CTA selection"


def protein_gene_legend(ax):
    handles = [
        Line2D([], [], marker="o", ls="", color=color, label=label, markersize=6)
        for color, label in [
            (TEAL, "One gene"),
            (ORANGE, "Multiple identical-protein genes"),
        ]
    ]
    ax.legend(
        handles=handles,
        fontsize=7,
        loc="upper right",
        title="Genes per protein entry",
        title_fontsize=8,
    )


def label_coverage_points(
    ax, data, column, *, top_n=10, extra_keys=(), x_column="protein_length_aa"
):
    keys = set(data.nlargest(top_n, column)[KEY]) | set(extra_keys)
    texts = [
        ax.text(row[x_column], row[column] * 100, row.Symbol, fontsize=8)
        for _, row in data[data[KEY].isin(keys)].iterrows()
    ]
    adjust_text(
        texts,
        x=data[x_column].to_numpy(),
        y=data[column].to_numpy() * 100,
        ax=ax,
        expand=(1.2, 1.4),
        iter_lim=300,
        arrowprops={"arrowstyle": "-", "color": "#60717B", "lw": 0.6},
    )


def plot_proteoform_collapse(universe, changes):
    """Explain the individual-gene comparison without implying a different cohort policy."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), layout="constrained")
    merged = universe[universe.n_member_genes.gt(1)].sort_values("n_member_genes")
    axes[0].barh(merged.Symbol, merged.n_member_genes, color=TEAL)
    axes[0].set(
        xlabel="Number of genes encoding the same protein", title="Genes combined into one entry"
    )
    xs = np.arange(4)
    individual, combined = [], []
    for f, p in FOCUSED:
        sub = changes[changes.scenario.eq(f"prevalence_gt{f}_transcriptome_p{p}")]
        individual.append(int(sub.selected_as_any_gene_before.sum()))
        combined.append(int(sub.selected_as_proteoform_now.sum()))
    axes[1].bar(
        xs - 0.18,
        individual,
        width=0.36,
        color=GRAY,
        label="Genes assessed separately (comparison)",
    )
    axes[1].bar(
        xs + 0.18,
        combined,
        width=0.36,
        color=TEAL,
        label="Same-protein genes combined (this report)",
    )
    for x, a, b in zip(xs, individual, combined):
        axes[1].text(x - 0.18, a + 0.5, str(a), ha="center")
        axes[1].text(x + 0.18, b + 0.5, str(b), ha="center")
    axes[1].set(
        xticks=xs,
        xticklabels=[f">{f}%\np{p}" for f, p in FOCUSED],
        ylabel="Selected protein identities",
        title="Effect of combining identical-protein genes",
        ylim=(0, max(individual + combined) * 1.3),
    )
    axes[1].legend(fontsize=8, loc="upper right")
    fig.supxlabel(
        "For each patient, sum expression across genes encoding the same protein (e.g. CTAG1A + CTAG1B).\n"
        "Both methods use the same eligible cohorts (>=20 patient groups); transcriptome cutoffs use the corresponding gene or combined entries.",
        fontsize=9,
    )
    return fig


def heatmap_cells(data, keys, cohort_codes):
    """Keep columns with any observed >10% prevalence, then retain ALL cell values."""
    sub = data[data[KEY].isin(keys) & data.cancer_code.isin(cohort_codes)].copy()
    qualifying = sub.available_sum_complete_measurement & (
        100 * sub.n_expressing > 10 * sub.n_patients
    )
    shown = set(sub.loc[qualifying, "cancer_code"])
    return sub[sub.cancer_code.isin(shown)].copy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "outputs/cta_proteoform_report_20260917")
    args = parser.parse_args()
    out = args.out
    (out / "plots").mkdir(exist_ok=True)
    universe = pd.read_csv(out / "proteoform_universe.csv").set_index(KEY, drop=False)
    cohorts = pd.read_csv(out / "cohort_audit.csv")
    cohorts = cohorts[cohorts.n_patients.ge(20)].copy()
    codes = cohorts.cancer_code.tolist()
    cohort_index = cohorts.set_index("cancer_code")
    groups = pd.read_csv(out / "cohort_overlap_groups.csv").set_index("cancer_code")
    metrics = pd.read_csv(out / "all_proteoform_cohort_metrics.csv.gz")
    metrics = metrics[metrics.cancer_code.isin(codes)].copy()
    summary = pd.read_csv(out / "selection_summary.csv")
    coverage = pd.read_csv(out / "proteoform_coverage_gt10.csv")
    changes = pd.read_csv(out / "gene_vs_proteoform_selection.csv")
    cases = pd.read_csv(out / "case_study_proteoforms.csv")
    patients = pd.read_csv(out / "case_study_patient_expression.csv.gz", dtype={"patient_id": str})
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "savefig.facecolor": "white",
        }
    )
    index = []

    def save(fig, name, title, category="overview"):
        fig.savefig(out / "plots" / f"{name}.png", dpi=165, bbox_inches="tight")
        fig.savefig(out / "plots" / f"{name}.svg", bbox_inches="tight")
        plt.close(fig)
        index.append({"name": name, "title": title, "category": category})

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
    grid = summary.pivot(
        index="prevalence_gt_pct", columns="transcriptome_percentile", values="n_proteoforms"
    ).loc[[50, 70], [70, 90]]
    axes[0].imshow(grid, cmap="GnBu", vmin=0, vmax=grid.to_numpy().max(), aspect="auto")
    for i in range(2):
        for j in range(2):
            n = grid.iloc[i, j]
            axes[0].text(
                j,
                i,
                str(n),
                ha="center",
                va="center",
                fontsize=26,
                color="white" if n > 40 else "#203746",
            )
    axes[0].set(
        xticks=[0, 1],
        xticklabels=["p70", "p90"],
        yticks=[0, 1],
        yticklabels=[">50%", ">70%"],
        title="Selected proteoforms",
        xlabel="Patient transcriptome cutoff",
        ylabel="Required cohort prevalence",
    )
    vals = [293, 300, len(universe), len(cohorts), groups.loc[codes].cancer_type_group.nunique()]
    labels = [
        "Original CTA genes",
        "Identical-protein member loci",
        "Unique CTA proteoforms",
        "Eligible cohort views (n >=20)",
        "Cancer-type groups",
    ]
    axes[1].barh(labels[::-1], vals[::-1], color=[PURPLE, TEAL, TEAL, GRAY, GRAY])
    for i, n in enumerate(vals[::-1]):
        axes[1].text(n + 3, i, str(n), va="center")
    axes[1].set(xlim=(0, 345), title="One focused analysis")
    save(fig, "selection_overview", "Focused proteoform selections and cohort denominator")

    save(plot_proteoform_collapse(universe, changes), "proteoform_collapse", COLLAPSE_TITLE)

    fig, axes = plt.subplots(1, 2, figsize=(14, 9), layout="constrained")
    for ax, p in zip(axes, [70, 90]):
        sub = (
            coverage[coverage.transcriptome_percentile.eq(p)]
            .nlargest(25, "n_cancer_type_groups_gt10")
            .iloc[::-1]
        )
        yy = np.arange(len(sub))
        ax.barh(
            yy,
            sub.fraction_cohort_views_gt10 * 100,
            color=GRAY,
            height=0.75,
            label=f"Cohort views (/{len(cohorts)})",
        )
        ax.barh(
            yy,
            sub.fraction_cancer_type_groups_gt10 * 100,
            color=TEAL,
            height=0.42,
            label=f"Cancer-type groups (/{int(sub.n_eligible_cancer_type_groups.iloc[0])})",
        )
        ax.set(
            yticks=yy,
            yticklabels=sub.Symbol,
            xlim=(0, 105),
            xlabel="Coverage with >10% positive patients (%)",
            title=f"Broad coverage at p{p}",
        )
        ax.legend(fontsize=8, loc="lower right")
    save(fig, "broad_coverage_gt10", "Coverage above 10% patient prevalence across cancer types")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    axes[0].hist(universe.protein_length_aa, bins=np.geomspace(50, 8000, 22), color=TEAL)
    axes[0].set_xscale("log")
    axes[0].set(
        xlabel="Protein length (aa; log scale)",
        ylabel="Proteoform entries",
        title="Lengths of all CTA protein identities",
    )
    broad = coverage[coverage.transcriptome_percentile.eq(90)]
    axes[1].scatter(
        broad.protein_length_aa,
        broad.fraction_cancer_type_groups_gt10 * 100,
        c=np.where(broad.n_member_genes.gt(1), ORANGE, TEAL),
        alpha=0.75,
    )
    axes[1].set_xscale("log")
    axes[1].set(
        xlabel="Protein length (aa)",
        ylabel="Cancer-type groups with >10% p90 positives (%)",
        title="Length versus cancer coverage",
    )
    axes[1].set_ylim(-3, broad.fraction_cancer_type_groups_gt10.max() * 100 + 8)
    protein_gene_legend(axes[1])
    label_coverage_points(axes[1], broad, "fraction_cancer_type_groups_gt10", extra_keys=cases[KEY])
    save(fig, "protein_length_coverage", "Protein lengths and broad p90 coverage")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    ordered = cohorts.sort_values("n_patients", ascending=False)
    axes[0].hist(
        cohorts.n_patients, bins=[20, 50, 100, 250, 500, 1200], color=TEAL, edgecolor="white"
    )
    axes[0].set(
        xlabel="Patient groups per eligible cohort",
        ylabel="Cohort views",
        title="Eligible cohort sizes",
    )
    completeness = (
        metrics[metrics.transcriptome_percentile.eq(90)]
        .groupby("cancer_code")
        .complete_measurement.sum()
    )
    axes[1].hist(completeness, bins=15, color=PURPLE)
    axes[1].set(
        xlabel=f"Number of CTA entries (out of {len(universe)})",
        ylabel="Number of cohorts",
        title="CTA entries with all\nmember-gene data available",
    )
    fig.supxlabel(
        "An entry requires expression values for every member gene in every patient group in the cohort.\n"
        "Recorded zeros count as available data; missing values do not. This measures data availability, not expression prevalence.",
        fontsize=9,
    )
    save(fig, "cohort_sizes_completeness", "Cohort sizes and CTA data availability")
    for page, start in enumerate(range(0, len(ordered), 35), 1):
        sub = ordered.iloc[start : start + 35].iloc[::-1]
        fig, ax = plt.subplots(figsize=(12, 11), layout="constrained")
        ax.barh(
            np.arange(len(sub)),
            sub.n_patients,
            color=np.where(sub.linear_tpm_comparable, TEAL, ORANGE),
        )
        for i, n in enumerate(sub.n_patients):
            ax.text(n * 1.04, i, str(n), va="center", fontsize=8)
        ax.set_xscale("log")
        ax.set(
            yticks=np.arange(len(sub)),
            yticklabels=[f"{r.cancer_code}: {r.cancer_name}" for _, r in sub.iterrows()],
            xlabel="Patient groups (log scale)",
            xlim=(17, 1600),
            title=f"Eligible cohorts ({page}/3): orange = non-TPM scale",
        )
        ax.tick_params(axis="y", labelsize=8)
        save(fig, f"cohort_sizes_{page}", f"All eligible cohort sizes {page}", "appendix")

    all_heatmap_cells = []
    for f, p in FOCUSED:
        scenario = f"prevalence_gt{f}_transcriptome_p{p}"
        selected = pd.read_csv(out / "selections/min20" / scenario / "proteoforms_ranked.csv")
        data = metrics[metrics.transcriptome_percentile.eq(p)]
        chosen = data[
            data[KEY].isin(selected[KEY])
            & data.complete_measurement
            & (100 * data.n_expressing > f * data.n_patients)
        ]
        fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout="constrained")
        top = selected.head(18).iloc[::-1]
        axes[0, 0].barh(top.Symbol, top.n_passing_cancer_type_groups, color=TEAL)
        axes[0, 0].set(
            xlabel="Cancer-type groups passing selection",
            title="A. Breadth at the selection threshold",
        )
        axes[0, 0].tick_params(axis="y", labelsize=8)
        axes[0, 1].scatter(
            selected.protein_length_aa,
            selected[f"p{p}_fraction_cancer_type_groups_gt10"] * 100,
            s=np.sqrt(selected.largest_passing_cohort_size) * 8,
            c=np.where(selected.n_member_genes.gt(1), ORANGE, TEAL),
            alpha=0.7,
        )
        axes[0, 1].set_xscale("log")
        axes[0, 1].set(
            xlabel="Protein length (aa)",
            ylabel="Cancer-type groups with >10% positives (%)",
            title="B. Broad coverage; point size = largest cohort",
        )
        coverage_column = f"p{p}_fraction_cancer_type_groups_gt10"
        axes[0, 1].set_ylim(-3, selected[coverage_column].max() * 100 + 12)
        protein_gene_legend(axes[0, 1])
        label_coverage_points(axes[0, 1], selected, coverage_column, top_n=6)
        linear = data[
            data.linear_tpm_comparable & data.complete_measurement & data.n_expressing.gt(0)
        ]
        axes[1, 0].scatter(
            linear.fraction_expressing * 100,
            linear.mean_expressing_tpm,
            color=GRAY,
            s=7,
            alpha=0.25,
        )
        hits = chosen[chosen.linear_tpm_comparable]
        axes[1, 0].scatter(
            hits.fraction_expressing * 100, hits.mean_expressing_tpm, color=TEAL, s=16, alpha=0.8
        )
        axes[1, 0].axvline(f, color=ORANGE, ls="--")
        axes[1, 0].set_yscale("log")
        axes[1, 0].set(
            xlabel="Positive patient groups (%)",
            ylabel="Mean TPM among positives",
            title="C. TPM means only from linear-TPM sources",
        )
        axes[1, 1].scatter(
            selected.normal_rna_max_member_somatic_ntpm,
            selected[f"p{p}_fraction_cancer_type_groups_gt10"] * 100,
            c=TEAL,
            s=40,
        )
        axes[1, 1].set_xscale("symlog", linthresh=0.1)
        axes[1, 1].set_xlim(left=0)
        axes[1, 1].set(
            xlabel="Maximum single-member normal somatic RNA (nTPM)",
            ylabel="Cancer-type groups with >10% positives (%)",
            title="D. Normal-tissue evidence for review",
        )
        normal_column = "normal_rna_max_member_somatic_ntpm"
        label_coverage_points(
            axes[1, 1],
            selected,
            f"p{p}_fraction_cancer_type_groups_gt10",
            top_n=len(selected) if len(selected) <= 13 else 6,
            extra_keys=selected.nlargest(5, normal_column)[KEY],
            x_column=normal_column,
        )
        fig.suptitle(
            f">{f}% prevalence above p{p}: {len(selected)} proteoforms",
            fontsize=17,
            fontweight="bold",
        )
        fig.supxlabel(
            "Panel D is a member-level HPA annotation, not a summed protein-level safety score. All cohorts have n >=20.",
            fontsize=9,
        )
        save(fig, f"selection_{scenario}", f"Selection process: >{f}% above p{p}", "selection")

        cells = heatmap_cells(data, selected[KEY], codes)
        cells["scenario"] = scenario
        cells["passes_selection"] = cells.complete_measurement & (
            100 * cells.n_expressing > f * cells.n_patients
        )
        all_heatmap_cells.append(cells)
        shown_codes = sorted(
            cells.cancer_code.unique(),
            key=lambda c: (
                groups.loc[c, "cancer_type_group"],
                -cohort_index.loc[c, "n_patients"],
                c,
            ),
        )
        keys = selected[KEY].tolist()
        for ri, rs in enumerate(range(0, len(keys), 24), 1):
            row_keys = keys[rs : rs + 24]
            for ci, cs in enumerate(range(0, len(shown_codes), 26), 1):
                column_codes = shown_codes[cs : cs + 26]
                vals = cells.pivot(
                    index=KEY, columns="cancer_code", values="fraction_expressing"
                ).reindex(index=row_keys, columns=column_codes)
                valid = (
                    cells.pivot(
                        index=KEY,
                        columns="cancer_code",
                        values="available_sum_complete_measurement",
                    )
                    .reindex(index=row_keys, columns=column_codes)
                    .fillna(False)
                )
                complete = (
                    cells.pivot(index=KEY, columns="cancer_code", values="complete_measurement")
                    .reindex(index=row_keys, columns=column_codes)
                    .fillna(False)
                )
                fig, ax = plt.subplots(
                    figsize=(
                        max(9, 0.44 * len(column_codes) + 3),
                        max(4, 0.35 * len(row_keys) + 3),
                    ),
                    layout="constrained",
                )
                cmap = plt.get_cmap("GnBu").copy()
                cmap.set_bad("#d1d1d1")
                cmap.set_under("white")
                im = ax.imshow(
                    vals.where(valid), cmap=cmap, vmin=np.nextafter(0.0, 1.0), vmax=1, aspect="auto"
                )
                for i in range(len(row_keys)):
                    for j in range(len(column_codes)):
                        value = vals.iloc[i, j]
                        if not valid.iloc[i, j] or pd.isna(value):
                            label = "NA"
                        else:
                            label = "<1" if 0 < value < 0.005 else f"{value * 100:.0f}"
                            if not complete.iloc[i, j]:
                                label += "†"
                            elif value > f / 100:
                                label += "*"
                                ax.add_patch(
                                    Rectangle(
                                        (j - 0.48, i - 0.48),
                                        0.96,
                                        0.96,
                                        fill=False,
                                        edgecolor=ORANGE,
                                        lw=0.65,
                                    )
                                )
                        ax.text(
                            j,
                            i,
                            label,
                            ha="center",
                            va="center",
                            fontsize=7,
                            color="white" if valid.iloc[i, j] and value > 0.62 else "#203746",
                        )
                ax.set(
                    xticks=np.arange(len(column_codes)),
                    xticklabels=[
                        f"{c} ({int(cohort_index.loc[c, 'n_patients'])})"
                        + (" [proxy]" if not cohort_index.loc[c, "linear_tpm_comparable"] else "")
                        for c in column_codes
                    ],
                    yticks=np.arange(len(row_keys)),
                    yticklabels=[
                        f"{universe.loc[k, 'Symbol']} | {universe.loc[k, 'protein_length_aa']} aa"
                        for k in row_keys
                    ],
                )
                ax.tick_params(axis="x", rotation=90, labelsize=8)
                ax.tick_params(axis="y", labelsize=8)
                fig.colorbar(
                    im, ax=ax, fraction=0.025, pad=0.01, label=f"Fraction above patient p{p}"
                )
                ax.set_title(
                    f"Selected at >{f}% / p{p}: full frequencies across cancer cohorts\nProteoforms {rs + 1}-{rs + len(row_keys)}/{len(keys)}; cohorts {cs + 1}-{cs + len(column_codes)}/{len(shown_codes)}",
                    fontsize=13,
                )
                fig.supxlabel(
                    "Column rule: any shortlisted proteoform >10%. All frequencies retained, including nonzero values below selection.\nCells are percentages; * = passes selection. † = incomplete member coverage (observed sum); gray/NA = unavailable; white = measured zero.",
                    fontsize=8.5,
                )
                save(
                    fig,
                    f"heatmap_{scenario}_r{ri}_c{ci}",
                    f">{f}% / p{p}: full frequencies, rows {ri}, cohorts {ci}",
                    "heatmaps",
                )
    pd.concat(all_heatmap_cells, ignore_index=True).to_csv(
        out / "shortlist_heatmap_all_displayed_cells.csv.gz", index=False
    )

    # Patient distributions for the five original protein identities, plus any new stringent hits.
    p90 = metrics[metrics.transcriptome_percentile.eq(90)].set_index([KEY, "cancer_code"])
    for _, case in cases.iterrows():
        key, label = case[KEY], case.Symbol
        vals = patients[patients[KEY].eq(key)]
        supported = p90.loc[key]
        passing_codes = (
            supported[
                supported.complete_measurement
                & (100 * supported.n_expressing > 70 * supported.n_patients)
            ]
            .sort_values("n_patients", ascending=False)
            .index.tolist()
        )
        collections = (
            [("passing", passing_codes, "Cohorts passing >70% / p90")] if passing_codes else []
        )
        ordered_codes = sorted(
            codes,
            key=lambda c: (
                groups.loc[c, "cancer_type_group"],
                -cohort_index.loc[c, "n_patients"],
                c,
            ),
        )
        collections += [
            (f"all_{i + 1}", ordered_codes[s : s + 30], f"All eligible cohorts ({i + 1}/4)")
            for i, s in enumerate(range(0, len(codes), 30))
        ]
        for suffix, subcodes, subtitle in collections:
            rng = np.random.default_rng(20260917)
            fig, axes = plt.subplots(
                1,
                2,
                figsize=(14, max(4.5, len(subcodes) * 0.30 + 2.2)),
                layout="constrained",
                sharey=True,
            )
            for i, code in enumerate(subcodes):
                sub = vals[vals.cancer_code.eq(code)]
                yy = i + rng.uniform(-0.18, 0.18, len(sub))
                for hit, color in [(False, GRAY), (True, TEAL)]:
                    mask = sub.positive_p90.eq(hit)
                    axes[0].scatter(
                        sub.loc[mask, "expression_to_p90_ratio"],
                        yy[mask],
                        c=color,
                        s=9,
                        alpha=0.5,
                        edgecolors="none",
                        rasterized=True,
                    )
                    if cohort_index.loc[code, "linear_tpm_comparable"]:
                        axes[1].scatter(
                            sub.loc[mask, "expression"],
                            yy[mask],
                            c=color,
                            s=9,
                            alpha=0.5,
                            edgecolors="none",
                            rasterized=True,
                        )
                positive = sub[sub.positive_p90]
                if len(positive) and cohort_index.loc[code, "linear_tpm_comparable"]:
                    axes[1].scatter([positive.expression.mean()], [i], marker="D", c="black", s=22)
                if not cohort_index.loc[code, "linear_tpm_comparable"]:
                    axes[1].text(
                        0.03,
                        i,
                        "proxy: TPM unavailable",
                        transform=axes[1].get_yaxis_transform(),
                        fontsize=8,
                        color=ORANGE,
                    )
            axes[0].axvline(1, c=ORANGE, ls="--", lw=1)
            axes[0].set(
                yticks=np.arange(len(subcodes)),
                yticklabels=[
                    f"{c} | n={int(cohort_index.loc[c, 'n_patients'])}"
                    + (" †" if not p90.loc[(key, c), "complete_member_coverage"] else "")
                    for c in subcodes
                ],
                xlabel="Summed expression / patient's p90 cutoff",
                title="Every eligible patient group",
            )
            axes[0].invert_yaxis()
            axes[1].set(
                xlabel="Summed clean TPM", title="Linear-TPM sources only; diamond = positive mean"
            )
            for ax in axes:
                ax.set_xscale("symlog", linthresh=0.1)
                ax.set_xlim(left=0)
            axes[0].tick_params(axis="y", labelsize=8)
            fig.suptitle(
                f"{label} | {case.protein_length_aa} aa | {subtitle}",
                fontsize=15,
                fontweight="bold",
            )
            fig.supxlabel(
                "Teal = above p90; gray = below or tied. † = incomplete gene-member coverage; blank rows = unavailable.\nThe transcriptome is collapsed before thresholding. Non-TPM cohorts contribute prevalence, but no TPM mean.",
                fontsize=9,
            )
            safe = label.replace("/", "_").replace(" ", "_")
            save(fig, f"patients_{safe}_{suffix}", f"{label}: {subtitle}", "patients")
    old_index = (
        pd.read_csv(out / "plot_index.csv") if (out / "plot_index.csv").exists() else pd.DataFrame()
    )
    retained = (
        old_index[old_index.category.isin(["mortality", "hpa"])] if len(old_index) else old_index
    )
    pd.concat([pd.DataFrame(index), retained], ignore_index=True).to_csv(
        out / "plot_index.csv", index=False
    )
    print(f"Created {len(index)} figures", flush=True)


if __name__ == "__main__":
    main()

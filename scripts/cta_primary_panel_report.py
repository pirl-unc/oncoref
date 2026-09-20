#!/usr/bin/env python3
"""A summary-first presentation of the threshold-selected CTA shortlist."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from cta_report_common import seal_stage, verify_analysis

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/cta_proteoform_report_20260917"
KEY = "proteoform_key"
SETTINGS = [(50, 70), (70, 70), (50, 90), (70, 90)]
COHORTS = [
    "LUAD",
    "LUSC",
    "COAD",
    "READ",
    "LIHC",
    "BRCA",
    "STAD",
    "PAAD",
    "ESCA",
    "HNSC",
    "PRAD",
    "CESC",
]
NAMES = [
    "Lung adeno.",
    "Lung squamous",
    "Colon",
    "Rectum",
    "Liver",
    "Breast",
    "Stomach",
    "Pancreas",
    "Esophagus",
    "Head/neck",
    "Prostate",
    "Cervix",
]


def prepare(out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from plot_cta_proteoform_report import label_coverage_points, protein_gene_legend

    verify_analysis(out)
    dest = out / "primary_panel"
    dest.mkdir(exist_ok=True)
    lists = {
        (f, p): pd.read_csv(
            out / f"selections/min20/prevalence_gt{f}_transcriptome_p{p}/proteoforms_ranked.csv"
        )
        for f, p in SETTINGS
    }
    primary = lists[50, 90]
    primary_n = len(primary)
    alternative_n = len(lists[50, 70])
    denominator = int(primary.p90_n_eligible_cancer_type_groups.iloc[0]) if len(primary) else 0
    primary.to_csv(dest / "primary_proteoforms.csv", index=False)
    (dest / "primary_proteoforms.txt").write_text("\n".join(primary.Symbol) + "\n")
    metrics = pd.read_csv(out / "all_proteoform_cohort_metrics.csv.gz")
    records = []
    for code in COHORTS:
        vals = pd.read_parquet(out / f"checkpoints/{code}_patients.parquet")
        cuts = (
            pd.read_csv(out / f"checkpoints/{code}_cutoffs.csv")
            .set_index("patient_id")
            .loc[vals.columns]
        )
        for (f, p), selected in lists.items():
            keys = selected[KEY]
            mm = (
                metrics[metrics.cancer_code.eq(code) & metrics.transcriptome_percentile.eq(90)]
                .set_index(KEY)
                .loc[keys]
            )
            assert mm.complete_measurement.all(), code
            values = vals.loc[keys].to_numpy()
            assert np.isfinite(values).all()
            for evaluation in (70, 90):
                positive = (values > cuts[f"p{evaluation}_cutoff"].to_numpy()).any(axis=0)
                records.append(
                    {
                        "cancer_code": code,
                        "selection_prevalence": f,
                        "selection_percentile": p,
                        "n_proteins": len(keys),
                        "evaluation_percentile": evaluation,
                        "n_patients": len(positive),
                        "n_expressing_any": int(positive.sum()),
                        "fraction_expressing_any": float(positive.mean()),
                    }
                )
    coverage = pd.DataFrame(records)
    coverage.to_csv(dest / "same_threshold_panel_comparison.csv", index=False)
    selected_metrics = metrics[
        metrics.cancer_code.isin(COHORTS)
        & metrics.transcriptome_percentile.eq(90)
        & metrics[KEY].isin(primary[KEY])
    ]
    selected_metrics.to_csv(dest / "primary_common_cancer_heatmap.csv", index=False)
    cohort_sizes = coverage.drop_duplicates("cancer_code").set_index("cancer_code").n_patients
    rows90 = (
        coverage[coverage.evaluation_percentile.eq(90)]
        .pivot(
            index="cancer_code",
            columns=["selection_prevalence", "selection_percentile"],
            values="fraction_expressing_any",
        )
        .loc[COHORTS]
    )
    wins = int(rows90[50, 90].gt(rows90[50, 70]).sum())
    manifest = {
        "primary_selection": ">50% prevalence above each patient's p90 in at least one eligible cohort",
        "n_primary_proteins": primary_n,
        "n_benchmark_cohorts": len(COHORTS),
        "n_cohorts_primary_beats_alternative_at_p90": wins,
        "benchmark_cohorts": COHORTS,
        "all_panel_members_measured_in_benchmark": True,
        "comparison_rule": "Patient-level OR at identical p70 or p90 evaluation thresholds; no pooled cohort or subtype maxima",
    }
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    np.random.seed(20260920)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "svg.hashsalt": "oncoref-cta-report",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_under("white")
    cmap.set_bad("#D4D4D4")

    def save(fig, name):
        fig.savefig(dest / f"{name}.png", dpi=165, bbox_inches="tight", metadata={"Date": None})
        fig.savefig(dest / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
        plt.close(fig)

    labels = [f"{name}\n{code} (n={cohort_sizes[code]})" for name, code in zip(NAMES, COHORTS)]
    data = selected_metrics.pivot(
        index=KEY, columns="cancer_code", values="fraction_expressing"
    ).reindex(index=primary[KEY], columns=COHORTS)
    matrix = np.vstack([rows90[50, 90].to_numpy(), data.to_numpy()])
    fig, ax = plt.subplots(figsize=(16, 10), layout="constrained")
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=np.nextafter(0.0, 1.0), vmax=1)
    for i in range(len(matrix)):
        for j, value in enumerate(matrix[i]):
            label = "<1%" if 0 < value < 0.005 else f"{value:.0%}"
            ax.text(
                j,
                i,
                label,
                ha="center",
                va="center",
                fontsize=10,
                color="white" if value > 0.55 else "#203746",
            )
    ax.axhline(0.5, color="#C8772B", lw=3)
    ax.set(
        xticks=range(12),
        xticklabels=labels,
        yticks=range(primary_n + 1),
        yticklabels=[
            f"ANY of the {primary_n} proteins",
            *[f"{r.Symbol} ({r.protein_length_aa} aa)" for r in primary.itertuples()],
        ],
        title=f"Primary {primary_n}-protein panel: expression above each patient's p90",
    )
    ax.tick_params(axis="x", rotation=40)
    plt.setp(ax.get_xticklabels(), ha="right")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="Fraction of patients expressing")
    fig.supxlabel(
        "Broad cohorts for common cancers, without selecting a favorable molecular subtype. All listed proteins have complete member-gene data here.\nANY counts each patient once if one or more listed proteins exceed p90. White = exactly zero; low nonzero frequencies remain colored.",
        fontsize=10,
    )
    save(fig, "primary_common_cancer_heatmap")

    fig, axes = plt.subplots(2, 1, figsize=(16, 9), layout="constrained")
    for ax, eval_p in zip(axes, [70, 90]):
        arr = []
        for f, p in SETTINGS:
            sub = coverage[
                coverage.selection_prevalence.eq(f)
                & coverage.selection_percentile.eq(p)
                & coverage.evaluation_percentile.eq(eval_p)
            ].set_index("cancer_code")
            arr.append(sub.loc[COHORTS].fraction_expressing_any.to_numpy())
        arr = np.array(arr)
        im = ax.imshow(arr, aspect="auto", cmap=cmap, vmin=np.nextafter(0.0, 1.0), vmax=1)
        for i in range(4):
            for j, value in enumerate(arr[i]):
                ax.text(
                    j,
                    i,
                    "<1%" if 0 < value < 0.005 else f"{value:.0%}",
                    ha="center",
                    va="center",
                    color="white" if value > 0.55 else "#203746",
                )
        ax.add_patch(Rectangle((-0.49, 1.51), 11.98, 0.98, fill=False, ec="#C8772B", lw=2.5))
        ax.set(
            xticks=range(12),
            xticklabels=NAMES,
            yticks=range(4),
            yticklabels=[f"{len(lists[f, p])} proteins: >{f}% / p{p}" for f, p in SETTINGS],
            title=f"Same evaluation threshold for every list: above patient p{eval_p}",
        )
        ax.tick_params(axis="x", rotation=25)
    fig.colorbar(
        im, ax=axes, fraction=0.018, pad=0.02, label="Patients expressing ANY listed protein"
    )
    fig.suptitle("Compare list size at the same expression threshold", fontsize=17)
    fig.supxlabel(
        f"Orange box = >50% / p90 panel ({primary_n} proteins). It exceeds the >50% / p70 panel ({alternative_n} proteins)\nat p90 in {wins}/{len(COHORTS)} broad cohorts. Lists are selected independently. White = reported zero.",
        fontsize=10,
    )
    save(fig, "list_size_coverage_comparison")

    fig, ax = plt.subplots(figsize=(12, 7), layout="constrained")
    x = "normal_rna_max_summed_somatic_ntpm"
    y = "p90_fraction_cancer_type_groups_gt10"
    ax.scatter(
        primary[x],
        primary[y] * 100,
        c=np.where(primary.n_member_genes.gt(1), "#C8772B", "#127C80"),
        s=65,
    )
    ax.set_xscale("symlog", linthresh=0.1)
    ax.set(
        xlim=(-0.015, primary[x].max() * 1.6),
        ylim=(-3, 52),
        xlabel="Maximum within-tissue sum across member genes (HPA nTPM)",
        ylabel="Cancer-type groups with >10% p90-positive patients (%)",
        title="Primary panel: cancer coverage and normal-tissue RNA evidence",
    )
    protein_gene_legend(ax)
    label_coverage_points(ax, primary, y, top_n=primary_n, x_column=x)
    fig.supxlabel(
        f"All {primary_n} proteins are labeled; breadth denominator = {denominator} cancer-type groups.\nNormal-tissue RNA sums the same loci as tumor RNA. Zero is rounded below 0.05 nTPM per locus; RNA does not establish protein absence.",
        fontsize=10,
    )
    save(fig, "primary_normal_tissue_scatter")
    seal_stage(
        out,
        "primary",
        [
            Path(__file__),
            out / "analysis_receipt.json",
            *[
                out / f"selections/min20/prevalence_gt{f}_transcriptome_p{p}/proteoforms_ranked.csv"
                for f, p in SETTINGS
            ],
        ],
        [
            *[
                dest / name
                for name in [
                    "primary_proteoforms.csv",
                    "same_threshold_panel_comparison.csv",
                    "primary_common_cancer_heatmap.csv",
                ]
            ],
            *sorted(dest.glob("*.png")),
            dest / "manifest.json",
        ],
    )
    print(
        f"Prepared {primary_n}-protein panel; higher p90 coverage in {wins}/{len(COHORTS)} cohorts"
    )


def render(out):
    from cta_decision_guide import build_guide

    build_guide(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["prepare", "render"])
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    (prepare if args.phase == "prepare" else render)(args.out)

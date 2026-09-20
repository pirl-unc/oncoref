#!/usr/bin/env python3
"""Compare the compact CTA lists across oncoref's leading world-mortality categories."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from cta_report_common import seal_stage, verify_analysis

import oncoref as od

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

os.environ.setdefault("SOURCE_DATE_EPOCH", "0")

ROOT = Path(__file__).resolve().parents[1]
KEY = "proteoform_key"
COMPACT = (50, 70)
TEAL = "#127C80"
ORANGE = "#C8772B"
LABELS = {
    "lung": "Lung",
    "colorectal": "Colorectal",
    "liver": "Liver",
    "breast": "Breast",
    "stomach": "Stomach",
    "pancreas": "Pancreas",
    "esophagus": "Esophagus",
    "head_and_neck": "Head and neck",
    "prostate": "Prostate",
    "cervix": "Cervix",
}
EXTERNAL_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC13343830/"
# Already retrieved during this conversation. These external values are audit-only;
# the ranking and all mortality figures come directly from oncoref's public API.
EXTERNAL_2024_DEATHS = {
    "lung": 1861839,
    "colorectal": 917895,
    "liver": 732489,
    "breast": 693660,
    "stomach": 641554,
    "pancreas": 490786,
    "esophagus": 441994,
    "prostate": 419849,
    "cervix": 279583,
    "leukemia": 292234,
    "head_and_neck": 194108 + 51823 + 49748 + 62060 + 96619,
}
EXTERNAL_TOTAL = 9762507


def panel_patient_coverage(values, cutoffs, fully_measured):
    """OR positive patients across fully measured proteins; preserve missing-panel bounds."""
    valid = np.asarray(fully_measured, dtype=bool)
    n_selected, n_patients = values.shape
    if not valid.any():
        return {
            "n_panel_proteins_available": 0,
            "n_selected_proteins": n_selected,
            "n_patients": n_patients,
            "n_expressing_any": np.nan,
            "fraction_expressing_any": np.nan,
            "coverage_is_lower_bound": True,
            "fraction_expressing_any_upper_bound": 1.0,
        }
    observed = np.any(values[valid] > np.asarray(cutoffs)[None, :], axis=0)
    fraction = float(observed.mean())
    return {
        "n_panel_proteins_available": int(valid.sum()),
        "n_selected_proteins": n_selected,
        "n_patients": n_patients,
        "n_expressing_any": int(observed.sum()),
        "fraction_expressing_any": fraction,
        "coverage_is_lower_bound": not valid.all(),
        "fraction_expressing_any_upper_bound": fraction if valid.all() else 1.0,
    }


def top_mortality_categories(reference):
    """Rank named categories once; a residual bucket is not a specific cancer type."""
    rows = reference[reference.derivation_basis.ne("residual")].copy()
    rows = rows.sort_values(
        ["world_mortality_pct", "burden_category"], ascending=[False, True]
    ).head(10)
    rows.insert(0, "mortality_rank", range(1, len(rows) + 1))
    rows["cancer_label"] = rows.burden_category.map(LABELS)
    assert rows.cancer_label.notna().all()
    return rows


BROAD_COHORTS = {
    "lung": {"LUAD", "LUSC"},
    "colorectal": {"COAD", "READ"},
    "liver": {"LIHC"},
    "breast": {"BRCA"},
    "stomach": {"STAD"},
    "pancreas": {"PAAD"},
    "esophagus": {"ESCA"},
    "head_and_neck": {"HNSC"},
    "prostate": {"PRAD"},
    "cervix": {"CESC"},
}


def report_burden_category(code):
    # Thymoma has no dedicated category in the current burden reference.
    return "other_and_unknown_primary" if code == "THYM" else od.burden_category(code)


def broad_cohort_mask(frame):
    return pd.Series(
        [r.cancer_code in BROAD_COHORTS.get(r.burden_category, set()) for r in frame.itertuples()],
        index=frame.index,
        dtype=bool,
    )


def best_observed(frame, value):
    candidates = frame[frame[value].notna() & broad_cohort_mask(frame)]
    if candidates.empty:
        return None
    return candidates.sort_values(
        [value, "n_patients", "cancer_code"], ascending=[False, False, True]
    ).iloc[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "outputs/cta_proteoform_report_20260917")
    out = parser.parse_args().out.resolve()
    verify_analysis(out)
    dest = out / "mortality"
    dest.mkdir(exist_ok=True)
    reference = od.cancer_burden_df()
    reference.to_csv(dest / "oncoref_mortality_reference.csv", index=False)
    top = top_mortality_categories(reference)
    top.to_csv(dest / "top10_world_mortality.csv", index=False)
    cohorts = pd.read_csv(out / "cohort_audit.csv")
    cohorts["burden_category"] = cohorts.cancer_code.map(report_burden_category)
    assert cohorts.burden_category.notna().all()
    mapping = cohorts[["cancer_code", "cancer_name", "n_patients", "burden_category"]].copy()
    mapping["eligible_n20"] = mapping.n_patients.ge(20)
    mapping["in_top10"] = mapping.burden_category.isin(top.burden_category)
    mapping.to_csv(dest / "cohort_burden_mapping.csv", index=False)
    mapping["broad_histology_match"] = broad_cohort_mask(mapping)
    mapping.to_csv(dest / "cohort_burden_mapping.csv", index=False)
    chosen_cohorts = mapping[
        mapping.eligible_n20 & mapping.in_top10 & mapping.broad_histology_match
    ].copy()
    chosen_cohorts["mortality_rank"] = chosen_cohorts.burden_category.map(
        top.set_index("burden_category").mortality_rank
    )
    chosen_cohorts = chosen_cohorts.sort_values(["mortality_rank", "cancer_code"])
    metrics = pd.read_csv(out / "all_proteoform_cohort_metrics.csv.gz")
    metrics = metrics[metrics.transcriptome_percentile.eq(90)]
    by_pair = metrics.set_index([KEY, "cancer_code"])
    lists = {
        f: pd.read_csv(
            out / f"selections/min20/prevalence_gt{f}_transcriptome_p90/proteoforms_ranked.csv"
        )
        for f in COMPACT
    }
    union_rows, protein_rows = [], []
    for f, selected in lists.items():
        keys = selected[KEY].tolist()
        for cohort in chosen_cohorts.itertuples():
            code = cohort.cancer_code
            vals = pd.read_parquet(out / f"checkpoints/{code}_patients.parquet").loc[keys]
            cutoffs = pd.read_csv(
                out / f"checkpoints/{code}_cutoffs.csv", dtype={"patient_id": str}
            )
            cuts = cutoffs.set_index("patient_id").loc[vals.columns, "p90_cutoff"].to_numpy()
            mm = metrics[metrics.cancer_code.eq(code)].set_index(KEY).loc[keys]
            result = panel_patient_coverage(vals.to_numpy(), cuts, mm.complete_measurement)
            result.update(
                {
                    "prevalence_gt_pct": f,
                    "transcriptome_percentile": 90,
                    "cancer_code": code,
                    "cancer_name": cohort.cancer_name,
                    "burden_category": cohort.burden_category,
                }
            )
            union_rows.append(result)
            for key in keys:
                r = by_pair.loc[key, code]
                # Verify patient-level alignment against the completed original report.
                if r.complete_measurement:
                    assert int((vals.loc[key].to_numpy() > cuts).sum()) == r.n_expressing
                protein_rows.append(
                    {
                        "prevalence_gt_pct": f,
                        KEY: key,
                        "Symbol": selected.set_index(KEY).loc[key, "Symbol"],
                        "cancer_code": code,
                        "burden_category": cohort.burden_category,
                        "n_patients": int(r.n_patients),
                        "complete_measurement": bool(r.complete_measurement),
                        "n_expressing": int(r.n_expressing) if r.complete_measurement else np.nan,
                        "fraction_expressing": r.fraction_expressing
                        if r.complete_measurement
                        else np.nan,
                    }
                )
    union = pd.DataFrame(union_rows)
    proteins = pd.DataFrame(protein_rows)
    union.to_csv(dest / "panel_patient_coverage_by_cohort.csv", index=False)
    proteins.to_csv(dest / "protein_coverage_by_cohort.csv", index=False)
    summary_rows, heatmap_rows = [], []
    for f, selected in lists.items():
        for category in top.itertuples():
            cat = category.burden_category
            sub = proteins[proteins.prevalence_gt_pct.eq(f) & proteins.burden_category.eq(cat)]
            counts = sub.groupby(KEY).fraction_expressing.max()
            u = union[union.prevalence_gt_pct.eq(f) & union.burden_category.eq(cat)]
            best = best_observed(u, "fraction_expressing_any")
            summary_rows.append(
                {
                    "mortality_rank": category.mortality_rank,
                    "burden_category": cat,
                    "cancer_label": category.cancer_label,
                    "world_mortality_pct": category.world_mortality_pct,
                    "prevalence_gt_pct": f,
                    "transcriptome_percentile": 90,
                    "n_selected_proteins": len(selected),
                    "n_matched_cohort_views": len(u),
                    "n_proteins_gt10_in_any_cohort": int(counts.gt(0.1).sum()),
                    "n_proteins_gt50_in_any_cohort": int(counts.gt(0.5).sum()),
                    "n_proteins_gt70_in_any_cohort": int(counts.gt(0.7).sum()),
                    "coverage_scope": "Best observed broad histology cohort; not category-wide or worldwide patient prevalence",
                    "best_panel_fraction": best.fraction_expressing_any
                    if best is not None
                    else np.nan,
                    "best_panel_cohort": best.cancer_code if best is not None else "",
                    "best_panel_cohort_n": best.n_patients if best is not None else np.nan,
                    "best_panel_is_lower_bound": best.coverage_is_lower_bound
                    if best is not None
                    else True,
                    "best_panel_available_proteins": best.n_panel_proteins_available
                    if best is not None
                    else 0,
                }
            )
            for r in selected.itertuples():
                best_p = best_observed(sub[sub[KEY].eq(getattr(r, KEY))], "fraction_expressing")
                heatmap_rows.append(
                    {
                        "prevalence_gt_pct": f,
                        KEY: getattr(r, KEY),
                        "Symbol": r.Symbol,
                        "burden_category": cat,
                        "best_fraction": best_p.fraction_expressing
                        if best_p is not None
                        else np.nan,
                        "best_cohort": best_p.cancer_code if best_p is not None else "",
                        "best_cohort_n": best_p.n_patients if best_p is not None else np.nan,
                    }
                )
    summary = pd.DataFrame(summary_rows)
    best_cells = pd.DataFrame(heatmap_rows)
    summary.to_csv(dest / "compact_list_coverage_summary.csv", index=False)
    best_cells.to_csv(dest / "top10_heatmap_cells.csv", index=False)

    # Audit the internal reference against the external figures already retrieved.
    comparison = []
    for r in top.itertuples():
        external_pct = 100 * EXTERNAL_2024_DEATHS[r.burden_category] / EXTERNAL_TOTAL
        comparison.append(
            {
                "burden_category": r.burden_category,
                "internal_world_mortality_pct": r.world_mortality_pct,
                "internal_source": r.source,
                "internal_world_counts_available": pd.notna(r.world_mortality_count),
                "external_year": 2024,
                "external_deaths": EXTERNAL_2024_DEATHS[r.burden_category],
                "external_world_total_deaths": EXTERNAL_TOTAL,
                "external_world_mortality_pct": external_pct,
                "difference_percentage_points": external_pct - r.world_mortality_pct,
                "comparison_note": "Different years; five-site head/neck approximation; salivary inclusion unresolved"
                if r.burden_category == "head_and_neck"
                else "Different reference years; same broad site",
                "external_source_url": EXTERNAL_URL,
            }
        )
    pd.DataFrame(comparison).to_csv(dest / "mortality_reference_comparison.csv", index=False)

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
    plot_rows = []

    def save(fig, name, title):
        fig.savefig(out / "plots" / f"{name}.png", dpi=165, bbox_inches="tight")
        fig.savefig(out / "plots" / f"{name}.svg", bbox_inches="tight")
        plt.close(fig)
        plot_rows.append({"name": name, "title": title, "category": "mortality"})

    fig, axes = plt.subplots(1, 2, figsize=(15, 8), layout="constrained")
    yy = np.arange(10)
    axes[0].barh(yy, top.world_mortality_pct, color=TEAL)
    for i, value in enumerate(top.world_mortality_pct):
        axes[0].text(value + 0.15, i, f"{value:g}%", va="center")
    axes[0].set(
        yticks=yy,
        yticklabels=[f"{r.mortality_rank}. {r.cancer_label}" for r in top.itertuples()],
        xlabel="Share of global cancer deaths (%)",
        title="Top 10 named oncoref burden categories",
        xlim=(0, 21),
    )
    axes[0].invert_yaxis()
    for f, shift, color in [(50, -0.18, TEAL), (70, 0.18, ORANGE)]:
        data = summary[summary.prevalence_gt_pct.eq(f)].sort_values("mortality_rank")
        axes[1].barh(
            yy + shift,
            data.best_panel_fraction * 100,
            height=0.33,
            color=color,
            label=f"{len(lists[f])} proteins (>{f}% / p90)",
        )
        for i, r in enumerate(data.itertuples()):
            label = f"{r.best_panel_fraction:.0%}{'+' if r.best_panel_is_lower_bound else ''}  {r.best_panel_cohort}"
            axes[1].text(r.best_panel_fraction * 100 + 1, i + shift, label, va="center", fontsize=8)
    axes[1].set(
        yticks=yy,
        yticklabels=top.cancer_label,
        xlim=(0, 137),
        xlabel="Patients expressing at least one listed protein (%)",
        title="Best observed cohort for each compact list",
    )
    axes[1].invert_yaxis()
    axes[1].legend(loc="lower right", fontsize=9)
    fig.supxlabel(
        "Mortality: oncoref, GLOBOCAN 2022 shares. Residual other/unknown excluded.\n"
        "Coverage is the maximum observed among matching broad histology cohorts; it is not worldwide patient prevalence. '+' = lower bound from an incomplete panel.\n"
        "Head/neck uses HNSC and pancreas uses PAAD; narrow subtype results do not stand in for the broad mortality category.",
        fontsize=9,
    )
    save(
        fig,
        "mortality_top10_overview",
        "Top 10 world mortality categories and compact-list patient coverage",
    )

    for f, selected in lists.items():
        keys = selected[KEY].tolist()
        cats = top.burden_category.tolist()
        b = best_cells[best_cells.prevalence_gt_pct.eq(f)]
        matrix = b.pivot(index=KEY, columns="burden_category", values="best_fraction").reindex(
            index=keys, columns=cats
        )
        s = summary[summary.prevalence_gt_pct.eq(f)].set_index("burden_category").loc[cats]
        grid = np.vstack([s.best_panel_fraction.to_numpy(), matrix.to_numpy()])
        fig, ax = plt.subplots(
            figsize=(19, max(6, 0.65 * (len(keys) + 1) + 3.3)), layout="constrained"
        )
        cmap = plt.get_cmap("YlGnBu").copy()
        cmap.set_bad("#D6DADF")
        cmap.set_under("white")
        im = ax.imshow(grid, vmin=np.nextafter(0.0, 1.0), vmax=1, cmap=cmap, aspect="auto")
        for j, cat in enumerate(cats):
            for i in range(len(keys) + 1):
                value = grid[i, j]
                if i == 0:
                    row = s.loc[cat]
                    code = row.best_panel_cohort
                    lower = row.best_panel_is_lower_bound
                else:
                    row = b[b[KEY].eq(keys[i - 1]) & b.burden_category.eq(cat)].iloc[0]
                    code, lower = row.best_cohort, False
                text = "NA" if pd.isna(value) else f"{value:.1%}{'+' if lower else ''}\n{code}"
                ax.text(
                    j,
                    i,
                    text,
                    ha="center",
                    va="center",
                    fontsize=7.5,
                    color="white" if pd.notna(value) and value > 0.55 else "#203746",
                )
                if i and pd.notna(value) and value > f / 100:
                    ax.add_patch(
                        Rectangle((j - 0.49, i - 0.49), 0.98, 0.98, fill=False, ec=ORANGE, lw=1.7)
                    )
        ax.axhline(0.5, color="white", lw=4)
        ax.set(
            xticks=range(10),
            xticklabels=[
                f"{r.mortality_rank}. {r.cancer_label}\n{r.world_mortality_pct:g}% of deaths"
                for r in top.itertuples()
            ],
            yticks=range(len(keys) + 1),
            yticklabels=["ANY listed protein"]
            + [f"{r.Symbol} ({r.protein_length_aa} aa)" for r in selected.itertuples()],
            title=f"{len(selected)}-protein list (>{f}% / p90): best observed prevalence by mortality category",
        )
        ax.tick_params(axis="x", rotation=35)
        plt.setp(ax.get_xticklabels(), ha="right")
        fig.colorbar(
            im, ax=ax, fraction=0.022, pad=0.015, label="Fraction above each patient's p90"
        )
        fig.supxlabel(
            "Each cell names its supporting cohort. All frequencies are shown, including zero; NA means no fully measured entry.\n"
            "ANY is a patient-level OR within one cohort, then the best cohort is selected. '+' marks an incomplete-panel lower bound.\n"
            "Orange border: that protein passes the list's prevalence cutoff in the named cohort. Maxima use broad histology cohorts; cohorts are never added or pooled.\n"
            "Head/neck uses HNSC; pancreas uses PAAD. Thymoma and neuroendocrine pancreas are excluded from these site estimates.",
            fontsize=9,
        )
        save(
            fig,
            f"mortality_top10_gt{f}_p90",
            f"Top mortality heatmap: {len(selected)}-protein p90 list",
        )

        # Every eligible matching cohort remains visible, even at zero prevalence.
        codes = chosen_cohorts.cancer_code.tolist()
        for page, start in enumerate(range(0, len(codes), 20), 1):
            page_codes = codes[start : start + 20]
            pp = proteins[proteins.prevalence_gt_pct.eq(f)]
            mat = pp.pivot(index=KEY, columns="cancer_code", values="fraction_expressing").reindex(
                index=keys, columns=page_codes
            )
            uu = union[union.prevalence_gt_pct.eq(f)].set_index("cancer_code").loc[page_codes]
            grid = np.vstack([uu.fraction_expressing_any.to_numpy(), mat.to_numpy()])
            fig, ax = plt.subplots(
                figsize=(18, max(6, 0.42 * (len(keys) + 1) + 3.5)), layout="constrained"
            )
            im = ax.imshow(grid, vmin=np.nextafter(0.0, 1.0), vmax=1, cmap=cmap, aspect="auto")
            for i in range(grid.shape[0]):
                for j in range(grid.shape[1]):
                    value = grid[i, j]
                    txt = (
                        "NA" if pd.isna(value) else ("<1%" if 0 < value < 0.005 else f"{value:.0%}")
                    )
                    if i == 0 and uu.iloc[j].coverage_is_lower_bound:
                        txt += "+"
                    ax.text(
                        j,
                        i,
                        txt,
                        ha="center",
                        va="center",
                        fontsize=8,
                        color="white" if pd.notna(value) and value > 0.55 else "#203746",
                    )
            cc = chosen_cohorts.set_index("cancer_code", drop=False).loc[page_codes]
            ax.set(
                xticks=range(len(page_codes)),
                xticklabels=[
                    f"{r.cancer_code} (n={r.n_patients})\n{LABELS[r.burden_category]}"
                    for r in cc.itertuples()
                ],
                yticks=range(len(keys) + 1),
                yticklabels=["ANY listed protein", *selected.Symbol.tolist()],
                title=f"{len(selected)}-protein list (>{f}% / p90): individual cohorts, page {page}",
            )
            ax.axhline(0.5, color="white", lw=4)
            ax.tick_params(axis="x", rotation=90)
            fig.colorbar(im, ax=ax, fraction=0.02, pad=0.015, label="Fraction above patient p90")
            fig.supxlabel(
                "Each column is a separate cohort view with n >=20 patient groups; overlapping parent/subtype views are not pooled.\n"
                "NA = incomplete member-gene data. '+' = ANY coverage is a lower bound using the fully measured listed proteins. All low and zero frequencies remain visible.",
                fontsize=9,
            )
            save(
                fig,
                f"mortality_cohorts_gt{f}_p90_{page}",
                f"Mortality cohorts: {len(selected)}-protein list, page {page}",
            )

    index = pd.read_csv(out / "plot_index.csv")
    index = pd.concat(
        [index[index.category.ne("mortality")], pd.DataFrame(plot_rows)], ignore_index=True
    )
    index.to_csv(out / "plot_index.csv", index=False)
    lines = [
        "# Compact CTA lists across oncoref's top 10 world-mortality categories",
        "",
        "Mortality ranking uses oncoref.cancer_burden_df() and world_mortality_pct (curated GLOBOCAN 2022 shares).",
        "Cohorts are mapped with oncoref.burden_category(). The residual other/unknown bucket is excluded from named-cancer ranking.",
        "No absolute death counts are inferred: oncoref's count/total fields are blank.",
        "",
        "Coverage uses p90 within each patient's collapsed transcriptome and n >=20 patient groups. Each site summary is the maximum observed over matched cohorts, restricted to broad histologies in BROAD_COHORTS; the supporting cohort is named. Overlapping views are never added or pooled. This is not population-weighted or worldwide patient coverage.",
        "The ANY result counts each patient once if at least one fully measured shortlisted protein passes p90. Missing panel members yield a lower bound (+), not assumed negatives. Per-protein NA is distinct from zero.",
        "",
        f"| Rank | Cancer category | Global mortality share | {len(lists[50])}-protein ANY coverage, broad cohort | {len(lists[70])}-protein ANY coverage, broad cohort |",
        "|---:|---|---:|---|---|",
    ]
    for r in top.itertuples():
        cells = []
        for f in COMPACT:
            q = summary[
                summary.prevalence_gt_pct.eq(f) & summary.burden_category.eq(r.burden_category)
            ].iloc[0]
            cells.append(
                f"{q.best_panel_fraction:.1%}{'+' if q.best_panel_is_lower_bound else ''}; {q.best_panel_cohort}, n={int(q.best_panel_cohort_n)}"
            )
        lines.append(
            f"| {r.mortality_rank} | {r.cancer_label} | {r.world_mortality_pct:g}% | {' | '.join(cells)} |"
        )
    lines += [
        "",
        "## Compatibility with the previously retrieved external figures",
        "",
        "The external IARC table is GLOBOCAN 2024, published 8 July 2026; it has the same top-five order but different shares and site grouping. It is an audit comparator only, not the source of this report's ranking.",
        "Oncoref combines oral cavity, pharynx, larynx and nasopharynx as head/neck, whereas the external table separates these. Oncoref splits AML and other leukemia; the external table combines leukemia. Therefore the top-ten lists are not directly interchangeable.",
        "For comparability, the external head/neck audit sums lip/oral cavity, oropharynx, hypopharynx, nasopharynx and larynx. Differences between years are documented, not treated as proven reference errors.",
        "Oncoref's liver burden includes intrahepatic bile duct, but its CHOL cohort currently maps to gallbladder_biliary. The report follows that mapping and flags the scope mismatch; no reference mappings were silently changed.",
        "Head/neck coverage uses HNSC. Thymoma and salivary cohorts are not used as proxies for that burden category. Pancreas coverage uses PAAD, excluding pancreatic neuroendocrine tumors. Each site estimate remains an observed cohort result, not global prevalence.",
        "Reference follow-ups: https://github.com/pirl-unc/oncoref/issues/542 (mortality refresh and counts) and https://github.com/pirl-unc/oncoref/issues/543 (CHOL/salivary category scope).",
        f"External audit source: [Sung et al., Table 1, DOI 10.3322/caac.70090]({EXTERNAL_URL}).",
        "",
        "## Files",
        "",
        "- top10_world_mortality.csv: internal ranking and provenance.",
        "- compact_list_coverage_summary.csv: per-list coverage and supporting cohorts.",
        "- top10_heatmap_cells.csv: every category/protein cell.",
        "- protein_coverage_by_cohort.csv and panel_patient_coverage_by_cohort.csv: complete cohort-level results.",
        "- mortality_reference_comparison.csv: internal shares alongside the already retrieved 2024 reference.",
        "- cohort_burden_mapping.csv: all included and excluded cohort mappings.",
    ]
    (dest / "report.md").write_text("\n".join(lines) + "\n")
    (dest / "report.txt").write_text("\n".join(line.lstrip("# ") for line in lines) + "\n")
    manifest = {
        "oncoref_version": od.__version__,
        "ranking_api": "oncoref.cancer_burden_df()",
        "mapping_api": "oncoref.burden_category()",
        "ranking_source_year": 2022,
        "external_audit_year": 2024,
        "external_audit_source": EXTERNAL_URL,
        "residual_category_excluded": "other_and_unknown_primary",
        "n_matched_cohort_views": len(chosen_cohorts),
        "n_figures": len(plot_rows),
        "n_patient_panel_checks": len(union),
        "internal_reference_sha256": hashlib.sha256(
            (ROOT / "oncoref/data/cancer-incidence-mortality.csv").read_bytes()
        ).hexdigest(),
        "coverage_definition": "Maximum observed cohort prevalence; ANY uses patient-level OR within each cohort; incomplete panel gives lower bound",
    }
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    # A compact list is nested in the larger one, so observed union coverage cannot increase.
    pivot = union.pivot(
        index="cancer_code", columns="prevalence_gt_pct", values="fraction_expressing_any"
    )
    assert (pivot[70] <= pivot[50] + 1e-12).all()
    assert union.n_expressing_any.le(union.n_patients).all()
    seal_stage(
        out,
        "mortality",
        [Path(__file__), out / "analysis_receipt.json"],
        sorted(dest.glob("*.csv")),
    )
    print(
        summary[
            ["cancer_label", "prevalence_gt_pct", "best_panel_fraction", "best_panel_cohort"]
        ].to_string(index=False)
    )
    print(
        f"Created mortality report with {len(chosen_cohorts)} cohort views and {len(plot_rows)} figures"
    )


if __name__ == "__main__":
    main()

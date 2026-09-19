#!/usr/bin/env python3
"""Extend the one-off CTA report with cohort-size, overlap and patient views.

Uses the saved, QC-filtered analytical snapshot. Minimum-size selection is a
separate, explicit axis: the default is >=10 patient groups, while the requested
larger-cohort analysis uses >=20 and linear-TPM-comparable sources only.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from oncoref.expression import per_sample_expression
from oncoref.source_matrices import CACHE_DIR_ENV_VAR, source_sample_namespace

ROOT = Path(__file__).resolve().parents[1]
ID = "Ensembl_Gene_ID"
PERCENTILES = (70, 90)
PREVALENCES = (50, 70)
SIZE_CUTS = (1, 5, 10, 20, 25, 50, 100, 250, 500)


def csvout(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, float_format="%.12g")


def passes(frame, prevalence, minimum=1, linear_only=False):
    mask = frame.complete_measurement & (100 * frame.n_expressing > prevalence * frame.n_patients)
    mask &= frame.n_patients >= minimum
    if linear_only:
        mask &= frame.linear_tpm_comparable
    return mask


def identity_key(row):
    patient = str(row.patient_id)
    if patient.startswith(("TCGA-", "TARGET-")):
        namespace = patient.split("-")[0]
    elif patient.startswith(("THR", "TH")) and row.patient_id_basis == "Treehouse_donor_prefix":
        namespace = "TREEHOUSE"
    else:
        namespace = source_sample_namespace(str(row.source_cohort))
    return f"{namespace}|{patient}"


def overlap_groups(audit, cohorts):
    """Connected components of known patient overlap across ALL cohort views.

    A broad cohort connects its subtype views even when the subtype views do not
    overlap each other. Counting one component is a deliberately conservative
    breadth metric; it is not a count of independent studies.
    """
    included = audit[audit.included_in_report].copy()
    included["identity_key"] = included.apply(identity_key, axis=1)
    sets = included.groupby("cancer_code").identity_key.agg(set).to_dict()
    codes = sorted(sets)
    parent = {c: c for c in codes}

    def root(c):
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    edges = []
    for left, right in itertools.combinations(codes, 2):
        shared = len(sets[left] & sets[right])
        if shared:
            parent[root(right)] = root(left)
            edges.append(
                {
                    "cohort_a": left,
                    "cohort_b": right,
                    "n_shared_patient_keys": shared,
                    "fraction_of_smaller_cohort": shared / min(len(sets[left]), len(sets[right])),
                }
            )
    components = {}
    for c in codes:
        components.setdefault(root(c), []).append(c)
    sizes = cohorts.set_index("cancer_code").n_patients.to_dict()
    rows = []
    for members in components.values():
        anchor = sorted(members, key=lambda c: (-sizes[c], c))[0]
        for c in members:
            rows.append(
                {
                    "cancer_code": c,
                    "overlap_group": anchor,
                    "group_members": ";".join(sorted(members)),
                    "n_cohort_views_in_group": len(members),
                    "n_unique_known_patient_keys_in_group": len(
                        set.union(*(sets[x] for x in members))
                    ),
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(
        edges,
        columns=["cohort_a", "cohort_b", "n_shared_patient_keys", "fraction_of_smaller_cohort"],
    )


def summarize_genes(selected, groups, min_size=10, large_min_size=20):
    columns = [
        ID,
        "Symbol",
        "n_passing_cohort_views",
        "n_overlap_adjusted_groups",
        "passing_cohorts",
        "passing_cancer_types",
        "passing_overlap_groups",
        "largest_passing_cohort_size",
        "largest_passing_cohorts",
        "largest_passing_cancer_types",
        "largest_linear_passing_cohort_size",
        "passes_min_cohort_size",
        "passes_large_linear_filter",
        "largest_cohort_n_expressing",
        "largest_cohort_fraction_expressing",
        "largest_cohort_mean_expressing",
        "largest_cohort_expression_unit",
    ]
    if selected.empty:
        return pd.DataFrame(columns=["breadth_rank", *columns])
    group_map = groups.set_index("cancer_code").overlap_group.to_dict()
    rows = []
    for gid, sub in selected.groupby(ID, sort=False):
        sub = sub.sort_values(["n_patients", "cancer_code"], ascending=[False, True])
        best = sub.iloc[0]
        largest = sub[sub.n_patients.eq(best.n_patients)]
        linear = sub[sub.linear_tpm_comparable]
        blocks = sorted({group_map[c] for c in sub.cancer_code})
        rows.append(
            {
                ID: gid,
                "Symbol": best.Symbol,
                "n_passing_cohort_views": len(sub),
                "n_overlap_adjusted_groups": len(blocks),
                "passing_cohorts": ";".join(sub.cancer_code),
                "passing_cancer_types": ";".join(sub.cancer_name),
                "passing_overlap_groups": ";".join(blocks),
                "largest_passing_cohort_size": int(best.n_patients),
                "largest_passing_cohorts": ";".join(largest.cancer_code),
                "largest_passing_cancer_types": ";".join(largest.cancer_name),
                "largest_linear_passing_cohort_size": int(linear.n_patients.max())
                if len(linear)
                else np.nan,
                "passes_min_cohort_size": bool(best.n_patients >= min_size),
                "passes_large_linear_filter": bool((linear.n_patients >= large_min_size).any()),
                "largest_cohort_n_expressing": int(best.n_expressing),
                "largest_cohort_fraction_expressing": best.fraction_expressing,
                "largest_cohort_mean_expressing": best.mean_expression_expressing,
                "largest_cohort_expression_unit": best.expression_unit,
            }
        )
    result = (
        pd.DataFrame(rows, columns=columns)
        .sort_values(
            ["n_overlap_adjusted_groups", "largest_passing_cohort_size", "Symbol"],
            ascending=[False, False, True],
        )
        .reset_index(drop=True)
    )
    # Equal breadth scores share a rank; largest n only orders ties for display.
    result.insert(
        0,
        "breadth_rank",
        result.n_overlap_adjusted_groups.rank(method="dense", ascending=False).astype(int),
    )
    return result


def extract_patient_values(out, metrics, audit, cutoffs, reuse=False):
    dest = out / "six_genes_patient_values.parquet"
    if reuse and dest.exists():
        return pd.read_parquet(dest)
    strict = metrics[(metrics.transcriptome_percentile == 90) & passes(metrics, 75)]
    genes = strict[[ID, "Symbol"]].drop_duplicates().sort_values("Symbol")
    pieces = []
    included = audit[audit.included_in_report]
    for i, (code, sub) in enumerate(included.groupby("cancer_code"), 1):
        print(f"Patient values [{i}/142] {code}", flush=True)
        matrix = per_sample_expression(code, normalize="tpm_clean", auto_fetch=False)
        by_id = matrix.set_index(ID)
        by_id = by_id.reindex(genes[ID])
        mapping = dict(zip(sub.sample_id, sub.patient_id))
        panel = by_id[list(mapping)].rename(columns=mapping).T.groupby(level=0).mean().T
        for (_, gene), values in zip(genes.iterrows(), panel.to_numpy()):
            part = pd.DataFrame(
                {
                    "cancer_code": code,
                    ID: gene[ID],
                    "Symbol": gene.Symbol,
                    "patient_id": panel.columns,
                    "expression": values,
                }
            )
            part = part.merge(
                cutoffs[cutoffs.cancer_code.eq(code)],
                on=["cancer_code", "patient_id"],
                validate="many_to_one",
            )
            part["positive_p90"] = (
                part.expression.gt(part.p90_expression_cutoff) & part.expression.notna()
            )
            part["expression_to_p90_ratio"] = np.divide(
                part.expression,
                part.p90_expression_cutoff,
                out=np.full(len(part), np.nan),
                where=part.p90_expression_cutoff > 0,
            )
            pieces.append(part)
    values = pd.concat(pieces, ignore_index=True)
    values.to_parquet(dest, index=False)
    csvout(values, out / "six_genes_patient_values.csv.gz")
    return values


def verify_patient_values(values, metrics):
    checked = 0
    for (gene, code), sub in values.groupby([ID, "cancer_code"]):
        truth = metrics[
            (metrics[ID] == gene)
            & metrics.cancer_code.eq(code)
            & metrics.transcriptome_percentile.eq(90)
        ].iloc[0]
        assert len(sub) == truth.n_patients
        if truth.complete_measurement:
            positives = sub[sub.positive_p90]
            assert len(positives) == truth.n_expressing
            if len(positives):
                assert np.isclose(
                    positives.expression.mean(), truth.mean_expression_expressing, rtol=1e-9
                )
            else:
                assert pd.isna(truth.mean_expression_expressing)
            checked += 1
    return checked


def build_tables(out, metrics, cohorts, audit, min_size, large_min_size):
    if (out / "selections").exists():
        shutil.rmtree(out / "selections")
    names = (
        pd.read_csv(ROOT / "oncoref/data/cancer-type-registry.csv")
        .set_index("code")["name"]
        .to_dict()
    )
    metrics = metrics.copy()
    metrics["cancer_name"] = metrics.cancer_code.map(names).fillna(metrics.cancer_code)
    groups, edges = overlap_groups(audit, cohorts)
    csvout(groups, out / "cohort_overlap_groups.csv")
    csvout(edges, out / "cohort_overlap_edges.csv")
    c = cohorts[
        [
            "cancer_code",
            "n_patients",
            "linear_tpm_comparable",
            "all_patient_ids_verified",
            "n_raw_samples",
            "n_qc_excluded",
            "source_cohort",
        ]
    ].copy()
    c["cancer_name"] = c.cancer_code.map(names)
    c["passes_default_size"] = c.n_patients >= min_size
    c["passes_large_linear_filter"] = (c.n_patients >= large_min_size) & c.linear_tpm_comparable
    csvout(c.sort_values("n_patients", ascending=False), out / "cohort_sizes.csv")
    availability = pd.DataFrame(
        [
            {
                "min_cohort_size": n,
                "n_cohorts": int((c.n_patients >= n).sum()),
                "n_linear_cohorts": int(((c.n_patients >= n) & c.linear_tpm_comparable).sum()),
            }
            for n in sorted({*SIZE_CUTS, min_size, large_min_size})
        ]
    )
    csvout(availability, out / "cohort_size_availability.csv")
    summary, sweep, gmts = [], [], {"all": [], "min_size": [], "large_linear": []}
    for p in PERCENTILES:
        data = metrics[metrics.transcriptome_percentile == p]
        for f in PREVALENCES:
            key = f"prevalence_gt{f}_transcriptome_p{p}"
            passing = data[passes(data, f)].copy()
            passing["passes_min_cohort_size"] = passing.n_patients >= min_size
            passing["passes_large_linear_filter"] = (
                passing.n_patients >= large_min_size
            ) & passing.linear_tpm_comparable
            passing["overlap_group"] = passing.cancer_code.map(
                groups.set_index("cancer_code").overlap_group
            )
            policies = {
                "all": passing,
                "min_size": passing[passing.passes_min_cohort_size],
                "large_linear": passing[passing.passes_large_linear_filter],
            }
            stats = {
                "scenario": key,
                "prevalence_gt_pct": f,
                "transcriptome_percentile": p,
                "min_cohort_size": min_size,
                "large_min_cohort_size": large_min_size,
            }
            for policy, selected in policies.items():
                gene_summary = summarize_genes(selected, groups, min_size, large_min_size)
                dest = out / "selections" / policy / key
                csvout(gene_summary, dest / "genes_ranked.csv")
                csvout(
                    selected.sort_values(["Symbol", "n_patients"], ascending=[True, False]),
                    dest / "passing_gene_cohorts.csv",
                )
                (dest / "symbols.txt").write_text(
                    "".join(s + "\n" for s in sorted(gene_summary.Symbol))
                )
                (dest / "ensembl_ids.txt").write_text(
                    "".join(s + "\n" for s in sorted(gene_summary[ID]))
                )
                stats[f"n_genes_{policy}"] = len(gene_summary)
                stats[f"n_cohorts_{policy}"] = selected.cancer_code.nunique()
                gmts[policy].append(
                    "\t".join(
                        [
                            key,
                            f"{policy}: prevalence >{f}%, transcriptome p{p}",
                            *sorted(gene_summary.Symbol),
                        ]
                    )
                )
            summary.append(stats)
            for n in sorted({*SIZE_CUTS, min_size, large_min_size}):
                for linear in (False, True):
                    selected = data[passes(data, f, n, linear)]
                    sweep.append(
                        {
                            "scenario": key,
                            "prevalence_gt_pct": f,
                            "transcriptome_percentile": p,
                            "min_cohort_size": n,
                            "linear_only": linear,
                            "n_genes": selected[ID].nunique(),
                            "n_passing_cohort_views": selected.cancer_code.nunique(),
                        }
                    )
    for policy, lines in gmts.items():
        (out / f"gene_sets_{policy}.gmt").write_text("\n".join(lines) + "\n")
    summary, sweep = pd.DataFrame(summary), pd.DataFrame(sweep)
    csvout(summary, out / "selection_summary.csv")
    csvout(sweep, out / "cohort_size_sensitivity.csv")
    strict = metrics[(metrics.transcriptome_percentile == 90) & passes(metrics, 75)].copy()
    strict["passes_default_size"] = strict.n_patients >= min_size
    strict["passes_large_linear_filter"] = (
        strict.n_patients.ge(large_min_size) & strict.linear_tpm_comparable
    )
    strict["overlap_group"] = strict.cancer_code.map(groups.set_index("cancer_code").overlap_group)
    csvout(
        strict.sort_values(["Symbol", "n_patients"], ascending=[True, False]),
        out / "strictest_passing_cancer_types.csv",
    )
    csvout(
        summarize_genes(strict, groups, min_size, large_min_size),
        out / "strictest_genes_ranked.csv",
    )
    return metrics, c, groups, edges, summary, sweep


def validate_outputs(out, summary, sweep, groups, min_size, large_min_size):
    for _, s in summary.iterrows():
        sets = {}
        for policy in ("all", "min_size", "large_linear"):
            folder = out / "selections" / policy / s.scenario
            genes = pd.read_csv(folder / "genes_ranked.csv")
            selected = pd.read_csv(folder / "passing_gene_cohorts.csv")
            sets[policy] = set(genes[ID])
            assert len(genes) == s[f"n_genes_{policy}"]
            assert set(selected[ID]) == sets[policy]
            for _, row in genes.iterrows():
                sub = selected[selected[ID].eq(row[ID])]
                assert sub.n_patients.max() == row.largest_passing_cohort_size
                assert sub.overlap_group.nunique() == row.n_overlap_adjusted_groups
            if policy == "min_size":
                assert selected.n_patients.ge(min_size).all()
            if policy == "large_linear":
                assert (
                    selected.n_patients.ge(large_min_size).all()
                    and selected.linear_tpm_comparable.all()
                )
        assert sets["min_size"] <= sets["all"] and sets["large_linear"] <= sets["all"]
    for _, sub in sweep.groupby(["scenario", "linear_only"]):
        assert (np.diff(sub.sort_values("min_cohort_size").n_genes) <= 0).all()
    for policy in ("all", "min_size", "large_linear"):
        sets = {
            (int(s.prevalence_gt_pct), int(s.transcriptome_percentile)): set(
                pd.read_csv(out / "selections" / policy / s.scenario / "genes_ranked.csv")[ID]
            )
            for _, s in summary.iterrows()
        }
        for p in PERCENTILES:
            for low, high in zip(PREVALENCES, PREVALENCES[1:]):
                assert sets[high, p] <= sets[low, p]
        for f in PREVALENCES:
            for low, high in zip(PERCENTILES, PERCENTILES[1:]):
                assert sets[f, high] <= sets[f, low]
    assert groups.cancer_code.is_unique
    return {
        "status": "passed",
        "checks": [
            "12 focused gene sets match supporting cohort rows",
            "minimum size applied to passing cohorts only",
            "largest passing cohort sizes verified",
            "linear-only eligibility verified",
            "overlap-adjusted breadth verified",
            "nested selections under all three thresholds",
        ],
    }


def make_plots(
    out, metrics, cohorts, groups, edges, summary, sweep, values, min_size, large_min_size
):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

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
    teal, orange, gray, purple = "#127C80", "#C8772B", "#BCC6CE", "#725C9E"
    dest = out / "plots"
    dest.mkdir(exist_ok=True)
    chart_index = []

    def save(fig, name, category, title):
        fig.savefig(dest / f"{name}.png", dpi=170, bbox_inches="tight")
        fig.savefig(dest / f"{name}.svg", bbox_inches="tight")
        plt.close(fig)
        chart_index.append({"name": name, "category": category, "title": title})

    # Size availability: empirical sizes after QC and donor grouping.
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    bins = [0, 9, 19, 49, 99, 249, 499, np.inf]
    labels = ["1-9", "10-19", "20-49", "50-99", "100-249", "250-499", "500+"]
    by_bin = pd.cut(cohorts.n_patients, bins=bins, labels=labels).value_counts(sort=False)
    axes[0].bar(labels, by_bin.values, color=teal)
    for i, n in enumerate(by_bin):
        axes[0].text(i, n + 0.5, str(n), ha="center", fontweight="bold")
    axes[0].set(
        xlabel="Patient groups per cohort",
        ylabel="Number of cohorts",
        title="A. Available cohort sizes",
    )
    availability = pd.read_csv(out / "cohort_size_availability.csv")
    for col, color, label in [
        ("n_cohorts", teal, "All scales"),
        ("n_linear_cohorts", orange, "Linear TPM only"),
    ]:
        axes[1].plot(
            availability.min_cohort_size, availability[col], marker="o", color=color, label=label
        )
    axes[1].set_xscale("log")
    axes[1].set(
        xticks=availability.min_cohort_size,
        xticklabels=availability.min_cohort_size.astype(str),
        xlabel="Minimum cohort size",
        ylabel="Eligible cohort views",
        title="B. Cohorts retained as the size cutoff increases",
    )
    axes[1].axvline(min_size, color=gray, ls="--")
    axes[1].axvline(large_min_size, color=orange, ls=":")
    axes[1].legend(frameon=False)
    fig.suptitle(
        f"142 cohorts | n={int(cohorts.n_patients.min())}-{int(cohorts.n_patients.max()):,} | median n={cohorts.n_patients.median():g}",
        fontsize=18,
        fontweight="bold",
    )
    save(fig, "cohort_size_overview", "main", "Available cohort sizes")

    ordered = cohorts.sort_values(["n_patients", "cancer_code"], ascending=[False, True])
    for page, start in enumerate(range(0, len(ordered), 48), 1):
        sub = ordered.iloc[start : start + 48].iloc[::-1]
        fig, ax = plt.subplots(figsize=(12, 12), layout="constrained")
        ax.barh(
            range(len(sub)), sub.n_patients, color=np.where(sub.linear_tpm_comparable, teal, orange)
        )
        for i, n in enumerate(sub.n_patients):
            ax.text(n * 1.08, i, str(n), va="center", fontsize=8)
        ax.set_xscale("log")
        ax.axvline(min_size, ls="--", color="black", lw=1, label=f"Default minimum: {min_size}")
        ax.axvline(
            large_min_size,
            ls=":",
            color=purple,
            lw=1,
            label=f"Larger-cohort minimum: {large_min_size}",
        )
        ax.set(
            yticks=range(len(sub)),
            yticklabels=[f"{r.cancer_code}: {r.cancer_name}" for _, r in sub.iterrows()],
            xlabel="Patient groups (log scale)",
            xlim=(0.8, 1600),
            title=f"All cohort sizes ({page}/3) | orange = proxy scale",
        )
        ax.tick_params(axis="y", labelsize=7.5)
        ax.legend(loc="lower right", frameon=False, fontsize=8)
        save(fig, f"cohort_sizes_{page}", "appendix", f"All cohort sizes {page}")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), layout="constrained")
    max_selected = summary.n_genes_all.max()
    for ax, policy, label in zip(
        axes,
        ("all", "min_size", "large_linear"),
        (
            "All passing cohorts",
            f"Minimum n={min_size}, all scales",
            f"Minimum n={large_min_size}, linear TPM",
        ),
    ):
        grid = summary.pivot(
            index="prevalence_gt_pct",
            columns="transcriptome_percentile",
            values=f"n_genes_{policy}",
        ).reindex(index=PREVALENCES, columns=PERCENTILES)
        ax.imshow(grid, cmap="GnBu", vmin=0, vmax=max_selected, aspect="auto")
        for i in range(len(PREVALENCES)):
            for j in range(len(PERCENTILES)):
                n = grid.iloc[i, j]
                ax.text(
                    j,
                    i,
                    str(n),
                    ha="center",
                    va="center",
                    fontsize=18,
                    color="white" if n > 0.65 * max_selected else "#203746",
                )
        ax.set(
            xticks=range(len(PERCENTILES)),
            xticklabels=[f"p{p}" for p in PERCENTILES],
            yticks=range(len(PREVALENCES)),
            yticklabels=[f">{f}%" for f in PREVALENCES],
            title=label,
            xlabel="Transcriptome expression cutoff",
            ylabel="Required cohort prevalence",
        )
    fig.suptitle("Selected CTA genes under three cohort policies", fontsize=19, fontweight="bold")
    save(fig, "selection_policy_comparison", "main", "Three cohort policies")

    fig, axes = plt.subplots(2, 2, figsize=(13, 8), layout="constrained", sharex=True, sharey=True)
    for i, f in enumerate(PREVALENCES):
        for j, p in enumerate(PERCENTILES):
            ax = axes[i, j]
            for linear, color, label in [(False, teal, "All scales"), (True, orange, "Linear TPM")]:
                sub = sweep[
                    sweep.prevalence_gt_pct.eq(f)
                    & sweep.transcriptome_percentile.eq(p)
                    & sweep.linear_only.eq(linear)
                ]
                ax.plot(sub.min_cohort_size, sub.n_genes, marker=".", color=color, label=label)
            ax.set_xscale("log")
            ax.axvline(min_size, color=gray, ls="--", lw=1)
            ax.axvline(large_min_size, color=orange, ls=":", lw=1)
            ax.set(title=f">{f}% above p{p}", ylim=(-2, 100))
            if i == len(PREVALENCES) - 1:
                ax.set_xlabel("Minimum cohort size")
            if j == 0:
                ax.set_ylabel("Selected genes")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle(
        "How every gene set changes with cohort-size requirements", fontsize=18, fontweight="bold"
    )
    save(fig, "cohort_size_sensitivity", "main", "Sensitivity to minimum cohort size")

    connected = (
        groups[groups.n_cohort_views_in_group.gt(1)]
        .drop_duplicates("overlap_group")
        .sort_values("n_cohort_views_in_group")
    )
    fig, ax = plt.subplots(figsize=(12, 6), layout="constrained")
    ax.barh(
        range(len(connected)),
        connected.n_cohort_views_in_group,
        color=gray,
        label="Original cohort views",
    )
    ax.barh(
        range(len(connected)),
        np.ones(len(connected)),
        color=teal,
        label="One overlap-adjusted group",
    )
    ax.set(
        yticks=range(len(connected)),
        yticklabels=connected.overlap_group,
        xlabel="Cohort views counted for breadth",
        title="142 cohort views form 108 groups after known-overlap adjustment",
    )
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(frameon=False, loc="lower right")
    fig.supxlabel(
        "Groups connect cohorts sharing known patient IDs, including subtype views connected through a parent.\nCounts are conservative support groups, not independent studies; unrecognized cross-source overlap can remain.",
        fontsize=9,
    )
    save(fig, "overlap_adjustment", "main", "Accounting for overlapping cohorts")

    strict = pd.read_csv(out / "strictest_passing_cancer_types.csv")
    six = pd.read_csv(out / "strictest_genes_ranked.csv").sort_values("largest_passing_cohort_size")
    fig, ax = plt.subplots(figsize=(12, 5), layout="constrained")
    ax.barh(
        six.Symbol,
        six.largest_passing_cohort_size,
        color=np.where(six.passes_large_linear_filter, teal, orange),
    )
    for i, (_, r) in enumerate(six.iterrows()):
        ax.text(
            r.largest_passing_cohort_size * 1.1,
            i,
            f"n={r.largest_passing_cohort_size} | {r.largest_passing_cohorts}",
            va="center",
            fontsize=10,
        )
    ax.set_xscale("log")
    ax.axvline(min_size, color="black", ls="--", label=f"Default minimum n={min_size}")
    ax.axvline(large_min_size, color=purple, ls=":", label=f"Larger minimum n={large_min_size}")
    ax.set(
        xlim=(0.8, 3000),
        xlabel="Largest passing cohort size (log scale)",
        title="Strictest list: the largest cohort that actually passes >75% / p90",
    )
    ax.legend(frameon=False, loc="lower right")
    fig.supxlabel(
        f"Teal: supported in a linear-TPM cohort with n >= {large_min_size}. Orange: only smaller or proxy-scale support.",
        fontsize=9,
    )
    save(fig, "strictest_largest_cohorts", "main", "Largest passing cohort for the six genes")

    gene_order = sorted(six.Symbol)
    codes = (
        strict.groupby("cancer_code").n_patients.max().sort_values(ascending=False).index.tolist()
    )
    cohort_lookup = cohorts.set_index("cancer_code")
    fig, ax = plt.subplots(figsize=(13, 11), layout="constrained")
    mat = strict.pivot(index="cancer_code", columns="Symbol", values="fraction_expressing").reindex(
        index=codes, columns=gene_order
    )
    im = ax.imshow(
        np.ma.masked_invalid(mat.to_numpy()), cmap="GnBu", vmin=0.75, vmax=1, aspect="auto"
    )
    for i, code in enumerate(codes):
        for j, gene in enumerate(gene_order):
            row = strict[(strict.cancer_code == code) & (strict.Symbol == gene)]
            if len(row):
                r = row.iloc[0]
                ax.text(
                    j,
                    i,
                    f"{int(r.n_expressing)}/{int(r.n_patients)}\n{r.fraction_expressing:.0%}",
                    ha="center",
                    va="center",
                    fontsize=9,
                    color="white" if r.fraction_expressing > 0.90 else "#203746",
                )
    labels = []
    for code in codes:
        r = cohort_lookup.loc[code]
        flags = (" [proxy]" if not r.linear_tpm_comparable else "") + (
            f" [n<{large_min_size}]" if r.n_patients < large_min_size else ""
        )
        labels.append(f"{r.cancer_name} ({code}){flags}")
    ax.set(
        xticks=range(6),
        xticklabels=gene_order,
        yticks=range(len(codes)),
        yticklabels=labels,
        title="All cancer cohorts passing >75% prevalence above p90",
    )
    ax.tick_params(axis="y", labelsize=8)
    fig.colorbar(im, ax=ax, shrink=0.6, label="Fraction expressing")
    fig.supxlabel(
        "Each filled cell is a passing gene/cohort pair. Cells show positive patient groups / cohort size. Blank = does not pass.\nProxy and small cohorts are shown explicitly; parent and subtype rows can overlap.",
        fontsize=9,
    )
    save(
        fig,
        "strictest_passing_cancer_types",
        "main",
        "Which cancer types pass for each of the six genes",
    )

    # One dashboard per expression/prevalence pair; every selected gene remains in CSV.
    for _, s in summary.iterrows():
        f, p = int(s.prevalence_gt_pct), int(s.transcriptome_percentile)
        all_genes = pd.read_csv(out / "selections/all" / s.scenario / "genes_ranked.csv")
        ranked = (
            pd.read_csv(out / "selections/large_linear" / s.scenario / "genes_ranked.csv")
            .head(15)
            .iloc[::-1]
        )
        eligible = metrics[
            metrics.transcriptome_percentile.eq(p)
            & metrics.linear_tpm_comparable
            & metrics.n_patients.ge(large_min_size)
            & metrics.complete_measurement
            & metrics.n_expressing.gt(0)
        ]
        fig, axes = plt.subplots(2, 2, figsize=(13, 9.5), layout="constrained")
        counts = [s.n_genes_all, s.n_genes_min_size, s.n_genes_large_linear]
        axes[0, 0].barh(range(3), counts, color=[gray, purple, teal])
        for i, n in enumerate(counts):
            axes[0, 0].text(n + 2, i, str(n), va="center", fontweight="bold")
        axes[0, 0].set(
            yticks=range(3),
            yticklabels=[
                "All passing cohorts",
                f"n >= {min_size}; all scales",
                f"n >= {large_min_size}; linear TPM",
            ],
            xlim=(0, max(counts) * 1.15 + 3),
            xlabel="Selected genes",
            title="A. Add the cohort-size and scale requirements",
        )
        axes[0, 0].invert_yaxis()
        axes[0, 1].scatter(
            all_genes.largest_passing_cohort_size,
            all_genes.n_overlap_adjusted_groups,
            c=np.where(all_genes.passes_large_linear_filter, teal, orange),
            s=20,
            alpha=0.6,
        )
        axes[0, 1].set_xscale("log")
        axes[0, 1].axvline(min_size, color="black", ls="--", lw=1)
        axes[0, 1].axvline(large_min_size, color=purple, ls=":", lw=1)
        axes[0, 1].set(
            xlabel="Largest passing cohort (patient groups)",
            ylabel="Passing overlap-adjusted groups",
            title="B. Largest support versus breadth",
        )
        axes[0, 1].yaxis.set_major_locator(MaxNLocator(integer=True))
        axes[0, 1].text(
            0.02,
            0.96,
            f"Teal = passes n >= {large_min_size} and linear TPM",
            transform=axes[0, 1].transAxes,
            va="top",
            fontsize=8,
        )
        axes[1, 0].barh(
            range(len(ranked)), ranked.n_passing_cohort_views, color=gray, label="Cohort views"
        )
        axes[1, 0].barh(
            range(len(ranked)),
            ranked.n_overlap_adjusted_groups,
            color=teal,
            label="Overlap-adjusted groups",
        )
        axes[1, 0].set(
            yticks=range(len(ranked)),
            yticklabels=ranked.Symbol,
            xlabel="Passing cohort views / groups",
            title=f"C. Broadest genes: n >= {large_min_size}, linear TPM",
        )
        axes[1, 0].xaxis.set_major_locator(MaxNLocator(integer=True))
        axes[1, 0].legend(frameon=False, fontsize=8, loc="lower right")
        positive = passes(eligible, f, large_min_size, True)
        for flag, color, alpha in [(False, gray, 0.25), (True, teal, 0.75)]:
            sub = eligible[positive.eq(flag)]
            axes[1, 1].scatter(
                100 * sub.fraction_expressing,
                sub.mean_expression_expressing,
                s=8 if flag else 4,
                c=color,
                alpha=alpha,
                rasterized=True,
            )
        axes[1, 1].axvline(f, color=orange, ls="--")
        axes[1, 1].set_yscale("log")
        axes[1, 1].set(
            xlabel="Positive patient groups (%)",
            ylabel="Mean positive-patient clean TPM",
            title=f"D. Expressing subset: n >= {large_min_size}, linear TPM",
            xlim=(0, 102),
        )
        fig.suptitle(
            f">{f}% above p{p} | {int(s.n_genes_large_linear)} genes with n >= {large_min_size}, linear TPM",
            fontsize=18,
            fontweight="bold",
        )
        save(fig, f"selection_{s.scenario}", "main", f"Selection and ranking: >{f}% / p{p}")

    # All-cohort prevalence heatmaps of the six genes, including failures.
    six_metrics = metrics[metrics.Symbol.isin(gene_order) & metrics.transcriptome_percentile.eq(90)]
    for page, start in enumerate(range(0, len(cohorts), 48), 1):
        subcodes = sorted(cohorts.cancer_code)[start : start + 48]
        mat = six_metrics.pivot(
            index="cancer_code", columns="Symbol", values="fraction_expressing"
        ).reindex(index=subcodes, columns=gene_order)
        fig, ax = plt.subplots(figsize=(11, 12), layout="constrained")
        im = ax.imshow(
            np.ma.masked_invalid(mat.to_numpy()), vmin=0, vmax=1, cmap="GnBu", aspect="auto"
        )
        for i, code in enumerate(subcodes):
            for j, gene in enumerate(gene_order):
                v = mat.loc[code, gene]
                if pd.notna(v):
                    ax.text(
                        j,
                        i,
                        f"{100 * v:.0f}%" + (" +" if v > 0.75 else ""),
                        ha="center",
                        va="center",
                        fontsize=7,
                        color="white" if v > 0.65 else "#203746",
                    )
        ax.set(
            xticks=range(6),
            xticklabels=gene_order,
            yticks=range(len(subcodes)),
            yticklabels=[f"{c} (n={int(cohort_lookup.loc[c, 'n_patients'])})" for c in subcodes],
            title=f"Six genes across all cohorts: fraction above p90 ({page}/3)",
        )
        ax.tick_params(axis="y", labelsize=8)
        fig.colorbar(im, ax=ax, shrink=0.6, label="Fraction expressing")
        fig.supxlabel(
            "+ marks prevalence >75%. Percentages are rounded only for display; unrounded values determine selection. Blank = unavailable.",
            fontsize=8,
        )
        save(
            fig,
            f"six_genes_all_cohorts_{page}",
            "main",
            f"Six-gene prevalence across all cohorts {page}",
        )

    # Patient-level plots preserve every group. X=expression / that patient's p90
    # is directly interpretable; raw expression retains its source-specific units.
    def distribution_plot(gene, codes, name, category, title):
        rng = np.random.default_rng(20260917)
        fig, axes = plt.subplots(
            1, 2, figsize=(14, max(4.5, 0.29 * len(codes) + 2)), layout="constrained", sharey=True
        )
        for i, code in enumerate(codes):
            sub = values[values.Symbol.eq(gene) & values.cancer_code.eq(code)]
            y = i + rng.uniform(-0.18, 0.18, len(sub))
            for positive, color in [(False, gray), (True, teal)]:
                m = sub.positive_p90.eq(positive).to_numpy()
                axes[0].scatter(
                    sub.loc[m, "expression_to_p90_ratio"],
                    y[m],
                    s=9,
                    alpha=0.45 if len(sub) > 100 else 0.75,
                    c=color,
                    edgecolors="none",
                    rasterized=True,
                )
                axes[1].scatter(
                    sub.loc[m, "expression"],
                    y[m],
                    s=9,
                    alpha=0.45 if len(sub) > 100 else 0.75,
                    c=color,
                    edgecolors="none",
                    rasterized=True,
                )
            pos = sub[sub.positive_p90]
            if len(pos):
                axes[1].scatter([pos.expression.mean()], [i], marker="D", c="black", s=23, zorder=5)
        labels = []
        for code in codes:
            r = cohort_lookup.loc[code]
            labels.append(
                f"{code} | n={int(r.n_patients)}"
                + (" | proxy" if not r.linear_tpm_comparable else "")
            )
        axes[0].set(
            yticks=range(len(codes)),
            yticklabels=labels,
            xlabel="Expression / patient's own p90 cutoff",
            title="Within-patient threshold",
        )
        axes[0].invert_yaxis()
        axes[0].set_xscale("symlog", linthresh=0.1)
        axes[0].axvline(1, color=orange, ls="--", lw=1)
        axes[0].set_xlim(left=0)
        axes[1].set_xscale("symlog", linthresh=0.1)
        axes[1].set_xlim(left=0)
        axes[1].set(
            xlabel="Clean TPM (proxy rows use source-normalized units)",
            title="Expression in every patient group",
        )
        axes[0].tick_params(axis="y", labelsize=8)
        fig.suptitle(title, fontsize=16, fontweight="bold")
        fig.supxlabel(
            "Every dot is one patient group: teal = above p90, gray = below or tied. Black diamond = mean among positive groups.\nBlank rows = gene unavailable. Symlog axes include zero. Proxy-scale absolute values are not comparable with TPM.",
            fontsize=9,
        )
        save(fig, name, category, title)

    for gene in gene_order:
        passing_codes = (
            strict[strict.Symbol.eq(gene)]
            .sort_values("n_patients", ascending=False)
            .cancer_code.tolist()
        )
        distribution_plot(
            gene,
            passing_codes,
            f"patients_{gene}_passing",
            "main",
            f"{gene}: patient-level expression in cohorts passing >75% / p90",
        )
        for page, start in enumerate(range(0, len(cohorts), 40), 1):
            subcodes = sorted(cohorts.cancer_code)[start : start + 40]
            distribution_plot(
                gene,
                subcodes,
                f"patients_{gene}_all_{page}",
                "appendix",
                f"{gene}: every cohort, including nonpassing cohorts ({page}/4)",
            )
    csvout(pd.DataFrame(chart_index), out / "plot_index.csv")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / "outputs/cta_threshold_report_20260916")
    parser.add_argument("--out", type=Path, default=ROOT / "outputs/cta_threshold_report_20260917")
    parser.add_argument(
        "--source-cache", type=Path, default=ROOT / "tmp/cta_threshold_report/source-matrices"
    )
    parser.add_argument("--min-cohort-size", type=int, default=10)
    parser.add_argument("--large-min-cohort-size", type=int, default=20)
    parser.add_argument("--reuse-patient-values", action="store_true")
    parser.add_argument("--tables-only", action="store_true")
    args = parser.parse_args()
    if min(args.min_cohort_size, args.large_min_cohort_size) < 1:
        parser.error("minimum cohort sizes must be positive")
    os.environ[CACHE_DIR_ENV_VAR] = str(args.source_cache.resolve())
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    metrics = pd.read_csv(args.base / "all_gene_cohort_metrics.csv.gz")
    cohorts = pd.read_csv(args.base / "cohort_audit.csv")
    audit = pd.read_csv(args.base / "sample_patient_audit.csv.gz")
    cutoffs = pd.read_csv(args.base / "patient_transcriptome_cutoffs.csv.gz")
    metrics, c, groups, edges, summary, sweep = build_tables(
        out, metrics, cohorts, audit, args.min_cohort_size, args.large_min_cohort_size
    )
    if args.tables_only:
        return
    values = extract_patient_values(out, metrics, audit, cutoffs, args.reuse_patient_values)
    checked = verify_patient_values(values, metrics)
    make_plots(
        out,
        metrics,
        c,
        groups,
        edges,
        summary,
        sweep,
        values,
        args.min_cohort_size,
        args.large_min_cohort_size,
    )
    validation = validate_outputs(
        out, summary, sweep, groups, args.min_cohort_size, args.large_min_cohort_size
    )
    validation["patient_expression_gene_cohort_checks"] = checked
    (out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (out / "run_manifest.json").write_text(
        json.dumps(
            {
                "analysis_date": "2026-09-17",
                "base_report": str(args.base.resolve()),
                "min_cohort_size": args.min_cohort_size,
                "large_min_cohort_size": args.large_min_cohort_size,
                "large_cohort_linear_only": True,
                "n_overlap_groups": groups.overlap_group.nunique(),
                "overlap_policy": "Connected components of known shared patient keys across all retained cohort views; one vote per component. Conservative, not independent-study count.",
                "unknown_overlap_policy": "Different physical sources are not assumed disjoint when donor mappings are unavailable; only observed shared identities can be removed.",
                "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "base_metrics_sha256": hashlib.sha256(
                    (args.base / "all_gene_cohort_metrics.csv.gz").read_bytes()
                ).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()

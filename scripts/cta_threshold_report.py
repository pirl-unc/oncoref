#!/usr/bin/env python3
"""One-off CTA prevalence / transcriptome-threshold analysis.

Run from this checkout with its Python environment. No package reference data is
modified. Per-cohort results are checkpointed; --resume reuses those checkpoints.
The report records its explicit analytical definitions in methodology.md.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from cta_report_common import (
    MIN_RANKED_PATIENTS,
    background_values,
    checkpoint_payload,
    checkpoint_valid,
    common_background,
    fingerprint,
    implementation_hash,
    seal_stage,
    write_checkpoint,
)
from cta_report_common import write_csv as reproducible_csv

from oncoref import __version__, source_matrices
from oncoref.cta import cta_gene_id_to_name
from oncoref.expression import (
    cancer_reference_expression_source_metadata,
    per_sample_expression,
    sample_expression_qc_from_matrix,
)
from oncoref.expression_engine import sample_columns
from oncoref.gene_families import clean_tpm_censored_gene_ids
from oncoref.gene_ids import canonical_gene_symbol, resolve_ensembl_id
from oncoref.normalization import clean_tpm
from oncoref.samples import sample_manifest, sample_molecular_provenance

PERCENTILES = (30, 50, 70, 90)
PREVALENCES = (25, 50, 75)
ID = "Ensembl_Gene_ID"
ROOT = Path(__file__).resolve().parents[1]


def write_csv(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    reproducible_csv(df, path, float_format="%.12g")


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def cta_universe():
    rows = []
    for gid, symbol in sorted(cta_gene_id_to_name().items()):
        canonical = resolve_ensembl_id(gid)
        rows.append(
            {
                ID: canonical,
                "Symbol": canonical_gene_symbol(canonical) or symbol,
                "cta_input_gene_id": gid,
                "cta_input_symbol": symbol,
            }
        )
    return pd.DataFrame(rows).drop_duplicates(ID).sort_values("Symbol").reset_index(drop=True)


def donor_maps(mtc_path=None):
    """Use physical-source manifests; never match unrelated cohorts by sample ID."""
    result = {}
    m = sample_manifest()
    m = m[m.included.astype(str).str.lower().isin(["true", "1"])]
    for _, r in m.dropna(subset=["case_id", "sample_id"]).iterrows():
        result[(str(r.source_cohort), str(r.sample_id))] = (str(r.case_id), "case_manifest")
    for _, r in sample_molecular_provenance().dropna(subset=["donor_id"]).iterrows():
        for col in ["sample_id", "library_id"]:
            if pd.notna(r[col]):
                result[(str(r.source_cohort), str(r[col]))] = (str(r.donor_id), "donor_manifest")
    # The cached MTC source titles identify repeat hybridizations of the same case.
    mtc = Path(mtc_path) if mtc_path else None
    if mtc is not None:
        if not mtc.is_file():
            raise FileNotFoundError(mtc)
        fields = {}
        with gzip.open(mtc, "rt") as f:
            for line in f:
                if line.startswith(("!Sample_geo_accession\t", "!Sample_title\t")):
                    row = next(csv.reader([line], delimiter="\t"))
                    fields[row[0]] = row[1:]
                if line.startswith("!series_matrix_table_begin"):
                    break
        for sample, title in zip(
            fields.get("!Sample_geo_accession", []), fields.get("!Sample_title", [])
        ):
            match = re.search(r"\[(MTC\d+)", title)
            if match:
                result[("GSE32662_PRINGLE_2012_MTC", sample)] = (match[1], "GEO_case_title")
    return result


def donor_identity(sample, source, mappings):
    if (source, sample) in mappings:
        return mappings[source, sample]
    if re.match(r"^(TCGA|TARGET)-[^-]+-[^-]+-", sample):
        return "-".join(sample.split("-")[:3]), "case_barcode"
    if re.match(r"^THR?\d+_[^_]+_S\d+", sample):
        return re.sub(r"_S\d+.*$", "", sample), "Treehouse_donor_prefix"
    return sample, "source_sample_unverified_patient"


def group_patients(values, samples, source, mappings):
    audit = pd.DataFrame(
        [
            {
                "sample_id": s,
                "patient_id": donor_identity(s, source, mappings)[0],
                "patient_id_basis": donor_identity(s, source, mappings)[1],
            }
            for s in samples
        ]
    )
    # Keep one TCGA specimen class per donor, preferring the primary tumor.
    # Metastatic and primary specimens are distinct biological states.
    audit["sample_type"] = audit.sample_id.map(
        lambda sample: (
            sample.split("-")[3][:2]
            if sample.startswith("TCGA-") and len(sample.split("-")) > 3
            else "unknown"
        )
    )
    audit["included_in_patient_group"] = True
    for _, group in audit.groupby("patient_id"):
        known = sorted(set(group.sample_type) - {"unknown"})
        if len(known) > 1:
            selected_type = "01" if "01" in known else known[0]
            audit.loc[group.index, "included_in_patient_group"] = group.sample_type.eq(
                selected_type
            )
    mask = audit.included_in_patient_group.to_numpy()
    frame = pd.DataFrame(values[:, mask], columns=audit.loc[mask, "patient_id"].tolist())
    grouped = frame.T.groupby(level=0, sort=True).mean().T
    return grouped.to_numpy(dtype=float), grouped.columns.tolist(), audit


def threshold_metrics(values, gene_rows, percentiles=PERCENTILES, background=None):
    """All measured biological genes x patient groups; CTA row indices separately.

    The percentile threshold is an interpolated expression quantile, not a
    percentile across patients. Strict > rejects ties at the quantile and zeros.
    Missing gene measurements remain missing, not zero.
    """
    if values.shape[1] == 0:
        raise ValueError("No evaluable patient groups")
    if np.any(values[np.isfinite(values)] < 0) or np.isinf(values).any():
        raise ValueError("Expression must be finite nonnegative values or missing")
    background = values if background is None else background
    cutoffs = np.nanquantile(background, np.asarray(percentiles) / 100, axis=0, method="linear")
    panel = values[gene_rows]
    # Average-rank percentiles are descriptive only. Selection uses value > quantile.
    ranks = pd.DataFrame(values).rank(axis=0, pct=True, method="average").to_numpy()[gene_rows]
    stats = []
    for p, cutoff in zip(percentiles, cutoffs):
        hits = np.isfinite(panel) & (panel > cutoff)
        n = hits.sum(axis=1)
        mean = np.divide(
            np.where(hits, panel, 0).sum(axis=1), n, out=np.full(len(panel), np.nan), where=n > 0
        )
        rank_mean = np.divide(
            np.where(hits, ranks, 0).sum(axis=1), n, out=np.full(len(panel), np.nan), where=n > 0
        )
        stats.append(
            pd.DataFrame(
                {
                    "transcriptome_percentile": p,
                    "n_patients": values.shape[1],
                    "n_measured_patients": np.isfinite(panel).sum(axis=1),
                    "n_expressing": n,
                    "fraction_expressing": n / values.shape[1],
                    "mean_expression_expressing": mean,
                    "mean_transcriptome_rank_expressing": rank_mean,
                    "mean_expression_all_patients": np.nanmean(panel, axis=1),
                    "complete_measurement": np.isfinite(panel).all(axis=1),
                }
            )
        )
    return stats, cutoffs


def analyze_cohort(code, universe, mappings, out, background_ids=None):
    info = source_matrices.cohort_info(code)
    meta = cancer_reference_expression_source_metadata(code)
    path = source_matrices.local_path(code)
    if not path.exists():
        raise FileNotFoundError(path)
    raw = per_sample_expression(code, normalize="tpm_raw", auto_fetch=False)
    samples = sample_columns(raw)
    if code == "MTC" and any((info["source_cohort"], sample) not in mappings for sample in samples):
        raise ValueError("MTC requires a donor title mapping for every source sample")
    qc = sample_expression_qc_from_matrix(raw, cancer_type=code)
    invalid = ((raw[samples] < 0) | np.isinf(raw[samples])).any(axis=0)
    invalid_samples = set(invalid.index[invalid])
    bad = qc.sample_id.isin(invalid_samples)
    qc.loc[bad, "sample_qc_status"] = "fail"
    qc.loc[bad, "passes_expression_qc"] = False
    qc.loc[bad, "sample_qc_reasons"] = (
        qc.loc[bad, "sample_qc_reasons"].fillna("").str.rstrip(";")
        + ";invalid_negative_or_infinite_expression"
    ).str.lstrip(";")
    qc["patient_id"] = [donor_identity(s, info["source_cohort"], mappings)[0] for s in qc.sample_id]
    qc["patient_id_basis"] = [
        donor_identity(s, info["source_cohort"], mappings)[1] for s in qc.sample_id
    ]
    qc["included_in_report"] = qc.sample_qc_status.isin(["pass", "warn"])
    kept = qc.loc[qc.included_in_report, "sample_id"].tolist()
    base = {
        "cancer_code": code,
        **info,
        **meta,
        "source_matrix_version": source_matrices.source_matrix_version(code),
        "source_matrix_url": source_matrices.release_url(code),
        "source_matrix_path": str(path.resolve()),
        "source_matrix_sha256": sha256(path),
        "n_raw_samples": len(samples),
        "n_qc_samples": len(kept),
        "n_qc_excluded": len(samples) - len(kept),
        "n_invalid_expression_samples": len(invalid_samples),
    }
    write_csv(qc, out / "checkpoints" / f"{code}_samples.csv")
    if not kept:
        base.update(status="no_QC_eligible_samples", n_patients=0)
        return pd.DataFrame(), base, pd.DataFrame()
    clean = clean_tpm(raw[kept], gene_table=raw[[ID, "Symbol"]])
    biological = ~raw[ID].astype(str).isin(clean_tpm_censored_gene_ids())
    ids = raw.loc[biological, ID].astype(str).reset_index(drop=True)
    values, patients, audit = group_patients(
        clean.loc[biological].to_numpy(), kept, info["source_cohort"], mappings
    )
    selected_samples = set(audit.loc[audit.included_in_patient_group, "sample_id"])
    qc["included_in_report"] &= qc.sample_id.isin(selected_samples)
    qc["specimen_policy_excluded"] = qc.sample_id.isin(set(kept) - selected_samples)
    write_csv(qc, out / "checkpoints" / f"{code}_samples.csv")
    base.update(
        n_specimen_type_excluded=len(kept) - len(selected_samples),
        n_patients=len(patients),
        n_repeat_samples_merged=len(selected_samples) - len(patients),
        n_biological_genes=len(ids),
        all_patient_ids_verified=not audit.patient_id_basis.eq(
            "source_sample_unverified_patient"
        ).any(),
        n_unverified_patient_groups=audit.loc[
            audit.patient_id_basis.eq("source_sample_unverified_patient"), "patient_id"
        ].nunique(),
        status="analyzed",
    )
    index = pd.Index(ids).get_indexer(universe[ID])
    found = index >= 0
    background = None
    if background_ids is not None:
        background = background_values(
            pd.DataFrame(values, index=ids, columns=patients).rename_axis(ID).reset_index(),
            background_ids,
            patients,
        )
        base["n_common_background_genes"] = len(background_ids)
    metrics, cutoffs = threshold_metrics(values, index[found], background=background)
    frames = []
    for p, stats in zip(PERCENTILES, metrics):
        observed = pd.concat(
            [universe.loc[found, [ID, "Symbol"]].reset_index(drop=True), stats], axis=1
        )
        full = universe[[ID, "Symbol"]].merge(observed, on=[ID, "Symbol"], how="left")
        full["gene_measured"] = found
        full["transcriptome_percentile"] = p
        full["n_patients"] = len(patients)
        full["complete_measurement"] = full.complete_measurement.fillna(False).astype(bool)
        for col in [
            "cancer_code",
            "source_cohort",
            "source_scale_class",
            "linear_tpm_comparable",
            "all_patient_ids_verified",
        ]:
            full[col] = base[col]
        full["expression_unit"] = (
            "biological_clean_TPM"
            if meta["linear_tpm_comparable"]
            else "clean_normalized_source_proxy"
        )
        # Only actual TPM-compatible means occupy the TPM column.
        full["mean_expressing_clean_tpm"] = full.mean_expression_expressing.where(
            full.linear_tpm_comparable
        )
        frames.append(full)
    thresholds = pd.DataFrame(cutoffs.T, columns=[f"p{p}_expression_cutoff" for p in PERCENTILES])
    thresholds.insert(0, "patient_id", patients)
    thresholds.insert(0, "cancer_code", code)
    thresholds["expression_unit"] = (
        "biological_clean_TPM" if meta["linear_tpm_comparable"] else "clean_normalized_source_proxy"
    )
    return pd.concat(frames, ignore_index=True), base, thresholds


def qualifies(frame, prevalence, minimum_patients=MIN_RANKED_PATIENTS):
    # Integer arithmetic avoids rounding a boundary case into the selected set.
    return (
        frame.complete_measurement
        & frame.n_patients.ge(minimum_patients)
        & (100 * frame.n_expressing > prevalence * frame.n_patients)
    )


def selection_outputs(metrics, universe, out):
    summary, memberships, gmt = [], [], []
    for p in PERCENTILES:
        data = metrics[metrics.transcriptome_percentile == p]
        for f in PREVALENCES:
            key = f"prevalence_gt{f}_transcriptome_p{p}"
            dest = out / "gene_sets" / key
            dest.mkdir(parents=True, exist_ok=True)
            selected = data[qualifies(data, f)].sort_values(["Symbol", "cancer_code"])
            gene_rows = []
            for gid, group in selected.groupby(ID, sort=False):
                best = group.sort_values(
                    ["fraction_expressing", "n_patients", "cancer_code"],
                    ascending=[False, False, True],
                ).iloc[0]
                robust = group[(group.n_patients >= 10) & group.linear_tpm_comparable]
                gene_rows.append(
                    {
                        ID: gid,
                        "Symbol": best.Symbol,
                        "n_qualifying_cohorts": group.overlap_group.nunique(),
                        "n_qualifying_cohort_views": len(group),
                        "qualifying_cohorts": ";".join(group.cancer_code),
                        "supported_in_n10_linear_cohort": not robust.empty,
                        "supported_in_verified_patient_cohort": bool(
                            group.all_patient_ids_verified.any()
                        ),
                        "best_cohort": best.cancer_code,
                        "best_cohort_n_patients": best.n_patients,
                        "best_cohort_n_expressing": best.n_expressing,
                        "best_cohort_fraction_expressing": best.fraction_expressing,
                        "best_cohort_mean_expressing": best.mean_expression_expressing,
                        "best_cohort_expression_unit": best.expression_unit,
                        "best_cohort_mean_positive_rank": best.mean_transcriptome_rank_expressing,
                    }
                )
            genes = pd.DataFrame(
                gene_rows,
                columns=[
                    ID,
                    "Symbol",
                    "n_qualifying_cohorts",
                    "n_qualifying_cohort_views",
                    "qualifying_cohorts",
                    "supported_in_n10_linear_cohort",
                    "supported_in_verified_patient_cohort",
                    "best_cohort",
                    "best_cohort_n_patients",
                    "best_cohort_n_expressing",
                    "best_cohort_fraction_expressing",
                    "best_cohort_mean_expressing",
                    "best_cohort_expression_unit",
                    "best_cohort_mean_positive_rank",
                ],
            )
            write_csv(genes, dest / "genes.csv")
            write_csv(selected, dest / "qualifying_gene_cohorts.csv")
            symbols = sorted(genes.Symbol.tolist())
            (dest / "symbols.txt").write_text("".join(s + "\n" for s in symbols))
            (dest / "ensembl_ids.txt").write_text("".join(s + "\n" for s in sorted(genes[ID])))
            gmt.append(
                "\t".join(
                    [key, f"CTA >{f}% patients above individual transcriptome p{p}", *symbols]
                )
            )
            robust = selected[(selected.n_patients >= 10) & selected.linear_tpm_comparable]
            verified = selected[selected.all_patient_ids_verified]
            summary.append(
                {
                    "scenario": key,
                    "prevalence_gt_pct": f,
                    "transcriptome_percentile": p,
                    "n_genes": len(genes),
                    "n_qualifying_cohorts": selected.overlap_group.nunique(),
                    "n_qualifying_cohort_views": selected.cancer_code.nunique(),
                    "interpretation": "exploratory low-expression screen; cutoff may be zero"
                    if p < 70
                    else "ranked selection",
                    "n_qualifying_gene_cohort_pairs": len(selected),
                    "n_genes_n10_linear": robust[ID].nunique(),
                    "n_genes_verified_patients": verified[ID].nunique(),
                }
            )
            member = universe[[ID, "Symbol"]].copy()
            member["scenario"] = key
            member["selected"] = member[ID].isin(genes[ID])
            member["selected_n10_linear"] = member[ID].isin(robust[ID])
            memberships.append(member)
    summary = pd.DataFrame(summary)
    write_csv(summary, out / "scenario_summary.csv")
    write_csv(pd.concat(memberships), out / "gene_set_membership.csv")
    (out / "all_gene_sets.gmt").write_text("\n".join(gmt) + "\n")
    return summary


def plots(metrics, universe, cohorts, summary, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

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
    teal, gray, orange = "#127C80", "#BBC2CC", "#D27B35"
    dest = out / "plots"
    dest.mkdir(exist_ok=True)

    def save(fig, name):
        fig.savefig(dest / (name + ".png"), dpi=180, bbox_inches="tight")
        fig.savefig(dest / (name + ".svg"), bbox_inches="tight")
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    for ax, metric, title in zip(
        axes,
        ["n_genes", "n_genes_n10_linear"],
        ["Any eligible cohort", "Sensitivity: linear TPM, at least 10 patient groups"],
    ):
        matrix = summary.pivot(
            index="prevalence_gt_pct", columns="transcriptome_percentile", values=metric
        ).reindex(index=PREVALENCES, columns=PERCENTILES)
        ax.imshow(matrix, cmap="GnBu", vmin=0, vmax=len(universe), aspect="auto")
        for i in range(3):
            for j in range(4):
                ax.text(
                    j,
                    i,
                    str(matrix.iloc[i, j]),
                    ha="center",
                    va="center",
                    fontsize=22,
                    color="white" if matrix.iloc[i, j] > len(universe) * 0.65 else "#162D40",
                )
        ax.set(
            xticks=range(4),
            xticklabels=[f"p{p}" for p in PERCENTILES],
            yticks=range(3),
            yticklabels=[f">{f}%" for f in PREVALENCES],
            xlabel="Patient's transcriptome expression cutoff",
            ylabel="Cohort prevalence required",
            title=title,
        )
    fig.suptitle("Number of selected CTA genes", fontsize=18, fontweight="bold")
    save(fig, "threshold_grid")

    cutoffs = pd.read_csv(out / "patient_transcriptome_cutoffs.csv.gz")
    linear_cutoffs = cutoffs[cutoffs.expression_unit.eq("biological_clean_TPM")]
    cutoff_cols = [f"p{p}_expression_cutoff" for p in PERCENTILES]
    cohort_medians = linear_cutoffs.groupby("cancer_code")[cutoff_cols].median()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
    rng = np.random.default_rng(20260916)
    for j, col in enumerate(cutoff_cols):
        vals = cohort_medians[col].to_numpy()
        axes[0].scatter(j + rng.uniform(-0.16, 0.16, len(vals)), vals, s=12, alpha=0.45, color=teal)
        axes[0].plot([j - 0.2, j + 0.2], [np.median(vals)] * 2, color=orange, lw=3)
    axes[0].set_yscale("symlog", linthresh=0.1)
    axes[0].set(
        xticks=range(4),
        xticklabels=[f"p{p}" for p in PERCENTILES],
        ylabel="Expression cutoff (clean TPM; symlog scale)",
        title="A. Median cutoff within each linear-TPM cohort",
    )
    zero = [linear_cutoffs[c].eq(0).mean() * 100 for c in cutoff_cols]
    axes[1].bar(range(4), zero, color=teal, width=0.6)
    for j, pct in enumerate(zero):
        axes[1].text(j, pct + 2, f"{pct:.1f}%", ha="center", fontweight="bold")
    axes[1].set(
        xticks=range(4),
        xticklabels=[f"p{p}" for p in PERCENTILES],
        ylim=(0, 105),
        ylabel="Cohort/patient entries with zero cutoff (%)",
        title="B. A zero cutoff means any positive expression",
    )
    fig.suptitle(
        "What the transcriptome percentiles mean in expression units",
        fontsize=16,
        fontweight="bold",
    )
    fig.supxlabel(
        "Dots: cohort medians. Orange: median across cohorts. Cohort/patient entries can overlap across subtype views.",
        fontsize=9,
    )
    save(fig, "expression_cutoff_diagnostics")

    for _, scenario in summary.iterrows():
        f, p = int(scenario.prevalence_gt_pct), int(scenario.transcriptome_percentile)
        data = metrics[metrics.transcriptome_percentile == p]
        selected = data[qualifies(data, f)]
        measured = data[data.complete_measurement]
        fig = plt.figure(figsize=(13, 10), layout="constrained")
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], width_ratios=[1.05, 1])
        ax = fig.add_subplot(gs[0, 0])
        counts = [
            len(universe),
            measured[ID].nunique(),
            measured.loc[measured.n_expressing > 0, ID].nunique(),
            selected[ID].nunique(),
        ]
        labels = [
            "Canonical CTA input",
            "Measured in >=1 cohort",
            f">p{p} in >=1 patient group",
            f">{f}% in >=1 cohort",
        ]
        ax.barh(range(4), counts, color=[gray, "#8DA6B7", "#589EA8", teal])
        for i, n in enumerate(counts):
            ax.text(n + 3, i, str(n), va="center", fontweight="bold")
        ax.set(
            yticks=range(4),
            yticklabels=labels,
            xlim=(0, len(universe) * 1.15),
            title="A. Selection funnel",
            xlabel="Unique genes",
        )
        ax.invert_yaxis()

        ax = fig.add_subplot(gs[0, 1])
        best = measured.sort_values(
            ["fraction_expressing", "n_patients", "cancer_code"], ascending=[False, False, True]
        ).drop_duplicates(ID)
        passing = qualifies(best, f)
        # All genes, including zero-positive genes (placed at y=0, explicitly stated).
        y = best.mean_transcriptome_rank_expressing.fillna(0)
        ax.scatter(
            best.fraction_expressing * 100,
            y * 100,
            c=np.where(passing, teal, gray),
            s=24,
            alpha=0.65,
            edgecolors="none",
        )
        ax.axvline(f, color=orange, ls="--", lw=1.5)
        ax.axhline(p, color=orange, ls=":", lw=1)
        ax.set(
            xlim=(-2, 103),
            ylim=(-3, 103),
            xlabel="Maximum cohort prevalence (%)",
            ylabel="Mean positive-patient gene rank (%)",
            title="B. Where each CTA passes",
        )
        ax.text(
            0.02,
            0.02,
            "One point per gene; strongest-prevalence cohort.\nNo positives: plotted at rank 0. Teal = selected.",
            transform=ax.transAxes,
            fontsize=8,
            va="bottom",
        )

        ax = fig.add_subplot(gs[1, 0])
        sizes = (
            selected.groupby("cancer_code")[ID]
            .nunique()
            .sort_values(ascending=False)
            .head(18)
            .sort_values()
        )
        n_lookup = cohorts.set_index("cancer_code").n_patients
        ax.barh(range(len(sizes)), sizes.values, color=teal)
        ax.set(
            yticks=range(len(sizes)),
            yticklabels=[f"{c} (n={int(n_lookup[c])})" for c in sizes.index],
            xlabel="Selected genes in this cohort",
            title="C. Cohorts contributing the most genes",
        )
        ax.tick_params(axis="y", labelsize=8)

        ax = fig.add_subplot(gs[1, 1])
        linear = measured[measured.linear_tpm_comparable & (measured.n_expressing > 0)]
        hits = qualifies(linear, f)
        for mask, color, label in [
            (~hits, gray, "Below prevalence cutoff"),
            (hits, teal, "Qualifying gene-cohort pair"),
        ]:
            sub = linear[mask]
            ax.scatter(
                sub.fraction_expressing * 100,
                sub.mean_expression_expressing,
                s=12 if label == "Qualifying gene-cohort pair" else 5,
                alpha=0.75 if label == "Qualifying gene-cohort pair" else 0.25,
                c=color,
                label=label,
                rasterized=True,
                edgecolors="none",
            )
        ax.axvline(f, color=orange, ls="--", lw=1.5)
        ax.set_yscale("log")
        ax.set(
            xlim=(-2, 103),
            xlabel="Patients above their transcriptome cutoff (%)",
            ylabel="Mean clean TPM among positive patients",
            title="D. Expression in the expressing subset",
        )
        ax.legend(fontsize=7, loc="lower right", frameon=False)
        ax.text(
            0.02,
            0.98,
            "Linear TPM cohorts only; each dot = gene x cohort",
            transform=ax.transAxes,
            fontsize=8,
            va="top",
        )
        fig.suptitle(
            f">{f}% of a cohort above transcriptome p{p}  |  {int(scenario.n_genes)} genes",
            fontsize=17,
            fontweight="bold",
        )
        save(fig, scenario.scenario)

    # A compact overview keeps every cohort visible instead of only the top 18.
    codes = cohorts.loc[cohorts.status.eq("analyzed"), "cancer_code"].tolist()
    matrix = pd.DataFrame(0, index=codes, columns=summary.scenario)
    for _, s in summary.iterrows():
        data = metrics[metrics.transcriptome_percentile == s.transcriptome_percentile]
        cnt = data[qualifies(data, s.prevalence_gt_pct)].groupby("cancer_code")[ID].nunique()
        matrix[s.scenario] = cnt.reindex(codes).fillna(0).astype(int)
    write_csv(matrix.rename_axis("cancer_code").reset_index(), out / "cohort_selection_counts.csv")
    for page, start in enumerate(range(0, len(codes), 48), 1):
        subset = matrix.iloc[start : start + 48]
        fig, ax = plt.subplots(figsize=(11, 12), layout="constrained")
        im = ax.imshow(subset, cmap="GnBu", vmin=0, vmax=matrix.to_numpy().max(), aspect="auto")
        ax.set(
            xticks=range(12),
            xticklabels=[
                f"p{int(s.transcriptome_percentile)}\n>{int(s.prevalence_gt_pct)}%"
                for _, s in summary.iterrows()
            ],
            yticks=range(len(subset)),
            yticklabels=[
                f"{c}  n={int(cohorts.set_index('cancer_code').loc[c, 'n_patients'])}"
                for c in subset.index
            ],
            title=f"Selected CTA count in every cohort ({page})",
        )
        ax.tick_params(labelsize=8)
        for i in range(len(subset)):
            for j in range(12):
                ax.text(
                    j,
                    i,
                    str(subset.iloc[i, j]),
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="white"
                    if subset.iloc[i, j] > matrix.to_numpy().max() * 0.6
                    else "#173247",
                )
        fig.colorbar(im, ax=ax, shrink=0.5, label="Selected genes")
        save(fig, f"all_cohort_counts_{page}")


def validate(metrics, summary, out):
    valid = metrics[metrics.complete_measurement]
    assert not metrics.duplicated([ID, "cancer_code", "transcriptome_percentile"]).any()
    assert valid.n_expressing.between(0, valid.n_patients).all()
    assert np.allclose(valid.fraction_expressing, valid.n_expressing / valid.n_patients)
    assert valid.loc[valid.n_expressing.eq(0), "mean_expression_expressing"].isna().all()
    assert (valid.loc[valid.n_expressing.gt(0), "mean_expression_expressing"] > 0).all()
    piv = valid.pivot(
        index=[ID, "cancer_code"], columns="transcriptome_percentile", values="n_expressing"
    )
    assert (np.diff(piv.to_numpy(), axis=1) <= 0).all()
    selected_sets = {}
    for _, s in summary.iterrows():
        frame = pd.read_csv(out / "gene_sets" / s.scenario / "genes.csv")
        selected_sets[s.prevalence_gt_pct, s.transcriptome_percentile] = set(frame[ID])
        assert len(frame) == s.n_genes
    for f in PREVALENCES:
        for low, high in zip(PERCENTILES, PERCENTILES[1:]):
            assert selected_sets[f, high] <= selected_sets[f, low]
    for p in PERCENTILES:
        for low, high in zip(PREVALENCES, PREVALENCES[1:]):
            assert selected_sets[high, p] <= selected_sets[low, p]
    (out / "validation.json").write_text(
        json.dumps(
            {
                "status": "passed",
                "checks": [
                    "unique gene/cohort/percentile keys",
                    "numerator <= denominator",
                    "fraction arithmetic",
                    "empty positive subsets have missing means",
                    "positive subset means > 0",
                    "positive counts nonincreasing with percentile",
                    "12 exports match summary",
                    "gene sets nested under both threshold axes",
                ],
            },
            indent=2,
        )
        + "\n"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "outputs/cta_threshold_report_20260916")
    parser.add_argument("--source-cache", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--mtc-donor-matrix",
        type=Path,
        help="GSE32662 series matrix containing MTC donor titles; required when MTC is included",
    )
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    if args.source_cache:
        os.environ[source_matrices.CACHE_DIR_ENV_VAR] = str(args.source_cache.resolve())
    out = args.out.resolve()
    (out / "checkpoints").mkdir(parents=True, exist_ok=True)
    universe = cta_universe()
    write_csv(universe, out / "cta_input_universe.csv")
    codes = source_matrices.available_cohorts()
    if "MTC" in codes and args.mtc_donor_matrix is None:
        parser.error(
            "--mtc-donor-matrix is required: MTC repeat hybridizations need verified donor mapping"
        )
    script_hash = sha256(Path(__file__))
    common_hash = sha256(ROOT / "scripts/cta_report_common.py")
    mappings = donor_maps(args.mtc_donor_matrix)
    background_ids, source_hashes = common_background(codes)
    write_csv(pd.DataFrame({ID: background_ids}), out / "common_background.csv")
    policy = fingerprint(
        {
            "implementation": implementation_hash(),
            "background": background_ids,
            "universe": universe.to_dict("records"),
            "donors": sorted((list(k), v) for k, v in mappings.items()),
        }
    )
    all_metrics, cohort_rows, all_cutoffs = [], [], []
    for i, code in enumerate(source_matrices.available_cohorts(), 1):
        cp = out / "checkpoints" / code
        expected = fingerprint([policy, code, source_hashes[code]])
        if args.resume and checkpoint_valid(
            cp, expected, [".parquet", "_cutoffs.csv", "_samples.csv"]
        ):
            cohort = checkpoint_payload(cp)
            data = pd.read_parquet(cp.with_suffix(".parquet"))
            cutoffs = (
                pd.read_csv(cp.with_name(code + "_cutoffs.csv"))
                if cp.with_name(code + "_cutoffs.csv").exists()
                else pd.DataFrame()
            )
        else:
            print(f"[{i}/142] {code}: computing", flush=True)
            data, cohort, cutoffs = analyze_cohort(code, universe, mappings, out, background_ids)
            data.to_parquet(cp.with_suffix(".parquet"), index=False)
            cp.with_suffix(".json").write_text(json.dumps(cohort, indent=2, default=str) + "\n")
            write_csv(cutoffs, cp.with_name(code + "_cutoffs.csv"))
            write_checkpoint(cp, cohort, expected, [".parquet", "_cutoffs.csv", "_samples.csv"])
        all_metrics.append(data)
        cohort_rows.append(cohort)
        all_cutoffs.append(cutoffs)
        print(
            f"[{i}/142] {code}: {cohort.get('n_patients', 0)} patient groups, {cohort['status']}",
            flush=True,
        )
    metrics = pd.concat(all_metrics, ignore_index=True)
    cohorts = pd.DataFrame(cohort_rows)
    cohorts["n_invalid_expression_samples"] = cohorts.n_invalid_expression_samples.fillna(0).astype(
        int
    )
    write_csv(metrics, out / "all_gene_cohort_metrics.csv.gz")
    write_csv(cohorts, out / "cohort_audit.csv")
    write_csv(pd.concat(all_cutoffs), out / "patient_transcriptome_cutoffs.csv.gz")
    write_csv(
        pd.concat([pd.read_csv(p) for p in sorted((out / "checkpoints").glob("*_samples.csv"))]),
        out / "sample_patient_audit.csv.gz",
    )
    from cta_threshold_report_details import overlap_groups

    patient_audit = pd.read_csv(out / "sample_patient_audit.csv.gz")
    groups, edges = overlap_groups(patient_audit, cohorts)
    write_csv(groups, out / "cohort_overlap_groups.csv")
    write_csv(edges, out / "cohort_overlap_edges.csv")
    metrics["overlap_group"] = metrics.cancer_code.map(
        groups.set_index("cancer_code").overlap_group
    )
    write_csv(metrics, out / "all_gene_cohort_metrics.csv.gz")
    summary = selection_outputs(metrics, universe, out)
    validate(metrics, summary, out)
    if not args.no_plots:
        plots(metrics, universe, cohorts, summary, out)
    # Provenance is best-effort: the report must still run from a source
    # tarball or an installed copy, where there is no git checkout to ask.
    git_commit = None
    with contextlib.suppress(OSError, subprocess.CalledProcessError):
        git_commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            cwd=ROOT,
        ).stdout.strip()
    if (
        sha256(Path(__file__)) != script_hash
        or sha256(ROOT / "scripts/cta_report_common.py") != common_hash
    ):
        raise RuntimeError("Analysis code changed during this run; rebuild before publishing")
    manifest = {
        "oncoref_version": __version__,
        "git_commit": git_commit,
        "script_sha256": script_hash,
        "common_script_sha256": common_hash,
        "percentiles": PERCENTILES,
        "minimum_patients_for_selection": MIN_RANKED_PATIENTS,
        "exploratory_percentiles": [30, 50],
        "common_background_genes": len(background_ids),
        "common_background_sha256": sha256(out / "common_background.csv"),
        "donor_mapping_inputs": {str(args.mtc_donor_matrix): sha256(args.mtc_donor_matrix)}
        if args.mtc_donor_matrix
        else {},
        "prevalence_gt_pct": PREVALENCES,
        "cta_genes": len(universe),
        "cohorts_registered": len(cohorts),
        "cohorts_analyzed": int(cohorts.status.eq("analyzed").sum()),
        "qc_policy": "pass_or_warn; proxy-scale warning retained and labeled",
        "invalid_expression_policy": "Exclude samples with negative or infinite raw expression",
        "definition": "expression > linear-interpolated quantile of a fixed common biological gene background, within patient group; p30/p50 are exploratory and can have zero cutoffs",
        "repeat_sample_policy": "one TCGA specimen class per donor, preferring primary tumor; arithmetic mean within retained class",
        "mean_definition": "arithmetic linear-scale mean among positive patient groups only",
        "prevalence_denominator": "all QC-eligible patient groups within each cohort; unknown patient IDs use source samples",
        "missing_policy": "missing gene/cohort measurements never imputed zero; only complete CTA measurements qualify",
        "input_table_sha256": {
            n: sha256(ROOT / "oncoref/data" / n)
            for n in [
                "cancer-testis-antigens.csv",
                "cta-specificity-audit.csv",
                "source-matrices.csv",
                "expression_sources.yaml",
            ]
        },
    }
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "methodology.md").write_text(
        "# CTA gene report methodology\n\n" + json.dumps(manifest, indent=2) + "\n"
    )
    seal_stage(
        out,
        "analysis",
        [
            Path(__file__),
            ROOT / "scripts/cta_report_common.py",
            *[source_matrices.local_path(code) for code in codes],
            *([args.mtc_donor_matrix] if args.mtc_donor_matrix else []),
        ],
        [
            out / n
            for n in [
                "validation.json",
                "run_manifest.json",
                "common_background.csv",
                "all_gene_cohort_metrics.csv.gz",
                "cohort_audit.csv",
                "sample_patient_audit.csv.gz",
                "scenario_summary.csv",
                "patient_transcriptome_cutoffs.csv.gz",
                "cta_input_universe.csv",
                "cohort_overlap_groups.csv",
            ]
        ]
        + sorted((out / "gene_sets").rglob("*.csv")),
    )
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()

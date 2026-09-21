#!/usr/bin/env python3
"""Recompute the CTA screen in genome-wide identical-protein expression space."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from cta_report_common import (
    MIN_RANKED_PATIENTS,
    RANKED_COHORT_POLICY,
    background_values,
    checkpoint_payload,
    checkpoint_valid,
    fingerprint,
    implementation_hash,
    ranked_selection_dir,
    seal_stage,
    verify_analysis,
    write_checkpoint,
    write_csv,
)

from oncoref import __version__
from oncoref.expression import per_sample_expression
from oncoref.gene_families import clean_tpm_censored_gene_ids
from oncoref.proteoforms import collapse_to_proteoforms, gene_to_proteoform_id, proteoform_groups
from oncoref.source_matrices import CACHE_DIR_ENV_VAR, local_path

ROOT = Path(__file__).resolve().parents[1]
ID = "Ensembl_Gene_ID"
KEY = "proteoform_key"
PERCENTILES = (70, 90)
FOCUSED = [(50, 70), (70, 70), (50, 90), (70, 90)]
POLICIES = {RANKED_COHORT_POLICY: (MIN_RANKED_PATIENTS, False)}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def csv(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    write_csv(frame, path)


def universe_tables(base, out, *, fasta=None):
    # Importing here lets focused numerical tests load this module by file path.
    from cta_report_protein_lengths import build_annotations

    seeds = pd.read_csv(base / "cta_input_universe.csv")[[ID, "Symbol"]]
    registry = proteoform_groups(scope="genome")
    labels = set(registry.loc[registry.member_gene_id.isin(seeds[ID]), "proteoform_id"])
    grouped = registry[registry.proteoform_id.isin(labels)]
    expanded = pd.concat(
        [
            seeds,
            grouped[["member_gene_id", "member_symbol"]].rename(
                columns={"member_gene_id": ID, "member_symbol": "Symbol"}
            ),
        ]
    ).drop_duplicates(ID)
    fasta = fasta or (
        Path.home() / "Library/Caches/pyensembl/GRCh38/ensembl112/Homo_sapiens.GRCh38.pep.all.fa.gz"
    )
    proteins, isoforms = build_annotations(expanded, fasta)
    mapping = gene_to_proteoform_id(expanded[ID], scope="genome")
    expanded[KEY] = expanded[ID].map(mapping)
    expanded["in_original_cta_set"] = expanded[ID].isin(seeds[ID])
    hpa = pd.read_csv(ROOT / "oncoref/data/cancer-testis-antigens.csv").drop_duplicates(ID)
    hpa_cols = [
        ID,
        "restriction",
        "restriction_confidence",
        "rna_restriction_level",
        "rna_max_somatic_ntpm",
        "rna_max_somatic_tissue",
        "safety_flags",
        "protein_restriction",
    ]
    member_evidence = expanded.merge(proteins.drop(columns="Symbol"), on=ID, validate="one_to_one")
    member_evidence = member_evidence.merge(hpa[hpa_cols], on=ID, how="left", validate="one_to_one")
    from oncoref import hpa
    from oncoref.cta_tissues import NON_SOMATIC_TISSUES

    normal_atlas = hpa._read_hpa("hpa_rna_consensus", "v23")
    normal_atlas = normal_atlas[~normal_atlas.Tissue.isin(NON_SOMATIC_TISSUES)]
    normal_wide = normal_atlas.pivot(index="Gene", columns="Tissue", values="nTPM")
    rows = []
    for key, members in expanded.groupby(KEY, sort=True):
        pp = proteins[proteins[ID].isin(members[ID])]
        # A longest isoform can tie; select a common longest sequence, not merely equal lengths.
        candidates = isoforms[isoforms[ID].isin(members[ID])].copy()
        candidates = candidates[
            candidates.protein_length_aa.eq(candidates[ID].map(pp.set_index(ID).protein_length_aa))
        ]
        common = set.intersection(
            *(set(s.protein_sequence_sha256) for _, s in candidates.groupby(ID))
        )
        assert common, f"Registry members lack an identical longest sequence: {key}"
        sequence_hash = sorted(common)[0]
        chosen = candidates[candidates.protein_sequence_sha256.eq(sequence_hash)].sort_values(
            "protein_id"
        )
        display = key if len(members) > 1 else members.Symbol.iloc[0]
        if len(display) > 30:
            display = f"CT47A ({len(members)} loci)" if display.startswith("CT47A") else display
        evidence = member_evidence[member_evidence[ID].isin(members[ID])]
        normal = pd.to_numeric(evidence.rna_max_somatic_ntpm, errors="coerce")
        top_normal = evidence.loc[normal.idxmax()] if normal.notna().any() else None
        normal_sum = normal_wide.reindex(members[ID]).sum(axis=0, min_count=len(members))
        normal_max = normal_sum.max()
        rows.append(
            {
                KEY: key,
                "normal_rna_max_summed_somatic_ntpm": normal_max,
                "normal_rna_max_summed_status": "unavailable"
                if pd.isna(normal_max)
                else "reported_zero"
                if normal_max == 0
                else "positive_estimate",
                "normal_rna_zero_upper_bound_ntpm": len(members) * 0.05
                if normal_max == 0
                else np.nan,
                "normal_rna_summed_max_tissues": ";".join(
                    sorted(normal_sum.index[normal_sum.eq(normal_max)])
                ),
                "normal_rna_summed_tissues_measured": int(normal_sum.notna().sum()),
                "Symbol": display,
                "member_symbols": ";".join(sorted(members.Symbol)),
                "member_gene_ids": ";".join(sorted(members[ID])),
                "n_member_genes": len(members),
                "n_original_cta_members": int(members.in_original_cta_set.sum()),
                "additional_member_symbols": ";".join(
                    sorted(members.loc[~members.in_original_cta_set, "Symbol"])
                ),
                "protein_length_aa": int(chosen.protein_length_aa.iloc[0]),
                "protein_ids": ";".join(sorted(chosen.protein_id.unique())),
                "protein_sequence_sha256": sequence_hash,
                "max_member_isoform_count": int(pp.protein_isoform_count.max()),
                "member_restriction_annotations": ";".join(
                    f"{r.Symbol}:{r.restriction}"
                    for _, r in evidence.iterrows()
                    if pd.notna(r.restriction)
                ),
                "member_restriction_confidence": ";".join(
                    sorted(evidence.restriction_confidence.dropna().unique())
                ),
                "normal_rna_max_member_somatic_ntpm": normal.max(),
                "normal_rna_max_member_and_tissue": f"{top_normal.Symbol}:{top_normal.rna_max_somatic_tissue}"
                if top_normal is not None
                else "",
                "normal_rna_members_with_evidence": int(normal.notna().sum()),
                "member_annotation_flags": ";".join(
                    sorted(evidence.safety_flags.dropna().astype(str).unique())
                ),
            }
        )
    universe = pd.DataFrame(rows).sort_values("Symbol").reset_index(drop=True)
    assert len(universe) == len(set(gene_to_proteoform_id(seeds[ID], scope="genome").values()))
    assert universe[KEY].is_unique
    csv(universe, out / "proteoform_universe.csv")
    csv(member_evidence, out / "proteoform_member_annotations.csv")
    csv(isoforms, out / "member_protein_isoforms.csv")
    csv(registry, out / "genome_identical_protein_registry.csv")
    csv(proteoform_groups(scope="cta"), out / "cta_identical_protein_registry.csv")
    csv(expanded, out / "gene_to_proteoform_mapping.csv")
    (out / "annotation_manifest.json").write_text(
        json.dumps(
            {
                "ensembl_release": 112,
                "assembly": "GRCh38",
                "protein_fasta_sha256": sha(fasta),
                "protein_fasta_path": str(fasta.resolve()),
                "protein_source": "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz",
                "genome_registry_sha256": sha(ROOT / "oncoref/data/proteoform-groups-genome.csv"),
                "normal_tissue_annotation_sha256": sha(
                    ROOT / "oncoref/data/cancer-testis-antigens.csv"
                ),
                "cancer_type_registry_sha256": sha(ROOT / "oncoref/data/cancer-type-registry.csv"),
                "n_original_cta_genes": len(seeds),
                "n_member_genes": len(expanded),
                "n_proteoforms": len(universe),
                "n_multigene_proteoforms": int(universe.n_member_genes.gt(1).sum()),
                "sequence_identity_validated": True,
                "identity_definition": "Identical longest Ensembl 112 protein sequence; not an isoform-resolved or post-translational proteomics measurement.",
                "normal_tissue_definition": "Primary RNA axis is maximum within-tissue sum across identical-protein loci in pinned HPA v23, requiring all members measured. NON_SOMATIC_TISSUES excluded. Member maximum retained as diagnostic. Reported zero is rounded; upper bound 0.05 nTPM per member. Neither axis is a safety score.",
            },
            indent=2,
        )
        + "\n"
    )
    return universe


def panel_metrics(background, panel, full_members, percentiles=PERCENTILES):
    """Recompute per-patient quantiles after collapse; incomplete sums cannot qualify."""
    if np.isinf(background).any() or np.any(background[np.isfinite(background)] < 0):
        raise ValueError("Invalid expression background")
    cuts = np.nanquantile(background, np.asarray(percentiles) / 100, axis=0, method="linear")
    available = np.isfinite(panel).all(axis=1)
    frames = []
    for p, cutoff in zip(percentiles, cuts):
        positive = np.isfinite(panel) & (panel > cutoff)
        n = positive.sum(axis=1)
        means = np.divide(
            np.where(positive, panel, 0).sum(axis=1),
            n,
            out=np.full(len(panel), np.nan),
            where=n > 0,
        )
        medians = np.array(
            [np.median(row[hit]) if hit.any() else np.nan for row, hit in zip(panel, positive)]
        )
        all_means = np.divide(
            np.nansum(panel, axis=1),
            np.isfinite(panel).sum(axis=1),
            out=np.full(len(panel), np.nan),
            where=np.isfinite(panel).sum(axis=1) > 0,
        )
        frames.append(
            pd.DataFrame(
                {
                    "transcriptome_percentile": p,
                    "n_patients": panel.shape[1],
                    "n_measured_patients": np.isfinite(panel).sum(axis=1),
                    "n_expressing": n,
                    "fraction_expressing": n / panel.shape[1],
                    "mean_expression_expressing": means,
                    "median_expression_expressing": medians,
                    "mean_expression_all_patients": all_means,
                    "available_sum_complete_measurement": available,
                    "complete_member_coverage": full_members,
                    "complete_measurement": available & full_members,
                }
            )
        )
    return frames, cuts


def analyze_cohort(code, universe, audit, cohort, out, background_ids=None):
    path = local_path(code)
    assert sha(path) == cohort.source_matrix_sha256, f"Changed source matrix: {code}"
    keep = audit[audit.included_in_report].copy()
    clean = per_sample_expression(code, normalize="tpm_clean", auto_fetch=False)
    samples = keep.sample_id.tolist()
    clean = clean[~clean[ID].isin(clean_tpm_censored_gene_ids())].reset_index(drop=True)
    assert set(samples) <= set(clean.columns)
    assert clean[ID].is_unique
    grouped = clean[samples].T
    grouped.index = keep.patient_id.to_numpy()
    grouped = grouped.groupby(level=0, sort=True).mean().T
    patients = grouped.columns.tolist()
    assert len(patients) == cohort.n_patients
    genes = clean[[ID, "Symbol"]].join(grouped)
    collapsed = collapse_to_proteoforms(genes, scope="genome", sample_cols=patients)
    assert collapsed[KEY].is_unique
    original_total = genes[patients].sum(axis=0).to_numpy()
    assert np.allclose(original_total, collapsed[patients].sum(axis=0).to_numpy(), rtol=1e-12)
    index = collapsed.set_index(KEY)
    panel = index.reindex(universe[KEY])[patients].to_numpy(dtype=float)
    measured = set(genes[ID])
    gene_values = genes.set_index(ID)[patients]
    full, member_counts, missing = [], [], []
    for members in universe.member_gene_ids:
        ids = members.split(";")
        found = [g for g in ids if g in measured]
        member_counts.append(len(found))
        missing.append(";".join(g for g in ids if g not in measured))
        full.append(len(found) == len(ids) and np.isfinite(gene_values.loc[found].to_numpy()).all())
    full = np.asarray(full)
    background = collapsed[patients].to_numpy()
    if background_ids is not None:
        # Collapse the same loci in every cohort; do not bring extra paralogs
        # from larger source matrices into the percentile background.
        background_genes = genes.set_index(ID).loc[background_ids].reset_index()
        fixed = collapse_to_proteoforms(background_genes, scope="genome", sample_cols=patients)
        background = background_values(fixed, sorted(fixed[KEY]), patients, KEY)
    frames, cuts = panel_metrics(background, panel, full)
    result = []
    for stats in frames:
        frame = pd.concat([universe[[KEY, "Symbol"]], stats], axis=1)
        frame["cancer_code"] = code
        frame["n_member_genes_observed"] = member_counts
        frame["missing_member_gene_ids"] = missing
        frame["linear_tpm_comparable"] = bool(cohort.linear_tpm_comparable)
        frame["expression_unit"] = (
            "biological_clean_TPM" if cohort.linear_tpm_comparable else "source_proxy"
        )
        result.append(frame)
    cp = out / "checkpoints" / code
    pd.DataFrame(panel, index=universe[KEY], columns=patients).rename_axis(KEY).to_parquet(
        cp.with_name(code + "_patients.parquet")
    )
    cut_frame = pd.DataFrame(cuts.T, columns=[f"p{p}_cutoff" for p in PERCENTILES])
    cut_frame.insert(0, "patient_id", patients)
    cut_frame.insert(0, "cancer_code", code)
    csv(cut_frame, cp.with_name(code + "_cutoffs.csv"))
    meta = {
        "cancer_code": code,
        "n_patients": len(patients),
        "n_gene_background_rows": len(genes),
        "n_proteoform_background_rows": len(background),
        "n_common_background_genes": len(background_ids)
        if background_ids is not None
        else len(genes),
        "n_complete_cta_proteoforms": int(full.sum()),
        "n_partial_cta_proteoforms": int(((np.asarray(member_counts) > 0) & ~full).sum()),
        "mass_conservation_verified": True,
        "source_matrix_sha256": cohort.source_matrix_sha256,
    }
    pd.concat(result, ignore_index=True).to_parquet(cp.with_suffix(".parquet"), index=False)
    cp.with_suffix(".json").write_text(json.dumps(meta, indent=2) + "\n")
    return pd.concat(result, ignore_index=True), meta


def eligible(cohorts, policy):
    minimum, linear = POLICIES[policy]
    return cohorts[
        (cohorts.n_patients >= minimum) & (cohorts.linear_tpm_comparable | ~np.bool_(linear))
    ]


def passing(frame, prevalence, allow_partial=False):
    measured = (
        frame.available_sum_complete_measurement if allow_partial else frame.complete_measurement
    )
    return (
        measured
        & frame.n_patients.ge(MIN_RANKED_PATIENTS)
        & (100 * frame.n_expressing > prevalence * frame.n_patients)
    )


def cancer_type_groups(groups, registry, cohorts):
    """Merge known overlap plus molecular/risk subtype views of one cancer type."""
    lookup = registry.fillna("").set_index("code")
    codes = groups.cancer_code.tolist()
    parents = {c: c for c in codes}

    def root(c):
        while parents[c] != c:
            parents[c] = parents[parents[c]]
            c = parents[c]
        return c

    def type_code(c):
        seen = set()
        while c in lookup.index and c not in seen:
            seen.add(c)
            row = lookup.loc[c]
            if (
                row.ontology_level not in ("molecular_subtype", "evidence_scope")
                or not row.parent_code
            ):
                break
            c = row.parent_code
        return c

    types = {c: type_code(c) for c in codes}
    type_sets = {}
    for c, t in types.items():
        type_sets.setdefault(t, []).append(c)
    for members in [*groups.groupby("overlap_group").cancer_code.agg(list), *type_sets.values()]:
        for c in members[1:]:
            parents[root(c)] = root(members[0])
    sizes = cohorts.set_index("cancer_code").n_patients
    blocks = {}
    for c in codes:
        blocks.setdefault(root(c), []).append(c)
    assignment, members_by_code = {}, {}
    for members in blocks.values():
        anchor = sorted(members, key=lambda c: (-sizes[c], c))[0]
        for c in members:
            assignment[c] = types[anchor]
            members_by_code[c] = ";".join(sorted(members))
    result = groups.copy()
    result["registry_parent_type"] = result.cancer_code.map(types)
    result["cancer_type_group"] = result.cancer_code.map(assignment)
    result["cancer_type_group_members"] = result.cancer_code.map(members_by_code)
    return result


def breadth_rows(data, cohorts, groups, prevalence=10):
    """Fixed eligible-cohort denominators; unavailable measurements stay in denominator."""
    mapping = groups.set_index("cancer_code").overlap_group
    types = groups.set_index("cancer_code").cancer_type_group
    n_cohorts = len(cohorts)
    n_groups = mapping.loc[cohorts.cancer_code].nunique()
    rows = []
    for key, sub in data[data.cancer_code.isin(cohorts.cancer_code)].groupby(KEY):
        good = sub[passing(sub, prevalence)]
        partial = sub[passing(sub, prevalence, allow_partial=True)]
        complete = sub[sub.complete_measurement]
        rows.append(
            {
                KEY: key,
                "n_eligible_cohort_views": n_cohorts,
                "n_eligible_overlap_groups": n_groups,
                "n_eligible_cancer_type_groups": types.loc[cohorts.cancer_code].nunique(),
                "n_fully_measured_cohort_views": len(complete),
                "n_fully_measured_overlap_groups": complete.cancer_code.map(mapping).nunique(),
                "n_cohort_views_gt10": len(good),
                "fraction_cohort_views_gt10": len(good) / n_cohorts,
                "n_overlap_groups_gt10": good.cancer_code.map(mapping).nunique(),
                "fraction_overlap_groups_gt10": good.cancer_code.map(mapping).nunique() / n_groups,
                "n_cancer_type_groups_gt10": good.cancer_code.map(types).nunique(),
                "fraction_cancer_type_groups_gt10": good.cancer_code.map(types).nunique()
                / types.loc[cohorts.cancer_code].nunique(),
                "cancer_type_groups_gt10": ";".join(sorted(good.cancer_code.map(types).unique())),
                "fraction_measured_cohort_views_gt10": len(good) / len(complete)
                if len(complete)
                else np.nan,
                "n_cohort_views_gt10_allow_partial": len(partial),
                "cohort_views_gt10": ";".join(sorted(good.cancer_code)),
                "overlap_groups_gt10": ";".join(sorted(good.cancer_code.map(mapping).unique())),
            }
        )
    return pd.DataFrame(rows)


def summarize(selected, universe, group_map):
    rows = []
    for key, sub in selected.groupby(KEY):
        sub = sub.sort_values(["n_patients", "cancer_code"], ascending=[False, True])
        best = sub.iloc[0]
        linear = sub[sub.linear_tpm_comparable]
        maxn = best.n_patients
        rows.append(
            {
                KEY: key,
                "n_passing_cohort_views": len(sub),
                "n_passing_overlap_groups": sub.cancer_code.map(group_map).nunique(),
                "n_passing_cancer_type_groups": sub.cancer_type_group.nunique(),
                "largest_passing_cohort_size": int(maxn),
                "largest_passing_cohorts": ";".join(
                    sub.loc[sub.n_patients.eq(maxn), "cancer_code"]
                ),
                "largest_cohort_prevalence": best.fraction_expressing,
                "largest_cohort_n_expressing": int(best.n_expressing),
                "largest_cohort_positive_mean": best.mean_expression_expressing,
                "largest_cohort_positive_median": best.median_expression_expressing,
                "largest_cohort_expression_unit": best.expression_unit,
                "largest_cohort_positive_mean_tpm": best.mean_expression_expressing
                if best.linear_tpm_comparable
                else np.nan,
                "largest_linear_passing_cohort": linear.cancer_code.iloc[0] if len(linear) else "",
                "largest_linear_passing_cohort_size": int(linear.n_patients.iloc[0])
                if len(linear)
                else 0,
                "mean_positive_tpm_in_largest_linear_cohort": linear.mean_expression_expressing.iloc[
                    0
                ]
                if len(linear)
                else np.nan,
                "max_cohort_prevalence": sub.fraction_expressing.max(),
                "passing_cohorts": ";".join(sub.cancer_code),
            }
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return universe.iloc[:0].copy()
    frame = universe.merge(frame, on=KEY, validate="one_to_one")
    frame["breadth_rank"] = frame.n_passing_cancer_type_groups.rank(
        method="dense", ascending=False
    ).astype(int)
    return frame.sort_values(
        ["n_passing_cancer_type_groups", "largest_passing_cohort_size", "Symbol"],
        ascending=[False, False, True],
    )


def patient_positive_status(expression, cutoff, complete_members):
    """A missing or partial member sum has no binary positivity verdict."""
    measured = np.isfinite(expression) & np.isfinite(cutoff) & complete_members
    return (expression > cutoff).astype("boolean").where(measured)


def build_exports(out, universe, metrics, cohorts, groups, base):
    if (out / "selections").exists():
        shutil.rmtree(out / "selections")
    for old_gmt in out.glob("proteoform_sets_*.gmt"):
        old_gmt.unlink()
    mapping = groups.set_index("cancer_code").overlap_group
    names = cohorts.set_index("cancer_code").cancer_name
    metrics["cancer_name"] = metrics.cancer_code.map(names)
    metrics["analysis_role"] = metrics.cancer_code.map(
        cohorts.set_index("cancer_code").analysis_role
    )
    metrics["overlap_group"] = metrics.cancer_code.map(mapping)
    metrics["cancer_type_group"] = metrics.cancer_code.map(
        groups.set_index("cancer_code").cancer_type_group
    )
    metrics["mean_expressing_tpm"] = metrics.mean_expression_expressing.where(
        metrics.linear_tpm_comparable
    )
    metrics["median_expressing_tpm"] = metrics.median_expression_expressing.where(
        metrics.linear_tpm_comparable
    )
    csv(metrics, out / "all_proteoform_cohort_metrics.csv.gz")
    coverages = []
    for policy in POLICIES:
        for p in (70, 90):
            breadth = breadth_rows(
                metrics[metrics.transcriptome_percentile.eq(p)], eligible(cohorts, policy), groups
            )
            breadth["cohort_policy"] = policy
            breadth["transcriptome_percentile"] = p
            coverages.append(breadth)
    coverage = pd.concat(coverages, ignore_index=True)
    csv(
        universe.merge(coverage, on=KEY, validate="one_to_many"),
        out / "proteoform_coverage_gt10.csv",
    )
    summary = []
    for policy in POLICIES:
        broad = universe[[KEY]].copy()
        for p in (70, 90):
            b = coverage[
                coverage.cohort_policy.eq(policy) & coverage.transcriptome_percentile.eq(p)
            ].drop(columns=["cohort_policy", "transcriptome_percentile"])
            broad = broad.merge(
                b.rename(columns={c: f"p{p}_{c}" for c in b if c != KEY}),
                on=KEY,
                validate="one_to_one",
            )
        gmts = []
        for f in (50, 70):
            for p in PERCENTILES:
                data = metrics[
                    metrics.transcriptome_percentile.eq(p)
                    & metrics.cancer_code.isin(eligible(cohorts, policy).cancer_code)
                ]
                chosen = data[passing(data, f)]
                ranked = summarize(chosen, universe, mapping).merge(
                    broad, on=KEY, validate="one_to_one"
                )
                scenario = f"prevalence_gt{f}_transcriptome_p{p}"
                dest = out / "selections" / policy / scenario
                csv(ranked, dest / "proteoforms_ranked.csv")
                csv(
                    chosen.merge(universe.drop(columns="Symbol"), on=KEY),
                    dest / "passing_proteoform_cohorts.csv",
                )
                (dest / "proteoform_keys.txt").write_text("\n".join(sorted(ranked[KEY])) + "\n")
                gmts.append(
                    "\t".join(
                        [scenario, policy + " identical-protein entries", *sorted(ranked[KEY])]
                    )
                )
                summary.append(
                    {
                        "cohort_policy": policy,
                        "prevalence_gt_pct": f,
                        "transcriptome_percentile": p,
                        "scenario": scenario,
                        "n_proteoforms": len(ranked),
                        "n_passing_cohort_views": chosen.cancer_code.nunique(),
                        "n_selected_allow_partial": data.loc[
                            passing(data, f, allow_partial=True), KEY
                        ].nunique(),
                        "n_eligible_cohort_views": len(eligible(cohorts, policy)),
                        "n_eligible_overlap_groups": eligible(cohorts, policy)
                        .cancer_code.map(mapping)
                        .nunique(),
                        "n_eligible_cancer_type_groups": eligible(cohorts, policy)
                        .cancer_code.map(groups.set_index("cancer_code").cancer_type_group)
                        .nunique(),
                    }
                )
        (out / f"proteoform_sets_{policy}.gmt").write_text("\n".join(gmts) + "\n")
    summary = pd.DataFrame(summary)
    csv(summary, out / "selection_summary.csv")
    csv(
        summary[
            [
                (r.prevalence_gt_pct, r.transcriptome_percentile) in FOCUSED
                for _, r in summary.iterrows()
            ]
        ],
        out / "focused_selection_summary.csv",
    )
    # The original six genes are now five protein identities; retain the same case study.
    original_six = ["CTAG1A", "CTAG1B", "DPPA5", "EGFL6", "POTEF", "PRAME"]
    original_keys = universe.loc[
        universe.member_symbols.map(lambda x: bool(set(x.split(";")) & set(original_six))), KEY
    ]
    stringent = pd.read_csv(ranked_selection_dir(out, 70, 90) / "proteoforms_ranked.csv")
    cases = set(original_keys) | set(stringent[KEY])
    csv(universe[universe[KEY].isin(cases)], out / "case_study_proteoforms.csv")
    patient_frames = []
    for code in eligible(cohorts, RANKED_COHORT_POLICY).cancer_code:
        cp = out / "checkpoints" / code
        values = pd.read_parquet(cp.with_name(code + "_patients.parquet")).loc[sorted(cases)]
        long = (
            values.rename_axis(KEY)
            .reset_index()
            .melt(id_vars=KEY, var_name="patient_id", value_name="expression")
        )
        cutoffs = pd.read_csv(cp.with_name(code + "_cutoffs.csv"), dtype={"patient_id": str})
        cutoffs = cutoffs[["cancer_code", "patient_id", "p70_cutoff", "p90_cutoff"]]
        long = long.merge(cutoffs, on="patient_id", validate="many_to_one")
        long["expression_measured"] = np.isfinite(long.expression)
        long["complete_member_coverage"] = long[KEY].map(
            metrics[metrics.cancer_code.eq(code) & metrics.transcriptome_percentile.eq(90)]
            .set_index(KEY)
            .complete_member_coverage
        )
        long["positive_p90"] = patient_positive_status(
            long.expression, long.p90_cutoff, long.complete_member_coverage
        )
        long["expression_to_p90_ratio"] = long.expression / long.p90_cutoff
        assert (long.p90_cutoff > 0).all()
        patient_frames.append(long)
    csv(pd.concat(patient_frames, ignore_index=True), out / "case_study_patient_expression.csv.gz")
    # Comparison to prior gene-level focused sets, mapped to protein identities.
    gene_map = pd.read_csv(out / "gene_to_proteoform_mapping.csv").set_index(ID)[KEY]
    changes = []
    old_metrics = pd.read_csv(base / "all_gene_cohort_metrics.csv.gz")
    for f, p in FOCUSED:
        scenario = f"prevalence_gt{f}_transcriptome_p{p}"
        for policy in POLICIES:
            old = old_metrics[
                old_metrics.transcriptome_percentile.eq(p)
                & old_metrics.cancer_code.isin(eligible(cohorts, policy).cancer_code)
                & old_metrics.complete_measurement
                & (100 * old_metrics.n_expressing > f * old_metrics.n_patients)
            ][[ID, "Symbol"]].drop_duplicates()
            old_keys = set(old[ID].map(gene_map))
            new = pd.read_csv(out / "selections" / policy / scenario / "proteoforms_ranked.csv")
            for key in sorted(old_keys | set(new[KEY])):
                changes.append(
                    {
                        KEY: key,
                        "scenario": scenario,
                        "cohort_policy": policy,
                        "selected_as_any_gene_before": key in old_keys,
                        "selected_as_proteoform_now": key in set(new[KEY]),
                        "old_selected_member_symbols": ";".join(
                            old.loc[old[ID].map(gene_map).eq(key), "Symbol"]
                        ),
                    }
                )
    csv(universe.merge(pd.DataFrame(changes), on=KEY), out / "gene_vs_proteoform_selection.csv")
    return summary, coverage


def validate(out, universe, metrics, summary, groups):
    assert len(metrics) == len(universe) * metrics.cancer_code.nunique() * len(PERCENTILES)
    assert not metrics.duplicated([KEY, "cancer_code", "transcriptome_percentile"]).any()
    assert (metrics.n_expressing <= metrics.n_patients).all()
    assert np.allclose(metrics.fraction_expressing, metrics.n_expressing / metrics.n_patients)
    assert metrics.loc[metrics.n_expressing.eq(0), "mean_expression_expressing"].isna().all()
    sets = {}
    for _, s in summary.iterrows():
        folder = out / "selections" / s.cohort_policy / s.scenario
        ranked = pd.read_csv(folder / "proteoforms_ranked.csv")
        passed = pd.read_csv(folder / "passing_proteoform_cohorts.csv")
        assert len(ranked) == s.n_proteoforms
        assert set(ranked[KEY]) == set(passed[KEY])
        assert passed.complete_measurement.all()
        assert passed.n_patients.ge(MIN_RANKED_PATIENTS).all()
        assert passing(passed, s.prevalence_gt_pct).all()
        minimum, linear = POLICIES[s.cohort_policy]
        assert passed.n_patients.ge(minimum).all()
        if linear:
            assert passed.linear_tpm_comparable.all()
        for _, r in ranked.iterrows():
            sub = passed[passed[KEY].eq(r[KEY])]
            assert sub.n_patients.max() == r.largest_passing_cohort_size
            assert sub.overlap_group.nunique() == r.n_passing_overlap_groups
        sets[s.cohort_policy, int(s.prevalence_gt_pct), int(s.transcriptome_percentile)] = set(
            ranked[KEY]
        )
    for policy in POLICIES:
        for p in PERCENTILES:
            for lo, hi in [(50, 70)]:
                assert sets[policy, hi, p] <= sets[policy, lo, p]
        for f in (50, 70):
            for lo, hi in zip(PERCENTILES, PERCENTILES[1:]):
                assert sets[policy, f, hi] <= sets[policy, f, lo]
    patient_values = pd.read_csv(out / "case_study_patient_expression.csv.gz")
    indexed = metrics[metrics.transcriptome_percentile.eq(90)].set_index([KEY, "cancer_code"])
    n_patient_checks = 0
    for (key, code), sub in patient_values.groupby([KEY, "cancer_code"]):
        truth = indexed.loc[(key, code)]
        positive = sub[sub.positive_p90.fillna(False).astype(bool)]
        assert len(sub) == truth.n_patients
        if not truth.complete_measurement:
            assert sub.positive_p90.isna().all()
            continue
        assert len(positive) == truth.n_expressing
        if len(positive):
            assert np.isclose(
                positive.expression.mean(), truth.mean_expression_expressing, rtol=1e-12
            )
        n_patient_checks += 1
    assert metrics.loc[~metrics.linear_tpm_comparable, "mean_expressing_tpm"].isna().all()
    validation = {
        "status": "passed",
        "n_proteoforms": len(universe),
        "n_proteoform_set_exports": len(summary),
        "n_patient_distribution_cohort_checks": n_patient_checks,
        "checks": [
            "exact longest-protein sequence identity",
            "per-patient TPM mass conservation in all 142 cohorts",
            "source matrix hashes match original analysis",
            "same QC and patient groups",
            "strict prevalence boundaries",
            "complete member coverage for selection",
            "nested thresholds and cohort policies",
            "largest passing cohort and overlap-adjusted breadth",
        ],
        "n_overlap_groups": groups.overlap_group.nunique(),
    }
    (out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / "outputs/cta_threshold_report_20260916")
    parser.add_argument("--out", type=Path, default=ROOT / "outputs/cta_proteoform_report_20260917")
    parser.add_argument(
        "--source-cache", type=Path, default=ROOT / "tmp/cta_threshold_report/source-matrices"
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--fasta",
        type=Path,
        help="Checksum-pinned Ensembl 112 protein FASTA; defaults to the pyensembl macOS cache",
    )
    args = parser.parse_args()
    os.environ[CACHE_DIR_ENV_VAR] = str(args.source_cache.resolve())
    out = args.out.resolve()
    (out / "checkpoints").mkdir(parents=True, exist_ok=True)
    script_hash = sha(__file__)
    verify_analysis(args.base)
    universe = universe_tables(args.base, out, fasta=args.fasta)
    background_ids = pd.read_csv(args.base / "common_background.csv")[ID].tolist()
    policy = fingerprint(
        {
            "implementation": implementation_hash(),
            "background": background_ids,
            "universe": universe.to_dict("records"),
            "sample_audit": sha(args.base / "sample_patient_audit.csv.gz"),
        }
    )
    audit = pd.read_csv(args.base / "sample_patient_audit.csv.gz")
    cohorts = pd.read_csv(args.base / "cohort_audit.csv")
    registry = pd.read_csv(ROOT / "oncoref/data/cancer-type-registry.csv")
    cohorts["cancer_name"] = cohorts.cancer_code.map(registry.set_index("code")["name"])
    from cta_threshold_report_details import overlap_groups

    groups, edges = overlap_groups(audit, cohorts)
    groups = cancer_type_groups(groups, registry, cohorts)
    csv(groups, out / "cohort_overlap_groups.csv")
    csv(edges, out / "cohort_overlap_edges.csv")
    frames, metas = [], []
    for i, (_, cohort) in enumerate(cohorts.iterrows(), 1):
        code = cohort.cancer_code
        cp = out / "checkpoints" / code
        expected = fingerprint([policy, code, sha(local_path(code))])
        suffixes = [".parquet", "_patients.parquet", "_cutoffs.csv"]
        if args.resume and checkpoint_valid(cp, expected, suffixes):
            data = pd.read_parquet(cp.with_suffix(".parquet"))
            meta = checkpoint_payload(cp)
        else:
            data, meta = analyze_cohort(
                code, universe, audit[audit.cancer_code.eq(code)], cohort, out, background_ids
            )
        write_checkpoint(cp, meta, expected, suffixes)
        frames.append(data)
        metas.append(meta)
        print(
            f"[{i}/{len(cohorts)}] {code}: n={cohort.n_patients}, complete proteoforms={meta['n_complete_cta_proteoforms']}",
            flush=True,
        )
    metrics = pd.concat(frames, ignore_index=True)
    metrics = metrics[metrics.transcriptome_percentile.isin(PERCENTILES)].reset_index(drop=True)
    cohorts = cohorts.merge(
        pd.DataFrame(metas).drop(
            columns=["n_patients", "source_matrix_sha256", "n_common_background_genes"]
        ),
        on="cancer_code",
        validate="one_to_one",
    )
    cohorts["analysis_role"] = np.where(
        cohorts.n_patients.ge(MIN_RANKED_PATIENTS), "ranked", "exploratory_only"
    )
    csv(cohorts, out / "cohort_audit.csv")
    csv(
        cohorts.loc[
            cohorts.analysis_role.eq("exploratory_only"),
            ["cancer_code", "cancer_name", "n_patients", "analysis_role"],
        ],
        out / "exploratory_cohort_audit.csv",
    )
    cutoffs = pd.concat(
        [pd.read_csv(out / "checkpoints" / f"{code}_cutoffs.csv") for code in cohorts.cancer_code],
        ignore_index=True,
    )
    csv(
        cutoffs[["cancer_code", "patient_id", "p70_cutoff", "p90_cutoff"]],
        out / "patient_transcriptome_cutoffs.csv.gz",
    )
    summary, _ = build_exports(out, universe, metrics, cohorts, groups, args.base)
    validate(out, universe, metrics, summary, groups)
    shutil.copy2(args.base / "sample_patient_audit.csv.gz", out / "sample_patient_audit.csv.gz")
    if sha(__file__) != script_hash:
        raise RuntimeError("Analysis code changed during this run; rebuild before publishing")
    (out / "run_manifest.json").write_text(
        json.dumps(
            {
                "common_background_genes": len(background_ids),
                "common_background_sha256": sha(args.base / "common_background.csv"),
                "oncoref_version": __version__,
                "base_gene_report": str(args.base.resolve()),
                "expression_level": "proteoform",
                "collapse_scope": "genome",
                "percentiles": PERCENTILES,
                "prevalence_cutoffs": [50, 70],
                "focused_combinations": FOCUSED,
                "coverage_prevalence_gt_percent": 10,
                "cohort_policies": POLICIES,
                "minimum_patients_for_selection": MIN_RANKED_PATIENTS,
                "small_cohort_policy": f"Fewer than {MIN_RANKED_PATIENTS} patients: exploratory only, excluded from rankings and coverage denominators; retained in all-cohort metrics and the exploratory audit.",
                "partial_member_policy": "Primary selection requires every registered member locus measured; available-member sums are exported as a sensitivity analysis.",
                "coverage_denominator": "All eligible cohort views or known-overlap components, including unavailable entries; also report evaluable-only view fraction.",
                "script_sha256": script_hash,
                "base_sample_audit_sha256": sha(args.base / "sample_patient_audit.csv.gz"),
            },
            indent=2,
        )
        + "\n"
    )
    seal_stage(
        out,
        "analysis",
        [
            Path(__file__),
            ROOT / "scripts/cta_report_common.py",
            args.base / "analysis_receipt.json",
            args.base / "common_background.csv",
            ROOT / "scripts/cta_report_protein_lengths.py",
            Path(json.loads((out / "annotation_manifest.json").read_text())["protein_fasta_path"]),
            *json.loads((args.base / "analysis_receipt.json").read_text())["inputs"],
            *json.loads((args.base / "analysis_receipt.json").read_text())["outputs"],
        ],
        [
            out / n
            for n in [
                "validation.json",
                "run_manifest.json",
                "annotation_manifest.json",
                "proteoform_universe.csv",
                "all_proteoform_cohort_metrics.csv.gz",
                "cohort_audit.csv",
                "exploratory_cohort_audit.csv",
                "cohort_overlap_groups.csv",
                "selection_summary.csv",
                "proteoform_member_annotations.csv",
                "member_protein_isoforms.csv",
                "gene_to_proteoform_mapping.csv",
                "genome_identical_protein_registry.csv",
                "cta_identical_protein_registry.csv",
                "proteoform_coverage_gt10.csv",
                "case_study_proteoforms.csv",
                "case_study_patient_expression.csv.gz",
                "gene_vs_proteoform_selection.csv",
                "sample_patient_audit.csv.gz",
                "patient_transcriptome_cutoffs.csv.gz",
                "cohort_overlap_edges.csv",
            ]
        ]
        + sorted((out / "selections").rglob("*.csv"))
        + sorted((out / "checkpoints").glob("*_patients.parquet"))
        + sorted((out / "checkpoints").glob("*_cutoffs.csv")),
    )
    print(
        summary[
            (summary.prevalence_gt_pct.isin([50, 70]))
            & summary.transcriptome_percentile.isin([70, 90])
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()

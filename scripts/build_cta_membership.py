#!/usr/bin/env python3
"""Reprocess canonical genes into a source-backed biological CTA/oncofetal set."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import shutil
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import (  # noqa: E402
    cta_membership as policy,
)
from oncoref import (  # noqa: E402
    cta_policy,
    cta_provenance,
    cta_research_policy,
    cta_tissues,
    hpa,
    hpa_cancer,
    normal_expression_breadth,
    reference_data,
)

CORE = {"testis", "ovary", "placenta"}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(frame, path):
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "mtime": 0} if path.name.endswith(".gz") else None,
    )


def join(values):
    return "|".join(sorted(set(values.dropna().astype(str)) - {""}))


def expression_facts(base, normal, membership):
    """Keep tumor pTPM and normal nTPM as separate measurement systems."""
    d = base.rename(
        columns={
            "Target": "Previous member",
            "Category": "Previous category",
            "Detail": "Previous selection detail",
            "Reason": "Previous selection reason",
            "Policy": "Previous selection policy",
        }
    ).copy()
    normal = normal.loc[normal.Gene.isin(d["Gene ID"])].copy()
    if normal.duplicated(["Gene", "Tissue"]).any():
        raise ValueError("Resolve duplicate normal observations before selection")
    if not np.isfinite(normal.nTPM).all() or normal.nTPM.lt(0).any():
        raise ValueError("Invalid HPA consensus measurements")
    normal = normal.pivot(index="Gene", columns="Tissue", values="nTPM")
    if len(normal.columns) != 50 or not set(normal.columns) >= CORE:
        raise ValueError("HPA v23 panel changed; review the named tissue universe")
    core = normal[sorted(CORE)]
    somatic = normal.drop(columns=sorted(CORE))
    other_reproductive = sorted(cta_tissues.EXTENDED_REPRODUCTIVE_TISSUES - CORE)
    nonreproductive = normal.drop(columns=sorted(cta_tissues.EXTENDED_REPRODUCTIVE_TISSUES))
    for label, series in {
        "Core RNA (nTPM)": core.max(axis=1),
        "Core tissues measured": core.notna().sum(axis=1),
        "Core comparator maximum (nTPM)": somatic.max(axis=1),
        "Core comparator tissues measured": somatic.notna().sum(axis=1),
        "Testis RNA (nTPM)": normal.testis,
        "Ovary RNA (nTPM)": normal.ovary,
        "Placenta RNA (nTPM)": normal.placenta,
        "Other reproductive RNA (nTPM)": normal[other_reproductive].max(axis=1),
        "Other reproductive tissues measured": normal[other_reproductive].notna().sum(axis=1),
        "Nonreproductive maximum (nTPM)": nonreproductive.max(axis=1),
        "Nonreproductive tissues measured": nonreproductive.notna().sum(axis=1),
    }.items():
        d[label] = d["Gene ID"].map(series)
    for label in [
        "Core tissues measured",
        "Core comparator tissues measured",
        "Other reproductive tissues measured",
        "Nonreproductive tissues measured",
    ]:
        d[label] = d[label].fillna(0).astype(int)
    # Enrichment is a positive class-association measure, not a safety score.
    # Other reproductive organs remain in this deliberately conservative
    # 47-tissue comparator; the independent anatomical breadth model is unchanged.
    d["Core / comparator maximum"] = d["Core RNA (nTPM)"] / d[
        "Core comparator maximum (nTPM)"
    ].replace(0, np.nan)
    d["Core comparator maximum is zero"] = d["Core comparator maximum (nTPM)"].eq(0)
    for category in ["Cancer-germline", "Cancer-placental"]:
        eligible = membership.loc[
            membership["Counts toward consensus"] & membership.Class.eq(category)
        ]
        counts = eligible.groupby("Ensembl_Gene_ID").doi.nunique()
        d[category + " papers"] = d["Gene ID"].map(counts).fillna(0).astype(int)
    d["Publication classes"] = (
        d["Gene ID"].map(membership.groupby("Ensembl_Gene_ID").Class.agg(join)).fillna("")
    )
    d["Nomination DOIs"] = (
        d["Gene ID"].map(membership.groupby("Ensembl_Gene_ID").doi.agg(join)).fillna("")
    )
    d["Nomination source tags"] = (
        d["Gene ID"].map(membership.groupby("Ensembl_Gene_ID").source_tag.agg(join)).fillna("")
    )
    return d


def tumor_facts(facts, rna, ihc, thresholds, *, rna_floor=1):
    d = facts.copy()
    suffix = str(rna_floor).replace(".", "_")
    count = "expressed_samples_ptpm_ge_" + suffix
    prevalence = "prevalence_ptpm_ge_" + suffix
    rna = rna.copy()
    rna["Positive"] = rna[count].ge(thresholds.tumor_samples) & rna[prevalence].ge(
        thresholds.tumor_prevalence
    )
    rna["Cohort key"] = rna.cancer_code + ":" + rna.cohort
    ihc = ihc.copy()
    ihc["Positive"] = ihc.medium_high.ge(thresholds.tumor_ihc_cases) & ihc.measurement_status.eq(
        "measured"
    )
    for label, series in {
        "Tumor RNA positive": rna.groupby("gene_id").Positive.any(),
        "Tumor IHC positive": ihc.groupby("gene_id").Positive.any(),
    }.items():
        d[label] = d["Gene ID"].map(series).eq(True).fillna(False).astype(bool)
    d["Tumor RNA cohorts measured"] = (
        d["Gene ID"].map(rna.loc[rna.samples.gt(0)].groupby("gene_id").size()).fillna(0).astype(int)
    )
    d["Tumor IHC groups measured"] = (
        d["Gene ID"]
        .map(
            ihc.loc[ihc.measurement_status.eq("measured") & ihc.total.gt(0)]
            .groupby("gene_id")
            .size()
        )
        .fillna(0)
        .astype(int)
    )
    d["Positive RNA cohorts"] = (
        d["Gene ID"]
        .map(rna.loc[rna.Positive].groupby("gene_id")["Cohort key"].agg(join))
        .fillna("")
    )
    d["Positive IHC groups"] = (
        d["Gene ID"].map(ihc.loc[ihc.Positive].groupby("gene_id").cancer.agg(join)).fillna("")
    )
    d["Tumor RNA support only TGCT"] = d["Positive RNA cohorts"].eq("TGCT:TCGA")
    d["Tumor RNA threshold (pTPM)"] = rna_floor
    d["Tumor minimum samples"] = thresholds.tumor_samples
    d["Tumor minimum prevalence"] = thresholds.tumor_prevalence
    d["Tumor IHC minimum cases"] = thresholds.tumor_ihc_cases
    return d, rna, ihc


def breadth(normal, gene_ids):
    """The existing fixed anatomical score for every processed gene, with gaps."""
    panel = normal_expression_breadth.anatomical_panel()
    profiles = pd.MultiIndex.from_product(
        [sorted(gene_ids), panel.compartment], names=["gene_id", "compartment"]
    ).to_frame(index=False)
    source = normal.rename(columns={"Gene": "gene_id", "Tissue": "compartment", "nTPM": "value"})
    profiles = profiles.merge(
        source[["gene_id", "compartment", "value"]],
        how="left",
        on=["gene_id", "compartment"],
        validate="one_to_one",
    )
    ceilings = {key: values[0] for key, values in cta_policy.manifest()["rna_ceiling_ntpm"].items()}
    profiles["ceiling"] = profiles.compartment.map(panel.set_index("compartment").category).map(
        ceilings
    )
    profiles["source"], profiles["version"], profiles["taxonomy_id"] = (
        "hpa_rna_consensus",
        "v23",
        9606,
    )
    profiles["unit"], profiles["cohort"], profiles["summary"] = "nTPM", "HPA consensus", "consensus"
    _, scores = normal_expression_breadth.burden(profiles, panel)
    return scores


def nonselection_reasons(decisions, thresholds):
    """Explain failure of positive membership evidence, separately from risk."""
    d = decisions.loc[~decisions.Member].copy()
    publication = d[["Cancer-germline papers", "Cancer-placental papers"]].max(axis=1).ge(1)
    tumor = d["Tumor RNA positive"] | d["Tumor IHC positive"] | publication
    core_coverage = d["Core tissues measured"].eq(3) & d["Core comparator tissues measured"].eq(47)
    other_coverage = (
        d["Other reproductive tissues measured"].eq(8)
        & d["Nonreproductive tissues measured"].eq(39)
        & publication
    )
    core_floor = d["Core RNA (nTPM)"].ge(thresholds.reproductive_ntpm)
    other_floor = d["Other reproductive RNA (nTPM)"].ge(thresholds.reproductive_ntpm)
    d["Nonselection reason"] = np.select(
        [
            ~tumor,
            ~(core_coverage | other_coverage),
            ~((core_coverage & core_floor) | (other_coverage & other_floor)),
        ],
        [
            "No qualifying published or measured tumor support",
            "Incomplete required reproductive/comparator RNA coverage",
            "Reproductive RNA below admission floor in eligible measured scope",
        ],
        default="Reproductive enrichment below required fold in eligible measured scope",
    )
    columns = [
        "Gene ID",
        "Gene",
        "Nominated",
        "Nonselection reason",
        "Expression category",
        "Expression detail",
        "Cancer-germline papers",
        "Cancer-placental papers",
        "Core RNA (nTPM)",
        "Core comparator maximum (nTPM)",
        "Core tissues measured",
        "Core comparator tissues measured",
        "Other reproductive RNA (nTPM)",
        "Nonreproductive maximum (nTPM)",
        "Other reproductive tissues measured",
        "Nonreproductive tissues measured",
        "Tumor RNA positive",
        "Tumor IHC positive",
        "Reviewed admission role",
    ]
    return d[columns]


def build(out, *, publish=False):
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    sources = []
    for name, version in [
        ("hpa_rna_consensus", "v23"),
        ("hpa_cancer_rna", "v25.1-1"),
        ("hpa_cancer_ihc", "v25.1-1"),
    ]:
        record = reference_data.provenance(name, version, verify_content=True)
        if record["verification_state"] != "verified":
            raise ValueError(f"Unverified source {name}: {record['verification_state']}")
        sources.append({k: record[k] for k in ["name", "version", "url", "bytes", "sha256"]})
    base = cta_research_policy.results()
    roles = pd.read_csv(policy.DATA / "cta-membership-source-roles.csv", keep_default_na=False)
    membership = cta_provenance.paper_membership().merge(
        roles, on="source_tag", how="left", validate="many_to_one"
    )
    if membership["Counts toward consensus"].isna().any():
        raise ValueError("Unreviewed nomination source role")
    membership = membership.loc[
        membership.protein_candidate_eligible & membership.Ensembl_Gene_ID.isin(base["Gene ID"])
    ]
    normal = hpa._read_hpa("hpa_rna_consensus", "v23")
    rna = hpa_cancer.hpa_cancer_rna_prevalence()
    ihc = hpa_cancer.hpa_cancer_ihc_prevalence()
    facts = expression_facts(base, normal, membership)
    thresholds = policy.Thresholds()
    facts, rna_assessed, ihc_assessed = tumor_facts(facts, rna, ihc, thresholds)
    reviewed = policy.evidence()
    decisions = policy.classify(facts, reviewed, thresholds=thresholds)
    decisions["Change from 215"] = np.select(
        [
            decisions.Member & ~decisions["Previous member"],
            ~decisions.Member & decisions["Previous member"],
        ],
        ["Added", "Removed"],
        default="Unchanged",
    )
    cell = pd.read_csv(ROOT / "docs/audits/cta-cell-thresholds-20261026/gene-assessments.csv")
    write(cell, out / "prior-cell-rna-assessments.csv")
    for field in ["cta_excluded", "reference_review", "no_source_records"]:
        decisions["Prior cell audit: " + field] = decisions["Gene ID"].map(
            cell.set_index("Gene ID")[field]
        )
    decisions["Cell RNA audit coverage"] = np.where(
        decisions["Gene ID"].isin(cell["Gene ID"]),
        "Prior tissue/IHC audit; see linked rows",
        "Not assessed in prior cell-RNA audit",
    )
    decisions["Risk assessment"] = "Retained expression findings; not clinical safety clearance"
    flags = [
        "Safety veto",
        "Normal IHC positive",
        "Normal RNA failure",
        "High normal RNA",
        "Variable",
        "Antibody concern",
        "Cell RNA concern",
        "Accepted normal protein",
        "Unresolved normal evidence",
        "Expression forbidden",
    ]
    decisions["Normal expression flags"] = decisions.apply(
        lambda row: "|".join(field for field in flags if row[field]), axis=1
    )
    decisions["Developmental atlas coverage"] = (
        "Not imported; oncofetal admission uses cited gene-specific evidence"
    )
    donor_path = ROOT / "docs/audits/normal-breadth-20261026/donor-profiles.csv.gz"
    donor_genes = set(pd.read_csv(donor_path).gene_id)
    decisions["Donor uncertainty"] = np.where(
        decisions["Gene ID"].isin(donor_genes),
        "Limited imported donors; see separate donor summaries and intervals",
        "Individual donors unassessed; consensus RNA is not a donor uncertainty estimate",
    )
    scores = breadth(normal, set(decisions["Gene ID"]))
    write(scores, out / "normal-breadth.csv.gz")
    for category in sorted(normal_expression_breadth.TIERS):
        subset = scores.loc[scores.category.eq(category)].set_index("gene_id")
        for value in ["lower", "upper", "missing_weight"]:
            decisions[f"Normal breadth {category} {value}"] = decisions["Gene ID"].map(
                subset[value]
            )
    exclusions = nonselection_reasons(decisions, thresholds)
    write(exclusions, out / "nonselection-reasons.csv.gz")
    reason_counts = []
    for scope, rows in [
        ("All canonical nonmembers", exclusions),
        ("Previously nominated nonmembers", exclusions.loc[exclusions.Nominated]),
    ]:
        reason_counts.extend(
            {"Scope": scope, "Reason": reason, "Genes": int(count)}
            for reason, count in rows["Nonselection reason"].value_counts().items()
        )
    write(pd.DataFrame(reason_counts), out / "nonselection-counts.csv")
    risk_counts = []
    for scope, rows in [
        ("All canonical genes", decisions),
        ("Prior nominated genes", decisions.loc[decisions.Nominated]),
        ("Selected biological members", decisions.loc[decisions.Member]),
        ("Unselected canonical genes", decisions.loc[~decisions.Member]),
        ("Selected prior nominees", decisions.loc[decisions.Nominated & decisions.Member]),
        ("Unselected prior nominees", decisions.loc[decisions.Nominated & ~decisions.Member]),
    ]:
        risk_counts.append(
            {
                "Scope": scope,
                "Genes": len(rows),
                **{
                    category: int(rows["Expression category"].eq(category).sum())
                    for category in ["CTA", "Unsafe", "Neither"]
                },
            }
        )
    write(pd.DataFrame(risk_counts), out / "expression-status-counts.csv")
    named = [
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
        "MAGEA1",
        "MAGEA3",
        "MAGEA6",
        "MAGEA10",
        "MAGEC1",
        "MAGEC2",
    ]
    write(
        decisions.loc[decisions.Gene.isin(named)].sort_values("Gene"),
        out / "named-gene-results.csv",
    )
    write(facts, out / "decision-facts.csv.gz")
    write(decisions, out / "gene-results.csv.gz")
    write(decisions.loc[decisions.Member], out / "cta-genes.csv")
    write(decisions.loc[decisions["Change from 215"].ne("Unchanged")], out / "changes.csv")
    write(reviewed, out / "reviewed-membership-evidence.csv")
    write(membership, out / "nomination-membership.csv.gz")
    write(
        rna_assessed.loc[rna_assessed.gene_id.isin(decisions.loc[decisions.Member, "Gene ID"])],
        out / "tumor-rna-evidence.csv.gz",
    )
    write(
        ihc_assessed.loc[ihc_assessed.gene_id.isin(decisions.loc[decisions.Member, "Gene ID"])],
        out / "tumor-ihc-evidence.csv.gz",
    )
    benchmark = policy.controls().merge(
        decisions.drop(columns="Gene"), on="Gene ID", validate="one_to_one"
    )
    benchmark["Benchmark passed"] = benchmark.Member.eq(benchmark["Expected member"])
    write(benchmark, out / "control-results.csv")
    if not benchmark["Benchmark passed"].all():
        raise ValueError("Biological controls failed; inspect control-results.csv")
    sensitivity = []
    # One-axis variations make the admission effect legible; these controls are
    # calibration benchmarks, not held-out validation or a sensitivity estimate.
    variants = [("baseline", thresholds, 1)]
    for field, values in {
        "reproductive_ntpm": [0.1, 2.0, 5.0],
        "reproductive_fold": [2.0, 8.0],
        "tumor_prevalence": [0.01, 0.1],
        "tumor_samples": [1, 5],
        "tumor_ihc_cases": [1, 3],
        "papers": [1, 3],
    }.items():
        variants.extend(
            (f"{field}={value}", policy.Thresholds(**{**asdict(thresholds), field: value}), 1)
            for value in values
        )
    variants.extend((f"tumor_rna_floor={value}", thresholds, value) for value in [0.1, 2, 5])
    for label, variant, floor in variants:
        inputs, _, _ = tumor_facts(facts, rna, ihc, variant, rna_floor=floor)
        result = policy.classify(inputs, reviewed, thresholds=variant)
        controls = policy.controls().merge(
            result[["Gene ID", "Member"]], on="Gene ID", validate="one_to_one"
        )
        sensitivity.append(
            {
                "Scenario": label,
                **asdict(variant),
                "tumor_rna_floor_ptpm": floor,
                "Selected": int(result.Member.sum()),
                "Positive controls selected": int(
                    controls.loc[controls["Expected member"], "Member"].sum()
                ),
                "Negative controls selected": int(
                    controls.loc[~controls["Expected member"], "Member"].sum()
                ),
            }
        )
    write(pd.DataFrame(sensitivity), out / "threshold-sensitivity.csv")
    source_sets = {doi: set(group.Ensembl_Gene_ID) for doi, group in membership.groupby("doi")}
    selected = set(decisions.loc[decisions.Member, "Gene ID"])
    write(
        pd.DataFrame(
            [
                {"DOI": doi, "Mapped nominees": len(ids), "Selected": len(ids & selected)}
                for doi, ids in source_sets.items()
            ]
        ),
        out / "source-counts.csv",
    )
    write(
        pd.DataFrame(
            [
                {
                    "Source A": a,
                    "Source B": b,
                    "Overlap": len(source_sets[a] & source_sets[b]),
                    "Selected overlap": len(source_sets[a] & source_sets[b] & selected),
                }
                for a, b in itertools.product(source_sets, repeat=2)
            ]
        ),
        out / "source-overlaps.csv",
    )
    summary = {
        "policy": policy.POLICY_ID,
        "canonical_genes_processed": len(decisions),
        "prior_nominated_genes": int(base.Nominated.sum()),
        "prior_selection": int(base.Target.sum()),
        "selected_genes": int(decisions.Member.sum()),
        "selected_evidence_tiers": decisions.loc[decisions.Member, "Evidence tier"]
        .value_counts()
        .to_dict(),
        "selected_classes_nonexclusive": {
            c: int(
                decisions.loc[decisions.Member, "Classes"]
                .str.split("|")
                .map(lambda x, c=c: c in x)
                .sum()
            )
            for c in sorted(policy.CLASSES)
        },
        "changes": decisions["Change from 215"].value_counts().to_dict(),
        "selected_outside_prior_nomination": int((decisions.Member & ~decisions.Nominated).sum()),
        "selected_with_expression_forbidden_flag": int(
            (decisions.Member & decisions["Expression forbidden"]).sum()
        ),
        "selected_with_only_tgct_rna_support": int(
            (decisions.Member & decisions["Tumor RNA support only TGCT"]).sum()
        ),
        "positive_controls": int(benchmark["Expected member"].sum()),
        "negative_controls": int((~benchmark["Expected member"]).sum()),
        "all_controls_pass": bool(benchmark["Benchmark passed"].all()),
        "context_only_reviewed_genes": sorted(
            reviewed.loc[reviewed["Admission role"].eq("Context only"), "Gene"].unique()
        ),
        "expression_status_counts": risk_counts,
        "nonselection_counts": reason_counts,
        "thresholds": {**asdict(thresholds), "tumor_rna_floor_ptpm": 1},
        "calibration": "User-requested biological benchmarks; no independent validation or safety calibration",
        "normal_risk_threshold_changes": 0,
        "forbidden_sequence_reference_changes": 0,
        "developmental_atlas_measurements_imported": False,
        "comprehensive_individual_literature_review_complete": False,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    inputs = [
        "oncoref/cta_membership.py",
        "scripts/build_cta_membership.py",
        "oncoref/cta_provenance.py",
        "oncoref/cta_sources.py",
        "oncoref/hpa_cancer.py",
        "oncoref/normal_expression_breadth.py",
        "oncoref/cta_tissues.py",
        "oncoref/data/cta-policy.json",
        "oncoref/data/cta-membership-evidence.csv",
        "oncoref/data/cta-membership-controls.csv",
        "oncoref/data/cta-membership-paper-receipts.json",
        "oncoref/data/cta-membership-source-roles.csv",
        "oncoref/data/cta-research-gene-sets.csv.gz",
        "oncoref/data/cta-publication-membership.csv",
        "oncoref/data/cta-publication-sources.csv",
        "docs/audits/cta-cell-thresholds-20261026/gene-assessments.csv",
        "docs/audits/normal-breadth-20261026/donor-profiles.csv.gz",
    ]
    manifest = {
        "policy": policy.POLICY_ID,
        "inputs": {name: digest(ROOT / name) for name in inputs},
        "external_sources": sources,
        "outputs": {p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file()},
    }
    (out / "run-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if publish:
        shutil.copyfile(
            out / "gene-results.csv.gz", policy.DATA / "cta-membership-gene-sets.csv.gz"
        )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--publish", action="store_true")
    args = parser.parse_args()
    build(args.output, publish=args.publish)

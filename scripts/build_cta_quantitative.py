#!/usr/bin/env python3
"""Rebuild a quantitative research candidate set without literature admission."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import shutil
import sys
from dataclasses import asdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import cta_membership, cta_tissues, hpa, hpa_cancer, reference_data  # noqa: E402
from oncoref import cta_quantitative as policy  # noqa: E402


def digest(path):
    """Return the content SHA-256 used to pin an audit input or output."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(frame: pd.DataFrame, path):
    """Write a data table with deterministic compression when requested."""
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "mtime": 0} if path.name.endswith(".gz") else None,
    )


def build(prepared, out: Path, *, publish=False) -> None:
    """Build the versioned audit tables, report and source provenance."""
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    out.mkdir(parents=True)
    inputs = {}

    def read(name):
        path = prepared / name
        inputs[name] = digest(path)
        return pd.read_csv(path, low_memory=False)

    # The prior table supplies the fixed identity universe and OUTPUT annotations only.
    old = cta_membership.results()
    facts = old[["Gene ID", "Gene", "Taxonomy ID", "Reference"]].rename(
        columns={
            "Gene ID": "gene_id",
            "Gene": "symbol",
            "Taxonomy ID": "taxonomy_id",
            "Reference": "reference",
        }
    )
    normal = hpa._read_hpa("hpa_rna_consensus", "v23")
    normal = normal[normal.Gene.isin(facts.gene_id)].pivot(
        index="Gene", columns="Tissue", values="nTPM"
    )
    if len(normal.columns) != 50:
        raise ValueError("Expected pinned 50-tissue HPA bulk panel")
    core = normal[sorted(cta_tissues.CORE_REPRODUCTIVE_TISSUES)]
    adult = normal.drop(columns=sorted(cta_tissues.EXTENDED_REPRODUCTIVE_TISSUES))
    for name, s in {
        "reproductive_ntpm": core.max(axis=1),
        "adult_bulk_max_ntpm": adult.max(axis=1),
        "normal_bulk_tissues": normal.notna().sum(axis=1),
    }.items():
        facts[name] = facts.gene_id.map(s)
    dev = read("developmental-maxima.csv").set_index("Gene ID")
    for col, name in {
        "Lineage": "developmental_lineage",
        "Mean TPM": "developmental_mean_tpm",
        "Cells": "developmental_cells",
        "Positive cells TPM>=1": "developmental_positive_cells",
        "Fraction TPM>=1": "developmental_fraction",
        "Embryos": "developmental_embryos",
    }.items():
        facts[name] = facts.gene_id.map(dev[col])
    lineages = read("developmental-lineages.csv.gz")
    facts["developmental_gene_measured"] = facts.gene_id.isin(lineages["Gene ID"])
    tumor = hpa_cancer.hpa_cancer_rna_prevalence()
    tumor["qualifies"] = tumor.expressed_samples_ptpm_ge_1.ge(3) & tumor.prevalence_ptpm_ge_1.ge(
        0.05
    )
    facts["tumor_positive"] = facts.gene_id.map(tumor.groupby("gene_id").qualifies.any()).eq(True)
    facts["tumor_cohorts_measured"] = facts.gene_id.map(
        tumor[tumor.samples.gt(0)].groupby("gene_id").size()
    )
    facts["tumor_cohorts_positive"] = facts.gene_id.map(
        tumor[tumor.qualifies]
        .groupby("gene_id")
        .cancer_code.agg(lambda s: "|".join(sorted(set(s))))
    )
    facts["tumor_floor_ptpm"], facts["tumor_minimum_samples"], facts["tumor_minimum_prevalence"] = (
        1,
        3,
        0.05,
    )
    cells = (
        read("adult-cell-summary.csv").drop(columns="Gene").rename(columns={"Gene ID": "gene_id"})
    )
    facts = facts.merge(cells, on="gene_id", how="left", validate="one_to_one")
    ihc = read("normal-ihc.csv.gz")
    ihc = ihc[~ihc.Category.isin(["reproductive", "disposable"])].copy()
    # HPA contains placeholder rows for absent assays; these are not coverage.
    ihc = ihc.dropna(subset=["Tissue", "Cell type", "Level"])
    ihc["detected"] = ihc.Level.isin(["Low", "Medium", "High", "Ascending", "Descending"])
    ihc["reliable"] = ihc.Reliability.isin(["Approved", "Supported", "Enhanced"])
    facts["ihc_detected"] = facts.gene_id.isin(ihc.loc[ihc.detected, "Gene"])
    facts["ihc_reliable_detected"] = facts.gene_id.isin(
        ihc.loc[ihc.detected & ihc.reliable, "Gene"]
    )
    facts["ihc_uncertain_detected"] = facts.gene_id.isin(
        ihc.loc[ihc.detected & ~ihc.reliable, "Gene"]
    )
    facts["ihc_records"] = facts.gene_id.map(ihc.groupby("Gene").size())
    facts["ihc_reliable_records"] = facts.gene_id.map(ihc[ihc.reliable].groupby("Gene").size())
    facts["ihc_missing"] = facts.ihc_records.isna()
    facts = facts.sort_values("gene_id").reset_index(drop=True)
    write(facts, out / "quantitative-facts.csv.gz")
    result = policy.classify(facts)
    annotations = [
        "Expression category",
        "Expression detail",
        "Expression forbidden",
        "Safety veto",
        "Accepted normal protein",
        "Antibody concern",
        "Donor uncertainty",
        "Normal expression flags",
        "Risk assessment",
    ]
    annotations += [c for c in old if c.startswith("Normal breadth ")]
    result = result.merge(
        old[["Gene ID", *annotations]].rename(columns={"Gene ID": "gene_id"}),
        on="gene_id",
        how="left",
        validate="one_to_one",
    )
    write(result, out / "gene-results.csv.gz")
    write(result[result.selected], out / "selected-genes.csv")
    named = [
        "PRAME",
        "MAGEA4",
        "CTAG1A",
        "CTAG1B",
        "CTAG2",
        "XAGE1A",
        "XAGE1B",
        "SSX2",
        "CT83",
        "TBXT",
        "GPC3",
        "AFP",
        "CLDN6",
        "SALL4",
        "MAGEA3",
        "MAGEA12",
        "SSX1",
        "CTCFL",
        "PLAC1",
        "WT1",
        "BIRC5",
        "EGFR",
        "ERBB2",
        "GAGE1",
        "SSX4",
        "CT47A1",
        "TERT",
        "VGLL1",
        "CT45A1",
    ]
    write(result[result.symbol.isin(named)], out / "named-gene-results.csv")
    counts = (
        result.groupby(["selected", "Expression category"], dropna=False)
        .size()
        .rename("genes")
        .reset_index()
    )
    write(counts, out / "expression-status-counts.csv")
    failures = result.groupby("first_failure").size().rename("genes").reset_index()
    write(failures, out / "selection-counts.csv")
    rows = []
    benchmarks = {
        "PRAME",
        "MAGEA4",
        "CTAG1A",
        "CTAG1B",
        "CTAG2",
        "XAGE1A",
        "XAGE1B",
        "SSX2",
        "CT83",
        "TBXT",
    }
    for bulk, cell, dev_floor in itertools.product(
        [3, 5, 6, 10], [5, 10, 20, 25, 30, 50], [5, 10, 20]
    ):
        t = policy.Thresholds(
            adult_bulk_ntpm=bulk, adult_cell_ntpm=cell, developmental_tpm=dev_floor
        )
        r = policy.classify(facts, thresholds=t)
        selected = r[r.selected]
        rows.append(
            {
                **asdict(t),
                "selected": len(selected),
                "benchmarks_retained": len(benchmarks & set(selected.symbol)),
                "benchmarks_lost": "|".join(sorted(benchmarks - set(selected.symbol))),
            }
        )
    write(pd.DataFrame(rows), out / "threshold-sensitivity.csv")
    s = result[result.selected]
    summary = {
        "policy": policy.POLICY_ID,
        "thresholds": asdict(policy.Thresholds()),
        "genes": len(result),
        "selected": len(s),
        "reproductive_only": int((s.reproductive_route & ~s.developmental_route).sum()),
        "developmental_only": int((s.developmental_route & ~s.reproductive_route).sum()),
        "both": int((s.reproductive_route & s.developmental_route).sum()),
        "selected_missing_ihc": int(s.ihc_missing.sum()),
        "selected_expression_categories": s["Expression category"].value_counts().to_dict(),
        "developmental_genes_measured": int(facts.developmental_gene_measured.sum()),
        "developmental_embryos": 1,
        "reproductive_tissues": sorted(cta_tissues.CORE_REPRODUCTIVE_TISSUES),
        "adult_bulk_exempt_reproductive_tissues": sorted(cta_tissues.EXTENDED_REPRODUCTIVE_TISSUES),
        "first_failure_counts": dict(zip(failures.first_failure, failures.genes, strict=True)),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    sources = json.loads((prepared / "cell-source-receipts.json").read_text())
    for source, version in [
        ("hpa_rna_consensus", "v23"),
        ("hpa_normal_tissue", "v23"),
        ("hpa_cancer_rna", "v25.1-1"),
    ]:
        r = reference_data.provenance(source, version, verify_content=True)
        sources.append({k: r[k] for k in ["name", "version", "url", "sha256"]})
    manifest = {
        "policy": policy.POLICY_ID,
        "prepared_inputs": inputs,
        "sources": sources,
        "developmental_source": json.loads(
            (ROOT / "docs/audits/cta-quantitative-20261026/developmental-source.json").read_text()
        ),
        "reference": facts.reference.iloc[0],
        "annotation_source_sha256": digest(ROOT / "oncoref/data/cta-membership-gene-sets.csv.gz"),
        "builder_sha256": digest(__file__),
        "preparer_sha256": digest(ROOT / "scripts/prepare_cta_quantitative.py"),
        "classifier_sha256": digest(ROOT / "oncoref/cta_quantitative.py"),
        "outputs": {p.name: digest(p) for p in sorted(out.iterdir())},
    }
    (out / "run-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if publish:
        shutil.copyfile(
            out / "gene-results.csv.gz", ROOT / "oncoref/data/cta-quantitative-gene-sets.csv.gz"
        )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--publish", action="store_true")
    args = parser.parse_args()
    build(args.prepared, args.out, publish=args.publish)

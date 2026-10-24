#!/usr/bin/env python3
"""Join the frozen CTA decisions to cancer occurrence and cited biology.

No classifications are changed. Every nomination gets a row; missing source
identities remain missing. RNA/IHC cohorts, units and denominators stay separate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from oncoref import hpa_cancer, protein_annotations, reference_data

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/audits/cta-neither-review-20261024"
OUT = ROOT / "docs/audits/cta-atlas-20261024"


def joined(values):
    return " | ".join(dict.fromkeys(str(v) for v in values if pd.notna(v) and str(v)))


def normal_summary(contexts, unit):
    rows = contexts.loc[contexts.unit.eq(unit) & ~contexts.permitted].copy()
    rows = rows.sort_values(
        ["maximum", "source", "tissue", "cell_type"], ascending=[False, True, True, True]
    )
    rows = rows.drop_duplicates("gene_id")
    return {
        r.gene_id: f"{r.maximum:.4g} {unit}; {r.tissue or 'unsited'}/{r.cell_type or 'bulk'}; {r.source}; {r.category or 'unresolved'}"
        for r in rows.itertuples()
    }


def top_rna(rows):
    measured = rows.loc[rows.samples.ge(20)]
    if measured.empty:
        return "No matched RNA cohort with >=20 measured samples"
    best = measured.sort_values(
        ["prevalence_ptpm_ge_5", "expressed_samples_ptpm_ge_5", "cancer"],
        ascending=[False, False, True],
    ).head(3)
    return "; ".join(
        f"{r.cancer}: {int(r.expressed_samples_ptpm_ge_5)}/{int(r.samples)} ({r.prevalence_ptpm_ge_5:.1%})"
        for r in best.itertuples()
    )


def top_ihc(rows):
    measured = rows.loc[rows.total.ge(5)]
    if measured.empty:
        return "No matched IHC group with >=5 scored patients"
    best = measured.sort_values(
        ["prevalence_medium_high", "medium_high", "cancer"], ascending=[False, False, True]
    ).head(2)
    return "; ".join(
        f"{r.cancer}: {int(r.medium_high)}/{int(r.total)} ({r.prevalence_medium_high:.1%})"
        for r in best.itertuples()
    )


def boundary_flags(genes, tissue):
    """Review flags only; never imply that changing one cutoff admits a gene."""
    flags = {}
    for row in genes.to_dict("records"):
        gid = row["Gene ID"]
        notes = []
        maxima = [row["Normal RNA maximum (TPM)"], row["Normal RNA maximum (nTPM)"]]
        known = [float(v) for v in maxima if pd.notna(v) and v != ""]
        if row["Category"] == "CTA" and any(24 <= v <= 30 for v in known):
            notes.append("Within 20% below the 30 ceiling in at least one native RNA unit")
        if row["Detail"] == "High RNA" and known and 30 < max(known) <= 33:
            notes.append(
                "All audited maxima are within 10% above 30; other gates still require review"
            )
        if row["Detail"] in {"RNA", "Selectivity"}:
            t = tissue.loc[tissue.Ensembl_Gene_ID.eq(gid)]
            ratios = []
            for numerator, denominator in [
                ("rna_ntpm", "rna_ceiling_ntpm"),
                ("rna_ntpm", "continuous_selectivity_ceiling_ntpm"),
                ("rna_ntpm", "reproductive_rna_max_ntpm"),
            ]:
                z = t[numerator] / t[denominator].replace(0, np.nan)
                ratios.extend(z.loc[z.gt(1)].tolist())
            if ratios and max(ratios) <= 1.2:
                notes.append(
                    "Failed bulk numeric comparisons are within 20% of their cutoffs; not an automatic rescue"
                )
        flags[gid] = "; ".join(notes)
    return flags


def build(out=OUT):
    out.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(AUDIT / "gene-results.csv.gz", low_memory=False)
    genes = genes.loc[genes["Previous CTA screen"].ne("Not nominated")].copy()
    assert len(genes) == 2540 and not genes["Gene ID"].duplicated().any()
    genes = genes.sort_values(["Gene", "Gene ID"], na_position="last").reset_index(drop=True)
    ids = set(genes["Gene ID"])
    contexts = pd.read_csv(AUDIT / "rna-context-maxima.csv.gz", keep_default_na=False)
    contexts = contexts.loc[contexts.gene_id.isin(ids)]
    rna = hpa_cancer.hpa_cancer_rna_prevalence(ids)
    ihc = hpa_cancer.hpa_cancer_ihc_prevalence(ids)
    annotations = protein_annotations.gene_annotations(ids)
    annotations.to_csv(
        out / "biology.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )
    rna.to_csv(out / "cancer-rna.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    ihc.to_csv(out / "cancer-ihc.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    contexts.to_csv(
        out / "normal-contexts.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )
    genes["Highest normal TPM context"] = genes["Gene ID"].map(normal_summary(contexts, "TPM"))
    genes["Highest normal nTPM context"] = genes["Gene ID"].map(normal_summary(contexts, "nTPM"))
    genes["Cancer RNA at >=5 pTPM"] = (
        genes["Gene ID"]
        .map({gid: top_rna(g) for gid, g in rna.groupby("gene_id")})
        .fillna("No exact-ID RNA coverage in the HPA cancer snapshot")
    )
    genes["Cancer IHC Medium/High"] = (
        genes["Gene ID"]
        .map({gid: top_ihc(g) for gid, g in ihc.groupby("gene_id")})
        .fillna("No exact-ID IHC coverage in the HPA cancer snapshot")
    )
    genes["Cancer RNA cohorts"] = (
        genes["Gene ID"].map(rna.groupby("gene_id").size()).fillna(0).astype(int)
    )
    genes["Cancer IHC groups"] = (
        genes["Gene ID"].map(ihc.groupby("gene_id").size()).fillna(0).astype(int)
    )
    for source, dest in [
        ("accession", "UniProt entries"),
        ("url", "Biology sources"),
        ("function", "Biological function"),
        ("tissue_specificity", "Curated tissue biology"),
        ("disease", "Disease biology"),
        ("caution", "Biology caveats"),
        ("location", "Protein location"),
        ("function_pubmed_ids", "Function PubMed IDs"),
        ("tissue_specificity_pubmed_ids", "Tissue PubMed IDs"),
    ]:
        genes[dest] = (
            genes["Gene ID"].map(annotations.groupby("gene_id")[source].agg(joined)).fillna("")
        )
    genes["Shared UniProt entry"] = (
        genes["Gene ID"].map(annotations.groupby("gene_id").shared_entry.any()).eq(True)
    )
    genes["Biology coverage"] = np.where(
        genes["UniProt entries"].ne(""),
        "Reviewed UniProt annotation; not a new primary-paper review",
        "No exact Ensembl mapping to a reviewed UniProt entry",
    )
    trials = pd.read_csv(ROOT / "oncoref/data/cta-trial-evidence.csv", keep_default_na=False)
    names = {
        "registered_intervention": "Trial registration only",
        "approved_therapy": "Approved product",
        "shared_epitope_clinical": "Clinical shared-epitope target",
        "family_label_clinical": "Clinical family label; locus unresolved",
        "primary_clinical_report": "Clinical publication",
        "institutional_clinical_report": "Institutional clinical report",
        "clinical_safety_context": "Clinical safety concern",
    }
    trials["Evidence"] = trials.evidence_kind.map(names)
    assert trials.Evidence.notna().all()
    genes["Trial evidence"] = (
        genes["Gene ID"]
        .map(trials.groupby("Ensembl_Gene_ID").Evidence.agg(joined))
        .fillna("No trial link in the curated intake; not an exhaustive negative search")
    )
    genes["Trial sources"] = (
        genes["Gene ID"].map(trials.groupby("Ensembl_Gene_ID").source_url.agg(joined)).fillna("")
    )
    reviews = yaml.safe_load((OUT / "literature.yml").read_text())["records"]
    review_rows = []
    for review in reviews:
        for gene in review["genes"]:
            review_rows.append({"Gene": gene, **{k: v for k, v in review.items() if k != "genes"}})
    review_frame = pd.DataFrame(review_rows)
    review_frame.to_csv(out / "literature-reviews.csv", index=False)
    genes["Reviewed literature"] = genes.Gene.map(
        review_frame.groupby("Gene").finding.agg(joined)
    ).fillna(
        "No individual clinical-outcome review in this report; no inference of safety or ineffectiveness"
    )
    genes["Reviewed literature sources"] = genes.Gene.map(
        review_frame.groupby("Gene").source.agg(joined)
    ).fillna("")
    genes["Reviewed evidence types"] = genes.Gene.map(
        review_frame.groupby("Gene").evidence.agg(joined)
    ).fillna("")
    membership = pd.read_csv(
        ROOT / "oncoref/data/cta-publication-membership.csv", keep_default_na=False
    )
    genes["Nomination sources"] = (
        genes["Gene ID"]
        .map(membership.groupby("Ensembl_Gene_ID").source_tag.agg(joined))
        .fillna("")
    )
    genes["Nomination scope"] = (
        genes["Gene ID"]
        .map(membership.groupby("Ensembl_Gene_ID").validation_scope.agg(joined))
        .fillna("")
    )
    followup = pd.read_csv(AUDIT / "follow-up.csv", keep_default_na=False)
    for col in ["Review question", "Reconsideration evidence", "Limitations"]:
        genes[col] = genes.Gene.map(followup.set_index("Gene")[col]).fillna("")
    tissue = pd.read_csv(
        ROOT / "docs/audits/cta-screen-20261023/tissue-evidence.csv.gz", low_memory=False
    )
    genes["Boundary review"] = genes["Gene ID"].map(boundary_flags(genes, tissue))
    conflicts = yaml.safe_load((OUT / "retained-literature-conflicts.yml").read_text())["records"]
    conflict_rows = [
        {"Gene": symbol, **{k: v for k, v in row.items() if k != "genes"}}
        for row in conflicts
        for symbol in row["genes"]
    ]
    conflict_frame = pd.DataFrame(conflict_rows)
    assert not conflict_frame.Gene.duplicated().any()
    assert set(conflict_frame.Gene) <= set(genes.loc[genes.Category.eq("CTA"), "Gene"])
    for source, column in [
        ("priority", "Literature conflict priority"),
        ("finding", "Literature conflict"),
        ("action", "Literature follow-up"),
        ("source", "Literature conflict source"),
    ]:
        genes[column] = genes.Gene.map(conflict_frame.set_index("Gene")[source]).fillna("")
    genes.loc[genes["Literature conflict"].ne("")].to_csv(
        out / "retained-literature-conflicts.csv", index=False
    )
    # Preserve source-specific margins and all other vetoes for adjudication.
    tissue = tissue.loc[tissue.Ensembl_Gene_ID.isin(ids)]
    tissue.to_csv(
        out / "bulk-tissue-rules.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )
    genes.to_csv(out / "all-genes.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    genes[
        [
            "Gene",
            "Gene ID",
            "Category",
            "Detail",
            "Reproductive RNA (nTPM)",
            "Highest normal nTPM context",
            "Highest normal TPM context",
            "Cancer RNA at >=5 pTPM",
            "Cancer IHC Medium/High",
            "Boundary review",
            "Literature conflict",
            "Literature follow-up",
            "Reviewed evidence types",
            "Reviewed literature sources",
            "Review question",
            "Reconsideration evidence",
        ]
    ].to_csv(out / "quick-table.csv", index=False)
    for category, filename in [("CTA", "retained"), ("Unsafe", "excluded"), ("Neither", "neither")]:
        genes.loc[genes.Category.eq(category)].to_csv(out / (filename + ".csv"), index=False)
    genes.loc[genes["Boundary review"].ne("")].to_csv(out / "boundaries.csv", index=False)
    genes.loc[genes["Trial sources"].ne("")].to_csv(out / "clinical-targets.csv", index=False)
    genes.loc[genes["Review question"].ne("")].to_csv(out / "reconsideration.csv", index=False)
    unsafe = genes.Category.eq("Unsafe")
    assert not genes.loc[~unsafe, "High normal RNA"].any()
    assert (
        not genes.loc[~unsafe, ["Normal RNA maximum (TPM)", "Normal RNA maximum (nTPM)"]]
        .gt(30)
        .any()
        .any()
    )
    summary = {
        "policy": "oncoref.target-sets.2026-10-24.1",
        "date": "2026-10-24",
        "counts": genes.Category.value_counts().to_dict(),
        "details": genes.groupby(["Category", "Detail"])
        .size()
        .rename("Genes")
        .reset_index()
        .to_dict("records"),
        "coverage": {
            "UniProt exact-ID genes": annotations.gene_id.nunique(),
            "UniProt function genes": annotations.loc[
                annotations.function.ne("")
            ].gene_id.nunique(),
            "Cancer RNA genes": rna.gene_id.nunique(),
            "Cancer IHC genes": ihc.gene_id.nunique(),
            "Targeted literature review genes": int(genes["Reviewed evidence types"].ne("").sum()),
            "Trial-linked genes": int(genes["Trial sources"].ne("").sum()),
            "Boundary-flagged genes": int(genes["Boundary review"].ne("").sum()),
            "Retained genes with literature conflicts": len(conflict_frame),
        },
        "cancer_rna_rows": len(rna),
        "cancer_ihc_rows": len(ihc),
        "scope": "All 2540 mapped nominations. BAGE/GAGE2B unmapped intake shown separately; not assigned false-negative measurements.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    inputs = [
        AUDIT / "gene-results.csv.gz",
        AUDIT / "rna-context-maxima.csv.gz",
        AUDIT / "follow-up.csv",
        OUT / "literature.yml",
        OUT / "retained-literature-conflicts.yml",
        ROOT / "oncoref/data/cta-trial-evidence.csv",
        ROOT / "oncoref/data/cta-publication-membership.csv",
        ROOT / "docs/audits/cta-screen-20261023/tissue-evidence.csv.gz",
        ROOT / "oncoref/protein_annotations.py",
        ROOT / "oncoref/hpa_cancer.py",
        ROOT / "oncoref/data/protein-annotation-sources.json",
        ROOT / "scripts/build_cta_atlas.py",
        ROOT / "scripts/render_cta_atlas.py",
        ROOT / "scripts/publish_cta_atlas_files.py",
        ROOT / "scripts/templates/cta-atlas-explorer.html",
    ]
    provenance = {
        "policy": summary["policy"],
        "created": summary["date"],
        "inputs": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs
        },
        "references": [
            reference_data.provenance(n, verify_content=True)
            for n in ["hpa_cancer_rna", "hpa_cancer_ihc", "uniprot_reviewed_human"]
        ],
        "limitations": [
            "Not a systematic literature review of every gene; reviewed clinical outcomes and UniProt annotations are labelled separately.",
            "HPA cancer pTPM is not normal-tissue nTPM; cohorts and modalities are not paired.",
            "Gene-level cancer occurrence does not establish isoform expression, peptide presentation, malignant-cell localization or therapeutic eligibility.",
            "Native Ensembl IDs are joined exactly across releases; absent IDs and unreviewed genes remain coverage gaps.",
            "Earlier Unsafe assignments were not all re-audited against every new normal source; missing maxima are not zero.",
            "Gene counts are not unique protein counts; shared exact sequences use Unsafe > Neither > CTA.",
        ],
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    build(args.output)

#!/usr/bin/env python3
"""Package current CTA decisions and publish a species-scoped source inventory."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_report_common import sha256, write_csv  # noqa: E402

from oncoref import canine, human  # noqa: E402
from oncoref.cta_selection import RESULT_COLUMNS  # noqa: E402
from oncoref.expression_registry import expression_sources  # noqa: E402

AUDIT = ROOT / "docs/audits/cta-release-1.21.0"
DATA = ROOT / "oncoref/data"


def source_inventory() -> pd.DataFrame:
    """List biological roles separately from native availability and selection use."""
    rows = []

    def add(
        taxon: int,
        source: str,
        name: str,
        categories: list[str],
        status: str,
        unit: str,
        stage: str,
        url: str,
        limitation: str,
    ) -> None:
        for category in categories:
            rows.append(
                {
                    "taxon_id": taxon,
                    "source_id": source,
                    "name": name,
                    "evidence_category": category,
                    "status": status,
                    "native_unit": unit,
                    "stage": stage,
                    "source_url": url,
                    "limitations": limitation,
                }
            )

    add(
        9606,
        "hpa_rna_consensus",
        "HPA v23 / GTEx consensus",
        ["normal_bulk_rna", "reproductive_expression"],
        "used_in_selection",
        "nTPM",
        "adult; placenta handled separately",
        "https://v23.proteinatlas.org/about/download",
        "Bulk source summaries; not exhaustive pediatric, pregnancy or rare-cell coverage",
    )
    add(
        9606,
        "tyser_E-MTAB-9388",
        "Tyser human gastrulation atlas",
        ["developmental_expression"],
        "used_in_selection",
        "TPM",
        "Carnegie stage 7",
        "https://doi.org/10.1038/s41586-021-04158-y",
        "One embryo; cell support does not establish independent-donor replication",
    )
    add(
        9606,
        "oncoref_tumor_union",
        "OncoRef source-mapped tumor union",
        ["tumor_expression"],
        "used_in_selection",
        "TPM",
        "cohort-specific",
        "https://github.com/pirl-unc/oncoref/tree/v1.21.0/docs/audits/cta-expanded-funnel-20261027/tumor-union",
        "132 compatible cohorts; source namespaces deduplicate views, not patients across studies",
    )
    add(
        9606,
        "hpa_cancer",
        "HPA cancer sample measurements",
        ["tumor_expression"],
        "used_in_selection",
        "pTPM",
        "cohort-specific",
        "https://v23.proteinatlas.org/about/download",
        "Protein-coding RNA denominator; no numerical protein or peptide abundance",
    )
    add(
        9606,
        "hpa_single_cell_tissue",
        "HPA tissue/cell clusters",
        ["normal_cell_rna"],
        "used_in_selection",
        "nTPM",
        "source-labelled adult; reproductive locations explicit",
        "https://v23.proteinatlas.org/about/download",
        "At least 20 reads; missing cell types and disease states remain uncovered",
    )
    for source, name, unit in [
        ("hpa_immune_cell", "HPA immune summaries", "nTPM"),
        ("hpa_immune_monaco", "Monaco immune summaries", "TPM"),
    ]:
        add(
            9606,
            source,
            name,
            ["immune_rna"],
            "used_in_selection",
            unit,
            "adult",
            "https://v23.proteinatlas.org/about/download",
            "Sorted populations retain native units; no matched negative-IHC relaxation",
        )
    add(
        9606,
        "dice",
        "DICE immune distributions",
        ["immune_rna"],
        "used_in_selection",
        "TPM",
        "adult",
        "https://dice-database.org/",
        "15 distributions; third-highest value; unpaired locus distributions retain bounds; optional coverage",
    )
    add(
        9606,
        "GSE155683_adult",
        "Collin adult corneal atlas",
        ["normal_cell_rna"],
        "used_in_selection",
        "CP10K",
        "adult",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155683",
        "Two supported sample groups; optional coverage; CP10K is not length-normalized TPM",
    )
    add(
        9606,
        "hpa_normal_tissue",
        "HPA normal-tissue IHC",
        ["normal_protein"],
        "used_in_selection",
        "HPA categorical level and reliability",
        "adult",
        "https://v23.proteinatlas.org/about/assays%2Bannotation",
        "Reliability describes evidence support, not a calibrated false-negative probability",
    )
    add(
        9606,
        "external_ihc",
        "Curated external staining observations",
        ["normal_protein"],
        "used_in_selection",
        "native categorical findings",
        "source-specific",
        "https://github.com/pirl-unc/oncoref/tree/v1.21.0/docs/audits/external-ihc-20261027",
        "Reagent, tissue, cell type and validation scope remain explicit; no inherited HPA grade",
    )
    add(
        9606,
        "cta_policy_findings",
        "Reviewed normal protein and clinical exclusions",
        ["normal_protein", "clinical_safety"],
        "used_in_selection",
        "categorical reviewed findings",
        "source-specific",
        "https://github.com/pirl-unc/oncoref/blob/v1.21.0/oncoref/data/cta-policy-findings.csv",
        "Qualitative reported Low is distinct from HPA Low and numerical intensity; clinical evidence is product-specific",
    )
    papers = pd.read_csv(DATA / "cta-publication-sources.csv").drop_duplicates("source_tag")
    for row in papers.itertuples():
        add(
            9606,
            row.source_tag,
            row.citation,
            ["nomination"],
            "annotation",
            "gene-list membership",
            "publication-specific",
            f"https://doi.org/{row.doi}" if pd.notna(row.doi) else row.source_url,
            "Comparison annotation only; no membership bypass of the quantitative funnel",
        )
    dog_roles = {
        "bgee-srp009687": ["normal_bulk_rna", "reproductive_expression"],
        "bgee-srp181218": ["normal_bulk_rna", "reproductive_expression"],
        "barkbase": ["normal_bulk_rna", "reproductive_expression", "developmental_expression"],
        "pfizer-dog-protein": ["normal_protein"],
        "doga": ["normal_bulk_rna", "reproductive_expression", "developmental_expression"],
        "pfizer-dog-rna": ["normal_bulk_rna"],
        "gse304066-osteosarcoma": ["tumor_expression"],
    }
    urls = {
        "bgee-srp009687": "https://www.ncbi.nlm.nih.gov/sra/?term=SRP009687",
        "bgee-srp181218": "https://www.ncbi.nlm.nih.gov/sra/?term=SRP181218",
        "barkbase": "https://pubmed.ncbi.nlm.nih.gov/31181663/",
        "pfizer-dog-protein": "https://doi.org/10.1038/s41597-026-07512-9",
        "pfizer-dog-rna": "https://doi.org/10.1038/s41597-026-07512-9",
        "doga": "https://doi.org/10.1038/s41467-024-52798-1",
        "gse304066-osteosarcoma": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304066",
    }
    for source in expression_sources(species="dog"):
        status = {"metadata_only": "metadata_only", "blocked_provenance": "blocked_provenance"}.get(
            source.acquisition_status, "native_matrix"
        )
        add(
            9615,
            source.id,
            source.id,
            dog_roles[source.id],
            status,
            source.unit or "unresolved",
            "sample metadata; adult/embryonic must remain separate",
            urls[source.id],
            source.expression_contract.provenance_notes,
        )
    native = canine.cta_snapshot().as_dict()
    parents = {record["id"]: record for record in native["evidence"]}
    for source in native["sources"]:
        facts = [
            fact
            for fact in native["facts"]
            if parents[fact["evidence_id"]]["source_id"] == source["id"]
        ]
        categories = {"nomination"}
        for fact in facts:
            if fact["observation"]["status"] == "not_reported":
                continue
            role = fact["tissue_role"]
            protein = fact["measurement"]["feature_level"] in {"protein", "protein_group"}
            if role == "reproductive":
                categories.add("reproductive_expression")
            elif role == "tumor":
                categories.add("tumor_expression")
            elif role == "somatic":
                categories.add("normal_protein" if protein else "normal_bulk_rna")
        add(
            9615,
            source["id"],
            source["id"],
            sorted(categories),
            "native_observations",
            "native PCR, categorical protein and reported study summaries",
            "reported study scope",
            source["url"],
            "14 nominees across two studies; family assays, ambiguous identities and conflicting counterevidence remain explicit; no quantitative CTA policy",
        )
    for source, name, categories, unit, stage, url, limitation in [
        (
            "tabula_muris",
            "Tabula Muris / Senis",
            ["normal_cell_rna", "immune_rna"],
            "platform-specific counts and processed expression",
            "adult and ageing",
            "https://doi.org/10.1038/s41586-020-2496-1",
            "Literature identified; matrices, strain/sex, stage and reference reconciliation are not integrated",
        ),
        (
            "moca",
            "Mouse Organogenesis Cell Atlas",
            ["developmental_expression"],
            "UMI counts",
            "E9.5-E13.5",
            "https://doi.org/10.1038/s41586-019-0969-x",
            "Literature identified; UMI counts cannot use human TPM cutoffs without a declared measurement policy",
        ),
        (
            "mouse_gastrulation",
            "Mouse gastrulation and early organogenesis",
            ["developmental_expression"],
            "native single-cell RNA counts",
            "embryonic",
            "https://doi.org/10.1038/s41586-019-0933-9",
            "Literature identified; embryo, stage and strain provenance need ingestion",
        ),
        (
            "gxd",
            "MGI Gene Expression Database",
            ["developmental_expression", "normal_protein", "normal_bulk_rna"],
            "assay-specific RNA and protein observations",
            "stage-specific",
            "https://www.informatics.jax.org/mgihome/GXD/aboutGXD.shtml",
            "Literature identified; distinguish measured not-detected from not-analyzed search results; no complete normal panel inferred",
        ),
    ]:
        add(10090, source, name, categories, "not_integrated", unit, stage, url, limitation)
    return pd.DataFrame(rows)


def build() -> None:
    """Export a small wheel-backed table, preserving every group and member identity."""
    source = AUDIT / "protein-group-classifications.csv.gz"
    decisions = pd.read_csv(source, low_memory=False)
    source_identity = "sha256:" + sha256(source)
    reference = human.cta_snapshot().as_dict()["references"][0]["identity"]
    rows = []
    flag_columns = [
        "Antibody concern",
        "external_review_flag",
        "low_normal_protein_flag",
        "limited_protein_review_flag",
        "stem_progenitor_expression_flag",
    ]
    for row in decisions.to_dict("records"):
        identities = [
            {
                "gene_id": member["canonical_gene_id"],
                "symbol": member["canonical_symbol"],
                "namespace": "Ensembl gene",
                "historical_gene_ids": json.loads(member["historical_ensembl_gene_ids"])
                if member["historical_ensembl_gene_ids"] is not None
                else None,
                "alternate_gene_ids": json.loads(member["alternate_ensembl_gene_ids"])
                if member["alternate_ensembl_gene_ids"] is not None
                else None,
                "synonyms": json.loads(member["gene_name_synonyms"])
                if member["gene_name_synonyms"] is not None
                else None,
                "identity_status": member["identity_status"],
                "history_status": "curated_snapshot"
                if member["identity_status"] == "matched"
                else "outside_identity_catalog",
            }
            for member in json.loads(row["member_identities"])
        ]
        rows.append(
            {
                "taxon_id": 9606,
                "species": "human",
                "entity_id": row["protein_equivalence_group_id"],
                "entity_kind": "protein_equivalence_group",
                "label": row["protein_group_label"],
                "gene_ids": row["member_gene_ids"],
                "gene_identities": json.dumps(identities, separators=(",", ":")),
                "identity_status": "resolved"
                if all(member["identity_status"] == "matched" for member in identities)
                else "partial",
                "reference_id": reference["reference_id"],
                "category": row["category"],
                "assessment_status": "screened",
                "policy_id": row["policy"],
                "reasons": json.dumps(
                    row["exclusion_reasons"].split(" | ")
                    if pd.notna(row["exclusion_reasons"])
                    else []
                ),
                "evidence_flags": json.dumps([column for column in flag_columns if row[column]]),
                "source_snapshot_id": source_identity,
            }
        )
    write_csv(pd.DataFrame(rows, columns=RESULT_COLUMNS), DATA / "cta-selection-human.csv.gz")
    inventory = source_inventory()
    write_csv(inventory, DATA / "cta-selection-sources.csv")
    write_csv(inventory, AUDIT / "species-source-inventory.csv")
    metadata = {
        "schema": "oncoref.cta-selection",
        "schema_version": 1,
        "release": "1.21.0",
        "data_sha256": {
            name: sha256(DATA / name)
            for name in ["cta-selection-human.csv.gz", "cta-selection-sources.csv"]
        },
        "species": {
            "human": {
                "selection_status": "available",
                "evidence_status": "partial",
                "policy_id": decisions.policy.iloc[0],
                "universe": "19863 protein-equivalence groups covering 20097 gene loci",
                "reference": reference,
                "limitations": [
                    "Research prioritization calibrated using PRAME; no clinical safety or pediatric validation",
                    "Only one embryonic atlas embryo; later fetal atlases are not integrated",
                    "Native RNA units remain separate and nominal thresholds lack assay-specific toxicity calibration",
                    "Low normal protein is a tolerated flag; author descriptions and HPA categories retain different scales",
                    "Representative-protein equivalence is not exhaustive isoform identity",
                ],
            },
            "dog": {
                "selection_status": "unavailable",
                "evidence_status": "partial",
                "policy_id": None,
                "universe": "14 curated native nominations; not a genome-wide classification",
                "reference": canine.cta_candidates()[0]["reference"],
                "limitations": [
                    "No validated species-specific quantitative selection policy",
                    "Native normal matrices need reference reconciliation and source-specific scale handling",
                    "Tumor cohort integration is blocked by unit/reference provenance",
                    "Unassessed nominations are Neither; human expression and safety evidence never transfer by orthology",
                ],
            },
            "mouse": {
                "selection_status": "unavailable",
                "evidence_status": "unavailable",
                "policy_id": None,
                "universe": "No curated mouse nomination or selection snapshot",
                "reference": None,
                "limitations": [
                    "Source inventory identifies relevant resources but their data are not integrated",
                    "Zero returned rows is unavailable curation, not zero biological CTA genes",
                    "A native tumor cohort, reference-bound gene mapping and normal-panel policy are required",
                ],
            },
        },
        "inputs": {
            str(source.relative_to(ROOT)): sha256(source),
            "report_manifest": sha256(AUDIT / "manifest.json"),
        },
    }
    (DATA / "cta-selection.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(
        json.dumps(
            {
                "rows": len(rows),
                "categories": decisions.category.value_counts().to_dict(),
                "source_category_rows": len(inventory),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    build()

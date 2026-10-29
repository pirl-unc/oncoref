#!/usr/bin/env python3
"""Import the declared Sasaroli vascular and Hirbo placental nomination sets.

The source lists annotate the quantitative selection; membership never supplies
an expression value or changes an admission threshold. Original names, exact
source locations, noncoding rows and ambiguous identities remain inspectable.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref.gene_matching import gene_identity_catalog, match_reference_genes  # noqa: E402

DATA = ROOT / "oncoref/data"
AUDIT = ROOT / "docs/audits/cta-source-walkthrough-20261028/publication-intake"
SASAROLI = "Sasaroli2011_vascular"
HIRBO = "Hirbo2015_placenta"
# The original paper explicitly expands each old label in Results, paragraph
# beginning 'More than 30 genes'. SLC11A1 is natural resistance-associated
# macrophage protein 1, not the distinct SLC11A2 locus.
VASCULAR_NAMES = [
    ("ADAM12", "ADAM12"),
    ("BLAME", "SLAMF8"),
    ("CDCP1", "CDCP1"),
    ("CSPG2", "VCAN"),
    ("DR6", "TNFRSF21"),
    ("EGFL6", "EGFL6"),
    ("ESM1", "ESM1"),
    ("FLJ46072", "FAM83H"),
    ("LGALS3BP", "LGALS3BP"),
    ("NRAMP", "SLC11A1"),
    ("SDC1", "SDC1"),
    ("ST14", "ST14"),
    ("TSG6", "TNFAIP6"),
]


def build() -> None:
    """Read pinned published rows and update the shared publication registry."""
    workbook = AUDIT / "hirbo-s1.xlsx"
    table = pd.read_excel(workbook, header=None)
    if (
        table.iloc[2, 9] != "Placenta"
        or table.iloc[1, 28] != "Fast evolving placentally enriched genes"
    ):
        raise ValueError("Unexpected Hirbo supplementary-table layout")
    placenta = table.iloc[3:, 9].dropna()
    fast = set(table.iloc[3:, 28].dropna())
    if len(placenta) != 86 or len(fast) != 30 or not fast <= set(placenta):
        raise ValueError("Published Hirbo set sizes do not reconcile")
    original = []
    for old, name in VASCULAR_NAMES:
        original.append(
            {
                "source_tag": SASAROLI,
                "source_name": old,
                "match_name": name,
                "source_location": "Results: 13 promising genes; Figure 2B",
                "subset": "13 detailed tumor-vascular candidates",
            }
        )
    for row, name in placenta.items():
        classes = [
            str(table.iloc[2, c]) for c in [29, 30, 31] if name in set(table.iloc[3:, c].dropna())
        ]
        original.append(
            {
                "source_tag": HIRBO,
                "source_name": name,
                "match_name": name,
                "source_location": f"S1 Table, Sheet1!J{row + 1}",
                "subset": "placenta_enriched"
                + (";fast_evolving:" + ",".join(classes) if classes else ""),
            }
        )
    source = pd.DataFrame(original)
    identity = match_reference_genes(
        source,
        id_column=None,
        symbol_column="match_name",
        source="published_related_lists",
        source_version="2011/2015",
    )
    reviewed = pd.concat([source, identity], axis=1)
    reviewed.to_csv(AUDIT / "mapped-source-rows.csv", index=False)
    catalog = gene_identity_catalog().set_index("canonical_gene_id")
    intake_path = DATA / "cta-publication-membership.csv"
    intake = pd.read_csv(intake_path).fillna("")
    intake = intake[~intake.source_tag.isin([SASAROLI, HIRBO])]
    additions = []
    for r in reviewed.itertuples():
        matched = r.identity_status == "matched"
        gene = catalog.loc[r.canonical_gene_id] if matched else None
        protein = matched and gene.biotype == "protein_coding"
        additions.append(
            {
                "source_tag": r.source_tag,
                "source_gene_id": "",
                "source_symbol": r.source_name,
                "source_biotype": "",
                "evidence_type": "tumor_vascular_nomination"
                if r.source_tag == SASAROLI
                else "placental_RNA_enrichment",
                "validation_scope": "Published nomination; no quantitative admission bypass",
                "Ensembl_Gene_ID": r.canonical_gene_id,
                "Symbol": gene.canonical_symbol if matched else "",
                "biotype": gene.biotype if matched else "",
                "mapping_status": "mapped" if matched else r.identity_status,
                "protein_candidate_eligible": protein,
                "source_id_type": "symbol",
                "source_file": "sasaroli-transcribed-13.csv"
                if r.source_tag == SASAROLI
                else "hirbo-s1.xlsx",
                "source_sheet": "Results/Figure 2B" if r.source_tag == SASAROLI else "Sheet1",
                "source_rows": r.source_location,
                "source_subset": r.subset,
                "assay_resolution": "published_gene_nomination",
                "mapping_method": r.identity_method,
            }
        )
    pd.concat([intake, pd.DataFrame(additions)], ignore_index=True).to_csv(intake_path, index=False)
    transcribed = source[source.source_tag.eq(SASAROLI)]
    transcribed.to_csv(AUDIT / "sasaroli-transcribed-13.csv", index=False)
    registry_path = DATA / "cta-publication-sources.csv"
    registry = pd.read_csv(registry_path)
    template = registry.iloc[0]
    registry = registry[~registry.source_tag.isin([SASAROLI, HIRBO])]
    entries = []
    for tag, author, title, doi, url, filename, count, location, scope in [
        (
            SASAROLI,
            "Sasaroli et al. (2011)",
            "Novel surface targets and serum biomarkers from the ovarian cancer vasculature",
            "10.4161/cbt.12.3.16260",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC3230481/",
            "sasaroli-transcribed-13.csv",
            13,
            "Results: 13 promising genes; Figure 2B",
            "Tumor-vascular candidate list; includes secreted and transmembrane proteins; not a CTA-restricted list",
        ),
        (
            HIRBO,
            "Hirbo et al. (2015)",
            "Integrating Diverse Types of Genomic Data to Identify Genes that Underlie Adverse Pregnancy Phenotypes",
            "10.1371/journal.pone.0144155",
            "https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0144155.s002&type=supplementary",
            "hirbo-s1.xlsx",
            86,
            "S1 Table: J4:J89; AC4:AF33 tracks the 30 fast-evolving subset",
            "HPA v13 placenta-enriched nominations; 30/86 fast-evolving subset retained as metadata; not tumor validation",
        ),
    ]:
        entries.append(
            {
                "source_tag": tag,
                "citation": author,
                "title": title,
                "doi": doi,
                "source_url": url,
                "source_table": location,
                "published_rows": count,
                "input_sha256": hashlib.sha256((AUDIT / filename).read_bytes()).hexdigest(),
                "hash_scope": "Original supplementary XLSX bytes"
                if tag == HIRBO
                else "Transcribed published 13-member list with original names and explicit paper aliases, UTF-8 CSV bytes",
                "evidence_scope": scope,
                "canonical_reference_sha256": template.canonical_reference_sha256,
                "canonical_ensembl_release": int(template.canonical_ensembl_release),
                "retrieved_date": "2026-10-28",
                "paper_reported_genes": count,
            }
        )
    pd.concat([registry, pd.DataFrame(entries)], ignore_index=True).to_csv(
        registry_path, index=False
    )
    (AUDIT / "manifest.json").write_text(
        json.dumps(
            {
                "identity_catalog_sha256": catalog.identity_catalog_sha256.iloc[0],
                "sources": entries,
                "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "outputs": {
                    p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in AUDIT.iterdir()
                    if p.suffix == ".csv"
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(reviewed.groupby(["source_tag", "identity_status"]).size().to_string())


if __name__ == "__main__":
    build()

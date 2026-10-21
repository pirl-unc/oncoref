#!/usr/bin/env python3
"""Import the primary BAGE/SPAG9 nominations and reviewed NCBI intake identities.

Inputs are NCBI PubMed efetch XML (PMIDs 7895173 and 17332284) and Gene
esummary JSON (IDs 574 and 645037). No absent Ensembl locus is invented.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref.cta_landscape import resolve_landscape_membership  # noqa: E402
from oncoref.cta_regen import regenerate_cta_columns  # noqa: E402
from oncoref.cta_sources import add_publication_candidates  # noqa: E402
from oncoref.gene_ids import canonical_gene_space  # noqa: E402


def build(pubmed_xml, ncbi_json, retrieved_date, out):
    original = ROOT / "oncoref/data"
    membership = pd.read_csv(original / "cta-publication-membership.csv", keep_default_na=False)
    sources = pd.read_csv(original / "cta-publication-sources.csv", keep_default_na=False)
    articles = {
        article.findtext(".//PMID"): article
        for article in ET.fromstring(pubmed_xml.read_bytes()).findall(".//PubmedArticle")
    }
    genes = json.loads(ncbi_json.read_text())["result"]
    for gid, name in [("574", "BAGE"), ("645037", "GAGE2B")]:
        assert genes[gid]["nomenclaturesymbol"] == name
        assert genes[gid]["organism"]["taxid"] == 9606
    space = canonical_gene_space()
    canonical_hash = hashlib.sha256(space.to_csv(index=False).encode()).hexdigest()
    article_hash = hashlib.sha256(pubmed_xml.read_bytes()).hexdigest()
    gene_hash = hashlib.sha256(ncbi_json.read_bytes()).hexdigest()
    definitions = [
        (
            "Boel1995_BAGE",
            "7895173",
            "BAGE",
            "574",
            "Boël et al., Immunity (1995)",
            "10.1016/s1074-7613(95)80053-0",
            "Original BAGE tumor-antigen nomination and autologous CTL recognition; contemporary paralog resolution is not inferred",
            "Original BAGE target label; not every BAGE-family locus",
        ),
        (
            "Garg2007_SPAG9",
            "17332284",
            "SPAG9",
            "9043",
            "Garg et al., Clinical Cancer Research (2007)",
            "10.1158/1078-0432.CCR-06-2340",
            "Ovarian-cancer RNA/protein expression and patient antibody responses; no HLA presentation or T-cell validation claimed",
            "Gene-labelled RT-PCR, RNA in situ hybridization, IHC and patient serology; no isoform-specific inference",
        ),
    ]
    additions, metadata = [], []
    for tag, pmid, symbol, ncbi_id, citation, doi, scope, assay in definitions:
        article = articles[pmid]
        article_dois = {
            n.text.lower()
            for n in article.findall(".//ArticleId")
            if n.attrib.get("IdType") == "doi"
        }
        assert doi.lower() in article_dois
        title = "".join(article.find(".//ArticleTitle").itertext())
        additions.append(
            {
                "source_tag": tag,
                "source_gene_id": ncbi_id,
                "source_symbol": symbol,
                "source_biotype": "protein target in original paper",
                "evidence_type": "published_CTA_nomination",
                "validation_scope": scope,
                "source_id_type": "entrez",
                "source_file": pubmed_xml.name,
                "source_sheet": "PubMed abstract",
                "source_rows": f"PMID:{pmid}",
                "source_subset": "targeted_primary_nomination",
                "assay_resolution": assay,
            }
        )
        metadata.append(
            {
                "source_tag": tag,
                "citation": citation,
                "title": title,
                "doi": doi,
                "source_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "source_table": f"Targeted primary nomination: {symbol}; PMID {pmid}",
                "published_rows": 1,
                "input_sha256": article_hash,
                "hash_scope": "Original NCBI PubMed efetch XML bytes for PMIDs 17332284 and 7895173",
                "evidence_scope": scope,
                "canonical_reference_sha256": canonical_hash,
                "canonical_ensembl_release": ";".join(
                    sorted(space.ensembl_release.astype(str).unique())
                ),
                "retrieved_date": retrieved_date,
            }
        )
    tags = {r[0] for r in definitions}
    mapped = resolve_landscape_membership(pd.DataFrame(additions))
    membership = pd.concat(
        [membership[~membership.source_tag.isin(tags)], mapped], ignore_index=True
    )
    sources = pd.concat(
        [sources[~sources.source_tag.isin(tags)], pd.DataFrame(metadata)], ignore_index=True
    )
    reviews = []
    for name, gid in [("BAGE", "574"), ("GAGE2B", "645037")]:
        assert not space.symbol.fillna("").eq(name).any()
        reviews.append(
            {
                "source_symbol": name,
                "authority_id": f"NCBIGene:{gid}",
                "authority_symbol": name,
                "authority_url": f"https://www.ncbi.nlm.nih.gov/gene/{gid}",
                "review_status": "official_identity_without_canonical_ensembl_mapping",
                "reviewed_date": retrieved_date,
                "source_sha256": gene_hash,
                "canonical_reference_sha256": canonical_hash,
                "rationale": "Official human NCBI identity confirmed; no matching symbol in pinned Ensembl 115 canonical space. Retain in intake; do not substitute a paralog or impute normal expression.",
                "assay_caveat": (
                    "Original BAGE target and later panel label do not identify every BAGE-family locus"
                    if name == "BAGE"
                    else "Jamin microarray probes are shared across GAGE loci; nomination is not locus-specific expression validation"
                ),
            }
        )
    core = pd.read_csv(original / "cancer-testis-antigens.csv")
    extended = pd.read_csv(original / "cancer-testis-antigens-extended.csv")
    updated = add_publication_candidates(core, membership)
    new = updated.loc[~updated.Ensembl_Gene_ID.isin(core.Ensembl_Gene_ID)]
    if len(new):
        extended = pd.concat(
            [extended, regenerate_cta_columns(new, tissue_scope="extended")], ignore_index=True
        )
    out.mkdir(parents=True, exist_ok=True)
    for name, frame in {
        "cta-publication-membership.csv": membership,
        "cta-publication-sources.csv": sources,
        "cta-intake-identity-reviews.csv": pd.DataFrame(reviews),
        "cancer-testis-antigens.csv": updated,
        "cancer-testis-antigens-extended.csv": extended,
    }.items():
        frame.to_csv(out / name, index=False)
        print(name, len(frame))
    (out / "cta-intake-identity-snapshot.json").write_bytes(ncbi_json.read_bytes())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pubmed-xml", type=Path, required=True)
    parser.add_argument("--ncbi-json", type=Path, required=True)
    parser.add_argument("--retrieved-date", type=date.fromisoformat, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    build(args.pubmed_xml, args.ncbi_json, args.retrieved_date.isoformat(), args.out)

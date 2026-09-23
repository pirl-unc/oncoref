#!/usr/bin/env python3
"""Import Gong 2021 S5/S6 and Bradley 2020 Fig. 3 into the CTA evidence table.

Download Gong's original supplementary workbook from the URL below, then run:
  python scripts/import_cta_publications.py --gong-xlsx workbook.xlsx
Only new, currently protein-coding genes receive freshly computed HPA columns.
Existing evidence is preserved; neither paper bypasses the default filters.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd

from oncoref.cta_sources import (
    BRADLEY,
    BRADLEY_GENES,
    GONG_NC,
    GONG_PC,
    add_gene_evidence_tags,
    add_publication_candidates,
    extract_publications,
)
from oncoref.gene_ids import canonical_gene_space
from oncoref.load_dataset import get_data

GONG_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41467-021-22695-y/MediaObjects/41467_2021_22695_MOESM2_ESM.xlsx"
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gong-xlsx", type=Path, required=True)
    parser.add_argument(
        "--out-dir", type=Path, default=Path(__file__).resolve().parents[1] / "oncoref/data"
    )
    args = parser.parse_args()
    membership = extract_publications(args.gong_xlsx)
    original = get_data("cancer-testis-antigens")
    candidates = add_gene_evidence_tags(add_publication_candidates(original, membership))
    canonical_hash = hashlib.sha256(canonical_gene_space().to_csv(index=False).encode()).hexdigest()
    workbook_hash = hashlib.sha256(args.gong_xlsx.read_bytes()).hexdigest()
    bradley_hash = hashlib.sha256(("\n".join(BRADLEY_GENES) + "\n").encode()).hexdigest()
    sources = []
    for tag, table, count in (
        (GONG_PC, "Supplementary Data 5: placenta rows", 71),
        (GONG_NC, "Supplementary Data 6: placenta rows", 74),
        (BRADLEY, "Figure 3: ten CPA candidates", 10),
    ):
        gong = tag != BRADLEY
        sources.append(
            {
                "source_tag": tag,
                "citation": "Gong et al., Nature Communications (2021)"
                if gong
                else "Bradley et al., Nature Communications (2020)",
                "title": (
                    "The RNA landscape of the human placenta in health and disease"
                    if gong
                    else "Vestigial-like 1 is a shared targetable cancer-placenta antigen expressed by pancreatic and basal-like breast cancers"
                ),
                "doi": "10.1038/s41467-021-22695-y" if gong else "10.1038/s41467-020-19141-w",
                "source_url": GONG_URL
                if gong
                else "https://www.nature.com/articles/s41467-020-19141-w/figures/3",
                "source_table": table,
                "published_rows": count,
                "input_sha256": workbook_hash if gong else bradley_hash,
                "hash_scope": "Original XLSX bytes"
                if gong
                else "Transcribed Figure 3 symbols, listed order, UTF-8 one per line with final newline",
                "evidence_scope": "Placental RNA enrichment; no tumor-antigen validation"
                if gong
                else "CPA expression nominations; HLA peptide and antigen-specific T-cell validation for VGLL1",
                "canonical_reference_sha256": canonical_hash,
                "canonical_ensembl_release": ";".join(
                    sorted(canonical_gene_space().ensembl_release.astype(str).unique())
                ),
                "retrieved_date": "2026-09-23",
            }
        )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    membership.to_csv(args.out_dir / "cta-publication-membership.csv", index=False)
    pd.DataFrame(sources).to_csv(args.out_dir / "cta-publication-sources.csv", index=False)
    candidates.to_csv(args.out_dir / "cancer-testis-antigens.csv", index=False)
    print(
        f"{len(membership)} publication memberships; {len(original)} -> {len(candidates)} candidate genes"
    )
    print(
        membership.groupby("source_tag")
        .agg(published=("source_symbol", "size"), coding=("protein_candidate_eligible", "sum"))
        .to_string()
    )


if __name__ == "__main__":
    main()

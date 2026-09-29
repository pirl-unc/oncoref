#!/usr/bin/env python3
"""Import Gong 2021 S5/S6 and Bradley 2020 Fig. 3 into the CTA evidence table.

Download Gong's original supplementary workbook from the URL below, then run:
  python scripts/import_cta_publications.py --gong-xlsx workbook.xlsx --retrieved-date YYYY-MM-DD
Add --apply to overwrite the bundled tables.
Only new, currently protein-coding genes receive freshly computed core HPA
columns. Existing core evidence is preserved; the paired extended snapshot is
regenerated from the same candidates. Neither paper bypasses the default filters.

Safe by default: writes ``*.import.csv`` sidecars and prints the intake report,
but does NOT overwrite the bundled CSVs unless ``--apply`` is passed.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from datetime import date
from pathlib import Path

import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))

from oncoref.cta_regen import regenerate_cta_columns  # noqa: E402
from oncoref.cta_sources import (  # noqa: E402
    BRADLEY,
    BRADLEY_GENES,
    GONG_NC,
    GONG_PC,
    add_gene_evidence_tags,
    add_publication_candidates,
    extract_publications,
    publication_membership,
    publication_sources,
)
from oncoref.gene_ids import canonical_gene_space  # noqa: E402
from oncoref.load_dataset import get_data  # noqa: E402

GONG_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41467-021-22695-y/MediaObjects/41467_2021_22695_MOESM2_ESM.xlsx"
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gong-xlsx", type=Path, required=True)
    parser.add_argument(
        "--retrieved-date",
        type=date.fromisoformat,
        required=True,
        help="Actual source retrieval/verification date (YYYY-MM-DD), not the import date.",
    )
    parser.add_argument("--out-dir", type=Path, default=_REPO_ROOT / "oncoref/data")
    parser.add_argument("--apply", action="store_true", help="Overwrite the bundled CSVs in place.")
    args = parser.parse_args()
    imported = extract_publications(args.gong_xlsx)
    tags = {GONG_PC, GONG_NC, BRADLEY}
    prior = publication_membership()
    membership = pd.concat([imported, prior[~prior.source_tag.isin(tags)]], ignore_index=True)
    original = get_data("cancer-testis-antigens")
    candidates = add_gene_evidence_tags(add_publication_candidates(original, membership))
    extended = regenerate_cta_columns(candidates, tissue_scope="extended")
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
                "retrieved_date": args.retrieved_date.isoformat(),
            }
        )
    prior_sources = publication_sources()
    sources = pd.concat(
        [pd.DataFrame(sources), prior_sources[~prior_sources.source_tag.isin(tags)]],
        ignore_index=True,
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    outputs = (
        ("cta-publication-membership.csv", membership),
        ("cta-publication-sources.csv", pd.DataFrame(sources)),
        ("cancer-testis-antigens.csv", candidates),
        ("cancer-testis-antigens-extended.csv", extended),
    )
    for name, frame in outputs:
        shipped = args.out_dir / name
        dest = shipped if args.apply else shipped.with_suffix(".import.csv")
        frame.to_csv(dest, index=False)
        print(f"{'wrote' if args.apply else '(dry run) wrote'} {len(frame)} rows to {dest}")
    print(
        f"{len(membership)} publication memberships; {len(original)} -> {len(candidates)} candidate genes"
    )
    print(
        membership.groupby("source_tag")
        .agg(published=("source_symbol", "size"), coding=("protein_candidate_eligible", "sum"))
        .to_string()
    )
    if not args.apply:
        print("\nShipped CSVs untouched. Re-run with --apply to overwrite them.")


if __name__ == "__main__":
    main()

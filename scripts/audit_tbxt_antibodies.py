#!/usr/bin/env python3
"""Extract the pinned evidence behind the TBXT antibody discrepancy.

This is a source comparison, not a classification override. Whole-organ negatives,
cell-type RNA, and exact IHC compartments remain separately identifiable.
"""

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import hpa, normal_cell_data, normal_cell_screen, reference_data  # noqa: E402

GENE = "ENSG00000164458"
OUT = ROOT / "docs/audits/brachyury-antibody-20261027"


def main() -> None:
    """Run the command-line audit with explicit input and output paths."""
    OUT.mkdir(parents=True, exist_ok=True)
    receipts = []
    sources = [
        "hpa_antibody_xml",
        "hpa_rna_consensus",
        "hpa_normal_tissue",
        "hpa_single_cell_tissue",
        "hpa_immune_cell",
        "hpa_immune_monaco",
    ]
    for source in sources:
        if source in {"hpa_rna_consensus", "hpa_normal_tissue"}:
            data = hpa._read_hpa(source, "v23")
            data = data[data.Gene.eq(GENE)]
        else:
            data = normal_cell_data.table(source, genes=[GENE])
        if source == "hpa_antibody_xml":
            data = normal_cell_screen.antibody_findings(data)
            data = data.drop(columns=["Antigen sequence", "Image URLs"], errors="ignore")
        data.to_csv(OUT / f"{source}.csv", index=False)
        p = reference_data.provenance(source, verify_content=True)
        receipts.append({k: p[k] for k in ["name", "version", "url", "sha256"]})
        print(source, len(data), flush=True)
    peaks = pd.read_csv(
        ROOT / "docs/audits/cta-expanded-funnel-20261027/normal/normal-peak-evidence.csv.gz"
    )
    peaks[peaks.gene_id.eq(GENE)].to_csv(OUT / "funnel-rna-peaks.csv", index=False)
    manifest = {
        "gene_id": GENE,
        "symbol": "TBXT",
        "review_date": "2026-10-27",
        "classification_override": False,
        "sources": receipts,
        "outputs": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob("*.csv")},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

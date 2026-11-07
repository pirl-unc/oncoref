#!/usr/bin/env python3
"""Extract native localization annotations from the pinned Swiss-Prot download.

Usage: .venv/bin/python scripts/import_cta_localization.py /path/to/pinned.json.gz
The report rebuild uses the committed extract and needs no network or cache.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def extract(path):
    source = json.loads((ROOT / "oncoref/data/protein-annotation-sources.json").read_text())[
        "artifacts"
    ]["uniprot_reviewed_human"]
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != source["sha256"]:
        raise ValueError("Input does not match the pinned UniProt snapshot")
    guide = json.loads(
        gzip.decompress(
            (ROOT / "docs/audits/cta-evidence-guide-20261106/data.json.gz").read_bytes()
        )
    )
    wanted = {g for r in guide["rows"] if r["category"] == "CTA" for g in r["gene_ids"]}
    entries = []
    for entry in json.loads(gzip.decompress(raw))["results"]:
        ids = sorted(
            {
                p["value"].split(".")[0]
                for x in entry.get("uniProtKBCrossReferences", [])
                if x["database"] == "Ensembl"
                for p in x.get("properties", [])
                if p["key"] == "GeneId" and p["value"].startswith("ENSG")
            }
        )
        if not wanted.intersection(ids):
            continue
        description = entry.get("proteinDescription", {})
        name = description.get("recommendedName", description.get("submissionNames", [{}])[0])
        entries.append(
            {
                "accession": entry["primaryAccession"],
                "gene_ids": ids,
                "protein_name": name.get("fullName", {}).get("value", ""),
                "protein_existence": entry.get("proteinExistence", ""),
                "comments": [
                    c
                    for c in entry.get("comments", [])
                    if c["commentType"] in {"SUBCELLULAR LOCATION", "CAUTION"}
                ],
            }
        )
    return {
        "schema": "oncoref.cta-localization-input",
        "version": 1,
        "source": source,
        "mapping": "Explicit Ensembl GeneId cross-references only; no symbol or paralog transfer.",
        "entries": sorted(entries, key=lambda e: e["accession"]),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    args = parser.parse_args()
    result = extract(args.input)
    output = ROOT / "scripts/cta_evidence_assets/cta-localization-input.json.gz"
    output.write_bytes(
        gzip.compress((json.dumps(result, ensure_ascii=False, indent=2) + "\n").encode(), mtime=0)
    )
    print(f"{len(result['entries'])} mapped entries → {output}")

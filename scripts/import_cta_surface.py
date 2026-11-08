#!/usr/bin/env python3
"""Pin already-retrieved UniProt features and local HPA prevalence for 37 CTAs.

The normal report build uses this compact input offline. Acquisition of the
UniProt JSON/HTTP receipt and the verified sequence bundle is explicit.
"""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from oncoref import cta_literature, hpa_cancer, sequence_bundle

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uniprot-json", type=Path, required=True)
    parser.add_argument("--uniprot-receipt", type=Path, required=True)
    parser.add_argument("--sequence-bundle", type=Path, required=True)
    args = parser.parse_args()
    guide = json.loads(
        gzip.decompress(
            (ROOT / "docs/audits/cta-evidence-guide-20261106/data.json.gz").read_bytes()
        )
    )
    locations = cta_literature.localization().to_dict("records")
    locations = [r for r in locations if "cell_surface" in r["facets"]]
    ids = {r["entity_id"] for r in locations}
    rows = [r for r in guide["rows"] if r["entity_id"] in ids]
    hashes = {h for r in rows for h in r["representative_sequence_sha256s"]}
    bundle = sequence_bundle.load(args.sequence_bundle)
    sequences = {}
    for record in bundle.records():
        if record["sequence_sha256"] in hashes:
            sequences[record["sequence_sha256"]] = record["sequence"]
    if set(sequences) != hashes:
        raise ValueError("Missing selected CTA sequence")
    genes = sorted({g for r in rows for g in r["gene_ids"]})
    native = json.loads(args.uniprot_json.read_bytes())
    receipt = json.loads(args.uniprot_receipt.read_bytes())
    native_sha = hashlib.sha256(args.uniprot_json.read_bytes()).hexdigest()
    if receipt["sha256"] != native_sha or set(receipt["accessions"]) != {
        r["primaryAccession"] for r in native["results"]
    }:
        raise ValueError("UniProt response does not match its retrieval receipt")
    if len(native["results"]) != 37:
        raise ValueError("Expected the 37 mapped UniProt entries")
    entries = []
    for entry in native["results"]:
        entries.append(
            {
                k: entry.get(k)
                for k in [
                    "primaryAccession",
                    "entryAudit",
                    "genes",
                    "sequence",
                    "features",
                    "comments",
                ]
            }
        )
    data = {
        "schema": "oncoref.cta-surface-input",
        "version": 1,
        "uniprot_receipt": receipt,
        "uniprot_native_sha256": native_sha,
        "uniprot_entries": entries,
        "sequence_bundle_id": bundle.identity,
        "reference": bundle.manifest["reference"],
        "selected_sequences": sequences,
        "hpa_source": hpa_cancer.hpa_cancer_sources(),
        "rna": json.loads(hpa_cancer.hpa_cancer_rna_prevalence(genes).to_json(orient="records")),
        "ihc": json.loads(hpa_cancer.hpa_cancer_ihc_prevalence(genes).to_json(orient="records")),
    }
    # pandas floating missing values must not become non-standard JSON NaNs.
    payload = json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    target = ROOT / "scripts/cta_evidence_assets/cta-surface-input.json.gz"
    target.write_bytes(gzip.compress(payload, mtime=0))
    print(
        f"Wrote {target}: {len(entries)} entries, {len(data['rna'])} RNA and {len(data['ihc'])} IHC rows"
    )


if __name__ == "__main__":
    main()

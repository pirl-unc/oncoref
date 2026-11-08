#!/usr/bin/env python3
"""Build source-bound CDS qualification from verified Ensembl GTF and peptide FASTA."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref.translation_qualification import POLICY_ID, qualify  # noqa: E402


def digest(path):
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def build(source):
    meta = json.loads((source / "reference.json").read_text())
    reference, paths = meta["reference"], meta["paths"]
    for role in ("gtf", "protein"):
        if digest(Path(paths[role])) != reference["content_provenance"][role]:
            raise ValueError(f"Source hash mismatch: {role}")
    transcripts = {}
    with gzip.open(paths["gtf"], "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[2] not in {"transcript", "CDS", "start_codon", "stop_codon"}:
                continue
            pairs = re.findall(r'(\w+) "([^"]*)";', fields[8])
            attrs = dict(pairs)
            tid = attrs["transcript_id"]
            record = transcripts.setdefault(
                tid,
                {
                    "transcript_id": tid,
                    "gene_id": attrs["gene_id"],
                    "gene_version": attrs.get("gene_version", ""),
                    "transcript_version": attrs.get("transcript_version", ""),
                    "symbol": attrs.get("gene_name", ""),
                    "contig": fields[0],
                    "strand": fields[6],
                    "gene_biotype": attrs.get("gene_biotype", ""),
                    "transcript_biotype": attrs.get("transcript_biotype", ""),
                    "tags": [],
                    "cds_bases": 0,
                    "start_codon_bases": 0,
                    "stop_codon_bases": 0,
                },
            )
            record["tags"] = sorted(set(record["tags"]) | {v for k, v in pairs if k == "tag"})
            if fields[2] == "transcript":
                record.update(locus_start=int(fields[3]), locus_end=int(fields[4]))
                for key in ("transcript_support_level", "ccds_id", "transcript_source"):
                    record[key] = attrs.get(key, "")
            else:
                key = {
                    "CDS": "cds_bases",
                    "start_codon": "start_codon_bases",
                    "stop_codon": "stop_codon_bases",
                }[fields[2]]
                record[key] += int(fields[4]) - int(fields[3]) + 1
    rows = []

    def add(header, pieces):
        fields = dict(re.findall(r"(\S+?):(\S+)", header))
        tid = fields["transcript"].split(".")[0]
        record = dict(transcripts[tid])
        sequence = "".join(pieces)
        record.update(
            protein_id=header.split()[0],
            transcript_accession=fields["transcript"],
            gene_accession=fields["gene"],
            length=len(sequence),
            sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
            reference_id=reference["reference_id"],
        )
        record.update(qualify(sequence, record))
        rows.append(record)

    with gzip.open(paths["protein"], "rt") as handle:
        header, pieces = None, []
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    add(header, pieces)
                header, pieces = line[1:].strip(), []
            else:
                pieces.append(line.strip())
        if header is not None:
            add(header, pieces)
    table = pd.DataFrame(rows).sort_values("protein_id")
    if table.protein_id.duplicated().any():
        raise ValueError("Duplicate protein accession")
    for column in ("tags", "qualification_reasons"):
        table[column] = table[column].map(json.dumps)
    data = ROOT / "oncoref/data"
    output = data / "translation-qualification.csv.gz"
    output.write_bytes(
        gzip.compress(table.to_csv(index=False, lineterminator="\n").encode(), mtime=0)
    )
    receipt = {
        "schema": "oncoref.translation-qualification",
        "schema_version": 1,
        "policy_id": POLICY_ID,
        "reference": reference,
        "orfdata_identity": meta["orfdata_identity"],
        "records_sha256": digest(output),
        "records": len(table),
        "qualified_records": int(table.target_qualified.sum()),
        "reason_counts": dict(
            Counter(r for value in table.qualification_reasons for r in json.loads(value))
        ),
        "inputs": {
            str(p.relative_to(ROOT)): digest(p)
            for p in [Path(__file__), ROOT / "oncoref/translation_qualification.py"]
        },
        "limitations": [
            "Annotated translation competence, not experimental protein production",
            "ncORF support is not transferred to a different reference CDS",
            "Alternate loci retained but cannot rescue primary candidates",
        ],
    }
    (data / "translation-qualification.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({k: receipt[k] for k in ["records", "qualified_records", "reason_counts"]}))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    build(p.parse_args().source)

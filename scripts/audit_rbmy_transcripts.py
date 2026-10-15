#!/usr/bin/env python3
"""Extract the primary annotation evidence for the two CTA corrections in #576.

Read existing Ensembl 112/115 GTF and peptide downloads without pyensembl or
network access. This audit does not select transcripts by protein length and
does not modify the CTA tables.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path

GENES = {"ENSG00000169800": "RBMY1F", "ENSG00000226941": "RBMY1J"}
RELEASES = (112, 115)


def source_record(path: Path, release: int, kind: str) -> dict:
    directory = "gtf/homo_sapiens" if kind == "gtf" else "fasta/homo_sapiens/pep"
    return {
        "filename": path.name,
        "url": f"https://ftp.ensembl.org/pub/release-{release}/{directory}/{path.name}",
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def extract_reference(root: Path, release: int) -> dict:
    directory = root / f"ensembl{release}"
    gtf = directory / f"Homo_sapiens.GRCh38.{release}.chr_patch_hapl_scaff.gtf.gz"
    fasta = directory / "Homo_sapiens.GRCh38.pep.all.fa.gz"
    transcripts = {}
    with gzip.open(gtf, "rt") as stream:
        for line in stream:
            if line.startswith("#") or not any(gene in line for gene in GENES):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "transcript":
                continue
            attributes = re.findall(r'(\w+) "([^"]*)";', fields[8])
            attrs = dict(attributes)
            if attrs["gene_id"] not in GENES:
                continue
            tid = attrs["transcript_id"]
            if tid in transcripts:
                raise ValueError(f"Duplicate transcript {tid} in {gtf.name}")
            transcripts[tid] = {
                "gene_id": attrs["gene_id"],
                "symbol": attrs["gene_name"],
                "transcript_id": tid,
                "transcript_version": attrs["transcript_version"],
                "tags": [value for key, value in attributes if key == "tag"],
                "gtf_line": line.rstrip("\n"),
                "peptide": None,
            }

    def retain_peptide(header, sequence):
        fields = dict(re.findall(r"(\w+):([^\s]+)", header))
        tid = fields.get("transcript", "").split(".")[0]
        if tid not in transcripts:
            return
        row = transcripts[tid]
        if fields["gene"].split(".")[0] != row["gene_id"] or row["peptide"] is not None:
            raise ValueError(f"Inconsistent peptide association for {tid}")
        protein = "".join(sequence)
        row["peptide"] = {
            "protein_id": header.split()[0][1:],
            "fasta_header": header,
            "length": len(protein),
            "sequence": protein,
            "sha256": hashlib.sha256(protein.encode("ascii")).hexdigest(),
        }

    header, sequence = "", []
    with gzip.open(fasta, "rt") as stream:
        for line in stream:
            if line.startswith(">"):
                retain_peptide(header, sequence)
                header, sequence = line.rstrip("\n"), []
            else:
                sequence.append(line.strip())
        retain_peptide(header, sequence)
    for gid in GENES:
        canonical = [
            row
            for row in transcripts.values()
            if row["gene_id"] == gid and "Ensembl_canonical" in row["tags"]
        ]
        if len(canonical) != 1 or canonical[0]["peptide"] is None:
            raise ValueError(f"Expected one translated Ensembl canonical transcript for {gid}")
    return {
        "ensembl_release": release,
        "sources": {
            "gtf": source_record(gtf, release, "gtf"),
            "peptide_fasta": source_record(fasta, release, "peptide"),
        },
        "transcripts": [transcripts[tid] for tid in sorted(transcripts)],
    }


def build_audit(root: Path) -> dict:
    return {
        "issue": "https://github.com/pirl-unc/oncoref/issues/576",
        "taxon_id": 9606,
        "assembly": "GRCh38",
        "selection_rule": "Ensembl_canonical transcript tag; MANE_Select is corroborating evidence. No longest-protein selection or expression allocation.",
        "previous_cta_transcript_ids": {
            "ENSG00000169800": "ENST00000454978",
            "ENSG00000226941": None,
        },
        "references": [extract_reference(root, release) for release in RELEASES],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(build_audit(args.reference_root), indent=2) + "\n")

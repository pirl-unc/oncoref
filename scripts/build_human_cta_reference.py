#!/usr/bin/env python3
"""Build human CTA annotation mappings from pinned Ensembl GTF/FASTA bytes.

No network or annotation indexes. This maps annotation identities, not historical
assay targets or measured isoform abundance. See docs/species-cta-evidence.md.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

from oncoref.reference_identity import ReferenceIdentity, _canonical_json, _identity
from oncoref.species_assets import registered_references


def _attributes(text):
    return dict(re.findall(r'(\w+) "([^"]*)"', text))


def _versioned(attrs, kind):
    return f"{attrs[kind + '_id']}.{attrs[kind + '_version']}"


def annotation_records(path, gene_ids):
    """Retain all loci and transcript/protein alternatives for exact gene IDs."""
    genes = defaultdict(lambda: {"loci": set(), "transcripts": {}, "symbols": set()})
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] not in {"gene", "transcript", "CDS"}:
                continue
            attrs = _attributes(fields[8])
            gid = attrs["gene_id"]
            if gid not in gene_ids:
                continue
            gene = genes[gid]
            gene["symbols"].add(attrs.get("gene_name", gid))
            if fields[2] == "gene":
                gene["loci"].add(
                    (
                        _versioned(attrs, "gene"),
                        fields[0],
                        int(fields[3]),
                        int(fields[4]),
                        fields[6],
                    )
                )
                continue
            tid = _versioned(attrs, "transcript")
            tx = gene["transcripts"].setdefault(tid, {"genes": set(), "proteins": set()})
            tx["genes"].add(_versioned(attrs, "gene"))
            if fields[2] == "CDS":
                tx["proteins"].add(_versioned(attrs, "protein"))
    return dict(genes)


def fasta_records(paths, identifiers):
    """Verify complete versioned IDs and metadata; conflicting duplicates fail."""
    records = {}

    def retain(header, chunks):
        if header is None or header.split()[0] not in identifiers:
            return
        key = header.split()[0]
        attrs = dict(part.split(":", 1) for part in header.split()[1:] if ":" in part)
        sequence = "".join(chunks)
        row = {
            "gene": attrs.get("gene"),
            "transcript": attrs.get("transcript"),
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "length": len(sequence),
        }
        if not sequence or (key in records and records[key] != row):
            raise ValueError(f"Empty or conflicting FASTA identity: {key}")
        records[key] = row

    for path in paths:
        with gzip.open(path, "rt") as handle:
            header, chunks = None, []
            for line in handle:
                if line.startswith(">"):
                    retain(header, chunks)
                    header, chunks = line[1:].strip(), []
                elif header is not None and header.split()[0] in identifiers:
                    chunks.append(line.strip())
            retain(header, chunks)
    return records


def build(reference, paths, gene_ids):
    if (
        reference.taxon_id != 9606
        or reference.annotation_source != "Ensembl"
        or reference.annotation_release != "115"
        or reference.assembly_accession != "GCA_000001405.29"
        or set(paths) != {"gtf", "rna", "ncrna", "protein"}
        or set(paths) != set(reference.content_provenance)
    ):
        raise ValueError("Human mapping builder requires the complete Ensembl 115 reference")
    for role, path in paths.items():
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        if digest.hexdigest() != reference.content_provenance[role]:
            raise ValueError(f"Reference checksum mismatch: {role}")
    genes = annotation_records(paths["gtf"], set(gene_ids))
    tids = {tid for gene in genes.values() for tid in gene["transcripts"]}
    pids = {
        pid
        for gene in genes.values()
        for tx in gene["transcripts"].values()
        for pid in tx["proteins"]
    }
    rna = fasta_records([paths["rna"], paths["ncrna"]], tids)
    protein = fasta_records([paths["protein"]], pids)
    mappings, sequence_audit = {}, {}
    for gid in sorted(gene_ids):
        gene = genes.get(gid)
        pairs, missing = [], []
        if gene:
            for tid, tx in sorted(gene["transcripts"].items()):
                native_rna = rna.get(tid)
                if native_rna is None:
                    missing.append(tid)
                elif native_rna["gene"] not in tx["genes"]:
                    raise ValueError(f"RNA/GTF gene identity conflict: {tid}")
                for pid in sorted(tx["proteins"]) or [None]:
                    native_protein = protein.get(pid)
                    if pid is not None:
                        if native_protein is None:
                            missing.append(pid)
                        elif (
                            native_protein["gene"] not in tx["genes"]
                            or native_protein["transcript"] != tid
                        ):
                            raise ValueError(f"Protein/GTF identity conflict: {pid}")
                    # Keep annotation IDs even when bytes are absent. Such a
                    # mapping is explicitly incomplete and cannot be admitted
                    # as a unique complete reference mapping by a consumer.
                    pairs.append({"transcript_id": tid, "protein_id": pid})
                    sequence_audit[tid] = native_rna
                    if pid:
                        sequence_audit[pid] = native_protein
        loci = sorted(gene["loci"]) if gene else []
        complete = bool(loci and pairs) and not missing
        status = "unresolved" if not gene else "unique" if len(loci) == 1 else "ambiguous"
        mappings[gid] = {
            "reference_id": reference.reference_id,
            "gene_namespace": "Ensembl gene",
            "status": status,
            "method": "Exact Ensembl gene ID; all versioned GTF transcript/protein pairs verified against pinned RNA/peptide FASTA headers and bytes",
            "complete": complete,
            "targets": [{"id": gid, "symbol": sorted(gene["symbols"])[0], "sequences": pairs}]
            if gene
            else [],
            "notes": "Annotation inventory, not measured isoform abundance. No preferred isoform is selected. "
            "Loci are [versioned gene ID, contig, 1-based inclusive start, end, strand]: "
            + _canonical_json(loci)
            + ". Missing sequence assets: "
            + _canonical_json(sorted(set(missing)))
            + ". Null protein IDs mean no annotated CDS, not an unmeasured protein inferred from RNA.",
        }
    payload = {
        "schema": "oncoref.cta-reference-mappings",
        "schema_version": 1,
        "reference": reference.as_dict(),
        "mappings": mappings,
        "sequence_audit": sequence_audit,
    }
    return {**payload, "mapping_id": _identity("cta-reference-mappings", 1, payload)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "paths", type=Path, help="JSON mapping gtf/rna/ncrna/protein to local files"
    )
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    from oncoref.cta import cta_evidence

    registration = next(r for r in registered_references("human") if r["id"] == "ensembl115")
    result = build(
        ReferenceIdentity.from_dict(registration["identity"]),
        json.loads(args.paths.read_text()),
        set(cta_evidence().Ensembl_Gene_ID),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(gzip.compress(_canonical_json(result).encode("ascii"), mtime=0))
    print(
        json.dumps(
            {
                "mapping_id": result["mapping_id"],
                "genes": len(result["mappings"]),
                "sequence_records": len(result["sequence_audit"]),
                "bytes": args.output.stat().st_size,
            }
        )
    )


if __name__ == "__main__":
    main()

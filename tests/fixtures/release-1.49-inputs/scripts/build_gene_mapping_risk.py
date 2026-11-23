#!/usr/bin/env python3
"""Annotate genes whose sparse RNA counts can come from another locus.

A count matrix cannot show where a read aligned. Two reference properties make
a mis-assigned molecule plausible, and both are computed here from Ensembl:

- shared sequence: sense-strand k-mers of the gene's transcripts that also occur
  in transcripts of a gene, pseudogene or noncoding RNA at another locus;
- host overlap: an exon lying inside another gene on the same strand, so the
  host's unspliced RNA is credited to the nested gene when introns are counted.

Members of one protein group are not each other's confounders.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import pandas as pd
from pyensembl import EnsemblRelease

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "oncoref/data"
RELEASE = 115
K = 31
TAIL = 1000
#: A tenth of the transcript, or of its 3' kilobase, shared with another locus.
SHARED_FRACTION = 0.10
PRIMARY = {*map(str, range(1, 23)), "X", "Y", "MT"}


def kmers(sequence: str, k: int = K):
    return {sequence[i : i + k] for i in range(len(sequence) - k + 1)}


def build(proteins: pd.DataFrame) -> pd.DataFrame:
    genome = EnsemblRelease(RELEASE)
    group = {g: p.entity_id for p in proteins.itertuples() for g in json.loads(p.member_gene_ids)}
    owners: dict[str, set[str]] = {}
    own, tail, locus = {}, {}, {}
    for gene_id in group:
        try:
            gene = genome.gene_by_id(gene_id)
        except ValueError:
            continue
        own[gene_id], tail[gene_id] = set(), set()
        locus[gene_id] = (gene.contig, gene.strand, gene.start, gene.end)
        for transcript in gene.transcripts:
            sequence = transcript.sequence
            if not sequence:
                continue
            own[gene_id] |= kmers(sequence)
            tail[gene_id] |= kmers(sequence[-TAIL:])
        for kmer in own[gene_id]:
            owners.setdefault(kmer, set()).add(gene_id)
    shared = {g: set() for g in own}
    partners = {g: Counter() for g in own}
    for transcript in genome.transcripts():
        # Alternate scaffolds and patches repeat primary loci; they are not paralogs.
        if transcript.contig not in PRIMARY:
            continue
        sequence = transcript.sequence
        if not sequence:
            continue
        other = transcript.gene_id
        hits = Counter()
        for i in range(len(sequence) - K + 1):
            targets = owners.get(sequence[i : i + K])
            if targets:
                for target in targets:
                    contig, strand, start, end = locus[target]
                    # Another record of the same locus is not a paralog; see host overlap.
                    same_locus = (
                        transcript.contig == contig
                        and transcript.strand == strand
                        and transcript.start <= end
                        and start <= transcript.end
                    )
                    if not same_locus and group.get(other) != group[target]:
                        shared[target].add(sequence[i : i + K])
                        hits[target] += 1
        for target, n in hits.items():
            label = f"{transcript.gene_name or other} ({transcript.biotype})"
            partners[target][label] = max(partners[target][label], n)
    rows = []
    for gene_id in group:
        if gene_id not in own:
            rows.append({"gene_id": gene_id, "annotated": False})
            continue
        gene = genome.gene_by_id(gene_id)
        exons = {(e.start, e.end) for t in gene.transcripts for e in t.exons}
        hosts = []
        for other in genome.genes_at_locus(gene.contig, gene.start, gene.end, strand=gene.strand):
            if other.gene_id == gene_id or group.get(other.gene_id) == group[gene_id]:
                continue
            if other.gene_name and other.gene_name == gene.gene_name:
                continue  # A duplicate record of this gene is not a host.
            if any(other.start <= start and end <= other.end for start, end in exons):
                hosts.append(f"{other.gene_name or other.gene_id} ({other.biotype})")
        n_own, n_tail = len(own[gene_id]), len(tail[gene_id])
        fraction = len(shared[gene_id]) / n_own if n_own else None
        tail_fraction = len(shared[gene_id] & tail[gene_id]) / n_tail if n_tail else None
        ambiguous = bool(n_own) and max(fraction, tail_fraction) >= SHARED_FRACTION
        rows.append(
            {
                "gene_id": gene_id,
                "annotated": True,
                "symbol": gene.gene_name,
                "transcript_kmers": n_own,
                "shared_kmer_fraction": fraction,
                "shared_3prime_kmer_fraction": tail_fraction,
                "sequence_partners": "; ".join(n for n, _ in partners[gene_id].most_common(5)),
                "same_strand_hosts": "; ".join(sorted(hosts)),
                "paralog_ambiguous": ambiguous,
                "host_overlapped": bool(hosts),
                "mapping_risk": ambiguous or bool(hosts),
            }
        )
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", type=Path, required=True)
    args = parser.parse_args()
    proteins = pd.read_csv(args.proteins)
    table = build(proteins).sort_values("gene_id")
    path = DATA / "gene-mapping-risk.csv.gz"
    table.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    provenance = {
        "schema": "oncoref.gene-mapping-risk",
        "ensembl_release": RELEASE,
        "k": K,
        "three_prime_window": TAIL,
        "shared_fraction_flag": SHARED_FRACTION,
        "strand": "sense k-mers only; the assays are stranded",
        "genes": len(table),
        "unannotated": int((~table.annotated).sum()),
        "paralog_ambiguous": int(table.paralog_ambiguous.fillna(False).sum()),
        "host_overlapped": int(table.host_overlapped.fillna(False).sum()),
        "mapping_risk": int(table.mapping_risk.fillna(False).sum()),
        "limitations": [
            "Exact shared k-mers; diverged or unannotated paralogous sequence is not found",
            "A flag marks a possible confounder, not a demonstrated mis-assignment",
            "An unflagged gene can still collect intronic or genomic-DNA reads",
        ],
        "sha256": digest,
    }
    (DATA / "gene-mapping-risk.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({k: v for k, v in provenance.items() if isinstance(v, int)}))


if __name__ == "__main__":
    main()

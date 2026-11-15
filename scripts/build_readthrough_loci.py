"""Find read-through loci among the catalog's proteins from the gene annotation.

A locus is a read-through when every one of its protein-coding transcripts shares
exons with two other protein-coding genes on the same strand that do not overlap each
other. Requiring every transcript leaves out ordinary genes that merely have one
read-through isoform. Requiring disjoint parents leaves out clustered families that
share constant exons (protocadherins, UGT1A): their members all overlap one another.
A read-through protein is its parents' sequence joined; it has no expression or
biology of its own, and reads assigned to it belong to the parents.

Writes oncoref/data/readthrough-loci.csv, which the catalog's sequence check reads.
Usage: PYTHONPATH=. python scripts/build_readthrough_loci.py
"""

import json
from pathlib import Path

import pandas as pd
from pyensembl import EnsemblRelease

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz"
OUTPUT = ROOT / "oncoref/data/readthrough-loci.csv"
ENSEMBL_RELEASE = 115
#: Loci named as read-throughs by NCBI Gene that the structural rule does not reach.
NAMED = {"ZNF816-ZNF321P": "NCBI Gene name: ZNF816-ZNF321P readthrough"}


def exons(transcripts):
    return {(e.start, e.end) for t in transcripts for e in t.exons}


def parents(genome, gene, transcripts):
    """Other same-strand protein-coding genes whose exons these transcripts share."""
    mine, found = exons(transcripts), {}
    for other in genome.genes_at_locus(gene.contig, gene.start, gene.end, strand=gene.strand):
        if other.gene_id == gene.gene_id or other.biotype != "protein_coding":
            continue
        theirs = exons(t for t in other.transcripts if t.biotype == "protein_coding")
        shared = sum(any(a[0] <= b[1] and b[0] <= a[1] for b in theirs) for a in mine)
        if shared:
            found[other.gene_id] = (
                other.gene_name or other.gene_id,
                shared,
                other.start,
                other.end,
            )
    return found


def joins_two_genes(found) -> bool:
    spans = [(s, e) for _, _, s, e in found.values()]
    return any(a[1] < b[0] or b[1] < a[0] for a in spans for b in spans)


def build():
    genome = EnsemblRelease(ENSEMBL_RELEASE)
    proteins = pd.read_csv(
        BASE, usecols=["protein_equivalence_group_id", "protein_group_label", "member_gene_ids"]
    )
    rows = []
    for r in proteins.itertuples():
        for gene_id in json.loads(r.member_gene_ids):
            if r.protein_group_label in NAMED:
                rows.append(
                    (
                        r.protein_equivalence_group_id,
                        r.protein_group_label,
                        gene_id,
                        "",
                        NAMED[r.protein_group_label],
                    )
                )
                continue
            try:
                gene = genome.gene_by_id(gene_id)
            except ValueError:
                continue
            coding = [t for t in gene.transcripts if t.biotype == "protein_coding"]
            # A gene with one read-through isoform among ordinary ones is not a read-through locus.
            if not coding or not all(joins_two_genes(parents(genome, gene, [t])) for t in coding):
                continue
            found = parents(genome, gene, coding)
            text = "; ".join(f"{n} ({k} exons)" for n, k, _, _ in sorted(found.values()))
            basis = (
                f"Ensembl {ENSEMBL_RELEASE}: every coding transcript joins two separate "
                f"protein-coding genes on the same strand; {len(found)} parents in all"
            )
            rows.append(
                (r.protein_equivalence_group_id, r.protein_group_label, gene_id, text, basis)
            )
    columns = ["entity_id", "protein_group_label", "gene_id", "parent_genes", "basis"]
    out = pd.DataFrame(rows, columns=columns).sort_values(["protein_group_label", "gene_id"])
    out.to_csv(OUTPUT, index=False)
    return out


if __name__ == "__main__":
    print(f"{len(build())} read-through loci")

#!/usr/bin/env python3
"""Build the CTA protein universe: coding genes plus established ncORFs.

Every classified protein must be translated from a coding transcript or from an
established non-canonical ORF. Established means at least two independent
Ribo-seq studies or any mass-spectrometry peptide support in orfdata's catalog.
Genes whose only reference CDS comes from NMD or CDS-not-defined transcripts,
and ORFs supported by one Ribo-seq study, remain in the table as excluded.
Reference-CDS microproteins (under 50 residues) nested inside a same-strand
protein-coding gene stay in the universe, tagged ``nested_microprotein``: their RNA is
indistinguishable from the host's and they are not established ncORFs, so the
selection policy classifies them Other:Nested and never CTA.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "oncoref/data/cta-protein-universe.csv"
CODING = {
    "protein_coding", "protein_coding_LoF",
    "IG_V_gene", "IG_C_gene", "IG_J_gene", "IG_D_gene",
    "TR_V_gene", "TR_C_gene", "TR_J_gene", "TR_D_gene",
}  # fmt: skip
MICROPROTEIN_MAX_AA = 49  # shorter than 50 residues
MS_EVIDENCE = ("immunopeptidomics", "proteomics", "ms_unclassified")


def nested_microproteins(genome, genes) -> set[str]:
    """Genes whose longest protein is a microprotein inside a protein-coding host."""
    result = set()
    for gene_id in genes:
        gene = genome.gene_by_id(gene_id)
        longest = max((len(t.protein_sequence or "") for t in gene.transcripts), default=0)
        if not 0 < longest <= MICROPROTEIN_MAX_AA:
            continue
        for host in genome.genes_at_locus(gene.contig, gene.start, gene.end, gene.strand):
            if (
                host.gene_id != gene_id
                and host.biotype == "protein_coding"
                and host.start <= gene.start
                and gene.end <= host.end
                and host.end - host.start > gene.end - gene.start
            ):
                result.add(gene_id)
                break
    return result


def build(release: int = 115, nc_orf_sources: str = "bundled") -> pd.DataFrame:
    from orfdata import ReferenceTranscriptome
    from pyensembl import EnsemblRelease

    genome = EnsemblRelease(release, species="human")
    proteome = ReferenceTranscriptome.from_genome(genome).proteome(nc_orf_sources=nc_orf_sources)
    coding, translated, established, observed = set(), set(), {}, {}
    for protein in proteome.by_accession().values():
        gene = protein.gene_id
        if protein.protein_source == "nc_orf":
            detail = protein.evidence_detail or {}
            studies = len(detail.get("ribosome_studies") or [])
            mass_spec = any(x in protein.evidence for x in MS_EVIDENCE) or bool(
                detail.get("peptides")
            )
            observed[gene] = max(observed.get(gene, 0), studies)
            if studies >= 2 or mass_spec:
                established[gene] = "mass_spectrometry" if mass_spec else "ribo_seq_2_studies"
        else:
            translated.add(gene)
            if protein.transcript_biotype in CODING:
                coding.add(gene)
    micro = nested_microproteins(genome, coding - set(established))
    rows = []
    for gene in sorted(translated | set(observed)):
        if gene in micro:
            rows.append((gene, "coding", "reference_cds", "nested_microprotein"))
        elif gene in coding:
            rows.append((gene, "coding", "reference_cds", "coding_transcript"))
        elif gene in established:
            rows.append((gene, "coding", "nc_orf", established[gene]))
        elif gene in observed:
            rows.append((gene, "excluded", "nc_orf", "single_study_ribo_seq"))
        else:
            rows.append((gene, "excluded", "reference_cds", "nmd_or_cds_not_defined_only"))
    return pd.DataFrame(rows, columns=["gene_id", "universe", "protein_source", "evidence"])


if __name__ == "__main__":
    table = build()
    table.to_csv(OUTPUT, index=False)
    print(table.groupby(["universe", "protein_source", "evidence"]).size().to_string())
    sys.exit(0)

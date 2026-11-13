"""Rebuild the packaged curation tables from the literature sweep and the annotation.

Reads docs/audits/target-filters-20261112/normal-function-sweep.csv (every finding,
verdict and source from the literature sweep) and oncoref/data/target-curation.yaml,
and writes oncoref/data/target-normal-function-findings.csv and
oncoref/data/target-readthrough-loci.csv.
Usage: PYTHONPATH=. python scripts/build_target_curation.py
"""

import json
from pathlib import Path

import pandas as pd
from pyensembl import EnsemblRelease

from oncoref import target_curation, target_selection

ROOT = Path(__file__).resolve().parents[1]
SWEEP = ROOT / "docs/audits/target-filters-20261112/normal-function-sweep.csv"
DATA = ROOT / "oncoref/data"
ENSEMBL_RELEASE = 115
#: Loci named as read-throughs by NCBI Gene that the structural rule does not reach.
NAMED_READTHROUGH = {"ZNF816-ZNF321P": "NCBI Gene name: ZNF816-ZNF321P readthrough"}


def build_findings():
    sweep = pd.read_csv(SWEEP)
    assigned = target_curation.registry()["finding_compartments"]
    tier = target_curation.compartments().set_index("compartment").tier
    confirmed = (
        sweep.verdict.eq("confirmed")
        & sweep.human_adult_normal.eq(True)
        & sweep.gene_specific.eq(True)
    )
    keys = sweep.gene.map(lambda g: assigned.get(g, []))
    sweep["compartments"] = keys.where(confirmed, keys.map(lambda _: [])).map("; ".join)
    sweep["tier"] = sweep.compartments.map(
        lambda text: target_curation.most_protective(tier[k] for k in text.split("; ") if k)
    )
    sweep["acts"] = confirmed & sweep.tier.ne("")
    columns = [
        "gene",
        "acts",
        "tier",
        "compartments",
        "tissue",
        "cell_type",
        "evidence_type",
        "claim",
        "source",
        "best_source",
        "verdict",
        "human_adult_normal",
        "gene_specific",
        "verify_reason",
    ]
    out = sweep.sort_values(["acts", "gene"], ascending=[False, True])[columns]
    out.to_csv(DATA / "target-normal-function-findings.csv", index=False)
    return out


def build_readthrough():
    genome = EnsemblRelease(ENSEMBL_RELEASE)
    table = target_selection.results()
    table = table[table.is_target_candidate.astype(bool)]

    def exons(gene):
        return {
            (e.start, e.end)
            for t in gene.transcripts
            if t.biotype == "protein_coding"
            for e in t.exons
        }

    rows = []
    for r in table.itertuples():
        members = r.member_gene_ids
        for gene_id in json.loads(members) if isinstance(members, str) else list(members):
            if r.protein_group_label in NAMED_READTHROUGH:
                rows.append(
                    (r.protein_group_label, gene_id, "", NAMED_READTHROUGH[r.protein_group_label])
                )
                continue
            try:
                gene = genome.gene_by_id(gene_id)
            except ValueError:
                continue
            mine = exons(gene)
            parents = {}
            for other in genome.genes_at_locus(
                gene.contig, gene.start, gene.end, strand=gene.strand
            ):
                if other.gene_id == gene_id or other.biotype != "protein_coding":
                    continue
                theirs = exons(other)
                shared = sum(any(a[0] <= b[1] and b[0] <= a[1] for b in theirs) for a in mine)
                if shared:
                    parents[other.gene_name or other.gene_id] = shared
            if len(parents) >= 2:
                basis = (
                    f"Ensembl {ENSEMBL_RELEASE}: exons overlap {len(parents)} other "
                    "protein-coding genes on the same strand"
                )
                text = "; ".join(f"{k} ({v} exons)" for k, v in sorted(parents.items()))
                rows.append((r.protein_group_label, gene_id, text, basis))
    out = pd.DataFrame(rows, columns=["protein_group_label", "gene_id", "parent_genes", "basis"])
    out.to_csv(DATA / "target-readthrough-loci.csv", index=False)
    return out


if __name__ == "__main__":
    loci = build_readthrough()
    found = build_findings()
    print(f"{len(loci)} read-through loci; {len(found)} findings, {int(found.acts.sum())} acting")

#!/usr/bin/env python3
"""Build disjoint gene and exact-protein exports from pinned expression and Ensembl assets."""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_scope_report import sha256, write_csv  # noqa: E402

from oncoref.gene_identity import resolve_gene_identity  # noqa: E402
from oncoref.gene_ids import canonical_gene_space  # noqa: E402
from oncoref.normal_sequence_sets import (  # noqa: E402
    GROUPS,
    SET_POLICY_ID,
    partition_genes,
    partition_proteins,
)
from oncoref.reference_identity import ReferenceIdentity  # noqa: E402
from oncoref.sequence_groups import normalize_protein_sequence  # noqa: E402

ASSET_ROOT = "https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/"
PINNED_HASHES = {
    "gtf": "02df5a2a8c7a41e513209e0780c56cf197941b6b35ef8bab3e9a6317eecd17ef",
    "proteins": "7991ab37d09986eefdb5c77032b15370ab0e99918b78f235f3230c0302851c96",
}

ATTR = re.compile(r'(\w+) "([^"]*)";')


def read_gtf(path):
    """Read authoritative gene loci and explicit translation-completeness features."""
    genes, transcripts, header = {}, defaultdict(lambda: {"start": 0, "stop": 0}), {}
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.startswith("#!"):
                key, value = line[2:].strip().split(" ", 1)
                header[key] = value
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            kind = fields[2]
            if kind not in {"gene", "transcript", "start_codon", "stop_codon"}:
                continue
            pairs = ATTR.findall(fields[8])
            attrs = dict(pairs)
            if kind == "gene":
                genes[attrs["gene_id"]] = SimpleNamespace(
                    gene_name=attrs.get("gene_name", ""),
                    contig=fields[0],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                )
            elif kind == "transcript":
                transcripts[attrs["transcript_id"]]["incomplete"] = any(
                    key == "tag" and value in {"cds_start_NF", "cds_end_NF"} for key, value in pairs
                )
            else:
                transcripts[attrs["transcript_id"]][
                    "start" if kind == "start_codon" else "stop"
                ] += int(fields[4]) - int(fields[3]) + 1
    return genes, transcripts, header


def fasta_records(path):
    header, parts = None, []
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
    if header is not None:
        yield header, "".join(parts)


class Annotation:
    species = SimpleNamespace(latin_name="homo_sapiens")
    reference_name = "GRCh38"
    release = 115

    def __init__(self, genes):
        self.genes = genes

    def gene_by_id(self, gene_id):
        if gene_id not in self.genes:
            raise ValueError(gene_id)
        return self.genes[gene_id]


def write_fasta(path, records):
    with (
        path.open("wb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream,
    ):
        for name, sequence in records:
            stream.write((">" + name + "\n").encode())
            for start in range(0, len(sequence), 80):
                stream.write((sequence[start : start + 80] + "\n").encode())


def build(out, categories, gtf, proteins):
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((categories / "run-manifest.json").read_text())
    for name, path in [("gtf", gtf), ("proteins", proteins)]:
        if sha256(path) != PINNED_HASHES[name]:
            raise ValueError(f"{name} does not match the pinned Ensembl 115 artifact")
    space = canonical_gene_space()
    space = space.loc[space.biotype.eq("protein_coding")]
    import hashlib

    if (
        hashlib.sha256(space.to_csv(index=False).encode()).hexdigest()
        != manifest["canonical_sha256"]
    ):
        raise ValueError("Expression categories use a different canonical gene reference")
    assert (
        sha256(categories / "gene-categories.csv.gz")
        == manifest["outputs"]["gene-categories.csv.gz"]
    )
    frame = pd.read_csv(categories / "gene-categories.csv.gz", keep_default_na=False)
    safety = pd.read_csv(ROOT / "oncoref/data/cta-external-safety-reviews.csv")
    genes = partition_genes(
        frame, normal_evidence_genes=safety.loc[safety.action.eq("hold"), "Ensembl_Gene_ID"]
    )
    loci, transcripts, header = read_gtf(gtf)
    assert header["genome-build"] == "GRCh38.p14"
    reference = ReferenceIdentity(
        taxon_id=9606,
        assembly_accession=header["genome-build-accession"],
        annotation_source="Ensembl",
        annotation_release="115",
        source_version="115",
        identifier_namespaces={"gene": "Ensembl", "transcript": "Ensembl", "protein": "Ensembl"},
        content_provenance={"annotation_gtf": sha256(gtf), "protein_fasta": sha256(proteins)},
    )
    annotation = Annotation(loci)
    gene_groups = genes.set_index("gene_id").group.to_dict()
    identities, occurrences, unresolved = {}, [], []
    for header, sequence in fasta_records(proteins):
        bits = header.split()
        meta = dict(bit.split(":", 1) for bit in bits[1:] if ":" in bit)
        gid = meta["gene"].split(".")[0]
        tid = meta["transcript"].split(".")[0]
        if gid not in identities:
            identities[gid] = resolve_gene_identity(gid, genome=annotation).as_dict()
        mapping = identities[gid]
        canonical = (
            mapping["canonical_gene_id"] if mapping["status"] in {"canonical", "alias"} else None
        )
        record = {
            "reference_id": reference.reference_id,
            "occurrence_id": bits[0] + ":" + meta["transcript"],
            "gene_id": canonical or gid,
            "transcript_id": meta["transcript"],
            "protein_id": bits[0],
            "source_gene_id": meta["gene"],
            "source_header": header,
            "source": "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz",
            "sequence": sequence,
            "complete": False,
            "gene_mapping": mapping["status"],
            "expression_group": gene_groups.get(canonical or gid, "Unassessed"),
            "exclude_by_default": gene_groups.get(canonical or gid, "Unassessed")
            in {"Normal", "Unassessed"},
        }
        tx = transcripts.get(tid, {})
        reasons = []
        if tx.get("start") != 3 or tx.get("stop") != 3 or tx.get("incomplete", True):
            reasons.append("Incomplete or unverified translation")
        try:
            normalize_protein_sequence(sequence)
        except ValueError:
            reasons.append("Nonstandard or ambiguous residues")
        if reasons:
            unresolved.append({**record, "reason": "; ".join(reasons)})
        else:
            record["complete"] = True
            occurrences.append(record)
    groups = partition_proteins(
        occurrences, genes, reference=reference, expression_reference_id=reference.reference_id
    )
    rows, members = [], []
    for group in groups:
        rows.append(
            {
                "Protein sequence": group["protein_sequence_id"],
                "Group": group["group"],
                "Genes": ";".join(group["member_genes"]),
                "Length": len(group["sequence"]),
                "Source records": len(group["occurrences"]),
                "Mixed gene groups": group["mixed_gene_groups"],
                "Observed normal expression": group["observed_normal_expression"],
                "Exclude by default": group["exclude_by_default"],
                "Reference": reference.reference_id,
            }
        )
        for record in group["occurrences"]:
            members.append(
                {k: v for k, v in record.items() if k != "sequence"}
                | {"protein_sequence_id": group["protein_sequence_id"], "group": group["group"]}
            )
    table = pd.DataFrame(rows)
    complete_counts = defaultdict(set)
    for group in groups:
        for gene in group["member_genes"]:
            complete_counts[gene].add(group["protein_sequence_id"])
    genes["complete_protein_sequences"] = genes.gene_id.map(lambda gene: len(complete_counts[gene]))
    readable = genes.rename(
        columns={
            "gene_id": "Gene ID",
            "symbol": "Gene",
            "group": "Group",
            "expression_group": "Expression group",
            "other_subgroup": "Other subgroup",
            "observed_normal_expression": "Observed normal expression",
            "reason": "Reason",
            "exclude_by_default": "Exclude by default",
            "complete_protein_sequences": "Complete protein sequences",
        }
    )
    write_csv(readable, out / "genes.csv.gz")
    write_csv(table, out / "proteins.csv.gz")
    write_csv(pd.DataFrame(members), out / "protein-occurrences.csv.gz")
    write_csv(
        pd.DataFrame([{k: v for k, v in row.items() if k != "sequence"} for row in unresolved]),
        out / "translations-for-review.csv.gz",
    )
    for group in GROUPS:
        name = group.lower()
        write_csv(readable.loc[readable.Group.eq(group)], out / f"{name}-genes.csv.gz")
        write_csv(table.loc[table.Group.eq(group)], out / f"{name}-proteins.csv.gz")
        write_fasta(
            out / f"{name}.fasta.gz",
            ((p["protein_sequence_id"], p["sequence"]) for p in groups if p["group"] == group),
        )
    write_fasta(
        out / "translations-for-review.fasta.gz",
        ((row["occurrence_id"], row["sequence"]) for row in unresolved),
    )
    write_fasta(
        out / "additional-background.fasta.gz",
        (
            (row["occurrence_id"], row["sequence"])
            for row in unresolved
            if row["exclude_by_default"]
        ),
    )
    observed = {row.gene_id for row in genes.itertuples() if row.observed_normal_expression}
    write_fasta(
        out / "observed-normal.fasta.gz",
        (
            (p["protein_sequence_id"], p["sequence"])
            for p in groups
            if any(gene in observed for gene in p["member_genes"])
        ),
    )
    write_csv(pd.DataFrame(identities.values()), out / "gene-mappings.csv.gz")
    counts = pd.DataFrame(
        {
            "Group": GROUPS,
            "Genes": [int(genes.group.eq(g).sum()) for g in GROUPS],
            "Protein sequences": [int(table.Group.eq(g).sum()) for g in GROUPS],
        }
    )
    write_csv(counts, out / "counts.csv")
    represented = {row["gene_id"] for row in occurrences}
    missing = readable.loc[~readable["Gene ID"].isin(represented)]
    write_csv(missing, out / "genes-without-complete-proteins.csv.gz")
    (out / "reference.json").write_text(json.dumps(reference.as_dict(), indent=2) + "\n")
    (out / "exclusion-basis.json").write_text(
        json.dumps(
            {
                "policy": SET_POLICY_ID,
                "reference_id": reference.reference_id,
                "target": ["reproductive.fasta.gz"],
                "exclude": [
                    "normal.fasta.gz",
                    "unassessed.fasta.gz",
                    "additional-background.fasta.gz",
                ],
                "optional_stricter_background": ["observed-normal.fasta.gz"],
                "other": ["other.fasta.gz"],
                "rule": "Enumerate contiguous AA20 k-mers within each record; never cross stops, ambiguous residues or record boundaries. Subtract the exclusion union from target k-mers.",
                "caveat": "Disjoint full proteins can share k-mers. Gene RNA does not establish isoform-specific expression. Missing translations leave a coverage gap.",
            },
            indent=2,
        )
        + "\n"
    )
    lines = "\n".join(
        f"| {r.Group} | {r.Genes:,} | {r._2:,} |" for r in counts.itertuples(index=False)
    )
    (out / "index.md").write_text(
        f"""# Gene and protein sets for k-mer comparisons

Each gene and each complete exact protein sequence appears in **one group only**.

| Group | Genes | Exact protein sequences |
|---|---:|---:|
{lines}

## Choose a download

| Group | Gene list | Protein list | Sequences |
|---|---|---|---|
"""
        + "\n".join(
            f"| {g} | [Genes]({g.lower()}-genes.csv.gz) | [Proteins]({g.lower()}-proteins.csv.gz) | [FASTA]({g.lower()}.fasta.gz) |"
            for g in GROUPS
        )
        + f"""

## How shared proteins are assigned

All complete annotated isoforms are compared by exact amino-acid sequence. If any encoding gene is **Normal**, the shared protein is Normal. Otherwise, any **Unassessed** gene makes it Unassessed; then Reproductive takes precedence over Other. Every original gene, transcript and protein record remains in the [occurrence table](protein-occurrences.csv.gz). Verified alternate-contig identities use the existing annotation-aware resolver; matching names alone do not merge loci.

Gene groups follow the approved expression screen. External gene-specific normal-tissue evidence adds MAGEA12 to Normal; its original HPA expression group remains visible. These are sequence-defined protein forms, without PTM-specific distinctions. All isoforms inherit gene-level evidence; isoform-specific expression has not been established.

## Build a k-mer exclusion background

1. Enumerate target k-mers from **Reproductive**.
2. Subtract every k-mer in **Normal**, **Unassessed**, and the [additional background FASTA](additional-background.fasta.gz) containing incomplete or unusual translations from Normal/Unassessed genes.
3. Keep record boundaries intact. For incomplete or unusual translations, enumerate only contiguous standard amino-acid windows; never join across an ambiguous residue or a stop.
4. Other is available separately. It contains Disposable/Silent below-threshold expression, and does not supply positive reproductive evidence.

Full-protein sets are disjoint; their k-mers can still overlap. The subtraction step remains necessary. The [machine-readable recipe](exclusion-basis.json) lists the exact files. [Complete proteins with any observed normal expression](observed-normal.fasta.gz) are an optional stricter background; this can overlap Reproductive because a passing RNA cutoff does not erase a measured trace signal. For incomplete translations, the review table preserves each gene's group; the optional complete-protein file does not establish exhaustive sequence coverage.

## Coverage and provenance

Pinned human Ensembl **115 / GRCh38.p14**, all peptide FASTA records, with transcript completeness checked against start/stop codons and CDS completeness flags in the GTF. The gene categories use HPA v23. No expression cutoffs are transferred to another species.

**{len(unresolved):,} translation records** are incomplete or contain nonstandard/ambiguous residues. They remain [listed with reasons and gene groups](translations-for-review.csv.gz). Only records from Normal/Unassessed genes enter the additional background; incompleteness alone does not turn a Reproductive gene into a normal-tissue exclusion. **{len(missing):,} classified genes** have [no complete usable protein](genes-without-complete-proteins.csv.gz); that is a coverage gap, not evidence of sequence absence.

- [All gene assignments](genes.csv.gz) · [All protein assignments](proteins.csv.gz) · [Counts](counts.csv)
- [Gene mapping decisions](gene-mappings.csv.gz) · [Reference identity](reference.json) · [Checksums](run-manifest.json)
- [CTA curation report](../cta-curation-20261021-r2/index.md) · [Expression groups](../normal-expression-categories-20261021-r2/index.md)

Rebuild with `python scripts/normal_sequence_sets_report.py --out <empty-directory> --categories <expression-report> --gtf <Ensembl-115-GTF.gz> --proteins <Ensembl-115-peptide-FASTA.gz>`.
"""
    )
    index_path = out / "index.md"
    index_path.write_text(
        re.sub(
            r"\]\(([^()]+)\)",
            lambda match: (
                "](" + ASSET_ROOT + match[1] + ")" if (out / match[1]).is_file() else match[0]
            ),
            index_path.read_text(),
        )
    )
    sources = {
        "gtf": {
            "url": "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.chr_patch_hapl_scaff.gtf.gz",
            "sha256": sha256(gtf),
        },
        "proteins": {
            "url": "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz",
            "sha256": sha256(proteins),
        },
    }
    result = {
        "schema": SET_POLICY_ID,
        "asset_base_url": ASSET_ROOT,
        "reference": reference.as_dict(),
        "sources": sources,
        "inputs": {
            str(p.resolve().relative_to(ROOT)): sha256(p)
            for p in [
                categories / "gene-categories.csv.gz",
                categories / "run-manifest.json",
                ROOT / "oncoref/data/cta-external-safety-reviews.csv",
            ]
        },
        "implementation_sha256": {
            name: sha256(ROOT / name)
            for name in [
                "scripts/normal_sequence_sets_report.py",
                "oncoref/normal_sequence_sets.py",
                "oncoref/sequence_groups.py",
                "oncoref/gene_identity.py",
                "oncoref/reference_identity.py",
            ]
        },
        "outputs": {p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()},
    }
    (out / "run-manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    print(counts.to_string(index=False), flush=True)
    print(
        "Translations requiring review:",
        len(unresolved),
        "Genes without complete protein:",
        len(missing),
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["out", "categories", "gtf", "proteins"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    build(args.out, args.categories, args.gtf, args.proteins)

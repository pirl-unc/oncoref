#!/usr/bin/env python3
"""Add version-pinned protein-length metadata to the one-off CTA reports."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ID = "Ensembl_Gene_ID"
SOURCE_URL = "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
SOURCE_SHA256 = "2868236131342d14afbfa8dbd8a5e3fb56ba85210fdca196b65db3e17c7ece83"
LENGTH_METHOD = (
    "Protein length is the amino-acid count of the longest annotated Ensembl 112 "
    "GRCh38 protein sequence for each gene, matching the repository's proteoform "
    "convention. This is a longest-isoform choice, not an Ensembl canonical-transcript "
    "designation or a tumor-specific isoform measurement. Terminal stop markers are "
    "excluded; signal peptides and propeptides are included. Equal-length choices are "
    "resolved by protein ID, then transcript ID. Protein and transcript IDs, all "
    "annotated isoform lengths, and source hashes accompany the exports."
)


def build_annotations(universe, fasta):
    observed = hashlib.sha256(Path(fasta).read_bytes()).hexdigest()
    if observed != SOURCE_SHA256:
        raise ValueError(f"Expected pinned Ensembl 112 FASTA {SOURCE_SHA256}, got {observed}")
    genes = set(universe[ID])
    records = []
    header, sequence = None, []

    def record():
        if header is None:
            return
        fields = dict(item.split(":", 1) for item in header.split() if ":" in item)
        gene_id = fields.get("gene", "").split(".")[0]
        if gene_id not in genes:
            return
        seq = "".join(sequence).rstrip("*")
        if not seq or "*" in seq:
            raise ValueError(f"Empty sequence or internal stop: {header}")
        records.append(
            {
                ID: gene_id,
                "protein_id": header.split()[0],
                "protein_transcript_id": fields["transcript"],
                "protein_length_aa": len(seq),
                "protein_sequence_sha256": hashlib.sha256(seq.encode()).hexdigest(),
            }
        )

    with gzip.open(fasta, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                record()
                header, sequence = line[1:].strip(), []
            else:
                sequence.append(line.strip())
        record()
    isoforms = pd.DataFrame(records)
    assert not isoforms.duplicated([ID, "protein_id"]).any()
    assert set(isoforms[ID]) == genes, "Some CTA genes lack protein sequences"
    selected = isoforms.sort_values(
        ["protein_length_aa", "protein_id", "protein_transcript_id"], ascending=[False, True, True]
    ).drop_duplicates(ID)
    bounds = isoforms.groupby(ID).protein_length_aa.agg(
        protein_length_min_aa="min", protein_length_max_aa="max", protein_isoform_count="size"
    )
    annotations = universe[[ID, "Symbol"]].merge(selected, on=ID, validate="one_to_one")
    annotations = annotations.merge(bounds, on=ID, validate="one_to_one")
    annotations["protein_length_definition"] = "longest_annotated_protein"
    annotations["protein_annotation_release"] = 112
    annotations["protein_annotation_source"] = SOURCE_URL
    assert annotations.protein_length_aa.eq(annotations.protein_length_max_aa).all()
    isoforms = isoforms.merge(universe[[ID, "Symbol"]], on=ID, validate="many_to_one")
    return annotations.sort_values("Symbol"), isoforms.sort_values(["Symbol", "protein_id"])


def annotate_report(out, annotations):
    """Enrich presentation exports; preserve the raw analytical metric snapshot."""
    paths = {out / "cta_input_universe.csv"}
    paths.update((out / "gene_sets").glob("*/genes.csv"))
    paths.update((out / "gene_sets").glob("*/qualifying_gene_cohorts.csv"))
    paths.update((out / "selections").glob("*/*/genes_ranked.csv"))
    paths.update((out / "selections").glob("*/*/passing_gene_cohorts.csv"))
    for name in [
        "gene_set_membership.csv",
        "strictest_genes_ranked.csv",
        "strictest_passing_cancer_types.csv",
    ]:
        if (out / name).exists():
            paths.add(out / name)
    metadata = annotations.drop(columns="Symbol")
    for path in sorted(paths):
        if not path.exists():
            continue
        frame = pd.read_csv(path)
        original = frame.drop(columns=[c for c in frame if c.startswith("protein_")])
        enriched = original.merge(metadata, on=ID, how="left", sort=False, validate="many_to_one")
        assert enriched.protein_length_aa.notna().all(), str(path)
        assert len(enriched) == len(original)
        pd.testing.assert_frame_equal(original, enriched[original.columns])
        if not frame.equals(enriched):
            enriched.to_csv(path, index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / "outputs/cta_threshold_report_20260916")
    parser.add_argument("--out", type=Path, default=ROOT / "outputs/cta_threshold_report_20260917")
    parser.add_argument(
        "--fasta",
        type=Path,
        default=Path.home()
        / "Library/Caches/pyensembl/GRCh38/ensembl112/Homo_sapiens.GRCh38.pep.all.fa.gz",
    )
    args = parser.parse_args()
    universe = pd.read_csv(args.base / "cta_input_universe.csv")
    annotations, isoforms = build_annotations(universe, args.fasta)
    shipped = pd.read_csv(ROOT / "oncoref/data/proteoform-groups.csv")
    crosscheck = annotations.merge(shipped, left_on=ID, right_on="member_gene_id")
    assert crosscheck.protein_length_aa.eq(crosscheck.protein_length).all()
    manifest = {
        "source_url": SOURCE_URL,
        "source_local_path": str(args.fasta.resolve()),
        "source_sha256": hashlib.sha256(args.fasta.read_bytes()).hexdigest(),
        "ensembl_release": 112,
        "assembly": "GRCh38",
        "method": LENGTH_METHOD,
        "n_genes": len(annotations),
        "n_annotated_proteins": len(isoforms),
        "n_genes_missing_lengths": 0,
        "n_shipped_proteoform_lengths_verified": len(crosscheck),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    for out in [args.base, args.out]:
        annotations.to_csv(out / "cta_protein_lengths.csv", index=False)
        isoforms.to_csv(out / "cta_protein_isoforms.csv", index=False)
        (out / "protein_length_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        annotate_report(out, annotations)
    print(json.dumps({k: v for k, v in manifest.items() if k.startswith("n_")}, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Independently verify packaged policy and exact sequence release assets."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import pandas as pd


def records(path):
    name, chunks = None, []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip().split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def verify(assets, package):
    manifest = json.loads((package / "cta-policy.json").read_text())
    for name, expected in manifest["data_sha256"].items():
        assert hashlib.sha256((package / name).read_bytes()).hexdigest() == expected, name
    for name, expected in manifest["sequence_artifacts"]["sha256"].items():
        assert hashlib.sha256((assets / name).read_bytes()).hexdigest() == expected, name
    genes = pd.concat(
        [
            pd.read_csv(package / "target-gene-sets.csv.gz", low_memory=False),
            pd.read_csv(package / "sequence-gene-rna-reviews.csv.gz", low_memory=False),
        ]
    )
    assert genes["Gene ID"].is_unique
    categories = genes.set_index("Gene ID").Category.to_dict()
    table = pd.read_csv(assets / "target-proteins.csv.gz")
    assert table["Protein sequence"].is_unique
    priority = {"CTA": 0, "Neither": 1, "Unsafe": 2}
    expected = table.Genes.map(
        lambda value: max(
            (categories.get(g, "Neither") for g in value.split(";")), key=priority.get
        )
    )
    assert expected.eq(table.Group).all()
    seen = set()
    for group in ["CTA", "Neither", "Unsafe"]:
        ids = set()
        for pid, sequence in records(assets / f"{group.lower()}.fasta.gz"):
            assert (
                pid
                == "openvax:protein-sequence:v1:sha256:"
                + hashlib.sha256(sequence.encode()).hexdigest()
            )
            assert set(sequence) <= set("ACDEFGHIKLMNPQRSTVWY")
            assert pid not in seen
            seen.add(pid)
            ids.add(pid)
        assert ids == set(table.loc[table.Group.eq(group), "Protein sequence"])
        assert len(ids) == manifest["counts"]["proteins"][group]
    assert len(seen) == 129985
    members = pd.read_csv(assets / "target-protein-occurrences.csv.gz")
    assert members.group.eq(
        members.protein_sequence_id.map(table.set_index("Protein sequence").Group)
    ).all()
    assert members.gene_category.eq(members.gene_id.map(categories).fillna("Neither")).all()
    translations = pd.read_csv(assets / "translations-for-review.csv.gz")
    assert translations.forbidden.eq(translations.gene_id.map(categories).eq("Unsafe")).all()
    extra = {name for name, _ in records(assets / "additional-forbidden.fasta.gz")}
    assert extra == set(translations.loc[translations.forbidden, "occurrence_id"])
    for name in ["XAGE1A", "XAGE1B", "CTAG1A", "CTAG1B"]:
        row = genes.loc[genes.Gene.eq(name)].iloc[0]
        assert row.Category == "CTA"
        assert members.loc[members.gene_id.eq(row["Gene ID"]), "group"].eq("CTA").all()
    recipe = json.loads((assets / "exclusion-basis.json").read_text())
    assert recipe["policy"] == manifest["policy"]
    assert recipe["reference_id"] == manifest["reference"]
    assert recipe["forbidden"] == ["unsafe.fasta.gz", "additional-forbidden.fasta.gz"]
    print(
        json.dumps(
            {
                "policy": manifest["policy"],
                "verified_complete_proteins": len(seen),
                "counts": manifest["counts"]["proteins"],
                "additional_forbidden_records": len(extra),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assets", type=Path, required=True)
    parser.add_argument(
        "--package-data", type=Path, default=Path(__file__).resolve().parents[1] / "oncoref/data"
    )
    args = parser.parse_args()
    verify(args.assets, args.package_data)

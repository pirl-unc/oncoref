#!/usr/bin/env python3
"""Verify curated peptide origins against the exact pinned reference FASTA.

Usage: python scripts/verify_tcr_trial_sequences.py /path/to/Ensembl115.pep.all.fa.gz
This validates listed links; it is not an exhaustive cross-reactivity search.
"""

import argparse
import gzip
import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import tcr_trials  # noqa: E402


def verify(path: Path) -> int:
    """Fail if a peptide, isoform, reference checksum or audited position changed."""
    evidence = tcr_trials.evidence()
    if hashlib.sha256(path.read_bytes()).hexdigest() != evidence["sequence_audit"]["sha256"]:
        raise ValueError("Reference FASTA does not match the curated sequence audit")
    expected = {}
    for product in evidence["products"]:
        for target in product["targets"]:
            for match in target["sequence_matches"]:
                key = (target["gene_id"], match["protein_id"], target["peptide"])
                if key in expected and expected[key] != match:
                    raise ValueError("Contradictory peptide origin")
                expected[key] = match
    seen = set()
    header, chunks = None, []

    def check():
        if header is None:
            return
        gene = re.search(r"gene:(ENSG\d+)", header)
        if gene is None:
            return
        protein = header.split()[0][1:]
        sequence = "".join(chunks)
        for key, match in expected.items():
            if key[:2] != (gene[1], protein):
                continue
            start = match["start_1based"] - 1
            if (
                len(sequence) != match["length"]
                or hashlib.sha256(sequence.encode()).hexdigest() != match["sequence_sha256"]
                or sequence[start : start + len(key[2])] != key[2]
            ):
                raise ValueError(f"Sequence evidence changed: {key}")
            seen.add(key)

    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                check()
                header, chunks = line.strip(), []
            else:
                chunks.append(line.strip())
        check()
    if seen != expected.keys():
        raise ValueError(f"Missing peptide origins: {expected.keys() - seen}")
    return len(seen)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path)
    args = parser.parse_args()
    print(f"Verified {verify(args.fasta)} peptide/isoform links")

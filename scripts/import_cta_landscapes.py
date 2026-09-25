#!/usr/bin/env python3
"""Import complete published CTA landscapes before HPA/default filtering.

Use the original files and directory names documented in cta-landscape-sources.md.
Requires the reports extra and xlrd for the original XLS supplements.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from oncoref.cta_landscape import INPUTS, PAPERS, extract_landscapes
from oncoref.cta_sources import (
    add_publication_candidates,
    publication_membership,
    publication_sources,
)
from oncoref.load_dataset import get_data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument(
        "--out-dir", type=Path, default=Path(__file__).resolve().parents[1] / "oncoref/data"
    )
    parser.add_argument("--memberships-only", action="store_true")
    args = parser.parse_args()
    refs, sources = extract_landscapes(args.input_dir)
    old_refs, old_sources = publication_membership(), publication_sources()
    refs = pd.concat([old_refs[~old_refs.source_tag.isin(PAPERS)], refs], ignore_index=True)
    sources = pd.concat(
        [old_sources[~old_sources.source_tag.isin(PAPERS)], sources], ignore_index=True
    )
    print(
        refs.groupby("source_tag")
        .agg(
            source_entries=("source_symbol", "size"),
            mapped=("mapping_status", lambda values: values.eq("mapped").sum()),
            coding=("protein_candidate_eligible", "sum"),
        )
        .to_string(),
        flush=True,
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if not args.memberships_only:
        original = get_data("cancer-testis-antigens")
        candidates = add_publication_candidates(original, refs)
        candidates.to_csv(args.out_dir / "cancer-testis-antigens.csv", index=False)
        print(f"Candidate universe: {len(original)} -> {len(candidates)} genes", flush=True)
    refs.to_csv(args.out_dir / "cta-publication-membership.csv", index=False)
    sources.to_csv(args.out_dir / "cta-publication-sources.csv", index=False)
    receipt = {
        "inputs": INPUTS,
        "outputs": {
            name: hashlib.sha256((args.out_dir / name).read_bytes()).hexdigest()
            for name in ["cta-publication-membership.csv", "cta-publication-sources.csv"]
        },
    }
    (args.out_dir / "cta-landscape-import-receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()

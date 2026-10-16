#!/usr/bin/env python3
"""Audit every bundled identifier for the two driver-evidence tables, without downloads."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from oncoref import drivers
from oncoref.gene_ids import canonical_gene_space, ensembl_id_aliases, entrez_gene_mappings
from oncoref.load_dataset import get_data

INPUTS = (
    "canonical-gene-space.csv.gz",
    "ncbi-symbol-synonyms.csv.gz",
    "ensembl-id-aliases.csv",
    "ncbi-entrez-gene-mappings.csv.gz",
    "driver-gene-evidence.csv",
    "driver-variant-evidence.csv",
)


def audit_table(kind: str) -> dict:
    frame = get_data(f"driver-{kind}-evidence")
    query = getattr(drivers, f"driver_{kind}_evidence")
    space = canonical_gene_space()
    symbols = space.set_index("ensembl_gene_id")["symbol"].to_dict()
    ids_by_symbol = space.groupby("symbol")["ensembl_gene_id"].agg(set).to_dict()
    expected = frame.groupby("ensembl_gene_id")["record_id"].agg(set).to_dict()
    probes: dict[str, set[tuple[str, str]]] = {
        name: set()
        for name in (
            "source_symbol",
            "current_symbol",
            "recorded_alias",
            "ensembl_id",
            "versioned_ensembl_id",
            "migrated_ensembl_id",
            "entrez_id",
            "ncbi_synonym",
        )
    }
    renamed = []
    for row in (
        frame[["gene_symbol", "aliases", "ensembl_gene_id"]]
        .drop_duplicates()
        .itertuples(index=False)
    ):
        gid = row.ensembl_gene_id
        current = symbols[gid]
        probes["source_symbol"].add((row.gene_symbol, gid))
        probes["current_symbol"].add((current, gid))
        probes["ensembl_id"].add((gid, gid))
        probes["versioned_ensembl_id"].add((gid + ".7", gid))
        aliases = (
            [] if pd.isna(row.aliases) else [a.strip() for a in row.aliases.split(";") if a.strip()]
        )
        probes["recorded_alias"].update((alias, gid) for alias in aliases)
        if row.gene_symbol != current:
            renamed.append(
                {
                    "source_symbol": row.gene_symbol,
                    "current_symbol": current,
                    "ensembl_gene_id": gid,
                    "aliases": aliases,
                }
            )
    probes["migrated_ensembl_id"].update(
        (alias, gid) for alias, gid in ensembl_id_aliases().items() if gid in expected
    )
    probes["entrez_id"].update(
        (str(row.entrez_id), row.canonical_ensembl_gene_id)
        for row in entrez_gene_mappings().itertuples(index=False)
        if row.canonical_ensembl_gene_id in expected
    )
    for row in get_data("ncbi-symbol-synonyms").itertuples(index=False):
        ids = ids_by_symbol.get(row.official_symbol, set())
        if ids.intersection(expected):
            if len(ids) != 1:
                raise ValueError(
                    f"Nonunique driver synonym target: {row.alias}/{row.official_symbol}"
                )
            probes["ncbi_synonym"].add((row.alias, next(iter(ids))))
    for category, pairs in probes.items():
        for token, gid in sorted(pairs):
            observed = set(query(token)["record_id"])
            if observed != expected[gid]:
                raise ValueError(
                    f"{kind}: {category} {token!r} does not return exactly the rows for {gid}"
                )
    return {
        "records": len(frame),
        "genes": len(expected),
        "verified_queries": {name: len(pairs) for name, pairs in probes.items()},
        "renamed_symbols": sorted(renamed, key=lambda row: row["source_symbol"]),
    }


def build_audit() -> dict:
    data = Path(drivers.__file__).parent / "data"
    return {
        "scope": "Every source/current symbol, recorded alias, Ensembl ID and shipped Ensembl/Entrez/NCBI synonym mapping for the driver evidence tables. Source rows remain unchanged.",
        "inputs_sha256": {
            name: hashlib.sha256((data / name).read_bytes()).hexdigest() for name in INPUTS
        },
        "tables": {f"driver-{kind}-evidence": audit_table(kind) for kind in ("gene", "variant")},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path("docs/audits/driver-identifier-mappings.json")
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(build_audit(), indent=2) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit(f"Stale driver identifier audit: {args.output}")
    else:
        args.output.write_text(rendered)
    print(f"Verified driver identifier mappings: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

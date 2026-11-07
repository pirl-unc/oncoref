#!/usr/bin/env python3
"""Retrieve individual literature leads for the 96-protein follow-up review.

Retrieval is separate from adjudication. Search results do not create nominations,
normal-tissue exclusions, or evidence that a protein is absent from the literature.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import re
from pathlib import Path

from audit_cta_literature import SearchClient

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "scripts/cta_evidence_assets"
OUT = ROOT / "docs/audits/cta-literature-20261107"


def review_population():
    report = json.loads(
        gzip.decompress(
            (ROOT / "docs/audits/cta-evidence-guide-20261106/data.json.gz").read_bytes()
        )
    )
    literature = json.loads((ASSETS / "literature-review.json").read_text())
    catalogs = json.loads((ASSETS / "catalog-review.json").read_text())
    previously_reviewed = {
        label
        for finding in literature["findings"]
        if "proteinatlas.org" not in finding["url"]
        for label in finding["proteins"]
    }
    previously_reviewed.update(
        match["protein"] for catalog in catalogs["catalogs"] for match in catalog["matches"]
    )
    rows = sorted(
        (
            row
            for row in report["rows"]
            if row["category"] == "CTA"
            and not row["papers"]
            and row["label"] not in previously_reviewed
        ),
        key=lambda row: row["label"],
    )
    if len(rows) != 96:
        raise ValueError("The pinned 1.31.0 follow-up population must contain 96 proteins")
    return rows


def queries(row, cutoff):
    symbols = {row["label"]}
    aliases = {
        name
        for identity in row["gene_identities"]
        for name in identity.get("synonyms", [])
        if len(name) >= 4 and re.fullmatch(r"[A-Za-z][A-Za-z0-9.\-]+", name)
    }
    # These common words/abbreviations do not identify the requested gene alone.
    aliases -= {"JUNO", "AMAC", "IIS", "SC13", "MSOG", "HMGPI"}
    terms = sorted(symbols | aliases)
    identity = " OR ".join(f'TITLE_ABS:"{term}"' for term in terms)
    fulltext = " OR ".join(f'"{term}"' for term in sorted(symbols | set(row["gene_ids"])))
    date = f"FIRST_PDATE:[1900-01-01 TO {cutoff}]"
    return {
        "gene_and_aliases": f"({identity}) AND SRC:MED AND {date}",
        "fulltext_symbol_or_locus": f"({fulltext}) AND OPEN_ACCESS:Y AND {date}",
    }, terms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=ROOT / "tmp/cta-literature-20261107/raw")
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--cutoff", default="2026-11-07")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    rows = review_population()
    client = SearchClient(args.cache, interval=0.4)
    seeds, jobs = [], []
    for row in rows:
        searches, terms = queries(row, args.cutoff)
        seeds.append(
            {
                k: row[k]
                for k in (
                    "entity_id",
                    "label",
                    "gene_ids",
                    "gene_identities",
                    "biology_categories",
                    "route",
                )
            }
            | {"search_terms": terms}
        )
        for kind, query in searches.items():
            jobs.append((row, kind, query))

    def run(job):
        row, kind, query = job
        receipt = {
            "protein": row["label"],
            "entity_id": row["entity_id"],
            "kind": kind,
            "query": query,
            "cutoff": args.cutoff,
        }
        try:
            record, path = client.fetch(query, 25 if kind == "gene_and_aliases" else 10)
            data = record["response"]
            papers = data["resultList"]["result"]
            receipt.update(
                status="retrieved",
                hits=data["hitCount"],
                retrieved=len(papers),
                truncated=data["hitCount"] > len(papers),
                paper_ids=[p.get("source", "MED") + ":" + p["id"] for p in papers],
                url=record["url"],
                retrieved_utc=record["retrieved_utc"],
                cache_file=path.name,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            )
        except Exception as exc:
            receipt.update(status="failed", error=str(exc))
        return receipt

    receipts = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for receipt in pool.map(run, jobs):
            receipts.append(receipt)
            if len(receipts) % 16 == 0:
                print(f"Retrieved {len(receipts)}/{len(jobs)} searches", flush=True)
            (args.out / "search-ledger.json").write_text(json.dumps(receipts, indent=2) + "\n")
    (args.out / "review-population.json").write_text(json.dumps(seeds, indent=2) + "\n")
    print(
        f"Finished: {len(rows)} proteins; {len(receipts)} searches; "
        f"{sum(r['status'] == 'failed' for r in receipts)} failures",
        flush=True,
    )


if __name__ == "__main__":
    main()

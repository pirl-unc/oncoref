#!/usr/bin/env python3
"""Expand selected literature queries with cursor receipts and alias disambiguation.

The 25 October audit stays immutable. This creates a new retrieval, not a claim
of paper review or a systematic search. Ambiguous original queries remain in the
old ledger as a sensitivity source rather than becoming negative evidence.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path

import pandas as pd
from audit_cta_literature import SearchClient

ROOT = Path(__file__).resolve().parents[1]
QUALIFIED_TERMS = {
    'TITLE_ABS:"IMPACT"': '(TITLE_ABS:"IMPACT gene" OR TITLE_ABS:"IMPACT protein" OR TITLE_ABS:"imprinted and ancient")',
    'TITLE_ABS:"SPATIAL"': '(TITLE_ABS:"SPATIAL gene" OR TITLE_ABS:"SPATIAL protein" OR TITLE_ABS:"stromal protein associated with thymocyte")',
    'TITLE_ABS:"FATE"': '(TITLE_ABS:"FATE gene" OR TITLE_ABS:"FATE protein" OR TITLE_ABS:"fetal and adult testis expressed")',
    'TITLE_ABS:"MRP8"': '(TITLE_ABS:"multidrug resistance protein 8" OR (TITLE_ABS:"MRP8" AND TITLE_ABS:"transporter"))',
}


def qualify(query):
    for term, replacement in QUALIFIED_TERMS.items():
        query = query.replace(term, replacement)
    return query


def retrieve(client, query, *, page_size=1000, max_pages=100):
    cursor, seen, papers, pages, total = "*", set(), {}, [], None
    for _ in range(max_pages):
        if cursor in seen:
            raise ValueError("Repeated cursor; cannot claim complete retrieval")
        seen.add(cursor)
        record, path = client.fetch(query, page_size, cursor=cursor, result_type="lite")
        data = record["response"]
        if total is not None and int(data["hitCount"]) != total:
            raise ValueError("Search index changed during pagination; retry as a new retrieval")
        total = int(data["hitCount"])
        pages.append(
            {
                "cache_file": path.name,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "url": record["url"],
                "retrieved_utc": record["retrieved_utc"],
            }
        )
        for paper in data["resultList"]["result"]:
            papers[(paper.get("source", "MED"), paper["id"])] = paper
        if len(papers) == total:
            return list(papers.values()), pages, total, True
        if len(papers) > total:
            raise ValueError("Retrieved identity count exceeds declared hit count")
        cursor = data.get("nextCursorMark")
        if not cursor or not data["resultList"]["result"]:
            raise ValueError("Incomplete cursor chain; cannot claim complete retrieval")
    return list(papers.values()), pages, total, False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ledger", type=Path, default=ROOT / "docs/audits/cta-literature-20261025/searches.csv"
    )
    parser.add_argument("--symbols", nargs="+", help="Omit to expand every truncated query")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--max-pages", type=int, default=100)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    args.cache.mkdir(parents=True, exist_ok=True)
    ledger = pd.read_csv(args.ledger, keep_default_na=False)
    selected = (
        ledger.loc[ledger.Gene.isin(args.symbols)] if args.symbols else ledger.loc[ledger.Truncated]
    )
    if args.symbols and set(args.symbols) - set(selected.Gene):
        raise ValueError("A requested symbol is outside the mapped search ledger")
    client = SearchClient(args.cache)

    def run(row):
        query = qualify(row["Query"])
        result = {
            "Gene": row["Gene"],
            "Gene ID": row["Gene ID"],
            "Search": row["Search"],
            "Original query": row["Query"],
            "Query": query,
            "Query changed": query != row["Query"],
            "Paper review complete": False,
        }
        try:
            papers, pages, total, complete = retrieve(client, query, max_pages=args.max_pages)
            result.update(
                {
                    "Status": "Retrieved" if complete else "Page cap reached",
                    "Hits": total,
                    "Retrieved": len(papers),
                    "Complete pagination": complete,
                    "Pages": json.dumps(pages),
                    "Error": "",
                }
            )
            leads = [
                {
                    "Gene": row["Gene"],
                    "Gene ID": row["Gene ID"],
                    "Search": row["Search"],
                    "PMID": p["id"],
                    "Title": p.get("title", ""),
                    "Year": p.get("pubYear", ""),
                    "URL": f"https://pubmed.ncbi.nlm.nih.gov/{p['id']}/",
                    "Review": "Retrieved lead; identity and findings unassessed",
                }
                for p in papers
            ]
        except Exception as exc:
            result.update({"Status": "Failed", "Complete pagination": False, "Error": str(exc)})
            leads = []
        return result, leads

    results, leads = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for i, (result, records) in enumerate(pool.map(run, selected.to_dict("records")), 1):
            results.append(result)
            leads.extend(records)
            print(
                f"{i}/{len(selected)} {result['Gene']} {result['Search']}: {result['Status']}",
                flush=True,
            )
            pd.DataFrame(results).to_csv(args.out / "searches.csv", index=False)
    pd.DataFrame(leads).drop_duplicates(["Gene ID", "Search", "PMID"]).to_csv(
        args.out / "papers.csv", index=False
    )
    summary = {
        "original_ledger_sha256": hashlib.sha256(args.ledger.read_bytes()).hexdigest(),
        "queries": len(results),
        "pagination_complete": sum(r["Complete pagination"] for r in results),
        "primary_paper_assessments": 0,
        "scope": "Selected queries; historical cutoff retained; disambiguated variants do not establish exhaustive literature coverage",
        "outputs": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in args.out.glob("*.csv")
        },
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()

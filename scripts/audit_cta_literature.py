#!/usr/bin/env python3
"""Search each mapped CTA nomination and preserve a resumable literature ledger.

Search hits are leads, not adjudicated findings. This does not change target sets.
Raw Europe PMC core responses stay in the configurable local cache; the public
ledger records queries, dates, hit counts, retrieved IDs and checksums. Published
normal-tissue PubMed references from the pinned UniProt annotations are retrieved
separately so relevance ranking cannot hide those existing citations.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import re
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
ENDPOINT = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EXTRA_ALIASES = {
    "CTAG1A": ["NY-ESO-1", "NYESO1"],
    "CTAG1B": ["NY-ESO-1", "NYESO1"],
    "CTAG2": ["LAGE-1", "LAGE1"],
    "CTCFL": ["BORIS"],
    "CT55": ["CXorf48"],
    "CT83": ["KK-LC-1", "KKLC1"],
    "MAGEC2": ["HCA587", "MAGE-C2"],
    "XAGE1A": ["XAGE-1", "XAGE1"],
    "XAGE1B": ["XAGE-1", "XAGE1"],
    "YBX2": ["MSY2"],
    "ABCC11": ["MRP8"],
    "TPTE2": ["TPIP"],
}


class SearchClient:
    def __init__(self, cache, interval=0.3):
        self.cache = cache
        self.interval = interval
        self.lock = threading.Lock()
        self.key_locks = {}
        self.last = 0.0

    def fetch(self, query, page_size):
        params = {"query": query, "format": "json", "resultType": "core", "pageSize": page_size}
        key = hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()
        with self.lock:
            key_lock = self.key_locks.setdefault(key, threading.Lock())
        # Shared-family citation queries can be identical. Serialize by cache key
        # so every referencing ledger row identifies the same immutable response.
        with key_lock:
            return self._fetch(params, key)

    def _fetch(self, params, key):
        dest = self.cache / f"{key}.json.gz"
        if dest.exists():
            with gzip.open(dest, "rt") as handle:
                return json.load(handle), dest
        for attempt in range(5):
            with self.lock:
                delay = self.interval - (time.monotonic() - self.last)
                if delay > 0:
                    time.sleep(delay)
                self.last = time.monotonic()
            try:
                response = requests.get(ENDPOINT, params=params, timeout=60)
                response.raise_for_status()
                data = response.json()
                if "hitCount" not in data or "resultList" not in data:
                    raise ValueError(f"Invalid search response: {str(data)[:200]}")
                record = {
                    "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                    "url": response.url,
                    "parameters": params,
                    "response": data,
                }
                staging = dest.with_suffix(".tmp")
                with gzip.open(staging, "wt") as handle:
                    json.dump(record, handle)
                staging.replace(dest)
                return record, dest
            except (requests.RequestException, ValueError):
                if attempt == 4:
                    raise
                time.sleep(2**attempt)


def gene_queries(row, cutoff):
    symbol = row["Gene"] or row["Gene ID"]
    names = [symbol, *EXTRA_ALIASES.get(symbol, [])]
    if re.fullmatch(r"(?:MAGE|SSX|GAGE)[A-Z]?\d+[A-Z]?", symbol):
        names.append(re.sub(r"^(MAGE|SSX|GAGE)", r"\1-", symbol))
    # Exact annotation-linked synonyms are supplied by prepare_aliases, never a paralog substitution.
    names += row.get("Search aliases", [])
    terms = " OR ".join(f'TITLE_ABS:"{name}"' for name in sorted(set(names)))
    date = f"FIRST_PDATE:[1900-01-01 TO {cutoff}]"
    common = f"({terms}) AND SRC:MED AND {date}"
    return {
        "normal": common
        + " AND (normal OR healthy OR tissue OR testis OR brain OR expression OR isoform)",
        "therapy": common
        + ' AND (immunotherapy OR vaccine OR "T cell" OR toxicity OR safety OR trial)',
    }


def prepare_aliases(genes):
    from oncoref import protein_annotations

    wanted = set(genes["UniProt entries"].str.split("|").explode()) - {""}
    aliases = {}
    for entry in protein_annotations.reviewed_human():
        accession = entry["primaryAccession"]
        if accession not in wanted:
            continue
        names = []
        for gene in entry.get("genes", []):
            names.extend(
                x["value"]
                for x in gene.get("synonyms", [])
                if len(x["value"]) >= 4 and re.fullmatch(r"[A-Za-z][A-Za-z0-9-]+", x["value"])
            )
        aliases[accession] = sorted(set(names))
    return aliases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--atlas", type=Path, default=ROOT / "docs/audits/cta-atlas-20261024/all-genes.csv"
    )
    parser.add_argument("--out", type=Path, default=ROOT / "docs/audits/cta-literature-20261025")
    parser.add_argument("--cache", type=Path, default=ROOT / "tmp/cta-literature-20261025/raw")
    parser.add_argument("--cutoff", default="2026-10-25")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--page-size", type=int, default=20)
    parser.add_argument("--limit", type=int, help="Debug only; zero/omitted processes every gene")
    args = parser.parse_args()
    genes = pd.read_csv(args.atlas, keep_default_na=False)
    if len(genes) != 2540 or not genes["Gene ID"].is_unique:
        raise ValueError("Expected the unique 2,540-gene reviewed nomination universe")
    genes["sort"] = genes.Category.map({"CTA": 0, "Neither": 1, "Unsafe": 2})
    genes = genes.sort_values(["sort", "Gene"])
    if args.limit:
        genes = genes.head(args.limit)
    args.cache.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    aliases = prepare_aliases(genes)
    client = SearchClient(args.cache)
    jobs = []
    for row in genes.to_dict("records"):
        row["Search aliases"] = sorted(
            {name for acc in row["UniProt entries"].split("|") for name in aliases.get(acc, [])}
        )
        for kind, query in gene_queries(row, args.cutoff).items():
            jobs.append((row, kind, query, args.page_size))
        pmids = sorted(set(re.findall(r"\d+", row["Tissue PubMed IDs"])))
        if pmids:
            jobs.append(
                (
                    row,
                    "annotation citations",
                    "SRC:MED AND (" + " OR ".join("EXT_ID:" + p for p in pmids) + ")",
                    max(25, len(pmids)),
                )
            )

    def run(job):
        row, kind, query, size = job
        result = {
            "Gene": row["Gene"],
            "Gene ID": row["Gene ID"],
            "Category": row["Category"],
            "Search": kind,
            "Query": query,
            "Identity caveat": "Search terms may match a family or ambiguous alias; every paper needs identity review",
        }
        try:
            record, dest = client.fetch(query, size)
            data = record["response"]
            papers = data["resultList"]["result"]
            result.update(
                {
                    "Status": "Retrieved",
                    "Hits": data["hitCount"],
                    "Retrieved": len(papers),
                    "Truncated": data["hitCount"] > len(papers),
                    "Retrieved UTC": record["retrieved_utc"],
                    "Paper IDs": "|".join(p["id"] for p in papers),
                    "URL": record["url"],
                    "Cache file": dest.name,
                    "SHA256": hashlib.sha256(dest.read_bytes()).hexdigest(),
                    "Error": "",
                }
            )
        except Exception as exc:
            result.update({"Status": "Failed", "Error": str(exc)})
        return result

    ledger = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        for result in executor.map(run, jobs):
            ledger.append(result)
            if len(ledger) % 50 == 0:
                staging = args.out / "searches.csv.tmp"
                pd.DataFrame(ledger).to_csv(staging, index=False)
                staging.replace(args.out / "searches.csv")
                print(
                    f"{len(ledger)}/{len(jobs)} searches; failures={sum(r['Status'] == 'Failed' for r in ledger)}; latest={result['Gene']}",
                    flush=True,
                )
    pd.DataFrame(ledger).to_csv(args.out / "searches.csv", index=False)
    metadata = {
        "universe": len(genes),
        "searches": len(ledger),
        "cutoff": args.cutoff,
        "source_sha256": hashlib.sha256(args.atlas.read_bytes()).hexdigest(),
        "failed": sum(r["Status"] == "Failed" for r in ledger),
        "status": "Literature retrieval only; not a completed primary-paper adjudication",
    }
    (args.out / "search-summary.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata), flush=True)


if __name__ == "__main__":
    main()

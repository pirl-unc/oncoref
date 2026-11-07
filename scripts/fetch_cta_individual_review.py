#!/usr/bin/env python3
"""Fetch annotation-linked and selected papers for the individual CTA review.

Raw abstracts and full texts remain in the local cache. Public receipts retain
identifiers and checksums; manually written findings retain their evidence scope.
"""

from __future__ import annotations

import concurrent.futures
import gzip
import hashlib
import html
import json
import re
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import requests
from audit_cta_literature import SearchClient

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "tmp/cta-literature-20261107"
OUT = ROOT / "docs/audits/cta-literature-20261107"


def main():
    client = SearchClient(CACHE / "raw", interval=0.4)
    ledger = json.loads((OUT / "search-ledger.json").read_text())
    selected = json.loads(
        (ROOT / "scripts/cta_evidence_assets/review-retrieval-plan.json").read_text()
    )
    pmids = set(selected["pmids"])
    notes = json.loads(
        (ROOT / "scripts/cta_evidence_assets/individual-review-notes.json").read_text()
    )
    pmids.update(p for r in notes["reviews"] for p in r["pmids"])
    # Annotation-linked normal/function evidence is retrieved independently of ranking.
    pmids.update(selected["annotation_pmids"])
    receipts = []
    for offset in range(0, len(pmids), 50):
        chunk = sorted(pmids)[offset : offset + 50]
        query = "SRC:MED AND (" + " OR ".join("EXT_ID:" + p for p in chunk) + ")"
        record, path = client.fetch(query, 100)
        receipts.append(
            {
                "kind": "selected_and_annotation_citations",
                "query": query,
                "url": record["url"],
                "cache_file": path.name,
                "retrieved_utc": record["retrieved_utc"],
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "requested_ids": chunk,
                "retrieved_ids": [p["id"] for p in record["response"]["resultList"]["result"]],
            }
        )
    for protein, query in selected.get("queries", {}).items():
        record, path = client.fetch(query + " AND FIRST_PDATE:[1900-01-01 TO 2026-11-07]", 30)
        receipts.append(
            {
                "protein": protein,
                "kind": "refined_identity_query",
                "query": query,
                "url": record["url"],
                "cache_file": path.name,
                "retrieved_utc": record["retrieved_utc"],
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "hits": record["response"]["hitCount"],
            }
        )
    (OUT / "followup-retrievals.json").write_text(json.dumps(receipts, indent=2) + "\n")
    papers = {}
    for receipt in [*ledger, *receipts]:
        if not receipt.get("cache_file"):
            continue
        with gzip.open(CACHE / "raw" / receipt["cache_file"], "rt") as handle:
            record = json.load(handle)
        for paper in record["response"]["resultList"]["result"]:
            if paper.get("pmid"):
                papers[paper["pmid"]] = paper
    (CACHE / "papers.json").write_text(json.dumps(papers, indent=2))
    bibliography = {}
    for pmid in sorted({p for r in notes["reviews"] for p in r["pmids"]}):
        if pmid not in papers:
            raise ValueError(f"Missing cited paper {pmid}")
        paper = papers[pmid]
        bibliography[pmid] = {
            "pmid": pmid,
            "pmcid": paper.get("pmcid"),
            "doi": paper.get("doi"),
            "title": re.sub(r"<[^>]*>", "", html.unescape(paper["title"])),
            "authors": paper.get("authorString", ""),
            "year": paper.get("pubYear"),
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            "first_publication_date": paper.get("firstPublicationDate"),
        }
    (ROOT / "scripts/cta_evidence_assets/review-bibliography.json").write_text(
        json.dumps(bibliography, indent=2) + "\n"
    )
    wanted = set(selected.get("fulltext_pmids", []))
    wanted.update(p for r in notes["reviews"] for p in r["pmids"])
    (CACHE / "fulltext").mkdir(exist_ok=True)

    def fetch(pmid):
        paper = papers.get(pmid, {})
        pmcid = paper.get("pmcid")
        if not pmcid:
            return {"pmid": pmid, "status": "no_pmc_fulltext"}
        dest = CACHE / "fulltext" / (pmcid + ".xml")
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
        if not dest.exists():
            for attempt in range(3):
                response = requests.get(url, timeout=45)
                if response.status_code == 200:
                    try:
                        ET.fromstring(response.content)
                    except ET.ParseError:
                        return {"pmid": pmid, "status": "invalid_xml", "url": url}
                    dest.write_bytes(response.content)
                    break
                if response.status_code in {403, 404}:
                    return {
                        "pmid": pmid,
                        "pmcid": pmcid,
                        "status": "fulltext_unavailable",
                        "url": url,
                    }
                time.sleep(2**attempt)
            else:
                return {"pmid": pmid, "status": "retrieval_failed", "url": url}
        return {
            "pmid": pmid,
            "pmcid": pmcid,
            "status": "retrieved",
            "url": url,
            "cache_file": dest.name,
            "sha256": hashlib.sha256(dest.read_bytes()).hexdigest(),
        }

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        texts = list(pool.map(fetch, sorted(wanted)))
    (OUT / "fulltext-retrievals.json").write_text(json.dumps(texts, indent=2) + "\n")
    print(
        f"Indexed {len(papers)} records; {sum(r['status'] == 'retrieved' for r in texts)}/{len(texts)} full texts"
    )


if __name__ == "__main__":
    main()

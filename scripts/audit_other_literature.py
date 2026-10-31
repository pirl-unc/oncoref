#!/usr/bin/env python3
"""Search every baseline Other member and retain source-linked research leads.

Retrieval is not adjudication. Immutable raw responses have checksums and dates;
accepted findings live separately in the gene-knowledge reference.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import re
from pathlib import Path

import pandas as pd
from audit_cta_literature import SearchClient

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-other-literature-20261030"


def prepare(out: Path):
    from oncoref import protein_annotations

    out.mkdir(parents=True, exist_ok=True)
    baseline = out / "baseline-proteins.csv.gz"
    if not baseline.exists():
        frame = pd.read_csv(
            ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz",
            low_memory=False,
        )
        frame = frame[frame.category.eq("Other")]
        if len(frame) != 600:
            raise ValueError("Freeze the 1.25.0 baseline before changing the policy")
        frame.to_csv(baseline, index=False, compression={"method": "gzip", "mtime": 0})
    frame = pd.read_csv(baseline, low_memory=False)
    rows = []
    for protein in frame.itertuples():
        for member in json.loads(protein.member_identities):
            aliases = json.loads(member["gene_name_synonyms"])
            rows.append(
                {
                    "gene_id": member["canonical_gene_id"],
                    "gene_symbol": member["canonical_symbol"] or member["canonical_gene_id"],
                    "aliases": json.dumps(aliases),
                    "protein_id": protein.protein_equivalence_group_id,
                    "protein_label": protein.protein_group_label,
                    "baseline_classification": protein.classification,
                }
            )
    genes = pd.DataFrame(rows).sort_values("gene_id")
    if genes.gene_id.duplicated().any() or len(genes) != 611:
        raise ValueError("Expected 611 distinct member genes")
    genes.to_csv(out / "genes.csv", index=False)
    annotations = out / "uniprot-annotations.csv"
    if not annotations.exists():
        protein_annotations.gene_annotations(genes.gene_id).to_csv(annotations, index=False)
    return genes, pd.read_csv(annotations, keep_default_na=False)


def search(out: Path, cache: Path, cutoff: str, workers: int):
    genes, annotations = prepare(out)
    cache.mkdir(parents=True, exist_ok=True)
    client = SearchClient(cache)
    jobs = []
    for gene in genes.to_dict("records"):
        names = {gene["gene_symbol"], gene["gene_id"]}
        names.update(
            n
            for n in json.loads(gene["aliases"])
            if len(n) >= 4 and re.fullmatch(r"[A-Za-z][A-Za-z0-9-]+", n)
        )
        terms = " OR ".join(f'TITLE_ABS:"{n}"' for n in sorted(names))
        common = f"({terms}) AND SRC:MED AND FIRST_PDATE:[1900-01-01 TO {cutoff}]"
        for topic, qualifier in [
            (
                "normal",
                "normal OR healthy OR expression OR tissue OR neuron OR immune OR testis OR ovary OR skin OR hair OR placenta OR embryo",
            ),
            (
                "cancer",
                "cancer OR tumor OR tumour OR antigen OR malignancy OR carcinoma OR leukemia",
            ),
        ]:
            jobs.append((gene, topic, common + " AND (" + qualifier + ")"))
        local = annotations[annotations.gene_id.eq(gene["gene_id"])]
        pmids = set()
        for col in [c for c in local if c.endswith("_pubmed_ids")]:
            pmids.update(p for s in local[col] for p in s.split("|") if p)
        if pmids:
            jobs.append(
                (
                    gene,
                    "annotation references",
                    "SRC:MED AND (" + " OR ".join("EXT_ID:" + p for p in sorted(pmids)) + ")",
                )
            )

    def run(job):
        gene, topic, query = job
        cursor, pages, papers, seen = "*", [], {}, set()
        result = {**gene, "topic": topic, "query": query}
        try:
            for _ in range(50):
                if cursor in seen:
                    raise ValueError("Repeated search cursor")
                seen.add(cursor)
                record, path = client.fetch(query, 100, cursor=cursor)
                data = record["response"]
                pages.append(
                    {
                        "file": path.name,
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "url": record["url"],
                        "retrieved_utc": record["retrieved_utc"],
                    }
                )
                for p in data["resultList"]["result"]:
                    papers[p["id"]] = p
                total = int(data["hitCount"])
                if len(papers) >= total:
                    break
                cursor = data.get("nextCursorMark")
                if not cursor:
                    break
            result.update(
                status="retrieved",
                hits=total,
                retrieved=len(papers),
                complete=len(papers) == total,
                pages=json.dumps(pages),
                error="",
            )
        except Exception as exc:
            result.update(status="failed", complete=False, pages=json.dumps(pages), error=str(exc))
        leads = []
        for p in papers.values():
            # Full responses remain in the private local cache. The published
            # ledger contains bibliography, not redistributed full abstracts.
            leads.append(
                {
                    "gene_id": gene["gene_id"],
                    "gene_symbol": gene["gene_symbol"],
                    "topic": topic,
                    "pmid": p["id"],
                    "title": re.sub(r"<[^>]+>", "", p.get("title", "")),
                    "year": p.get("pubYear", ""),
                    "doi": p.get("doi", ""),
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{p['id']}/",
                    "review_status": "search_lead",
                }
            )
        return result, leads

    ledger, papers = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for i, (receipt, leads) in enumerate(pool.map(run, jobs), 1):
            ledger.append(receipt)
            papers.extend(leads)
            if i % 25 == 0 or i == len(jobs):
                pd.DataFrame(ledger).to_csv(out / "searches.csv", index=False)
                pd.DataFrame(papers).drop_duplicates(["gene_id", "topic", "pmid"]).to_csv(
                    out / "search-leads.csv.gz",
                    index=False,
                    compression={"method": "gzip", "mtime": 0},
                )
                print(
                    f"{i}/{len(jobs)} queries; {sum(r['status'] == 'failed' for r in ledger)} failed",
                    flush=True,
                )
    (out / "search-manifest.json").write_text(
        json.dumps(
            {
                "scope": "611 member genes of 600 baseline Other proteins; retrieval is not paper adjudication",
                "cutoff": cutoff,
                "queries": len(ledger),
                "successful": sum(r["status"] == "retrieved" for r in ledger),
                "complete_pagination": sum(r["complete"] for r in ledger),
                "files": {
                    p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in out.iterdir()
                    if p.name
                    in {
                        "genes.csv",
                        "baseline-proteins.csv.gz",
                        "uniprot-annotations.csv",
                        "searches.csv",
                        "search-leads.csv.gz",
                    }
                },
            },
            indent=2,
        )
        + "\n"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--cache", type=Path, default=ROOT / "tmp/cta-other-literature/raw")
    parser.add_argument("--cutoff", default="2026-10-30")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    search(args.out, args.cache, args.cutoff, args.workers)


if __name__ == "__main__":
    main()

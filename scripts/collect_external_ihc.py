#!/usr/bin/env python3
"""Inventory public MSVA tissue panels without treating supplier claims as adjudication.

Downloads are cached outside the repository. Extracted observations are deliberately
coarse: organ-level negative images do not establish absence in every rare cell type.
Raw HTML, receipt hashes, unresolved target mappings and absent tables remain auditable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from html.parser import HTMLParser
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref.gene_matching import match_reference_genes  # noqa: E402
from oncoref.ihc_observations import parse_staining_statement  # noqa: E402

BASE = "https://ms-validatedantibodies.com"
DATE = "2026-10-27"


class Text(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts = []

    def handle_data(self, data):
        self.parts.append(data)


def plain(html: str) -> str:
    parser = Text()
    parser.feed(html)
    return " ".join(" ".join(parser.parts).split())


def parse_panel(html: str, url: str) -> tuple[dict, list[dict]]:
    """Extract explicit table and paragraph observations without inventing tissue matches."""
    title = re.search(r"<h1[^>]*>(.*?)</h1>", html, re.S | re.I)
    title = plain(title.group(1)) if title else ""
    section = re.search(
        r"<h[23][^>]*>\s*Staining Pattern in Normal Tissues\s*</h[23]>(.*?)"
        r"<h[23][^>]*>\s*Staining Pattern in Relevant Tumor Types",
        html,
        re.S | re.I,
    )
    if not section:
        return {"title": title, "url": url, "parse_status": "no_normal_table"}, []
    target_links = sorted(set(re.findall(r"proteinatlas\.org/(ENSG\d{11})", section.group(1))))
    # Only the normal panel's own HPA links can identify the target. Related
    # products and comparison paragraphs are excluded from target matching.
    record = {
        "title": title,
        "url": url,
        "source_gene_id": target_links[0] if len(target_links) == 1 else "",
        "target_link_candidates": json.dumps(target_links),
        "parse_status": "parsed",
        "source_class": "manufacturer_panel",
    }
    rows = []
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", section.group(1), re.S | re.I):
        cells = [plain(x) for x in re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S | re.I)]
        if len(cells) < 2:
            continue
        tissue, description = cells[-2:]
        if not tissue or tissue.lower() == "remarks":
            continue
        parsed = parse_staining_statement(description)
        parsed["reported_percentages"] = json.dumps(parsed["reported_percentages"])
        rows.append(
            {
                "source_url": url,
                "tissue": tissue,
                "cell_type": "unspecified",
                **parsed,
                "exact_cell_coverage": False,
                "assay": "IHC",
                "source_class": "manufacturer_panel",
                "review_status": "pending",
                "reagent_title": title,
                "source_locator": "Staining Pattern in Normal Tissues / " + tissue,
                "anatomy_mapping_status": "organ_only",
                "extraction_method": "table_text",
            }
        )
    if not rows:
        normal_text = plain(section.group(1))
        # Control recommendations and image captions are not independent specimens.
        normal_text = re.split(
            r"Positive control|Appropriate control|Suggested positive", normal_text, flags=re.I
        )[0]
        statements = re.split(r"(?<=[.!])\s+(?=[A-Z])", normal_text)
        for number, statement in enumerate(statements, 1):
            parsed = parse_staining_statement(statement)
            if parsed["result"] == "not_assessed":
                continue
            parsed["reported_percentages"] = json.dumps(parsed["reported_percentages"])
            rows.append(
                {
                    "source_url": url,
                    "tissue": "unresolved",
                    "cell_type": "unspecified",
                    **parsed,
                    "exact_cell_coverage": False,
                    "assay": "IHC",
                    "source_class": "manufacturer_panel",
                    "review_status": "pending",
                    "reagent_title": title,
                    "source_locator": f"Normal tissue description / sentence {number}",
                    "anatomy_mapping_status": "requires_context_mapping",
                    "extraction_method": "paragraph_text",
                }
            )
    record["normal_contexts"] = len(rows)
    record["negative_contexts"] = sum(r["result"] == "not_detected" for r in rows)
    record["parse_status"] = "categorical_observations" if rows else "no_explicit_observation"
    return record, rows


def collect(output: Path, cache: Path) -> None:
    """Collect the catalog, identity audit, categorical observations and source receipts."""
    output.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    receipts = []

    def fetch(url):
        key = hashlib.sha256(url.encode()).hexdigest()[:20]
        path = cache / (key + ".html")
        if not path.exists():
            response = requests.get(url, timeout=40)
            response.raise_for_status()
            path.write_bytes(response.content)
        data = path.read_bytes()
        return data.decode("utf-8"), {
            "url": url,
            "sha256": hashlib.sha256(data).hexdigest(),
            "cache_path": str(path),
            "retrieved": DATE,
        }

    first, receipt = fetch(BASE + "/products/")
    receipts.append(receipt)
    pages = {BASE + "/products/", *re.findall(r'href="(' + BASE + r'/products/page/\d+/)"', first)}
    links = set()
    for page in sorted(pages):
        content, receipt = fetch(page)
        receipts.append(receipt)
        links.update(re.findall(r'href="(' + BASE + r'/product/[^"?#]+/)"', content))
    records, observations, failures = [], [], []
    with ThreadPoolExecutor(max_workers=4) as pool:

        def one(url):
            try:
                html, receipt = fetch(url)
                record, rows = parse_panel(html, url)
                record["content_sha256"] = receipt["sha256"]
                return record, rows, receipt, None
            except Exception as error:
                return None, [], None, {"url": url, "error": str(error)}

        for i, (record, rows, receipt, error) in enumerate(pool.map(one, sorted(links)), 1):
            if error:
                failures.append(error)
            else:
                records.append(record)
                observations.extend(rows)
                receipts.append(receipt)
            if i % 20 == 0 or i == len(links):
                print(f"Catalog pages {i}/{len(links)}; failed {len(failures)}", flush=True)
    inventory = pd.DataFrame(records)
    inventory["source_gene_id"] = inventory.get("source_gene_id", "").fillna("")
    matches = match_reference_genes(
        inventory, id_column="source_gene_id", source="MSVA", source_version="retrieved:" + DATE
    )
    inventory = pd.concat(
        [
            inventory.reset_index(drop=True),
            matches.drop(columns="source_gene_id").reset_index(drop=True),
        ],
        axis=1,
    )
    facts = pd.read_csv(
        ROOT / "docs/audits/cta-expanded-funnel-20261027/gene-classifications.csv.gz"
    )
    inventory = inventory.merge(
        facts[
            [
                "gene_id",
                "symbol",
                "category",
                "ihc_missing",
                "ihc_reliable_records",
                "Antibody concern",
                "quantitative_candidate",
            ]
        ],
        left_on="canonical_gene_id",
        right_on="gene_id",
        how="left",
        validate="many_to_one",
    )
    inventory.to_csv(output / "external-panel-inventory.csv", index=False)
    obs = (
        pd.DataFrame(observations)
        .merge(
            inventory[["url", "canonical_gene_id", "identity_status", "content_sha256"]],
            left_on="source_url",
            right_on="url",
            how="left",
            validate="many_to_one",
        )
        .drop(columns="url")
    )
    obs.to_csv(
        output / "external-panel-observations.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.DataFrame(failures, columns=["url", "error"]).to_csv(
        output / "collection-failures.csv", index=False
    )
    priority = facts[facts.quantitative_candidate & facts.category.eq("Neither")].copy()
    priority["external_panel_links"] = priority.gene_id.map(
        inventory.groupby("canonical_gene_id").url.agg(lambda s: " | ".join(sorted(s)))
    )
    priority.to_csv(output / "antibody-only-review-queue.csv", index=False)
    summary = {
        "discovered_products": len(links),
        "downloaded_products": len(records),
        "failed_products": len(failures),
        "parsed_normal_panels": int(inventory.parse_status.eq("categorical_observations").sum()),
        "unmapped_anatomy_observations": int(
            obs.anatomy_mapping_status.eq("requires_context_mapping").sum()
        ),
        "result_counts": obs.result.value_counts().to_dict(),
        "unique_matched_genes": inventory.loc[
            inventory.identity_status.eq("matched"), "canonical_gene_id"
        ].nunique(),
        "normal_observations": len(obs),
        "antibody_only_priority_genes": len(priority),
        "priority_genes_with_panel": int(priority.external_panel_links.notna().sum()),
        "automatic_overrides": 0,
    }
    (output / "collection-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "source-receipts.json").write_text(json.dumps(receipts, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "docs/audits/external-ihc-20261027")
    parser.add_argument("--cache", type=Path, default=Path("/tmp/oncoref-external-ihc-msva"))
    args = parser.parse_args()
    collect(args.output, args.cache)

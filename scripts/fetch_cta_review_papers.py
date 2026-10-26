#!/usr/bin/env python3
"""Retrieve public primary-paper metadata and available OA XML for targeted review."""

import argparse
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree

import requests

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("pmids", nargs="+")
parser.add_argument("--out", type=Path, default=Path("tmp/cta-literature-20261025/targeted"))
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=True)
ledger_path = args.out / "retrieval.json"
previous = json.loads(ledger_path.read_text()) if ledger_path.exists() else []
records = {row["PMID"]: row for row in previous}
base = "https://www.ebi.ac.uk/europepmc/webservices/rest/"
for pmid in args.pmids:
    if not pmid.isdecimal():
        parser.error("Supply PubMed numeric IDs")
    metadata = args.out / f"{pmid}.json"
    if metadata.exists():
        record = json.loads(metadata.read_text())
    else:
        response = requests.get(
            base + "search",
            params={"query": f"EXT_ID:{pmid} AND SRC:MED", "resultType": "core", "format": "json"},
            timeout=60,
        )
        response.raise_for_status()
        data = response.json()["resultList"]["result"]
        if len(data) != 1 or data[0]["id"] != pmid:
            raise ValueError(f"No unique PubMed match: {pmid}")
        record = {
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "url": response.url,
            "paper": data[0],
        }
        metadata.write_text(json.dumps(record, indent=2))
        time.sleep(0.3)
    pmc = record["paper"].get("pmcid")
    status = "Metadata/abstract only"
    if pmc:
        xmlpath = args.out / f"{pmc}.xml"
        if not xmlpath.exists():
            response = requests.get(base + pmc + "/fullTextXML", timeout=60)
            if response.status_code == 200:
                ElementTree.fromstring(response.content)
                xmlpath.write_bytes(response.content)
            time.sleep(0.3)
        if xmlpath.exists():
            tree = ElementTree.fromstring(xmlpath.read_bytes())
            blocks = []
            for element in tree.iter():
                if element.tag in {"article-title", "title", "p", "table-wrap"}:
                    blocks.append(" ".join("".join(element.itertext()).split()))
            (args.out / f"{pmid}.txt").write_text("\n\n".join(blocks))
            status = "Full text XML retrieved"
    records[pmid] = {
        "PMID": pmid,
        "PMCID": pmc,
        "title": record["paper"]["title"],
        "status": status,
        "metadata_sha256": hashlib.sha256(metadata.read_bytes()).hexdigest(),
    }
    print(pmid, status, record["paper"]["title"], flush=True)
ledger_path.write_text(json.dumps(list(records.values()), indent=2))

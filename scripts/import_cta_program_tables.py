#!/usr/bin/env python3
"""Extract complete, pinned cancer–embryo program tables, preserving source rows.

Download the three linked workbooks into --inputs, then run this importer. The
standard-library reader avoids an Excel dependency in the release environment.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "scripts/cta_evidence_assets"
FILES = {
    "Chew2019-mmc2.xlsx": (
        "https://ars.els-cdn.com/content/image/1-s2.0-S1534580719305295-mmc2.xlsx",
        "b8ae23bd27e6039abb2efb4a1d987fcceabe244d4e7d77e61d9f0230d7aedd6f",
    ),
    "Chew2019-mmc4.xlsx": (
        "https://ars.els-cdn.com/content/image/1-s2.0-S1534580719305295-mmc4.xlsx",
        "4a176f519ae81388881b0a0c459725635774684eda0ba3b28de91c3c440cfb65",
    ),
    "Smith2023-mmc4.xlsx": (
        "https://ars.els-cdn.com/content/image/1-s2.0-S2211124723011257-mmc4.xlsx",
        "70fbf0881a6ea959d28ed96f64c78bbd8133790b51c341a362ce86cde36a0394",
    ),
}
NS = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def cells(path, sheet):
    """Read the text/numeric cells of one worksheet, indexed by Excel address."""
    with zipfile.ZipFile(path) as archive:
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            strings = [
                "".join(s.itertext()) for s in ET.fromstring(archive.read("xl/sharedStrings.xml"))
            ]
        root = ET.fromstring(archive.read(f"xl/worksheets/sheet{sheet}.xml"))
        out = {}
        for cell in root.findall(".//x:c", NS):
            value = cell.find("x:v", NS)
            if cell.get("t") == "inlineStr":
                out[cell.get("r")] = "".join(cell.find("x:is", NS).itertext())
            elif value is not None and value.text is not None:
                out[cell.get("r")] = (
                    strings[int(value.text)] if cell.get("t") == "s" else value.text
                )
        return out


def extract(inputs):
    receipts = []
    for name, (url, digest) in FILES.items():
        path = inputs / name
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != digest:
            raise ValueError(f"Source hash mismatch: {name}")
        receipts.append({"file": name, "url": url, "sha256": digest})
    rows = []
    c = cells(inputs / "Chew2019-mmc4.xlsx", 1)
    for n in range(4, 51):
        if f"B{n}" in c:
            rows.append(
                {
                    "source_id": "10.1016/j.devcel.2019.06.011",
                    "table": "S3",
                    "sheet": "Table S3",
                    "row": n,
                    "program": "DUX4 high-confidence targets",
                    "type": c.get(f"A{n}", ""),
                    "symbol": c[f"B{n}"],
                    "gene_id": c.get(f"C{n}", ""),
                }
            )
    c = cells(inputs / "Smith2023-mmc4.xlsx", 1)
    for n in range(5, 722):
        for col, program in [("A", "ZGA genes"), ("B", "8C genes")]:
            if c.get(f"{col}{n}"):
                rows.append(
                    {
                        "source_id": "10.1016/j.celrep.2023.113114",
                        "table": "S3A",
                        "sheet": "A",
                        "row": n,
                        "column": col,
                        "program": program,
                        "type": "gene",
                        "symbol": c[f"{col}{n}"],
                        "gene_id": "",
                        "upstream_source_id": "10.1016/j.stem.2022.01.014",
                    }
                )
    # S3B annotates the 8C program with DUX4-target status; it is not an independent nomination.
    c = cells(inputs / "Smith2023-mmc4.xlsx", 2)
    for n in range(2, 720):
        if c.get(f"A{n}"):
            rows.append(
                {
                    "source_id": "10.1016/j.celrep.2023.113114",
                    "table": "S3B",
                    "sheet": "B",
                    "row": n,
                    "program": "8C program: DUX4-target annotation",
                    "type": "gene",
                    "symbol": c[f"A{n}"],
                    "gene_id": "",
                    "direct_target_marker": c.get(f"B{n}", ""),
                }
            )
    c = cells(inputs / "Chew2019-mmc2.xlsx", 1)
    screen = [
        {"row": n, "gene_id": c[f"A{n}"], "symbol": c.get(f"B{n}", "")}
        for n in range(4, 12701)
        if c.get(f"A{n}")
    ]
    return {
        "schema": "oncoref.program-table-extraction",
        "version": 1,
        "source_files": receipts,
        "rows": rows,
        "screen_only": {
            "source_id": "10.1016/j.devcel.2019.06.011",
            "table": "S1",
            "role": "background_screen",
            "description": "All scored genes; membership does not nominate a CTA.",
            "rows": screen,
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=ASSETS / "program-source-rows.json.gz")
    args = parser.parse_args()
    result = extract(args.inputs)
    args.out.write_bytes(gzip.compress((json.dumps(result, indent=2) + "\n").encode(), mtime=0))
    print(
        f"Extracted {len(result['rows'])} program rows and {len(result['screen_only']['rows'])} screening rows"
    )

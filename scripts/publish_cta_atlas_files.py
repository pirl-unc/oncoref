#!/usr/bin/env python3
"""Assemble local static report assets; this does not upload anything."""

import hashlib
import json
import shutil
import zipfile
from pathlib import Path

import pandas as pd
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-atlas-20261024"
genes = pd.read_csv(OUT / "all-genes.csv.gz", low_memory=False).fillna("")
index = pd.read_csv(OUT / "pdf-gene-index.csv")
assert len(index) == len(genes) == 2540 and not index["Gene ID"].duplicated().any()
genes = genes.merge(index, on="Gene ID", validate="one_to_one")
template = (ROOT / "scripts/templates/cta-atlas-explorer.html").read_text()
payload = json.dumps(genes.to_dict("records"), ensure_ascii=False, allow_nan=False).replace(
    "</", "<\\/"
)
(OUT / "explore.html").write_text(template.replace("__DATA__", payload))
(OUT / "pdf").mkdir(exist_ok=True)
pages = {}
for path in (ROOT / "output/pdf").glob("cta-*.pdf"):
    shutil.copy2(path, OUT / "pdf" / path.name)
    pages[path.name] = len(PdfReader(path).pages)
shutil.copy2(ROOT / "tmp/pdfs/cancer-rna-heatmap.png", OUT / "cancer-rna-heatmap.png")
genes.to_csv(OUT / "all-genes.csv", index=False)
with zipfile.ZipFile(
    OUT / "cta-evidence-tables.zip", "w", compression=zipfile.ZIP_DEFLATED
) as archive:
    for path in sorted(OUT.iterdir()):
        if path.name != "artifacts.json" and (
            path.suffix in {".csv", ".json", ".yml"} or path.name.endswith(".csv.gz")
        ):
            archive.write(path, path.name)
manifest = {
    "policy": "oncoref.target-sets.2026-10-24.1",
    "date": "2026-10-24",
    "pdf_pages": pages,
    "files": {
        str(p.relative_to(OUT)): {
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            "bytes": p.stat().st_size,
        }
        for p in sorted(OUT.rglob("*"))
        if p.is_file() and p.name != "artifacts.json"
    },
}
(OUT / "artifacts.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps(pages, indent=2))

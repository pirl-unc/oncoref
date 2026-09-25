#!/usr/bin/env python3
"""Render all CTA source/curation figures and a combined vector PDF."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path

from pypdf import PdfReader, PdfWriter

from oncoref import cta_curation_plots


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = cta_curation_plots.render(args.out)
    combined = args.out / "oncoref-cta-landscape-figures.pdf"
    writer = PdfWriter()
    for key, path in result["paths"].items():
        writer.append(str(path.with_suffix(".pdf")), outline_item=key.replace("_", " ").title())
    writer.add_metadata({"/Title": "Complete CTA source intake and curation"})
    writer.write(str(combined))
    writer.close()
    if len(PdfReader(combined).pages) != len(result["paths"]):
        raise ValueError("Combined PDF does not contain every figure")
    manifest_path = args.out / "run-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["outputs"][combined.name] = hashlib.sha256(combined.read_bytes()).hexdigest()
    manifest["python"] = platform.python_version()
    manifest["packages"] = {
        name: importlib.metadata.version(name)
        for name in ("numpy", "pandas", "matplotlib", "matplotlib-venn", "pypdf")
    }
    manifest["implementation"] = {
        name: hashlib.sha256((Path(__file__).resolve().parents[1] / name).read_bytes()).hexdigest()
        for name in (
            "oncoref/cta.py",
            "oncoref/cta_curation_plots.py",
            "oncoref/cta_landscape.py",
            "oncoref/cta_sources.py",
            "oncoref/cta_regen.py",
            "oncoref/cta_tissues.py",
            "scripts/render_cta_curation.py",
        )
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    index_path = args.out / "index.md"
    index_path.write_text(
        index_path.read_text().replace(
            "# CTA source intake and curation\n",
            f"# CTA source intake and curation\n\n[All figures, vector PDF]({combined.name})\n",
            1,
        )
    )
    print(combined)


if __name__ == "__main__":
    main()

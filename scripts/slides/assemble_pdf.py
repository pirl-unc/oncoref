#!/usr/bin/env python3
"""Package native slide renders with bookmarks and clickable source links.

This intentionally produces an image-based PDF. The adjacent native PPTX is
editable. It avoids the missing standard-font resource in artifact-tool's PDF
export in the original runtime. No classification or rendering logic lives here.
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

from pypdf import PdfReader
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


def assemble(build_dir: Path, output: Path) -> None:
    manifest = json.loads((build_dir / "manifest.json").read_text())
    pages = manifest["pages"]
    width, height = (manifest["slideSize"][key] for key in ("width", "height"))
    if not pages or width <= 0 or height <= 0:
        raise ValueError("A nonempty deck and positive page dimensions are required")
    images = [build_dir / "rendered" / f"slide-{n:02d}.png" for n in range(1, len(pages) + 1)]
    for image in images:
        if not image.is_file():
            raise FileNotFoundError(image)
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(suffix=".pdf", dir=output.parent)
    os.close(fd)
    staging = Path(name)
    try:
        pdf = canvas.Canvas(str(staging), pagesize=(width, height), pageCompression=1)
        pdf.setTitle(pages[0]["title"])
        pdf.setAuthor("OncoRef")
        pdf.setSubject(f"Policy {manifest['policy']}; source {manifest['commit']}")
        for n, (page, image) in enumerate(zip(pages, images), 1):
            pdf.drawImage(ImageReader(str(image)), 0, 0, width=width, height=height)
            key = f"slide-{n}"
            pdf.bookmarkPage(key)
            pdf.addOutlineEntry(f"{n:02d} | {page['title']}", key, level=0)
            for link in page["links"]:
                x, y, w, h = (link[key] for key in ("x", "y", "w", "h"))
                pdf.linkURL(link["url"], (x, height - y - h, x + w, height - y), thickness=0)
            pdf.showPage()
        pdf.save()
        reader = PdfReader(staging)
        if len(reader.pages) != len(pages) or len(reader.outline) != len(pages):
            raise ValueError("Page/bookmark count does not match the manifest")
        if sum(len(p.get("/Annots", [])) for p in reader.pages) != sum(
            len(p["links"]) for p in pages
        ):
            raise ValueError("Source links were lost during PDF assembly")
        staging.replace(output)
    finally:
        staging.unlink(missing_ok=True)
    print(f"{output}: {len(pages)} slides; bookmarks and source links verified")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    assemble(args.build_dir, args.output)

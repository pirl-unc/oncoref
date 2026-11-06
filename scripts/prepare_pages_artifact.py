"""Serve large historical CSV downloads from the exact repository revision.

Only the generated site is changed. Report sources, manifests and data remain
in the repository; interactive JSON data and current explorer downloads stay
on Pages. Run after ``mkdocs build`` and before uploading the Pages artifact.
"""

import argparse
import hashlib
import html
import json
import posixpath
import re
from pathlib import Path
from urllib.parse import quote, urlsplit

CURRENT_DOWNLOADS = "audits/cta-evidence-guide-20261106/"
MINIMUM_BYTES = 1_000_000
MAXIMUM_SITE_BYTES = 1_000_000_000
HREF = re.compile(r"(\bhref\s*=\s*)([\"'])(.*?)(\2)", re.IGNORECASE)


def download_path(url, page):
    """Resolve local report links, including absolute links to this Pages site."""
    parsed = urlsplit(html.unescape(url))
    if parsed.netloc and parsed.netloc != "pirl-unc.github.io":
        return None
    if parsed.scheme and parsed.scheme not in {"http", "https"}:
        return None
    if not parsed.path:
        return None
    if parsed.path.startswith("/oncoref/"):
        return parsed.path.removeprefix("/oncoref/")
    if parsed.netloc or parsed.path.startswith("/"):
        return None
    return posixpath.normpath(posixpath.join(page.parent.as_posix(), parsed.path))


def prepare(site_dir, docs_dir, revision):
    site_dir, docs_dir = Path(site_dir).resolve(), Path(docs_dir).resolve()
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ValueError("Use a complete Git commit SHA for immutable downloads")
    if site_dir == docs_dir or site_dir in docs_dir.parents or docs_dir in site_dir.parents:
        raise ValueError("The generated site and source directories must be separate")
    if not site_dir.is_dir() or not docs_dir.is_dir():
        raise ValueError("Build the site before preparing the Pages artifact")

    downloads = {}
    for source in sorted((docs_dir / "audits").rglob("*.csv*")):
        relative = source.relative_to(docs_dir).as_posix()
        built = site_dir / relative
        if (
            not source.is_file()
            or not relative.endswith((".csv", ".csv.gz"))
            or relative.startswith(CURRENT_DOWNLOADS)
            or source.stat().st_size < MINIMUM_BYTES
            or not built.is_file()
        ):
            continue
        if (
            hashlib.sha256(source.read_bytes()).digest()
            != hashlib.sha256(built.read_bytes()).digest()
        ):
            raise ValueError(f"Generated download differs from its source: {relative}")
        downloads[relative] = (
            f"https://raw.githubusercontent.com/pirl-unc/oncoref/{revision}/docs/"
            + quote(relative, safe="/")
        )

    links = 0
    for page in sorted(site_dir.rglob("*.html")):
        relative = page.relative_to(site_dir)

        def rewrite(match, page_relative=relative):
            nonlocal links
            path = download_path(match[3], page_relative)
            if path not in downloads:
                return match[0]
            parsed = urlsplit(html.unescape(match[3]))
            url = downloads[path]
            if parsed.query:
                url += "?" + parsed.query
            if parsed.fragment:
                url += "#" + parsed.fragment
            links += 1
            return match[1] + match[2] + html.escape(url, quote=True) + match[4]

        original = page.read_text(encoding="utf-8")
        rewritten = HREF.sub(rewrite, original)
        if rewritten != original:
            page.write_text(rewritten, encoding="utf-8")

    removed_bytes = 0
    for relative in downloads:
        built = site_dir / relative
        removed_bytes += built.stat().st_size
        built.unlink()
    total_bytes = sum(p.stat().st_size for p in site_dir.rglob("*") if p.is_file())
    result = {
        "revision": revision,
        "external_downloads": len(downloads),
        "rewritten_links": links,
        "removed_bytes": removed_bytes,
        "site_bytes": total_bytes,
    }
    print(json.dumps(result, sort_keys=True))
    if total_bytes > MAXIMUM_SITE_BYTES:
        raise ValueError("Generated Pages site still exceeds the 1 GB size limit")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-dir", default="site")
    parser.add_argument("--docs-dir", default="docs")
    parser.add_argument("--revision", required=True)
    args = parser.parse_args()
    prepare(args.site_dir, args.docs_dir, args.revision)

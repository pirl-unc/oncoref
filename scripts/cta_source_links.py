"""Readable source links for reports; artifact URLs remain separate provenance."""

from __future__ import annotations

from urllib.parse import urlparse


def information_url(url: str) -> str:
    """Prefer explanatory release pages over an immediate database-file download."""
    parsed = urlparse(url)
    if (
        parsed.hostname
        and parsed.hostname.endswith("proteinatlas.org")
        and parsed.path.startswith("/download/")
    ):
        return f"{parsed.scheme}://{parsed.netloc}/about/download"
    return url


def markdown_link(label: str, url: str) -> str:
    """Escape labels while retaining a readable title rather than a bare URL."""
    if not url or urlparse(url).scheme not in {"http", "https"}:
        raise ValueError(f"An explicit web source URL is required: {url!r}")
    label = str(label).replace("[", "\\[").replace("]", "\\]").replace("|", " / ")
    return f"[{label}]({url.replace(' ', '%20')})"


def linked_sources(frame):
    """Link source names, citations and paper titles in a display-only copy.

    CSV/API values stay plain text with independent URLs. Names link to an
    informative landing page; the separate Source column retains the artifact.
    """
    display = frame.copy()
    if "source_url" not in display:
        return display
    for column in ["name", "citation", "title"]:
        if column in display:
            display[column] = [
                markdown_link(label, information_url(url))
                for label, url in zip(display[column], display.source_url)
            ]
    display["source_url"] = display.source_url.map(lambda url: markdown_link("Source", url))
    return display

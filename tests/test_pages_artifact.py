import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "prepare_pages_artifact", Path("scripts/prepare_pages_artifact.py")
)
pages = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pages)
REVISION = "a" * 40


def test_large_archived_downloads_keep_exact_bytes_and_immutable_links(tmp_path):
    docs, site = tmp_path / "docs", tmp_path / "site"
    archive = "audits/old/cohort.csv.gz"
    current = pages.CURRENT_DOWNLOADS + "normal-top-contexts.csv.gz"
    files = {
        archive: b"x" * pages.MINIMUM_BYTES,
        current: b"y" * pages.MINIMUM_BYTES,
        next(iter(pages.ARCHIVED_REPORT_DOWNLOADS)): b"p" * pages.MINIMUM_BYTES,
        "audits/old/small.csv": b"a,b\n1,2\n",
        "audits/old/data.json.gz": b"z" * pages.MINIMUM_BYTES,
    }
    for base in (docs, site):
        for name, body in files.items():
            path = base / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(body)
    page = site / "audits/old/index.html"
    page.write_text(
        '<a href="cohort.csv.gz?download=1&amp;a=2#data">relative</a>'
        "<a href='/oncoref/audits/old/cohort.csv.gz'>absolute</a>"
        '<a href="https://pirl-unc.github.io/oncoref/audits/old/cohort.csv.gz">site</a>'
        '<a href="https://example.org/cohort.csv.gz">external</a>'
        '<a href="small.csv">small</a><a href="#data">anchor</a>'
    )

    result = pages.prepare(site, docs, REVISION)
    expected = f"https://raw.githubusercontent.com/pirl-unc/oncoref/{REVISION}/docs/{archive}"
    assert result["external_downloads"] == 1
    assert result["rewritten_links"] == 3
    assert expected + "?download=1&amp;a=2#data" in page.read_text()
    assert page.read_text().count(expected) == 3
    assert 'href="https://example.org/cohort.csv.gz"' in page.read_text()
    assert 'href="small.csv"' in page.read_text()
    assert 'href="#data"' in page.read_text()
    assert not (site / archive).exists()
    for name, body in files.items():
        assert (docs / name).read_bytes() == body
        if name != archive:
            assert (site / name).read_bytes() == body


def test_rejects_mutable_revision_and_source_directory(tmp_path):
    docs, site = tmp_path / "docs", tmp_path / "site"
    docs.mkdir()
    site.mkdir()
    with pytest.raises(ValueError, match="complete Git commit"):
        pages.prepare(site, docs, "main")
    with pytest.raises(ValueError, match="must be separate"):
        pages.prepare(docs, docs, REVISION)


def test_rejects_changed_data_and_oversized_site(tmp_path, monkeypatch):
    docs, site = tmp_path / "docs", tmp_path / "site"
    for base in (docs, site):
        (base / "audits/old").mkdir(parents=True)
        (base / "audits/old/data.csv").write_bytes(b"x" * pages.MINIMUM_BYTES)
    (site / "audits/old/data.csv").write_bytes(b"y" * pages.MINIMUM_BYTES)
    with pytest.raises(ValueError, match="differs from its source"):
        pages.prepare(site, docs, REVISION)
    assert (site / "audits/old/data.csv").exists()
    (site / "audits/old/data.csv").write_bytes(b"x" * pages.MINIMUM_BYTES)
    (site / "index.html").write_text("Still too large")
    monkeypatch.setattr(pages, "MAXIMUM_SITE_BYTES", 1)
    with pytest.raises(ValueError, match="exceeds the 1 GB"):
        pages.prepare(site, docs, REVISION)

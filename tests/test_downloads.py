import gzip
import io
import logging
import sys
import urllib.error
import zipfile
from types import SimpleNamespace

import datacache.download
import pytest
import requests

from oncoref import _downloads, catalog, cli, expression_builders, expression_source_adapters


@pytest.mark.parametrize("kind", ["gzip", "zip", "query-zip"])
def test_preserves_archive_bytes_with_temporary_destination(tmp_path, http_downloads, kind):
    content = b"Gene\tTPM\nA\t10\n"
    if kind == "gzip":
        archive = gzip.compress(content)
        url = "https://example.test/matrix.tsv.gz?download=1"
    else:
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w") as zipped:
            zipped.writestr("matrix.tsv", content)
        archive = buffer.getvalue()
        url = "https://example.test/matrix.ZIP"
        if kind == "query-zip":
            url = "https://example.test/download.php?file_name=data.zip"
    http_downloads(lambda url: archive)
    destination = tmp_path / "matrix.tsv.part"
    _downloads.fetch_file(url, destination, label="test v1")
    assert destination.read_bytes() == archive
    assert list(tmp_path.iterdir()) == [destination]


def test_cache_reuse_is_offline_and_describes_destination(tmp_path, http_downloads, capsys):
    destination = tmp_path / "reference.tsv"
    destination.write_bytes(b"previous")
    http_downloads(lambda url: pytest.fail("cache reuse must not make a request"))
    _downloads.fetch_file("https://example.test/data", destination, label="HPA v23")
    output = capsys.readouterr()
    assert output.out == ""
    assert f"using cached HPA v23 at {destination}" in output.err


def test_retries_report_attempts_then_publish(tmp_path, monkeypatch, caplog, capsys):
    attempts = []
    delays = []

    def get(url, **kwargs):
        attempts.append(url)
        response = requests.Response()
        response.url = url
        response.status_code = 503 if len(attempts) == 1 else 200
        response.raw = io.BytesIO(b"complete")
        response.headers["Content-Length"] = "8"
        return response

    monkeypatch.setattr(requests, "get", get)
    monkeypatch.setattr(datacache.download.time, "sleep", delays.append)
    with caplog.at_level(logging.WARNING):
        path = _downloads.fetch_file(
            "https://example.test/data", tmp_path / "matrix", label="LUAD v1"
        )
    assert path.read_bytes() == b"complete"
    assert len(attempts) == 2
    assert delays == [1.0]
    assert "attempt 1/3" in caplog.text
    assert "retrying" in caplog.text
    output = capsys.readouterr()
    assert not output.out
    assert "downloaded LUAD v1" in output.err and "MiB/s" in output.err


@pytest.mark.parametrize("status,attempt_count", [(404, 1), (503, 3)])
def test_http_failure_preserves_previous_file_and_cleans_staging(
    tmp_path, monkeypatch, status, attempt_count
):
    destination = tmp_path / "reference.tsv"
    destination.write_bytes(b"previous")
    attempts = []

    def get(url, **kwargs):
        attempts.append(url)
        response = requests.Response()
        response.url = url
        response.status_code = status
        response.raw = io.BytesIO(b"failure")
        return response

    monkeypatch.setattr(requests, "get", get)
    monkeypatch.setattr(datacache.download.time, "sleep", lambda delay: None)
    with pytest.raises(urllib.error.HTTPError) as error:
        _downloads.fetch_file("https://example.test/data", destination, label="test", force=True)
    assert error.value.code == status
    assert len(attempts) == attempt_count
    assert destination.read_bytes() == b"previous"
    assert list(tmp_path.iterdir()) == [destination]


def test_progress_modes_and_quiet_restore_context(monkeypatch):
    monkeypatch.setenv("ONCOREF_DOWNLOAD_PROGRESS", "auto")
    monkeypatch.setattr(sys.stderr, "isatty", lambda: False)
    monkeypatch.setitem(sys.modules, "IPython", SimpleNamespace(get_ipython=lambda: None))
    assert not _downloads._show_progress(True)
    monkeypatch.setattr(sys.stderr, "isatty", lambda: True)
    assert _downloads._show_progress(True)
    with _downloads.download_options(progress="never"):
        assert not _downloads._show_progress(True)
    with _downloads.download_options(progress="always", verbose=False):
        assert not _downloads._show_progress(True)
        assert not _downloads.is_verbose()
    assert _downloads.is_verbose()
    assert not _downloads._show_progress(False)
    monkeypatch.setattr(sys.stderr, "isatty", lambda: False)
    monkeypatch.setitem(
        sys.modules, "IPython", SimpleNamespace(get_ipython=lambda: SimpleNamespace(kernel=True))
    )
    assert _downloads._show_progress(True)


def test_forced_bar_and_unknown_total(monkeypatch, tmp_path, capsys):
    def get(url, **kwargs):
        response = requests.Response()
        response.status_code = 200
        response.url = url
        response.raw = io.BytesIO(b"data")
        # Deliberately no Content-Length: the display must not invent a total.
        return response

    monkeypatch.setattr(requests, "get", get)
    with _downloads.download_options(progress="always"):
        _downloads.fetch_file("https://example.test/data", tmp_path / "data", label="test")
    output = capsys.readouterr()
    assert not output.out
    assert "Downloading" in output.err
    assert "100%" not in output.err


@pytest.mark.parametrize("module", [expression_builders, expression_source_adapters])
def test_build_source_downloads_preserve_gzip_and_repair_empty_cache(
    module, tmp_path, http_downloads
):
    content = gzip.compress(b"source matrix")
    http_downloads(lambda url: content)
    path = tmp_path / "matrix.tsv.gz"
    path.touch()
    assert module._download("https://example.test/matrix.tsv.gz", path).read_bytes() == content


def test_cli_quiet_keeps_results_and_progress_override(monkeypatch, capsys):
    monkeypatch.setenv("ONCOREF_DOWNLOAD_PROGRESS", "always")

    def fetch(*args, **kwargs):
        _downloads.report("downloading HPA")
        assert not _downloads._show_progress(True)
        return ["hpa_rna_consensus"]

    monkeypatch.setattr(catalog, "fetch", fetch)
    assert cli.main(["data", "fetch", "hpa", "--quiet"]) == 0
    output = capsys.readouterr()
    assert output.out.strip() == "Downloaded: hpa_rna_consensus"
    assert not output.err
    assert _downloads.is_verbose()

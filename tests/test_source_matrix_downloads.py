# Licensed under the Apache License, Version 2.0.
"""Network-free concurrency and interrupted-transfer regressions for #590."""

import io
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from http.client import IncompleteRead
from pathlib import Path
from threading import Barrier, local
from unittest.mock import patch

import pandas as pd
import pytest

from oncoref import source_matrices as sm


@pytest.fixture
def matrix_bytes():
    return pd.DataFrame({"sample": [1.0, 2.0]}, index=["gene1", "gene2"]).to_parquet()


def _synchronized_first_check(destination, barrier):
    """Make both callers observe a cold cache before either can download."""
    real_exists = Path.exists
    checks = local()

    def exists(path):
        present = real_exists(path)
        if path == destination and not getattr(checks, "checked", False):
            checks.checked = True
            barrier.wait(timeout=20)
        return present

    return exists


@pytest.fixture
def destination(monkeypatch, tmp_path):
    path = tmp_path / "TEST.parquet"
    monkeypatch.setattr(sm, "local_path", lambda code: path)
    monkeypatch.setattr(sm, "release_url", lambda code: "https://example.invalid/TEST.parquet")
    return path


@pytest.mark.parametrize("method", ["fetch", "ensure"])
def test_simultaneous_threads_download_once(monkeypatch, destination, matrix_bytes, method):
    calls = []

    def urlopen(url):
        calls.append(url)
        return io.BytesIO(matrix_bytes)

    monkeypatch.setattr(sm.urllib.request, "urlopen", urlopen)
    with monkeypatch.context() as racing:
        racing.setattr(Path, "exists", _synchronized_first_check(destination, Barrier(2)))
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(getattr(sm, method), "TEST") for _ in range(2)]
            results = [future.result(timeout=25) for future in futures]

    assert results == [destination, destination]
    assert len(calls) == 1
    assert destination.read_bytes() == matrix_bytes
    # Both public paths must reuse the completed artifact without any HTTP.
    assert sm.fetch("TEST") == sm.ensure("TEST") == destination
    assert len(calls) == 1
    assert not list(destination.parent.glob("*.part"))


def _process_fetch(destination, barrier, calls, results, payload, method):
    """Spawn-safe worker: only HTTP is faked; locking and file I/O are real."""
    destination = Path(destination)

    def urlopen(url):
        calls.put(url)
        return io.BytesIO(payload)

    try:
        with (
            patch.object(sm, "local_path", return_value=destination),
            patch.object(sm, "release_url", return_value="https://example.invalid/TEST.parquet"),
            patch.object(sm.urllib.request, "urlopen", urlopen),
            patch.object(Path, "exists", _synchronized_first_check(destination, barrier)),
        ):
            result = getattr(sm, method)("TEST")
            results.put(("ok", str(result)))
    except Exception as exc:
        results.put(("error", repr(exc)))


@pytest.mark.parametrize("method", ["fetch", "ensure"])
def test_simultaneous_processes_download_once(tmp_path, matrix_bytes, method):
    destination = tmp_path / "TEST.parquet"
    ctx = multiprocessing.get_context("spawn")
    barrier = ctx.Barrier(2)
    calls = ctx.Queue()
    results = ctx.Queue()
    processes = [
        ctx.Process(
            target=_process_fetch,
            args=(str(destination), barrier, calls, results, matrix_bytes, method),
        )
        for _ in range(2)
    ]
    try:
        for process in processes:
            process.start()
        outcomes = [results.get(timeout=30) for _ in processes]
        for process in processes:
            process.join(timeout=10)
            assert process.exitcode == 0
        assert outcomes == [("ok", str(destination))] * 2
        assert calls.get(timeout=5) == "https://example.invalid/TEST.parquet"
        assert calls.empty()  # joined workers have flushed their queue feeders
        assert destination.read_bytes() == matrix_bytes
        assert not list(tmp_path.glob("*.part"))
    finally:
        for process in processes:
            if process.is_alive():
                process.terminate()
                process.join(timeout=10)
        calls.close()
        results.close()


@pytest.mark.parametrize("separate_versions", [False, True])
def test_different_destinations_can_download_concurrently(
    monkeypatch, tmp_path, matrix_bytes, separate_versions
):
    def local_path(code):
        path = (
            tmp_path / code / "TEST.parquet" if separate_versions else tmp_path / f"{code}.parquet"
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    monkeypatch.setattr(sm, "local_path", local_path)
    monkeypatch.setattr(sm, "release_url", lambda code: f"https://example.invalid/{code}")
    downloading = Barrier(2)

    def urlopen(url):
        downloading.wait(timeout=20)
        return io.BytesIO(matrix_bytes)

    monkeypatch.setattr(sm.urllib.request, "urlopen", urlopen)
    with ThreadPoolExecutor(max_workers=2) as pool:
        paths = list(pool.map(sm.fetch, ["FIRST", "SECOND"]))
    assert paths == [local_path("FIRST"), local_path("SECOND")]
    assert all(path.read_bytes() == matrix_bytes for path in paths)


def test_concurrent_force_calls_each_refresh(monkeypatch, destination, matrix_bytes):
    destination.write_bytes(matrix_bytes)
    calls = []

    def urlopen(url):
        calls.append(url)
        return io.BytesIO(matrix_bytes)

    monkeypatch.setattr(sm.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(Path, "exists", _synchronized_first_check(destination, Barrier(2)))
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(sm.fetch, "TEST", force=True) for _ in range(2)]
        assert [future.result(timeout=25) for future in futures] == [destination, destination]
    assert len(calls) == 2
    assert destination.read_bytes() == matrix_bytes
    assert not list(destination.parent.glob("*.part"))


def test_force_refresh_uses_unique_staging_and_atomic_replacement(
    monkeypatch, destination, matrix_bytes
):
    destination.write_bytes(b"previous matrix")
    stale = destination.with_suffix(".parquet.part")
    stale.write_bytes(b"unrelated interrupted writer")
    stages = []

    class Response(io.BytesIO):
        def read(self, size=-1):
            assert destination.read_bytes() == b"previous matrix"
            active = set(destination.parent.glob("*.part")) - {stale}
            assert len(active) == 1
            stage = active.pop()
            if stage not in stages:
                stages.append(stage)
            return super().read(size)

    monkeypatch.setattr(sm.urllib.request, "urlopen", lambda url: Response(matrix_bytes))
    for _ in range(2):
        destination.write_bytes(b"previous matrix")
        assert sm.fetch("TEST", force=True) == destination
        assert destination.read_bytes() == matrix_bytes
    assert len(set(stages)) == 2
    assert set(destination.parent.glob("*.part")) == {stale}
    assert stale.read_bytes() == b"unrelated interrupted writer"


@pytest.mark.parametrize("existing", [False, True])
@pytest.mark.parametrize("error_type", [OSError, IncompleteRead, KeyboardInterrupt])
def test_interrupted_transfer_cleans_staging_and_preserves_cache(
    monkeypatch, destination, matrix_bytes, existing, error_type
):
    if existing:
        destination.write_bytes(b"previous matrix")

    class InterruptedResponse(io.BytesIO):
        def read(self, size=-1):
            if self.tell():
                raise error_type(b"interrupted")
            return super().read(8)

    monkeypatch.setattr(sm.urllib.request, "urlopen", lambda url: InterruptedResponse(matrix_bytes))
    expected = KeyboardInterrupt if error_type is KeyboardInterrupt else sm.SourceMatrixError
    with pytest.raises(expected):
        sm.fetch("TEST", force=existing)
    assert not list(destination.parent.glob("*.part"))
    if existing:
        assert destination.read_bytes() == b"previous matrix"
    else:
        assert not destination.exists()

    # A failed writer releases its lock so the next attempt can complete.
    monkeypatch.setattr(sm.urllib.request, "urlopen", lambda url: io.BytesIO(matrix_bytes))
    assert sm.fetch("TEST", force=existing) == destination
    assert destination.read_bytes() == matrix_bytes


def test_failed_promotion_preserves_cache_and_cleans_staging(
    monkeypatch, destination, matrix_bytes
):
    destination.write_bytes(b"previous matrix")
    monkeypatch.setattr(sm.urllib.request, "urlopen", lambda url: io.BytesIO(matrix_bytes))

    def fail_replace(source, target):
        raise OSError("cannot promote")

    monkeypatch.setattr(Path, "replace", fail_replace)
    with pytest.raises(sm.SourceMatrixError, match="cannot promote"):
        sm.fetch("TEST", force=True)
    assert destination.read_bytes() == b"previous matrix"
    assert not list(destination.parent.glob("*.part"))

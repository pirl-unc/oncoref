"""Shared raw-file transfers and user-facing download status.

Source owners retain their own validation, extraction and cache contracts.
Datacache handles streaming, bounded retries and atomic transfer publication.
"""

from __future__ import annotations

import os
import stat
import sys
import tempfile
import time
import urllib.error
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path

import datacache
import requests

_OPTIONS: ContextVar[tuple[str | None, bool]] = ContextVar(
    "oncoref_download_options", default=(None, True)
)


@contextmanager
def download_options(*, progress: str | None = None, verbose: bool = True):
    """Temporarily configure CLI output without changing process-wide logging."""
    token = _OPTIONS.set((progress, verbose))
    try:
        yield
    finally:
        _OPTIONS.reset(token)


def is_verbose(verbose: bool = True) -> bool:
    return verbose and _OPTIONS.get()[1]


def report(message: str, *, verbose: bool = True) -> None:
    if is_verbose(verbose):
        print(f"oncoref: {message}", file=sys.stderr, flush=True)


def report_cached(label: str, path: Path, *, verbose: bool = True) -> None:
    report(f"using cached {label} at {path}", verbose=verbose)


def _show_progress(verbose: bool) -> bool:
    if not is_verbose(verbose):
        return False
    mode = _OPTIONS.get()[0] or os.environ.get("ONCOREF_DOWNLOAD_PROGRESS", "auto")
    if mode not in {"auto", "always", "never"}:
        raise ValueError("ONCOREF_DOWNLOAD_PROGRESS must be auto, always, or never")
    if mode != "auto":
        return mode == "always"
    if getattr(sys.stderr, "isatty", lambda: False)():
        return True
    # Detect an existing notebook without importing or starting IPython.
    ipython = sys.modules.get("IPython")
    shell = ipython.get_ipython() if ipython is not None else None
    return shell is not None and getattr(shell, "kernel", None) is not None


def fetch_file(
    url: str,
    destination: Path,
    *,
    label: str,
    force: bool = False,
    verbose: bool = True,
    timeout: float = 180,
    display_path: Path | None = None,
) -> Path:
    """Fetch raw bytes; leave an existing destination intact on transfer failure.

    Use datacache's inferred filename in a private sibling directory: an explicit
    ``.part``/``.tsv`` destination would implicitly decompress a ZIP/gzip source,
    even with ``decompress=False``. Oncoref must verify the original archive bytes
    and select archive members itself. No temporary cache survives the operation.
    HTTP exceptions retain the urllib interface used by the existing source APIs.
    """
    destination = Path(destination)
    if destination.is_file() and not force:
        report_cached(label, destination, verbose=verbose)
        return destination
    show_progress = _show_progress(verbose)
    destination.parent.mkdir(parents=True, exist_ok=True)
    report(
        f"downloading {label}\n  from {url}\n  to   {display_path or destination}",
        verbose=verbose,
    )
    started = time.monotonic()
    try:
        with tempfile.TemporaryDirectory(prefix=".download-", dir=destination.parent) as staging:
            path = Path(
                datacache.fetch_file(
                    url,
                    cache_root=staging,
                    decompress=False,
                    timeout=timeout,
                    max_retries=2,
                    show_progress=show_progress,
                )
            )
            # Preserve the permissions of a replaced cache file.
            if destination.exists():
                path.chmod(stat.S_IMODE(destination.stat().st_mode))
            path.replace(destination)
    except requests.HTTPError as error:
        response = error.response
        if response is not None:
            raise urllib.error.HTTPError(
                url, response.status_code, str(error), response.headers, None
            ) from error
        raise urllib.error.URLError(str(error)) from error
    except requests.RequestException as error:
        raise urllib.error.URLError(str(error)) from error
    elapsed = max(time.monotonic() - started, 0.001)
    size = destination.stat().st_size / (1024 * 1024)
    report(
        f"downloaded {label}: {size:.2f} MiB in {elapsed:.1f}s ({size / elapsed:.2f} MiB/s)",
        verbose=verbose,
    )
    return destination

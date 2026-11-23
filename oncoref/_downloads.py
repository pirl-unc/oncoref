"""Shared cache-root resolution, raw-file transfers, and download status (#747).

All OncoRef cache subsystems resolve their root through ``oncoref_cache_root()``
backed by ``datacache.get_cache_root``, honoring ``OPENVAX_DATA_CACHE`` and
``ONCOREF_CACHE_DIR`` while preserving existing ``~/.cache/oncoref`` caches
without forced re-downloads.
"""

from __future__ import annotations

import os
import stat
import sys
import time
import urllib.error
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path

import datacache
import requests

OPENVAX_CACHE_ENV_VAR = "OPENVAX_DATA_CACHE"
ONCOREF_CACHE_ENV_VAR = "ONCOREF_CACHE_DIR"
LEGACY_DATA_DIR_ENV_VAR = "CANCERDATA_DATA_DIR"
OFFLINE_ENV_VAR = "OPENVAX_OFFLINE"
ONCOREF_OFFLINE_ENV_VAR = "ONCOREF_OFFLINE"

_TRUTHY_ENV = frozenset({"1", "true", "yes", "on"})

_OPTIONS: ContextVar[tuple[str | None, bool, bool, Path | None]] = ContextVar(
    "oncoref_download_options",
    default=(None, True, False, None),
)


@contextmanager
def download_options(
    *,
    progress: str | None = None,
    verbose: bool = True,
    offline: bool | None = None,
    cache_root: str | os.PathLike | None = None,
):
    """Temporarily configure CLI download, offline, and cache-root options."""
    prev_progress, _prev_verbose, prev_offline, prev_root = _OPTIONS.get()
    resolved_root = Path(cache_root).expanduser() if cache_root is not None else prev_root
    token = _OPTIONS.set(
        (
            progress if progress is not None else prev_progress,
            verbose,
            prev_offline if offline is None else bool(offline),
            resolved_root,
        )
    )
    try:
        yield
    finally:
        _OPTIONS.reset(token)


def is_verbose(verbose: bool = True) -> bool:
    return verbose and _OPTIONS.get()[1]


def is_offline(offline: bool = False) -> bool:
    """Return True when offline mode is requested via flag, context, or environment."""
    if offline or _OPTIONS.get()[2]:
        return True
    for key in (OFFLINE_ENV_VAR, ONCOREF_OFFLINE_ENV_VAR):
        if os.environ.get(key, "").strip().lower() in _TRUTHY_ENV:
            return True
    return False


def active_cache_root_override() -> Path | None:
    """Return any scoped ``--cache-root`` override set via :func:`download_options`."""
    return _OPTIONS.get()[3]


def has_explicit_cache_override(*, envkeys: tuple[str, ...] = ()) -> bool:
    """Return True when a caller override or cache-root environment variable is active."""
    if active_cache_root_override() is not None:
        return True
    for key in (*envkeys, LEGACY_DATA_DIR_ENV_VAR, ONCOREF_CACHE_ENV_VAR, OPENVAX_CACHE_ENV_VAR):
        if os.environ.get(key, "").strip():
            return True
    return False


def resolve_cache_subdir(
    subpath: str,
    *,
    override: str | os.PathLike | None = None,
    envkeys: tuple[str, ...] = (),
) -> Path:
    """Resolve ``<cache_root>/<subpath>``, reusing ``~/.cache/oncoref/<subpath>`` on split-brain hosts."""
    if override is not None:
        return Path(override).expanduser()
    for key in envkeys:
        val = os.environ.get(key, "").strip()
        if val:
            return Path(val).expanduser()
    base = oncoref_cache_root()
    primary = base / subpath
    if not has_explicit_cache_override():
        legacy_sub = (Path.home() / ".cache" / "oncoref" / subpath).expanduser()
        if legacy_sub != primary and not primary.exists() and legacy_sub.exists():
            return legacy_sub
    return primary


def oncoref_cache_root(
    *,
    override: str | os.PathLike | None = None,
    envkeys: tuple[str, ...] = (),
) -> Path:
    """Resolve OncoRef's unified base cache directory without creating it.

    Precedence:
    1. Explicit ``override`` argument or active ``download_options(cache_root=...)``
    2. Subsystem-specific ``envkeys`` (when supplied)
    3. ``CANCERDATA_DATA_DIR`` or ``ONCOREF_CACHE_DIR`` (direct OncoRef root)
    4. ``OPENVAX_DATA_CACHE`` (shared OpenVax workspace root -> ``<root>/oncoref``)
    5. ``datacache.get_cache_root("oncoref", legacy=(~/.cache/oncoref,))``
       (platform cache dir, reusing an existing ``~/.cache/oncoref`` when populated)
    """
    effective_override = override if override is not None else active_cache_root_override()
    if effective_override is not None:
        return Path(effective_override).expanduser()
    for key in (*envkeys, LEGACY_DATA_DIR_ENV_VAR, ONCOREF_CACHE_ENV_VAR):
        val = os.environ.get(key, "").strip()
        if val:
            return Path(val).expanduser()
    openvax_root = os.environ.get(OPENVAX_CACHE_ENV_VAR, "").strip()
    if openvax_root:
        return Path(openvax_root).expanduser() / "oncoref"
    legacy_dot_cache = (Path.home() / ".cache" / "oncoref").expanduser()
    legacy_dirs = (str(legacy_dot_cache),) if legacy_dot_cache.is_absolute() else ()
    return Path(datacache.get_cache_root("oncoref", legacy=legacy_dirs))


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
    """Fetch raw bytes directly via ``datacache.fetch_file(..., raw=True)``.

    Preserves archive bytes regardless of output suffix, leaves an existing
    destination intact on transfer failure, honors offline mode, and converts
    HTTP errors to ``urllib.error`` for compatibility with source callers.
    """
    destination = Path(destination)
    if destination.is_file() and not force:
        report_cached(label, destination, verbose=verbose)
        return destination
    if is_offline():
        raise FileNotFoundError(f"Offline mode is enabled; refusing to download {label} from {url}")
    show_progress = _show_progress(verbose)
    destination.parent.mkdir(parents=True, exist_ok=True)
    prior_mode = stat.S_IMODE(destination.stat().st_mode) if destination.exists() else None
    report(
        f"downloading {label}\n  from {url}\n  to   {display_path or destination}",
        verbose=verbose,
    )
    started = time.monotonic()
    try:
        datacache.fetch_file(
            url,
            destination=destination,
            raw=True,
            force=True,
            timeout=timeout,
            max_retries=2,
            show_progress=show_progress,
        )
        if prior_mode is not None:
            destination.chmod(prior_mode)
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

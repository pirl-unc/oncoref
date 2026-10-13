# Licensed under the Apache License, Version 2.0.
"""Content-aware CSV/TSV reads for local source and cache files."""

from __future__ import annotations

import bz2
import gzip
import lzma
import tarfile
from pathlib import Path

import pandas as pd


def _is_tar_header(header):
    try:
        tarfile.TarInfo.frombuf(header, encoding="utf-8", errors="surrogateescape")
    except tarfile.HeaderError:
        return False
    return True


def _compression(path):
    with path.open("rb") as stream:
        header = stream.read(512)
    # Check inside compressed streams as well: pandas treats a single-member
    # tar.gz/bz2/xz as an archive, rather than parsing its tar header as table data.
    for magic, method, opener in (
        (b"\x1f\x8b", "gzip", gzip.open),
        (b"BZh", "bz2", bz2.open),
        (b"\xfd7zXZ\x00", "xz", lzma.open),
    ):
        if header.startswith(magic):
            with opener(path, "rb") as stream:
                return "tar" if _is_tar_header(stream.read(512)) else method
    if header.startswith((b"PK\x03\x04", b"PK\x05\x06", b"PK\x07\x08")):
        return "zip"
    if header.startswith(b"\x28\xb5\x2f\xfd"):
        return "zstd"
    return "tar" if _is_tar_header(header) else None


def read_csv(path, **kwargs):
    """Read a local table using its bytes, independent of its name or symlink.

    Supports plain text, gzip, bzip2, xz, Zstandard (with pandas' optional
    dependency), and pandas' single-member zip/tar archives. Parsing options and
    chunking are passed through unchanged. Pandas owns its file handles, so a
    returned chunk reader remains usable after this function returns.

    Callers remain responsible for any required asset checksum verification.
    """
    path = Path(path)
    return pd.read_csv(path, compression=_compression(path), **kwargs)

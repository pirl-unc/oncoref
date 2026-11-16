"""Exact bytes of files that dated receipts bound, after the live file moved on.

A receipt names a repository path and a hash. If the live file still has that hash
it is used. Otherwise the bytes come from a frozen fixture set that records the same
path with the same hash. Nothing is accepted on the strength of its name alone.
"""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = [
    ROOT / "tests/fixtures/release-1.39-inputs",
    ROOT / "tests/fixtures/release-1.38-inputs",
    ROOT / "tests/fixtures/release-1.37-inputs",
    ROOT / "tests/fixtures/cta-1.31-inputs",
    ROOT / "tests/fixtures/cta-1.21-inputs",
]


def _recorded(folder: Path) -> dict:
    receipt = json.loads((folder / "manifest.json").read_text())
    if "files" in receipt:
        return {name: entry["sha256"] for name, entry in receipt["files"].items()}
    return receipt["sha256"]


def path(name: str, expected: str) -> Path:
    """Live path if it still has the bound hash, else the frozen copy that does."""
    current = ROOT / name
    if current.is_file() and hashlib.sha256(current.read_bytes()).hexdigest() == expected:
        return current
    for folder in FIXTURES:
        if _recorded(folder).get(name) == expected:
            frozen = folder / name
            assert hashlib.sha256(frozen.read_bytes()).hexdigest() == expected, name
            return frozen
    raise AssertionError(f"{name} no longer has the bound hash and no frozen copy records it")


def read(name: str, expected: str) -> bytes:
    return path(name, expected).read_bytes()

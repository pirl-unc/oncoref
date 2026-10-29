"""Read historical curated data from verified wheels without installing old code.

Only the requested data member is extracted. A package release identifies its
curation, not the upstream release of each study represented in the table.
"""

from __future__ import annotations

import json
import re
import tempfile
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

from filelock import FileLock

from . import _downloads, reference_data
from .version import __version__


def resolve_version(version: str | None = None) -> str:
    """Concrete package release; latest is the installed registered release."""
    if version in (None, "latest"):
        return f"oncoref-{__version__}"
    if not re.fullmatch(r"oncoref-[0-9]+\.[0-9]+\.[0-9]+", version):
        raise ValueError(f"Invalid packaged-data version {version!r}; use oncoref-X.Y.Z")
    return version


def local_path(relative: str, version: str) -> Path:
    """Version-isolated path for a curated data member, never executable code."""
    member = Path(relative)
    if member.is_absolute() or ".." in member.parts or not member.parts:
        raise ValueError("Packaged data member must stay within oncoref/data")
    return reference_data.cache_dir() / "packaged" / resolve_version(version) / member


def provenance(relative: str, version: str) -> dict:
    """Cached wheel provenance; missing historical releases have no fabricated hash."""
    dest = local_path(relative, version)
    receipt = dest.with_name(dest.name + ".provenance.json")
    try:
        return json.loads(receipt.read_text()) if receipt.exists() else {}
    except (OSError, json.JSONDecodeError):
        return {}


def ensure(relative: str, version: str) -> Path:
    """Fetch one data member from an exact PyPI wheel, checking the wheel digest.

    A missing release or member fails explicitly. No code is imported, installed
    or executed; older data never reuse the installed package's table.
    """
    concrete = resolve_version(version)
    dest = local_path(relative, concrete)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with FileLock(str(dest.with_name(dest.name + ".lock"))):
        record = provenance(relative, concrete)
        if dest.is_file() and record.get("sha256") == reference_data._sha256(dest):
            return dest
        release = concrete.removeprefix("oncoref-")
        metadata_url = f"https://pypi.org/pypi/oncoref/{release}/json"
        try:
            with urllib.request.urlopen(metadata_url, timeout=60) as response:
                metadata = json.load(response)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                raise ValueError(f"No published OncoRef package {release}") from error
            raise
        wheels = [row for row in metadata["urls"] if row["filename"].endswith("-py3-none-any.whl")]
        if len(wheels) != 1:
            raise ValueError(f"Expected one platform-independent wheel for {concrete}")
        wheel = wheels[0]
        with tempfile.TemporaryDirectory(dir=dest.parent) as directory:
            archive = Path(directory) / "source.whl"
            _downloads.fetch_file(wheel["url"], archive, label=concrete, force=True)
            expected = wheel["digests"]["sha256"]
            if reference_data._sha256(archive) != expected:
                raise ValueError(f"Wheel checksum mismatch for {concrete}")
            member = "oncoref/data/" + Path(relative).as_posix()
            with zipfile.ZipFile(archive) as zipped:
                if member not in zipped.namelist():
                    raise FileNotFoundError(f"{relative!r} is absent from {concrete}")
                staged = Path(directory) / "member"
                staged.write_bytes(zipped.read(member))
            record = {
                "version": concrete,
                "url": wheel["url"],
                "archive_sha256": expected,
                "sha256": reference_data._sha256(staged),
                "member": member,
            }
            receipt = Path(directory) / "receipt.json"
            receipt.write_text(json.dumps(record, indent=2) + "\n")
            staged.replace(dest)
            receipt.replace(dest.with_name(dest.name + ".provenance.json"))
    return dest

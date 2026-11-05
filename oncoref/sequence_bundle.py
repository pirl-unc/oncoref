"""Verified, offline sequence records for the explicitly versioned CTA policy.

Acquisition is opt-in. Original records, every occurrence and interrupted valid
segments remain available; exact sequence grouping never discards an exclusion.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import re
import tempfile
import zipfile
from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path

SCHEMA = "oncoref.sequence-bundle"
PAYLOADS = frozenset(
    {"records.jsonl.gz", "coverage.json", "decisions.csv.gz", "orfdata-proteome.json.gz"}
)
ROLE_ORDER = {"allow": 0, "review": 1, "exclude": 2}


def digest(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def content_id(manifest: dict) -> str:
    body = {k: v for k, v in manifest.items() if k != "bundle_id"}
    return "oncoref:sequence-bundle:v1:sha256:" + digest(
        json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    )


def descriptor() -> dict:
    """Packaged current asset pin; reading it performs no download or cache write."""
    return json.loads(files("oncoref").joinpath("data/sequence-bundle.json").read_text())


def _decode_records(payload):
    with gzip.GzipFile(fileobj=io.BytesIO(payload)) as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


@dataclass(frozen=True)
class SequenceBundle:
    """Immutable verified bytes; public objects are detached on every access."""

    _manifest: bytes
    _records: bytes
    _coverage: bytes
    _orfdata: bytes

    @property
    def identity(self):
        return json.loads(self._manifest)["bundle_id"]

    @property
    def manifest(self):
        return json.loads(self._manifest)

    @property
    def coverage(self):
        return json.loads(self._coverage)

    def records(self):
        """Yield full originals with qualification, roles and all source occurrences."""
        yield from _decode_records(self._records)

    def segments(self, *, role=None, min_length=1):
        """Yield coordinate-bearing runs; never concatenate across an interruption."""
        if role is not None and role not in ROLE_ORDER:
            raise ValueError("Unknown sequence role")
        if type(min_length) is not int or min_length < 1:
            raise ValueError("min_length must be a positive integer")
        for record in self.records():
            if role is not None and record["shared_sequence"] != role:
                continue
            for segment in record["valid_segments"]:
                if segment["end"] - segment["start"] >= min_length:
                    yield {
                        **segment,
                        "record_id": record["record_id"],
                        "shared_sequence": record["shared_sequence"],
                        "occurrences": deepcopy(record["occurrences"]),
                        "reference_id": record["reference_id"],
                    }

    def groups(self):
        """Exact original-sequence groups; retain records/occurrences and role precedence.

        Full originals containing ambiguity are grouped only by their exact raw
        bytes; use segments for matching. Group identity does not prove production.
        """
        grouped = defaultdict(list)
        for record in self.records():
            grouped[record["sequence_sha256"]].append(record)
        for key, records in sorted(grouped.items()):
            yield {
                "sequence_sha256": key,
                "sequence": records[0]["sequence"],
                "shared_sequence": max((r["shared_sequence"] for r in records), key=ROLE_ORDER.get),
                "records": records,
                "record_count": len(records),
                "occurrence_count": sum(len(r["occurrences"]) for r in records),
            }

    def as_orfdata(self):
        """Load the preserved native ORFdata proteome with its original identity.

        Requires the optional quantitative extra. Policy roles remain in this
        bundle; the returned raw Proteome does not itself apply CTA exclusions.
        """
        try:
            from orfdata import Proteome
        except ImportError as error:
            raise ImportError("Install oncoref[quantitative] to use the ORFdata adapter") from error
        with tempfile.TemporaryDirectory(prefix="oncoref-orfdata-") as folder:
            path = Path(folder) / "proteome.json.gz"
            path.write_bytes(self._orfdata)
            proteome = Proteome.load(path)
        if proteome.identity != self.manifest["orfdata_identity"]:
            raise ValueError("Native ORFdata content identity mismatch")
        return proteome


def load(path=None, *, expected_id=None, require_current=True) -> SequenceBundle:
    """Verify an unpacked bundle without network access; default to the current pin.

    For a named historical bundle, explicitly pass its trusted expected_id and
    require_current=False. A manifest from the bundle itself is not a trusted pin.
    Missing coverage cannot be represented by an empty record list.
    """
    pin = descriptor() if require_current else None
    if not require_current and not expected_id:
        raise ValueError("Historical loading requires an explicit trusted bundle identity")
    if pin and expected_id and expected_id != pin["bundle_id"]:
        raise ValueError("Expected bundle is not the current package policy")
    expected_id = pin["bundle_id"] if pin else expected_id
    folder = Path(path) if path is not None else ensure(download=False)
    manifest_bytes = (folder / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    if (manifest.get("schema"), manifest.get("schema_version")) != (SCHEMA, 1):
        raise ValueError("Unsupported sequence bundle schema")
    if manifest.get("bundle_id") != expected_id or content_id(manifest) != expected_id:
        raise ValueError("Sequence bundle identity mismatch")
    if set(manifest["payloads"]) != PAYLOADS:
        raise ValueError("Missing or unexpected sequence bundle payload")
    if pin:
        for name, expected in manifest["package_dependencies"].items():
            if Path(name).name != name:
                raise ValueError("Invalid package dependency name")
            if digest(files("oncoref").joinpath("data", name).read_bytes()) != expected:
                raise ValueError(f"Mixed package/bundle policy: {name}")
    payloads = {}
    for name, entry in manifest["payloads"].items():
        payload = (folder / name).read_bytes()
        if len(payload) != entry["bytes"] or digest(payload) != entry["sha256"]:
            raise ValueError(f"Sequence bundle payload checksum mismatch: {name}")
        payloads[name] = payload
    coverage = json.loads(payloads["coverage.json"])
    if coverage.get("reference_id") != manifest["reference"]["reference_id"]:
        raise ValueError("Coverage reference mismatch")
    ids, occurrences = set(), set()
    for record in _decode_records(payloads["records.jsonl.gz"]):
        rid = record["record_id"]
        if rid in ids or record["reference_id"] != coverage["reference_id"]:
            raise ValueError("Repeated record or mixed reference")
        ids.add(rid)
        if not record["occurrences"] or record["shared_sequence"] not in ROLE_ORDER:
            raise ValueError("Record lacks an occurrence or explicit role")
        sequence = record["sequence"]
        if digest(sequence.encode()) != record["sequence_sha256"]:
            raise ValueError("Original sequence checksum mismatch")
        previous = 0
        for segment in record["valid_segments"]:
            start, end = segment["start"], segment["end"]
            if not 0 <= previous <= start < end <= len(sequence):
                raise ValueError("Invalid segment coordinates")
            value = sequence[start:end]
            if (
                value != segment["sequence"]
                or set(value) - set("ACDEFGHIKLMNPQRSTVWY")
                or digest(value.encode()) != segment["sha256"]
            ):
                raise ValueError("Invalid segment sequence")
            previous = end
        expected_segments = [
            (m.start(), m.end()) for m in re.finditer(r"[ACDEFGHIKLMNPQRSTVWY]+", sequence)
        ]
        if [(s["start"], s["end"]) for s in record["valid_segments"]] != expected_segments:
            raise ValueError("Valid segments do not preserve every complete AA20 run")
        for occurrence in record["occurrences"]:
            oid = occurrence["occurrence_id"]
            if oid in occurrences:
                raise ValueError("Duplicate occurrence identity")
            occurrences.add(oid)
            if occurrence.get("shared_sequence") not in ROLE_ORDER:
                raise ValueError("Occurrence lacks an explicit sequence role")
        role = max((o["shared_sequence"] for o in record["occurrences"]), key=ROLE_ORDER.get)
        if record["shared_sequence"] != role:
            raise ValueError("Record lost occurrence role precedence")
    if len(ids) != manifest["record_count"] or len(occurrences) != manifest["occurrence_count"]:
        raise ValueError("Sequence bundle count mismatch")
    return SequenceBundle(
        manifest_bytes,
        payloads["records.jsonl.gz"],
        payloads["coverage.json"],
        payloads["orfdata-proteome.json.gz"],
    )


def ensure(*, cache_dir=None, download=False) -> Path:
    """Acquire only the one pinned ZIP, bounded by its declared byte limits.

    No source discovery, gene downloads or automatic repairs. Corrupt caches fail;
    download=True explicitly permits replacing them with a verified pinned asset.
    """
    from filelock import FileLock

    pin = descriptor()
    root = Path(
        cache_dir
        or os.environ.get(
            "ONCOREF_SEQUENCE_CACHE_DIR", Path.home() / ".cache/oncoref/sequence-bundles"
        )
    )
    destination = root / pin["bundle_id"].rsplit(":", 1)[-1]
    if destination.exists():
        try:
            load(destination)
            return destination
        except (ValueError, OSError):
            if not download:
                raise
    if not download:
        raise FileNotFoundError(
            "Pinned sequence bundle is not cached; use ensure(download=True) explicitly"
        )
    root.mkdir(parents=True, exist_ok=True)
    with FileLock(str(root / (destination.name + ".lock"))):
        import shutil

        import requests

        with tempfile.TemporaryDirectory(prefix=".sequence-", dir=root) as temporary:
            temp = Path(temporary)
            archive = temp / "asset.zip"
            size = 0
            with requests.get(pin["url"], stream=True, timeout=(20, 120)) as response:
                response.raise_for_status()
                with archive.open("wb") as output:
                    for block in response.iter_content(1024 * 1024):
                        size += len(block)
                        if size > pin["archive_bytes"]:
                            raise ValueError("Sequence download exceeds declared bound")
                        output.write(block)
            if (
                size != pin["archive_bytes"]
                or digest(archive.read_bytes()) != pin["archive_sha256"]
            ):
                raise ValueError("Sequence archive checksum mismatch")
            unpacked = temp / "bundle"
            unpacked.mkdir()
            with zipfile.ZipFile(archive) as zipped:
                names = zipped.namelist()
                if set(names) != PAYLOADS | {"manifest.json"} or len(names) != len(set(names)):
                    raise ValueError("Unexpected sequence archive members")
                if sum(i.file_size for i in zipped.infolist()) > pin["unpacked_bytes"]:
                    raise ValueError("Sequence archive exceeds expansion bound")
                for name in names:
                    (unpacked / name).write_bytes(zipped.read(name))
            load(unpacked)
            if destination.exists():
                shutil.rmtree(destination)
            unpacked.rename(destination)
    return destination

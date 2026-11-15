# Licensed under the Apache License, Version 2.0.
"""Reconstruct the released bytes of packaged data files that were later renamed.

Released receipts and the pinned sequence bundle bind data files by name and hash.
When such a file is renamed or reworded without a change of content,
``data-file-renames.json`` records the exact text substitutions. :func:`as_released`
reverses them and checks the result against the released hash, so a bound file is
still verified byte for byte. An edit that is not a listed substitution fails here.
"""

from __future__ import annotations

import hashlib
import json
from functools import cache
from importlib.resources import files

_DATA = files("oncoref") / "data"


@cache
def renames() -> dict:
    """Released file name -> current name, released hash and substitutions."""
    doc = json.loads((_DATA / "data-file-renames.json").read_text())
    if (doc.get("schema"), doc.get("schema_version")) != ("oncoref.data-file-renames", 1):
        raise ValueError("Unsupported data-file rename schema")
    return doc["files"]


def current_name(name: str) -> str:
    """Name under which a released data file is packaged now."""
    return renames().get(name, {}).get("current_name", name)


def as_released(name: str) -> bytes:
    """Bytes of the data file as released under ``name``; current bytes if never renamed."""
    record = renames().get(name)
    if record is None:
        return (_DATA / name).read_bytes()
    text = (_DATA / record["current_name"]).read_bytes().decode()
    for item in record["substitutions"]:
        if text.count(item["current"]) != item["occurrences"]:
            raise ValueError(f"{record['current_name']} no longer matches its rename record")
        text = text.replace(item["current"], item["as_built"])
    payload = text.encode()
    if hashlib.sha256(payload).hexdigest() != record["as_built_sha256"]:
        raise ValueError(f"{record['current_name']} differs from released {name}")
    return payload

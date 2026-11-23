#!/usr/bin/env python3
"""Snapshot live code/script files bound by dated docs/audits manifests (#753).

Copies any repository ``oncoref/`` (excluding ``oncoref/data/``) or ``scripts/``
file whose SHA-256 is pinned in ``docs/audits/**/*manifest*.json`` and not yet
archived in an earlier ``tests/fixtures/*-inputs`` directory into a frozen
release-inputs directory (default ``tests/fixtures/release-1.49-inputs``).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST_SECTIONS = (
    "implementation_sha256",
    "inputs",
    "inputs_sha256",
    "implementation_and_inputs",
    "implementation",
)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _existing_recorded(fixtures_dir: Path, exclude_folder: Path) -> dict[str, set[str]]:
    recorded: dict[str, set[str]] = {}
    for mpath in sorted(fixtures_dir.glob("*-inputs/manifest.json")):
        if mpath.parent.resolve() == exclude_folder.resolve():
            continue
        doc = json.loads(mpath.read_text(encoding="utf-8"))
        mapping = (
            {k: v["sha256"] for k, v in doc["files"].items()}
            if "files" in doc
            else doc.get("sha256", {})
        )
        for rel, sha in mapping.items():
            recorded.setdefault(rel, set()).add(sha)
    return recorded


def snapshot(
    root: Path = ROOT,
    *,
    dest: Path | None = None,
    source_ref: str = "v1.49.0",
) -> dict[str, dict[str, int | str]]:
    if dest is None:
        dest = root / "tests/fixtures/release-1.49-inputs"
    prior = _existing_recorded(root / "tests/fixtures", dest)

    bound: dict[str, tuple[str, Path]] = {}
    for mpath in sorted((root / "docs/audits").rglob("*manifest*.json")):
        try:
            doc = json.loads(mpath.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        for sec_name in MANIFEST_SECTIONS:
            sec = doc.get(sec_name)
            if not isinstance(sec, dict):
                continue
            for rel, val in sec.items():
                sha = val.get("sha256") if isinstance(val, dict) else val
                if not (isinstance(sha, str) and len(sha) == 64):
                    continue
                if not (rel.startswith("oncoref/") or rel.startswith("scripts/")):
                    continue
                if rel.startswith("oncoref/data/"):
                    continue
                if sha in prior.get(rel, set()):
                    continue
                cur = root / rel
                if cur.is_file() and _sha256(cur) == sha:
                    bound[rel] = (sha, cur)

    dest.mkdir(parents=True, exist_ok=True)
    files_manifest: dict[str, dict[str, int | str]] = {}
    for rel in sorted(bound):
        sha, src = bound[rel]
        out_file = dest / rel
        out_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, out_file)
        files_manifest[rel] = {"sha256": sha, "bytes": out_file.stat().st_size}

    payload = {
        "source_ref": source_ref,
        "scope": "Frozen v1.49.0 code and script inputs bound by dated docs/audits manifests (#753).",
        "files": files_manifest,
    }
    (dest / "manifest.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return files_manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dest",
        type=Path,
        default=ROOT / "tests/fixtures/release-1.49-inputs",
    )
    parser.add_argument("--source-ref", default="v1.49.0")
    args = parser.parse_args(argv)
    files = snapshot(ROOT, dest=args.dest, source_ref=args.source_ref)
    total = sum(int(v["bytes"]) for v in files.values())
    print(f"Snapshotted {len(files)} files ({total:,} bytes) -> {args.dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

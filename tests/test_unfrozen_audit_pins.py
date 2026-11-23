"""Tests verifying that dated audit manifests stay decoupled from live source edits (#753)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from tests import released_inputs

ROOT = Path(__file__).resolve().parents[1]


def test_release_149_fixture_manifest_covers_target_selection_rebuild_sources():
    rebuild_manifest = json.loads(
        (ROOT / "docs/audits/target-selection-20261110/rebuild/manifest.json").read_text()
    )
    fixture_dir = ROOT / "tests/fixtures/release-1.49-inputs"
    fixture_manifest = json.loads((fixture_dir / "manifest.json").read_text())["files"]

    for rel, expected_sha in rebuild_manifest["implementation_sha256"].items():
        # Either recorded in release-1.49-inputs or an earlier fixture folder
        resolved = released_inputs.path(rel, expected_sha)
        assert hashlib.sha256(resolved.read_bytes()).hexdigest() == expected_sha
        if rel in fixture_manifest:
            frozen = fixture_dir / rel
            assert frozen.is_file()
            assert hashlib.sha256(frozen.read_bytes()).hexdigest() == expected_sha


def test_released_inputs_falls_back_when_live_runtime_module_evolves(monkeypatch, tmp_path):
    rebuild_manifest = json.loads(
        (ROOT / "docs/audits/target-selection-20261110/rebuild/manifest.json").read_text()
    )
    fake_root = tmp_path / "repo"
    (fake_root / "oncoref").mkdir(parents=True)
    # Simulate refactored live files whose digests differ from the 2026-11-10 audit
    for rel in ("oncoref/target_selection.py", "oncoref/host_context.py", "oncoref/cta.py"):
        (fake_root / rel).write_text("# refactored runtime implementation\n", encoding="utf-8")

    monkeypatch.setattr(released_inputs, "ROOT", fake_root)
    for rel in ("oncoref/target_selection.py", "oncoref/host_context.py"):
        expected = rebuild_manifest["implementation_sha256"][rel]
        content = released_inputs.read(rel, expected)
        assert hashlib.sha256(content).hexdigest() == expected

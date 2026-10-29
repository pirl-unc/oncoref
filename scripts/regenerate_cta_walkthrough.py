#!/usr/bin/env python3
"""Regenerate current decisions, packaged APIs, report and scientific figures.

Run the slide builder afterward. Frozen release artifacts are never overwritten.
Use --reuse-contexts only when the full-context replay receipt still verifies.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import audit_cta_source_coverage
import build_cta_selection_audit
import build_cta_selection_release
import render_cta_source_walkthrough
from cta_report_common import input_hashes, sha256
from replay_cta_selection import OUTPUT, ROOT, build


def receipts() -> None:
    """Record the measurement replay and current exports without self-hashing."""
    inputs = [
        ROOT / "oncoref/cta_selection_policy.py",
        ROOT / "oncoref/gene_matching.py",
        OUTPUT / "alternative-evidence/manifest.json",
        OUTPUT / "tumor-identity-repair/manifest.json",
        ROOT / "oncoref/data/cta-ihc-compartments.csv",
        ROOT / "scripts/replay_cta_selection.py",
        ROOT / "docs/audits/cta-release-1.21.0/manifest.json",
        ROOT / "docs/audits/cta-quantitative-20261026/developmental-lineages.csv.gz",
        ROOT
        / "docs/audits/cta-harmonized-20261027/group-normal/developmental-lineage-bounds.csv.gz",
        ROOT / "docs/audits/external-ihc-20261027/assessed-evidence.csv.gz",
    ]
    (OUTPUT / "manifest.json").write_text(
        json.dumps(
            {
                "inputs": input_hashes(inputs),
                "outputs": {
                    p.name: sha256(p)
                    for p in sorted(OUTPUT.iterdir())
                    if p.is_file()
                    and p.name not in {"manifest.json", "species-source-inventory.csv"}
                },
            },
            indent=2,
        )
        + "\n"
    )


def walkthrough_receipt() -> None:
    """Bind the final reports/figures/slides to their code and packaged evidence."""
    destination = build_cta_selection_audit.OUTPUT
    inputs = [
        OUTPUT / "manifest.json",
        ROOT / "oncoref/data/cta-selection.json",
        ROOT / "oncoref/data/cta-selection-audit.json",
        ROOT / "scripts/build_cta_selection_audit.py",
        ROOT / "scripts/build_cta_selection_release.py",
        ROOT / "scripts/render_cta_source_walkthrough.py",
        ROOT / "scripts/build_cta_source_slides.mjs",
        ROOT / "oncoref/data/reference-releases.json",
        ROOT / "docs/data-versions.md",
        ROOT / "scripts/audit_cta_embryonic_coverage.py",
        ROOT / "scripts/audit_cta_source_coverage.py",
        ROOT / "scripts/build_cta_alternative_evidence.py",
        ROOT / "scripts/import_cta_related_sources.py",
        ROOT / "oncoref/data/cta-publication-membership.csv",
        ROOT / "oncoref/data/cta-publication-context.csv",
        Path(__file__),
    ]
    receipt = {
        "inputs": input_hashes(inputs),
        "outputs": {
            str(path.relative_to(destination)): sha256(path)
            for path in sorted(destination.rglob("*"))
            if path.is_file() and path.name != "manifest.json"
        },
    }
    (destination / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")


def main() -> None:
    """Run dependent generation steps in evidence-to-presentation order."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse-contexts", action="store_true")
    parser.add_argument(
        "--verify-replay",
        action="store_true",
        help="Compare complete and cached replay outputs byte for byte",
    )
    parser.add_argument(
        "--finalize-receipts",
        action="store_true",
        help="Refresh the walkthrough receipt after finalizing the slides",
    )
    args = parser.parse_args()
    if args.finalize_receipts:
        walkthrough_receipt()
        return
    if args.verify_replay and args.reuse_contexts:
        parser.error("--verify-replay requires a complete first replay")
    build(reuse_contexts=args.reuse_contexts)
    compared = [
        "protein-group-classifications.csv.gz",
        "filtering-funnel.csv",
        "known-cta-decisions.csv.gz",
        "summary.json",
    ]
    if args.verify_replay:
        before = {name: sha256(OUTPUT / name) for name in compared}
        build(reuse_contexts=True)
        after = {name: sha256(OUTPUT / name) for name in compared}
        if before != after:
            raise ValueError("Complete and cached replay outputs differ")
        (OUTPUT / "reproducibility.json").write_text(
            json.dumps({"complete_vs_cached_replay": "byte_identical", "outputs": after}, indent=2)
            + "\n"
        )
        print("Complete and cached replay outputs are byte-identical", flush=True)
    else:
        check = OUTPUT / "reproducibility.json"
        if check.exists() and json.loads(check.read_text())["outputs"] != {
            name: sha256(OUTPUT / name) for name in compared
        }:
            check.unlink()
    (OUTPUT / "index.md").write_text(
        "# Current CTA selection replay\n\n[Human and dog report, figures and sources](../cta-source-walkthrough-20261028/index.md).\n\nThis directory contains the versioned full-context replay, expression-profile evidence, classifications, source receipts and gene/group identity crosswalks.\n"
    )
    receipts()
    build_cta_selection_release.build(OUTPUT)
    audit_cta_source_coverage.build()
    build_cta_selection_audit.build()
    render_cta_source_walkthrough.build()
    walkthrough_receipt()


if __name__ == "__main__":
    main()

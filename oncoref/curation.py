"""Run and verify the curation dependency graph from a source checkout."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from importlib.resources import files
from pathlib import Path

from . import catalog, reference_data, target_selection

RECEIPT = "oncoref/data/curation-build-receipt.json"
TUMOR_MANIFEST = "docs/audits/cta-selection-20261028/tumor-identity-repair/manifest.json"
#: Summary bulk/brain/immune references read by the curation build outside reference_inputs.
SUMMARY_REFERENCE_INPUTS = (
    ("hpa_rna_consensus", "v23"),
    ("rna_tissue_hpa", "v23"),
    ("rna_tissue_gtex", "v23"),
    ("rna_brain_gtex", "v23"),
    ("rna_brain_hpa", "v23"),
    ("gtex_tissue_rna", "v11"),
    ("hpa_immune_cell_samples", "v23"),
    ("hpa_immune_monaco", "v23"),
)


def manifest():
    return json.loads(files("oncoref").joinpath("data/curation-build.json").read_text())


def reference_input(role, *, import_path=None):
    """Resolve a declared source/version to checksum-verified, content-addressed bytes."""
    source = manifest()["reference_inputs"][role]
    return reference_data.pinned_path(source["source"], source["version"], import_path=import_path)


def required_reference_inputs() -> list[tuple[str, str]]:
    """Every external reference (source, version) required by the curation build."""
    spec = manifest()
    items: list[tuple[str, str]] = [
        (info["source"], info["version"]) for info in spec["reference_inputs"].values()
    ]
    items.extend(SUMMARY_REFERENCE_INPUTS)
    for name in sorted(reference_data.REFERENCE_SOURCES):
        if name.startswith("dice_"):
            items.append((name, reference_data.resolve_version(name)))
    for panel in spec.get("reference_panels", {}).values():
        source = panel["source"]
        prefix = panel["release"] + "-"
        pins = reference_data.REFERENCE_SOURCES[source].get("pins", {})
        for version in sorted(v for v in pins if v.startswith(prefix)):
            items.append((source, version))
    return list(dict.fromkeys(items))


def required_tumor_matrices(root: Path) -> list[tuple[str, str, str]]:
    """Pinned per-sample tumor cohorts as (cancer_code, version, sha256) tuples."""
    raw = json.loads((Path(root) / TUMOR_MANIFEST).read_text())
    out = []
    for entry, expected_sha256 in sorted(raw["matrices"].items()):
        path = Path(entry)
        out.append((path.stem, path.parent.name.removeprefix("v"), expected_sha256))
    return out


def fetch_inputs(root: Path | None = None, *, include_tumors: bool = False) -> dict:
    """Materialize and verify every external dataset required to rebuild the catalog."""
    fetched_refs = 0
    for source, version in required_reference_inputs():
        pins = reference_data.REFERENCE_SOURCES[source].get("pins", {})
        if version in pins and "sha256" in pins[version]:
            reference_data.pinned_path(source, version)
        else:
            reference_data.ensure(source, version)
        fetched_refs += 1
    fetched_tumors = 0
    if include_tumors:
        root = Path(__file__).resolve().parents[1] if root is None else Path(root)
        for code, version, expected_sha256 in required_tumor_matrices(root):
            path = catalog.ensure(f"per-sample:{code}", version=version)
            if digest(path) != expected_sha256:
                raise ValueError(f"Checksum mismatch for tumor matrix {code}@{version}: {path}")
            fetched_tumors += 1
    return {
        "reference_inputs_verified": fetched_refs,
        "tumor_matrices_verified": fetched_tumors,
    }


def ordered_steps(steps):
    """Topological order; missing dependencies and cycles fail before any mutation."""
    ordered, visiting = [], set()

    def visit(name):
        if name in ordered:
            return
        if name not in steps or name in visiting:
            raise ValueError(f"Missing or cyclic curation dependency: {name}")
        visiting.add(name)
        for dependency in steps[name]["needs"]:
            visit(dependency)
        visiting.remove(name)
        ordered.append(name)

    for name in steps:
        visit(name)
    return ordered


def digest(path):
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def write_receipt(root):
    """Bind outputs and the implementation used to generate them, without code copies."""
    spec = manifest()
    outputs = {root / name for name in spec["outputs"]}
    for name in spec["output_roots"]:
        outputs.update(p for p in (root / name).rglob("*") if p.is_file())
    # The report already enumerates every implementation dependency it binds.
    report = json.loads(
        (root / "docs/audits/target-selection-20261110/rebuild/manifest.json").read_text()
    )
    implementation = {
        "oncoref/curation.py",
        "oncoref/antigen_data.py",
        "oncoref/catalog.py",
        "scripts/cta_evidence_assets/antigen-whitepaper.md",
        "docs/audits/normal-cell-atlas-20261110/rare-cell/findings.csv",
        "oncoref/data/cta-trial-evidence.csv",
        "oncoref/data/curation-build.json",
        "oncoref/data/antigen-source-inventory.csv",
        "oncoref/data/cta-literature-benchmarks.csv",
        "oncoref/data/marcinkowski2019-ct83-expression.csv",
        "oncoref/target_filters.py",
        "oncoref/host_context.py",
        "scripts/cta_evidence_assets/target-selection-guide.md",
        "scripts/cta_evidence_assets/target-filter-guide.md",
        "oncoref/antigen_sets.py",
        "oncoref/antigen_cli.py",
        "oncoref/reference_data.py",
        "scripts/gtex_sample_evidence.py",
        "oncoref/data/reference-releases.json",
        "oncoref/data/normal-cell-sources.json",
        "oncoref/data/cta-literature-followup.json",
        "oncoref/tissue_tiers.py",
        "oncoref/normal_tissues.py",
        "scripts/replay_cta_selection.py",
    }
    implementation.update(step["argv"][0] for step in spec["steps"].values())
    implementation.update(report["implementation_sha256"])
    receipt = {
        "schema": "oncoref.curation-build-receipt.v1",
        "policy_id": target_selection.POLICY_ID,
        "steps": ordered_steps(spec["steps"]),
        "implementation": {name: digest(root / name) for name in sorted(implementation)},
        "outputs": {str(p.relative_to(root)): digest(p) for p in sorted(outputs)},
    }
    (root / RECEIPT).write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def check(root):
    """Verify the stored build's policy, graph, implementation and output hashes."""
    receipt = json.loads((root / RECEIPT).read_text())
    if receipt["policy_id"] != target_selection.POLICY_ID:
        raise ValueError("Curation receipt belongs to another policy")
    if receipt["steps"] != ordered_steps(manifest()["steps"]):
        raise ValueError("Curation dependency graph changed")
    errors = []
    for section in ("implementation", "outputs"):
        for name, expected in receipt[section].items():
            path = root / name
            if not path.is_file() or digest(path) != expected:
                errors.append(name)
    required = set(manifest().get("outputs", []))
    errors.extend(sorted(required - receipt["outputs"].keys()))
    for stem, field in (
        ("target-selection", "data_sha256"),
        ("rare-cell-evidence", "sha256"),
        ("host-exposure", "sha256"),
    ):
        metadata = root / f"oncoref/data/{stem}.json"
        data = root / f"oncoref/data/{stem}.csv.gz"
        if (
            str(metadata.relative_to(root)) in receipt["outputs"]
            and metadata.is_file()
            and data.is_file()
            and json.loads(metadata.read_text())[field] != digest(data)
        ):
            errors.append(str(metadata.relative_to(root)))
    report_path = root / "docs/audits/target-selection-20261110/rebuild/manifest.json"
    if report_path.is_file():
        report = json.loads(report_path.read_text())
        for name, expected in report["implementation_sha256"].items():
            if not (root / name).is_file() or digest(root / name) != expected:
                errors.append(name)
        for name, expected in report["outputs_sha256"].items():
            path = report_path.parent / name
            if not path.is_file() or digest(path) != expected:
                errors.append(str(path.relative_to(root)))
    if errors:
        raise ValueError("Stale or missing curation artifacts: " + ", ".join(sorted(set(errors))))
    return {"policy_id": receipt["policy_id"], "verified_outputs": len(receipt["outputs"])}


def category_diff(previous, current=None):
    """Compare data-only catalog snapshots by stable protein identity."""
    from .target_filters import POSITIVE_CONTROLS

    current = target_selection.results() if current is None else current
    for frame in (previous, current):
        if not frame.entity_id.is_unique:
            raise ValueError("Catalog snapshots must have unique protein identities")
    before = previous[["entity_id", "category"]].rename(columns={"category": "previous_category"})
    after = current[["entity_id", "protein_group_label", "category"]]
    result = before.merge(after, on="entity_id", how="outer", validate="one_to_one")
    result["positive_control"] = result.protein_group_label.isin(POSITIVE_CONTROLS)
    return (
        result[result.previous_category.fillna("").ne(result.category.fillna(""))]
        .sort_values(["positive_control", "protein_group_label"], ascending=[False, True])
        .reset_index(drop=True)
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["build", "check", "diff", "fetch"])
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--check", action="store_true", help="Verify the existing build without rebuilding"
    )
    parser.add_argument(
        "--reuse-tumor-cache",
        action="store_true",
        help="Use the builder's checksum-verified tumor cache",
    )
    parser.add_argument(
        "--include-tumors",
        action="store_true",
        help="For fetch: also download and verify all 132 pinned per-sample tumor matrices",
    )
    parser.add_argument(
        "--against",
        type=Path,
        help="Prior catalog CSV, optionally gzipped; no frozen Python code required",
    )
    parser.add_argument("--format", choices=["table", "csv", "json"], default="table")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    if args.command == "fetch":
        print(json.dumps(fetch_inputs(root, include_tumors=args.include_tumors)))
        return
    if args.command == "diff":
        import pandas as pd

        from .antigen_cli import emit

        if args.against is None:
            parser.error("diff requires --against PATH to a prior catalog snapshot")
        emit(category_diff(pd.read_csv(args.against)), args.format)
        return
    if args.check or args.command == "check":
        print(json.dumps(check(root)))
        return
    spec = manifest()
    order = ordered_steps(spec["steps"])
    for name in order:
        if not (root / spec["steps"][name]["argv"][0]).is_file():
            parser.error("Curation builds require a complete source checkout")
    env = os.environ | {"PYTHONPATH": str(root)}
    for name in order:
        command = [sys.executable, *spec["steps"][name]["argv"]]
        if name == "catalog" and args.reuse_tumor_cache:
            command.append("--reuse-tumor-cache")
        print(f"Curation step: {name}", flush=True)
        subprocess.run(command, cwd=root, env=env, check=True)
    write_receipt(root)
    print(json.dumps(check(root)))


if __name__ == "__main__":
    main()

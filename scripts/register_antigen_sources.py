"""Publish source selectors and checksums for the installed catalog from a committed data snapshot."""

import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd

from oncoref import curation, reference_data, source_matrices

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "oncoref/data"
BUILD = Path("docs/audits/target-selection-20261110/rebuild")
ATLAS = Path("docs/audits/normal-cell-atlas-20261110")
EVIDENCE = Path("docs/audits/antigen-evidence-review-20261118")
EXTENSIONS = (".csv.gz", ".csv", ".json.gz", ".json", ".parquet", ".yaml")


def register(revision):
    revision = subprocess.check_output(
        ["git", "rev-parse", "--verify", revision + "^{commit}"], cwd=ROOT, text=True
    ).strip()
    releases_path = DATA / "reference-releases.json"
    releases = json.loads(releases_path.read_text())
    rows = {}

    def record(name, version, role, digest, url="", location=""):
        key = (name, version)
        if key in rows and rows[key]["sha256"] != digest:
            raise ValueError(f"Conflicting input hashes for {key}")
        rows[key] = {
            "dataset": name,
            "version": version,
            "role": role,
            "sha256": digest,
            "url": url,
            "source_file": location,
        }

    def ref(name, version, role, expected=None):
        spec = reference_data.REFERENCE_SOURCES[name]
        pin = spec.get("pins", {}).get(version, {})
        actual = pin.get("sha256", expected)
        if expected and actual != expected:
            raise ValueError(f"Wrong upstream pin for {name}/{version}")
        record(name, version, role, actual, spec["urls"][version])

    def artifact(path, role, expected=None):
        path = Path(path)
        digest = curation.digest(ROOT / path)
        if expected and digest != expected:
            raise ValueError(f"Wrong source bytes for {path}")
        # Pin only bytes already public at this commit (or about to be pushed with it).
        committed = subprocess.check_output(
            ["git", "rev-parse", f"{revision}:{path.as_posix()}"], cwd=ROOT, text=True
        ).strip()
        local = subprocess.check_output(["git", "hash-object", str(ROOT / path)], text=True).strip()
        if committed != local:
            raise ValueError(f"Commit the input before registering it: {path}")
        url = f"https://raw.githubusercontent.com/pirl-unc/oncoref/{revision}/{path.as_posix()}"
        if path.parent == Path("oncoref/data"):
            record(path.name, "", role, digest, url, str(path))
            return
        stem = next(path.name.removesuffix(ext) for ext in EXTENSIONS if path.name.endswith(ext))
        name = "antigen-" + stem
        if stem == "manifest":
            name += "-" + path.parent.name
        version = revision[:12]
        spec = releases["sources"].setdefault(
            name,
            {
                "description": "Antigen curation " + role + " — " + stem.replace("-", " "),
                "filename": path.name,
                "archive_format": "file",
                "kind": "reference",
                "versions": {},
            },
        )
        if spec["filename"] != path.name:
            raise ValueError(f"Ambiguous source name: {path}")
        # Unchanged source bytes retain the same immutable download and selector.
        for existing, previous in spec["versions"].items():
            if previous["sha256"] == digest and previous["url"].endswith("/" + path.as_posix()):
                spec["default_version"] = existing
                record(name, existing, role, digest, previous["url"], str(path))
                return
        spec["default_version"] = version
        spec["versions"][version] = {
            "url": url,
            "sha256": digest,
            "archive_sha256": digest,
            "bytes": (ROOT / path).stat().st_size,
            "taxon_id": 9606,
        }
        record(name, version, role, digest, url, str(path))

    # Include supplemental native sources even when the assay contributes context only.
    for name in reference_data._normal_manifest["artifacts"]:
        ref(name, reference_data.resolve_version(name), "native assay or supplement")
    for spec in curation.manifest()["reference_inputs"].values():
        ref(spec["source"], spec["version"], "expression or sample design")
    receipts = [
        json.loads((ROOT / BUILD / "inputs.json").read_text()),
        json.loads((ROOT / ATLAS / "rare-cell/regrade-inputs.json").read_text())["inputs"],
        json.loads((DATA / "host-exposure.json").read_text())["inputs"],
    ]
    for receipt in receipts:
        for key, digest in receipt.items():
            if key.startswith("reference:"):
                _, name, version = key.split(":", 2)
                ref(name, version, "curation input", digest)
            elif not key.endswith(".py"):
                artifact(key, "curation input", digest)

    artifact(
        Path("docs/audits/cta-selection-20261028/widely-studied-contexts.csv.gz"),
        "assay observations",
    )

    matrices = json.loads(
        (
            ROOT / "docs/audits/cta-selection-20261028/tumor-identity-repair/manifest.json"
        ).read_text()
    )["matrices"]
    for path, digest in matrices.items():
        path = Path(path)
        code, version = path.stem, path.parent.name.removeprefix("v")
        record(
            "per-sample:" + code,
            version,
            "tumor expression",
            digest,
            source_matrices.release_url(code, version),
        )

    for filename in [
        "normal-comparison.csv.gz",
        "reproductive-comparison.csv.gz",
        "tumor-witnesses.csv.gz",
        "library-qc.csv.gz",
        "tumor-samples.csv.gz",
        "gtex-library-qc.csv.gz",
        "normal-library-qc.csv",
        "normal-rna/immune-library-qc.csv",
        "tumor-source-coverage.csv",
    ]:
        artifact(BUILD / filename, "comparison or sample metadata")
    for filename in [
        "source-registry.csv",
        "dataset-versions.csv",
        "rare-cell/panel.csv",
        "siletti-granule-gene-detection.csv.gz",
    ]:
        artifact(ATLAS / filename, "cell atlas")
    for filename in [
        "hpa-assay-measurements.csv.gz",
        "protein-staining.csv.gz",
        "donor-evidence.csv.gz",
    ]:
        artifact(EVIDENCE / filename, "assay observations")
    for filename in [
        "cta-literature-benchmarks.csv",
        "cta-trial-evidence.csv",
        "normal-cell-publication-evidence.csv",
        "marcinkowski2019-ct83-expression.csv",
        "cta-publication-sources.csv",
        "cta-publication-membership.csv",
        "host-exposure.csv.gz",
        "target-selection.csv.gz",
    ]:
        artifact(Path("oncoref/data") / filename, "curated data")
    releases_path.write_text(json.dumps(releases, indent=2) + "\n")
    pd.DataFrame(rows.values()).sort_values(["role", "dataset", "version"]).to_csv(
        DATA / "antigen-source-inventory.csv", index=False
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("revision", help="Commit containing the exact input and comparison files")
    register(parser.parse_args().revision)

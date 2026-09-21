"""Shared, testable provenance and background contracts for CTA reports."""

from __future__ import annotations

import gzip
import hashlib
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ID = "Ensembl_Gene_ID"
ROOT = Path(__file__).resolve().parents[1]
MIN_RANKED_PATIENTS = 10
RANKED_COHORT_POLICY = f"min{MIN_RANKED_PATIENTS}"


def ranked_selection_dir(out, prevalence, percentile):
    return (
        Path(out)
        / "selections"
        / RANKED_COHORT_POLICY
        / f"prevalence_gt{prevalence}_transcriptome_p{percentile}"
    )


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def implementation_hash():
    paths = sorted((ROOT / "scripts").glob("*cta*.py"))
    paths += sorted((ROOT / "oncoref").glob("*.py"))
    paths += sorted((ROOT / "oncoref/data").glob("*.csv"))
    paths += sorted((ROOT / "oncoref/data").glob("*.yaml"))
    return fingerprint({str(p.relative_to(ROOT)): sha256(p) for p in paths})


def common_background(cohorts):
    """Canonical biological loci present in every declared source matrix.

    Identity intersection is determined before looking at candidate expression.
    The analysis refuses incomplete patient backgrounds rather than silently
    taking a quantile over a different set of genes.
    """
    from oncoref import source_matrices
    from oncoref.gene_families import clean_tpm_censored_gene_ids
    from oncoref.gene_ids import resolve_ensembl_id

    shared = None
    sources = {}
    for code in sorted(cohorts):
        path = source_matrices.local_path(code)
        frame = pd.read_parquet(path, columns=[ID])
        ids = set(frame[ID].astype(str).map(resolve_ensembl_id))
        ids -= clean_tpm_censored_gene_ids()
        shared = ids if shared is None else shared & ids
        sources[code] = sha256(path)
    if not shared:
        raise ValueError("No common biological gene background")
    return sorted(shared), sources


def background_values(frame, identifiers, columns, id_column=ID):
    selected = frame.set_index(id_column).reindex(identifiers)[columns].to_numpy(dtype=float)
    if not np.isfinite(selected).all() or (selected < 0).any():
        raise ValueError("Common background contains missing or invalid measurements")
    return selected


def checkpoint_valid(prefix, expected, suffixes):
    prefix = Path(prefix)
    metadata = prefix.with_suffix(".json")
    if not metadata.exists():
        return False
    try:
        stored = json.loads(metadata.read_text())
        return (
            stored.get("checkpoint_fingerprint") == expected
            and all(prefix.with_name(prefix.name + suffix).is_file() for suffix in suffixes)
            and all(
                sha256(prefix.with_name(prefix.name + suffix)) == digest
                for suffix, digest in stored.get("checkpoint_outputs", {}).items()
            )
            and set(stored.get("checkpoint_outputs", {})) == set(suffixes)
        )
    except (ValueError, OSError):
        return False


def checkpoint_payload(prefix):
    """Return scientific metadata without cache bookkeeping in exported tables."""
    stored = json.loads(Path(prefix).with_suffix(".json").read_text())
    return {key: value for key, value in stored.items() if not key.startswith("checkpoint_")}


def write_checkpoint(prefix, metadata, expected, suffixes):
    prefix = Path(prefix)
    metadata = dict(
        metadata,
        checkpoint_fingerprint=expected,
        checkpoint_outputs={s: sha256(prefix.with_name(prefix.name + s)) for s in suffixes},
    )
    prefix.with_suffix(".json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")


def write_csv(frame, path, **kwargs):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        if path.suffix == ".gz":
            # gzip otherwise derives its header filename even from a file object.
            with gzip.GzipFile(fileobj=handle, mode="wb", filename="", mtime=0) as compressed:
                frame.to_csv(compressed, index=False, **kwargs)
        else:
            frame.to_csv(handle, index=False, **kwargs)


def seal_stage(out, name, inputs, outputs):
    out = Path(out)
    receipt = {
        "inputs": {str(Path(p).resolve()): sha256(p) for p in inputs},
        "outputs": {str(Path(p).resolve()): sha256(p) for p in outputs},
    }
    (out / f"{name}_receipt.json").write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")


def verify_stage(out, name):
    receipt = Path(out) / f"{name}_receipt.json"
    if not receipt.is_file():
        raise ValueError(f"Missing {name} provenance receipt; rebuild the report")
    for paths in json.loads(receipt.read_text()).values():
        for filename, digest in paths.items():
            if not Path(filename).is_file() or sha256(filename) != digest:
                raise ValueError(f"Stale {name} input/output: {filename}; rebuild the report")


def verify_analysis(out):
    out = Path(out)
    if json.loads((out / "validation.json").read_text()).get("status") != "passed":
        raise ValueError("Analysis validation did not pass")
    verify_stage(out, "analysis")


def deterministic_zip(out, destination):
    out, destination = Path(out), Path(destination)
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(out.rglob("*")):
            if not path.is_file() or path == destination or path.suffix in {".zip", ".log"}:
                continue
            info = zipfile.ZipInfo(str(path.relative_to(out)), date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            archive.writestr(info, path.read_bytes())

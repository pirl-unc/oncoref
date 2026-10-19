#!/usr/bin/env python3
"""Build genome-wide HPA cancer RNA summaries from checksum-pinned sample rows.

Acquire the raw files listed in oncoref/data/hpa-cancer-sources.json (upstream
URLs or the archived copies), then run with --sources DIR --output DIR.
The builder never fetches mutable upstream data or filters to CTA genes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import zipfile
from pathlib import Path

import pandas as pd

from oncoref.hpa_cancer import (
    _summarize_ihc,
    _summarize_rna_rows,
    hpa_cancer_crosswalk,
    hpa_cancer_rna_cohorts,
    hpa_cancer_sources,
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_cohort_metadata(sources: Path):
    """Require the shipped RNA vocabulary/counts to match the pinned raw metadata."""
    frames = []
    for cohort, name in [("TCGA", "tcga"), ("validation", "validation")]:
        frame = pd.read_csv(sources / f"cancer_rna_{name}_cancers.tsv.zip", sep="\t")
        frame = frame.rename(
            columns={
                "Cancer": "cancer_type",
                "Abbriviation": "cancer_code",
                "Cancer group": "cancer",
                "Organ": "organ",
                "Num samples": "nominal_samples",
            }
        )
        for col in ["cancer_type", "cancer_code", "cancer", "organ"]:
            frame[col] = frame[col].str.strip()
        frame["cancer"] = frame["cancer"].str.lower()
        frame["cohort"] = cohort
        frames.append(frame)
    actual = pd.concat(frames, ignore_index=True)
    expected = hpa_cancer_rna_cohorts()
    keys = ["cohort", "cancer_code"]
    pd.testing.assert_frame_equal(
        actual.set_index(keys).sort_index(),
        expected.set_index(keys).sort_index(),
        check_dtype=False,
    )
    return expected


def complete_gene_chunks(path, *, chunksize=500_000):
    """Yield complete contiguous gene blocks and reject a repeated gene block.

    Carry the trailing gene across parser chunks so duplicate sample identities
    straddling a chunk boundary are checked together. The pinned HPA source is
    gene-contiguous; an incompatible future ordering fails rather than creating
    duplicate or incomplete summary rows.
    """
    seen = set()
    carry = None
    reader = pd.read_csv(path, sep="\t", chunksize=chunksize)
    for chunk in reader:
        if carry is not None:
            chunk = pd.concat([carry, chunk], ignore_index=True)
        if chunk["Gene"].isna().any():
            raise ValueError("missing RNA gene identity")
        last = chunk["Gene"].iloc[-1]
        # Only the trailing contiguous block is held back (not earlier copies).
        starts = chunk.index[chunk["Gene"].ne(chunk["Gene"].shift())]
        boundary = starts[-1]
        carry = chunk.iloc[boundary:].copy()
        batch = chunk.iloc[:boundary]
        if not batch.empty:
            genes = batch.loc[batch["Gene"].ne(batch["Gene"].shift()), "Gene"]
            if genes.duplicated().any() or seen.intersection(genes):
                raise ValueError("RNA source is not gene-contiguous")
            seen.update(genes)
            yield batch
        if last in seen:
            raise ValueError("RNA source is not gene-contiguous")
    if carry is not None and not carry.empty:
        if carry["Gene"].iloc[0] in seen:
            raise ValueError("RNA source is not gene-contiguous")
        yield carry


def build(sources: Path, output: Path, *, chunksize=500_000):
    manifest = hpa_cancer_sources()
    output.mkdir(parents=True, exist_ok=True)
    for name, spec in manifest["raw_sources"].items():
        path = sources / name
        if path.stat().st_size != spec["bytes"] or sha256(path) != spec["sha256"]:
            raise ValueError(f"source checksum mismatch: {name}")
    ihc = _summarize_ihc(pd.read_csv(sources / "cancer_data.tsv.zip", sep="\t"))
    if set(ihc["cancer"]) != set(hpa_cancer_crosswalk()["cancer"]):
        raise ValueError("crosswalk does not cover the source IHC groups exactly")
    cohorts = verify_cohort_metadata(sources)
    symbols = ihc[["gene_id", "hpa_gene_name"]].drop_duplicates()
    if symbols["gene_id"].duplicated().any():
        raise ValueError("inconsistent HPA gene names")
    path = output / "hpa_cancer_rna_summary.tsv"
    tmp = path.with_suffix(".tsv.part")
    rows = source_rows = genes = 0
    cohort_rows = {}
    missing_measurements = 0
    with tmp.open("w") as stream:
        for index, batch in enumerate(
            complete_gene_chunks(sources / "rna_cancer_sample.tsv.gz", chunksize=chunksize)
        ):
            summary = _summarize_rna_rows(batch, cohorts)
            if summary.duplicated(["gene_id", "cancer_code", "cohort"]).any():
                raise ValueError("duplicate normalized RNA gene/type/cohort")
            summary = summary.merge(symbols, on="gene_id", how="left", validate="many_to_one")
            summary.to_csv(stream, sep="\t", index=False, header=(index == 0), float_format="%.12g")
            rows += len(summary)
            source_rows += len(batch)
            genes += batch["Gene"].nunique()
            missing_measurements += int(summary["missing_samples"].sum())
            for cohort, count in summary["cohort"].value_counts().items():
                cohort_rows[cohort] = cohort_rows.get(cohort, 0) + int(count)
            if index % 20 == 0:
                print(
                    f"{source_rows:,} sample rows; {genes:,} genes; {rows:,} summary rows",
                    flush=True,
                )
    tmp.replace(path)
    archive = path.with_suffix(".tsv.zip")
    # Stable ZIP metadata makes repeat builds byte-reproducible.
    info = zipfile.ZipInfo(path.name, date_time=(1980, 1, 1, 0, 0, 0))
    info.compress_type = zipfile.ZIP_DEFLATED
    with zipfile.ZipFile(archive, "w") as zf, path.open("rb") as src, zf.open(info, "w") as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024)
    artifact = {
        "description": "HPA 25.1 cancer RNA native-pTPM prevalence (genome-wide; TCGA/validation separate)",
        "filename": path.name,
        "url": f"https://github.com/pirl-unc/oncoref/releases/download/hpa-cancer-v25.1-1/{archive.name}",
        "archive_sha256": sha256(archive),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }
    audit = {
        "scope": "genome-wide",
        "hpa_release": manifest["hpa_release"],
        "rna_sample_rows": source_rows,
        "rna_genes": genes,
        "rna_summary_rows": rows,
        "rna_rows_by_cohort": cohort_rows,
        "explicit_missing_rna_measurements": missing_measurements,
        "ihc_rows": len(ihc),
        "ihc_genes": ihc["gene_id"].nunique(),
        "ihc_groups": ihc["cancer"].nunique(),
        "ihc_status_counts": {
            k: int(v) for k, v in ihc["measurement_status"].value_counts().items()
        },
        "artifact": artifact,
    }
    (output / "build-audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2), flush=True)
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chunksize", type=int, default=500_000)
    args = parser.parse_args()
    build(args.sources, args.output, chunksize=args.chunksize)


if __name__ == "__main__":
    main()

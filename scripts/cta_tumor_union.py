#!/usr/bin/env python3
"""Audit maximum observed tumor RNA across all compatible OncoRef matrices.

Keep native linear TPM separate from HPA pTPM. No clean-TPM, microarray,
CPM-proxy, medoid or cohort-median substitutes for an observed sample.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import sys
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_report_common import input_hashes  # noqa: E402
from cta_report_common import sha256 as digest  # noqa: E402
from cta_report_common import write_csv as write  # noqa: E402

from oncoref import data_bundle, expression, source_matrices  # noqa: E402
from oncoref.gene_matching import (  # noqa: E402
    attach_gene_identities,
    harmonize_reference_expression,
)
from oncoref.protein_equivalence import equivalence_groups  # noqa: E402


def summarize_matrix(frame: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    """One actual winning sample per gene, with missingness retained."""
    values = frame[samples].to_numpy(dtype=float)
    if np.isinf(values).any() or (values < 0).any():
        raise ValueError("Expected nonnegative finite linear TPM or missing values")
    measured = np.isfinite(values).sum(axis=1)
    best = np.where(np.isfinite(values), values, -np.inf).argmax(axis=1)
    maximum = values[np.arange(len(values)), best]
    return pd.DataFrame(
        {
            "gene_id": frame.Ensembl_Gene_ID.to_numpy(),
            "maximum_tpm": np.where(measured > 0, maximum, np.nan),
            "maximum_sample": np.where(measured > 0, np.asarray(samples)[best], ""),
            "measured_samples": measured,
            "samples_ge2": (values >= 2).sum(axis=1),
        }
    )


def winning_observations(paths: Iterable[Path], *, entity_column: str) -> pd.DataFrame:
    """Select actual winning rows with stable ties and bounded merge memory.

    Each cohort retains its complete identity bundle. Only current winners are
    carried between merges, avoiding an in-memory copy of every cohort's
    repeated identity annotations. Measured values take precedence over missing
    values; ties use cancer code and then sample identity.
    """
    best = pd.DataFrame()
    for path in paths:
        part = pd.read_csv(path)
        best = (
            pd.concat([best, part], ignore_index=True)
            .sort_values(
                ["maximum_tpm", "cancer_code", "maximum_sample"], ascending=[False, True, True]
            )
            .drop_duplicates(entity_column)
        )
    return best.reset_index(drop=True)


def build(out: Path, download_cache: Path, *, cancer_codes: set[str] | None = None) -> None:
    """Recompute per-sample maxima and auditable identities across compatible cohorts."""
    out.mkdir(parents=True, exist_ok=True)
    shards = out / "cohorts"
    shards.mkdir(exist_ok=True)
    identities = out / "identity"
    identities.mkdir(exist_ok=True)
    group_shards = out / "protein-groups"
    group_shards.mkdir(exist_ok=True)
    universe = pd.read_csv(
        ROOT / "docs/audits/cta-geometric-thresholds-20261027/gene-comparison.csv.gz",
        usecols=["gene_id"],
    )
    groups = equivalence_groups(universe.gene_id)
    groups = groups[groups.member_count.gt(1)]
    members = {
        gid: r.protein_equivalence_group_id
        for r in groups.itertuples()
        for gid in json.loads(r.member_gene_ids)
    }
    registry = source_matrices.registry()
    if cancer_codes is not None:
        unknown = cancer_codes - set(registry.cancer_code)
        if unknown:
            raise ValueError(f"Unknown cancer codes: {sorted(unknown)}")
        registry = registry[registry.cancer_code.isin(cancer_codes)].copy()
    availability_path = ROOT / "oncoref/data/cancer-reference-expression-availability.csv"
    availability = pd.read_csv(availability_path).drop_duplicates(["cancer_code", "source_cohort"])
    inventory = registry.merge(
        availability, on=["cancer_code", "source_cohort"], how="left", validate="one_to_one"
    )
    inventory["eligible_scale"] = inventory.linear_tpm_comparable.eq(True)
    existing = {code: source_matrices.local_path(code) for code in registry.cancer_code}
    os.environ[source_matrices.CACHE_DIR_ENV_VAR] = str(download_cache)
    eligible = inventory[inventory.eligible_scale]

    def acquire(code):
        path = existing[code]
        if not path.exists():
            path = source_matrices.ensure(code)
        return code, path

    paths = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(acquire, code) for code in eligible.cancer_code]
        for future in concurrent.futures.as_completed(futures):
            code, path = future.result()
            paths[code] = path
            print(f"Available {len(paths)}/{len(eligible)}: {code}", flush=True)
    qc_path = data_bundle.cache_dir() / "source-matrix-sample-qc.csv"
    qc = pd.read_csv(qc_path, low_memory=False)
    records, samples, hashes = [], [], {}
    for r in inventory.itertuples():
        record = {
            "cancer_code": r.cancer_code,
            "source_cohort": r.source_cohort,
            "source_scale_class": r.source_scale_class,
            "source_matrix_version": r.source_matrix_version,
            "eligible_scale": r.eligible_scale,
            "registered_samples": r.n_samples,
        }
        if not r.eligible_scale:
            record.update(status="Excluded: incompatible scale", samples_used=0)
            records.append(record)
            continue
        path = paths[r.cancer_code]
        hashes[str(path)] = digest(path)
        q = qc[qc.cancer_code.eq(r.cancer_code) & qc.source_cohort.eq(r.source_cohort)]
        if q.empty or q.sample_id.duplicated().any():
            raise ValueError(f"Missing/ambiguous sample QC for {r.cancer_code}")
        use = q[q.sample_qc_status.eq("pass") & q.recommended_for_absolute_tpm_floor.eq(True)]
        raw = pd.read_parquet(path)
        all_samples = expression.sample_columns(raw)
        if set(all_samples) != set(q.sample_id.astype(str)):
            raise ValueError(f"QC does not describe exact matrix samples: {r.cancer_code}")
        sample_ids = sorted(set(use.sample_id.astype(str)))
        record.update(
            status="Used" if sample_ids else "No passing sample QC",
            samples_used=len(sample_ids),
            matrix_sha256=hashes[str(path)],
        )
        records.append(record)
        if not sample_ids:
            continue
        matrix, identity = harmonize_reference_expression(
            raw[[*expression.id_columns(raw), *sample_ids]],
            sample_cols=sample_ids,
            source=r.source_cohort,
            source_version=str(r.source_matrix_version),
            symbol_column="Symbol" if "Symbol" in raw else None,
            duplicate_policy="quarantine",
        )
        write(identity, identities / f"{r.cancer_code}.csv.gz")
        result = summarize_matrix(matrix, sample_ids)
        result["cancer_code"], result["source_cohort"] = r.cancer_code, r.source_cohort
        result["unit"] = "TPM"
        result = attach_gene_identities(result)
        part = matrix[matrix.Ensembl_Gene_ID.isin(members)].copy()
        part["group_id"] = part.Ensembl_Gene_ID.map(members)
        summed = part.groupby("group_id")[sample_ids].sum(min_count=1)
        measured = part.groupby("group_id")[sample_ids].count()
        group_result = summarize_matrix(
            summed.reset_index().rename(columns={"group_id": "Ensembl_Gene_ID"}), sample_ids
        )
        group_result = group_result.rename(columns={"gene_id": "protein_equivalence_group_id"})
        group_result["measured_members_at_maximum"] = [
            int(measured.loc[r.protein_equivalence_group_id, r.maximum_sample])
            if r.maximum_sample
            else 0
            for r in group_result.itertuples()
        ]
        group_result["cancer_code"], group_result["source_cohort"] = r.cancer_code, r.source_cohort
        group_result["unit"] = "TPM"
        write(
            group_result.merge(groups, on="protein_equivalence_group_id", validate="one_to_one"),
            group_shards / f"{r.cancer_code}.csv.gz",
        )
        result.to_csv(
            shards / f"{r.cancer_code}.csv.gz",
            index=False,
            compression={"method": "gzip", "mtime": 0},
        )
        namespace = source_matrices.source_sample_namespace(r.source_cohort)
        samples.extend(
            {
                "cancer_code": r.cancer_code,
                "source_cohort": r.source_cohort,
                "namespace": namespace,
                "sample_id": s,
            }
            for s in sample_ids
        )
        print(f"Summarized {r.cancer_code}: {len(sample_ids)} samples", flush=True)
        del raw, matrix, result
    counts = pd.DataFrame(records)
    counts.to_csv(out / "cohort-availability.csv", index=False)
    sample_table = pd.DataFrame(samples)
    sample_table.to_csv(
        out / "sample-identity.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )
    used_codes = counts.loc[counts.samples_used.gt(0), "cancer_code"]
    best = winning_observations(
        (shards / f"{code}.csv.gz" for code in used_codes), entity_column="gene_id"
    )
    best.to_csv(
        out / "tumor-maxima.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )
    group_best = winning_observations(
        (group_shards / f"{code}.csv.gz" for code in used_codes),
        entity_column="protein_equivalence_group_id",
    )
    write(group_best, out / "protein-group-tumor-maxima.csv.gz")
    summary = {
        "registered_cohorts": len(registry),
        "compatible_cohorts": len(eligible),
        "used_cohorts": int(counts.samples_used.gt(0).sum()),
        "sample_occurrences": len(sample_table),
        "namespace_deduplicated_samples": len(
            sample_table.drop_duplicates(["namespace", "sample_id"])
        ),
        "source_namespaces": sample_table.namespace.nunique(),
        "note": "Namespaces deduplicate shared source views, not patients across independent studies or HPA.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    inputs = input_hashes(
        [
            availability_path,
            qc_path,
            ROOT / "oncoref/data/source-matrices.csv",
            ROOT / "oncoref/gene_matching.py",
            ROOT / "oncoref/protein_equivalence.py",
            ROOT / "oncoref/gene_ids.py",
            ROOT / "oncoref/data/canonical-gene-space.csv.gz",
            ROOT / "oncoref/data/ensembl-id-aliases.csv",
            ROOT / "oncoref/data/ncbi-symbol-synonyms.csv.gz",
            ROOT / "oncoref/data/proteoform-groups-genome.csv",
            Path(__file__),
        ]
    )
    (out / "manifest.json").write_text(
        json.dumps(
            {
                "inputs": inputs,
                "matrices": hashes,
                "outputs": {
                    str(p.relative_to(out)): digest(p)
                    for p in out.rglob("*")
                    if p.is_file() and p.name != "manifest.json"
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--download-cache", type=Path, default=Path("/tmp/oncoref-cta-expanded-matrices")
    )
    parser.add_argument("--cancer-codes", nargs="+", help="Rebuild specified cohort views only")
    args = parser.parse_args()
    build(
        args.output,
        args.download_cache,
        cancer_codes=set(args.cancer_codes) if args.cancer_codes else None,
    )

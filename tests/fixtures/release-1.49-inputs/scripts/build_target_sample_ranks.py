#!/usr/bin/env python3
"""Rank candidates against the complete measured protein inventory per tumor."""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import target_selection  # noqa: E402
from oncoref.gene_matching import harmonize_reference_expression  # noqa: E402

BASE = ROOT / "docs/audits/cta-selection-20261028/tumor-identity-repair"
REBUILD = ROOT / "docs/audits/target-selection-20261109/rebuild"
CACHE = ROOT / "tmp/target-threshold-cache"


def sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def build():
    CACHE.mkdir(parents=True, exist_ok=True)
    catalog = target_selection.results()
    gene_to_entity = {
        g: r.entity_id for r in catalog.itertuples() for g in json.loads(r.member_gene_ids)
    }
    expected = catalog.set_index("entity_id").member_gene_ids.map(lambda v: len(json.loads(v)))
    matrix_path = ROOT / "tmp/target-catalog-cache/tumor-comparison.parquet"
    assert (
        sha(matrix_path) == target_selection.metadata()["comparison"]["tumor_matrix_cache_sha256"]
    )
    matrix = pd.read_parquet(matrix_path)
    samples = pd.read_csv(REBUILD / "tumor-samples.csv.gz", keep_default_na=False).set_index(
        "sample_key"
    )
    paths = {
        Path(p).stem: (Path(p), digest)
        for p, digest in json.loads((BASE / "manifest.json").read_text())["matrices"].items()
    }
    original_samples = pd.read_csv(BASE / "sample-identity.csv.gz")
    blocks, cutoffs, used = [], [], {}
    for code, meta in samples.groupby("cancer_code", sort=True):
        path, digest = paths[code]
        if sha(path) != digest:
            raise ValueError(f"Changed source: {path}")
        used[str(path)] = digest
        columns = (
            original_samples.loc[original_samples.cancer_code.eq(code), "sample_id"]
            .astype(str)
            .tolist()
        )
        raw = pq.read_table(path, columns=["Ensembl_Gene_ID", "Symbol", *columns]).to_pandas()
        raw, _ = harmonize_reference_expression(
            raw,
            sample_cols=columns,
            symbol_column="Symbol",
            source="within-sample-rank-audit",
            source_version="5.22.10",
            duplicate_policy="quarantine",
        )
        raw = raw[raw.Ensembl_Gene_ID.isin(gene_to_entity)]
        d = raw[meta.sample_id.tolist()].copy()
        if np.isinf(d.to_numpy()).any() or d.lt(0).any().any():
            raise ValueError(f"Invalid full-inventory values in {code}")
        d.index = raw.Ensembl_Gene_ID.map(gene_to_entity)
        values = d.groupby(level=0).sum(min_count=1)
        counts = d.groupby(level=0).count()
        values = values.where(counts.eq(expected.reindex(values.index), axis=0))
        values.columns = meta.index
        # The admitted subset must reproduce the independently built tumor panel.
        pd.testing.assert_frame_equal(
            values.reindex(matrix.columns).T,
            matrix.loc[meta.index],
            check_names=False,
            check_dtype=False,
        )
        ranks = values.rank(axis=0, method="average", pct=True)
        blocks.append(ranks.reindex(matrix.columns).T)
        q = values.quantile([0.8, 0.9, 0.95], axis=0, interpolation="lower").T
        q.columns = ["p80_tpm", "p90_tpm", "p95_tpm"]
        q["measured_proteins"] = values.notna().sum()
        q["measured_fraction"] = q.measured_proteins / len(catalog)
        cutoffs.append(q)
        print(
            f"{code}: {len(meta)} tumors; median {q.measured_proteins.median():g} measured proteins",
            flush=True,
        )
    ranks = pd.concat(blocks).reindex(matrix.index)
    assert ranks.index.is_unique and ranks.notna().equals(matrix.notna())
    ranks.to_parquet(CACHE / "sample-ranks.parquet", compression="zstd")
    pd.concat(cutoffs).join(samples).to_csv(
        CACHE / "sample-rank-cutoffs.csv.gz", compression={"method": "gzip", "mtime": 0}
    )
    receipt = {
        "rank_universe": "All measured current-catalog protein records, including Unsafe and Other; members summed only within one sample; missing members leave the record unmeasured.",
        "rank_ties": "average",
        "candidate_matrix_sha256": sha(matrix_path),
        "catalog_sha256": sha(ROOT / "oncoref/data/target-selection.csv.gz"),
        "samples_sha256": sha(REBUILD / "tumor-samples.csv.gz"),
        "source_matrices": used,
        "implementation_sha256": {
            p: sha(ROOT / p)
            for p in ["scripts/build_target_sample_ranks.py", "oncoref/gene_matching.py"]
        },
        "outputs_sha256": {
            p: sha(CACHE / p) for p in ["sample-ranks.parquet", "sample-rank-cutoffs.csv.gz"]
        },
    }
    (CACHE / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    build()

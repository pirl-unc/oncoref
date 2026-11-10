#!/usr/bin/env python3
"""Stream donor-resolved counts for proteins needing normal-tissue reassessment.

Uses a pinned inventory; keeps all ages and assays separate. Results are evidence
for review, not automatic Unsafe decisions. No imputed/integrated expression is
used. The separate >=2 count column is a sensitivity analysis, not a new gate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import cellxgene_census as cc
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import tiledbsoma as soma
from inventory_normal_census import FILTER, GROUPS

ROOT = Path(__file__).resolve().parents[1]
PROFILES = ROOT / "docs/audits/target-selection-20261109/category-profiles"


def age_scope(stage: str) -> str:
    """Conservative age parsing. Unresolved ranges never become adult controls."""
    if any(
        x in stage
        for x in [
            "post-fertilization",
            "Carnegie",
            "LMP month",
            "embryonic",
            "organogenesis",
            "prenatal",
        ]
    ):
        return "prenatal"
    match = re.fullmatch(r"(\d+)-year-old stage", stage)
    if match:
        return "adult" if int(match[1]) >= 18 else "child"
    if stage in {
        "late adult stage",
        "middle aged stage",
        "60-79 year-old stage",
        "80 year-old and over stage",
        "90 year-old and over stage",
    }:
        return "adult"
    if stage in {
        f"{d} decade stage"
        for d in ["third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth"]
    }:
        return "adult"
    if re.fullmatch(r"\d+-month-old stage", stage) or stage in {
        "newborn stage (0-28 days)",
        "juvenile stage (5-14 yo)",
        "child stage (1-4 yo)",
        "infant stage",
    }:
        return "child"
    return "unresolved"


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def make_groups(inventory, output):
    groups = pd.read_csv(inventory / "populations.csv.gz")
    groups["population_id"] = np.arange(len(groups), dtype=np.int32)
    # HsapDv adult/prime/young-adult terms start at 15, not 18. They remain
    # unresolved for an >=18 analysis unless the original donor age is recovered.
    groups["age_scope"] = groups.development_stage.map(age_scope)
    datasets = pd.read_csv(inventory / "datasets.csv")
    groups = groups.merge(
        datasets[
            [
                "dataset_id",
                "collection_id",
                "collection_name",
                "collection_doi",
                "dataset_version_id",
            ]
        ],
        on="dataset_id",
        validate="many_to_one",
    )
    groups["donor_identity_known"] = groups.donor_id.notna() & ~groups.donor_id.str.lower().isin(
        ["unknown", "na", "nan", "pooled", ""]
    )
    groups["count_unit"] = np.where(
        groups.assay.str.startswith("10x ") & ~groups.assay.str.contains("flex"),
        "UMI_protocol",
        "assay_review",
    )
    groups["source_qc_review"] = "pending"
    groups.to_parquet(output / "populations.parquet", index=False)
    return groups


def make_row_map(inventory, groups, output):
    """Bounded-memory metadata join, independent of file/row ordering."""
    files = sorted(inventory.glob("observations-*.parquet"))
    max_id = 0
    for path in files:
        pf = pq.ParquetFile(path)
        for batch in pf.iter_batches(columns=["soma_joinid"], batch_size=2_000_000):
            max_id = max(max_id, int(batch.column(0).to_numpy().max()))
    mapping = np.full(max_id + 1, -1, dtype=np.int32)
    group_table = pa.Table.from_pandas(groups[[*GROUPS, "population_id"]], preserve_index=False)
    seen = 0
    for path in files:
        for batch in pq.ParquetFile(path).iter_batches(
            columns=["soma_joinid", *GROUPS], batch_size=1_000_000
        ):
            table = pa.Table.from_batches([batch])
            # Dictionary encodings differ across files, so join actual labels.
            table = table.cast(
                pa.schema(
                    [
                        (name, pa.string() if name != "soma_joinid" else pa.int64())
                        for name in table.column_names
                    ]
                )
            )
            joined = table.join(group_table, keys=GROUPS, join_type="left outer")
            if joined["population_id"].null_count or len(joined) != len(table):
                raise ValueError("Unmapped or duplicated inventory population")
            ids = joined["soma_joinid"].to_numpy()
            if (mapping[ids] >= 0).any() or len(np.unique(ids)) != len(ids):
                raise ValueError("Duplicated Census observation identity")
            mapping[ids] = joined["population_id"].to_numpy()
            seen += len(ids)
    if seen != int(groups.n_cells.sum()):
        raise ValueError("Population denominator mismatch")
    np.save(output / "observation-population-map.npy", mapping)
    print(f"Mapped {seen:,} observations to {len(groups):,} populations", flush=True)
    return mapping


def summarize_block(census, features, mapping, groups):
    n_genes = len(features)
    shape = (len(groups), n_genes)
    n_positive = np.zeros(shape, np.int32)
    n_ge2 = np.zeros(shape, np.int32)
    count_sum = np.zeros(shape, np.float64)
    count_max = np.zeros(shape, np.float32)
    n_fractional = np.zeros(shape, np.int32)
    var_ids = features.soma_joinid.to_numpy()
    var_map = np.full(int(var_ids.max()) + 1, -1, np.int32)
    var_map[var_ids] = np.arange(n_genes)
    with census["census_data"]["homo_sapiens"].axis_query(
        "RNA",
        obs_query=soma.AxisQuery(value_filter=FILTER),
        var_query=soma.AxisQuery(coords=(var_ids,)),
    ) as query:
        for table in query.X("raw").tables():
            rows = table["soma_dim_0"].to_numpy()
            cols = var_map[table["soma_dim_1"].to_numpy()]
            values = table["soma_data"].to_numpy()
            if np.any(rows >= len(mapping)) or np.any(mapping[rows] < 0) or np.any(cols < 0):
                raise ValueError("Expression does not match pinned inventory")
            if not np.all(np.isfinite(values) & (values > 0)):
                raise ValueError("Invalid sparse counts")
            idx = (mapping[rows], cols)
            np.add.at(n_positive, idx, 1)
            np.add.at(n_ge2, idx, values >= 2)
            np.add.at(n_fractional, idx, values != np.floor(values))
            np.add.at(count_sum, idx, values)
            np.maximum.at(count_max, idx, values)
            print(f"  streamed {len(table):,} nonzero counts", flush=True)
    row, col = np.nonzero(n_positive)
    out = pd.DataFrame(
        {
            "population_id": row,
            "gene_id": features.feature_id.to_numpy()[col],
            "n_positive": n_positive[row, col],
            "n_ge2": n_ge2[row, col],
            "count_sum": count_sum[row, col],
            "count_max": count_max[row, col],
            "n_fractional": n_fractional[row, col],
        }
    )
    if (out.n_positive.to_numpy() > groups.n_cells.to_numpy()[row]).any():
        raise ValueError("Positive count exceeds population denominator")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--block-size", type=int, default=100)
    parser.add_argument("--max-blocks", type=int)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    inventory_manifest = json.loads((args.inventory / "manifest.json").read_text())
    if inventory_manifest["filter"] != FILTER:
        raise ValueError("Inventory filter differs")
    profiles = pd.read_csv(PROFILES / "all-protein-profiles.csv.gz")
    familiar = set(pd.read_csv(PROFILES / "familiar-cta-profiles.csv").entity_id)
    selected = profiles[
        profiles.category.ne("Unsafe")
        | profiles.Unsafe_reasons.eq("Rare vital-tissue cell expression")
        | profiles.entity_id.isin(familiar)
    ].copy()
    genes = {g for ids in selected.member_gene_ids for g in json.loads(ids)}
    features = pd.read_parquet(args.inventory / "features.parquet")
    features = features[features.feature_id.isin(genes)].copy()
    # Examine the quoted genes first; all selected genes still receive identical logic.
    focus = {"PLAC1", "CBY3", "DMRT1", "CAGE1", "TMEM269", "NUTM1", "SYCP1", "MAGEB10", "PRAME"}
    features["priority"] = ~features.feature_name.isin(focus)
    features = features.sort_values(["priority", "feature_id"]).reset_index(drop=True)
    selected.to_csv(
        args.output / "selected-proteins.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    features.to_csv(args.output / "selected-features.csv", index=False)
    manifest = {
        "census_release": inventory_manifest["census_release"],
        "filter": FILTER,
        "inventory_sha256": digest(args.inventory / "manifest.json"),
        "profiles_sha256": digest(PROFILES / "all-protein-profiles.csv.gz"),
        "n_proteins": len(selected),
        "requested_genes": sorted(genes),
        "unavailable_genes": sorted(genes - set(features.feature_id)),
        "block_size": args.block_size,
        "classification_changed": False,
        "method_version": 1,
    }
    manifest_path = args.output / "manifest.json"
    if manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        raise ValueError("Resume identity mismatch; use a new output directory")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    groups = make_groups(args.inventory, args.output)
    map_path = args.output / "observation-population-map.npy"
    mapping = (
        np.load(map_path, mmap_mode="r")
        if map_path.exists()
        else make_row_map(args.inventory, groups, args.output)
    )
    ctx = cc.get_default_soma_context(
        {"soma.init_buffer_bytes": 32 * 1024**2, "py.init_buffer_bytes": 32 * 1024**2}
    )
    with cc.open_soma(census_version=manifest["census_release"], context=ctx) as census:
        if not (args.output / "dataset-feature-presence.parquet").exists():
            presence = cc.get_presence_matrix(census, organism="Homo sapiens")
            ds = pd.read_csv(args.inventory / "datasets.csv")
            ds = ds[ds.dataset_id.isin(groups.dataset_id)]
            sub = presence[ds.soma_joinid.to_numpy(), :][:, features.soma_joinid.to_numpy()].tocoo()
            pd.DataFrame(
                {
                    "dataset_id": ds.dataset_id.to_numpy()[sub.row],
                    "gene_id": features.feature_id.to_numpy()[sub.col],
                }
            ).to_parquet(args.output / "dataset-feature-presence.parquet", index=False)
        for start in range(0, len(features), args.block_size):
            block_id = start // args.block_size
            if args.max_blocks is not None and block_id >= args.max_blocks:
                break
            path = args.output / f"counts-{block_id:04d}.parquet"
            if path.exists():
                continue
            block = features.iloc[start : start + args.block_size]
            print(f"Block {block_id}: {len(block)} genes", flush=True)
            result = summarize_block(census, block, mapping, groups)
            temp = path.with_suffix(".partial")
            result.to_parquet(temp, index=False)
            temp.replace(path)
            print(f"Saved {len(result):,} positive gene/population rows", flush=True)
    print("Requested blocks complete", flush=True)


if __name__ == "__main__":
    main()

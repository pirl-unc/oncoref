#!/usr/bin/env python3
"""Count Target-candidate genes in adult normal nose, skin, scalp and tongue cells.

Reads the pinned CELLxGENE Census used by the whole-body atlas, restricted to the
adult, UMI-protocol populations of ``populations.parquet``. Writes per-population
positive cells, multi-molecule cells and count sums for ``build_nose_skin_reference``.

    python scripts/extract_nose_skin_cell_counts.py \\
        tmp/normal-atlas/wholebody-2025-11-08 tmp/normal-atlas/nose-skin-2025-11-08

Needs ``cellxgene_census`` (``tmp/normal-atlas-venv``) and network access.
"""

import json
import sys
from pathlib import Path

import cellxgene_census as cc
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ATLAS, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
CONTROLS = ["KRT5", "KRT14", "KRT1", "KRT10", "SCGB1A1", "BPIFA1", "MUC5AC", "FOXJ1", "KRT85"]
TISSUES = ["nose", "skin of body", "scalp", "tongue"]

pops = pd.read_parquet(ATLAS / "populations.parquet")
keep = pops[
    pops.tissue_general.isin(TISSUES)
    & (pops.age_scope == "adult")
    & (pops.count_unit == "UMI_protocol")
]
mapping = np.load(ATLAS / "observation-population-map.npy", mmap_mode="r")
catalog = pd.read_csv(ROOT / "oncoref/data/target-selection.csv.gz", low_memory=False)
genes = sorted(
    {g for ids in catalog[catalog.is_target_candidate].member_gene_ids for g in json.loads(ids)}
)
controls = CONTROLS

with cc.open_soma(census_version="2025-11-08") as census:
    human = census["census_data"]["homo_sapiens"]
    var = human.ms["RNA"].var.read().concat().to_pandas()
    var = var[var.feature_id.isin(genes) | var.feature_name.isin(controls)]
    var_ids = var.soma_joinid.to_numpy()
    col = np.full(int(var_ids.max()) + 1, -1, np.int32)
    col[var_ids] = np.arange(len(var_ids))
    obs = (
        human.obs.read(
            value_filter=(
                "is_primary_data == True and disease == 'normal' and tissue_type == 'tissue' "
                f"and tissue_general in {TISSUES}"
            ),
            column_names=["soma_joinid", "raw_sum"],
        )
        .concat()
        .to_pandas()
    )
    obs["population_id"] = mapping[obs.soma_joinid.to_numpy()]
    obs = obs[obs.population_id.isin(keep.population_id)]
    print(len(obs), "cells", obs.population_id.nunique(), "populations", flush=True)
    pid = np.sort(obs.population_id.unique())
    row_of = pd.Series(np.arange(len(pid)), index=pid)
    shape = (len(pid), len(var_ids))
    n_pos = np.zeros(shape, np.int32)
    n_ge2 = np.zeros(shape, np.int32)
    total = np.zeros(shape, np.float64)
    ids = np.sort(obs.soma_joinid.to_numpy())
    for start in range(0, len(ids), 50_000):
        chunk = ids[start : start + 50_000]
        for table in human.ms["RNA"].X["raw"].read(coords=(chunk, var_ids)).tables():
            r = row_of.reindex(mapping[table["soma_dim_0"].to_numpy()]).to_numpy()
            c = col[table["soma_dim_1"].to_numpy()]
            v = table["soma_data"].to_numpy()
            np.add.at(n_pos, (r, c), 1)
            np.add.at(n_ge2, (r, c), v >= 2)
            np.add.at(total, (r, c), v)
        print(start + len(chunk), flush=True)

depth = obs.groupby("population_id").raw_sum.sum().reindex(pid)
cells = obs.groupby("population_id").size().reindex(pid)
r, c = np.nonzero(n_pos)
long = pd.DataFrame(
    {
        "population_id": pid[r],
        "gene_id": var.feature_id.to_numpy()[c],
        "gene_name": var.feature_name.to_numpy()[c],
        "n_positive": n_pos[r, c],
        "n_ge2": n_ge2[r, c],
        "count_sum": total[r, c],
    }
)
long.to_parquet(OUT / "counts.parquet")
meta = keep.merge(
    pd.DataFrame({"population_id": pid, "cells": cells.to_numpy(), "raw_sum": depth.to_numpy()}),
    on="population_id",
)
meta.to_parquet(OUT / "populations.parquet")
var.to_csv(OUT / "var.csv", index=False)
print("done", len(long), len(meta))

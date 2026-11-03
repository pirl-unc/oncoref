# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Rare expression in vital-tissue cell types, as a filter on every protein group.

A gene that reads near zero in bulk and cluster panels can still be expressed by a
small population of an essential tissue (MAGEA12 in rare midbrain dopaminergic
neurons). ``scripts/build_cta_rare_cell_screen.py`` records, for every gene, the
brain, spinal cord, heart and lung (alveolar, endothelial, airway basal) cell types with the highest fraction
of expressing cells in CELLxGENE normal samples. Immune lineages, stroma, adipose
and developmental cell types are ignored, so rare expression in spleen, adipose
or immune subsets is allowed.

A protein group is flagged when any member gene has at least 1% of a vital cell type
expressing it. Detection frequency is not abundance and can include ambient RNA,
doublets or paralog multimapping; the filter is an automated screen. Unscreened
groups (genes absent from the service) are reported as such and never flagged.
"""

from __future__ import annotations

import hashlib
import json
from importlib.resources import files

import pandas as pd

TABLE, PROVENANCE = "cta-rare-cell-screen.csv.gz", "cta-rare-cell-screen.json"
FLAG_PCT = 1.0
#: Positive expression in an essential tissue is Unsafe evidence, like a reviewed
#: unresolved Protected-tissue finding. Set False to make it a review hold instead.
UNSAFE = True


def provenance() -> dict:
    return json.loads(files("oncoref").joinpath("data", PROVENANCE).read_text())


def screen_table() -> pd.DataFrame:
    """The packaged per-gene screen, verified against its provenance checksum."""
    path = files("oncoref").joinpath("data", TABLE)
    if hashlib.sha256(path.read_bytes()).hexdigest() != provenance()["sha256"]:
        raise ValueError("Rare-cell screen table does not match its provenance")
    with path.open("rb") as handle:
        return pd.read_csv(handle, compression="gzip")


def annotate(facts: pd.DataFrame, *, member_column: str = "member_gene_ids") -> pd.DataFrame:
    """Add rare-cell columns for each protein group from its member genes."""
    table = screen_table()
    screened = set(provenance()["screened_gene_ids"])
    table = table[table.tier.eq("vital")]  # Manageable cell types are recorded, never flagged
    top = (
        table.sort_values("pct_expressing", ascending=False)
        .drop_duplicates("gene_id")
        .set_index("gene_id")
    )
    members = facts[member_column].map(json.loads)

    def best(ids):
        hits = [g for g in ids if g in top.index]
        return max(hits, key=lambda g: top.pct_expressing[g], default=None)

    chosen = members.map(best)
    has = chosen.map(lambda g: isinstance(g, str))
    out = facts.copy()
    out["rare_cell_screened"] = members.map(lambda ids: any(g in screened for g in ids)).astype(
        bool
    )
    out["rare_cell_max_pct"] = chosen.map(
        lambda g: top.pct_expressing[g] if isinstance(g, str) else float("nan")
    )
    out["rare_cell_type"] = chosen.map(
        lambda g: f"{top.tissue[g]}: {top.cell_type[g]}" if isinstance(g, str) else ""
    )
    out["rare_cell_detail"] = chosen.map(
        lambda g: (
            f"{int(top.n_expressing[g]):,} of {int(top.n_cells[g]):,} cells"
            if isinstance(g, str)
            else ""
        )
    )
    out["rare_cell_flag"] = (has & out.rare_cell_max_pct.ge(FLAG_PCT)).astype(bool)
    return out

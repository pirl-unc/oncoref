#!/usr/bin/env python3
"""Tabulate every filter's status for every protein group, including Unsafe ones.

Every layer evaluates every protein, so a protein that is already Unsafe still has
a pass/hold/unsafe/fail status for each other filter. This table exposes them.

Statuses
- target gates (reproductive association, cancer expression, coverage, concentration):
  ``pass`` or ``fail``. A failure blocks CTA only.
- safety layers: ``pass``, ``hold`` (blocks CTA) or ``unsafe`` (positive exclusion).
- Unsafe sequence overlap: ``fail`` at 95% or more of 9-mers inside Unsafe proteins.
  It removes a would-be CTA only. It is ``not_applicable`` for Unsafe proteins, which are
  inside the Unsafe set by definition.
- Rare vital-tissue cell expression: ``unsafe`` at 1% or more of a brain, spinal cord,
  heart or lung alveolar, endothelial or airway basal cell type, ``pass``, or ``no_data`` for genes absent from the
  single-cell service (never flagged).
- Universe: out-of-universe groups have every other filter ``not_evaluated``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
AUDIT = ROOT / "docs/audits/cta-selection-20261028"

from oncoref.cta_selection_policy import CANDIDATE_ONLY_LAYERS, CONTAINED_FLOOR  # noqa: E402


def build(audit: Path = AUDIT) -> tuple[pd.DataFrame, ...]:
    facts = pd.read_csv(audit / "protein-group-classifications.csv.gz", low_memory=False)
    layers = [c[5:] for c in facts if c.startswith("hold:")]
    table = facts[["protein_equivalence_group_id", "protein_group_label", "member_gene_ids"]].copy()
    table = table.rename(columns={"protein_group_label": "label", "member_gene_ids": "gene_ids"})
    table["universe"] = "in_universe"
    table["protein_source"] = facts.protein_source
    unsafe_flags = []
    for layer in layers:
        hold, unsafe = facts["hold:" + layer].astype(bool), facts["unsafe:" + layer].astype(bool)
        if layer in CANDIDATE_ONLY_LAYERS:
            status = pd.Series("pass", index=facts.index).mask(hold | unsafe, "fail")
        else:
            status = pd.Series("pass", index=facts.index).mask(hold, "hold").mask(unsafe, "unsafe")
            unsafe_flags.append(unsafe.rename(layer))
        table["filter: " + layer] = status
    # An Unsafe protein lies wholly inside the Unsafe set by definition, so the rule
    # is not applicable to it; the overlap is meaningful for every other protein.
    is_unsafe = facts.classification.eq("Unsafe")
    overlap = facts.unsafe_sequence_overlap_pct.mask(is_unsafe)
    table["filter: Unsafe sequence overlap"] = (
        pd.Series("pass", index=facts.index)
        .mask(overlap.ge(100 * CONTAINED_FLOOR), "fail")
        .mask(overlap.isna(), "no_9mer")
        .mask(is_unsafe, "not_applicable")
    )
    table["unsafe_sequence_overlap_pct"] = overlap
    # Every protein is screened for rare vital-tissue cell expression; genes absent from
    # the single-cell service have no data and are never flagged.
    rare = "filter: Rare vital-tissue cell expression"
    table.loc[~facts.rare_cell_screened.astype(bool), rare] = "no_data"
    table["rare_cell_max_pct"] = facts.rare_cell_max_pct
    table["rare_cell_type"] = facts.rare_cell_type.fillna("")
    unsafe_matrix = pd.concat(unsafe_flags, axis=1)
    table["n_unsafe_layers"] = unsafe_matrix.sum(axis=1)
    table["sole_unsafe_layer"] = unsafe_matrix.apply(
        lambda row: row.index[row][0] if row.sum() == 1 else ("multiple" if row.sum() else ""),
        axis=1,
    )
    blockers = pd.DataFrame(
        {
            c[8:]: table[c].isin(["fail", "hold", "unsafe"])
            for c in table
            if c.startswith("filter: ")
        }
    )
    table["blocking_filters"] = blockers.apply(lambda r: "; ".join(r.index[r]), axis=1)
    table["n_blocking_filters"] = blockers.sum(axis=1)
    table["sole_blocking_filter"] = blockers.apply(
        lambda r: r.index[r][0] if r.sum() == 1 else "", axis=1
    )
    table["classification"] = facts.classification
    table["assessment_status"] = facts.assessment_status
    # Out-of-universe groups: listed, every filter not evaluated.
    out = pd.read_csv(audit / "out-of-universe.csv.gz")
    rows = out.rename(columns={"protein_group_label": "label"}).copy()
    rows["universe"] = "out_of_universe"
    rows["classification"] = "out_of_universe"
    rows["assessment_status"] = rows.reason
    for column in table:
        if column.startswith("filter: "):
            rows[column] = "not_evaluated"
    full = pd.concat([table, rows[[c for c in table if c in rows]]], ignore_index=True)
    full["gene_ids"] = full.gene_ids.fillna("[]")
    # Per-filter summary.
    filters = [c for c in full if c.startswith("filter: ")]
    summary = (
        full.melt(
            id_vars="classification", value_vars=filters, var_name="filter", value_name="status"
        )
        .assign(filter=lambda d: d["filter"].str.removeprefix("filter: "))
        .groupby(["filter", "status"])
        .size()
        .rename("proteins")
        .reset_index()
    )
    sole = (
        full[full.n_blocking_filters.eq(1)]
        .sole_blocking_filter.value_counts()
        .rename_axis("filter")
        .rename("proteins_blocked_only_by_this_filter")
        .reset_index()
    )
    unsafe_only = (
        full[full.classification.eq("Unsafe")]
        .sole_unsafe_layer.value_counts()
        .rename_axis("layer")
        .rename("unsafe_proteins_with_this_as_sole_layer")
        .reset_index()
    )
    return full, summary, unsafe_only, sole


if __name__ == "__main__":
    full, summary, unsafe_only, sole = build()
    full.to_csv(
        AUDIT / "filter-status-by-protein.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    summary.to_csv(AUDIT / "filter-status-summary.csv", index=False)
    sole.to_csv(AUDIT / "filter-status-sole-blocker.csv", index=False)
    unsafe_only.to_csv(AUDIT / "filter-status-unsafe-sole-layer.csv", index=False)
    genes = full.assign(gene_id=full.gene_ids.map(json.loads)).explode("gene_id")
    genes = genes[genes.gene_id.notna()][
        ["gene_id", "protein_equivalence_group_id", "label", "classification"]
    ]
    genes.to_csv(
        AUDIT / "filter-status-gene-map.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    print(len(full), "protein groups;", genes.gene_id.nunique(), "genes")

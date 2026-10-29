"""Shared pinned RNA profile loading for selection and diffuseness reports."""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from build_cta_alternative_evidence import GTEX_TISSUES

from oncoref import expression_diffuseness as metric
from oncoref import gtex, hpa

ROOT = Path(__file__).resolve().parents[1]
CURRENT = ROOT / "docs/audits/cta-selection-20261028"
ID = "protein_equivalence_group_id"


def hpa_profiles(facts, current=CURRENT):
    """Sum identical-protein members within each native tissue, with full coverage."""
    members = {
        g: r.protein_equivalence_group_id
        for r in facts.itertuples()
        for g in json.loads(r.member_gene_ids)
    }
    counts = facts.set_index(ID).member_count
    raw = hpa.hpa_rna_consensus("v23").copy()
    ids = pd.read_csv(current / "hpa_rna_consensus-accepted-identities.csv.gz")
    raw["gene_id"] = raw.Gene.map(ids.set_index("source_gene_id").canonical_gene_id)
    raw = raw[raw.gene_id.isin(members)].copy()
    raw["entity_id"] = raw.gene_id.map(members)
    grouped = (
        raw.groupby(["entity_id", "Tissue"])
        .agg(value=("nTPM", "sum"), members=("gene_id", "nunique"))
        .reset_index()
    )
    grouped.loc[grouped.members.ne(grouped.entity_id.map(counts)), "value"] = np.nan
    hp = grouped.pivot(index="entity_id", columns="Tissue", values="value").reindex(facts[ID])
    hp.index.name = "entity_id"
    return hp


def profiles(facts):
    """Load independent HPA and GTEx protein profiles without pooling assays."""
    hp = hpa_profiles(facts)
    members = {
        g: r.protein_equivalence_group_id
        for r in facts.itertuples()
        for g in json.loads(r.member_gene_ids)
    }
    counts = facts.set_index(ID).member_count
    native = gtex.tissue_rna("v11")
    tissue_map = {
        re.sub(r"[ ()-]+", "_", name.replace("c-1", "cMINUS1"))
        .strip("_")
        .replace("cMINUS1", "c-1"): tissue
        for name, tissue in GTEX_TISSUES.items()
    }
    scope = json.loads((CURRENT / "alternative-evidence/gtex-scope.json").read_text())
    if set(tissue_map) != set(scope["included_bulk_contexts"]):
        raise ValueError("GTEx core does not match the pinned adapter scope")
    accepted = pd.read_csv(CURRENT / "alternative-evidence/gtex_tissue_rna-row-decisions.csv.gz")
    accepted = accepted[accepted.row_position.notna() & accepted.gene_id.isin(members)]
    gp = native.iloc[accepted.row_position.astype(int)][list(tissue_map)].copy()
    gp.index = accepted.gene_id.map(members).to_numpy()
    present = gp.groupby(level=0).count()
    gp = gp.groupby(level=0).sum(min_count=1)
    gp = gp.where(present.eq(counts.reindex(gp.index), axis=0)).reindex(facts[ID])
    gp.index.name = "entity_id"
    return hp, gp, tissue_map


def summarize(matrix, reference_columns, panel, source, version, unit, statistic):
    """Use the shared continuous metric without borrowing references across assays."""
    meta = {
        "source": source,
        "version": version,
        "taxonomy_id": 9606,
        "unit": unit,
        "cohort": source,
        "summary": statistic,
    }
    data = (
        matrix[list(panel.compartment)]
        .reset_index()
        .melt(id_vars="entity_id", var_name="compartment", value_name="value")
        .assign(**meta)
    )
    reference = (
        pd.DataFrame(
            {
                "reference_value": matrix[reference_columns].max(axis=1),
                "reference_complete": matrix[reference_columns].notna().all(axis=1),
            }
        )
        .reset_index()
        .assign(**meta)
    )
    _, summary = metric.assess(data, reference, panel)
    summary["core_complete"] = summary.entity_id.map(
        matrix[[*reference_columns, *panel.compartment]].notna().all(axis=1)
    )
    return summary


def hpa_summary(facts, current=CURRENT):
    """Selection's fixed HPA v23 panel; thymus is never a positive reference."""
    matrix = hpa_profiles(facts, current)
    refs = sorted(t for ts in metric.HPA_REPRODUCTIVE_BRANCHES.values() for t in ts)
    if set(matrix) != set(metric.hpa_bulk_roles().compartment):
        raise ValueError("Unexpected HPA v23 tissue vocabulary")
    return summarize(
        matrix, refs, metric.hpa_bulk_panel(), "hpa_rna_consensus", "v23", "nTPM", "consensus"
    )

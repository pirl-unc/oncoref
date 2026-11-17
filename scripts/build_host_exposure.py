#!/usr/bin/env python3
"""Build host exposure observations independently of positive RA qualification."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import host_context, reference_data, target_selection  # noqa: E402

OUT = ROOT / "docs/audits/target-selection-20261110/rebuild"


def summarize(observations, *, compartment):
    """Retain finite identity-resolved measurements, including zeros and sparse positives.

    RA's positive-cell/donor criterion must never gate a negative reference.
    Embryo scope remains incomplete regardless of the measured maximum.
    """
    good = observations.identity_resolved & observations.value.notna()
    if "library_qc_review" in observations:
        good &= observations.library_qc_review.fillna("").eq("")
    measured = observations.loc[good]
    values = measured.value
    if not np.isfinite(values).all() or values.lt(0).any():
        raise ValueError("Host expression must be finite and nonnegative")
    return {
        "coverage": "measured" if len(measured) else "missing",
        "max_tpm": float(values.max()) if len(measured) else None,
        "measured_observations": len(measured),
        "unresolved_observations": int((~good).sum()),
        "scope_complete": bool(len(measured) and compartment != "embryo"),
        "scope": "early embryo only; later fetal/developmental stages unavailable" if compartment == "embryo" else "sampled adult tissue/state; no universal donor or cell-type coverage claim",
        "sources": sorted(measured.source_id.unique()),
        "max_ntpm": None,
        "protein_status": "unassessed",
    }


def build():
    inputs = {}

    def read(path, **kwargs):
        path = Path(path)
        inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        return pd.read_csv(path, **kwargs)

    catalog = target_selection.results()
    observations = read(OUT / "reproductive-comparison.csv.gz", low_memory=False)
    evidence = {}
    for entity, rows in observations.groupby("entity_id", sort=True):
        evidence[entity] = {
            c: summarize(rows[rows.compartment.eq(c)], compartment=c)
            for c in host_context.COMPARTMENTS
        }
    # Reuse the catalog's source-to-canonical identity adjudication, not symbol guesses.
    identity = read(ROOT / "docs/audits/cta-selection-20261028/hpa_rna_consensus-accepted-identities.csv.gz")
    lookup = identity.set_index("source_gene_id").canonical_gene_id
    members = {r.entity_id: json.loads(r.member_gene_ids) for r in catalog.itertuples()}
    gene_entity = {g: e for e, genes in members.items() for g in genes}
    sizes = {e: len(genes) for e, genes in members.items()}
    for source in ("hpa_single_cell_tissue", "hpa_normal_tissue"):
        path = reference_data.local_path(source, "v23")
        data = read(path, sep="\t")
        data["entity_id"] = data.Gene.map(lookup).map(gene_entity)
        data = data[data.entity_id.isin(evidence)]
        for (entity, tissue), rows in data.groupby(["entity_id", "Tissue"], sort=True):
            compartments = catalog_compartments(tissue)
            for compartment in compartments:
                item = evidence[entity][compartment]
                item.setdefault("orthogonal_sources", []).append(source + ":v23")
                if source == "hpa_normal_tissue":
                    positive = rows.Level.isin(["Low", "Medium", "High"]) & rows.Reliability.isin(["Approved", "Supported", "Enhanced"])
                    item["protein_status"] = "unresolved" if positive.any() else "not_detected" if rows.Level.eq("Not detected").all() else "unassessed"
                    item["protein_findings"] = rows.loc[positive, ["Gene", "Cell type", "Level", "Reliability"]].to_dict("records")
                    # These are observations requiring attribution review, not new Unsafe calls.
                else:
                    rows = rows.copy()
                    keys = [c for c in ("Cluster", "Cell type") if c in rows]
                    grouped = rows.groupby(keys).agg(total=("nTPM", "sum"), measured=("Gene", "nunique"))
                    complete = grouped.measured.eq(sizes[entity])
                    if complete.any():
                        item["max_ntpm"] = float(grouped.loc[complete, "total"].max())
                    if not complete.all():
                        item["scope_complete"] = False
                    item["cell_rna_scope"] = "unlocalized mammary state" if tissue == "breast" else "named tissue cell clusters"
    rows = [{"entity_id": e, "exposures": json.dumps(v, allow_nan=False, sort_keys=True)} for e, v in sorted(evidence.items())]
    target = ROOT / "oncoref/data/host-exposure.csv.gz"
    pd.DataFrame(rows).to_csv(target, index=False, compression={"method": "gzip", "mtime": 0})
    metadata = {
        "schema": "oncoref.host-exposure.v1", "policy_id": host_context.POLICY_ID,
        "catalog_policy_id": target_selection.POLICY_ID,
        "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
        "inputs": inputs,
        "limitations": ["Native TPM and nTPM are screened separately, never converted.", "Positive HPA IHC remains unresolved host-specific protein evidence until attribution is reviewed.", "Early embryo observations cannot clear all stages of pregnancy.", "Mammary state is expression-inferred; state-unspecified orthogonal findings apply conservatively to both states.", "Unmeasured cell populations and treatment-induced states are not excluded by this research screen."],
    }
    (ROOT / "oncoref/data/host-exposure.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"host_exposure_groups": len(rows), "policy_id": host_context.POLICY_ID}))


def catalog_compartments(tissue):
    if tissue == "breast":
        return ("lactating_breast", "resting_breast")
    return (tissue,) if tissue in host_context.COMPARTMENTS else ()


if __name__ == "__main__":
    build()

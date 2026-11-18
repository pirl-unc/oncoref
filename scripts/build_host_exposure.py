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
from oncoref import curation, host_context, reference_data, target_selection  # noqa: E402

OUT = ROOT / "docs/audits/target-selection-20261110/rebuild"


def summarize(observations, *, compartment):
    """Retain finite identity-resolved measurements, including zeros and sparse positives.

    RA's positive-cell/donor criterion must never gate a negative reference.
    Embryo scope remains incomplete regardless of the measured maximum.
    """
    good = observations.identity_resolved & observations.value.notna()
    qc_review = (
        int(observations.library_qc_review.fillna("").ne("").sum())
        if "library_qc_review" in observations
        else 0
    )
    measured = observations.loc[good]
    values = measured.value
    if not np.isfinite(values).all() or values.lt(0).any():
        raise ValueError("Host expression must be finite and nonnegative")
    return {
        "coverage": "measured" if len(measured) else "missing",
        "max_tpm": float(values.max()) if len(measured) else None,
        "measured_observations": len(measured),
        "unresolved_observations": int((~good).sum()),
        "qc_review_observations": qc_review,
        "scope_complete": bool(len(measured) and compartment != "embryo" and not qc_review),
        "scope": "early embryo only; later fetal/developmental stages unavailable"
        if compartment == "embryo"
        else "sampled adult tissue/state; no universal donor or cell-type coverage claim",
        "sources": sorted(measured.source_id.unique()),
        "max_ntpm": None,
        "protein_status": "unassessed",
    }


def build():
    inputs = {}

    def read(path, *, source_key=None, pin=None, **kwargs):
        path = Path(path)
        with path.open("rb") as handle:
            checksum = hashlib.sha256()
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                checksum.update(chunk)
            digest = checksum.hexdigest()
        if pin and (digest != pin["sha256"] or path.stat().st_size != pin["bytes"]):
            raise ValueError(f"Host source does not match its immutable pin: {source_key}")
        key = source_key or (
            str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        )
        inputs[key] = digest
        return pd.read_csv(path, **kwargs)

    catalog = target_selection.results()
    observations = read(OUT / "reproductive-comparison.csv.gz", low_memory=False)
    evidence = {}
    for entity, rows in observations.groupby("entity_id", sort=True):
        evidence[entity] = {}
        for c in host_context.COMPARTMENTS:
            unknown_state = (
                rows.compartment.eq("breast")
                if c in target_selection.MAMMARY_COMPARTMENTS
                else pd.Series(False, index=rows.index)
            )
            evidence[entity][c] = summarize(
                rows[rows.compartment.eq(c) | unknown_state], compartment=c
            )
            if unknown_state.any():
                evidence[entity][c]["scope_complete"] = False
                evidence[entity][c]["mammary_state_unknown"] = True
    # Reuse the catalog's source-to-canonical identity adjudication, not symbol guesses.
    identity = read(
        ROOT / "docs/audits/cta-selection-20261028/hpa_rna_consensus-accepted-identities.csv.gz"
    )
    lookup = identity.set_index("source_gene_id").canonical_gene_id
    members = {r.entity_id: json.loads(r.member_gene_ids) for r in catalog.itertuples()}
    gene_entity = {g: e for e, genes in members.items() for g in genes}
    sizes = {e: len(genes) for e, genes in members.items()}
    references = curation.manifest()["reference_inputs"]
    for role in ("host_cell_rna", "host_protein"):
        source, version = references[role]["source"], references[role]["version"]
        path = reference_data.local_path(source, version)
        pin = reference_data.REFERENCE_SOURCES[source]["pins"][version]
        data = read(path, sep="\t", source_key=f"reference:{source}:{version}", pin=pin)
        data["entity_id"] = data.Gene.map(lookup).map(gene_entity)
        data = data[data.entity_id.isin(evidence)]
        data["compartment"] = data.Tissue.map(catalog_compartments)
        data = data.explode("compartment").dropna(subset=["compartment"])
        # Combine every source label for a compartment before grading. A negative
        # alias must never overwrite a positive finding from another label.
        for (entity, compartment), rows in data.groupby(["entity_id", "compartment"], sort=True):
            item = evidence[entity][compartment]
            item.setdefault("orthogonal_sources", []).append(source + ":" + version)
            if source == "hpa_normal_tissue":
                item.update(protein_evidence(rows))
            else:
                keys = [c for c in ("Tissue", "Cluster", "Cell type") if c in rows]
                grouped = rows.groupby(keys).agg(
                    total=("nTPM", "sum"), measured=("Gene", "nunique")
                )
                complete = grouped.measured.eq(sizes[entity])
                if complete.any():
                    item["max_ntpm"] = float(grouped.loc[complete, "total"].max())
                if not complete.all():
                    item["scope_complete"] = False
                item["cell_rna_scope"] = (
                    "unlocalized mammary state"
                    if rows.Tissue.eq("breast").any()
                    else "named tissue cell clusters"
                )
    rows = [
        {"entity_id": e, "exposures": json.dumps(v, allow_nan=False, sort_keys=True)}
        for e, v in sorted(evidence.items())
    ]
    target = ROOT / "oncoref/data/host-exposure.csv.gz"
    pd.DataFrame(rows).to_csv(target, index=False, compression={"method": "gzip", "mtime": 0})
    metadata = {
        "schema": "oncoref.host-exposure.v1",
        "policy_id": host_context.POLICY_ID,
        "catalog_policy_id": target_selection.POLICY_ID,
        "catalog_data_sha256": target_selection.metadata()["data_sha256"],
        "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
        "inputs": inputs,
        "reference_inputs": {role: references[role] for role in ("host_cell_rna", "host_protein")},
        "limitations": [
            "Native TPM and nTPM are screened separately, never converted.",
            "Positive HPA IHC remains unresolved host-specific protein evidence until attribution is reviewed.",
            "Early embryo observations cannot clear all stages of pregnancy.",
            "Mammary state is expression-inferred; state-unspecified orthogonal findings apply conservatively to both states.",
            "Unmeasured cell populations and treatment-induced states are not excluded by this research screen.",
        ],
    }
    (ROOT / "oncoref/data/host-exposure.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"host_exposure_groups": len(rows), "policy_id": host_context.POLICY_ID}))


def catalog_compartments(tissue):
    # HPA distinguishes numbered endometrial panels and a lactating breast label.
    # Retain their native label in findings while mapping to catalog compartments.
    tissue = {
        "endometrium 1": "endometrium",
        "endometrium 2": "endometrium",
        "lactating breast": "lactating_breast",
    }.get(tissue, tissue)
    if tissue == "breast":
        return ("lactating_breast", "resting_breast")
    return (tissue,) if tissue in host_context.COMPARTMENTS else ()


def protein_evidence(rows):
    """Pool source labels without allowing negative panels to erase positive ones."""
    positive = rows.Level.isin(["Low", "Medium", "High"]) & rows.Reliability.isin(
        ["Approved", "Supported", "Enhanced"]
    )
    return {
        "protein_status": "unresolved"
        if positive.any()
        else "not_detected"
        if rows.Level.eq("Not detected").all()
        else "unassessed",
        "protein_findings": rows.loc[
            positive, ["Gene", "Tissue", "Cell type", "Level", "Reliability"]
        ].to_dict("records"),
    }


if __name__ == "__main__":
    build()

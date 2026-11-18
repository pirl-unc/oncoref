"""Pinned individual tissue RNA observations for the symmetric antigen screen."""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from cta_report_common import sha256, write_csv

ROOT = Path(__file__).resolve().parents[1]
HPA_SOURCE = {
    "source_id": "hpa_tissue_samples:v23",
    "url": "https://v23.proteinatlas.org/download/transcript_rna_tissue.tsv.zip",
    "sha256": "4f4790de54380f4d47800b16efaed5a0d6bc364a90ddfdd0582f8f2d13e4f9e1",
    "reference": "HPA v23, Ensembl 109",
    "unit": "TPM",
    "statistic": "individual tissue sample; summed transcript TPM",
    "participant_identity": "Not provided; sample identifiers are not assumed to be donor IDs",
}
LACTATION_MARKERS = ("LALBA", "CSN1S1", "CSN2", "CSN3")
LACTATION_MARKER_TPM = 100.0
LACTATION_MIN_MARKERS = 3


def mammary_sample_states(genes, accepted_ids) -> dict:
    """Expression-inferred state; neither signature absence nor missing data proves rest."""
    from oncoref.gene_matching import gene_identity_catalog

    breast_cols = [c for c in genes.columns if c.removeprefix("TPM.").rsplit(".", 1)[0] == "breast"]
    if not breast_cols:
        return {}
    lookup = accepted_ids.set_index("source_gene_id").canonical_gene_id
    mapped = genes.copy()
    mapped.index = mapped.index.map(lookup)
    mapped = mapped[mapped.index.notna()]
    catalog = gene_identity_catalog()
    marker_ids = catalog.loc[
        catalog.canonical_symbol.isin(LACTATION_MARKERS), "canonical_gene_id"
    ].tolist()
    if len(marker_ids) != len(LACTATION_MARKERS):
        raise ValueError("Missing canonical gene identities for lactation markers")
    panel = mapped.reindex(marker_ids)[breast_cols]
    states = {}
    for column in panel:
        values = panel[column]
        complete = values.notna().all() and np.isfinite(values).all() and values.ge(0).all()
        state = (
            "unknown"
            if not complete
            else "secretory_signature"
            if values.ge(LACTATION_MARKER_TPM).sum() >= LACTATION_MIN_MARKERS
            else "no_high_secretory_signature"
        )
        states[column.removeprefix("TPM.")] = {
            "state": state,
            "inference": "milk-program-v1",
            "provider_confirmed": False,
            "marker_tpm": {
                gene: float(value) if pd.notna(value) and np.isfinite(value) else None
                for gene, value in values.items()
            },
        }
    return states


def lactating_mammary_samples(genes, accepted_ids) -> frozenset[str]:
    """Compatibility name: samples with an inferred secretory/lactation-like signature."""
    return frozenset(
        sample
        for sample, record in mammary_sample_states(genes, accepted_ids).items()
        if record["state"] == "secretory_signature"
    )


def sum_transcript_chunks(chunks):
    """Sum transcripts within each source gene and sample, preserving missingness."""
    pieces, incomplete = [], []
    for chunk in chunks:
        values = chunk.drop(columns="ensgid").apply(pd.to_numeric, errors="raise")
        if (values.lt(0) | (np.isinf(values))).any().any():
            raise ValueError("Transcript TPM must be finite and nonnegative")
        values.index = chunk.ensgid
        pieces.append(values.groupby(level=0, sort=False).sum(min_count=1))
        incomplete.append(values.isna().groupby(level=0, sort=False).any())
    summed = pd.concat(pieces).groupby(level=0, sort=False).sum(min_count=1)
    missing = pd.concat(incomplete).groupby(level=0, sort=False).any()
    return summed.mask(missing)


def tissue_sample_genes(archive=None):
    """Verify the provider archive and extract gene TPM without renormalizing it."""
    from oncoref.curation import reference_input

    archive = reference_input("hpa_tissue_samples") if archive is None else archive
    if sha256(Path(archive)) != HPA_SOURCE["sha256"]:
        raise ValueError("HPA tissue-sample source checksum mismatch")
    chunks = pd.read_csv(
        archive,
        sep="\t",
        chunksize=5000,
        usecols=lambda column: column == "ensgid" or column.startswith("TPM."),
    )
    genes = sum_transcript_chunks(chunks)
    totals = genes.sum()
    if not totals.between(999_990, 1_000_010).all():
        raise ValueError("HPA sample TPM totals do not reproduce the source scale")
    return genes


def tissue_samples(facts, accepted_ids, *, genes=None):
    """Group exact member loci in each real sample; missing members never become zero."""
    raw_genes = tissue_sample_genes() if genes is None else genes.copy()
    mammary_states = mammary_sample_states(raw_genes, accepted_ids)
    lactating = {s for s, r in mammary_states.items() if r["state"] == "secretory_signature"}
    genes = raw_genes
    lookup = accepted_ids.set_index("source_gene_id").canonical_gene_id
    genes.index = genes.index.map(lookup)
    genes = genes[genes.index.notna()]
    if genes.index.duplicated().any():
        raise ValueError("Multiple HPA source genes map to one canonical gene")
    available = set(genes.index)
    rows = []
    for protein in facts.itertuples():
        members = json.loads(protein.member_gene_ids)
        present = genes.reindex(members)
        values = present.sum(min_count=len(members))
        for column, value in values.items():
            sample = column.removeprefix("TPM.")
            tissue, identifier = sample.rsplit(".", 1)
            if not identifier.isdigit():
                raise ValueError(f"Unexpected HPA sample identifier: {sample}")
            compartment = (
                "lactating_breast"
                if sample in lactating
                else (
                    "resting_breast"
                    if tissue == "breast" and mammary_states[sample]["state"] != "unknown"
                    else tissue
                )
            )
            rows.append(
                {
                    "entity_id": protein.entity_id,
                    "source_id": HPA_SOURCE["source_id"],
                    "sample_id": sample,
                    "participant_id": "",
                    "tissue": tissue,
                    "compartment": compartment,
                    "mammary_state_inference": json.dumps(
                        mammary_states.get(sample, {}), allow_nan=False
                    ),
                    "source_context": sample,
                    "cell_type": "bulk",
                    "statistic": HPA_SOURCE["statistic"],
                    "observation_kind": "individual_tissue_sample",
                    "unit": "TPM",
                    "value": value,
                    "qc_pass": pd.notna(value),
                    "identity_resolved": set(members) <= available,
                }
            )
    return pd.DataFrame(rows)


def tissue_library_qc(genes, accepted_ids):
    """Audit within tissue, with reproductive programs recorded as notes only.

    These are tumor-screen signatures, not validated normal-library exclusion
    criteria. Flagged normals remain in the denominator pending adjudication.
    """
    from oncoref import library_qc
    from oncoref.gene_matching import gene_identity_catalog

    lookup = accepted_ids.set_index("source_gene_id").canonical_gene_id
    values = genes.copy()
    values.index = values.index.map(lookup)
    values = values[values.index.notna()]
    if not values.index.is_unique:
        raise ValueError("Multiple HPA source genes map to one canonical gene")
    labels = gene_identity_catalog().set_index("canonical_gene_id").canonical_symbol
    symbols = labels.reindex(values.index).fillna("")
    rows = []
    columns = pd.Series(values.columns, index=values.columns)
    for tissue, samples in columns.groupby(lambda c: c.removeprefix("TPM.").rsplit(".", 1)[0]):
        audit = library_qc.signatures(
            values, list(samples), gene_ids=pd.Series(values.index), symbols=symbols
        )
        audit["tissue"] = tissue
        audit["library_qc_review"] = audit.apply(
            lambda r: ";".join(library_qc.reasons(r)),
            axis=1,
        )
        audit["library_qc_notes"] = audit.apply(lambda r: ";".join(library_qc.notes(r)), axis=1)
        rows.append(audit)
    audit = pd.concat(rows).reset_index()
    audit["sample_id"] = audit.sample_id.str.removeprefix("TPM.")
    return audit


def developmental_samples(facts):
    """TPM means within one embryo and lineage, never across biological donors.

    A well-sampled population needs >=20 cells, >=5 positive and >=10% positive.
    Smaller populations need >=5 cells, >=50% positive and two independent embryos
    at the same stage and lineage reaching the fold. Cells are not independent witnesses.
    """
    from scipy.io import mmread

    from oncoref.curation import reference_input
    from oncoref.gene_matching import match_reference_genes, resolve_measurement_rows
    from oncoref.target_selection import developmental_donor_requirement

    sources = [
        ("petropoulos_E-MTAB-3929", "petropoulos_tpm", "petropoulos_design"),
        ("tyser_E-MTAB-9388", "tyser_tpm", "tyser_design"),
    ]
    rows, provenance = [], []
    for source, archive_role, design_role in sources:
        archive_path, design_path = reference_input(archive_role), reference_input(design_role)
        for role, path in [(archive_role, archive_path), (design_role, design_path)]:
            provenance.append(
                {"role": role, "path": str(path), "sha256": sha256(path), "source_id": source}
            )
        with zipfile.ZipFile(archive_path) as archive:
            name = next(n for n in archive.namelist() if n.endswith(".mtx"))
            genes = pd.read_csv(
                archive.open(name + "_rows"), sep="\t", header=None, names=["gene_id", "label"]
            )
            cells = pd.read_csv(archive.open(name + "_cols"), sep="\t", header=None)[0]
            matrix = mmread(archive.open(name)).tocsr()
        if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
            raise ValueError("Developmental TPM must be finite and nonnegative")
        identity = match_reference_genes(
            genes,
            id_column="gene_id",
            symbol_column="label",
            source=source,
            source_version="Ensembl106_filtered_TPM",
        )
        decisions = resolve_measurement_rows(
            identity, equivalent=lambda a, b, values=matrix: (values[a] != values[b]).nnz == 0
        )
        write_csv(
            decisions,
            ROOT
            / "docs/audits/target-selection-20261110/rebuild/reproductive-identities"
            / (source + ".csv.gz"),
        )
        chosen = decisions[decisions.row_position.notna()].reset_index(drop=True)
        matrix = matrix[chosen.row_position.astype(int).to_numpy()]
        if not cells.is_unique:
            raise ValueError("Duplicated developmental cell IDs")
        meta = pd.read_csv(design_path, sep="\t").set_index("Assay").loc[cells]
        if source.startswith("petropoulos"):
            meta["embryo"] = meta["Sample Characteristic[individual]"]
            meta["stage"] = meta["Sample Characteristic[age]"]
            meta["lineage"] = (
                meta["Sample Characteristic[age]"]
                + " / "
                + meta["Factor Value[inferred lineage]"].fillna("unspecified lineage")
            )
        else:
            meta["embryo"] = "single-embryo"
            meta["stage"] = "Carnegie stage 7"
            meta["lineage"] = (
                "Carnegie stage 7 / "
                + meta["Sample Characteristic[inferred cell type - authors labels]"]
            )
        if meta[["embryo", "lineage"]].isna().any().any():
            raise ValueError("Missing developmental donor or context")
        positions = dict(zip(chosen.gene_id, range(len(chosen))))
        groups = list(meta.groupby(["embryo", "lineage"]).indices.items())
        for protein in facts.itertuples():
            members = json.loads(protein.member_gene_ids)
            complete = all(g in positions for g in members)
            values = (
                np.asarray(matrix[[positions[g] for g in members]].sum(axis=0)).ravel()
                if complete
                else None
            )
            for (embryo, lineage), indices in groups:
                stage = meta.iloc[indices].stage.unique()
                if len(stage) != 1:
                    raise ValueError("Developmental lineage spans multiple ages")
                stage = stage[0]
                n_positive = int((values[indices] >= 1).sum()) if complete else 0
                donors_required = developmental_donor_requirement(len(indices), n_positive)
                supported = complete and donors_required > 0
                rows.append(
                    {
                        "entity_id": protein.entity_id,
                        "source_id": source,
                        "sample_id": str(embryo) + " / " + lineage,
                        "participant_id": source + ":" + str(embryo),
                        "tissue": "embryo",
                        "compartment": "embryo",
                        "source_context": lineage,
                        "cell_type": lineage,
                        "scope": "developmental",
                        "statistic": "mean native TPM within one embryo and lineage",
                        "observation_kind": "embryo_lineage_mean",
                        "unit": "TPM",
                        "value": float(values[indices].mean()) if complete else np.nan,
                        "qc_pass": bool(supported),
                        "identity_resolved": complete,
                        "cells": len(indices),
                        "positive_cells": n_positive,
                        "developmental_stage": stage,
                        "required_embryos": donors_required,
                    }
                )
    return pd.DataFrame(rows), provenance

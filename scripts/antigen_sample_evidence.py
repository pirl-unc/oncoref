"""Pinned individual tissue RNA observations for the symmetric antigen screen."""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from cta_report_common import sha256, write_csv

ROOT = Path(__file__).resolve().parents[1]
HPA_ARCHIVE = ROOT / "tmp/harmonize-antigens/transcript_rna_tissue.v23.tsv.zip"
HPA_SOURCE = {
    "source_id": "hpa_tissue_samples:v23",
    "url": "https://v23.proteinatlas.org/download/transcript_rna_tissue.tsv.zip",
    "sha256": "4f4790de54380f4d47800b16efaed5a0d6bc364a90ddfdd0582f8f2d13e4f9e1",
    "reference": "HPA v23, Ensembl 109",
    "unit": "TPM",
    "statistic": "individual tissue sample; summed transcript TPM",
    "participant_identity": "Not provided; sample identifiers are not assumed to be donor IDs",
}


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


def tissue_sample_genes(archive=HPA_ARCHIVE):
    """Verify the provider archive and extract gene TPM without renormalizing it."""
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
    genes = tissue_sample_genes() if genes is None else genes.copy()
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
            rows.append(
                {
                    "entity_id": protein.entity_id,
                    "source_id": HPA_SOURCE["source_id"],
                    "sample_id": sample,
                    "participant_id": "",
                    "tissue": tissue,
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


def developmental_samples(facts):
    """TPM means within one embryo and lineage, never across biological donors.

    A positive observation needs >=20 measured cells, >=5 cells with >=1 TPM,
    and >=10% positive cells. Cell counts are support, not independent witnesses.
    """
    from scipy.io import mmread

    from oncoref import reference_data
    from oncoref.gene_matching import match_reference_genes, resolve_measurement_rows

    directory = ROOT / "tmp/cta-source-consolidation"
    tyser = ROOT / "tmp/cta-quantitative"
    receipt = json.loads(
        (ROOT / "docs/audits/cta-selection-20261028/developmental-source.json").read_text()
    )
    sources = [
        (
            "petropoulos_E-MTAB-3929",
            directory / "petropoulos-tpm.zip",
            directory / "petropoulos-design.tsv",
            reference_data.REFERENCE_SOURCES["petropoulos_embryonic_tpm"]["pins"][
                "E-MTAB-3929-20261028"
            ]["sha256"],
            reference_data.REFERENCE_SOURCES["petropoulos_embryonic_design"]["pins"][
                "E-MTAB-3929-20261028"
            ]["sha256"],
        ),
        (
            "tyser_E-MTAB-9388",
            tyser / "atlas_tpm",
            tyser / "atlas_design",
            receipt["files"]["atlas_tpm"]["sha256"],
            receipt["files"]["atlas_design"]["sha256"],
        ),
    ]
    rows, provenance = [], []
    for source, archive_path, design_path, archive_sha, design_sha in sources:
        for path, digest in [(archive_path, archive_sha), (design_path, design_sha)]:
            if sha256(path) != digest:
                raise ValueError(f"Changed developmental input: {path}")
            provenance.append(
                {"path": str(path.relative_to(ROOT)), "sha256": digest, "source_id": source}
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
            meta["lineage"] = (
                meta["Sample Characteristic[age]"]
                + " / "
                + meta["Factor Value[inferred lineage]"].fillna("unspecified lineage")
            )
        else:
            meta["embryo"] = "single-embryo"
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
                n_positive = int((values[indices] >= 1).sum()) if complete else 0
                supported = (
                    complete
                    and len(indices) >= 20
                    and n_positive >= 5
                    and n_positive / len(indices) >= 0.1
                )
                rows.append(
                    {
                        "entity_id": protein.entity_id,
                        "source_id": source,
                        "sample_id": str(embryo) + " / " + lineage,
                        "participant_id": source + ":" + str(embryo),
                        "tissue": "embryo",
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
                    }
                )
    return pd.DataFrame(rows), provenance

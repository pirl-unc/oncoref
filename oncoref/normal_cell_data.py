# Licensed under the Apache License, Version 2.0.
"""Versioned normal-cell evidence, alongside the original HPA reference tables.

``sources``, ``table``, ``iter_tables`` and ``observations`` share one registry
and cache. Native counts, TPM, nTPM and antibody staining are never conflated.
Large matrices/XML have streaming readers; derived tables are checksum-keyed.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from filelock import FileLock

from . import reference_data
from .gene_ids import canonical_gene_space, resolve_symbol
from .load_dataset import get_data

ADAPTER_VERSION = "1"
ORIGINAL_SOURCES = ("hpa_rna_consensus", "hpa_normal_tissue", "hpa_single_cell")
OBSERVATION_COLUMNS = [
    "gene_id",
    "gene_name",
    "source",
    "version",
    "taxonomy_id",
    "tissue",
    "cell_type",
    "cluster",
    "donor",
    "sample",
    "state",
    "stage",
    "modality",
    "unit",
    "value",
    "level",
    "reliability",
    "antibody",
    "read_count",
]
IHC_COLUMNS = (
    "Gene",
    "Gene name",
    "Tissue",
    "Cell type",
    "Level",
    "Antibody",
    "Profile reliability",
    "Profile notes",
    "Antibody validation",
    "RNA consistency",
    "Literature consistency",
    "Antigen sequence",
    "Donor IDs",
    "Image URLs",
    "Staining summary",
    "Source URL",
)


def publication_evidence():
    """Reviewed human findings, citations and limitations; not gene allow/deny overrides."""
    return get_data("normal-cell-publication-evidence").copy()


def sources() -> pd.DataFrame:
    """Available datasets and their exact version, URL, checksum and native scope."""
    names = [*ORIGINAL_SOURCES, *reference_data._normal_manifest["artifacts"]]
    return pd.DataFrame(
        [
            {
                "source": n,
                "version": reference_data.resolve_version(n),
                "description": reference_data.REFERENCE_SOURCES[n]["description"],
                "url": reference_data.REFERENCE_SOURCES[n]["urls"][
                    reference_data.resolve_version(n)
                ],
                "taxonomy_id": 9606,
                "sha256": reference_data.provenance(n)["sha256"],
                "publication_url": reference_data.REFERENCE_SOURCES[n].get("publication_url", ""),
                "cached": reference_data.is_cached(n),
            }
            for n in names
        ]
    )


def _gene_filter(frame, genes):
    if genes is None:
        return frame
    selected = set(genes)
    if "Feature_name" in frame:
        ids = frame.Feature_name.str.split(".").str[0]
        symbols = frame.Additional_annotations.str.split(";").str[0]
        return frame.loc[ids.isin(selected) | symbols.isin(selected)]
    keys = [c for c in ("Gene", "Gene name", "ENSG ID", "Ensembl", "gene_name") if c in frame]
    if not keys:
        # Corneal counts are native wide matrices, indexed by the original symbol.
        return frame.loc[frame.index.isin(selected)]
    mask = pd.Series(False, index=frame.index)
    for key in keys:
        mask |= frame[key].isin(selected)
    return frame.loc[mask]


def _is_count_matrix(source):
    return source.startswith("cornea_collin_") and (
        source.endswith("_counts") or "_additional" in source
    )


def iter_tables(source: str, *, genes=None, version=None, chunksize=512):
    """Yield native table chunks; genes may be Ensembl IDs or source gene symbols.

    Corneal matrices use symbols as their row index and cell barcodes as columns.
    Unknown symbols are not substituted with paralogs. XML yields antibody-level
    tissue/cell aggregates, NOT invented patient-specific staining measurements.
    """
    if source not in set(sources().source):
        raise ValueError(f"Not a normal-cell evidence source: {source}")
    selected = None if genes is None else set(genes)
    if source == "hpa_antibody_xml":
        p = _antibody_cache(version)
        profile_source = reference_data.ensure(
            "hpa_gene_summary", reference_data.resolve_version(source, version)
        )
        profiles = pd.read_csv(
            profile_source, sep="\t", usecols=["Ensembl", "Reliability (IH)"], keep_default_na=False
        ).set_index("Ensembl")["Reliability (IH)"]
        for batch in pq.ParquetFile(p).iter_batches(batch_size=chunksize):
            frame = _gene_filter(batch.to_pandas(), selected)
            if len(frame):
                yield _ihc_with_summary(frame, profiles)
        return
    p = reference_data.ensure(source, version)
    counts = _is_count_matrix(source)
    for frame in pd.read_csv(
        p,
        sep="," if counts or source.startswith("dice_") else "\t",
        index_col=0 if counts else None,
        chunksize=chunksize,
        keep_default_na=False,
    ):
        frame = _gene_filter(frame, selected)
        if len(frame):
            yield frame


def table(source: str, *, genes=None, version=None) -> pd.DataFrame:
    """Read native data; use ``iter_tables`` for an entire large cell matrix."""
    parts = list(
        iter_tables(
            source,
            genes=genes,
            version=version,
            chunksize=50000 if not _is_count_matrix(source) else 128,
        )
    )
    if not parts:
        if source == "hpa_antibody_xml":
            return pd.DataFrame(columns=IHC_COLUMNS)
        counts = _is_count_matrix(source)
        return pd.read_csv(
            reference_data.ensure(source, version),
            sep="," if counts or source.startswith("dice_") else "\t",
            index_col=0 if counts else None,
            nrows=0,
        )
    result = pd.concat(parts, ignore_index=not _is_count_matrix(source))
    result.attrs["source"] = reference_data.provenance(source, version)
    return result


def _cache_path(source, inputs, suffix):
    fingerprint = hashlib.sha256(
        json.dumps({"adapter": ADAPTER_VERSION, "inputs": inputs}, sort_keys=True).encode()
    ).hexdigest()[:20]
    return source.parent / "derived" / f"{suffix}-{fingerprint}.parquet"


def _ihc_with_summary(frame, profiles):
    """The XML omits some uncertain curated profiles; the summary preserves them."""
    frame = frame.copy()
    frame["XML profile reliability"] = frame["Profile reliability"]
    supplement = frame.Gene.map(profiles).fillna("")
    missing = frame["Profile reliability"].isin(["Missing", ""])
    fill = missing & supplement.ne("")
    frame.loc[fill, "Profile reliability"] = supplement[fill]
    frame["Profile reliability source"] = np.where(fill, "hpa_gene_summary", "hpa_antibody_xml")
    return frame


def _entry_ihc(entry):
    identifier = entry.find('identifier[@db="Ensembl"]')
    if identifier is None:
        return []
    profile = entry.find('tissueExpression[@technology="IHC"][@assayType="tissue"]')
    reliability = profile.find('verification[@type="reliability"]') if profile is not None else None
    profile_score = (reliability.text or "").title() if reliability is not None else "Missing"
    notes = reliability.get("description", "") if reliability is not None else ""
    rows = []
    for antibody in entry.findall("antibody"):
        for assay in antibody.findall('tissueExpression[@technology="IHC"][@assayType="tissue"]'):
            for data in assay.findall("data"):
                donors = sorted({p.findtext("patientId", "") for p in data.findall("patient")})
                images = sorted({i.text for i in data.findall(".//imageUrl") if i.text})
                for cell in data.findall("tissueCell"):
                    rows.append(
                        {
                            "Gene": identifier.get("id", ""),
                            "Gene name": entry.findtext("name", ""),
                            "Tissue": data.findtext("tissue", ""),
                            "Cell type": cell.findtext("cellType", ""),
                            "Level": cell.findtext('level[@type="staining"]', "").capitalize(),
                            "Antibody": antibody.get("id", ""),
                            "Profile reliability": profile_score,
                            "Profile notes": notes,
                            "Antibody validation": assay.findtext(
                                'verification[@type="validation"]', ""
                            ).title(),
                            "RNA consistency": assay.findtext(
                                'validation[@type="RNAConsistency"]', ""
                            ),
                            "Literature consistency": assay.findtext(
                                'validation[@type="literatureConformity"]', ""
                            ),
                            "Antigen sequence": antibody.findtext("antigenSequence", ""),
                            "Donor IDs": ";".join(donors),
                            "Image URLs": ";".join(images),
                            "Staining summary": assay.findtext('summary[@type="tissue"]', ""),
                            "Source URL": entry.get("url", "").replace("http:", "https:")
                            + "/tissue/primary+data",
                        }
                    )
    return rows


def _antibody_cache(version=None):
    source = "hpa_antibody_xml"
    p = reference_data.ensure(source, version)
    provenance = reference_data.provenance(source, version)
    dest = _cache_path(p, [provenance["sha256"]], "primary-ihc")
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    with FileLock(str(dest) + ".lock"):
        if dest.exists():
            return dest
        schema = pa.schema([(c, pa.string()) for c in IHC_COLUMNS])
        fd, filename = tempfile.mkstemp(dir=dest.parent, suffix=".parquet.part")
        os.close(fd)
        try:
            with (
                pq.ParquetWriter(filename, schema, compression="zstd") as writer,
                gzip.open(p, "rb") as stream,
            ):
                context = ET.iterparse(stream, events=("start", "end"))
                _, root = next(context)
                rows = []
                for event, element in context:
                    if event == "end" and element.tag == "entry":
                        rows.extend(_entry_ihc(element))
                        element.clear()
                        root.clear()
                        if len(rows) >= 10000:
                            writer.write_table(pa.Table.from_pylist(rows, schema=schema))
                            rows.clear()
                if rows:
                    writer.write_table(pa.Table.from_pylist(rows, schema=schema))
            os.replace(filename, dest)
        finally:
            Path(filename).unlink(missing_ok=True)
    return dest


def observations(source: str, *, genes=None, version=None) -> pd.DataFrame:
    """Uniform long evidence with explicit units, stage, cell context and provenance.

    Metadata tables are accessed with ``table``. Corneal matrices return
    cluster/sample summaries (CP10K, not TPM), with raw detection/count columns.
    """
    if source == "cornea_collin_adult_counts":
        return cornea_summary(genes=genes)
    if source.startswith("dice_"):
        return dice_observations(source, genes=genes, version=version)
    native = table(source, genes=genes, version=version)
    if source.endswith("_cells") or source in ("hpa_gene_summary", "hpa_single_cell_clusters"):
        raise ValueError("Metadata source: use table(), not observations()")
    if _is_count_matrix(source):
        raise ValueError(
            "Use iter_tables() for developmental/culture count matrices; do not pool with adult normal tissue"
        )
    out = pd.DataFrame(index=native.index)
    out["gene_id"] = native.get("Gene", native.get("ENSG ID", ""))
    out["gene_name"] = native.get("Gene name", "")
    out["source"] = source
    out["version"] = reference_data.resolve_version(source, version)
    out["taxonomy_id"] = 9606
    out["tissue"] = native.get("Tissue", "blood" if "immune" in source else "")
    out["cell_type"] = native.get("Cell type", native.get("Immune cell", ""))
    out["cluster"] = native.get("Cluster", "")
    out["donor"] = native.get("Donor", "")
    out["sample"] = native.get("Sample ID", "")
    out["state"] = np.where(
        out.cell_type.str.contains("stimulated|activated", case=False), "activated", "unspecified"
    )
    out["stage"] = "adult"
    out["modality"] = "ihc" if "Level" in native else "rna"
    unit = "nTPM" if "nTPM" in native else "TPM" if "TPM" in native else "staining"
    out["unit"] = unit
    out["value"] = pd.to_numeric(native[unit], errors="coerce") if unit != "staining" else np.nan
    out["level"] = native.get("Level", "")
    out["reliability"] = native.get("Reliability", native.get("Profile reliability", ""))
    out["antibody"] = native.get("Antibody", "")
    out["read_count"] = pd.to_numeric(native.get("Read count", np.nan), errors="coerce")
    out.attrs["source"] = reference_data.provenance(source, version)
    return out


def dice_observations(source, *, genes=None, version=None):
    """DICE's per-gene specimen distributions in TPM, with no invented donor pairing.

    Public columns are rank/replicate slots and are NOT stable donor IDs across
    genes or cell types. They cannot support paired or cross-gene donor analyses.
    """
    spec = reference_data.REFERENCE_SOURCES[source]
    if not source.startswith("dice_"):
        raise ValueError("Expected an original DICE matrix")
    native = table(source, genes=genes, version=version)
    if native.empty:
        return pd.DataFrame(columns=[*OBSERVATION_COLUMNS, "replicate_slot"])
    native["gene_id"] = native.Feature_name.str.split(".").str[0]
    native["gene_name"] = native.Additional_annotations.str.split(";").str[0]
    slots = [c for c in native if c.isdigit()]
    out = native.melt(
        id_vars=["gene_id", "gene_name"],
        value_vars=slots,
        var_name="replicate_slot",
        value_name="value",
    )
    out["value"] = pd.to_numeric(out.value, errors="raise")
    if not np.isfinite(out.value).all() or out.value.lt(0).any():
        raise ValueError("Invalid DICE TPM")
    for key, value in {
        "source": source,
        "version": reference_data.resolve_version(source, version),
        "taxonomy_id": 9606,
        "tissue": "blood",
        "cell_type": spec["cell_type"],
        "state": spec["state"],
        "stage": "adult",
        "unit": "TPM",
        "modality": "rna",
        "donor": "",
        "sample": "",
        "cluster": "",
        "level": "",
        "antibody": "",
        "reliability": "",
    }.items():
        out[key] = value
    out["read_count"] = np.nan
    out.attrs["source"] = reference_data.provenance(source, version)
    out.attrs["donor_pairing"] = "Unavailable: columns are not stable donor IDs"
    return out


def cornea_cells() -> pd.DataFrame:
    """Published adult normal cell types; barcode suffixes are sample groups, not donor IDs."""
    frame = table("cornea_collin_adultcornea_cells").rename(
        columns={"cellName": "cell", "cluster": "cell_type"}
    )
    if frame.cell.duplicated().any():
        raise ValueError("Duplicate corneal cell annotation")
    frame["sample_group"] = frame.cell.str.rsplit("_", n=1).str[-1]
    frame["donor"] = ""  # The published annotation does not link donor identities.
    frame["tissue"] = "cornea/conjunctiva"
    frame["stage"] = "adult"
    return frame


def _count_rows(path):
    with gzip.open(path, "rt") as stream:
        cells = next(csv.reader([stream.readline()]))[1:]
        yield cells
        for line in stream:
            name, values = line.split(",", 1)
            counts = np.fromstring(values, sep=",", dtype=np.float64)
            if len(counts) != len(cells) or (counts < 0).any() or not np.isfinite(counts).all():
                raise ValueError("Invalid corneal count row")
            yield name.strip('"'), counts


def _read_cornea_summary(path, genes):
    """Resolve explicit source symbols against the shipped canonical identity map."""
    names = pq.read_table(path, columns=["gene_name"]).column(0).unique().to_pylist()
    space = canonical_gene_space()
    unique = space.loc[space.symbol.ne("") & ~space.symbol.duplicated(keep=False)]
    index = unique.set_index("symbol").ensembl_gene_id
    identities = pd.DataFrame({"gene_name": names})
    identities["canonical_symbol"] = identities.gene_name.map(resolve_symbol)
    identities["gene_id"] = identities.canonical_symbol.map(index).fillna("")
    ambiguous = identities.gene_id.ne("") & identities.gene_id.duplicated(keep=False)
    identities["gene_mapping"] = np.where(identities.gene_id.eq(""), "Unmapped", "Mapped")
    identities.loc[ambiguous, "gene_mapping"] = "Ambiguous"
    identities.loc[ambiguous, "gene_id"] = ""
    selected = identities
    if genes is not None:
        selected = identities.loc[
            identities.gene_name.isin(genes)
            | identities.gene_id.isin(genes)
            | identities.canonical_symbol.isin([resolve_symbol(g) for g in genes])
        ]
    if selected.empty:
        frame = pa.Table.from_batches([], schema=pq.ParquetFile(path).schema_arrow).to_pandas()
    else:
        frame = pd.read_parquet(path, filters=[("gene_name", "in", selected.gene_name.tolist())])
    frame = frame.drop(columns="gene_id").merge(
        selected, on="gene_name", how="left", validate="many_to_one"
    )
    return frame


def cornea_summary(*, genes=None) -> pd.DataFrame:
    """Adult normal cornea/conjunctiva RNA per cell type and barcode sample group.

    CP10K is the mean of count/whole-cell-library-size * 10,000, including zero
    cells. It is not TPM. Every source gene contributes to each library size,
    even when callers request a subset of genes. Zero is measured; absent genes
    have no row. Developmental and diseased/cultured samples are never pooled.
    """
    source = "cornea_collin_adult_counts"
    p = reference_data.ensure(source)
    meta = cornea_cells()
    keys = [
        reference_data.provenance(n)["sha256"] for n in (source, "cornea_collin_adultcornea_cells")
    ]
    dest = _cache_path(p, keys, "adult-cell-summary")
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        result = _read_cornea_summary(dest, genes)
        result.attrs["source"] = reference_data.provenance(source)
        return result
    with FileLock(str(dest) + ".lock"):
        if not dest.exists():
            rows = _count_rows(p)
            cells = next(rows)
            if set(cells) != set(meta.cell) or len(cells) != len(meta):
                raise ValueError("Corneal count cells do not match the published adult annotations")
            meta = meta.set_index("cell").loc[cells]
            libraries = np.zeros(len(cells), dtype=float)
            for _, counts in rows:
                libraries += counts
            if (libraries <= 0).any():
                raise ValueError("Empty corneal cell library")
            groups = list(meta.groupby(["cell_type", "sample_group"], sort=True).indices.items())
            counts_iter = _count_rows(p)
            next(counts_iter)
            fd, staged = tempfile.mkstemp(dir=dest.parent, suffix=".parquet.part")
            os.close(fd)
            writer = None
            try:
                result = []
                for gene, counts in counts_iter:
                    normalized = counts / libraries * 10000
                    for (cell_type, sample), idx in groups:
                        result.append(
                            {
                                "gene_id": "",
                                "gene_name": gene,
                                "source": source,
                                "version": reference_data.resolve_version(source),
                                "taxonomy_id": 9606,
                                "tissue": "cornea/conjunctiva",
                                "cell_type": cell_type,
                                "cluster": cell_type,
                                "sample": sample,
                                "donor": "",
                                "state": "normal",
                                "stage": "adult",
                                "modality": "rna",
                                "unit": "CP10K",
                                "value": float(normalized[idx].mean()),
                                "cells": len(idx),
                                "positive_cells": int((counts[idx] > 0).sum()),
                                "read_count": int(counts[idx].sum()),
                                "level": "",
                                "reliability": "",
                                "antibody": "",
                            }
                        )
                    if len(result) >= 10000:
                        batch = pa.Table.from_pylist(result)
                        if writer is None:
                            writer = pq.ParquetWriter(staged, batch.schema, compression="zstd")
                        writer.write_table(batch)
                        result.clear()
                if result:
                    batch = pa.Table.from_pylist(result)
                    if writer is None:
                        writer = pq.ParquetWriter(staged, batch.schema, compression="zstd")
                    writer.write_table(batch)
                if writer is None:
                    raise ValueError("Empty corneal counts")
                writer.close()
                writer = None
                os.replace(staged, dest)
            finally:
                if writer is not None:
                    writer.close()
                Path(staged).unlink(missing_ok=True)
    result = _read_cornea_summary(dest, genes)
    result.attrs["source"] = reference_data.provenance(source)
    return result

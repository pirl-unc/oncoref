# Licensed under the Apache License, Version 2.0.
"""Native expression matrices with species contracts and auditable sample metadata."""

from __future__ import annotations

import copy
import hashlib
import json
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from ._table_io import read_csv
from .expression_contract import ExpressionContract, validate_expression_join
from .expression_records import RECORD_KINDS, ExpressionRecords
from .expression_registry import (
    expression_source,
    expression_source_registry_entries,
    expression_source_registry_path,
)
from .reference_identity import _identity
from .species import resolve_species
from .species_assets import DEFAULT_MAX_BYTES, ensure_asset, species_asset


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _numbers(frame):
    """Blank/NA measurements stay missing; malformed or negative values fail."""
    values = frame.replace({"": np.nan, "NA": np.nan, "NaN": np.nan}).apply(
        pd.to_numeric, errors="raise"
    )
    array = values.to_numpy(dtype=float)
    if np.isinf(array).any() or (array < 0).any():
        raise ValueError("Expression values must be finite nonnegative numbers or missing")
    return values.astype(float)


class ExpressionDataset:
    """One native feature-by-library matrix; accessors return defensive copies.

    Source byte pins preserve the raw records. Metadata retains source IDs, tissue
    labels, all mapping alternatives, pooling and unresolved donor/reference facts.
    """

    def __init__(self, values, features, samples, contract, *, records=None):
        if not isinstance(contract, ExpressionContract):
            raise ValueError("ExpressionContract required")
        self._contract = contract
        if records is not None and not isinstance(records, ExpressionRecords):
            raise TypeError("records must be ExpressionRecords or None")
        self._records = records
        self._values = _numbers(values.copy())
        self._features = pd.DataFrame(
            copy.deepcopy(features.to_dict("list")), index=features.index.copy()
        )
        self._samples = pd.DataFrame(
            copy.deepcopy(samples.to_dict("list")), index=samples.index.copy()
        )
        for index in (values.index, values.columns):
            if not index.is_unique or any(not isinstance(x, str) or not x.strip() for x in index):
                raise ValueError("Unique nonempty source feature/sample identifiers required")
        if not features.index.equals(values.index) or not samples.index.equals(values.columns):
            raise ValueError("Matrix feature/sample IDs must exactly match their metadata")
        required = {
            "taxon_id",
            "library_id",
            "run_ids",
            "donor_id",
            "donor_scope",
            "pooling",
            "tissue",
        }
        if not required <= set(samples):
            raise ValueError(f"Missing sample metadata: {sorted(required - set(samples))}")
        if not samples.taxon_id.map(lambda x: type(x) is int and x == contract.taxon_id).all():
            raise ValueError("Matrix sample taxon mismatch")
        if not samples.pooling.isin(["single", "pooled", "unknown"]).all():
            raise ValueError("Explicit single/pooled/unknown sample pooling required")
        seen_runs = set()
        for runs in samples.run_ids:
            if not isinstance(runs, list) or any(
                not isinstance(run, str) or not run for run in runs
            ):
                raise ValueError("run_ids must be explicit lists (empty means unresolved)")
            if len(set(runs)) != len(runs) or seen_runs.intersection(runs):
                raise ValueError(
                    "Run reused across source libraries; reconcile lineage before summarizing"
                )
            seen_runs.update(runs)

    @property
    def contract(self):
        return self._contract

    @property
    def values(self):
        return self._values.copy()

    @property
    def features(self):
        return pd.DataFrame(
            copy.deepcopy(self._features.to_dict("list")), index=self._features.index.copy()
        )

    @property
    def samples(self):
        return pd.DataFrame(
            copy.deepcopy(self._samples.to_dict("list")), index=self._samples.index.copy()
        )

    @property
    def provenance(self):
        """Native dataset identity and processing, or None for a matrix-only snapshot."""
        return self._records.provenance if self._records is not None else None

    def native_records(self, kind):
        """Stream source records with their original schema, IDs and uncertainty.

        Kinds: samples, assays, measurements, reference_mappings, mappings,
        relations, pairings. An empty iterator does not establish biological absence.
        Historical matrix-only snapshots have no native record archive.
        """
        if kind not in RECORD_KINDS:
            raise ValueError(f"Unknown expression record kind {kind!r}")
        if self._records is not None:
            yield from self._records.records(kind)

    def tissue_annotations(self):
        """Cited current tissue-scope annotations; retain the frozen sample metadata.

        No matrix download, feature remapping or cross-species expression inference.
        Unknown source versions return unresolved mappings.
        """
        from .normal_tissues import resolve_normal_tissue

        def optional(value):
            return None if value is None or pd.isna(value) else str(value)

        records = [
            resolve_normal_tissue(
                optional(row["tissue"]),
                species=self.contract.taxon_id,
                cell_type=optional(row.get("cell_type")),
                stage=optional(row.get("stage")),
                sex=optional(row.get("sex")),
                physiological_state=optional(row.get("physiological_state")),
                source_id=self.contract.source_id,
                source_version=self.contract.source_version,
                context=self.contract.context,
            ).as_dict()
            for row in self._samples.to_dict("records")
        ]
        return pd.DataFrame(records, index=self._samples.index.copy())

    def tissue_summary(self):
        """Within-source native values; measured libraries and donor coverage remain separate.

        This does not allocate gene totals to isoforms or pool incompatible sources.
        Unknown tissues are retained. Developmental stages are grouped separately
        when supplied; embryonic and adult samples never share an organ summary.
        Pooled/unknown-pooling samples cannot supply
        independent-donor counts. Missing values never become zero observations.
        """
        rows = []
        group_columns = ["tissue", "stage"] if "stage" in self._samples else ["tissue"]
        for keys, metadata in self._samples.groupby(group_columns, dropna=False, sort=False):
            if not isinstance(keys, tuple):
                keys = (keys,)
            context = {
                col: None if pd.isna(value) else value for col, value in zip(group_columns, keys)
            }
            values = self._values.loc[:, metadata.index]
            measured = values.notna()
            individual = metadata[metadata.pooling.eq("single") & metadata.donor_id.notna()]
            # Count each known donor once per feature, regardless of technical replicates.
            donor_counts = pd.Series(0, index=values.index, dtype=int)
            for _, donor_samples in individual.groupby(["donor_scope", "donor_id"], dropna=False):
                donor_counts += measured.loc[:, donor_samples.index].any(axis=1).astype(int)
            summary = pd.DataFrame(
                {
                    "feature_id": values.index,
                    **context,
                    "native_mean": values.mean(axis=1).to_numpy(),
                    "native_max": values.max(axis=1).to_numpy(),
                    "measured_libraries": measured.sum(axis=1).to_numpy(),
                    "total_libraries": len(metadata),
                    "known_individual_donors_measured": donor_counts.to_numpy(),
                    "libraries_with_unresolved_or_pooled_donors": len(metadata) - len(individual),
                    "unit": self.contract.measurement.unit,
                }
            )
            rows.append(summary)
        result = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
        result.attrs["expression_contract"] = self.contract.as_dict()
        return result

    def save(self, destination):
        """Write a new portable artifact atomically; refuse to replace existing results."""
        destination = Path(destination)
        if destination.exists():
            raise FileExistsError(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=".oncoref-expression-", dir=destination.parent
        ) as temporary:
            staged = Path(temporary) / "artifact"
            staged.mkdir()
            self._values.to_parquet(staged / "values.parquet")
            version = 2 if self._records is not None else 1
            if version == 1:
                self._features.to_parquet(staged / "features.parquet")
                feature_file = "features.parquet"
            else:
                feature_file = "features.json"
                # Keep the index separate from arbitrary native feature columns.
                features = self._features.to_dict("tight")
                (staged / feature_file).write_text(json.dumps(features, allow_nan=False) + "\n")
            # JSON retains nested original sample facts and provenance without stringifying them.
            records = copy.deepcopy(
                self._samples.astype(object).where(self._samples.notna(), None).to_dict("records")
            )
            for sample_id, row in zip(self._samples.index, records):
                row["sample_id"] = sample_id
            (staged / "samples.json").write_text(json.dumps(records, allow_nan=False) + "\n")
            names = ["values.parquet", feature_file, "samples.json"]
            if self._records is not None:
                names.extend(self._records.write(staged))
            payload = {
                "schema": "oncoref.expression-dataset",
                "schema_version": version,
                "contract": self.contract.as_dict(),
                "files": {
                    name: {
                        "sha256": _sha256(staged / name),
                        "size_bytes": (staged / name).stat().st_size,
                    }
                    for name in names
                },
            }
            if version == 2:
                payload["provenance"] = self.provenance
            record = {**payload, "artifact_id": _identity("expression-dataset", version, payload)}
            (staged / "metadata.json").write_text(
                json.dumps(record, indent=2, sort_keys=True) + "\n"
            )
            os.rename(staged, destination)
        return destination

    @classmethod
    def load(cls, directory):
        """Offline replay validates the metadata identity and every artifact file."""
        directory = Path(directory)
        record = json.loads((directory / "metadata.json").read_text())
        version = record.get("schema_version")
        fields = {"schema", "schema_version", "contract", "files", "artifact_id"}
        if version == 2:
            fields.add("provenance")
        if (
            set(record) != fields
            or record["schema"] != "oncoref.expression-dataset"
            or type(record["schema_version"]) is not int
            or version not in (1, 2)
        ):
            raise ValueError("Unsupported expression artifact schema")
        payload = {key: value for key, value in record.items() if key != "artifact_id"}
        if record["artifact_id"] != _identity("expression-dataset", version, payload):
            raise ValueError("Expression artifact metadata digest mismatch")
        expected = {"values.parquet", "features.parquet", "samples.json"}
        if version == 2:
            expected = {"values.parquet", "features.json", "samples.json"}
            expected.update(f"{kind}.jsonl.gz" for kind in RECORD_KINDS)
        if set(record["files"]) != expected:
            raise ValueError("Unexpected expression artifact file inventory")
        for name, info in record["files"].items():
            path = directory / name
            if (
                set(info) != {"sha256", "size_bytes"}
                or path.stat().st_size != info["size_bytes"]
                or _sha256(path) != info["sha256"]
            ):
                raise ValueError(f"Expression artifact checksum mismatch: {name}")
        samples = pd.DataFrame(json.loads((directory / "samples.json").read_text())).set_index(
            "sample_id"
        )
        return cls(
            pd.read_parquet(directory / "values.parquet"),
            pd.read_parquet(directory / "features.parquet")
            if version == 1
            else pd.DataFrame.from_dict(
                json.loads((directory / "features.json").read_text()), orient="tight"
            ),
            samples,
            ExpressionContract.from_dict(record["contract"]),
            records=ExpressionRecords.from_directory(directory, record["provenance"])
            if version == 2
            else None,
        )


PROTEIN_SAMPLE_COLUMNS = (
    "Run",
    "Channel",
    "BioReplicate",
    "Condition",
    "TechRepMixture",
    "Mixture",
    "sample.type",
)


def protein_sample_id(record):
    return _identity("protein-sample", 1, {key: record[key] for key in PROTEIN_SAMPLE_COLUMNS})


def read_species_expression(
    source_id,
    *,
    species,
    offline=False,
    cache_root=None,
    local_overrides=None,
    max_bytes=DEFAULT_MAX_BYTES,
):
    """Load a registered source in its native space; never apply human canonicalization."""
    selected = resolve_species(species)
    source = expression_source(source_id, species=selected)
    if source is None:
        raise KeyError(f"Unknown {selected.name} expression source {source_id!r}")
    entry = next(
        row
        for row in expression_source_registry_entries(species=selected)
        if row["id"] == source_id
    )
    adapter = entry.get("adapter")
    if adapter not in {"bgee", "barkbase", "orfdata"}:
        raise ValueError(
            f"{selected.name}/{source_id} is unavailable ({source.acquisition_status or 'reader_not_registered'}): {source.notes or 'no native reader registered'}"
        )
    contract = source.expression_contract
    expected = {
        "bgee": ("gene", "TPM"),
        "barkbase": ("gene", "TMM-CPM"),
    }.get(adapter, (contract.measurement.feature_level, contract.measurement.unit))
    if (
        (contract.measurement.feature_level, contract.measurement.unit) != expected
        or contract.measurement.scale != "linear"
        or contract.measurement.is_proxy is not False
    ):
        raise ValueError("Native reader and declared measurement semantics disagree")
    filename = entry["sample_metadata"]
    if Path(filename).name != filename:
        raise ValueError("Sample metadata must be a package-local filename")
    metadata_path = expression_source_registry_path(species=selected).parent / filename
    if _sha256(metadata_path) != contract.content_provenance["samples"]:
        raise ValueError("Sample metadata checksum mismatch")
    metadata = json.loads(metadata_path.read_text())
    if (
        set(metadata) != {"schema", "schema_version", "taxon_id", "source_id", "samples"}
        or metadata["schema"] != "oncoref.expression-samples"
        or type(metadata["schema_version"]) is not int
        or metadata["schema_version"] != 1
        or metadata["taxon_id"] != selected.taxon_id
        or metadata["source_id"] != source_id
    ):
        raise ValueError("Expression sample metadata source/schema mismatch")
    samples = pd.DataFrame(metadata["samples"]).set_index("sample_id")
    if adapter == "orfdata":
        from .managed_expression import read_managed_expression

        return read_managed_expression(
            entry,
            samples,
            contract,
            offline=offline,
            cache_root=cache_root,
            local_overrides=local_overrides,
            max_bytes=max_bytes,
        )
    asset = species_asset(entry["matrix_asset"], species=selected)
    if asset.sha256 != contract.content_provenance["matrix"]:
        raise ValueError("Matrix asset and expression contract disagree")
    path = ensure_asset(
        asset,
        offline=offline,
        cache_root=cache_root,
        local_overrides=local_overrides,
        max_bytes=max_bytes,
    )
    if adapter == "bgee":
        chunks = read_csv(
            path,
            sep="\t",
            usecols=["Experiment ID", "Library ID", "Gene ID", "TPM"],
            dtype=str,
            keep_default_na=False,
            chunksize=100000,
        )
        frame = pd.concat(chunks, ignore_index=True)
        if set(frame["Experiment ID"]) != {source.accession}:
            raise ValueError("Bgee matrix contains an unexpected study")
        if frame.duplicated(["Gene ID", "Library ID"]).any():
            raise ValueError("Duplicate Bgee gene/library measurements")
        frame["TPM"] = _numbers(frame[["TPM"]])["TPM"]
        values = frame.pivot(index="Gene ID", columns="Library ID", values="TPM")
        features = pd.DataFrame({"source_feature_id": values.index}, index=values.index)
    else:
        frame = read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        required = {"StringTieID", "EnsemblGeneID", "GeneSymbol"}
        if not required <= set(frame):
            raise ValueError("BarkBase native feature columns are missing")
        frame = frame.set_index("StringTieID")
        # StringTie identity is primary: retain unmapped and multi-mapped feature rows.
        features = frame[["EnsemblGeneID", "GeneSymbol"]].copy()
        features.insert(0, "source_feature_id", features.index)
        values = _numbers(frame.drop(columns=["EnsemblGeneID", "GeneSymbol"]))
    if set(samples.index) != set(values.columns):
        raise ValueError("Quantitative libraries and sample metadata disagree")
    samples = samples.loc[values.columns]
    values.index.name = features.index.name = "feature_id"
    values.columns.name = samples.index.name = "sample_id"
    return ExpressionDataset(values, features, samples, contract)


def combine_expression(datasets):
    """Concatenate compatible native libraries without zero filling or double counting.

    Unknown run lineage cannot establish independent libraries. Source contracts
    must already carry verified reference identities; this is not a crosswalk.
    """
    datasets = tuple(datasets)
    if not datasets:
        raise ValueError("At least one expression dataset required")
    first = datasets[0]
    seen = set()
    for dataset in datasets:
        validate_expression_join(first.contract, dataset.contract, combine_values=True)
        for runs in dataset._samples.run_ids:
            if not runs:
                raise ValueError("Cannot combine libraries with unresolved run lineage")
            if seen.intersection(runs):
                raise ValueError("Shared runs would count reused libraries more than once")
            seen.update(runs)
    # Returning component contracts retains source provenance for all contributors.
    values = pd.concat([d._values for d in datasets], axis=1, verify_integrity=True)
    return values, tuple(d.contract.as_dict() for d in datasets)


def expression_lineage_overlaps(species):
    """Source inventory overlap evidence; independent-donor counts remain unresolved."""
    path = expression_source_registry_path(species=species).parent / "lineage_overlaps.json"
    return json.loads(path.read_text()) if path.exists() else None

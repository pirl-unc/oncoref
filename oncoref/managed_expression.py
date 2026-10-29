# Licensed under the Apache License, Version 2.0.
"""Adapt pinned ORFdata releases through its public, species-neutral interface."""

from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .expression_records import ExpressionRecords
from .expression_registry import expression_source_registry_entries, expression_source_registry_path


@dataclass(frozen=True)
class ManagedExpressionAsset:
    """One optional dataset exposed through the existing species data manager."""

    id: str
    taxon_id: int
    source_version: str
    spec: dict
    group: str = "expression"

    @property
    def size_bytes(self):
        download = self.spec["download"]
        return transfer_size(download)

    @property
    def sha256(self):
        return self.spec["download"]["assets"]["dataset.sqlite"]["sha256"]

    @property
    def description(self):
        return self.spec["manifest"]["title"] + "; managed ORFdata release"

    def as_dict(self):
        return {
            "id": self.id,
            "taxon_id": self.taxon_id,
            "group": self.group,
            "source_version": self.source_version,
            "provider": "orfdata",
            "size_bytes": self.size_bytes,
            "installed_size_bytes": self.spec["download"]["assets"]["dataset.sqlite"]["size"],
            "sha256": self.sha256,
            "manifest": self.spec["manifest"],
            "download": self.spec["download"],
        }


def managed_assets(species):
    result = []
    for entry in expression_source_registry_entries(species=species):
        if entry.get("adapter") == "orfdata":
            spec = specification(entry)
            result.append(
                ManagedExpressionAsset(
                    entry["id"], entry["taxon_id"], spec["manifest"]["version"], spec
                )
            )
        for profile in entry.get("evidence_sources", []):
            spec = evidence_specification(profile, entry["taxon_id"])
            result.append(
                ManagedExpressionAsset(
                    profile["id"], entry["taxon_id"], spec["manifest"]["version"], spec
                )
            )
    return tuple(result)


def inspect_managed(asset, *, cache_root=None, local_overrides=None, verify_files=True):
    from .species_assets import _overrides
    from .species_expression import _sha256

    spec = asset.spec
    name = spec["manifest"]["dataset_id"]
    path = _overrides(local_overrides).get(f"{asset.taxon_id}:{name}")
    status, error = "missing", None
    try:
        registry, name, version = _validated_registry(spec, cache_root)
        if path is not None:
            path = Path(path).expanduser()
            pin = spec["download"]["assets"]["dataset.sqlite"]
            if path.is_file():
                status = (
                    "available"
                    if path.stat().st_size == pin["size"]
                    and (not verify_files or _sha256(path) == pin["sha256"])
                    else "corrupt"
                )
        else:
            state = registry.inspect(name, version, verify_files=verify_files)
            status = state.status
            path = state.files["dataset.sqlite"].path if status == "available" else state.path
            error = str(state.error) if state.error else None
    except ImportError as exception:
        status, error = "dependency_missing", str(exception)
    return {
        "asset_id": asset.id,
        "taxon_id": asset.taxon_id,
        "group": asset.group,
        "provider": "orfdata",
        "path": str(path) if path is not None else None,
        "status": status,
        "error": error,
        "size_bytes": asset.size_bytes,
        "installed_size_bytes": spec["download"]["assets"]["dataset.sqlite"]["size"],
        "sha256": asset.sha256,
    }


def _digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    ).hexdigest()


def transfer_size(download):
    if "archive" in download:
        return download["archive"]["size"]
    if any(value.get("decompress") for value in download["assets"].values()):
        size = download["distribution"].get("download_size")
        if type(size) is not int or size < 1:
            raise ValueError("Compressed quantitative release needs a pinned download size")
        return size
    return sum(value["size"] for value in download["assets"].values())


def specification(entry):
    """Read a bundled pin without importing the optional quantitative dependency."""
    return _specification(
        entry["managed_dataset"],
        entry["taxon_id"],
        {"manifest", "download", "selection", "source_roles", "feature_fields"},
    )


def evidence_specification(profile, species):
    return _specification(profile["managed_dataset"], species, {"manifest", "download"})


def _specification(filename, species, keys):
    if Path(filename).name != filename:
        raise ValueError("Managed dataset specification must be package-local")
    path = expression_source_registry_path(species=species).parent / filename
    spec = json.loads(path.read_text())
    if set(spec) != keys:
        raise ValueError("Unsupported managed expression specification")
    if spec["manifest"].get("schema") != "orfdata.measurements.v1":
        raise ValueError("Unsupported quantitative manifest schema")
    return spec


def _registry(cache_root):
    if sys.version_info < (3, 11):
        raise ImportError(
            "Managed quantitative data requires Python 3.11+ and oncoref[quantitative]"
        )
    try:
        from orfdata import DatasetRegistry
    except ImportError as error:
        raise ImportError(
            "Install oncoref[quantitative] to read managed quantitative data"
        ) from error
    return DatasetRegistry(cache_root=cache_root)


def _validated_registry(spec, cache_root):
    registry = _registry(cache_root)
    manifest = spec["manifest"]
    name, version = manifest["dataset_id"], manifest["version"]
    if (
        registry.info(name, version).to_dict() != manifest
        or registry.download_metadata(name, version) != spec["download"]
        or not registry.published(name, version)
    ):
        raise ValueError("ORFdata release manifest/schema/artifact pins disagree with OncoRef")
    return registry, name, version


def open_managed(
    spec, *, taxon_id, offline=False, cache_root=None, local_overrides=None, max_bytes
):
    """Open verified managed bytes; download only the declared published artifact."""
    from .species_assets import _overrides
    from .species_expression import _sha256

    if type(max_bytes) is not int or max_bytes < 0:
        raise ValueError("max_bytes must be a nonnegative integer transfer budget")
    registry, name, version = _validated_registry(spec, cache_root)
    overrides = _overrides(local_overrides)
    key = f"{taxon_id}:{name}"
    if key in overrides:
        from orfdata import QuantitativeDataset

        path = Path(overrides[key]).expanduser()
        pin = spec["download"]["assets"]["dataset.sqlite"]
        if path.stat().st_size != pin["size"] or _sha256(path) != pin["sha256"]:
            raise ValueError("Managed expression override checksum mismatch")
        dataset = QuantitativeDataset(path)
    else:
        try:
            dataset = registry.open(name, version)
        except FileNotFoundError:
            if offline:
                raise
            download = spec["download"]
            size = transfer_size(download)
            if size > max_bytes:
                raise ValueError(
                    f"{name} requires {size} transfer bytes; budget is {max_bytes}"
                ) from None
            registry.download(name, version, timeout=60, max_retries=2)
            dataset = registry.open(name, version)
    if dataset.manifest.to_dict() != spec["manifest"]:
        raise ValueError("Installed quantitative manifest disagrees with pinned release")
    return dataset


def _validate_contract(spec, contract):
    manifest = spec["manifest"]
    references = manifest["references"]
    if len(references) != 1 or references[0]["taxon_id"] != contract.taxon_id:
        raise ValueError("Managed expression requires one matching reference context")
    if contract.content_provenance.get("quantitative_manifest") != _digest(manifest):
        raise ValueError("Quantitative manifest and expression contract disagree")
    pin = spec["download"]["assets"]["dataset.sqlite"]
    if contract.content_provenance.get("quantitative_artifact") != pin["sha256"]:
        raise ValueError("Quantitative artifact and expression contract disagree")
    sources = {row["asset_id"]: row for row in manifest["sources"]}
    for role, asset in spec["source_roles"].items():
        if contract.content_provenance.get(role) != sources[asset]["sha256"]:
            raise ValueError("Quantitative source and expression contract disagree")
    selection = spec["selection"]
    if (
        selection["kind"] != contract.measurement.feature_level
        or selection["unit"] != contract.measurement.unit
        or selection["scale"] != contract.measurement.scale
    ):
        raise ValueError("Quantitative selection and expression measurement disagree")


def read_managed_expression(entry, samples, contract, **options):
    """One-to-one native observations, never aggregation or isoform allocation."""
    from .species_expression import ExpressionDataset

    spec = specification(entry)
    _validate_contract(spec, contract)
    dataset = open_managed(spec, taxon_id=contract.taxon_id, **options)
    reference = dataset.manifest.references[0]
    native_samples = {row.sample_id: row.to_dict() for row in dataset.samples()}
    if set(native_samples) != set(samples.index):
        raise ValueError("Quantitative samples and curated metadata disagree")
    for sample_id, row in native_samples.items():
        if (
            row["reference"] != reference.to_dict()
            or row["tissue"] != samples.at[sample_id, "source_tissue"]
            or row["attributes"].get("source_attributes")
            != samples.at[sample_id, "source_attributes"]
        ):
            raise ValueError("Quantitative sample identity and curated source facts disagree")
    selection = spec["selection"]
    features, rows, seen = {}, [], set()
    for observation in dataset.measurements(
        **{key: selection[key] for key in ("quantity", "unit", "scale")}
    ):
        target = observation.target
        if (
            target.reference != reference
            or target.kind != selection["kind"]
            or target.namespace != selection["namespace"]
            or observation.sample_id not in native_samples
            or observation.status not in {"measured", "missing"}
            or observation.origin != "observed"
        ):
            raise ValueError("Quantitative observation identity/semantics disagree")
        feature_id = target.identifier
        key = feature_id, observation.sample_id
        if key in seen:
            raise ValueError("Duplicate feature/sample observations cannot be silently aggregated")
        seen.add(key)
        fields = {
            name: observation.provenance[field] for name, field in spec["feature_fields"].items()
        }
        feature = {"source_feature_id": feature_id, "source_target": target.to_dict(), **fields}
        if feature_id in features and features[feature_id] != feature:
            raise ValueError("Conflicting feature metadata in native observations")
        features[feature_id] = feature
        rows.append((feature_id, observation.sample_id, observation.value))
    if not rows:
        raise ValueError("No native measurements match the registered selection")
    values = pd.DataFrame(rows, columns=["feature_id", "sample_id", "value"]).pivot(
        index="feature_id", columns="sample_id", values="value"
    )
    values = values.reindex(columns=sorted(samples.index))
    frame = pd.DataFrame.from_dict(features, orient="index").loc[values.index]
    frame.index.name = "feature_id"
    provenance = {
        "producer": "orfdata",
        "manifest": spec["manifest"],
        "download": spec["download"],
        "selection": selection,
        "aggregation": "none; native feature/sample observations",
        "record_scope": "Complete source-native records from this quantitative release",
        "experiment": entry["accession"],
    }

    def records(kind):
        return (row.to_dict() for row in getattr(dataset, kind)())

    return ExpressionDataset(
        values,
        frame,
        samples.loc[values.columns],
        contract,
        records=ExpressionRecords(provenance, records),
    )

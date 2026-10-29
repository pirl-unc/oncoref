# Licensed under the Apache License, Version 2.0.
"""Small species manifests; datacache owns transfer, validation and publication."""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse

from .reference_identity import ReferenceIdentity, _label, _string_map
from .species import DATA_GROUPS, resolve_species

_MANIFEST_ROOT = Path(__file__).parent / "data" / "species"
DEFAULT_MAX_BYTES = 100 * 1024 * 1024


@dataclass(frozen=True)
class SpeciesAsset:
    id: str
    taxon_id: int
    group: str
    filename: str
    url: str
    sha256: str
    size_bytes: int
    source_version: str
    description: str
    license: str

    def __post_init__(self):
        for name in ("id", "filename"):
            value = getattr(self, name)
            if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
                raise ValueError(f"Asset {name} must be a plain identifier without path components")
        if type(self.taxon_id) is not int or self.taxon_id < 1:
            raise ValueError("Positive integer taxonomy identifier required")
        if self.group not in DATA_GROUPS:
            raise ValueError("Unknown asset data group")
        if urlparse(self.url).scheme not in {"https", "http"} or not urlparse(self.url).netloc:
            raise ValueError("Asset URL must be an absolute HTTP(S) URL")
        digest = _string_map({"asset": self.sha256}, "asset", digests=True)["asset"]
        object.__setattr__(self, "sha256", digest)
        if type(self.size_bytes) is not int or self.size_bytes < 1:
            raise ValueError("Asset requires a positive pinned byte size")
        for name in ("source_version", "description", "license"):
            _label(getattr(self, name), name)

    def as_dict(self):
        return {key: getattr(self, key) for key in self.__dataclass_fields__}


def _manifest(species):
    selected = resolve_species(species)
    path = _MANIFEST_ROOT / selected.name / "assets.json"
    if not path.exists():
        return {"assets": [], "references": []}
    record = json.loads(path.read_text())
    if set(record) != {"schema", "schema_version", "taxon_id", "assets", "references"} or (
        record["schema"] != "oncoref.species-assets"
        or type(record["schema_version"]) is not int
        or record["schema_version"] != 1
        or type(record["taxon_id"]) is not int
        or record["taxon_id"] != selected.taxon_id
    ):
        raise ValueError("Invalid species asset manifest schema or taxon")
    assets = [SpeciesAsset(**row) for row in record["assets"]]
    identities = {(a.id, a.source_version) for a in assets}
    if len(identities) != len(assets) or any(a.taxon_id != selected.taxon_id for a in assets):
        raise ValueError("Duplicate asset IDs or mixed species in manifest")
    seen = set()
    for row in record["references"]:
        if set(row) != {"id", "identity", "assets", "notes"} or row["id"] in seen:
            raise ValueError("Invalid or duplicate registered reference")
        _label(row["id"], "reference registration ID")
        _label(row["notes"], "reference notes")
        seen.add(row["id"])
        reference = ReferenceIdentity.from_dict(row["identity"])
        if reference.taxon_id != selected.taxon_id or set(row["assets"]) != set(
            reference.content_provenance
        ):
            raise ValueError("Registered reference taxon/asset roles disagree")
        for role, asset_id in row["assets"].items():
            asset = next(
                (
                    a
                    for a in assets
                    if a.id == asset_id and a.sha256 == reference.content_provenance[role]
                ),
                None,
            )
            if (
                asset is None
                or asset.group != "references"
                or asset.sha256 != reference.content_provenance[role]
            ):
                raise ValueError("Registered reference and asset digest disagree")
    return record


def species_assets(species, *, group=None) -> tuple[SpeciesAsset, ...]:
    """Newest registered artifact per stable ID (manifest order is release order)."""
    if group is not None and group not in DATA_GROUPS:
        raise ValueError(f"Unknown data group {group!r}")
    latest = {row["id"]: SpeciesAsset(**row) for row in _manifest(species)["assets"]}
    from .managed_expression import managed_assets

    latest.update((asset.id, asset) for asset in managed_assets(species))
    return tuple(asset for asset in latest.values() if group is None or asset.group == group)


def asset_versions(asset_id: str, *, species) -> tuple[str, ...]:
    """Concrete releases, newest first; adding a release preserves older rows."""
    from .managed_expression import managed_assets

    managed = next((a for a in managed_assets(species) if a.id == asset_id), None)
    if managed is not None:
        return (managed.source_version,)
    versions = tuple(
        row["source_version"]
        for row in reversed(_manifest(species)["assets"])
        if row["id"] == asset_id
    )
    if not versions:
        raise KeyError(f"Unknown {resolve_species(species).name} asset {asset_id!r}")
    return versions


def species_asset(asset_id: str, *, species, version: str | None = None) -> SpeciesAsset:
    """Resolve an exact source release; latest never becomes a cache identity."""
    available = asset_versions(asset_id, species=species)
    concrete = available[0] if version in (None, "latest") else version
    from .managed_expression import managed_assets

    managed = next((a for a in managed_assets(species) if a.id == asset_id), None)
    if managed is not None and concrete == managed.source_version:
        return managed
    for row in _manifest(species)["assets"]:
        if row["id"] == asset_id and row["source_version"] == concrete:
            return SpeciesAsset(**row)
    raise ValueError(f"{asset_id!r} has no version {version!r}; available: {', '.join(available)}")


def registered_references(species) -> tuple[dict, ...]:
    """Independent records; two registrations are never implicit assembly aliases."""
    return tuple(_manifest(species)["references"])


def asset_cache_dir(*, cache_root=None) -> Path:
    from datacache import Cache

    root = cache_root if cache_root is not None else os.environ.get("ONCOREF_ASSET_CACHE")
    return Path(Cache("oncoref", cache_root=root).cache_directory_path) / "species"


def _overrides(local_overrides):
    if local_overrides is None:
        file = os.environ.get("ONCOREF_ASSET_OVERRIDES")
        local_overrides = json.loads(Path(file).read_text()) if file else {}
    if not isinstance(local_overrides, dict) or any(
        not isinstance(key, str) or not isinstance(value, (str, Path))
        for key, value in local_overrides.items()
    ):
        raise ValueError("Asset overrides must map species-qualified asset IDs to local paths")
    return local_overrides


def asset_path(asset, *, cache_root=None, local_overrides=None) -> Path:
    """Compute a path without reads/writes. Overrides are verified before use."""
    from .managed_expression import ManagedExpressionAsset, inspect_managed

    if isinstance(asset, ManagedExpressionAsset):
        state = inspect_managed(
            asset, cache_root=cache_root, local_overrides=local_overrides, verify_files=False
        )
        if state["path"] is None:
            raise ImportError(state["error"])
        return Path(state["path"])
    overrides = _overrides(local_overrides)
    # Taxonomy prefixes prevent similarly named assets in different species colliding.
    key = f"{asset.taxon_id}:{asset.id}"
    if key in overrides:
        return Path(overrides[key]).expanduser()
    return (
        asset_cache_dir(cache_root=cache_root)
        / str(asset.taxon_id)
        / asset.group
        / asset.sha256
        / asset.filename
    )


def inspect_asset(asset, *, cache_root=None, local_overrides=None):
    """Read-only integrity check; never downloads, repairs or copies an override."""
    from datacache import inspect_file

    from .managed_expression import ManagedExpressionAsset, inspect_managed

    if isinstance(asset, ManagedExpressionAsset):
        return inspect_managed(asset, cache_root=cache_root, local_overrides=local_overrides)

    path = asset_path(asset, cache_root=cache_root, local_overrides=local_overrides)
    state = inspect_file(path, expected_sha256=asset.sha256, expected_size=asset.size_bytes)
    return {
        "asset_id": asset.id,
        "taxon_id": asset.taxon_id,
        "group": asset.group,
        "path": str(path),
        "status": state.status,
        "size_bytes": asset.size_bytes,
        "sha256": asset.sha256,
        "error": str(state.error) if state.error else None,
    }


def ensure_asset(
    asset, *, offline=False, cache_root=None, local_overrides=None, max_bytes=DEFAULT_MAX_BYTES
):
    """Verify/reuse local bytes or acquire one pinned asset within a transfer budget.

    An explicit override is authoritative: corrupt/missing override files fail,
    including online, and are never overwritten. Cached corruption also fails
    offline; online datacache can repair the managed cache atomically. On POSIX,
    HTTP resume bounds retained partial bytes by the manifest's exact expected
    size. Other platforms use a verified, non-resumable download.
    """
    from datacache import fetch_file

    from .managed_expression import ManagedExpressionAsset, open_managed

    if isinstance(asset, ManagedExpressionAsset):
        return open_managed(
            asset.spec,
            taxon_id=asset.taxon_id,
            offline=offline,
            cache_root=cache_root,
            local_overrides=local_overrides,
            max_bytes=max_bytes,
        ).path

    if type(max_bytes) is not int or max_bytes < 0:
        raise ValueError("max_bytes must be a nonnegative integer transfer budget")
    overrides = _overrides(local_overrides)
    state = inspect_asset(asset, cache_root=cache_root, local_overrides=overrides)
    if state["status"] == "available":
        return Path(state["path"])
    if offline or f"{asset.taxon_id}:{asset.id}" in overrides:
        raise ValueError(
            f"{asset.id}: {state['status']}: {state['error'] or 'verified local asset required'}"
        )
    if asset.size_bytes > max_bytes:
        raise ValueError(
            f"{asset.id} requires {asset.size_bytes} bytes; transfer budget is {max_bytes}"
        )
    return Path(
        fetch_file(
            asset.url,
            destination=state["path"],
            force=state["status"] == "corrupt",
            raw=True,
            resume=os.name == "posix",
            expected_sha256=asset.sha256,
            expected_size=asset.size_bytes,
            timeout=60,
            max_retries=2,
            chunk_size=64 * 1024,
            record_provenance=True,
        )
    )


def ensure_reference(
    reference_id,
    *,
    species,
    offline=False,
    cache_root=None,
    local_overrides=None,
    max_bytes=DEFAULT_MAX_BYTES,
):
    """Verify a complete declared asset inventory before returning its identity."""
    if type(max_bytes) is not int or max_bytes < 0:
        raise ValueError("max_bytes must be a nonnegative integer transfer budget")
    row = next((r for r in registered_references(species) if r["id"] == reference_id), None)
    if row is None:
        raise KeyError(f"Unknown registered reference {reference_id!r}")
    reference = ReferenceIdentity.from_dict(row["identity"])
    registered_assets = [SpeciesAsset(**a) for a in _manifest(species)["assets"]]
    assets = {
        role: next(
            a
            for a in registered_assets
            if a.id == asset_id and a.sha256 == reference.content_provenance[role]
        )
        for role, asset_id in row["assets"].items()
    }
    missing_bytes = sum(
        a.size_bytes
        for a in assets.values()
        if inspect_asset(a, cache_root=cache_root, local_overrides=local_overrides)["status"]
        != "available"
    )
    if not offline and missing_bytes > max_bytes:
        raise ValueError(f"Reference needs {missing_bytes} bytes; transfer budget is {max_bytes}")
    paths = {
        role: ensure_asset(
            asset,
            offline=offline,
            cache_root=cache_root,
            local_overrides=local_overrides,
            max_bytes=max_bytes,
        )
        for role, asset in assets.items()
    }
    return ReferenceIdentity.from_dict(row["identity"]), paths


def asset_inventory(species, *, group=None, cache_root=None, local_overrides=None):
    from .managed_expression import ManagedExpressionAsset, inspect_managed

    def available(asset):
        if isinstance(asset, ManagedExpressionAsset):
            return (
                inspect_managed(
                    asset,
                    cache_root=cache_root,
                    local_overrides=local_overrides,
                    verify_files=False,
                )["status"]
                == "available"
            )
        return asset_path(asset, cache_root=cache_root, local_overrides=local_overrides).is_file()

    selected = resolve_species(species)
    return [
        {
            "name": asset.id,
            "held": "asset",
            "category": asset.group,
            "group": asset.group,
            "available": available(asset),
            "integrity": "not_checked",
            "cohorts": None,
            "description": asset.description,
            "species": selected.name,
            "taxon_id": selected.taxon_id,
            "size_bytes": asset.size_bytes,
            "default_version": asset.source_version,
            "sha256": asset.sha256,
        }
        for asset in species_assets(selected, group=group)
    ]

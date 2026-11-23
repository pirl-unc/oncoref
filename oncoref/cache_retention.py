# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Bounded cache inventory and safe retention pruning across OncoRef cache families (#591).

Read-only inventory never creates directories, lockfiles, or network requests, and
deduplicates hard-linked inodes ``(st_dev, st_ino)`` so overlay bundles report both
logical bytes and true physical reclaimable bytes.
"""

from __future__ import annotations

import contextlib
import json
import os
import re
from collections.abc import Iterable
from pathlib import Path

from filelock import FileLock, Timeout

CACHE_TARGETS: tuple[str, ...] = ("bundle", "source", "reference", "sequence", "species")
_SEMVER_RE = re.compile(r"^v?(\d+)\.(\d+)\.(\d+)$")
_HEX64_RE = re.compile(r"^[0-9a-f]{64}$")


def _normalize_targets(targets: str | Iterable[str] | None) -> tuple[str, ...]:
    if targets is None:
        return CACHE_TARGETS
    raw = [targets] if isinstance(targets, str) else list(targets)
    if not raw or "all" in raw:
        return CACHE_TARGETS
    resolved: list[str] = []
    for item in raw:
        token = str(item).strip().lower()
        if token == "hpa":
            token = "reference"
        if token == "all":
            return CACHE_TARGETS
        if token not in CACHE_TARGETS:
            raise ValueError(
                f"Unknown cache prune target {item!r}; choose from ('all', 'hpa', {', '.join(repr(t) for t in CACHE_TARGETS)})"
            )
        if token not in resolved:
            resolved.append(token)
    return tuple(resolved)


def _semver_key(label: str) -> tuple[int, int, int] | None:
    match = _SEMVER_RE.match(str(label).strip())
    if not match:
        return None
    return int(match.group(1)), int(match.group(2)), int(match.group(3))


def _scan_inodes(path: Path) -> dict[tuple[int, int], int]:
    """Map ``(st_dev, st_ino) -> st_size`` for regular files under ``path`` without following symlinks."""
    inodes: dict[tuple[int, int], int] = {}
    try:
        st = os.lstat(path)
    except OSError:
        return inodes
    if os.path.islink(path):
        return inodes
    if path.is_file():
        inodes[(st.st_dev, st.st_ino)] = st.st_size
        return inodes
    if not path.is_dir():
        return inodes
    for root, dirs, files in os.walk(path, followlinks=False):
        root_path = Path(root)
        dirs[:] = [d for d in dirs if not (root_path / d).is_symlink()]
        for fname in files:
            fpath = root_path / fname
            try:
                fst = os.lstat(fpath)
            except OSError:
                continue
            if os.path.islink(fpath):
                continue
            inodes[(fst.st_dev, fst.st_ino)] = fst.st_size
    return inodes


def _active_lock_exists(lock_paths: Iterable[Path]) -> bool:
    """Return True only when an existing lockfile is actively held by another process/handle."""
    for lock_path in lock_paths:
        if not lock_path.exists():
            continue
        try:
            with FileLock(str(lock_path), timeout=0):
                pass
        except Timeout:
            return True
        except OSError:
            # Unwritable lock directory is conservatively treated as protected.
            return True
    return False


def _bundle_overlay_base_version(bundle_dir: Path) -> str | None:
    from . import data_bundle

    marker = data_bundle._read_completion_marker(bundle_dir)
    if isinstance(marker, dict):
        rel = marker.get("release_manifest")
        if isinstance(rel, dict) and rel.get("bundle_layout") == "overlay":
            base = rel.get("base_bundle")
            if isinstance(base, dict) and base.get("data_version"):
                return str(base["data_version"]).removeprefix("v")
    overlay_file = bundle_dir / "overlay_manifest.json"
    if overlay_file.is_file() and not overlay_file.is_symlink():
        with contextlib.suppress(OSError, json.JSONDecodeError):
            raw = json.loads(overlay_file.read_text())
            base = raw.get("base_bundle") or {}
            if isinstance(base, dict) and base.get("data_version"):
                return str(base["data_version"]).removeprefix("v")
            if raw.get("base_version"):
                return str(raw["base_version"]).removeprefix("v")
    return None


def _collect_bundle_entries(
    *,
    keep_current: bool,
    keep_set: set[str],
    include_legacy_pirlygenes: bool,
) -> list[dict]:
    from . import data_bundle

    primary_root = data_bundle.cache_root()
    current_store = data_bundle.store_dir()
    current_dir = data_bundle.cache_dir()
    roots: list[tuple[str, Path]] = [("oncoref", primary_root)]
    legacy_root = data_bundle._LEGACY_CACHE_PARENT
    if include_legacy_pirlygenes and legacy_root != primary_root:
        roots.append(("legacy_pirlygenes", legacy_root))

    discovered: list[tuple[str, Path, str]] = []
    by_version: dict[str, list[Path]] = {}
    for namespace, root in roots:
        if not root.exists() or root.is_symlink():
            continue
        for child in sorted(root.iterdir()):
            if child.name.startswith("."):
                continue
            ver = child.name.removeprefix("v") if _semver_key(child.name) else None
            discovered.append((namespace, child, ver or ""))
            if ver and child.is_dir() and not child.is_symlink():
                by_version.setdefault(ver, []).append(child)

    retained_versions: set[str] = set()
    if keep_current:
        retained_versions.add(data_bundle.DATA_VERSION)
    for token in keep_set:
        retained_versions.add(token.removeprefix("v"))

    # Transitively retain any overlay base versions required by retained bundles.
    overlay_bases: set[str] = set()
    queue = list(retained_versions)
    seen = set(queue)
    while queue:
        ver = queue.pop()
        for bdir in by_version.get(ver, ()):
            base_ver = _bundle_overlay_base_version(bdir)
            if base_ver and base_ver not in seen:
                seen.add(base_ver)
                overlay_bases.add(base_ver)
                queue.append(base_ver)

    entries: list[dict] = []
    for namespace, child, ver in discovered:
        lock_file = child.parent / f".{child.name}.oncoref.lock"
        if child.is_symlink():
            status, reason, prunable = "protected", "symlink_protected", False
        elif not ver or not child.is_dir():
            status, reason, prunable = "protected", "unmanaged_custom_entry", False
        elif _active_lock_exists([lock_file]):
            status, reason, prunable = "protected", "active_lock", False
        elif keep_current and child.resolve() in {current_store.resolve(), current_dir.resolve()}:
            status, reason, prunable = "retained", "active_pin", False
        elif ver in keep_set or f"v{ver}" in keep_set:
            status, reason, prunable = "retained", "explicit_keep_version", False
        elif ver in overlay_bases:
            status, reason, prunable = "retained", "overlay_base_dependency", False
        elif keep_current and ver == data_bundle.DATA_VERSION and namespace == "oncoref":
            status, reason, prunable = "retained", "active_pin", False
        else:
            status = "prunable"
            reason = (
                "legacy_pirlygenes_duplicate"
                if namespace == "legacy_pirlygenes"
                else "superseded_bundle_version"
            )
            prunable = True
        entries.append(
            {
                "target": "bundle",
                "namespace": namespace,
                "id": child.name,
                "version": f"v{ver}" if ver else child.name,
                "path": str(child),
                "status": status,
                "reason": reason,
                "prunable": prunable,
                "_path_obj": child,
            }
        )
    return entries


def _collect_source_entries(*, keep_current: bool, keep_set: set[str]) -> list[dict]:
    from . import source_matrices

    root = source_matrices.cache_root()
    if not root.exists() or root.is_symlink():
        return []

    pinned_by_code = {
        str(row.cancer_code): str(row.source_matrix_version)
        for row in source_matrices.registry().itertuples(index=False)
    }

    per_code_versions: dict[str, list[tuple[tuple[int, int, int], str]]] = {}
    raw_items: list[tuple[Path, str | None, str | None]] = []
    for vdir in sorted(root.iterdir()):
        if vdir.name.startswith("."):
            continue
        if vdir.is_symlink() or not vdir.is_dir() or _semver_key(vdir.name) is None:
            raw_items.append((vdir, None, None))
            continue
        ver = vdir.name.removeprefix("v")
        vkey = _semver_key(ver)
        assert vkey is not None
        for item in sorted(vdir.iterdir()):
            if item.name.startswith("."):
                continue
            if item.is_file() and not item.is_symlink() and item.suffix == ".parquet":
                code = item.stem
                raw_items.append((item, code, ver))
                per_code_versions.setdefault(code, []).append((vkey, ver))
            else:
                raw_items.append((item, None, ver))

    newest_local_by_code = {
        code: max(vers, key=lambda pair: pair[0])[1] for code, vers in per_code_versions.items()
    }

    entries: list[dict] = []
    for item, code, ver in raw_items:
        lock_candidates = [item.with_name(f".{item.name}.oncoref.lock")]
        if item.is_symlink():
            status, reason, prunable = "protected", "symlink_protected", False
        elif code is None or ver is None or code not in pinned_by_code:
            status, reason, prunable = "protected", "unmanaged_custom_entry", False
        elif _active_lock_exists(lock_candidates):
            status, reason, prunable = "protected", "active_lock", False
        elif ver in keep_set or f"v{ver}" in keep_set or f"{code}@{ver}" in keep_set:
            status, reason, prunable = "retained", "explicit_keep_version", False
        elif keep_current and pinned_by_code.get(code) == ver:
            status, reason, prunable = "retained", "active_cohort_pin", False
        elif keep_current and newest_local_by_code.get(code) == ver:
            status, reason, prunable = "retained", "newest_cached_cohort_version", False
        else:
            status, reason, prunable = "prunable", "superseded_cohort_matrix", True
        entries.append(
            {
                "target": "source",
                "namespace": "oncoref",
                "id": f"{code}@v{ver}" if code and ver else item.name,
                "version": f"v{ver}" if ver else "",
                "path": str(item),
                "status": status,
                "reason": reason,
                "prunable": prunable,
                "_path_obj": item,
            }
        )
    return entries


def _collect_reference_entries(*, keep_current: bool, keep_set: set[str]) -> list[dict]:
    from . import curation, reference_data

    base = reference_data.cache_base_dir()
    root = base / "sources"
    if not root.exists() or root.is_symlink():
        return []

    active_pairs: set[tuple[str, str]] = set()
    if keep_current:
        for name, spec in reference_data.REFERENCE_SOURCES.items():
            active_pairs.add((name, reference_data.resolve_version(name)))
            if reference_data.DEFAULT_HPA_VERSION in spec.get("urls", {}):
                active_pairs.add((name, reference_data.DEFAULT_HPA_VERSION))
        active_pairs.update(curation.required_reference_inputs())

    entries: list[dict] = []
    for source_dir in sorted(root.iterdir()):
        if source_dir.name.startswith(".") or source_dir.name == "manifest.json":
            continue
        if source_dir.is_symlink():
            entries.append(
                {
                    "target": "reference",
                    "namespace": "oncoref",
                    "id": source_dir.name,
                    "version": "",
                    "path": str(source_dir),
                    "status": "protected",
                    "reason": "symlink_protected",
                    "prunable": False,
                    "_path_obj": source_dir,
                }
            )
            continue
        if not source_dir.is_dir() or source_dir.name not in reference_data.REFERENCE_SOURCES:
            entries.append(
                {
                    "target": "reference",
                    "namespace": "oncoref",
                    "id": source_dir.name,
                    "version": "",
                    "path": str(source_dir),
                    "status": "protected",
                    "reason": "unmanaged_custom_entry",
                    "prunable": False,
                    "_path_obj": source_dir,
                }
            )
            continue
        name = source_dir.name
        known_versions = set(reference_data.REFERENCE_SOURCES[name].get("urls", ()))
        for vdir in sorted(source_dir.iterdir()):
            if vdir.name.startswith("."):
                continue
            ver = vdir.name
            lock_candidates = [
                *vdir.glob(".*.lock"),
                source_dir / f".{ver}.lock",
                vdir / ".oncoref.lock",
            ]
            if vdir.is_symlink():
                status, reason, prunable = "protected", "symlink_protected", False
            elif not vdir.is_dir():
                status, reason, prunable = "protected", "unmanaged_custom_entry", False
            elif _active_lock_exists(lock_candidates):
                status, reason, prunable = "protected", "active_lock", False
            elif ver in keep_set or f"{name}@{ver}" in keep_set:
                status, reason, prunable = "retained", "explicit_keep_version", False
            elif (name, ver) in active_pairs:
                status, reason, prunable = "retained", "active_reference_pin", False
            elif ver not in known_versions and not ver.startswith("v"):
                status, reason, prunable = "protected", "unmanaged_custom_entry", False
            else:
                status, reason, prunable = "prunable", "superseded_reference_version", True
            entries.append(
                {
                    "target": "reference",
                    "namespace": "oncoref",
                    "id": f"{name}@{ver}",
                    "version": ver,
                    "path": str(vdir),
                    "status": status,
                    "reason": reason,
                    "prunable": prunable,
                    "_path_obj": vdir,
                }
            )
    return entries


def _collect_sequence_entries(*, keep_current: bool, keep_set: set[str]) -> list[dict]:
    from . import sequence_bundle

    root = sequence_bundle.cache_root()
    if not root.exists() or root.is_symlink():
        return []

    pin_digest = sequence_bundle.descriptor()["bundle_id"].rsplit(":", 1)[-1]
    entries: list[dict] = []
    for child in sorted(root.iterdir()):
        if child.name.startswith(".") or child.name.endswith(".lock"):
            continue
        lock_file = root / f"{child.name}.lock"
        if child.is_symlink():
            status, reason, prunable = "protected", "symlink_protected", False
        elif not child.is_dir() or not _HEX64_RE.match(child.name):
            status, reason, prunable = "protected", "unmanaged_custom_entry", False
        elif _active_lock_exists([lock_file]):
            status, reason, prunable = "protected", "active_lock", False
        elif (
            child.name in keep_set or f"oncoref:sequence-bundle:v1:sha256:{child.name}" in keep_set
        ):
            status, reason, prunable = "retained", "explicit_keep_version", False
        elif keep_current and child.name == pin_digest:
            status, reason, prunable = "retained", "active_sequence_bundle_pin", False
        else:
            status, reason, prunable = "prunable", "superseded_sequence_bundle", True
        entries.append(
            {
                "target": "sequence",
                "namespace": "oncoref",
                "id": child.name,
                "version": child.name[:12],
                "path": str(child),
                "status": status,
                "reason": reason,
                "prunable": prunable,
                "_path_obj": child,
            }
        )
    return entries


def _collect_species_entries(*, keep_current: bool, keep_set: set[str]) -> list[dict]:
    from . import species, species_assets

    env_root = os.environ.get("ONCOREF_ASSET_CACHE")
    if env_root is not None:
        base = Path(env_root).expanduser()
        candidates = [base / "species", base / "oncoref" / "species"]
        root = next((c for c in candidates if c.exists()), candidates[0])
    else:
        root = species_assets.asset_cache_dir()
    if not root.exists() or root.is_symlink():
        return []

    pinned_digests: set[tuple[str, str, str]] = set()
    if keep_current:
        for sp in species.available_species():
            for asset in species_assets.species_assets(sp):
                pinned_digests.add((str(asset.taxon_id), asset.group, asset.sha256))

    entries: list[dict] = []
    for taxon_dir in sorted(root.iterdir()):
        if taxon_dir.name.startswith("."):
            continue
        if taxon_dir.is_symlink() or not taxon_dir.is_dir() or not taxon_dir.name.isdigit():
            entries.append(
                {
                    "target": "species",
                    "namespace": "oncoref",
                    "id": taxon_dir.name,
                    "version": "",
                    "path": str(taxon_dir),
                    "status": "protected",
                    "reason": "symlink_protected"
                    if taxon_dir.is_symlink()
                    else "unmanaged_custom_entry",
                    "prunable": False,
                    "_path_obj": taxon_dir,
                }
            )
            continue
        for group_dir in sorted(taxon_dir.iterdir()):
            if group_dir.name.startswith(".") or not group_dir.is_dir() or group_dir.is_symlink():
                continue
            for digest_dir in sorted(group_dir.iterdir()):
                if digest_dir.name.startswith("."):
                    continue
                lock_candidates = [
                    *digest_dir.glob("*.lock"),
                    *group_dir.glob(f".{digest_dir.name}*.lock"),
                ]
                key = (taxon_dir.name, group_dir.name, digest_dir.name)
                if digest_dir.is_symlink():
                    status, reason, prunable = "protected", "symlink_protected", False
                elif not digest_dir.is_dir() or not _HEX64_RE.match(digest_dir.name):
                    status, reason, prunable = "protected", "unmanaged_custom_entry", False
                elif _active_lock_exists(lock_candidates):
                    status, reason, prunable = "protected", "active_lock", False
                elif digest_dir.name in keep_set:
                    status, reason, prunable = "retained", "explicit_keep_version", False
                elif key in pinned_digests:
                    status, reason, prunable = "retained", "active_species_asset_pin", False
                else:
                    status, reason, prunable = "prunable", "superseded_species_asset", True
                entries.append(
                    {
                        "target": "species",
                        "namespace": "oncoref",
                        "id": f"{taxon_dir.name}:{group_dir.name}:{digest_dir.name[:12]}",
                        "version": digest_dir.name[:12],
                        "path": str(digest_dir),
                        "status": status,
                        "reason": reason,
                        "prunable": prunable,
                        "_path_obj": digest_dir,
                    }
                )
    return entries


def cache_inventory(
    *,
    targets: str | Iterable[str] | None = None,
    keep_current: bool = True,
    keep_versions: Iterable[str] = (),
    include_legacy_pirlygenes: bool = False,
) -> dict:
    """Read-only bounded cache inventory with hard-link-aware physical byte accounting."""
    selected = _normalize_targets(targets)
    keep_set = {str(v).strip() for v in keep_versions if str(v).strip()}

    raw_entries: list[dict] = []
    if "bundle" in selected:
        raw_entries.extend(
            _collect_bundle_entries(
                keep_current=keep_current,
                keep_set=keep_set,
                include_legacy_pirlygenes=include_legacy_pirlygenes,
            )
        )
    if "source" in selected:
        raw_entries.extend(_collect_source_entries(keep_current=keep_current, keep_set=keep_set))
    if "reference" in selected:
        raw_entries.extend(_collect_reference_entries(keep_current=keep_current, keep_set=keep_set))
    if "sequence" in selected:
        raw_entries.extend(_collect_sequence_entries(keep_current=keep_current, keep_set=keep_set))
    if "species" in selected:
        raw_entries.extend(_collect_species_entries(keep_current=keep_current, keep_set=keep_set))

    inode_maps = [_scan_inodes(entry["_path_obj"]) for entry in raw_entries]
    all_inodes: dict[tuple[int, int], int] = {}
    kept_inodes: dict[tuple[int, int], int] = {}
    prunable_inodes: dict[tuple[int, int], int] = {}

    for entry, inodes in zip(raw_entries, inode_maps):
        all_inodes.update(inodes)
        if entry["prunable"]:
            prunable_inodes.update(inodes)
        else:
            kept_inodes.update(inodes)

    public_entries: list[dict] = []
    for entry, inodes in zip(raw_entries, inode_maps):
        logical_bytes = sum(inodes.values())
        reclaimable_bytes = (
            sum(size for key, size in inodes.items() if key not in kept_inodes)
            if entry["prunable"]
            else 0
        )
        public_entries.append(
            {
                k: v
                for k, v in {
                    **entry,
                    "logical_bytes": logical_bytes,
                    "reclaimable_bytes": reclaimable_bytes,
                }.items()
                if k != "_path_obj"
            }
        )

    reclaimable_total = sum(size for key, size in prunable_inodes.items() if key not in kept_inodes)
    return {
        "targets": list(selected),
        "keep_current": bool(keep_current),
        "keep_versions": sorted(keep_set),
        "include_legacy_pirlygenes": bool(include_legacy_pirlygenes),
        "logical_bytes": sum(e["logical_bytes"] for e in public_entries),
        "unique_bytes": sum(all_inodes.values()),
        "retained_unique_bytes": sum(kept_inodes.values()),
        "reclaimable_bytes": reclaimable_total,
        "prunable_count": sum(1 for e in public_entries if e["prunable"]),
        "entries": public_entries,
    }


def _delete_path_safely(path: Path) -> bool:
    """Remove a non-symlink file or directory tree without following symlinks."""
    if not path.exists() or path.is_symlink():
        return False
    if path.is_file():
        try:
            path.unlink()
            return True
        except OSError:
            return False
    for root, dirs, files in os.walk(path, topdown=False, followlinks=False):
        root_path = Path(root)
        for fname in files:
            fpath = root_path / fname
            with contextlib.suppress(OSError):
                fpath.unlink()
        for dname in dirs:
            dpath = root_path / dname
            with contextlib.suppress(OSError):
                if dpath.is_symlink():
                    dpath.unlink()
                else:
                    dpath.rmdir()
    try:
        path.rmdir()
        return True
    except OSError:
        return False


def _cleanup_empty_version_dirs(deleted_paths: Iterable[Path]) -> None:
    parents = {p.parent for p in deleted_paths}
    for parent in sorted(parents, key=lambda p: len(p.parts), reverse=True):
        if not parent.exists() or parent.is_symlink() or not parent.is_dir():
            continue
        # Clean up empty source-matrix vX.Y.Z or reference <source> dirs when drained.
        remaining = [c for c in parent.iterdir() if not c.name.endswith(".oncoref.lock")]
        if not remaining:
            for lock in parent.glob("*.oncoref.lock"):
                with contextlib.suppress(OSError):
                    lock.unlink()
            with contextlib.suppress(OSError):
                parent.rmdir()


def prune_caches(
    *,
    targets: str | Iterable[str] | None = None,
    keep_current: bool = True,
    keep_versions: Iterable[str] = (),
    include_legacy_pirlygenes: bool = False,
    dry_run: bool = True,
) -> dict:
    """Inspect or delete obsolete unpinned cache entries across OncoRef cache roots.

    Defaults to ``dry_run=True`` so callers must opt into deletion explicitly.
    """
    before = cache_inventory(
        targets=targets,
        keep_current=keep_current,
        keep_versions=keep_versions,
        include_legacy_pirlygenes=include_legacy_pirlygenes,
    )
    planned = [entry for entry in before["entries"] if entry["prunable"]]
    if not dry_run and "bundle" in _normalize_targets(targets):
        import datacache.bundle_store as _bs

        for entry in before["entries"]:
            if entry["target"] == "bundle" and not entry["prunable"]:
                bpath = Path(entry["path"])
                with contextlib.suppress(*_bs.INVALID_STORE_ERRORS):
                    if _bs.store_at(bpath) is not None:
                        _bs.prune_bundles(bpath, keep=1)
    if dry_run or not planned:
        return {
            **before,
            "dry_run": bool(dry_run),
            "deleted_count": 0,
            "deleted": [],
            "reclaimed_bytes": 0,
        }

    deleted: list[dict] = []
    deleted_paths: list[Path] = []
    for entry in planned:
        path = Path(entry["path"])
        if _delete_path_safely(path):
            deleted.append(entry)
            deleted_paths.append(path)
    _cleanup_empty_version_dirs(deleted_paths)

    after = cache_inventory(
        targets=targets,
        keep_current=keep_current,
        keep_versions=keep_versions,
        include_legacy_pirlygenes=include_legacy_pirlygenes,
    )
    reclaimed = max(0, before["unique_bytes"] - after["unique_bytes"])
    return {
        **before,
        "dry_run": False,
        "deleted_count": len(deleted),
        "deleted": deleted,
        "reclaimed_bytes": reclaimed,
        "remaining_unique_bytes": after["unique_bytes"],
    }


__all__ = [
    "CACHE_TARGETS",
    "cache_inventory",
    "prune_caches",
]

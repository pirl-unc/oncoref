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

"""One catalog over every managed (fetchable) dataset.

oncoref holds two kinds of heavy/external data behind separate machinery: the
version-pinned per-cohort **expression bundle** (:mod:`oncoref.data_bundle`, a
single GitHub-release tarball) and versioned **reference** tables
(:mod:`oncoref.reference_data`: HPA, DICE and the human corneal atlas).
This module is the unified facade over both — one ``Dataset`` model and one
fetch/cache/status API — so a caller (or the CLI) manages everything uniformly
without knowing which backend a dataset lives in.

    from oncoref import catalog
    catalog.datasets()                       # every managed dataset
    catalog.status()                         # uniform present/size rows
    catalog.ensure("hpa_normal_tissue")      # download if needed -> Path
    catalog.ensure("cancer-reference-expression-percentiles")
    catalog.fetch("all")                     # materialize everything

Small curated tables also support path, version and provenance inspection.
They ship in the wheel and are listed by inventory(), not the fetch-all list.
All selectors accept version=None (newest registered) or an explicit pin;
historical curated tables can be read from a verified published wheel without
installing or executing its code.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from . import _downloads, data_bundle, data_manifest, packaged_data, reference_data, source_matrices
from .legacy import LEGACY_TYPED_ACCESSORS, legacy_dataset_dispositions
from .species import Species
from .version import __version__

_BUNDLE = "bundle"
_HPA = "hpa"
_REFERENCE = "reference"

#: prefix for addressing one cohort's raw per-sample matrix (``per-sample:LUAD``).
_PER_SAMPLE = "per-sample:"

#: Human-readable descriptions for the expression-bundle members.
_BUNDLE_DESCRIPTIONS = {
    "cancer-reference-expression": "per-cohort RNA-seq expression shards",
    "cancer-reference-expression-representatives": "per-cohort medoid/exemplar samples",
    "cancer-reference-expression-percentiles": "per-gene tail-weighted percentile vectors",
    "cancer-reference-expression-within-sample-top5": (
        "within-sample high-expression prevalence by cohort"
    ),
    "cancer-reference-expression-percentiles-proteoform-cta": (
        "CTA-scope proteoform-summed percentile vectors"
    ),
    "cancer-reference-expression-within-sample-top5-proteoform-cta": (
        "CTA-scope proteoform-summed within-sample prevalence"
    ),
    "source-matrix-sample-qc.csv": "per-sample QC manifest for rebuilt expression artifacts",
    "expression-artifact-build-metadata.csv": (
        "per-cohort build/QC provenance for rebuilt expression artifacts"
    ),
    "expression-artifact-build-metadata.json": (
        "bundle-level build/QC provenance for rebuilt expression artifacts"
    ),
    "pan-cancer-expression.csv": "pan-cancer HPA-tissue + TCGA expression matrix",
    "hpa-cell-type-expression.csv": "HPA cell-type nTPM matrix",
}


@dataclass(frozen=True)
class Dataset:
    """Dataset name, storage kind (bundle, hpa, reference or wheel), and description."""

    name: str
    kind: str
    description: str


def datasets() -> list[Dataset]:
    """Every managed dataset across both backends (bundle members + HPA sources)."""
    bundle_names = list(data_bundle.DOWNLOADABLE_PATHS)
    hpa_names = list(reference_data.REFERENCE_SOURCES)
    collisions = set(bundle_names) & set(hpa_names)
    if collisions:
        # The catalog routes by name; a name in both backends would silently
        # shadow one of them. Fail loudly instead of mis-routing.
        raise RuntimeError(f"dataset name collision across backends: {sorted(collisions)}")
    out = [
        Dataset(name, _BUNDLE, _BUNDLE_DESCRIPTIONS.get(name, "expression bundle artifact"))
        for name in bundle_names
    ]
    out += [
        Dataset(
            name,
            reference_data.REFERENCE_SOURCES[name].get("kind", _HPA),
            reference_data.REFERENCE_SOURCES[name]["description"],
        )
        for name in hpa_names
    ]
    return out


def _by_name() -> dict[str, Dataset]:
    return {d.name: d for d in datasets()}


def dataset(name: str) -> Dataset:
    """The :class:`Dataset` for ``name``; raises ``KeyError`` if unknown."""
    wheel = {**data_manifest.WHEEL, **data_manifest.CANCERDATA_ORIGINATED}
    if name in wheel and name not in _by_name():
        return Dataset(name, "wheel", wheel[name][1])
    try:
        return _by_name()[name]
    except KeyError:
        avail = ", ".join(sorted(_by_name()))
        raise KeyError(f"unknown dataset {name!r}; available: {avail}") from None


def _species_packaged(name: str, selected) -> tuple[Path, str] | None:
    """Installed native evidence snapshots have their own curated source version."""
    if name != "cta-evidence":
        return None
    import json

    dest = Path(__file__).parent / "data" / "species" / selected.name / "cta-evidence.json"
    if not dest.is_file():
        raise KeyError(f"No {selected.name} {name!r} dataset is registered")
    return dest, json.loads(dest.read_text())["source_version"]


def _packaged_versions(relative: str) -> tuple[str, ...]:
    """Installed package first, followed by locally cached historical wheels."""
    root = reference_data.cache_dir() / "packaged"
    cached = sorted(
        (p.relative_to(root).parts[0] for p in root.glob(f"oncoref-*/{relative}")), reverse=True
    )
    return tuple(dict.fromkeys((f"oncoref-{__version__}", *cached)))


def versions(name: str, *, species: Species | str | int = "human") -> tuple[str, ...]:
    """Known versions, newest registered default first; never queries the network.

    Upstream labels are retained (v11, study accession/date, etc.). Bundle and
    per-sample releases additionally accept explicit historical X.Y.Z versions;
    their availability and compatibility are checked when fetched.
    """
    from .species import resolve_species

    selected = resolve_species(species)
    if selected.taxon_id != 9606:
        packaged = _species_packaged(name, selected)
        if packaged:
            return _packaged_versions(f"species/{selected.name}/{packaged[0].name}")
        from .species_assets import asset_versions

        return asset_versions(name, species=selected)
    if name.startswith(_PER_SAMPLE):
        code = name[len(_PER_SAMPLE) :]
        current = source_matrices.resolve_version(code)
        cached = tuple(
            p.parent.name.removeprefix("v")
            for p in source_matrices.cache_dir().parent.glob(
                f"v*/{source_matrices._resolve(code)}.parquet"
            )
        )
        return tuple(dict.fromkeys((current, *sorted(cached, reverse=True))))
    d = dataset(name)
    if d.kind in (_HPA, _REFERENCE):
        return reference_data.versions(name)
    if d.kind == "wheel":
        return _packaged_versions(_wheel_path(name).name)
    cached = (row["version"].removeprefix("v") for row in data_bundle.list_cache_versions())
    return tuple(dict.fromkeys((data_bundle.DATA_VERSION, *cached)))


def resolve_version(
    name: str, version: str | None = None, *, species: Species | str | int = "human"
) -> str:
    """Return a concrete version; None/"latest" select the registered default.

    No numeric ordering is guessed across unrelated providers. The installed
    registry identifies the newest supported release for each source. An
    explicit unavailable pin never falls back to the default.
    """
    from .species import resolve_species

    selected = resolve_species(species)
    if selected.taxon_id != 9606:
        packaged = _species_packaged(name, selected)
        if packaged:
            return packaged_data.resolve_version(version)
        from .species_assets import species_asset

        return species_asset(name, species=selected, version=version).source_version
    if name.startswith(_PER_SAMPLE):
        return source_matrices.resolve_version(name[len(_PER_SAMPLE) :], version)
    d = dataset(name)
    if d.kind in (_HPA, _REFERENCE):
        return reference_data.resolve_version(name, version)
    if d.kind == "bundle":
        return data_bundle.resolve_version(version)
    return packaged_data.resolve_version(version)


def _wheel_path(name: str) -> Path:
    root = Path(__file__).parent / "data"
    candidates = [
        root / name,
        *(root / (name + suffix) for suffix in (".csv", ".csv.gz", ".json", ".yaml", ".yml")),
    ]
    matches = [path for path in candidates if path.is_file()]
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one packaged artifact for {name!r}, found {matches}")
    return matches[0]


def provenance(
    name: str,
    version: str | None = None,
    *,
    species: Species | str | int = "human",
    verify_content: bool = False,
) -> dict:
    """Source identity and concrete version with a content checksum when known.

    Version is the artifact release, not a gene-model version. Annotation,
    processing pipeline, units and upstream study identifiers remain separate
    metadata. Unverified/missing checksums are reported explicitly.
    """
    from .species import resolve_species

    selected = resolve_species(species)
    concrete = resolve_version(name, version, species=selected)
    if selected.taxon_id != 9606:
        packaged = _species_packaged(name, selected)
        if packaged:
            dest = path(name, concrete, species=selected)
            record = {
                "name": name,
                "version": concrete,
                "kind": "wheel",
                "taxon_id": selected.taxon_id,
                "path": str(dest) if dest else None,
                "exists": dest is not None,
                "sha256": None,
                "verification_state": "missing_file" if dest is None else "not_checked",
            }
            if concrete == f"oncoref-{__version__}":
                record.update(
                    sha256=reference_data._sha256(dest),
                    source_version=packaged[1],
                    verification_state="content_recorded",
                )
            else:
                relative = f"species/{selected.name}/{packaged[0].name}"
                record.update(packaged_data.provenance(relative, concrete))
                if dest and verify_content and record.get("sha256"):
                    matches = reference_data._sha256(dest) == record["sha256"]
                    record["checksum_matches"] = matches
                    record["verification_state"] = "verified" if matches else "checksum_mismatch"
                if dest:
                    import json

                    try:
                        record["source_version"] = json.loads(dest.read_text())["source_version"]
                    except (OSError, ValueError, KeyError):
                        record["source_version"] = None
            return record
        from .species_assets import asset_path, inspect_asset, species_asset

        asset = species_asset(name, species=selected, version=concrete)
        dest = asset_path(asset)
        state = inspect_asset(asset)["status"] if verify_content else "not_checked"
        return {
            "name": name,
            "version": concrete,
            "kind": "asset",
            "taxon_id": selected.taxon_id,
            "url": asset.url,
            "path": str(dest),
            "sha256": asset.sha256,
            "exists": dest.exists(),
            "verification_state": state,
        }
    if name in reference_data.REFERENCE_SOURCES:
        return {
            **reference_data.provenance(name, concrete, verify_content=verify_content),
            "kind": dataset(name).kind,
            "taxon_id": selected.taxon_id,
        }
    dest = path(name, version=concrete)
    kind = "source" if name.startswith(_PER_SAMPLE) else dataset(name).kind
    record = {
        "name": name,
        "version": concrete,
        "kind": kind,
        "taxon_id": selected.taxon_id,
        "path": str(dest) if dest else None,
        "exists": dest is not None,
        "sha256": None,
        "verification_state": "missing_file" if dest is None else "not_checked",
    }
    if kind == "source":
        record["url"] = source_matrices.release_url(name[len(_PER_SAMPLE) :], concrete)
    elif kind == "bundle":
        record["url"] = data_bundle.release_sources(concrete)[0]["url"]
        marker = data_bundle._read_completion_marker(data_bundle.cache_dir(concrete)) or {}
        record["archive_sha256"] = (
            (marker.get("release_manifest") or {}).get("tarball", {}).get("sha256")
        )
    if kind == "wheel" and concrete != f"oncoref-{__version__}":
        record.update(packaged_data.provenance(_wheel_path(name).name, concrete))
        if dest and verify_content and record.get("sha256"):
            matches = reference_data._sha256(dest) == record["sha256"]
            record["checksum_matches"] = matches
            record["verification_state"] = "verified" if matches else "checksum_mismatch"
    elif dest and dest.is_file() and (verify_content or kind == "wheel"):
        record["sha256"] = reference_data._sha256(dest)
        # Observing bytes is not the same as checking against an independent pin.
        record["verification_state"] = "content_recorded"
    return record


def path(
    name: str, version: str | None = None, *, species: Species | str | int = "human"
) -> Path | None:
    """Find exactly the selected artifact without downloading or substituting versions."""
    from .species import resolve_species

    selected = resolve_species(species)
    if selected.taxon_id != 9606:
        packaged = _species_packaged(name, selected)
        if packaged:
            concrete = packaged_data.resolve_version(version)
            if concrete == f"oncoref-{__version__}":
                return packaged[0]
            dest = packaged_data.local_path(f"species/{selected.name}/{packaged[0].name}", concrete)
            return dest if dest.is_file() else None
        from .species_assets import asset_path, species_asset

        dest = asset_path(species_asset(name, species=selected, version=version))
        return dest if dest.is_file() else None
    if name.startswith(_PER_SAMPLE):
        dest = source_matrices.local_path(name[len(_PER_SAMPLE) :], version)
        return dest if dest.exists() else None
    d = dataset(name)
    if d.kind in (_HPA, _REFERENCE):
        dest = reference_data.local_path(name, version)
        return dest if dest.exists() else None
    if d.kind == "wheel":
        concrete = resolve_version(name, version)
        if concrete == f"oncoref-{__version__}":
            return _wheel_path(name)
        dest = packaged_data.local_path(_wheel_path(name).name, concrete)
        return dest if dest.is_file() else None
    return data_bundle.find(name, **({"version": version} if version is not None else {}))


def _per_sample_status_row(name: str, version: str | None = None) -> dict:
    code = name[len(_PER_SAMPLE) :]
    try:
        info = source_matrices.cohort_info(code)
        p = path(name, version)
    except source_matrices.SourceMatrixError as e:
        raise KeyError(str(e)) from None
    resolved = str(info["cancer_code"])
    concrete = resolve_version(name, version)
    sample_note = (
        f"; n={int(info['n_samples'])}"
        if concrete == source_matrices.source_matrix_version(code)
        else "; historical build"
    )
    return {
        "name": f"{_PER_SAMPLE}{resolved}",
        "kind": "source",
        "version": resolve_version(name, version),
        "present": p is not None,
        "path": str(p) if p is not None else None,
        "size_bytes": _size_bytes(p),
        "cohorts": 1,
        "description": (
            f"raw per-sample TPM matrix for {resolved} ({info['source_cohort']}{sample_note})"
        ),
    }


def ensure(
    name: str, version: str | None = None, *, species: Species | str | int = "human"
) -> Path:
    """Local path to ``name``, downloading if absent. Returns the path.

    HPA sources fetch per-file; a bundle member triggers the one-time tarball
    download (all bundle members arrive together); ``per-sample:<code>`` fetches
    one cohort's raw matrix.
    """
    from .species import resolve_species

    selected = resolve_species(species)
    if selected.taxon_id != 9606:
        packaged = _species_packaged(name, selected)
        if packaged:
            concrete = packaged_data.resolve_version(version)
            return (
                packaged[0]
                if concrete == f"oncoref-{__version__}"
                else packaged_data.ensure(f"species/{selected.name}/{packaged[0].name}", concrete)
            )
        from .species_assets import ensure_asset, species_asset

        return ensure_asset(species_asset(name, species=selected, version=version))
    if name.startswith(_PER_SAMPLE):
        return source_matrices.ensure(
            name[len(_PER_SAMPLE) :], **({"version": version} if version is not None else {})
        )
    d = dataset(name)
    if d.kind in (_HPA, _REFERENCE):
        return reference_data.ensure(name, version)
    if d.kind == "wheel":
        concrete = resolve_version(name, version)
        return (
            _wheel_path(name)
            if concrete == f"oncoref-{__version__}"
            else packaged_data.ensure(_wheel_path(name).name, concrete)
        )
    data_bundle.ensure_local(**({"version": version} if version is not None else {}))
    resolved = data_bundle.find(name, **({"version": version} if version is not None else {}))
    if resolved is None:
        raise FileNotFoundError(
            f"{name!r} missing from the expression bundle after fetch ({data_bundle.cache_dir()})"
        )
    return resolved


def fetch(
    target: str = "all",
    *,
    version: str | None = None,
    species: Species | str | int = "human",
    force: bool = False,
) -> list[str]:
    """Materialize data, downloading what's missing. ``target`` selects what:

    - ``"all"`` — the expression bundle and every managed reference source;
    - ``"bundle"`` / ``"hpa"`` — just that kind (e.g. ``fetch("hpa")`` for HPA only);
    - ``"source"`` — every cohort's raw per-sample matrix (~21 GB; opt-in, not in ``"all"``);
    - ``"per-sample:<code>"`` — one cohort's raw matrix;
    - a dataset name — just that dataset.

    Returns the names actually downloaded (already-cached skipped unless ``force``).
    The bundle tarball is fetched at most once.
    """
    from .species import resolve_species

    selected = resolve_species(species)
    if selected.taxon_id != 9606:
        from .species_assets import species_assets

        if force:
            raise ValueError(
                "Pinned species assets are verified and reused; force refresh is unsupported"
            )
        names = [a.id for a in species_assets(selected)] if target == "all" else [target]
        if version is not None and target == "all":
            raise ValueError("Select one dataset when pinning a version")
        downloaded = []
        for name in names:
            before = path(name, version, species=selected)
            ensure(name, version, species=selected)
            if before is None:
                downloaded.append(name)
        return downloaded
    if version is not None and target in {"all", "source", "hpa", "reference"}:
        raise ValueError(
            "Select one dataset when pinning a version; providers have independent releases"
        )
    options = {"version": version} if version is not None else {}
    # one cohort's raw per-sample matrix
    if target.startswith(_PER_SAMPLE):
        code = target[len(_PER_SAMPLE) :]
        already = source_matrices.is_cached(code, **options)
        source_matrices.fetch(code, force=force, **options)
        return [target] if force or not already else []
    # every per-sample matrix (heavy; deliberately not part of "all")
    if target == "source":
        out = []
        for code in source_matrices.available_cohorts():
            already = source_matrices.is_cached(code, **options)
            source_matrices.fetch(code, force=force, **options)
            if force or not already:
                out.append(f"{_PER_SAMPLE}{code}")
        return out

    if target == "all":
        targets = [d.name for d in datasets()]
    elif target in (_BUNDLE, _HPA, _REFERENCE):
        targets = [d.name for d in datasets() if d.kind == target]
    else:
        targets = [dataset(target).name]

    if len(targets) == 1 and dataset(targets[0]).kind == "wheel":
        before = path(targets[0], version)
        ensure(targets[0], version)
        return [targets[0]] if before is None else []
    downloaded = []
    # The bundle is one tarball — fetch it a single time if any member is targeted.
    bundle_targets = [n for n in targets if dataset(n).kind == _BUNDLE]
    if bundle_targets and (force or not data_bundle.is_local(**options)):
        data_bundle.fetch(**options)
        downloaded += bundle_targets
    elif bundle_targets:
        _downloads.report_cached(
            f"expression bundle v{data_bundle.resolve_version(version)}",
            data_bundle.cache_dir(version),
        )
    for n in targets:
        if dataset(n).kind in (_HPA, _REFERENCE):
            dest = reference_data.local_path(n, version)
            concrete = reference_data.resolve_version(n, version)
            already = dest.exists() and reference_data._cached_file_ok(n, concrete, dest)
            reference_data.download(n, version, force=force)
            if force or not already:
                downloaded.append(n)
    return downloaded


def _size_bytes(p: Path | None) -> int:
    if p is None or not p.exists():
        return 0
    if p.is_dir():
        return sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
    return p.stat().st_size


def _cohort_count(p: Path | None) -> int | None:
    """Per-cohort file count for a directory dataset (one file per cancer cohort);
    ``None`` for a single-file dataset. This is the cohort scale held *inside* an
    expression artifact — e.g. ``cancer-reference-expression-percentiles`` is one
    catalog entry but ~118 cancer cohorts."""
    if p is None or not p.is_dir():
        return None
    return sum(1 for f in p.iterdir() if f.is_file() and not f.name.startswith("_"))


def inventory(*, species: Species | str | int = "human", group=None) -> list[dict]:
    """The complete oncoref-domain data inventory — the full picture behind the
    fetchable :func:`datasets`. One row per dataset with ``name``, ``held``
    (``wheel`` / ``bundle`` / ``hpa`` / ``source`` / ``planned``), ``category``,
    ``available`` (present locally / shipped), and ``description``. Legacy rows
    additionally expose ``current_status``, ``replacement_owner``,
    ``replacement_surface``, ``compatibility_policy``, and ``typed_accessor``;
    those fields are ``None`` for current datasets. Driven by
    :mod:`oncoref.data_manifest`; this is the current self-owned inventory.
    The separately versioned pirlygenes migration snapshot is historical audit
    evidence and does not expand when downstream adds a new dataset.
    """
    from .species import DATA_GROUPS, _human_data_group, resolve_species

    selected = resolve_species(species)
    if group is not None and group not in DATA_GROUPS:
        raise ValueError(f"Unknown data group {group!r}")
    if selected.taxon_id != 9606:
        from .species_antigens import antigen_inventory
        from .species_assets import asset_inventory

        rows = asset_inventory(selected, group=group)
        if group in (None, "antigens"):
            rows += antigen_inventory(selected)
        return rows
    bundle_member = {p.removesuffix(".csv"): p for p in data_bundle.DOWNLOADABLE_PATHS}
    rows: list[dict] = []
    legacy_rows = {
        str(row.dataset): row._asdict()
        for row in legacy_dataset_dispositions().itertuples(index=False)
    }

    def _add(name, held, category, description, available, cohorts=None):
        disposition = legacy_rows.get(name, {})
        rows.append(
            {
                "name": name,
                "held": held,
                "category": category,
                "available": available,
                "cohorts": cohorts,
                "description": description,
                "current_status": disposition.get("current_status"),
                "replacement_owner": disposition.get("replacement_owner"),
                "replacement_surface": disposition.get("replacement_surface"),
                "compatibility_policy": disposition.get("compatibility_policy"),
                "typed_accessor": LEGACY_TYPED_ACCESSORS.get(name),
            }
        )

    wheel = {**data_manifest.WHEEL, **data_manifest.CANCERDATA_ORIGINATED}
    for name, (cat, desc) in sorted(wheel.items()):
        _add(name, "wheel", cat, desc, True)  # ships in the wheel — always present
    for name, (cat, desc) in sorted(data_manifest.BUNDLE.items()):
        member = data_bundle.find(bundle_member[name])
        _add(name, "bundle", cat, desc, member is not None, _cohort_count(member))
    for name, (cat, desc) in sorted(data_manifest.HPA.items()):
        _add(name, "hpa", cat, desc, reference_data.local_path(name).exists())
    for name, spec in sorted(reference_data._normal_manifest["artifacts"].items()):
        _add(
            name,
            spec["kind"],
            "expression",
            spec["description"],
            reference_data.local_path(name).exists(),
        )
    for name, spec in sorted(reference_data._annotation_manifest["artifacts"].items()):
        _add(
            name,
            "reference",
            "protein-annotation",
            spec["description"],
            reference_data.local_path(name).exists(),
        )
    listed = {row["name"] for row in rows}
    for name, spec in sorted(reference_data.REFERENCE_SOURCES.items()):
        if name not in listed:
            _add(
                name,
                spec.get("kind", "hpa"),
                "expression",
                spec["description"],
                reference_data.local_path(name).exists(),
            )
    for name, (cat, desc) in sorted(data_manifest.SOURCE.items()):
        cohorts = source_matrices.available_cohorts()
        available = any(source_matrices.is_cached(c) for c in cohorts)
        _add(name, "source", cat, desc, available, len(cohorts))  # per-cohort fetchable
    for name, (cat, desc) in sorted(data_manifest.PLANNED.items()):
        _add(name, "planned", cat, desc, False)  # oncoref-domain, still to port
    for row in rows:
        catalog_name = bundle_member.get(row["name"], row["name"])
        if row["held"] in {"wheel", "bundle", "hpa", "reference"}:
            row["default_version"] = resolve_version(catalog_name)
        else:
            row["default_version"] = None
        row.update(taxon_id=selected.taxon_id, species=selected.name, group=_human_data_group(row))
    return [row for row in rows if group is None or row["group"] == group]


def status(
    name: str | None = None, version: str | None = None, *, species: Species | str | int = "human"
) -> list[dict]:
    """Uniform status rows over the catalog: ``name``, ``kind``, ``present``,
    ``path``, ``size_bytes``, ``cohorts`` (per-cohort file count for directory
    datasets, else ``None``), ``description``. One dataset if ``name`` given."""
    from .species import resolve_species

    selected = resolve_species(species)
    if selected.taxon_id != 9606:
        from .species_assets import species_assets

        names = [a.id for a in species_assets(selected)] if name in (None, "all") else [name]
        if version is not None and len(names) != 1:
            raise ValueError("Select one dataset when pinning a version")
        rows = []
        for n in names:
            record = provenance(n, version, species=selected)
            dest = path(n, version, species=selected)
            rows.append(
                {
                    "name": n,
                    "kind": record["kind"],
                    "version": record["version"],
                    "present": dest is not None,
                    "path": str(dest) if dest else None,
                    "size_bytes": _size_bytes(dest),
                    "cohorts": None,
                    "description": n,
                }
            )
        return rows
    if version is not None and name in (None, "all", "source", "hpa", "reference"):
        raise ValueError("Select one dataset when pinning a version")
    if name == "all":
        name = None
    if name is not None and name.startswith(_PER_SAMPLE):
        return [_per_sample_status_row(name, version)]
    if name == "source":
        return [
            _per_sample_status_row(f"{_PER_SAMPLE}{code}")
            for code in source_matrices.available_cohorts()
        ]
    if name in (_BUNDLE, _HPA, _REFERENCE):
        names = [d.name for d in datasets() if d.kind == name]
    else:
        names = [dataset(name).name] if name is not None else [d.name for d in datasets()]
    rows = []
    for n in names:
        d = dataset(n)
        p = path(n, version)
        rows.append(
            {
                "name": n,
                "kind": d.kind,
                "version": resolve_version(n, version),
                "present": p is not None,
                "path": str(p) if p is not None else None,
                "size_bytes": _size_bytes(p),
                "cohorts": _cohort_count(p),
                "description": d.description,
            }
        )
    return rows

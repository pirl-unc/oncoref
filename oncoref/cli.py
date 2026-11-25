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

"""``oncoref`` command-line interface.

Reference lookups over the bundled tables (cancer-type / TMB / burden / ICI) plus
the ``data`` manager for downloadable expression-bundle, HPA, and per-sample sources.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import textwrap

from . import (
    apd1,
    cancer_types,
    cta,
    data_bundle,
    expression_registry,
    ici,
    incidence,
    proteoforms,
    reference_data,
    samples,
    tmb,
)
from .version import __version__


def _fmt_bytes(n: int) -> str:
    size = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


def _listing_width():
    return max(20, min(120, shutil.get_terminal_size(fallback=(100, 24)).columns))


def _listing_text(value, *, indent=""):
    print(
        textwrap.fill(
            str(value), width=_listing_width(), initial_indent=indent, subsequent_indent=indent
        )
    )


def _listing_table(headers, records):
    """Wrap details; use labeled records when the table would exceed the terminal."""
    widths = [
        max(len(header), *(len(cells[i]) for cells, _ in records))
        for i, header in enumerate(headers)
    ]
    table_fits = sum(widths) + 2 * (len(headers) - 1) <= _listing_width()
    if table_fits:
        print("  ".join(header.ljust(width) for header, width in zip(headers, widths)).rstrip())
        print("  ".join("-" * width for width in widths))
    for cells, details in records:
        if table_fits:
            print("  ".join(cell.ljust(width) for cell, width in zip(cells, widths)).rstrip())
        else:
            _listing_text(cells[0])
            _listing_text(
                " | ".join(f"{header}: {cell}" for header, cell in zip(headers[1:], cells[1:])),
                indent="  ",
            )
        for detail in details:
            _listing_text(detail, indent="  ")


def _listing_title(species, title):
    name = "Canine" if species.taxon_id == 9615 else species.name.title()
    _listing_text(f"{name} {title}")
    _listing_text(f"{species.scientific_name} | NCBI taxon {species.taxon_id}")


def _species_asset_options(args):
    from pathlib import Path

    overrides = json.loads(Path(args.asset_overrides).read_text()) if args.asset_overrides else None
    if overrides is not None and not isinstance(overrides, dict):
        raise ValueError("Asset overrides must be a JSON object mapping taxon:asset IDs to paths")
    return {"cache_root": args.cache_root, "local_overrides": overrides}


def _cmd_data_list(args, selected):
    from . import catalog, source_matrices
    from .species import DATA_GROUPS

    if args.group is not None and args.group not in DATA_GROUPS:
        raise ValueError(f"Unknown data group {args.group!r}")
    if args.cache_root or args.asset_overrides:
        from . import species_assets as assets
        from .species_antigens import antigen_inventory

        rows = assets.asset_inventory(selected, group=args.group, **_species_asset_options(args))
        if selected.taxon_id == 9606:
            # Keep the historical inventory while applying explicit cache options
            # only to registered species assets.
            asset_names = {a.id for a in assets.species_assets(selected)}
            rows += [
                r
                for r in catalog.inventory(species=selected, group=args.group)
                if r["name"] not in asset_names
            ]
        elif args.group in (None, "antigens"):
            rows += antigen_inventory(selected)
    else:
        rows = catalog.inventory(species=selected, group=args.group)
    if args.name and args.name != "all":
        rows = [r for r in rows if args.name in (r["held"], r["category"], r["name"])]
        if not rows:
            print(f"Error: unknown data list filter {args.name!r}", file=sys.stderr)
            return 1
    if args.json:
        print(json.dumps(rows, indent=2))
        return 0
    _listing_title(selected, "data")
    _listing_text(
        f"{len(rows)} datasets | {sum(r['held'] == 'wheel' for r in rows)} included in the package"
    )
    if not rows:
        _listing_text("No datasets registered for this selection.")
        return 0
    for group in DATA_GROUPS:
        entries = sorted((r for r in rows if r["group"] == group), key=lambda r: r["name"])
        if not entries:
            continue
        print()
        _listing_text(f"{group.title()} ({len(entries)})")
        records = []
        for row in entries:
            local = (
                "bundled"
                if row["held"] == "wheel"
                else "present"
                if row["available"]
                else "missing"
            )
            details = [row["description"]]
            if row["held"] == "source":
                cohorts = source_matrices.available_cohorts()
                count = sum(source_matrices.is_cached(code) for code in cohorts)
                local = f"{count}/{len(cohorts)}"
                details.append(f"{count} of {len(cohorts)} cohort matrices present locally.")
            elif row["cohorts"] is not None:
                details.append(f"{row['cohorts']} cohort files present locally.")
            if row.get("size_bytes") is not None:
                details.append(f"Asset size: {_fmt_bytes(row['size_bytes'])}")
            records.append(((row["name"], row["held"], local), details))
        _listing_table(("Dataset", "Storage", "Local"), records)
    print()
    _listing_text(
        "Local reports file presence, not checksum verification or complete biological coverage."
    )
    _listing_text("Use --json for the full inventory records.")
    return 0


def _cmd_version(args: argparse.Namespace) -> int:
    print(f"oncoref v{__version__}")
    return 0


def _print_cache_prune(args: argparse.Namespace) -> int:
    from . import cache_retention

    report = cache_retention.prune_caches(
        targets=args.name or "bundle",
        keep_current=not args.include_current,
        keep_versions=args.keep_versions,
        include_legacy_pirlygenes=args.include_legacy_pirlygenes,
        dry_run=not args.yes,
    )
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    rows = [e for e in report["entries"] if e["prunable"]] if not args.yes else report["deleted"]
    if not rows:
        print("Nothing to prune.")
        return 0
    verb = "Would delete" if not args.yes else "Deleted"
    for entry in rows:
        print(
            f"{verb} [{entry['target']}] {entry['id']}  "
            f"(logical {_fmt_bytes(entry['logical_bytes'])}, "
            f"reclaimable {_fmt_bytes(entry['reclaimable_bytes'])})  {entry['path']}"
        )
    reclaimed = report["reclaimable_bytes"] if not args.yes else report["reclaimed_bytes"]
    print(
        f"\nTotal physical space {'reclaimable' if not args.yes else 'reclaimed'}: {_fmt_bytes(reclaimed)}"
    )
    if not args.yes:
        print("(dry run — pass --yes to delete)")
    return 0


def _cmd_species(args):
    from .species import available_species, for_species

    if args.name is None:
        rows = [
            {"species": s.name, "scientific_name": s.scientific_name, "taxon_id": s.taxon_id}
            for s in available_species()
        ]
    else:
        rows = [group.coverage() for group in for_species(args.name).groups]
    print(json.dumps(rows, indent=2))
    return 0


def _cmd_species_data(args):
    from . import species_assets as assets

    options = _species_asset_options(args)
    target = args.name or "all"
    if args.action in {"dir", "dirs"}:
        d = assets.asset_cache_dir(cache_root=args.cache_root)
        if args.json:
            print(json.dumps({"species": str(d)}, indent=2))
        else:
            print(d)
        return 0
    selected = assets.species_assets(args.species)
    if target != "all":
        selected = tuple(a for a in selected if a.id == target or a.group == target)
    if args.data_version is not None:
        if len(selected) != 1:
            raise ValueError("Select one asset when pinning a version")
        selected = (
            assets.species_asset(selected[0].id, species=args.species, version=args.data_version),
        )
    if not selected:
        print(
            f"Error: unavailable {args.species} asset/group {target!r}; no human fallback",
            file=sys.stderr,
        )
        return 1
    if args.action == "status":
        print(json.dumps([assets.inspect_asset(a, **options) for a in selected], indent=2))
    elif args.action in {"metadata", "info"}:
        rows = [
            {
                **a.as_dict(),
                "inspection": assets.inspect_asset(a, **options),
                "available_versions": list(assets.asset_versions(a.id, species=args.species)),
            }
            if args.action == "info"
            else a.as_dict()
            for a in selected
        ]
        print(json.dumps(rows[0] if args.action == "info" and len(rows) == 1 else rows, indent=2))
    elif args.action == "verify":
        inspected = [assets.inspect_asset(a, **options) for a in selected]
        ok = all(item["status"] == "available" for item in inspected)
        if args.json:
            print(json.dumps({"verified": ok, "assets": inspected}, indent=2))
        else:
            for item in inspected:
                print(f"{item['asset_id']}\t{item['status']}\t{item['path']}")
        return 0 if ok else 1
    elif args.action == "path":
        if len(selected) != 1:
            raise ValueError("Specify one asset for data path")
        state = assets.inspect_asset(selected[0], **options)
        if state["status"] != "available":
            raise ValueError(f"{selected[0].id}: {state['status']}")
        if args.json:
            print(json.dumps({"asset_id": selected[0].id, "path": state["path"]}, indent=2))
        else:
            print(state["path"])
    elif args.action == "fetch":
        if args.force:
            raise ValueError(
                "Pinned species assets are verified and reused; force refresh is unsupported"
            )
        missing = sum(
            a.size_bytes
            for a in selected
            if assets.inspect_asset(a, **options)["status"] != "available"
        )
        if not args.offline and missing > args.max_bytes:
            raise ValueError(f"Transfer requires {missing} bytes; budget is {args.max_bytes}")
        fetched_paths = [
            str(
                assets.ensure_asset(
                    asset, offline=args.offline, max_bytes=args.max_bytes, **options
                )
            )
            for asset in selected
        ]
        if args.json:
            print(json.dumps({"downloaded": fetched_paths}, indent=2))
        else:
            for p in fetched_paths:
                print(p)
    else:
        raise ValueError(f"data {args.action} is unavailable for {args.species}; no human fallback")
    return 0


def _resolve_cache_dirs() -> dict[str, str]:
    import os
    from pathlib import Path

    from . import _downloads, sequence_bundle, source_matrices, species_assets

    shard_env = os.environ.get("ONCOREF_SHARD_CACHE_DIR", "").strip()
    shard_dir = (
        Path(shard_env).expanduser()
        if shard_env
        else _downloads.oncoref_cache_root() / "shard_cache"
    )
    return {
        "root": str(_downloads.oncoref_cache_root()),
        "bundle": str(data_bundle.cache_dir()),
        "hpa": str(reference_data.cache_base_dir() / "sources"),
        "reference": str(reference_data.cache_base_dir() / "sources"),
        "source": str(source_matrices.cache_root() / f"v{source_matrices.SOURCE_MATRIX_VERSION}"),
        "sequence": str(sequence_bundle.cache_root()),
        "species": str(species_assets.asset_cache_dir()),
        "shard": str(shard_dir),
    }


def _cmd_data(args: argparse.Namespace) -> int:
    from ._downloads import download_options

    with download_options(
        progress=args.progress,
        verbose=not args.quiet,
        offline=True if args.offline else None,
        cache_root=args.cache_root,
    ):
        try:
            return _cmd_data_impl(args)
        except (ImportError, FileNotFoundError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1


def _cmd_data_impl(args: argparse.Namespace) -> int:
    """Unified OpenVax data management over every managed OncoRef dataset (#747)."""
    from . import _downloads, catalog
    from .species import resolve_species

    selected = resolve_species(args.species)
    if args.action in {"versions", "provenance"}:
        if not args.name:
            raise ValueError(f"data {args.action} requires a dataset name")
        result = (
            catalog.versions(args.name, species=selected)
            if args.action == "versions"
            else catalog.provenance(args.name, args.data_version, species=selected)
        )
        print(json.dumps(result, indent=2))
        return 0
    if args.data_version is not None and args.action not in {
        "status",
        "info",
        "verify",
        "fetch",
        "path",
    }:
        raise ValueError(
            "--version is supported for status, info, verify, fetch, path, and provenance"
        )
    from .species_assets import species_assets

    asset_target = any(args.name in (a.id, a.group) for a in species_assets(selected))
    if args.action not in {"list", "prune", "dirs"} and (selected.taxon_id != 9606 or asset_target):
        if args.action == "dir" and selected.taxon_id != 9606:
            return _cmd_species_data(args)
        return _cmd_species_data(args)
    if args.action == "list":
        return _cmd_data_list(args, selected)
    if args.asset_overrides and not asset_target and selected.taxon_id == 9606:
        raise ValueError("--asset-overrides requires a registered species asset")

    if args.action == "status":
        try:
            rows = catalog.status(
                args.name,
                **({"version": args.data_version} if args.data_version is not None else {}),
            )
        except KeyError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
        if args.json:
            print(json.dumps(rows, indent=2, default=str))
            return 0
        print(
            f"{'Dataset':<46} {'Kind':<7} {'Present':<8} {'Size':>10} {'Cohorts':>8}  Description"
        )
        print("-" * 110)
        for r in rows:
            present = "yes" if r["present"] else "no"
            size = _fmt_bytes(r["size_bytes"])
            cohorts = str(r["cohorts"]) if r["cohorts"] is not None else "-"
            print(
                f"{r['name']:<46} {r['kind']:<7} {present:<8} {size:>10} {cohorts:>8}  "
                f"{r['description']}"
            )
        return 0

    if args.action == "info":
        if args.name in (None, "bundle", "oncoref", "pirlygenes"):
            source = args.name if args.name in {"oncoref", "pirlygenes"} else "oncoref"
            try:
                payload = data_bundle.bundle_metadata(
                    source,
                    include_release_manifest=bool(getattr(args, "remote", False)),
                )
            except (ValueError, data_bundle.BundleIntegrityError, OSError) as e:
                print(f"Error: {e}", file=sys.stderr)
                return 1
            if args.json:
                print(json.dumps(payload, indent=2, sort_keys=True, default=str))
                return 0
            local_cache = payload["local_cache"]
            print(f"name: bundle ({source})")
            print(f"package_version: {payload['package_version']}")
            print(f"data_version: {payload['data_version']}")
            print(f"source_matrix_version: {payload['source_matrix_version']}")
            print(f"cache_dir: {local_cache['cache_dir']}")
            print(f"all_local: {local_cache['all_local']}")
            print(f"verified_sha256: {local_cache['completion_marker']['verified_sha256']}")
            print("\nUse --json (and --remote for release manifest) for complete metadata.")
            return 0
        try:
            prov = catalog.provenance(args.name, args.data_version, species=selected)
            vers = list(catalog.versions(args.name, species=selected))
        except (KeyError, ValueError, reference_data.ReferenceDataError) as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
        info_payload = {**prov, "available_versions": vers}
        if args.json:
            print(json.dumps(info_payload, indent=2, sort_keys=True, default=str))
            return 0
        for key, value in info_payload.items():
            print(f"{key}: {value}")
        return 0

    if args.action == "verify":
        target = args.name or "all"
        if target == "bundle":
            try:
                snap = data_bundle.verify_local()
            except (FileNotFoundError, RuntimeError) as e:
                if args.json:
                    print(
                        json.dumps(
                            {"target": "bundle", "verified": False, "error": str(e)}, indent=2
                        )
                    )
                else:
                    print(f"Error: {e}", file=sys.stderr)
                return 1
            if args.json:
                print(
                    json.dumps(
                        {
                            "target": "bundle",
                            "verified": True,
                            "data_version": snap["data_version"],
                            "cache_dir": snap["cache_dir"],
                        },
                        indent=2,
                    )
                )
            else:
                print(f"Verified bundle v{snap['data_version']} at {snap['cache_dir']}")
            return 0
        try:
            status_rows = catalog.status(
                None if target == "all" else target,
                **({"version": args.data_version} if args.data_version is not None else {}),
                species=selected,
            )
        except KeyError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
        checked = []
        all_ok = True
        for row in status_rows:
            if target == "all" and not row["present"]:
                continue
            if row["kind"] == "bundle":
                try:
                    data_bundle.verify_local()
                    vstate, ok = "verified", True
                except (FileNotFoundError, RuntimeError) as exc:
                    vstate, ok = f"invalid: {exc}", False
            else:
                prov = catalog.provenance(
                    row["name"],
                    args.data_version,
                    species=selected,
                    verify_content=True,
                )
                vstate = str(prov.get("verification_state") or "unknown")
                ok = vstate in {"verified", "content_recorded"} and bool(prov.get("exists"))
            all_ok = all_ok and ok
            checked.append(
                {
                    "name": row["name"],
                    "kind": row["kind"],
                    "version": row["version"],
                    "path": row["path"],
                    "verification_state": vstate,
                    "verified": ok,
                }
            )
        if args.json:
            print(json.dumps({"target": target, "verified": all_ok, "items": checked}, indent=2))
        elif not checked:
            print("No cached datasets present to verify.")
        else:
            for item in checked:
                print(f"{item['name']}\t{item['verification_state']}\t{item['path'] or '-'}")
        return 0 if all_ok else 1

    if args.action == "contract":
        print(json.dumps(data_bundle.bundle_contract(), indent=2, sort_keys=True, default=str))
        return 0

    if args.action == "metadata":
        source = args.name or "oncoref"
        try:
            metadata = data_bundle.bundle_metadata(source)
        except (ValueError, data_bundle.BundleIntegrityError, OSError) as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
        print(json.dumps(metadata, indent=2, sort_keys=True, default=str))
        return 0

    if args.action == "release-manifest":
        source = args.name or "oncoref"
        try:
            manifest = data_bundle.bundle_release_manifest(source)
        except (ValueError, data_bundle.BundleIntegrityError, OSError) as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
        print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
        return 0

    if args.action in {"dir", "dirs"}:
        target = args.name or "all"
        dirs = _resolve_cache_dirs()
        if target == "all":
            if args.json:
                print(json.dumps(dirs, indent=2, sort_keys=True))
            else:
                for name, path in dirs.items():
                    print(f"{name}\t{path}")
            return 0
        if target not in dirs:
            print(
                f"Error: data {args.action} target must be one of: all, {', '.join(dirs)}",
                file=sys.stderr,
            )
            return 1
        if args.json:
            print(json.dumps({target: dirs[target]}, indent=2))
        else:
            print(dirs[target])
        return 0

    if args.action == "fetch":
        try:
            downloaded = catalog.fetch(
                args.name or "all",
                force=args.force,
                **({"version": args.data_version} if args.data_version is not None else {}),
            )
        except KeyError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
        except Exception as e:  # network / offline / extract failure
            print(f"Error: {e}", file=sys.stderr)
            return 1
        if args.json:
            print(
                json.dumps(
                    {"target": args.name or "all", "downloaded": downloaded},
                    indent=2,
                )
            )
        else:
            print(f"Downloaded: {', '.join(downloaded)}" if downloaded else "Already present.")
        return 0

    if args.action == "path":
        if not args.name:
            print("Error: 'data path' requires a dataset name", file=sys.stderr)
            return 1
        try:
            if _downloads.is_offline():
                resolved_path = catalog.path(
                    args.name,
                    args.data_version,
                    species=selected,
                )
                if resolved_path is None:
                    raise FileNotFoundError(
                        f"Offline mode is enabled and {args.name!r} is not cached locally"
                    )
            else:
                resolved_path = catalog.ensure(
                    args.name,
                    **({"version": args.data_version} if args.data_version is not None else {}),
                )
        except KeyError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
        except Exception as e:  # network / offline / extract failure
            print(f"Error: {e}", file=sys.stderr)
            return 1
        if args.json:
            print(json.dumps({"name": args.name, "path": str(resolved_path)}, indent=2))
        else:
            print(resolved_path)
        return 0

    if args.action == "prune":
        return _print_cache_prune(args)

    return 1


def _cmd_cancer_type(args: argparse.Namespace) -> int:
    try:
        info = cancer_types.cancer_type_info(args.query)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    if info is None:
        print("Error: no cancer type given", file=sys.stderr)
        return 1
    print(json.dumps(info, indent=2, default=str))
    return 0


def _cmd_tmb(args: argparse.Namespace) -> int:
    if args.code is None:
        for code, value in sorted(tmb.cancer_tmb().items()):
            print(f"{code}\t{value:g}")
        return 0
    try:
        value = tmb.cancer_tmb(args.code, inherit=not args.no_inherit)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    if value is None:
        print(f"No TMB value for {args.code!r}", file=sys.stderr)
        return 1
    print(f"{value:g}")
    return 0


def _cmd_apd1(args: argparse.Namespace) -> int:
    if args.code is None:
        for code, value in sorted(apd1.cancer_apd1_response().items()):
            print(f"{code}\t{value:g}")
        return 0
    try:
        value = apd1.cancer_apd1_response(args.code, inherit=not args.no_inherit)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    if value is None:
        print(f"No anti-PD-1 ORR for {args.code!r}", file=sys.stderr)
        return 1
    print(f"{value:g}")
    return 0


def _cmd_ici(args: argparse.Namespace) -> int:
    regimen = args.regimen
    if args.code is None:
        mapping = ici.cancer_ici_response(regimen=regimen)
        for code, value in sorted(mapping.items()):
            print(f"{code}\t{value:g}")
        return 0
    try:
        if args.all_regimens:
            per = ici.cancer_ici_response(args.code, fallback=False, inherit=not args.no_inherit)
            if not per:
                print(
                    _no_ici_orr_message(args.code, inherit=not args.no_inherit),
                    file=sys.stderr,
                )
                return 1
            for reg in ici.ici_regimens():
                if reg in per:
                    print(f"{reg}\t{per[reg]:g}")
            return 0
        value = ici.cancer_ici_response(args.code, regimen=regimen, inherit=not args.no_inherit)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    if value is None:
        print(_no_ici_orr_message(args.code, inherit=not args.no_inherit), file=sys.stderr)
        return 1
    print(f"{value:g}")
    return 0


def _no_ici_orr_message(code: str, *, inherit: bool) -> str:
    """Explain *why* there is no ORR when the code is a reviewed gap.

    An audited gap and an uncurated code both resolve to no value; only the gap can
    say why, and the reason is the whole point of curating it."""
    resolved = ici.resolve_ici_response_source(code, inherit=inherit)
    if resolved.get("inheritance_kind") != "direct_missing":
        return f"No ICI ORR for {code!r}"
    reason = resolved.get("missing_reason") or "no_curated_representative_orr"
    message = f"No representative ICI ORR for {code!r}: {reason}"
    note = resolved.get("notes")
    return f"{message}\n{note}" if note else message


def _cmd_tcr_trials(args: argparse.Namespace) -> int:
    """Export the reviewed trial evidence without downloading source papers."""
    from . import tcr_trials

    if args.table == "all":
        print(json.dumps(tcr_trials.evidence(), indent=2))
    elif args.table == "benchmark":
        frame = tcr_trials.evaluate_against_gold_standard()
        print(frame.to_json(orient="records", indent=2))
    elif args.table == "benchmark_summary":
        print(json.dumps(tcr_trials.gold_standard_summary(), indent=2))
    else:
        frame = getattr(tcr_trials, args.table)()
        print(frame.to_json(orient="records", indent=2))
    return 0


def _cmd_cta(args: argparse.Namespace) -> int:
    print(
        "Legacy CTA snapshot; use 'oncoref antigens list --set CTA' for the current catalog.",
        file=sys.stderr,
    )
    if args.unfiltered:
        genes = cta.cta_unfiltered_gene_ids() if args.ids else cta.cta_unfiltered_gene_names()
    else:
        genes = cta.cta_gene_ids() if args.ids else cta.cta_gene_names()
    if args.count:
        print(len(genes))
    else:
        for g in sorted(genes):
            print(g)
    return 0


def _cmd_proteoforms(args: argparse.Namespace) -> int:
    if args.gene:
        symbol = proteoforms.proteoform_for_gene(args.gene)
        if symbol is None:
            print(f"{args.gene} is not in any proteoform group", file=sys.stderr)
            return 1
        members = proteoforms.proteoform_members_for_gene(args.gene) or ()
        print(f"{symbol}\t{', '.join(members)}")
        return 0
    symbol_map = proteoforms.proteoform_symbol_map()
    if args.count:
        print(len(symbol_map))
        return 0
    for label in sorted(symbol_map):
        print(f"{label}\t{', '.join(symbol_map[label])}")
    return 0


#: Plot families that write several figures into a directory rather than one PNG.
_DIRECTORY_PLOTS = frozenset(
    {"patient-coverage", "covering-set", "cta-curation", "expression-provenance"}
)

# Plot commands whose marks correspond one-to-one with registry cancer codes.
# Other commands use burden categories, genes, or multi-figure renderers and cannot
# truthfully promise that ``--highlight`` identifies a mark.
_HIGHLIGHTABLE_PLOTS = frozenset(
    {
        "apd1-vs-tmb",
        "apd1-orr-bars",
        "cta-addressable-burden",
        "cta-burden-vs-response",
        "cta-specific-9mer-load",
        "apd1-response-signature",
    }
)


def _cmd_plot(args: argparse.Namespace) -> int:
    from . import figure_output, figure_style, plots

    if args.highlight:
        if args.which not in _HIGHLIGHTABLE_PLOTS:
            print(
                f"Error: --highlight is not supported for {args.which}; this plot "
                "does not have one mark per cancer code",
                file=sys.stderr,
            )
            return 1
        try:
            args.highlight = cancer_types.resolve_cancer_type(args.highlight)
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1

    # Before any figure exists: the preset governs rcParams, page geometry and how
    # many rows survive, so it has to be chosen ahead of the first subplot.
    figure_style.use(args.preset)
    figure_style.set_highlight(args.highlight)

    # An explicit --out is used verbatim; omitting it opens a fresh timestamped run
    # directory under figures/, so a figure can always be traced to the run that
    # made it and no run overwrites another.
    name = args.which if args.preset == "print" else f"{args.which}-{args.preset}"
    if args.highlight:
        name = f"{name}-{args.highlight}"
    args.out = str(
        figure_output.resolve(
            args.out,
            name=name,
            suffix="" if args.which in _DIRECTORY_PLOTS else ".png",
        )
    )

    fns = {
        "apd1-vs-tmb": plots.apd1_vs_tmb,
        "apd1-orr-bars": plots.apd1_orr_bars,
        "incidence-vs-mortality": plots.incidence_vs_mortality,
        "cta-expression-heatmap": plots.cta_expression_heatmap,
        "cta-addressable-burden": plots.cta_addressable_burden,
        "cta-patient-heatmap": plots.cta_patient_count_heatmap,
        "cta-coverage-curves": plots.cta_coverage_curves,
        "cta-coverage-stacked": plots.cta_coverage_stacked_bars,
        "cta-burden-vs-response": plots.cta_burden_vs_response,
        "cta-specific-9mer-load": plots.cta_specific_9mer_load,
        "cta-covering-set": plots.cta_covering_set,
        "burden-category-bars": plots.burden_category_bars,
        "apd1-response-signature": plots.apd1_response_signature_scatter,
    }
    # Only pass threshold_tpm when the user explicitly set it, so each plot keeps its own
    # default (50 TPM for the TPM-based plots, within-sample p95 for the patient heatmap)
    # rather than being forced onto a single CLI-wide value.
    tpm = {} if args.threshold_tpm is None else {"threshold_tpm": args.threshold_tpm}
    merged_codes: list[str] = []
    if args.codes:
        merged_codes.extend(c.strip() for c in args.codes.split(",") if c.strip())
    if getattr(args, "cohort", None):
        merged_codes.extend(c.strip() for c in args.cohort if c and c.strip())
    selected_codes = list(dict.fromkeys(merged_codes)) if merged_codes else None

    try:
        if args.which == "patient-coverage":
            from . import coverage

            explicit_threshold = (
                args.threshold if args.threshold is not None else args.threshold_tpm
            )
            pc_kwargs: dict = {
                "cohorts": selected_codes,
                "out_dir": args.out,
            }
            if args.threshold_mode is not None:
                pc_kwargs["threshold_mode"] = args.threshold_mode
                if explicit_threshold is not None:
                    pc_kwargs["threshold"] = explicit_threshold
            else:
                pc_kwargs["threshold"] = (
                    explicit_threshold if explicit_threshold is not None else 25
                )
            result = coverage.render_patient_coverage(args.gene_set, **pc_kwargs)
            if result["n_cohorts"] == 0:
                print(
                    f"Error: no cohorts with cached per-sample data and coverage for "
                    f"gene set {result['label']!r}",
                    file=sys.stderr,
                )
                return 1
            thr_desc = result.get(
                "threshold_label",
                f"> {pc_kwargs.get('threshold', 25):g} TPM",
            )
            print(f"{result['label']}: {result['n_cohorts']} cohorts ({thr_desc})")
            for kind, path in result["paths"].items():
                print(f"  {kind}: {path}")
            return 0
        if args.which == "covering-set":
            from . import coverage

            cs_threshold = (
                args.threshold
                if args.threshold is not None
                else (
                    args.threshold_tpm
                    if args.threshold_tpm is not None
                    else coverage.ACTIONABLE_TPM
                )
            )
            result = coverage.render_covering_set(
                args.gene_set,
                stat=args.stat or "q3",
                threshold=cs_threshold,
                codes=selected_codes,
                n_genes=args.n_genes if args.n_genes is not None else 25,
                metric=args.metric,
                out_dir=args.out,
            )
            print(
                f"{result['label']}: {result['n_coverable']}/{result['n_cancer_types']} "
                f"cancer types coverable (> {cs_threshold:g} clean TPM at {result['stat']})"
            )
            for kind, path in result["paths"].items():
                print(f"  {kind}: {path}")
            return 0
        if args.which == "cta-curation":
            from . import cta_curation_plots

            kwargs = {"tissue_scope": args.tissue_scope} if args.tissue_scope != "core" else {}
            result = cta_curation_plots.render(out_dir=args.out, **kwargs)
            print(f"CTA curation figures ({result['n_genes']} evidence rows):")
            for label, remaining, dropped in result.get("stages", ()):
                suffix = f"  (-{dropped})" if dropped else ""
                print(f"  stage {label:<18} {remaining:>5}{suffix}")
            for kind, path in result["paths"].items():
                print(f"  {kind}: {path}")
            for kind, path in result.get("data_paths", {}).items():
                print(f"  audit {kind}: {path}")
            return 0
        if args.which == "expression-provenance":
            from . import expression_provenance_plots

            result = expression_provenance_plots.render(out_dir=args.out, region=args.region)
            summary, burden = result["summary"], result["burden"]
            print(
                f"Expression provenance figures: {summary['n_samples']:,} samples, "
                f"{summary['n_codes_covered']}/{summary['n_codes_total']} cancer types, "
                f"{summary['n_sources']} source projects"
            )
            inc = 100 * burden["incidence_covered"] / burden["incidence_total"]
            mor = 100 * burden["mortality_covered"] / burden["mortality_total"]
            print(
                f"  {args.region.upper()} burden covered: "
                f"{inc:.0f}% of incidence, {mor:.0f}% of mortality"
            )
            for kind, path in result["paths"].items():
                print(f"  {kind}: {path}")
            return 0
        if args.which in ("incidence-vs-mortality", "burden-category-bars"):
            kwargs = {"region": args.region}
        elif args.which == "apd1-response-signature":
            kwargs = {"signature": args.signature}
        elif args.which == "cta-expression-heatmap":
            kwargs = {"stat": args.stat or "median"}
        elif args.which == "cta-covering-set":
            kwargs = dict(tpm)
            if args.stat is not None:
                kwargs["stat"] = args.stat
            if args.n_genes is not None:
                kwargs["n_genes"] = args.n_genes
        elif args.which == "cta-addressable-burden":
            kwargs = {"source": args.source, **tpm}
        elif args.which == "cta-patient-heatmap":
            kwargs = {**tpm}
        elif args.which in ("cta-burden-vs-response", "cta-specific-9mer-load"):
            kwargs = {"against": args.against, **tpm}
        elif args.which in ("cta-coverage-curves", "cta-coverage-stacked"):
            if not selected_codes:
                print(f"Error: {args.which} needs --codes", file=sys.stderr)
                return 1
            kwargs = {
                "cancer_types": selected_codes,
                **tpm,
            }
        else:
            kwargs = {}
        fns[args.which](save=args.out, **kwargs)
    except (ValueError, ModuleNotFoundError, NotImplementedError, FileNotFoundError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    print(f"Wrote {args.out}")
    return 0


def _cmd_expression_sources(args: argparse.Namespace) -> int:
    from .species import resolve_species

    selected = resolve_species(args.species)
    srcs = expression_registry.expression_sources(species=selected, context=args.context)
    if args.code:
        srcs = [s for s in srcs if args.code in s.cancer_codes]
    if args.type:
        srcs = [s for s in srcs if s.source_type == args.type]
    srcs = sorted(srcs, key=lambda s: (s.context, s.source_type, s.id))
    if args.json:
        print(
            json.dumps(
                [{**vars(s), "contract": s.expression_contract.as_dict()} for s in srcs], indent=2
            )
        )
        return 0
    _listing_title(selected, "expression sources")
    _listing_text(f"{len(srcs)} sources")
    if not srcs:
        _listing_text("No matching expression sources.")
        return 0
    for context, kind in sorted({(s.context, s.source_type) for s in srcs}):
        entries = [s for s in srcs if (s.context, s.source_type) == (context, kind)]
        print()
        _listing_text(f"{context.replace('_', ' ').title()} / {kind} ({len(entries)})")
        records = []
        for source in entries:
            details = []
            if source.cancer_codes:
                details.append(
                    f"Cancer codes ({len(source.cancer_codes)}): {', '.join(source.cancer_codes)}"
                )
            if source.acquisition_status:
                details.append(f"Status: {source.acquisition_status.replace('_', ' ')}")
            if (
                source.acquisition_status in {"metadata_only", "blocked_provenance"}
                and source.notes
            ):
                details.append(f"Read unavailable: {source.notes}")
            records.append(((source.id, source.unit or "unknown"), details))
        _listing_table(("Source", "Native unit"), records)
    print()
    _listing_text("Units are source-native; unknown units remain unresolved.")
    _listing_text("Use --json for complete source metadata and measurement contracts.")
    return 0


def _cmd_samples(args: argparse.Namespace) -> int:
    if args.counts or (not args.code and not args.cohort):
        counts = samples.sample_counts_by_cancer_code(included_only=not args.all)
        for code, n in counts.items():
            print(f"{code}\t{n}")
        print(f"\n{int(counts.sum())} samples across {len(counts)} cancer codes", file=sys.stderr)
        return 0
    if args.code:
        df = samples.samples_for_cancer_code(args.code, included_only=not args.all)
    else:
        df = samples.samples_for_cohort(args.cohort, included_only=not args.all)
    cols = [c for c in ("sample_id", "source_cohort", "sample_type", "raw_unit") if c in df.columns]
    if df.empty:
        print("No matching samples.", file=sys.stderr)
        return 1
    print(df[cols].to_string(index=False))
    return 0


def _cmd_burden(args: argparse.Namespace) -> int:
    try:
        if args.category is None:
            for cat, pct in sorted(incidence.cancer_burden(metric=args.metric).items()):
                print(f"{cat}\t{pct:g}")
            return 0
        value = incidence.cancer_burden(args.category, metric=args.metric)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    if value is None:
        print(f"No burden value for category {args.category!r}", file=sys.stderr)
        return 1
    print(f"{value:g}")
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="oncoref",
        description="Curated cancer reference data: ontology, TMB, incidence/mortality, expression.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("version", help="Print the installed oncoref version").set_defaults(
        func=_cmd_version
    )

    # --- reference lookups ---
    p_ct = sub.add_parser(
        "cancer-type", help="Resolve a cancer type/alias/name and print its registry info"
    )
    p_ct.add_argument("query", help="Cancer code, alias, or display name (e.g. PRAD, prostate)")
    p_ct.set_defaults(func=_cmd_cancer_type)

    p_tmb = sub.add_parser("tmb", help="Median TMB (mut/Mb) for a code, or the full map")
    p_tmb.add_argument("code", nargs="?", default=None, help="Cancer code/alias (omit for all)")
    p_tmb.add_argument("--no-inherit", action="store_true", help="Do not inherit an ancestor's TMB")
    p_tmb.set_defaults(func=_cmd_tmb)

    p_apd1 = sub.add_parser(
        "apd1", help="Anti-PD-1 monotherapy ORR (%%) for a code, or the full map"
    )
    p_apd1.add_argument("code", nargs="?", default=None, help="Cancer code/alias (omit for all)")
    p_apd1.add_argument(
        "--no-inherit", action="store_true", help="Do not inherit an ancestor's ORR"
    )
    p_apd1.set_defaults(func=_cmd_apd1)

    p_ici = sub.add_parser(
        "ici", help="ICI ORR (%%) for a code by regimen (anti-PD-1/PD-L1/combo), or the full map"
    )
    p_ici.add_argument("code", nargs="?", default=None, help="Cancer code/alias (omit for all)")
    p_ici.add_argument(
        "--regimen",
        default=None,
        choices=list(ici.REGIMEN_FALLBACK),
        help="Pin a regimen (default: PD-1 → PD-L1 → PD-1+CTLA-4 fallback)",
    )
    p_ici.add_argument(
        "--all-regimens", action="store_true", help="Show every regimen present for the code"
    )
    p_ici.add_argument("--no-inherit", action="store_true", help="Do not inherit an ancestor's ORR")
    p_ici.set_defaults(func=_cmd_ici)

    from .antigen_cli import register as register_antigens

    register_antigens(sub)
    p_cta = sub.add_parser("cta", help="Legacy CTA snapshot (current catalog: antigens list)")
    p_cta.add_argument("--unfiltered", action="store_true", help="Full candidate universe")
    p_cta.add_argument("--ids", action="store_true", help="Ensembl gene IDs instead of symbols")
    p_cta.add_argument("--count", action="store_true", help="Print the count only")
    p_cta.set_defaults(func=_cmd_cta)

    p_trials = sub.add_parser(
        "tcr-trials",
        help="Reviewed TCR-T/ImmTAC cohorts, gold-standard toxicity benchmark, and events",
    )
    p_trials.add_argument(
        "--table",
        choices=[
            "all",
            "products",
            "cohorts",
            "events",
            "gold_standard",
            "benchmark",
            "benchmark_summary",
        ],
        default="products",
    )
    p_trials.set_defaults(func=_cmd_tcr_trials)

    p_proteoforms = sub.add_parser(
        "proteoforms",
        help="List identical-protein CGA groups whose TPM sums to proteoform level",
    )
    p_proteoforms.add_argument("--gene", help="Show the group for one gene (Ensembl ID or symbol)")
    p_proteoforms.add_argument("--count", action="store_true", help="Print the group count only")
    p_proteoforms.set_defaults(func=_cmd_proteoforms)

    p_plot = sub.add_parser("plot", help="Render a cancer-type reference plot to a PNG")
    p_plot.add_argument(
        "which",
        choices=[
            "apd1-vs-tmb",
            "apd1-orr-bars",
            "incidence-vs-mortality",
            "cta-expression-heatmap",
            "cta-addressable-burden",
            "cta-patient-heatmap",
            "cta-coverage-curves",
            "cta-coverage-stacked",
            "cta-burden-vs-response",
            "cta-specific-9mer-load",
            "cta-covering-set",
            "covering-set",
            "burden-category-bars",
            "apd1-response-signature",
            "patient-coverage",
            "cta-curation",
            "expression-provenance",
        ],
        help="Which plot to render",
    )
    p_plot.add_argument(
        "--out",
        default=None,
        help=(
            "Output PNG path, or output directory for the multi-figure families "
            "(patient-coverage, cta-curation, expression-provenance). Omit to write "
            "into a fresh figures/run_<YYYYMMDD-HHMMSS>/ directory."
        ),
    )
    p_plot.add_argument(
        "--highlight",
        default=None,
        help=(
            "Pick one cancer code out of the landscape (e.g. BRCA_Basal). Nothing is "
            "filtered or re-ranked: the same points stay in the same places, the rest "
            "fade to context."
        ),
    )
    p_plot.add_argument(
        "--preset",
        default="print",
        choices=["print", "slide"],
        help=(
            "Figure preset. 'print' is the dense manuscript figure; 'slide' uses large "
            "type, a 16:9 page, and keeps only the top rows so it reads from the back "
            "of a room."
        ),
    )
    p_plot.add_argument(
        "--signature",
        default="t_cell_inflamed",
        help="Response signature for apd1-response-signature (e.g. t_cell_inflamed, tgfb_exclusion)",
    )
    p_plot.add_argument(
        "--region",
        default="us",
        choices=["us", "world"],
        help="Region for incidence-vs-mortality and expression-provenance",
    )
    p_plot.add_argument(
        "--stat",
        default=None,
        choices=["q1", "median", "q3"],
        help=(
            "Statistic for cta-expression-heatmap and cta-covering-set. Omit to use "
            "each plot's own default (median for the heatmap, q3 for the covering set)."
        ),
    )
    p_plot.add_argument(
        "--tissue-scope",
        choices=["core", "extended", "both"],
        default="core",
        help="CTA curation panel(s) to render; core remains the default",
    )
    p_plot.add_argument(
        "--source",
        default="within_sample",
        choices=["within_sample", "per_sample"],
        help="Prevalence basis for cta-addressable-burden (per_sample needs matrices)",
    )
    p_plot.add_argument(
        "--against",
        default="apd1",
        choices=["apd1", "tmb"],
        help="Response metric for cta-burden-vs-response",
    )
    p_plot.add_argument(
        "--threshold-tpm",
        type=float,
        default=None,
        help=(
            "Clean-TPM 'expressed' cut for the per-sample CTA plots. Omit to use each "
            "plot's default (50 TPM, or within-sample p95 for the patient heatmap)."
        ),
    )
    p_plot.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="Coverage cutoff for patient-coverage or covering-set (matches pirlygenes).",
    )
    p_plot.add_argument(
        "--threshold-mode",
        choices=["auto", "tpm", "percentile"],
        default=None,
        help=(
            "Threshold contract for patient-coverage (auto uses clean TPM only when "
            "every selected source is linearly TPM-comparable, otherwise percentile)."
        ),
    )
    p_plot.add_argument(
        "--n-genes",
        type=int,
        default=None,
        help="Panel size to plot for covering-set / cta-covering-set (default: 25; 0 = full).",
    )
    p_plot.add_argument(
        "--metric",
        default="us_incidence_pct",
        choices=list(incidence._BURDEN_METRICS),
        help="Burden share column for covering-set patient weighting (default: us_incidence_pct).",
    )
    p_plot.add_argument(
        "--gene-set",
        default="cta",
        help=(
            "Gene set for patient-coverage or covering-set: cta, ca, ra, ca|ra, mito, "
            "housekeeping, family:<name>, therapy:<agent_class>, lineage:<code>, or CSV/TXT path"
        ),
    )
    p_plot.add_argument(
        "--codes",
        default=None,
        help="Comma-separated cancer codes (for cta coverage plots, covering-set, or patient-coverage)",
    )
    p_plot.add_argument(
        "--cohort",
        action="append",
        default=None,
        metavar="CODE",
        help="Restrict to specific cancer-type code(s); repeatable (matches pirlygenes).",
    )
    p_plot.set_defaults(func=_cmd_plot)

    p_burden = sub.add_parser(
        "burden", help="Incidence/mortality share for a burden category, or the full map"
    )
    p_burden.add_argument(
        "category", nargs="?", default=None, help="Burden category (omit for all)"
    )
    p_burden.add_argument(
        "--metric",
        default="us_incidence_pct",
        choices=list(incidence._BURDEN_METRICS),
        help="Which share to report (default: us_incidence_pct)",
    )
    p_burden.set_defaults(func=_cmd_burden)

    p_data = sub.add_parser(
        "data",
        help="Unified data management over every managed dataset (expression bundle + HPA)",
    )
    p_data.add_argument(
        "action",
        choices=[
            "list",
            "status",
            "info",
            "fetch",
            "path",
            "verify",
            "dirs",
            "prune",
            # Backward-compatible aliases:
            "dir",
            "contract",
            "metadata",
            "release-manifest",
            "versions",
            "provenance",
        ],
        help=(
            "list (catalog), status (local state), info (versions/provenance/contract), "
            "fetch (download), path (resolve local path), verify (integrity check), "
            "dirs (cache roots), prune (safe retention cleanup)"
        ),
    )
    p_data.add_argument(
        "name",
        nargs="?",
        default=None,
        help="Dataset name; or a fetch target: all / hpa / reference / bundle / source / "
        "per-sample:<CODE> (omit = all)",
    )
    p_data.add_argument(
        "--progress",
        choices=("auto", "always", "never"),
        default=None,
        help="Download bars: auto (terminal/notebook), always, or never; "
        "defaults to ONCOREF_DOWNLOAD_PROGRESS or auto",
    )
    p_data.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress download status and bars; retain results, errors, and retry warnings",
    )
    p_data.add_argument(
        "--version",
        dest="data_version",
        help="Concrete source release; default is newest registered (same as latest)",
    )
    p_data.add_argument("--force", action="store_true", help="Re-download even if cached")
    p_data.add_argument(
        "--yes",
        action="store_true",
        help="Actually delete stale bundle caches for 'data prune' (default: dry run)",
    )
    p_data.add_argument(
        "--include-current",
        action="store_true",
        help="For 'data prune', also delete current pinned cache entries",
    )
    p_data.add_argument(
        "--keep-version",
        dest="keep_versions",
        action="append",
        default=[],
        help="Additional version/digest pin to retain during 'data prune' (repeatable)",
    )
    p_data.add_argument(
        "--include-legacy-pirlygenes",
        action="store_true",
        help="For 'data prune', also inspect/prune obsolete ~/.cache/pirlygenes/bundled_data entries",
    )
    p_data.add_argument(
        "--species",
        default="human",
        help="Species name or NCBI taxonomy ID (default: human; dog/canine are aliases)",
    )
    p_data.add_argument(
        "--json",
        action="store_true",
        help="Emit machine-readable JSON (supported across all 'oncoref data' subcommands)",
    )
    p_data.add_argument("--group", default=None, help="Conceptual data group for list")
    p_data.add_argument(
        "--offline",
        action="store_true",
        help="Forbid network access (also honored via OPENVAX_OFFLINE=1)",
    )
    p_data.add_argument(
        "--remote",
        action="store_true",
        help="Include remote release manifest when running 'oncoref data info bundle'",
    )
    p_data.add_argument(
        "--asset-overrides", help="JSON mapping taxon:asset IDs to existing local files"
    )
    p_data.add_argument(
        "--cache-root", help="Override the unified OncoRef cache root for this command"
    )
    p_data.add_argument(
        "--max-bytes",
        type=int,
        default=100 * 1024 * 1024,
        help="Maximum total species-asset transfer bytes",
    )
    p_data.set_defaults(func=_cmd_data)

    p_species = sub.add_parser("species", help="List species or inspect a species' data groups")
    p_species.add_argument("name", nargs="?", default=None)
    p_species.set_defaults(func=_cmd_species)

    p_exprsrc = sub.add_parser(
        "expression-sources",
        help="List the cohort expression sources (Treehouse, GDC, GEO, recount3, …)",
    )
    p_exprsrc.add_argument("--code", default=None, help="Only sources feeding this cancer code")
    p_exprsrc.add_argument("--type", default=None, help="Only this source_type (e.g. gdc, geo)")
    p_exprsrc.add_argument(
        "--species",
        default="human",
        help="Species name or NCBI taxonomy ID (default: human; dog/canine are aliases)",
    )
    p_exprsrc.add_argument(
        "--json", action="store_true", help="Print complete source records and contracts as JSON"
    )
    p_exprsrc.add_argument("--context", choices=sorted(expression_registry.EXPRESSION_CONTEXTS))
    p_exprsrc.set_defaults(func=_cmd_expression_sources)

    p_samples = sub.add_parser(
        "samples", help="Per-sample curation manifest (counts, or rows for a code/cohort)"
    )
    p_samples.add_argument("--code", default=None, help="Samples assigned to this cancer code")
    p_samples.add_argument("--cohort", default=None, help="Samples from this source cohort")
    p_samples.add_argument("--counts", action="store_true", help="Per-cancer-code sample counts")
    p_samples.add_argument(
        "--all", action="store_true", help="Include excluded samples (default: included only)"
    )
    p_samples.set_defaults(func=_cmd_samples)

    return parser


def main(argv=None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

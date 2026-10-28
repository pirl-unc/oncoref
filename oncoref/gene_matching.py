"""Central, ambiguity-preserving matching for human reference data tables.

Registry matches are distinct from annotation-verified ``gene_identity`` records.
Stable IDs win over names, but a contradictory known name is quarantined. A name
cannot turn an unknown stable ID into a different locus. Original identifiers,
all alternatives, reference versions and reference hashes remain in the audit.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Sequence
from functools import lru_cache
from typing import Literal

import numpy as np
import pandas as pd

from .gene_ids import unversioned
from .load_dataset import _register_derived_cache, get_data


def _text(value: object) -> str:
    return "" if value is None or pd.isna(value) else str(value).strip()


def _json(values: Iterable[str]) -> str:
    return json.dumps(sorted(set(values)), separators=(",", ":"))


@lru_cache(maxsize=1)
def _catalog():
    canonical = get_data("canonical-gene-space").fillna("")
    aliases = get_data("ensembl-id-aliases").fillna("")
    synonyms = get_data("ncbi-symbol-synonyms").fillna("")
    ids = set(canonical.ensembl_gene_id)
    targets, methods = {}, {}
    for r in aliases.itertuples():
        a, b = unversioned(r.alt_haplotype_id), unversioned(r.primary_contig_id)
        targets.setdefault(a, set()).add(b)
        methods.setdefault(a, set()).add(r.source)

    def terminal(gid, visiting=frozenset()):
        if gid in visiting:
            raise ValueError(f"Cycle in registered Ensembl aliases: {gid}")
        if gid not in targets:
            return {gid} if gid in ids else set()
        result = {gid} if gid in ids else set()
        for target in targets[gid]:
            result |= terminal(target, visiting | {gid})
        return result

    resolved = {alias: terminal(alias) for alias in targets}
    by_name = {}
    for r in canonical.itertuples():
        if r.symbol:
            by_name.setdefault(r.symbol.upper(), set()).add(r.ensembl_gene_id)
    names = {gid: set() for gid in ids}
    for name, gene_ids in by_name.items():
        for gid in gene_ids:
            names[gid].add(name)
    synonym_index = {}
    for r in synonyms.itertuples():
        candidates = by_name.get(r.official_symbol.upper(), set())
        if r.alias and candidates:
            synonym_index.setdefault(r.alias.upper(), set()).update(candidates)
            for gid in candidates:
                names[gid].add(r.alias)
    history, alternate = {gid: set() for gid in ids}, {gid: set() for gid in ids}
    for alias, candidates in resolved.items():
        if len(candidates) != 1 or alias in ids:
            continue
        target = next(iter(candidates))
        dest = alternate if methods[alias] == {"ensembl_alt_contig"} else history
        dest[target].add(alias)
    payload = {
        "canonical": sorted(map(tuple, canonical.astype(str).to_numpy())),
        "aliases": sorted(map(tuple, aliases.astype(str).to_numpy())),
        "synonyms": sorted(map(tuple, synonyms.astype(str).to_numpy())),
    }
    signature = hashlib.sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()
    records = canonical.rename(
        columns={"ensembl_gene_id": "canonical_gene_id", "symbol": "canonical_symbol"}
    )
    records["historical_ensembl_gene_ids"] = records.canonical_gene_id.map(
        lambda gid: _json(history[gid])
    )
    records["alternate_ensembl_gene_ids"] = records.canonical_gene_id.map(
        lambda gid: _json(alternate[gid])
    )
    records["gene_name_synonyms"] = records.canonical_gene_id.map(lambda gid: _json(names[gid]))
    records["identity_catalog_sha256"] = signature
    records["taxonomy_id"] = 9606
    return records, ids, resolved, methods, by_name, synonym_index


_register_derived_cache(_catalog.cache_clear)


def gene_identity_catalog() -> pd.DataFrame:
    """One row per canonical ENSG, with historical IDs, synonyms and source hash.

    Empty lists mean no mapping is registered, not that no historical name exists.
    Alternate-contig IDs remain separate from historical replacement IDs.
    """
    return _catalog()[0].copy()


def match_reference_genes(
    frame: pd.DataFrame,
    *,
    id_column: str | None,
    symbol_column: str | None = None,
    source: str,
    source_version: str,
    taxonomy_id: int = 9606,
) -> pd.DataFrame:
    """Return one auditable match per input row; never select the first ambiguity.

    An explicit stable ID is never replaced by a symbol-only guess. Source labels
    that are unknown synonyms are retained as a name mismatch alongside a valid
    stable-ID match. Names resolving to another known locus quarantine the row.
    Source version is recorded and does not imply verification against its GTF.
    """
    if taxonomy_id != 9606:
        raise ValueError("Human registry cannot match another taxon")
    if not source or not source_version:
        raise ValueError("Explicit source and source version are required")
    catalog, ids, aliases, methods, official, synonyms = _catalog()
    pairs = pd.DataFrame(
        {
            "source_gene_id": frame[id_column].map(_text).to_numpy()
            if id_column
            else [""] * len(frame),
            "source_gene_name": frame[symbol_column].map(_text).to_numpy()
            if symbol_column
            else [""] * len(frame),
        }
    ).astype({"source_gene_id": "str", "source_gene_name": "str"})
    decisions = []
    for row in pairs.drop_duplicates().itertuples(index=False):
        raw, name = row
        gid = unversioned(raw).upper() if raw else ""
        # NaN is also a historical SCN11A synonym. In tabular reference input,
        # missing-value sentinels cannot supply or contradict gene identity.
        missing_name = name.casefold() in {"", "nan", "na", "n/a", "none", "null", "-"}
        name_candidates = (
            set() if missing_name else official.get(name.upper(), synonyms.get(name.upper(), set()))
        )
        candidates, status, method, chosen = set(), "unmapped", "", ""
        if raw:
            candidates = aliases.get(gid, {gid} if gid in ids else set())
            if len(candidates) > 1:
                status = "ambiguous_id"
            elif len(candidates) == 1:
                candidate = next(iter(candidates))
                if name_candidates and candidate not in name_candidates:
                    status = "id_name_conflict"
                    candidates = candidates | name_candidates
                else:
                    chosen = candidate
                    status = "matched"
                    method = (
                        "canonical_id"
                        if gid == chosen
                        else "+".join(sorted(methods.get(gid, {"registered_alias"})))
                    )
            else:
                status, candidates = "unmapped_id", name_candidates
        elif len(name_candidates) == 1:
            chosen = next(iter(name_candidates))
            status, method, candidates = (
                "matched",
                "unique_symbol" if name.upper() in official else "unique_synonym",
                name_candidates,
            )
        elif len(name_candidates) > 1:
            status, candidates = "ambiguous_name", name_candidates
        decisions.append(
            {
                "source_gene_id": raw,
                "source_gene_name": name,
                "canonical_gene_id": chosen,
                "identity_status": status,
                "identity_method": method,
                "candidate_gene_ids": _json(candidates),
                "source_name_unrecognized": bool(chosen and name and not name_candidates),
            }
        )
    schema = [
        "source_gene_id",
        "source_gene_name",
        "canonical_gene_id",
        "identity_status",
        "identity_method",
        "candidate_gene_ids",
        "source_name_unrecognized",
    ]
    result = pairs.merge(
        pd.DataFrame(decisions, columns=schema),
        on=["source_gene_id", "source_gene_name"],
        how="left",
        validate="many_to_one",
    )
    result = result.merge(catalog, on="canonical_gene_id", how="left", validate="many_to_one")
    result["source"], result["source_version"], result["taxonomy_id"] = (
        source,
        source_version,
        taxonomy_id,
    )
    result["identity_catalog_sha256"] = catalog.identity_catalog_sha256.iloc[0]
    result.index = frame.index
    return result


def attach_gene_identities(frame: pd.DataFrame, *, gene_column: str = "gene_id") -> pd.DataFrame:
    """Attach the complete identity bundle to a stage result without changing rows."""
    catalog = gene_identity_catalog().rename(columns={"canonical_gene_id": gene_column})
    columns = [
        gene_column,
        "canonical_symbol",
        "historical_ensembl_gene_ids",
        "alternate_ensembl_gene_ids",
        "gene_name_synonyms",
        "identity_catalog_sha256",
    ]
    d = frame.drop(columns=[c for c in columns[1:] if c in frame]).merge(
        catalog[columns], on=gene_column, how="left", validate="many_to_one"
    )
    if d.identity_catalog_sha256.isna().any():
        raise ValueError("Stage contains an identity outside the canonical registry")
    d.index = frame.index
    return d.copy()


def harmonize_reference_expression(
    frame: pd.DataFrame,
    *,
    sample_cols: Sequence[str],
    source: str,
    source_version: str,
    id_column: str = "Ensembl_Gene_ID",
    symbol_column: str | None = "Symbol",
    duplicate_policy: Literal["error", "equal_alias_rows", "quarantine"] = "error",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Resolve loci before sums, removing only explicitly allowed equal alias rows.

    Distinct protein-identical loci remain distinct here. Multiple rows resolving
    to one locus are never silently added. Unmatched rows stay in the returned
    audit; their measurements are not assigned to a guessed gene. ``quarantine``
    deduplicates equal alias rows but excludes every conflicting observation of
    a locus, retaining those source rows and the reason in the identity audit.
    """
    if duplicate_policy not in {"error", "equal_alias_rows", "quarantine"}:
        raise ValueError("Unknown duplicate policy")
    sample_cols = list(sample_cols)
    if not sample_cols or len(set(sample_cols)) != len(sample_cols):
        raise ValueError("Distinct explicit sample/context columns are required")
    matches = match_reference_genes(
        frame,
        id_column=id_column,
        symbol_column=symbol_column,
        source=source,
        source_version=source_version,
    )
    d = frame[list(sample_cols)].reset_index(drop=True).copy()
    values = d.to_numpy(dtype=float)
    if np.isinf(values).any() or (values < 0).any():
        raise ValueError("Nonnegative linear abundance or missing values required")
    d = pd.DataFrame(values, columns=sample_cols)
    d["Ensembl_Gene_ID"] = matches.canonical_gene_id.to_numpy()
    d["Symbol"] = matches.canonical_symbol.to_numpy()
    d = d[matches.identity_status.eq("matched").to_numpy()]
    keep = pd.Series(True, index=d.index)
    for _, group in d[d.Ensembl_Gene_ID.duplicated(keep=False)].groupby("Ensembl_Gene_ID"):
        a = group[list(sample_cols)].to_numpy(dtype=float)
        equal = np.all((a == a[0]) | (np.isnan(a) & np.isnan(a[0])))
        source_ids = matches.iloc[group.index].source_gene_id.map(unversioned)
        if duplicate_policy == "quarantine" and (not equal or not source_ids.is_unique):
            keep.loc[group.index] = False
            matches.iloc[group.index, matches.columns.get_loc("identity_status")] = (
                "conflicting_duplicate_measurements"
            )
            continue
        if duplicate_policy == "error" or not equal or not source_ids.is_unique:
            raise ValueError(
                f"Conflicting/repeated measurements for {group.Ensembl_Gene_ID.iloc[0]}"
            )
        keep.loc[group.index[1:]] = False
        matches.iloc[group.index[1:], matches.columns.get_loc("identity_status")] = (
            "duplicate_alias_not_added"
        )
    return d.loc[keep, ["Ensembl_Gene_ID", "Symbol", *sample_cols]].reset_index(drop=True), matches


def collapse_reference_expression(
    frame: pd.DataFrame,
    *,
    sample_cols: Sequence[str],
    source: str,
    source_version: str,
    unit: str,
    id_column: str = "Ensembl_Gene_ID",
    symbol_column: str | None = "Symbol",
    scope: str = "genome",
    duplicate_policy: Literal["error", "equal_alias_rows", "quarantine"] = "error",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compatibility import; use :mod:`oncoref.protein_equivalence` for aggregation."""
    from .protein_equivalence import collapse_reference_expression as collapse

    return collapse(
        frame,
        sample_cols=sample_cols,
        source=source,
        source_version=source_version,
        unit=unit,
        id_column=id_column,
        symbol_column=symbol_column,
        scope=scope,
        duplicate_policy=duplicate_policy,
    )

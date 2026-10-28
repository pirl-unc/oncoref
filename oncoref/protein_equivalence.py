"""Same-context expression sums for identical representative-protein groups.

A gene ID identifies a locus. A protein equivalence group identifies loci whose
chosen reference proteins have identical sequences in the pinned registry. It
is not a complete catalogue of splice isoforms or post-translational proteoforms.
Display aliases (NY-ESO-1, XAGE1A/B) never replace member gene identities.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Sequence
from typing import Literal, NamedTuple

import numpy as np
import pandas as pd

from .gene_matching import gene_identity_catalog, harmonize_reference_expression
from .proteoforms import proteoform_group_map, proteoform_symbol

LinearUnit = Literal["TPM", "nTPM", "pTPM", "CP10K", "counts"]


def _json(values: Iterable[str]) -> str:
    return json.dumps(sorted(set(values)), separators=(",", ":"))


def group_identifier(members: Iterable[str]) -> str:
    """Deterministic membership ID; deliberately not a protein sequence hash."""
    return "PEG:" + hashlib.sha256(_json(members).encode()).hexdigest()[:20]


def equivalence_groups(gene_ids: Iterable[str], *, scope: str = "genome") -> pd.DataFrame:
    """One row per group intersecting the requested loci, including singletons.

    All expected members remain listed, including members outside the requested
    universe. Their absence must not be interpreted as measured zero expression.
    """
    ids = set(gene_ids)
    catalog = gene_identity_catalog().set_index("canonical_gene_id")
    unknown = ids - set(catalog.index)
    if unknown:
        raise ValueError(f"Unknown canonical gene IDs: {sorted(unknown)}")
    registry = proteoform_group_map(scope=scope)
    registry_hash = hashlib.sha256(json.dumps(registry, sort_keys=True).encode()).hexdigest()
    groups = {
        proteoform_symbol(k): tuple(sorted(v)) for k, v in registry.items() if ids.intersection(v)
    }
    grouped = {gid for values in groups.values() for gid in values}
    groups.update({gid: (gid,) for gid in sorted(ids - grouped)})
    fields = [
        "canonical_symbol",
        "historical_ensembl_gene_ids",
        "alternate_ensembl_gene_ids",
        "gene_name_synonyms",
    ]
    bundles = catalog[fields].to_dict("index")
    catalog_hash = catalog.identity_catalog_sha256.iloc[0]
    rows = []
    for label, members in sorted(groups.items()):
        identities = [
            {
                "canonical_gene_id": gid,
                **bundles.get(gid, dict.fromkeys(fields)),
                "identity_status": "matched"
                if gid in bundles
                else "registry_member_outside_canonical_catalog",
            }
            for gid in members
        ]
        rows.append(
            {
                "protein_equivalence_group_id": group_identifier(members),
                "protein_group_label": label
                if len(members) > 1
                else bundles[members[0]]["canonical_symbol"],
                "proteoform_key": label,
                "representative_gene_id": members[0],
                "member_gene_ids": _json(members),
                "member_count": len(members),
                "member_identities": json.dumps(identities, separators=(",", ":")),
                "group_registry_sha256": registry_hash,
                "identity_catalog_sha256": catalog_hash,
            }
        )
    return pd.DataFrame(rows)


class CollapsedExpression(NamedTuple):
    """Observed sums, per-context member coverage, and source identity audit."""

    expression: pd.DataFrame
    coverage: pd.DataFrame
    identity_audit: pd.DataFrame


def sum_context_observations(
    observations: pd.DataFrame,
    groups: pd.DataFrame,
    *,
    context_columns: Sequence[str],
    value_columns: Sequence[str],
    unit: LinearUnit,
    gene_column: str = "gene_id",
) -> pd.DataFrame:
    """Sum distinct canonical loci within explicitly matched measurement contexts.

    Inputs must already pass identity matching. Each locus has at most one row
    per context. Unit and source/cohort/sample boundaries belong in the context
    keys; rank slots, gene maxima and categorical IHC are not valid inputs.
    ``complete_members`` is separate from a positive lower-bound measurement.
    """
    if unit not in {"TPM", "nTPM", "pTPM", "CP10K", "counts"}:
        raise ValueError("An explicit supported linear abundance unit is required")
    if not context_columns or not value_columns:
        raise ValueError("Explicit measurement contexts and values are required")
    if observations.duplicated([gene_column, *context_columns]).any():
        raise ValueError("Repeated locus/context observations cannot be added")
    values = observations[list(value_columns)].to_numpy(dtype=float)
    if np.isinf(values).any() or (values < 0).any():
        raise ValueError("Nonnegative linear abundance or missing values required")
    members = {
        gid: r.protein_equivalence_group_id
        for r in groups.itertuples()
        for gid in json.loads(r.member_gene_ids)
    }
    d = observations[observations[gene_column].isin(members)].copy()
    d["protein_equivalence_group_id"] = d[gene_column].map(members)
    keys = ["protein_equivalence_group_id", *context_columns]
    grouped = d.groupby(keys, dropna=False, sort=True)
    result = grouped[list(value_columns)].sum(min_count=1)
    for value in value_columns:
        result[f"measured_members:{value}"] = grouped[value].count()
    result = result.reset_index().merge(
        groups, on="protein_equivalence_group_id", validate="many_to_one"
    )
    result["complete_members"] = (
        result[[f"measured_members:{v}" for v in value_columns]]
        .eq(result.member_count, axis=0)
        .all(axis=1)
    )
    result["unit"] = unit
    return result


def unpaired_rank_bounds(
    values: Sequence[Sequence[float]], *, rank: int = 3
) -> tuple[float, float]:
    """Bounds on a ranked sum when genes lack paired specimen identities.

    Lower bound is the largest member's ranked value. Upper bound follows the
    union bound: allocate at most ``rank - 1`` exceedances among member ranks.
    It assumes a common specimen population; unequal/missing lengths give an
    unbounded upper limit. No sorted values are paired as invented donors.
    """
    arrays = [np.sort(np.asarray(v, dtype=float))[::-1] for v in values]
    if rank < 1 or not arrays or any(np.isinf(a).any() or (a < 0).any() for a in arrays):
        raise ValueError("Positive rank and finite nonnegative observations required")
    lower = max(
        (a[rank - 1] for a in arrays if len(a) >= rank and np.isfinite(a[rank - 1])), default=np.nan
    )
    if len({len(a) for a in arrays}) != 1 or any(
        len(a) < rank or np.isnan(a).any() for a in arrays
    ):
        return float(lower), float("inf")
    # For large groups dynamic programming avoids an exponential product.
    costs = np.full(rank, np.inf)
    costs[0] = 0.0
    for a in arrays:
        costs = np.array([min(costs[k - j] + a[j] for j in range(k + 1)) for k in range(rank)])
    return float(lower), float(costs.min())


def collapse_reference_expression(
    frame: pd.DataFrame,
    *,
    sample_cols: Sequence[str],
    source: str,
    source_version: str,
    unit: LinearUnit,
    id_column: str = "Ensembl_Gene_ID",
    symbol_column: str | None = "Symbol",
    scope: str = "genome",
    duplicate_policy: Literal["error", "equal_alias_rows", "quarantine"] = "error",
) -> CollapsedExpression:
    """Identity-audited, same-sample sums of distinct identical-protein loci.

    Returns (observed_sums, coverage, identity_audit). Values are nonnegative
    linear abundance in ONE explicitly named native unit. Columns must identify
    the SAME samples or measured contexts for every row; ranked donor slots,
    per-gene maxima, percentiles and IHC categories are not additive inputs.

    A partial sum is a lower bound. Missing loci are not measured zeros: the
    separate group-by-sample coverage frame has ``complete=False``. Partial sums
    can establish positive expression but must not establish normal absence.
    Groups identify identical chosen representative proteins, not every isoform.
    """
    from .proteoforms import (
        collapse_to_proteoforms,
        proteoform_group_map,
        proteoform_symbol,
    )

    if unit not in {"TPM", "nTPM", "pTPM", "CP10K", "counts"}:
        raise ValueError("An explicit supported linear abundance unit is required")
    sample_cols = list(sample_cols)
    canonical, audit = harmonize_reference_expression(
        frame,
        sample_cols=sample_cols,
        source=source,
        source_version=source_version,
        id_column=id_column,
        symbol_column=symbol_column,
        duplicate_policy=duplicate_policy,
    )
    if canonical.empty:
        raise ValueError("No uniquely matched expression rows")
    collapsed = collapse_to_proteoforms(canonical, scope=scope, sample_cols=sample_cols)
    groups = {
        proteoform_symbol(k): tuple(sorted(v)) for k, v in proteoform_group_map(scope=scope).items()
    }
    genes = canonical.set_index("Ensembl_Gene_ID")
    identities = gene_identity_catalog().set_index("canonical_gene_id")
    catalog = equivalence_groups(canonical.Ensembl_Gene_ID, scope=scope)
    bundles = catalog.set_index("proteoform_key").member_identities
    metadata, coverage = [], []
    for row in collapsed.itertuples():
        key = row.proteoform_key
        expected = groups.get(key, (row.Ensembl_Gene_ID,))
        observed = sorted(set(expected) & set(genes.index))
        missing = sorted(set(expected) - set(observed))
        n_measured = genes.loc[observed, sample_cols].notna().sum(axis=0)
        metadata.append(
            {
                "proteoform_key": key,
                "member_gene_ids": _json(expected),
                "observed_member_gene_ids": _json(observed),
                "missing_member_gene_ids": _json(missing),
                "member_identities": bundles[key],
            }
        )
        coverage.extend(
            {
                "proteoform_key": key,
                "sample": sample,
                "expected_members": len(expected),
                "measured_members": int(count),
                "complete": count == len(expected),
            }
            for sample, count in n_measured.items()
        )
    attrs = collapsed.attrs.copy()
    collapsed = collapsed.merge(pd.DataFrame(metadata), on="proteoform_key", validate="one_to_one")
    collapsed["source"], collapsed["source_version"], collapsed["unit"] = (
        source,
        source_version,
        unit,
    )
    collapsed["identity_catalog_sha256"] = identities.identity_catalog_sha256.iloc[0]
    collapsed.attrs.update(attrs)
    extra = [
        "protein_equivalence_group_id",
        "protein_group_label",
        "representative_gene_id",
        "member_count",
        "group_registry_sha256",
    ]
    collapsed = collapsed.merge(
        catalog[["proteoform_key", *extra]], on="proteoform_key", validate="one_to_one"
    )
    covered = pd.DataFrame(coverage).merge(
        catalog[["proteoform_key", "protein_equivalence_group_id"]],
        on="proteoform_key",
        validate="many_to_one",
    )
    return CollapsedExpression(collapsed, covered, audit)

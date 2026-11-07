"""Cited CTA literature reviews and distinct source roles.

This registry annotates the frozen 1.31.0 selection. A paper's placement in the
core or supplementary reading set never determines whether its findings apply.
Normal-expression and identity review findings are not yet classification edits.
"""

from __future__ import annotations

import gzip
import json
from importlib.resources import files

import pandas as pd


def registry() -> dict:
    """Return an independent copy of the packaged, versioned evidence registry."""
    payload = files("oncoref").joinpath("data/cta-literature.json.gz").read_bytes()
    return json.loads(gzip.decompress(payload))


def sources(*, collection: str | None = None, role: str | None = None) -> pd.DataFrame:
    """Source citations, coverage, overlap and retention reasons.

    `collection` is core or supplementary. `role` is one of registry()['roles'].
    All roles remain available regardless of whether a source adds unique names.
    """
    data = registry()
    if collection not in {None, "core", "supplementary"}:
        raise ValueError("collection must be core, supplementary or None")
    if role is not None and role not in data["roles"]:
        raise ValueError(f"Unknown literature role: {role}")
    result = pd.DataFrame(data["sources"])
    if collection is not None:
        result = result[result.collection.eq(collection)]
    if role is not None:
        result = result[result.roles.map(lambda roles: role in roles)]
    return result.reset_index(drop=True)


def reviews(gene: str | None = None) -> pd.DataFrame:
    """Individual reviews for the 41 unmatched and 55 broad-table-only proteins.

    An optional exact protein label or Ensembl gene ID selects its review. Source
    scope and unresolved assay identity are retained alongside every citation.
    """
    result = pd.DataFrame(registry()["reviews"])
    if gene is not None:
        result = result[result.protein.eq(gene) | result.gene_ids.map(lambda ids: gene in ids)]
    return result.reset_index(drop=True)


def memberships() -> pd.DataFrame:
    """All paper-to-protein links, including unresolved and noncoverage links."""
    return pd.DataFrame(registry()["memberships"])


def cta_coverage() -> pd.DataFrame:
    """Paper coverage for all 602 released CTAs; empty lists are coverage gaps."""
    return pd.DataFrame(registry()["cta_coverage"])


def _protein_rows(key: str, gene: str | None) -> pd.DataFrame:
    result = pd.DataFrame(registry()[key])
    if gene is not None:
        result = result[result.protein.eq(gene) | result.gene_ids.map(lambda ids: gene in ids)]
    return result.reset_index(drop=True)


def gap_assessments(gene: str | None = None) -> pd.DataFrame:
    """All 33 paper gaps: biology, released measurements and review recommendations."""
    return _protein_rows("gap_assessments", gene)


def concerns(gene: str | None = None) -> pd.DataFrame:
    """The 29 flagged reviews, separated by evidence type and next action."""
    return _protein_rows("concerns", gene)


def localization(gene: str | None = None) -> pd.DataFrame:
    """Gene-mapped locations for all 602 CTAs, with shared and partial scope.

    Facets are annotation summaries, not validated tumor-surface epitopes or
    selected-sequence localization. Native evidence, isoform notes and primary
    paper additions are retained in each row.
    """
    return _protein_rows("localization", gene)

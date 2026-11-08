"""Cited CTA literature reviews and distinct source roles.

This registry annotates the frozen 1.31.0 selection. A paper's placement in the
core or supplementary reading set never determines whether its findings apply.
This frozen registry does not itself edit classifications. Current transcript
decisions and reviewed coding holds are exposed by ``transcript_selection``.
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


def surface_domains(gene: str | None = None) -> pd.DataFrame:
    """37 membrane CTAs: inclusive extracellular coordinates and sequence attribution.

    UniProt canonical dimensions and selected-sequence projections are separate.
    Missing boundaries are unknown, not zero. Predicted topology is not measured
    tumor-surface accessibility.
    """
    return _protein_rows("surface_domains", gene)


def surface_cancer_rna(gene: str | None = None, *, cohort: str | None = None) -> pd.DataFrame:
    """HPA 25.1 native-pTPM sample counts for the 37 membrane CTAs.

    TCGA and validation cohorts remain separate. Numerators use >= thresholds;
    denominators count measured samples, not the population prevalence of cancer.
    """
    if cohort not in {None, "TCGA", "validation"}:
        raise ValueError("cohort must be TCGA, validation or None")
    result = _protein_rows("surface_cancer_rna", gene)
    if cohort is not None:
        result = result[result.cohort.eq(cohort)]
    return result.reset_index(drop=True)


def surface_cancer_ihc(gene: str | None = None) -> pd.DataFrame:
    """HPA scored-patient staining counts; retain unmeasured records explicitly.

    These counts do not establish antibody specificity, subcellular localization,
    the fraction of malignant cells positive, or selected-proteoform expression.
    """
    return _protein_rows("surface_cancer_ihc", gene)

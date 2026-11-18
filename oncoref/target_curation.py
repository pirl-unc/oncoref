# Licensed under the Apache License, Version 2.0.
"""Read the catalog's central curation records: literature, family rules, tissues.

This module owns no data. It is a view over the stores the catalog itself is built
from, so a filter can never act on a record the catalog has not seen:

- ``normal-expression-literature.json``: one reviewed entry per protein, the rubric
  that decides Unsafe, and the family rules.
- the tissue-tier registry and its updates, through
  :mod:`oncoref.tissue_tiers`.
- ``nose-skin-reference.csv``: normal nose and skin expression of each candidate.
- ``readthrough-loci.csv``: loci the catalog's sequence check moves to Other.

A literature entry either excludes (Unsafe in the catalog) or is recorded with its
tissue and tier. :func:`table` reports the recorded tier that
:mod:`oncoref.target_filters` may act on.
"""

from __future__ import annotations

import json
from functools import cache
from importlib.resources import files

import pandas as pd

from . import tissue_tiers

TIERS = tissue_tiers.TIERS  # most protective first
_REFERENCE_SOURCES = (
    "https://doi.org/10.1186/s12864-016-2960-3 | https://chanzuckerberg.github.io/cellxgene-census/"
)
_DATA = files("oncoref") / "data"


@cache
def _literature() -> dict:
    return json.loads((_DATA / "normal-expression-literature.json").read_text())


def literature() -> dict:
    """The literature review's scope, rubric, method and limits (entries left out)."""
    return {k: v for k, v in _literature().items() if k != "entries"}


def most_protective(tiers) -> str:
    """The most protective tier in ``tiers``, or "" when none is given."""
    given = set(tiers)
    return next((t for t in TIERS if t in given), "")


@cache
def _entries() -> pd.DataFrame:
    columns = [
        "entity_id",
        "protein",
        "review",
        "verdict",
        "effect",
        "rule",
        "decision_tissue",
        "tissue_tier",
        "recorded_tissue",
        "recorded_tier",
        "decision_reason",
    ]
    frame = pd.DataFrame(_literature()["entries"])
    frame["sources"] = frame.evidence.map(
        lambda items: " | ".join(dict.fromkeys(x["url"] for x in items if x.get("url")))
    )
    return frame[[*columns, "sources"]].fillna("")


def entries(protein: str | None = None) -> pd.DataFrame:
    """One row per reviewed protein: verdict, effect, deciding and recorded tissue."""
    d = _entries()
    return (d[d.protein == protein] if protein else d).reset_index(drop=True)


def findings(protein: str | None = None) -> pd.DataFrame:
    """Every evidence item behind the reviewed entries, one row each.

    ``quote`` is verbatim text retrieved from the source, or empty. Items from the
    candidate sweep also carry the second reader's verdict in
    ``second_reader_verdict`` and its note in ``limitations``.
    """
    rows = [
        {"protein": e["protein"], "effect": e["effect"], "rule": e["rule"], **item}
        for e in _literature()["entries"]
        if protein is None or e["protein"] == protein
        for item in e["evidence"]
    ]
    return pd.DataFrame(rows)


def families() -> pd.DataFrame:
    """Family rules: a name pattern, the tissue its members share, tier and sources.

    ``filter_acts`` marks the rules whose tier the candidate filters use. Those rules
    also name the home tissue for the catalog's Other:Unscreened tissue rule.
    """
    frame = pd.DataFrame(_literature()["family_rules"])
    frame["sources"] = frame.sources.map("; ".join)
    return frame


def compartments() -> pd.DataFrame:
    """Tissues added to the tissue-tier registry for these findings, with tier and references."""
    return tissue_tiers.updates()


@cache
def _readthrough() -> pd.DataFrame:
    return pd.read_csv(_DATA / "readthrough-loci.csv", keep_default_na=False)


def rna_ihc_calibration() -> dict:
    """How often a tissue stains at each bulk RNA band, and the grade cutoffs drawn from it."""
    return json.loads((_DATA / "rna-ihc-calibration.json").read_text())


def rna_grade(value: float) -> str:
    """Read an RNA value as "", "low", "medium" or "high" on the staining calibration."""
    floors = rna_ihc_calibration()["grade_floors_ntpm"]
    return next(
        (name for name in ("high", "medium", "low") if value >= floors[name]),
        "",
    )


def nose_skin(protein: str | None = None) -> pd.DataFrame:
    """Normal nose and skin expression of each Target candidate, from two references.

    Olfactory mucosa bulk RNA (four people) and single cells of adult nose, skin,
    scalp and tongue. ``recorded_tissue`` is set when several people or datasets
    agree; :func:`nose_skin_provenance` gives sources, rules and limits.
    """
    d = pd.read_csv(_DATA / "nose-skin-reference.csv", keep_default_na=False)
    return (d[d.protein == protein] if protein else d).reset_index(drop=True)


def nose_skin_provenance() -> dict:
    """Sources, thresholds, coverage and limits of the nose and skin reference."""
    return json.loads((_DATA / "nose-skin-reference.json").read_text())


def readthrough() -> pd.DataFrame:
    """Loci whose every coding transcript joins two separate genes, with the parents."""
    return _readthrough().copy()


def table(catalog: pd.DataFrame) -> pd.DataFrame:
    """Recorded tier per catalog row, from the literature columns the catalog carries.

    ``curated_tier``          most protective tier among the protein's recorded finding,
                              its nose and skin reference record and an acting family
                              rule; "" when none exists
    ``curated_compartments``  the tissues behind it
    ``curated_basis``         "finding", "reference", "family" or several
    ``curated_sources``       source links of the entry and the family rule
    """
    family = families().set_index("id")
    sources = _entries().set_index("entity_id").sources
    rows = []
    for r in catalog.itertuples():
        hits = []
        if isinstance(r.literature_recorded_tier, str) and r.literature_recorded_tier:
            hits.append(
                (
                    r.literature_recorded_tier,
                    r.literature_recorded_tissue,
                    "finding",
                    sources.get(r.entity_id, ""),
                )
            )
        if isinstance(r.reference_recorded_tier, str) and r.reference_recorded_tier:
            hits.append(
                (
                    r.reference_recorded_tier,
                    r.reference_recorded_tissue,
                    "reference",
                    _REFERENCE_SOURCES,
                )
            )
        rule = r.literature_family if isinstance(r.literature_family, str) else ""
        if rule and bool(family.loc[rule, "filter_acts"]):
            hits.append(
                (
                    r.literature_family_tier,
                    family.loc[rule, "tissue"],
                    "family",
                    family.loc[rule, "sources"],
                )
            )
        rows.append(
            {
                "curated_tier": most_protective(h[0] for h in hits),
                "curated_compartments": "; ".join(dict.fromkeys(h[1] for h in hits)),
                "curated_basis": "; ".join(dict.fromkeys(h[2] for h in hits)),
                "curated_sources": " | ".join(dict.fromkeys(h[3] for h in hits if h[3])),
            }
        )
    return pd.DataFrame(rows, index=catalog.index)

# Licensed under the Apache License, Version 2.0.
"""Literature curation of Target candidates: normal function and read-through loci.

Nothing here changes the catalog or the tissue hazard registry. The layer records,
with citations, where a candidate is reported to act in normal adult tissue outside
the expression panels, tiers that tissue (from :mod:`oncoref.tissue_hazards` where the
registry resolves it, otherwise from a cited proposal in ``target-curation.yaml``),
and lists read-through loci. :mod:`oncoref.target_filters` reads :func:`table`.
"""

from __future__ import annotations

import re
from functools import cache
from importlib.resources import files

import pandas as pd
import yaml

TIERS = ("vital", "protected", "manageable", "disposable")  # most protective first
_DATA = files("oncoref") / "data"


@cache
def registry() -> dict:
    """The curation file: compartments, proposed tiers with references, family rules."""
    return yaml.safe_load((_DATA / "target-curation.yaml").read_text())


@cache
def _compartments() -> pd.DataFrame:
    from . import tissue_hazards

    rows = []
    for key, c in registry()["registry_compartments"].items():
        cell = c.get("cell")
        found = (
            tissue_hazards.assess(c["tissue"], cell, resolution="cell")
            if cell
            else tissue_hazards.assess(c["tissue"])
        )
        if found.get("tier") not in TIERS:
            raise ValueError(f"registry does not resolve compartment {key}")
        source = f"tissue_hazards:{found['assessment_id']}"
        rows.append((key, c["label"], found["tier"], source, "", ""))
    for key, c in registry()["proposed_compartments"].items():
        cited = "; ".join(f"PMID {r['pmid']}" for r in c["references"])
        rows.append((key, c["label"], c["tier"], "proposed", c["confidence"], cited))
    columns = ["compartment", "label", "tier", "tier_source", "confidence", "tier_references"]
    return pd.DataFrame(rows, columns=columns)


def compartments() -> pd.DataFrame:
    """Every compartment a finding can point to, with its tier and where the tier comes from."""
    return _compartments().copy()


def most_protective(tiers) -> str:
    """The most protective tier in ``tiers``, or "" when none is given."""
    given = set(tiers)
    return next((t for t in TIERS if t in given), "")


@cache
def _findings() -> pd.DataFrame:
    return pd.read_csv(_DATA / "target-normal-function-findings.csv", keep_default_na=False)


def findings(gene: str | None = None) -> pd.DataFrame:
    """Every literature finding with its verdict, sources, compartments and tier.

    ``acts`` is true only for findings a second reader confirmed as human, adult,
    normal-tissue and gene-specific, and that name a tiered compartment.
    """
    d = _findings()
    return (d[d.gene == gene] if gene else d).reset_index(drop=True)


def families() -> pd.DataFrame:
    """Family rules, each with its pattern, compartment, tier, basis and sources."""
    tier = _compartments().set_index("compartment").tier
    rows = [
        {**f, "tier": tier[f["compartment"]], "sources": "; ".join(f["sources"])}
        for f in registry()["families"]
    ]
    return pd.DataFrame(rows)


@cache
def _readthrough() -> pd.DataFrame:
    return pd.read_csv(_DATA / "target-readthrough-loci.csv", keep_default_na=False)


def readthrough() -> pd.DataFrame:
    """Candidate loci that share exons with two or more other protein-coding genes."""
    return _readthrough().copy()


def table(labels) -> pd.DataFrame:
    """One row per protein label: the curated tier and what it rests on.

    ``curated_tier``          most protective tier among acting findings and family rules
    ``curated_compartments``  the compartments behind it
    ``curated_basis``         "finding", "family" or both
    ``curated_sources``       citations, kept verbatim from the curation records
    ``is_readthrough``        the locus is a read-through of other protein-coding genes
    """
    labels = pd.Index(pd.unique(pd.Series(labels, dtype=str)))
    tier = _compartments().set_index("compartment").tier
    records = {label: [] for label in labels}
    acting = _findings()
    acting = acting[acting.acts.astype(str).str.lower().eq("true")]
    for r in acting.itertuples():
        if r.gene in records:
            for key in r.compartments.split("; "):
                records[r.gene].append((tier[key], key, "finding", r.best_source or r.source))
    for f in registry()["families"]:
        if not f["acts"]:
            continue
        pattern = re.compile(f["pattern"])
        for label in labels:
            if pattern.match(label):
                sources = "; ".join(f["sources"])
                records[label].append((tier[f["compartment"]], f["compartment"], "family", sources))
    loci = set(_readthrough().protein_group_label)
    rows = []
    for label in labels:
        hits = records[label]
        rows.append(
            {
                "protein_group_label": label,
                "curated_tier": most_protective(h[0] for h in hits),
                "curated_compartments": "; ".join(dict.fromkeys(h[1] for h in hits)),
                "curated_basis": "; ".join(dict.fromkeys(h[2] for h in hits)),
                "curated_sources": " | ".join(dict.fromkeys(h[3] for h in hits)),
                "is_readthrough": label in loci,
            }
        )
    return pd.DataFrame(rows).set_index("protein_group_label")

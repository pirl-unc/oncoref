"""Portable explanations and source accounting for the current CTA selection.

All species use the same tables. An unavailable dog/mouse selection remains
explicitly unassessed; a nomination count is never a genome-wide funnel count.
The audit describes the released policy without recomputing classifications.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from dataclasses import dataclass
from importlib.resources import files

import pandas as pd

from .species import Species, resolve_species


@dataclass
class SelectionAudit:
    """Detached, species-scoped audit tables with their policy and source receipt.

    ``explanations`` contains rows per decision layer and entity, with a
    representative observation and a link to the complete evidence where
    available. It is not a complete inventory of every positive observation.
    ``memberships`` counts source annotations, not causal admissions. ``sources``
    counts published candidate-list annotations; ``cta_selection.sources()``
    inventories all measurement sources and their integration status.
    ``criteria`` separates safety, cancer presence, adult restriction and
    reproductive/developmental association. A failure of cancer, restriction
    or association admission is never itself an Unsafe observation.
    ``sources.primary_source`` identifies the compact source display. Redundant
    single-digit lists retain all memberships and their ``covered_by`` sources.
    """

    species: str
    status: str
    policy_id: str | None
    sources: pd.DataFrame
    memberships: pd.DataFrame
    overlaps: pd.DataFrame
    funnel: pd.DataFrame
    criteria: pd.DataFrame
    explanations: pd.DataFrame

    def explain(self, entity: str) -> pd.DataFrame:
        """Return all recorded reasons for an exact entity ID, gene ID or label.

        Display labels are matched case-insensitively. Ambiguous matches fail
        instead of selecting a preferred gene. Use ``entity_id`` to disambiguate.
        """
        frame = self.explanations
        if frame.empty:
            return frame.copy()
        match = frame.entity_id.eq(entity) | frame.label.str.casefold().eq(entity.casefold())
        match |= frame.gene_ids.map(lambda ids: entity in ids)
        found = frame.loc[match]
        if found.entity_id.nunique() > 1:
            raise ValueError(f"Ambiguous entity {entity!r}; use an entity_id")
        if found.empty:
            raise KeyError(f"No audit entity {entity!r} in {self.species}")
        return found.reset_index(drop=True).copy()


def audit(*, species: Species | str | int = "human") -> SelectionAudit:
    """Load the pinned explanation, provenance and funnel tables without network.

    Human counts use protein-equivalence groups. Dog and mouse counts use native assay or
    gene nominations, never a genome-wide denominator. Unknowns are
    retained as nulls, not represented as zero passing or failing genes.
    """
    selected = resolve_species(species)
    root = files("oncoref").joinpath("data")
    receipt = json.loads(root.joinpath("cta-selection-audit.json").read_text())
    if receipt.get("schema") != "oncoref.cta-selection-audit" or receipt.get("version") != 1:
        raise ValueError("Unsupported CTA selection audit schema")
    name = f"cta-selection-audit-{selected.name}.json.gz"
    content = root.joinpath(name).read_bytes()
    if hashlib.sha256(content).hexdigest() != receipt["files"][name]:
        raise ValueError(f"CTA audit checksum mismatch: {name}")
    data = json.loads(gzip.decompress(content))
    if data["taxon_id"] != selected.taxon_id:
        raise ValueError("CTA audit species mismatch")
    tables = {
        key: pd.DataFrame(data[key], columns=receipt["columns"][key]) for key in receipt["columns"]
    }
    for frame in tables.values():
        for column in ("gene_ids", "source_ids", "covered_by"):
            if column in frame:
                frame[column] = frame[column].map(tuple)
    return SelectionAudit(
        species=selected.name,
        status=data["status"],
        policy_id=data["policy_id"],
        **tables,
    )

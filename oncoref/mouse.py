# Licensed under the Apache License, Version 2.0.
"""Mouse data collections; uncurated datasets report unavailable coverage."""

from .species import for_species as _for_species

_data = _for_species(10090)
species = _data.species
references = _data.references
expression = _data.expression
cancer = _data.cancer
antigens = _data.antigens
normalization = _data.normalization
governance = _data.group("governance")
legacy = _data.group("legacy")
groups = (references, expression, cancer, antigens, normalization, governance, legacy)
inventory = _data.inventory
datasets = ()

cta_candidates = antigens.cta_candidates
cta_evidence = antigens.cta_evidence
cta_gene_ids = antigens.cta_gene_ids
cta_snapshot = antigens.cta_snapshot

__all__ = [
    "antigens",
    "cancer",
    "cta_candidates",
    "cta_evidence",
    "cta_gene_ids",
    "cta_snapshot",
    "datasets",
    "expression",
    "governance",
    "groups",
    "inventory",
    "legacy",
    "normalization",
    "references",
    "species",
]

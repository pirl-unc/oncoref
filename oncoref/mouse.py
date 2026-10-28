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

cta_selection = antigens.cta_selection
cta_selection_coverage = antigens.cta_selection_coverage
cta_sources = antigens.cta_sources
cta_candidates = antigens.cta_candidates
cta_evidence = antigens.cta_evidence
cta_gene_ids = antigens.cta_gene_ids
cta_snapshot = antigens.cta_snapshot
cta_tissue_annotations = antigens.cta_tissue_annotations
normal_tissue_scopes = antigens.normal_tissue_scopes
resolve_normal_tissue = antigens.resolve_normal_tissue

__all__ = [
    "antigens",
    "cancer",
    "cta_candidates",
    "cta_evidence",
    "cta_gene_ids",
    "cta_selection",
    "cta_selection_coverage",
    "cta_snapshot",
    "cta_sources",
    "cta_tissue_annotations",
    "datasets",
    "expression",
    "governance",
    "groups",
    "inventory",
    "legacy",
    "normal_tissue_scopes",
    "normalization",
    "references",
    "resolve_normal_tissue",
    "species",
]

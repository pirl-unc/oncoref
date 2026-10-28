# Licensed under the Apache License, Version 2.0.
"""Canine datasets and conceptual data groups, including explicit coverage gaps.

Named datasets are the same modules available directly under ``oncoref``.
Importing this module does not read matrices or download assets.
"""

from . import (
    barkbase,
    bgee_srp009687,
    bgee_srp181218,
    doga,
    gse304066,
    pfizer_canine_protein,
    pfizer_canine_rna,
)
from .species import for_species as _for_species

_data = _for_species(9615)
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

datasets = (
    bgee_srp009687,
    bgee_srp181218,
    barkbase,
    pfizer_canine_protein,
    doga,
    pfizer_canine_rna,
    gse304066,
)

__all__ = [
    "antigens",
    "barkbase",
    "bgee_srp009687",
    "bgee_srp181218",
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
    "doga",
    "expression",
    "governance",
    "groups",
    "gse304066",
    "inventory",
    "legacy",
    "normal_tissue_scopes",
    "normalization",
    "pfizer_canine_protein",
    "pfizer_canine_rna",
    "references",
    "resolve_normal_tissue",
    "species",
]

# Licensed under the Apache License, Version 2.0.
"""Human datasets and conceptual data groups.

Dataset modules are shared with their top-level imports, for example
``human.human_protein_atlas`` is ``oncoref.human_protein_atlas``.
"""

from . import human_protein_atlas, normal_cell_data
from .species import for_species as _for_species

_data = _for_species(9606)
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

# Named dataset modules exposed here; inventory() describes the complete catalog.
datasets = (human_protein_atlas, normal_cell_data)

cta_selection = antigens.cta_selection
cta_selection_coverage = antigens.cta_selection_coverage
cta_selection_audit = antigens.cta_selection_audit
cta_sources = antigens.cta_sources
cta_candidates = antigens.cta_candidates
cta_evidence = antigens.cta_evidence
cta_gene_ids = antigens.cta_gene_ids
cta_snapshot = antigens.cta_snapshot
cta_tissue_annotations = antigens.cta_tissue_annotations
evaluate_target_catalog = antigens.evaluate_target_catalog
host_compartments = antigens.host_compartments
normal_tissue_scopes = antigens.normal_tissue_scopes
resolve_normal_tissue = antigens.resolve_normal_tissue
target_selection = antigens.target_selection
target_selection_coverage = antigens.target_selection_coverage

__all__ = [
    "antigens",
    "cancer",
    "cta_candidates",
    "cta_evidence",
    "cta_gene_ids",
    "cta_selection",
    "cta_selection_audit",
    "cta_selection_coverage",
    "cta_snapshot",
    "cta_sources",
    "cta_tissue_annotations",
    "datasets",
    "evaluate_target_catalog",
    "expression",
    "governance",
    "groups",
    "host_compartments",
    "human_protein_atlas",
    "inventory",
    "legacy",
    "normal_cell_data",
    "normal_tissue_scopes",
    "normalization",
    "references",
    "resolve_normal_tissue",
    "species",
    "target_selection",
    "target_selection_coverage",
]

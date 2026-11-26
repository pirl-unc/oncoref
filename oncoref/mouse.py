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
cta_selection_audit = antigens.cta_selection_audit
cta_sources = antigens.cta_sources
cta_candidates = antigens.cta_candidates
cta_evidence = antigens.cta_evidence
cta_gene_ids = antigens.cta_gene_ids
cta_snapshot = antigens.cta_snapshot
cta_tissue_annotations = antigens.cta_tissue_annotations
evaluate_target_catalog = antigens.evaluate_target_catalog
host_compartments = antigens.host_compartments
host_exposure = antigens.host_exposure
normal_tissue_scopes = antigens.normal_tissue_scopes
resolve_normal_tissue = antigens.resolve_normal_tissue
screening_datasets = antigens.screening_datasets
screening_policy = antigens.screening_policy
screening_sources = antigens.screening_sources
sex_restricted_antigens = antigens.sex_restricted_antigens
target_selection = antigens.target_selection
target_selection_coverage = antigens.target_selection_coverage
tumor_cohorts = cancer.tumor_cohorts
cell_lines = cancer.cell_lines

__all__ = [
    "antigens",
    "cancer",
    "cell_lines",
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
    "host_exposure",
    "inventory",
    "legacy",
    "normal_tissue_scopes",
    "normalization",
    "references",
    "resolve_normal_tissue",
    "screening_datasets",
    "screening_policy",
    "screening_sources",
    "sex_restricted_antigens",
    "species",
    "target_selection",
    "target_selection_coverage",
    "tumor_cohorts",
]

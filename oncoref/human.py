# Licensed under the Apache License, Version 2.0.
"""Human datasets and conceptual data groups.

Dataset modules are shared with their top-level imports, for example
``human.human_protein_atlas`` is ``oncoref.human_protein_atlas``.
"""

from . import ca_only_review, gtex, human_protein_atlas, human_screening_inputs, normal_cell_data
from .ca_only_review import (
    ca_only_class_summary,
    ca_only_functional_review,
    ca_only_niche_adjudications,
)
from .human_screening_inputs import (
    build_human_screening_inputs as screening_inputs,
)
from .human_screening_inputs import (
    evaluate_human_catalog as evaluate_integrated_catalog,
)
from .human_screening_inputs import (
    next_policy_preview,
    next_policy_summary,
    upgraded_737_tumor_observations,
)
from .normal_tissues import normal_tissue_source_grid as tissue_source_grid
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
datasets = (gtex, human_protein_atlas, normal_cell_data)

cell_lines = antigens.cell_lines
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
screening_datasets = antigens.screening_datasets
screening_policy = antigens.screening_policy
sex_restricted_antigens = antigens.sex_restricted_antigens
target_selection = antigens.target_selection
target_selection_coverage = antigens.target_selection_coverage
tumor_cohorts = antigens.tumor_cohorts

__all__ = [
    "antigens",
    "ca_only_class_summary",
    "ca_only_functional_review",
    "ca_only_niche_adjudications",
    "ca_only_review",
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
    "evaluate_integrated_catalog",
    "evaluate_target_catalog",
    "expression",
    "governance",
    "groups",
    "gtex",
    "host_compartments",
    "human_protein_atlas",
    "human_screening_inputs",
    "inventory",
    "legacy",
    "next_policy_preview",
    "next_policy_summary",
    "normal_cell_data",
    "normal_tissue_scopes",
    "normalization",
    "references",
    "resolve_normal_tissue",
    "screening_datasets",
    "screening_inputs",
    "screening_policy",
    "sex_restricted_antigens",
    "species",
    "target_selection",
    "target_selection_coverage",
    "tissue_source_grid",
    "tumor_cohorts",
    "upgraded_737_tumor_observations",
]

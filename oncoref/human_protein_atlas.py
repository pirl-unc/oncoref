# Licensed under the Apache License, Version 2.0.
"""Human Protein Atlas normal-tissue, single-cell and cancer evidence.

These are the existing HPA loaders and caches, exposed together under the full
dataset name. Native units, assay limitations and historical imports are retained.
"""

from .hpa import (
    SAFETY_TISSUE_COVERAGE_LEVELS,
    SAFETY_TISSUE_COVERAGE_STATES,
    SAFETY_TISSUE_MAPPING_KINDS,
    SafetyTissueMapping,
    SafetyTissueResolution,
    SafetyTissueResolutionError,
    gene_cell_type_ntpm,
    gene_protein_tissues,
    gene_tissue_ntpm,
    resolve_safety_tissue_group,
    safety_tissue_mapping_table,
)
from .hpa import (
    hpa_cell_type_expression as cell_type_expression,
)
from .hpa import (
    hpa_normal_tissue as normal_tissue,
)
from .hpa import (
    hpa_normal_tissue_labels as normal_tissue_labels,
)
from .hpa import (
    hpa_rna_consensus as rna_consensus,
)
from .hpa import (
    hpa_single_cell as single_cell,
)
from .hpa import (
    resolve_hpa_normal_tissue_label as resolve_normal_tissue_label,
)
from .hpa_cancer import (
    hpa_cancer_assay_limitations as cancer_assay_limitations,
)
from .hpa_cancer import (
    hpa_cancer_crosswalk as cancer_crosswalk,
)
from .hpa_cancer import (
    hpa_cancer_ihc_prevalence as cancer_ihc_prevalence,
)
from .hpa_cancer import (
    hpa_cancer_rna_cohorts as cancer_rna_cohorts,
)
from .hpa_cancer import (
    hpa_cancer_rna_ihc_comparison as cancer_rna_ihc_comparison,
)
from .hpa_cancer import (
    hpa_cancer_rna_prevalence as cancer_rna_prevalence,
)
from .hpa_cancer import (
    hpa_cancer_sources as cancer_sources,
)

__all__ = [
    "SAFETY_TISSUE_COVERAGE_LEVELS",
    "SAFETY_TISSUE_COVERAGE_STATES",
    "SAFETY_TISSUE_MAPPING_KINDS",
    "SafetyTissueMapping",
    "SafetyTissueResolution",
    "SafetyTissueResolutionError",
    "cancer_assay_limitations",
    "cancer_crosswalk",
    "cancer_ihc_prevalence",
    "cancer_rna_cohorts",
    "cancer_rna_ihc_comparison",
    "cancer_rna_prevalence",
    "cancer_sources",
    "cell_type_expression",
    "gene_cell_type_ntpm",
    "gene_protein_tissues",
    "gene_tissue_ntpm",
    "normal_tissue",
    "normal_tissue_labels",
    "resolve_normal_tissue_label",
    "resolve_safety_tissue_group",
    "rna_consensus",
    "safety_tissue_mapping_table",
    "single_cell",
]

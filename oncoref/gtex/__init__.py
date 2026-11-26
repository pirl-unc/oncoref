# Licensed under the Apache License, Version 2.0.
"""Unified GTEx normal-tissue RNA package across versions, providers, QC, cleanup, & figures (#778).

Submodules:
- :mod:`oncoref.gtex.spec`: :class:`GtexSpec`, ``ALL_QC``, ``QC_ALIASES``, ``ALL_CLEANUP``,
  ``CLEANUP_ALIASES``, ``resolve_spec``, and deterministic Parquet cache helpers.
- :mod:`oncoref.gtex.deconvolution`: 52-shard GTEx v11 deconvolution, donor chromosomal sex,
  organ lineage/surgical-boundary topology, and packaged audit tables.
- :mod:`oncoref.gtex.loaders`: ``tissue_rna``, ``sample_tpm``, ``tissue_summary``, ``tissues``,
  ``versions``, ``donor_sex``, and ``coverage_grid``.
- :mod:`oncoref.gtex.plots`: publication-ready diagnostic and before/after evaluation figures.
"""

from __future__ import annotations

from . import deconvolution, loaders, plots, spec
from . import deconvolution as gtex_deconvolution
from .loaders import (
    apply_summary_spec,
    coverage_grid,
    donor_age,
    donor_sex,
    sample_tpm,
    tissue_rna,
    tissue_summary,
    tissues,
    versions,
)
from .spec import (
    ALL_CLEANUP,
    ALL_QC,
    CELL_LINE_COLUMNS,
    CLEANUP_ALIASES,
    GTEX_V11_SUBCOMPARTMENT_COLUMNS,
    HPA_BRAIN_GTEX_SOURCE_ID,
    HPA_TISSUE_GTEX_SOURCE_ID,
    HPA_UNIQUE_GTEX_TISSUES,
    QC_ALIASES,
    SHARD_SOURCE_ID,
    SOURCE_ID,
    GtexSpec,
    cached_specs,
    clear_derived_cache,
    derived_cache_dir,
    expand_cleanup_ops,
    expand_qc_ops,
    normalize_column_to_site,
    resolve_spec,
    spec_cache_key,
)

__all__ = [
    "ALL_CLEANUP",
    "ALL_QC",
    "CELL_LINE_COLUMNS",
    "CLEANUP_ALIASES",
    "GTEX_V11_SUBCOMPARTMENT_COLUMNS",
    "HPA_BRAIN_GTEX_SOURCE_ID",
    "HPA_TISSUE_GTEX_SOURCE_ID",
    "HPA_UNIQUE_GTEX_TISSUES",
    "QC_ALIASES",
    "SHARD_SOURCE_ID",
    "SOURCE_ID",
    "GtexSpec",
    "apply_summary_spec",
    "cached_specs",
    "clear_derived_cache",
    "coverage_grid",
    "deconvolution",
    "derived_cache_dir",
    "donor_age",
    "donor_sex",
    "expand_cleanup_ops",
    "expand_qc_ops",
    "gtex_deconvolution",
    "loaders",
    "normalize_column_to_site",
    "plots",
    "resolve_spec",
    "sample_tpm",
    "spec",
    "spec_cache_key",
    "tissue_rna",
    "tissue_summary",
    "tissues",
    "versions",
]

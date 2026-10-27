"""Experimental, uniform quantitative CTA/developmental candidate selection.

No publication, clinical-control, gene-family or gene-name admission routes.
This research view does not change the released expression/sequence policy.
Thresholds are exploratory and benchmark-informed, not validated safety limits.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

POLICY_ID = "oncoref.cta-quantitative.2026-10-26.1"


@dataclass(frozen=True)
class Thresholds:
    """Native-unit thresholds for the frozen 2026-10-26 comparison policy."""

    reproductive_ntpm: float = 5
    developmental_tpm: float = 5
    adult_bulk_ntpm: float = 10
    adult_cell_ntpm: float = 30
    immune_native: float = 5
    cornea_cp10k: float = 1

    def __post_init__(self):
        if any(isinstance(v, bool) or not np.isfinite(v) or v <= 0 for v in asdict(self).values()):
            raise ValueError("Thresholds must be positive and finite")


def classify(facts: pd.DataFrame, *, thresholds: Thresholds | None = None) -> pd.DataFrame:
    """Require quantitative association, tumor expression, coverage and adult gates.

    Bulk reproductive dominance compares nTPM within the same HPA assay only.
    Developmental TPM, tumor pTPM, immune assays and corneal CP10K are never
    divided by or averaged with bulk nTPM. Source-specific support statistics
    and exact frozen panel sizes are constructed by the reproducible builder.
    NaN supported maxima with observed source rows mean no supported signal;
    missing source rows fail coverage, never imply measured zero.
    """
    t = thresholds or Thresholds()
    d = facts.copy()
    measurements = [
        "reproductive_ntpm",
        "adult_bulk_max_ntpm",
        "developmental_mean_tpm",
        "developmental_cells",
        "developmental_positive_cells",
        "developmental_fraction",
        "normal_bulk_tissues",
        "cell_records",
        "cell_supported_records",
        "cell_max_ntpm",
        "hpa_immune_cell_records",
        "hpa_immune_cell_max_nTPM",
        "hpa_immune_monaco_records",
        "hpa_immune_monaco_max_TPM",
        "dice_sources",
        "dice_third_max_tpm",
        "cornea_records",
        "cornea_second_max_cp10k",
        "tumor_floor_ptpm",
        "tumor_minimum_samples",
        "tumor_minimum_prevalence",
    ]
    required = {
        *measurements,
        "gene_id",
        "taxonomy_id",
        "reference",
        "tumor_positive",
        "ihc_detected",
    }
    if not required <= set(d):
        raise ValueError(f"Missing quantitative facts: {sorted(required - set(d))}")
    if d.empty or d.gene_id.isna().any() or not d.gene_id.is_unique:
        raise ValueError("One row per canonical gene is required")
    if not d.taxonomy_id.eq(9606).all() or d.reference.isna().any() or d.reference.nunique() != 1:
        raise ValueError("A single explicit human reference is required")
    values = d[measurements].to_numpy(dtype=float)
    if np.isinf(values).any() or (values < 0).any():
        raise ValueError("Measurements must be nonnegative, finite or missing")
    for col in ["tumor_positive", "ihc_detected"]:
        if d[col].isna().any() or not pd.api.types.is_bool_dtype(d[col]):
            raise ValueError(f"Explicit boolean required: {col}")
    if not (
        d.tumor_floor_ptpm.eq(1)
        & d.tumor_minimum_samples.eq(3)
        & d.tumor_minimum_prevalence.eq(0.05)
    ).all():
        raise ValueError("Tumor summaries must use 1 pTPM, 3 samples and 5% within one cohort")
    d["reproductive_route"] = d.reproductive_ntpm.ge(t.reproductive_ntpm) & d.reproductive_ntpm.gt(
        d.adult_bulk_max_ntpm
    )
    d["developmental_route"] = (
        d.developmental_mean_tpm.ge(t.developmental_tpm)
        & d.developmental_cells.ge(20)
        & d.developmental_positive_cells.ge(5)
        & d.developmental_fraction.ge(0.1)
    )
    d["positive_context"] = d.reproductive_route | d.developmental_route
    d["rna_coverage_complete"] = (
        d.normal_bulk_tissues.eq(50)
        & d.cell_records.eq(426)
        & d.hpa_immune_cell_records.eq(19)
        & d.hpa_immune_monaco_records.eq(30)
        & d.dice_sources.eq(15)
        & d.cornea_records.eq(84)
        & d.reproductive_ntpm.notna()
        & d.adult_bulk_max_ntpm.notna()
        & d.hpa_immune_cell_max_nTPM.notna()
        & d.hpa_immune_monaco_max_TPM.notna()
        & d.dice_third_max_tpm.notna()
        & (d.cell_supported_records.eq(0) | d.cell_max_ntpm.notna())
    )
    d["adult_bulk_exclusion"] = d.adult_bulk_max_ntpm.gt(t.adult_bulk_ntpm)
    d["adult_cell_exclusion"] = d.cell_max_ntpm.ge(t.adult_cell_ntpm)
    d["immune_exclusion"] = (
        d.hpa_immune_cell_max_nTPM.ge(t.immune_native)
        | d.hpa_immune_monaco_max_TPM.ge(t.immune_native)
        | d.dice_third_max_tpm.ge(t.immune_native)
    )
    d["cornea_exclusion"] = d.cornea_second_max_cp10k.ge(t.cornea_cp10k)
    gates = [
        ("positive_context", True, "No qualifying reproductive/developmental expression"),
        ("tumor_positive", True, "No qualifying measured tumor expression"),
        ("adult_bulk_exclusion", False, "Adult bulk RNA above limit"),
        ("ihc_detected", False, "Detected adult protein; uncertain staining also held"),
        ("adult_cell_exclusion", False, "Supported adult-cell RNA at or above limit"),
        ("immune_exclusion", False, "Immune RNA at or above limit"),
        ("cornea_exclusion", False, "Recurrent adult corneal RNA at or above limit"),
        ("rna_coverage_complete", True, "Incomplete required RNA panels"),
    ]
    d["selected"] = True
    d["first_failure"] = "Selected"
    for col, expected, reason in gates:
        failure = ~d[col].eq(expected)
        d.loc[d.selected & failure, "first_failure"] = reason
        d["selected"] &= ~failure
    d["policy"] = POLICY_ID
    return d


def results(*, selected_only: bool = False) -> pd.DataFrame:
    """Frozen quantitative proposal, with earlier risk labels kept as annotations."""
    d = pd.read_csv(Path(__file__).parent / "data/cta-quantitative-gene-sets.csv.gz")
    return d.loc[d.selected].copy() if selected_only else d


def gene_ids() -> frozenset[str]:
    """Canonical gene loci selected by the frozen quantitative comparison."""
    return frozenset(results(selected_only=True).gene_id)

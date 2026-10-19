"""Explain the existing HPA CTA gates without changing candidate membership."""

from __future__ import annotations

import pandas as pd

from .cta import cta_gene_ids, cta_unfiltered_gene_ids, passes_filters_mask
from .cta_regen import _has_protein
from .cta_sources import BRADLEY_GENES
from .cta_tissues import HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS
from .load_dataset import get_data


def cta_gate_audit(table=None) -> pd.DataFrame:
    """One row per candidate with independent gate results and exact reasons.

    Protein reliability indicates whether IHC evidence is available, not its
    cancer specificity. Missing RNA fails the RNA gate. The recomputed HPA
    conjunction is checked against the stored gate; default membership remains
    owned by the public CTA helper and can add separate specificity decisions.
    """
    raw = get_data("cancer-testis-antigens") if table is None else table.copy()
    rows = raw.copy()
    rows["coding_gate_pass"] = rows.biotype.eq("protein_coding")
    rows["has_protein_evidence"] = rows.protein_reliability.map(_has_protein)
    rows["protein_gate_pass"] = ~rows.has_protein_evidence | rows.protein_reproductive.astype(
        str
    ).str.strip().str.lower().eq("true")
    tiers = (
        rows.protein_reliability.astype(str).str.strip().where(rows.has_protein_evidence, "Missing")
    )
    rows["required_rna_fraction"] = tiers.map(HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS).fillna(
        HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS["Missing"]
    )
    fraction = pd.to_numeric(rows.rna_deflated_reproductive_frac, errors="coerce")
    rows["rna_measured"] = fraction.notna()
    rows["rna_gate_pass"] = fraction.ge(rows.required_rna_fraction)
    rows["hpa_gate_recomputed"] = (
        rows.coding_gate_pass & rows.protein_gate_pass & rows.rna_gate_pass
    )
    rows["hpa_gate_recorded"] = passes_filters_mask(raw)
    if not rows.hpa_gate_recomputed.equals(rows.hpa_gate_recorded):
        raise ValueError("CTA gate audit does not reproduce the stored HPA filter")
    rows["family_eligible"] = rows.Ensembl_Gene_ID.isin(cta_unfiltered_gene_ids())
    rows["default_panel"] = rows.Ensembl_Gene_ID.isin(cta_gene_ids())

    def reasons(row):
        failed = []
        if not row.coding_gate_pass:
            failed.append("not_protein_coding")
        if not row.protein_gate_pass:
            failed.append("IHC_outside_allowed_reproductive_tissues")
        if not row.rna_measured:
            failed.append("missing_RNA_fraction")
        elif not row.rna_gate_pass:
            failed.append("RNA_fraction_below_protein_tier_threshold")
        if not row.family_eligible:
            failed.append("excluded_gene_family")
        if row.hpa_gate_recomputed and row.family_eligible and not row.default_panel:
            failed.append("default_expression_or_specificity_policy")
        return ";".join(failed) or "retained"

    rows["gate_reasons"] = rows.apply(reasons, axis=1)
    return rows


def bradley_candidate_audit() -> pd.DataFrame:
    """Bradley Figure 3 candidates, in the published order, through current gates."""
    all_rows = cta_gate_audit().set_index("Symbol", drop=False)
    missing = set(BRADLEY_GENES) - set(all_rows.index)
    if missing:
        raise ValueError(f"Bradley nominations absent from the source table: {sorted(missing)}")
    return all_rows.loc[list(BRADLEY_GENES)].reset_index(drop=True)

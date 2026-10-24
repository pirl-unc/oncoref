# Licensed under the Apache License, Version 2.0.
"""An absolute normal-RNA exclusion, independent of reproductive enrichment.

This is a conservative research-background rule, not evidence of clinical
toxicity. Unresolved anatomy cannot establish reproductive/disposable permission.
TPM and nTPM retain their units; counts and CP10K cannot enter this rule.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

POLICY_ID = "oncoref.target-sets.2026-10-24.1"
BASE_POLICY_ID = "oncoref.target-sets.2026-10-23.2"
NORMAL_RNA_CEILING = 30.0
PERMITTED_CATEGORIES = {"reproductive", "disposable"}


def followup_reviews():
    """Cited decisions and reconsideration questions for watched candidates.

    These records document review needs; they never override classification.
    """
    from .load_dataset import get_data

    return get_data("cta-followup-reviews").copy()


def sequence_gene_reviews():
    """Exact reference loci outside the canonical gene table, with RNA decisions."""
    from .load_dataset import get_data

    return get_data("sequence-gene-rna-reviews").copy()


def forbidden_gene_ids():
    """Unsafe canonical genes plus high-RNA reference loci used by protein exports."""
    from .target_sets import gene_ids

    extra = sequence_gene_reviews()
    return gene_ids("Unsafe") | set(extra.loc[extra.Category.eq("Unsafe"), "Gene ID"])


def findings(observations):
    """Annotate mapped adult-human RNA observations without averaging away peaks.

    A single observed value >30 is sufficient. Read counts, donor IDs and
    technical limitations remain in the evidence; they do not silently waive
    the exclusion. An observation exactly at 30 does not trigger this rule.
    """
    d = observations.copy()
    required = {
        "gene_id",
        "source",
        "version",
        "taxonomy_id",
        "modality",
        "stage",
        "unit",
        "value",
        "category",
        "mapping_status",
    }
    if not required <= set(d):
        raise ValueError(f"Missing RNA evidence columns: {sorted(required - set(d))}")
    if not d.taxonomy_id.eq(9606).all() or not d.stage.eq("adult").all():
        raise ValueError("Only adult human evidence supports this RNA ceiling")
    if not d.modality.eq("rna").all() or not d.unit.isin(["TPM", "nTPM"]).all():
        raise ValueError("Use RNA in explicit native TPM or nTPM units; never CP10K/counts")
    if d[["gene_id", "source", "version"]].isna().any().any() or (
        d[["gene_id", "source", "version"]].eq("").any().any()
    ):
        raise ValueError("Gene identity and source/version provenance are required")
    d["value"] = pd.to_numeric(d.value, errors="raise")
    if not np.isfinite(d.value).all() or d.value.lt(0).any():
        raise ValueError("Missing/invalid measurements must not become observed zero")
    d["permitted"] = d.mapping_status.eq("mapped") & d.category.isin(PERMITTED_CATEGORIES)
    d["above_ceiling"] = d.value.gt(NORMAL_RNA_CEILING)
    d["exclude"] = d.above_ceiling & ~d.permitted
    return d


def classify(base, summary):
    """Overlay a complete RNA audit on frozen v1.17 decisions; never rescue genes."""
    if not base.Policy.eq(BASE_POLICY_ID).all():
        raise ValueError("Supply the frozen v1.17 normal-cell decisions")
    if base["Gene ID"].duplicated().any() or summary["Gene ID"].duplicated().any():
        raise ValueError("Unique gene identities required")
    maxima = summary[["Normal RNA maximum (TPM)", "Normal RNA maximum (nTPM)"]].replace("", np.nan)
    maxima = maxima.apply(pd.to_numeric, errors="raise")
    if maxima.lt(0).any().any() or np.isinf(maxima.to_numpy()).any():
        raise ValueError("Invalid normal RNA maxima")
    if not maxima.gt(NORMAL_RNA_CEILING).any(axis=1).eq(summary["High normal RNA"]).all():
        raise ValueError("High RNA findings disagree with the recorded native-unit maxima")
    required_ids = set(base.loc[base.Category.ne("Unsafe"), "Gene ID"])
    assessed_ids = set(summary.loc[summary["Normal RNA assessed"].eq(True), "Gene ID"])
    if not required_ids <= assessed_ids:
        raise ValueError("Every CTA and Neither gene requires an explicit RNA audit")
    d = base.merge(summary, on="Gene ID", how="left", validate="one_to_one")
    d["Before cell review category"] = d["Previous category"]
    d["Previous category"], d["Previous detail"] = d.Category, d.Detail
    d["Previous policy"] = d.Policy
    d["Normal RNA assessed"] = d["Normal RNA assessed"].fillna(False).astype(bool)
    d["High normal RNA"] = d["High normal RNA"].fillna(False).astype(bool)
    d["High normal RNA findings"] = d["High normal RNA findings"].fillna("")
    # Retain existing protein/safety decisions. High RNA also overrides Variable.
    changed = d["High normal RNA"] & d.Category.ne("Unsafe")
    d.loc[changed, ["Category", "Detail", "Reason"]] = [
        "Unsafe",
        "High RNA",
        "Observed RNA above 30 TPM or nTPM outside explicitly permitted reproductive/disposable contexts; research exclusion, not proven clinical toxicity",
    ]
    d["Target"], d["Forbidden"] = d.Category.eq("CTA"), d.Category.eq("Unsafe")
    d["Policy"] = POLICY_ID
    d["Change"] = np.where(
        d.Category.eq(d["Previous category"]),
        "Unchanged",
        d["Previous category"] + " → " + d.Category,
    )
    if (d.Category.isin(["CTA", "Neither"]) & d["High normal RNA"]).any():
        raise AssertionError("High normal RNA cannot remain in CTA or Neither")
    return d.sort_values(["Category", "Gene ID"]).reset_index(drop=True)

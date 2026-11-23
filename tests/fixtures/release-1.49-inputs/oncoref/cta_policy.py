# Licensed under the Apache License, Version 2.0.
"""Current CTA decision hierarchy over explicit measurements and reviewed facts.

The numerical RNA/IHC gates are versioned inputs. This module does not count
papers, turn missing measurements into zero, or contain gene-name exceptions.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

POLICY_ID = "oncoref.target-sets.2026-10-25.2"
BASE_POLICY_ID = "oncoref.target-sets.2026-10-24.1"
SCHEMA = "oncoref.cta-policy.v1"
RNA_FLOOR = 5.0
NONPERMITTED = {"vital", "protected", "manageable"}
FINDINGS = {
    "Normal protein",
    "Normal RNA",
    "Reproductive protein",
    "Reproductive expression",
    "Tumor expression",
    "Tumor recognition",
}
FLAGS = [
    "Nominated",
    "Restriction passed",
    "Coverage adequate",
    "Safety veto",
    "Normal IHC positive",
    "Normal RNA failure",
    "High normal RNA",
    "Variable",
    "Antibody concern",
    "Cell RNA concern",
    "Reproductive protein",
]
REVIEW_COLUMNS = [
    "Review ID",
    "Gene ID",
    "Gene",
    "Taxonomy ID",
    "Reference",
    "Finding",
    "Status",
    "Tissue category",
    "Context",
    "Identity scope",
    "Assay",
    "Specificity",
    "Source",
    "Study",
    "Reviewed by",
    "Review date",
    "Review depth",
    "Rationale",
    "Limitations",
    "Issue",
]


def findings():
    """Source-specific Accepted/Unresolved/Rejected findings; not outcome overrides."""
    from .load_dataset import get_data

    return get_data("cta-policy-findings").fillna("").copy()


def manifest():
    """Active policy, exact guardrails, scopes and current artifact identities."""
    return json.loads((Path(__file__).parent / "data/cta-policy.json").read_text())


def adjudications():
    """Completed primary-antibody reviews, including unresolved biological conflicts."""
    from .load_dataset import get_data

    return get_data("cta-antibody-adjudications").fillna("").copy()


def _validate(facts, reviews):
    required = {"Gene ID", "Gene", "Taxonomy ID", "Reference", "Reproductive RNA (nTPM)", *FLAGS}
    if not required <= set(facts):
        raise ValueError(f"Missing decision facts: {sorted(required - set(facts))}")
    if (
        facts["Gene ID"].isna().any()
        or facts["Gene ID"].eq("").any()
        or facts["Gene ID"].duplicated().any()
    ):
        raise ValueError("Unique nonempty gene identities required")
    if not facts["Taxonomy ID"].eq(9606).all():
        raise ValueError(
            "This threshold policy requires human evidence; other species need calibration"
        )
    if (
        facts.Reference.isna().any()
        or facts.Reference.eq("").any()
        or facts.Reference.nunique() != 1
    ):
        raise ValueError("One explicit sequence/reference identity required")
    for flag in FLAGS:
        if (
            facts[flag].isna().any()
            or not facts[flag].map(lambda v: isinstance(v, (bool, np.bool_))).all()
        ):
            raise ValueError(f"{flag} requires explicit booleans, not missing or string truthiness")
    rna = pd.to_numeric(facts["Reproductive RNA (nTPM)"], errors="raise")
    if rna.lt(0).any() or np.isinf(rna).any():
        raise ValueError("Invalid reproductive RNA")
    if (facts["Restriction passed"] & (rna.isna() | rna.le(0))).any():
        raise ValueError("Quantitative restriction cannot pass without positive reproductive RNA")
    if not set(REVIEW_COLUMNS) <= set(reviews):
        raise ValueError("Incomplete reviewed-finding schema")
    if reviews.empty:
        return
    if (
        reviews["Review ID"].isna().any()
        or reviews["Review ID"].eq("").any()
        or reviews["Review ID"].duplicated().any()
    ):
        raise ValueError("Unique review identities required; do not duplicate evidence")
    if not reviews.Status.isin(["Accepted", "Unresolved", "Rejected"]).all():
        raise ValueError("Findings must be Accepted, Unresolved or Rejected")
    if not reviews.Finding.isin(FINDINGS).all():
        raise ValueError("Unknown finding; clinical nomination is not a protein observation")
    if not set(reviews["Gene ID"]) <= set(facts["Gene ID"]):
        raise ValueError("Review gene missing from the reference-bound universe")
    if (
        not reviews["Taxonomy ID"].eq(9606).all()
        or not reviews.Reference.eq(facts.Reference.iloc[0]).all()
    ):
        raise ValueError("Review species/reference does not match the measurements")
    for column in [
        "Context",
        "Assay",
        "Specificity",
        "Source",
        "Study",
        "Reviewed by",
        "Review date",
        "Review depth",
        "Rationale",
    ]:
        if reviews[column].isna().any() or reviews[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"Reviewed findings require {column}")
    if not reviews.Source.str.startswith(("https://", "http://")).all():
        raise ValueError("A source URL is required")
    pd.to_datetime(reviews["Review date"], format="%Y-%m-%d", errors="raise")
    accepted = reviews.loc[reviews.Status.eq("Accepted")]
    if accepted.Finding.eq("Normal RNA").any():
        raise ValueError(
            "Qualitative Normal RNA may hold or be rejected; quantitative RNA must use "
            "native-unit measurement gates and cannot become a protein veto"
        )
    if not accepted["Identity scope"].eq("gene").all():
        raise ValueError(
            "Accepted findings require resolved gene attribution; family or protein-only "
            "identity must be mapped and reviewed before this gene-level decision"
        )
    normal = accepted.loc[accepted.Finding.eq("Normal protein")]
    if not normal["Tissue category"].isin(NONPERMITTED).all():
        raise ValueError("Accepted normal-protein veto requires a non-permitted compartment")
    reproductive = accepted.loc[
        accepted.Finding.isin(["Reproductive protein", "Reproductive expression"])
    ]
    if not reproductive["Tissue category"].eq("reproductive").all():
        raise ValueError("Reproductive support requires reproductive anatomy, not Disposable")


def classify(facts, reviews):
    """Apply exclusion > unresolved/insufficient > qualified CTA, retaining reasons.

    ``facts`` contains the explicit pinned measurement-gate results, not a prior
    final category. Findings contribute only their accepted biological meaning.
    Unresolved normal-protein/identity concerns hold; rejected observations do
    not veto or supply admission evidence. New negative evidence cannot silently
    erase an accepted positive from another context.
    """
    _validate(facts, reviews)
    d = facts.copy().reset_index(drop=True)
    accepted = reviews.loc[reviews.Status.eq("Accepted")]

    def has(frame, kinds):
        return d["Gene ID"].isin(frame.loc[frame.Finding.isin(kinds), "Gene ID"])

    normal = has(accepted, {"Normal protein"})
    protein = d["Reproductive protein"] | has(accepted, {"Reproductive protein"})
    support = has(accepted, {"Reproductive protein", "Reproductive expression"}) & has(
        accepted, {"Tumor expression", "Tumor recognition"}
    )
    conflict = has(reviews.loc[reviews.Status.eq("Unresolved")], {"Normal protein", "Normal RNA"})
    rna = pd.to_numeric(d["Reproductive RNA (nTPM)"], errors="raise")
    admission = rna.ge(RNA_FLOOR) | protein | support
    normal_rna = d["Normal RNA failure"] & ~d.Variable
    unsafe = (
        d["Safety veto"] | normal | d["Normal IHC positive"] | d["High normal RNA"] | normal_rna
    )
    hold = conflict | d["Antibody concern"] | d["Cell RNA concern"] | d.Variable
    eligible = d.Nominated & d["Restriction passed"] & d["Coverage adequate"] & admission
    d["Category"] = np.select(
        [unsafe, ~unsafe & ~hold & eligible], ["Unsafe", "CTA"], default="Neither"
    )

    # All gates are evaluated. Precedence selects a readable principal reason.
    rules = [
        (
            d["Safety veto"],
            "Safety",
            "An applicable reviewed safety veto excludes unrestricted targeting",
        ),
        (normal, "Normal", "Accepted published protein in non-permitted normal cells"),
        (
            d["Normal IHC positive"],
            "Protein",
            "Accepted normal IHC is positive in a non-permitted compartment",
        ),
        (
            d["High normal RNA"],
            "High RNA",
            "Observed non-permitted or unsited normal RNA exceeds 30 TPM or nTPM",
        ),
        (normal_rna, "RNA", "Normal RNA exceeds its applicable local category/IHC ceiling"),
        (
            conflict,
            "Conflict",
            "Material normal-expression or identity evidence remains unresolved",
        ),
        (
            d.Variable,
            "Variable",
            "Reviewed vital-tissue donor variability with negative local IHC requires reconsideration",
        ),
        (d["Antibody concern"], "Antibody", "Uncertain normal antibody evidence requires review"),
        (
            d["Cell RNA concern"],
            "Cellular",
            "Normal-cell RNA crosses a source-specific review trigger",
        ),
        (
            ~d["Coverage adequate"],
            "Unassessed",
            "The declared measurement panel or reference coverage is incomplete",
        ),
        (
            ~admission,
            "Low",
            "Insufficient reproductive RNA/protein or reviewed affirmative support for admission",
        ),
        (
            d.Nominated & ~d["Restriction passed"],
            "Selectivity",
            "The pinned restriction/selectivity screen does not pass",
        ),
        (~d.Nominated, "Uncurated", "No source-backed cancer/reproductive candidate nomination"),
    ]
    d["Detail"] = np.select(
        [mask for mask, _, _ in rules], [name for _, name, _ in rules], default="Pass"
    )
    d["Reason"] = np.select(
        [mask for mask, _, _ in rules],
        [reason for _, _, reason in rules],
        default="Affirmative reproductive admission and all applicable normal gates pass",
    )
    supported = d.Category.eq("CTA") & rna.lt(RNA_FLOOR) & (protein | support)
    d.loc[supported, "Detail"] = "Supported"
    d.loc[supported, "Reason"] = (
        "Reproductive protein or reviewed reproductive/tumor evidence satisfies admission below 5 nTPM; normal gates pass"
    )
    d["All reasons"] = [
        "; ".join(reason for mask, _, reason in rules if bool(mask.iloc[i]))
        or "All applicable gates pass"
        for i in range(len(d))
    ]
    d["Accepted normal protein"] = normal
    d["Unresolved normal evidence"] = conflict
    d["Reproductive admission"] = admission
    d["Reviewed affirmative support"] = support
    for column, frame, value in [
        ("Review IDs", reviews, "Review ID"),
        ("Evidence sources", reviews, "Source"),
        ("Accepted evidence", accepted, "Review ID"),
        ("Unresolved evidence", reviews.loc[reviews.Status.eq("Unresolved")], "Review ID"),
        ("Rejected evidence", reviews.loc[reviews.Status.eq("Rejected")], "Review ID"),
    ]:
        values = frame.groupby("Gene ID")[value].agg(lambda v: " | ".join(sorted(set(v))))
        d[column] = d["Gene ID"].map(values).fillna("")
    d["Target"], d["Forbidden"] = d.Category.eq("CTA"), d.Category.eq("Unsafe")
    d["Policy"] = POLICY_ID
    if (d.Category.ne("Unsafe") & (d["High normal RNA"] | normal)).any():
        raise AssertionError("Decisive normal evidence escaped exclusion")
    return d

"""Evidence-weighted CTA research selection with explicit clinical qualification.

This policy admits a gene for target development, not every encoded peptide or
receptor. It preserves the stricter expression decisions as a separate column;
clinical qualification never rewrites observed normal expression as absent.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from . import cta_policy

POLICY_ID = "oncoref.cta-research-selection.2026-10-26.1"
BASE_POLICY_ID = cta_policy.POLICY_ID
BOOLEAN_FIELDS = [
    "Prospective treated cohort",
    "Repeated objective responses",
    "Target-specific evidence",
    "Normal selectivity assessed",
    "Unresolved target-mediated organ toxicity",
]
TEXT_FIELDS = [
    "Gene ID",
    "Gene",
    "Reference",
    "Product",
    "Evidence sources",
    "Normal selectivity sources",
    "Reviewed finding IDs",
    "Findings SHA256",
    "Normal finding scope",
    "Clinical outcome",
    "Adverse events",
    "Review depth",
    "Limitations",
    "Rationale",
    "Review date",
]


def clinical_evidence():
    """Reviewed positive and negative comparators; registration alone is insufficient."""
    return pd.read_csv(
        Path(__file__).parent / "data/cta-research-clinical-evidence.csv", keep_default_na=False
    )


def finding_digest(findings, gene_id):
    """Bind reconsideration to complete active finding contents, not just IDs."""
    selected = findings.loc[
        findings["Gene ID"].eq(gene_id)
        & findings.Finding.isin(["Normal protein", "Normal RNA"])
        & findings.Status.ne("Rejected"),
        cta_policy.REVIEW_COLUMNS,
    ].sort_values("Review ID")
    return hashlib.sha256(selected.to_json(orient="records").encode()).hexdigest()


def classify(facts, findings, clinical):
    """Apply the same evidence gates to all identities, without gene-name branches.

    The existing measured RNA/protein gates remain mandatory. A new route can
    qualify low-level or unresolved published protein findings when prospective
    target-specific clinical responses and independent normal-selectivity review
    support continued development. Every conflicting finding must be explicitly
    considered. Substantial normal protein, clinical organ toxicity, high RNA,
    donor variability and cell-RNA concerns cannot use this route.
    """
    base = cta_policy.classify(facts, findings)
    required = {*TEXT_FIELDS, *BOOLEAN_FIELDS, "Taxonomy ID"}
    if not required <= set(clinical):
        raise ValueError(f"Missing clinical evidence fields: {sorted(required - set(clinical))}")
    if clinical["Gene ID"].duplicated().any() or not set(clinical["Gene ID"]) <= set(
        base["Gene ID"]
    ):
        raise ValueError("Clinical records require unique, reference-bound gene identities")
    for col in BOOLEAN_FIELDS:
        if not clinical[col].map(lambda value: isinstance(value, (bool, np.bool_))).all():
            raise ValueError(f"{col} requires explicit booleans")
    for col in TEXT_FIELDS:
        if clinical[col].isna().any() or clinical[col].astype(str).str.strip().eq("").any():
            raise ValueError(f"Clinical evidence requires {col}")
    if (
        not clinical["Taxonomy ID"].eq(9606).all()
        or not clinical.Reference.eq(facts.Reference.iloc[0]).all()
    ):
        raise ValueError("Clinical species/reference must match the decision facts")
    if (
        not clinical["Normal finding scope"]
        .isin(["Limited", "Unresolved", "Substantial", "None"])
        .all()
    ):
        raise ValueError("Unknown normal-finding scope")
    for col in ["Evidence sources", "Normal selectivity sources"]:
        if not clinical[col].str.startswith("https://").all():
            raise ValueError("Primary source links are required")
    pd.to_datetime(clinical["Review date"], format="%Y-%m-%d", errors="raise")
    considered = {}
    rna_review_holds = set(
        findings.loc[findings.Finding.eq("Normal RNA") & findings.Status.ne("Rejected"), "Gene ID"]
    )
    for row in clinical.to_dict("records"):
        local = findings.loc[findings["Gene ID"].eq(row["Gene ID"])]
        expected = set(
            local.loc[
                local.Finding.isin(["Normal protein", "Normal RNA"]) & local.Status.ne("Rejected"),
                "Review ID",
            ]
        )
        declared = set(row["Reviewed finding IDs"].split("|")) - {"None"}
        if declared != expected:
            raise ValueError("Clinical review must account for every active normal finding")
        if row["Findings SHA256"] != finding_digest(findings, row["Gene ID"]):
            raise ValueError("Normal finding contents changed; clinical review must be refreshed")
        accepted = local.loc[local.Finding.eq("Normal protein") & local.Status.eq("Accepted")]
        if not accepted.empty and row["Normal finding scope"] not in {"Limited", "Substantial"}:
            raise ValueError("Accepted normal protein cannot be labelled absent or only unresolved")
        if expected and row["Normal finding scope"] == "None":
            raise ValueError("Active normal findings cannot be labelled absent")
        considered[row["Gene ID"]] = row
    d = base.copy()
    d["Expression category"] = d.Category
    d["Expression detail"] = d.Detail
    d["Expression reason"] = d.Reason
    d["Expression policy"] = d.Policy
    d["Clinical qualification"] = False
    d["Clinical evidence sources"] = ""
    d["Clinical qualification limits"] = ""
    d["Selection basis"] = np.where(d.Category.eq("CTA"), "Expression evidence", "Not selected")
    for index, row in d.iterrows():
        evidence = considered.get(row["Gene ID"])
        if evidence is None:
            continue
        positive = all(evidence[k] for k in BOOLEAN_FIELDS[:-1])
        no_veto = not evidence["Unresolved target-mediated organ toxicity"]
        measured = (
            row.Nominated
            and row["Restriction passed"]
            and row["Coverage adequate"]
            and not any(
                row[k]
                for k in [
                    "Safety veto",
                    "Normal IHC positive",
                    "Normal RNA failure",
                    "High normal RNA",
                    "Variable",
                    "Antibody concern",
                    "Cell RNA concern",
                ]
            )
        )
        qualified = (
            positive
            and no_veto
            and measured
            and row["Reproductive admission"]
            and row["Gene ID"] not in rna_review_holds
            and evidence["Normal finding scope"] != "Substantial"
        )
        d.at[index, "Clinical evidence sources"] = evidence["Evidence sources"]
        d.at[index, "Clinical qualification limits"] = evidence["Limitations"]
        if not no_veto and row.Category == "CTA":
            d.at[index, "Category"] = "Neither"
            d.at[index, "Detail"] = "Clinical concern"
            d.at[index, "Reason"] = "Unresolved target-mediated organ toxicity requires review"
            d.at[index, "Selection basis"] = "Clinical concern; not selected"
        if qualified:
            d.at[index, "Clinical qualification"] = True
            if row.Category != "CTA":
                d.at[index, "Category"] = "CTA"
                d.at[index, "Detail"] = "Clinically qualified"
                d.at[index, "Reason"] = evidence["Rationale"]
                d.at[index, "Selection basis"] = (
                    "Clinical qualification with normal-expression caveats"
                )
    d["Target"] = d.Category.eq("CTA")
    # A development shortlist must not silently regenerate a normal-protein
    # forbidden background or imply clearance of every whole-protein sequence.
    d = d.rename(columns={"Forbidden": "Expression forbidden", "All reasons": "Expression reasons"})
    d["Policy"] = POLICY_ID
    d["Change"] = np.where(
        d.Category.eq(d["Expression category"]),
        "Unchanged",
        d["Expression category"] + " → " + d.Category,
    )
    return d


def results():
    """All canonical decisions for the research selection, including caveats."""
    return pd.read_csv(
        Path(__file__).parent / "data/cta-research-gene-sets.csv.gz",
        keep_default_na=False,
        na_values={"Gene": [""], "Reproductive RNA (nTPM)": [""]},
    )


def gene_ids():
    """CTA research candidates; consult selection basis before product development."""
    d = results()
    return set(d.loc[d.Target, "Gene ID"])

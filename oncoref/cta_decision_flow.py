"""One public decision flow; detailed assay reasons remain attached to each row."""

from __future__ import annotations

from copy import deepcopy

import pandas as pd

FLOW_ID = "oncoref.cta-decision-flow.2026-11-08.1"
CATEGORIES = ("CTA", "Other", "Unsafe")
OTHER_SUBTYPES = {
    "Review": "Required evidence is missing, conflicting, or blocks admission without establishing Unsafe.",
    "Expression": "Cancer presence, reproductive/developmental association or the applicable concentration rule fails.",
    "Sequence": "A nested protein lacks independent expression evidence, or at least 95% of its 9-mers occur in Unsafe proteins.",
}
_STEPS = [
    {"order": 1, "question": "Normal or applicable clinical exclusion?", "if_yes": "Unsafe"},
    {
        "order": 2,
        "question": "Required evidence missing, conflicting or held?",
        "if_yes": "Other:Review",
    },
    {"order": 3, "question": "Nested or Unsafe-contained sequence?", "if_yes": "Other:Sequence"},
    {
        "order": 4,
        "question": "Cancer, reproductive/developmental or concentration requirement fails?",
        "if_yes": "Other:Expression",
    },
    {"order": 5, "question": "All applicable requirements pass", "if_yes": "CTA"},
]


def definition() -> dict:
    """Return a detached contract used by the classifier, API and reports."""
    return {
        "flow_id": FLOW_ID,
        "categories": list(CATEGORIES),
        "other_subtypes": dict(OTHER_SUBTYPES),
        "steps": deepcopy(_STEPS),
        "translation_rule": (
            "Qualifying coding reference or supported noncanonical ORF, plus a complete approved "
            "representative with conflicting coding evidence resolved. Experimental protein "
            "detection is not required for every annotated coding gene."
        ),
        "annotations_not_categories": [
            "paper coverage",
            "experimental protein detection",
            "biological program",
            "tumor prevalence",
            "surface localization",
            "antigen presentation",
        ],
        "scope": "Genes in the assessed human inventory; no transfer of human evidence to another species.",
    }


def apply(frame: pd.DataFrame) -> pd.DataFrame:
    """Standardize Other reasons without changing any top-level classification.

    The input must already have been assessed by the complete selection policy.
    This function cannot admit a protein or discard an Unsafe finding. Retaining
    the source label makes historical releases auditable under one display scheme.
    """
    required = {"category", "subtype", "classification", "assessment_status"}
    if not required <= set(frame):
        raise ValueError(f"Missing decision fields: {sorted(required - set(frame))}")
    if not frame.category.isin(CATEGORIES).all():
        raise ValueError("Uniform decisions require an assessed CTA, Other or Unsafe category")
    cta = frame.category.eq("CTA")
    if not frame.loc[cta, "assessment_status"].eq("screened").all():
        raise ValueError("CTA requires completed assessment")
    known = {
        "Held",
        "Translation",
        "Nested",
        "Contained",
        "Nonreproductive",
        "Reproductive",
        "Diffuse",
        *OTHER_SUBTYPES,
    }
    if not frame.loc[frame.category.eq("Other"), "subtype"].isin(known).all():
        raise ValueError("Unrecognized Other subtype")
    d = frame.copy()
    if "source_classification" not in d:
        d["source_classification"] = d.classification
    other = d.category.eq("Other")
    review = other & (
        ~d.assessment_status.eq("screened").fillna(False)
        | d.subtype.isin(["Held", "Translation", "Review"])
    )
    sequence = other & ~review & d.subtype.isin(["Nested", "Contained", "Sequence"])
    d["subtype"] = pd.Series("", index=d.index, dtype="string")
    d.loc[other, "subtype"] = "Expression"
    d.loc[sequence, "subtype"] = "Sequence"
    d.loc[review, "subtype"] = "Review"
    d["classification"] = d.category.astype("string")
    d.loc[other, "classification"] = "Other:" + d.loc[other, "subtype"]
    d["decision_flow_id"] = FLOW_ID
    return d

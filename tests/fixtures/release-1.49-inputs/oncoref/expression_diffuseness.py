"""Whole-profile expression contrast, independent of tissue safety and IHC.

This metric averages normal / (normal + reference) across a fixed
anatomical panel. The selection policy applies it to adult/placental Reproductive identity. Native source,
version, assay summary and units must match the positive reference. A review flag
is not a safety exclusion. Coverage bounds are not statistical confidence intervals.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

POLICY_ID = "oncoref.expression-diffuseness.2026-10-29.2"
MAX_DIFFUSENESS = 0.10
HPA_REPRODUCTIVE_BRANCHES = {
    "Testicular": ["testis"],
    "Other male": ["epididymis", "prostate", "seminal vesicle"],
    "Ovarian": ["ovary"],
    "Female tract": ["fallopian tube", "endometrium", "cervix", "vagina"],
    "Mammary": ["breast"],
    "Placental": ["placenta"],
}
PERMISSIVE_TISSUES = frozenset({"thymus"})
THYMUS_REFERENCE = "https://doi.org/10.1126/science.1075958"
KEYS = ["entity_id", "source", "version", "taxonomy_id", "unit", "cohort", "summary"]


def hpa_bulk_panel() -> pd.DataFrame:
    """Fixed 38-context / 29-anatomy HPA v23 comparison panel.

    Include Disposable contexts; exclude Permissive thymus from both the normal
    comparator and the Reproductive reference. Brain subdivisions share one weight; colon and
    rectum share one, as do duodenum and small intestine. Other contexts have equal
    anatomy weight. These are sampling weights, not tissue mass or safety weights.
    This panel describes available HPA bulk anatomy, not exhaustive body coverage.
    """
    anatomy = {
        "brain": [
            "amygdala",
            "basal ganglia",
            "cerebellum",
            "cerebral cortex",
            "choroid plexus",
            "hippocampal formation",
            "hypothalamus",
            "midbrain",
        ],
        "colon": ["colon", "rectum"],
        "small intestine": ["small intestine", "duodenum"],
        **{
            name: [name]
            for name in [
                "adipose tissue",
                "adrenal gland",
                "appendix",
                "bone marrow",
                "esophagus",
                "gallbladder",
                "heart muscle",
                "kidney",
                "liver",
                "lung",
                "lymph node",
                "pancreas",
                "parathyroid gland",
                "pituitary gland",
                "retina",
                "salivary gland",
                "skeletal muscle",
                "skin",
                "smooth muscle",
                "spinal cord",
                "spleen",
                "stomach",
                "thyroid gland",
                "tongue",
                "tonsil",
                "urinary bladder",
            ]
        },
    }
    return pd.DataFrame(
        [
            {"compartment": tissue, "anatomy": organ, "weight": 1 / len(anatomy) / len(tissues)}
            for organ, tissues in anatomy.items()
            for tissue in tissues
        ]
    )


def profile_role(annotation) -> str:
    """Specificity role orthogonal to safety, using resolved species/stage anatomy.

    Adult thymus and its mapped descendants are Permissive for expression identity,
    reflecting tissue-antigen expression during central tolerance. This is not a
    Disposable designation or a safety exemption. Explicit embryonic/fetal stage
    retains Reproductive precedence. Unknown anatomy stays unresolved.
    """
    from .normal_tissues import normal_tissue_taxonomy

    if annotation.mapping_status != "mapped":
        return "unresolved"
    if "reproductive" in annotation.scopes:
        return "reproductive"
    concepts = normal_tissue_taxonomy()["concepts"]
    pending = [annotation.concept_id]
    seen = set()
    while pending:
        concept = pending.pop()
        if concept == "thymus":
            return "permissive"
        if concept not in seen:
            seen.add(concept)
            pending.extend(concepts[concept]["parents"])
    return "normal"


def hpa_bulk_roles() -> pd.DataFrame:
    """All 50 HPA v23 contexts, including explicitly excluded thymus."""
    panel = hpa_bulk_panel().assign(role="normal")
    extra = [
        {"compartment": t, "anatomy": branch, "weight": 0.0, "role": "reproductive"}
        for branch, tissues in HPA_REPRODUCTIVE_BRANCHES.items()
        for t in tissues
    ] + [
        {"compartment": t, "anatomy": t, "weight": 0.0, "role": "permissive"}
        for t in sorted(PERMISSIVE_TISSUES)
    ]
    return pd.concat([panel, pd.DataFrame(extra)], ignore_index=True)


def assess(
    profiles: pd.DataFrame,
    references: pd.DataFrame,
    panel: pd.DataFrame,
    *,
    review_above: float = MAX_DIFFUSENESS,
    minimum_reference: float = 2.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return tissue contributions and source-specific diffuseness summaries.

    Profiles require KEYS, compartment, value. References require KEYS,
    reference_value and reference_complete. Reference values are positive-context
    maxima, aggregated across identical-protein loci only within the same context.
    Incomplete locus measurements must already be missing, not a low-value sum.

    A value of 0.5 means every normal context at the reference level in the uniform
    example; 1/11 means every context tenfold lower. One extreme context contributes
    at most its fixed weight. The score is neither a tissue percentage nor a
    whole-body RNA fraction. Zero/absent references do not imply concentration.

    The 10% identity boundary and 2-unit weak-reference guard are analysis
    settings, not assay detection limits or clinical thresholds. The selection gate
    additionally requires the adult/placental association floor. Tiny references
    receive a descriptive score but are not labelled concentrated/diffuse. An
    embryonic route needs a compatible developmental/adult reference, not a borrowed
    adult reproductive anchor or mixed-unit maximum.
    """
    if not 0 < review_above < 1 or not np.isfinite(minimum_reference) or minimum_reference <= 0:
        raise ValueError("Require review boundary in (0,1) and a positive reference floor")
    for frame, fields in [
        (profiles, {*KEYS, "compartment", "value"}),
        (references, {*KEYS, "reference_value", "reference_complete"}),
        (panel, {"compartment", "weight"}),
    ]:
        if not fields <= set(frame):
            raise ValueError(f"Missing fields: {sorted(fields - set(frame))}")
    p = panel.copy()
    if p.compartment.isna().any() or p.compartment.duplicated().any():
        raise ValueError("Unique panel compartments required")
    if not np.isfinite(p.weight).all() or p.weight.le(0).any() or not np.isclose(p.weight.sum(), 1):
        raise ValueError("Positive fixed panel weights must sum to one")
    d = profiles[[*KEYS, "compartment", "value"]].drop_duplicates().copy()
    r = references[[*KEYS, "reference_value", "reference_complete"]].copy()
    if d.duplicated([*KEYS, "compartment"]).any() or r.duplicated(KEYS).any():
        raise ValueError("Conflicting duplicate measurement or reference")
    if not set(d.compartment) <= set(p.compartment):
        raise ValueError("Unexpected compartment outside the fixed panel")
    for frame, value in [(d, "value"), (r, "reference_value")]:
        if frame[KEYS].isna().any().any() or frame[KEYS].eq("").any().any():
            raise ValueError("Explicit identity and source/assay metadata required")
        if not frame.unit.isin(["nTPM", "TPM", "pTPM", "CP10K"]).all():
            raise ValueError("An explicit linear abundance unit is required")
        frame[value] = pd.to_numeric(frame[value], errors="raise")
        if np.isinf(frame[value]).any() or frame[value].lt(0).any():
            raise ValueError("Finite nonnegative expression or missing values required")
    if r.reference_complete.isna().any() or r.reference_complete.dtype != bool:
        raise ValueError("Explicit boolean reference coverage required")
    full = (
        pd.concat([d[KEYS], r[KEYS]], ignore_index=True)
        .drop_duplicates()
        .merge(p, how="cross")
        .merge(d, on=[*KEYS, "compartment"], how="left", validate="one_to_one")
        .merge(r, on=KEYS, how="left", validate="many_to_one")
    )
    full["reference_complete"] = full.reference_complete.fillna(False).astype(bool)
    valid = full.reference_value.gt(0)
    full["similarity"] = (full.value / (full.value + full.reference_value)).where(valid)
    full["contribution"] = full.weight * full.similarity
    full["observed_weight"] = full.weight.where(full.value.notna(), 0)
    result = full.groupby(KEYS, as_index=False, sort=True).agg(
        reference_value=("reference_value", "first"),
        reference_complete=("reference_complete", "first"),
        observed_contribution=("contribution", "sum"),
        normal_coverage=("observed_weight", "sum"),
    )
    result["missing_weight"] = (1 - result.normal_coverage).clip(0, 1)
    positive = result.reference_value.gt(0)
    complete_reference = positive & result.reference_complete
    result["lower"] = result.observed_contribution.where(complete_reference, 0)
    result["upper"] = (result.observed_contribution + result.missing_weight).where(positive, 1)
    result["diffuseness"] = result.lower.where(
        complete_reference & np.isclose(result.normal_coverage, 1)
    )
    result["review"] = np.select(
        [
            ~positive,
            ~result.reference_complete,
            result.reference_value.lt(minimum_reference),
            result.lower.gt(review_above),
            result.upper.le(review_above),
        ],
        [
            "no_reference_signal",
            "incomplete_reference",
            "weak_reference",
            "diffuse",
            "concentrated",
        ],
        default="incomplete_normal_coverage",
    )
    result["review_above"] = review_above
    result["minimum_reference"] = minimum_reference
    result["policy"] = POLICY_ID
    return full, result

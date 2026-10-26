"""Descriptive normal-expression burden; never an eligibility or injury score.

Weights describe a fixed anatomical panel, not tissue mass. Studies, native
units and donor summaries remain separate. Missing weight is never normalized
away, and reproductive expression is not an input to the calculation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

POLICY_ID = "oncoref.normal-expression-breadth.2026-10-26.1"
GROUP = ["gene_id", "source", "version", "taxonomy_id", "unit", "cohort", "summary"]
TIERS = {"vital", "protected", "manageable"}


def anatomical_panel(*, brain_weight=1.0):
    """Fixed bulk-organ panel; brain regions share one organ's total weight.

    Retina and unmeasured non-retinal eye are disjoint components of one eye
    weight. Whole-organ RNA does not additionally count embedded cell types,
    vessels or subregions. Duodenum is not added to the small-intestine sample.
    The panel explicitly includes organs missing from the HPA consensus table;
    it does not claim exhaustive cell-state or developmental coverage.
    """
    if not np.isfinite(brain_weight) or brain_weight <= 0:
        raise ValueError("A positive fixed brain weight is required")
    organs = {
        "vital": {
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
            "eye": ["retina", "non-retinal eye"],
            **{
                x: [x]
                for x in [
                    "heart muscle",
                    "lung",
                    "spinal cord",
                    "liver",
                    "kidney",
                    "bone marrow",
                    "blood",
                    "adrenal gland",
                    "pituitary gland",
                    "pancreas",
                    "small intestine",
                    "skeletal muscle",
                    "smooth muscle",
                    "peripheral nerve",
                    "inner ear",
                ]
            },
        },
        "protected": {
            x: [x]
            for x in [
                "salivary gland",
                "parathyroid gland",
                "lymph node",
                "esophagus",
                "stomach",
                "colon",
                "rectum",
                "urinary bladder",
                "ureter",
                "cartilage",
                "oral mucosa",
                "nasopharynx",
                "tongue",
            ]
        },
        "manageable": {
            x: [x] for x in ["skin", "spleen", "thymus", "thyroid gland", "adipose tissue"]
        },
    }
    rows = []
    for category, groups in organs.items():
        total = sum(brain_weight if organ == "brain" else 1 for organ in groups)
        for organ, compartments in groups.items():
            organ_weight = brain_weight if organ == "brain" else 1
            for compartment in compartments:
                rows.append(
                    {
                        "compartment": compartment,
                        "organ": organ,
                        "category": category,
                        "weight": organ_weight / total / len(compartments),
                    }
                )
    return pd.DataFrame(rows)


def burden(profiles, panel):
    """Return compartment contributions and category coverage bounds.

    Input columns: GROUP plus compartment, value and ceiling. Exactly duplicated
    observations are collapsed, conflicting duplicates rejected. Different
    studies must use distinct source identifiers and are never pooled. The
    [lower, upper] range represents missing compartments, not a confidence
    interval. An all-missing profile has range [0, 1], not evidence of low burden.
    """
    required = {*GROUP, "compartment", "value", "ceiling"}
    if not required <= set(profiles):
        raise ValueError(f"Missing profile columns: {sorted(required - set(profiles))}")
    if not {"compartment", "category", "weight"} <= set(panel):
        raise ValueError("Panel requires compartment, category and weight")
    p = panel.copy()
    if p.compartment.isna().any() or p.compartment.eq("").any() or p.compartment.duplicated().any():
        raise ValueError("Panel compartments must be unique and nonempty")
    if not p.category.isin(TIERS).all():
        raise ValueError("Only screened normal categories belong in the panel")
    if not np.isfinite(p.weight).all() or p.weight.le(0).any():
        raise ValueError("Weights must be finite and positive")
    if not np.allclose(p.groupby("category").weight.sum(), 1, atol=1e-12, rtol=0):
        raise ValueError("Fixed weights must sum to one within each category")
    d = profiles[[*list(GROUP), "compartment", "value", "ceiling"]].copy()
    if (
        d[[*GROUP, "compartment"]].isna().any().any()
        or d[[*GROUP, "compartment"]].eq("").any().any()
    ):
        raise ValueError("Explicit gene, source, version, unit, summary and compartment required")
    # This implementation uses the human HPA v23 bulk ceilings only. Cell data,
    # other units, species and platforms require an independently specified policy.
    if (
        not d.unit.eq("nTPM").all()
        or not d.version.eq("v23").all()
        or not d.taxonomy_id.eq(9606).all()
    ):
        raise ValueError("Only human HPA v23 nTPM profiles are calibrated here")
    if not d.source.isin(["hpa_rna_consensus", "hpa_tissue_donors"]).all():
        raise ValueError("Unsupported source; do not pool or transfer ceilings")
    d.value = pd.to_numeric(d.value, errors="raise")
    d.ceiling = pd.to_numeric(d.ceiling, errors="raise")
    if np.isinf(d.value).any() or d.value.lt(0).any():
        raise ValueError("Invalid expression; missing is NaN, never zero")
    if not np.isfinite(d.ceiling).all() or d.ceiling.le(0).any():
        raise ValueError("Finite positive native-unit ceilings required")
    d = d.drop_duplicates()
    if d.duplicated([*GROUP, "compartment"]).any():
        raise ValueError("Conflicting duplicate observations; adjudicate before aggregation")
    if not set(d.compartment) <= set(p.compartment):
        raise ValueError("Unexpected compartment; use the fixed panel explicitly")
    grid = d[GROUP].drop_duplicates().merge(p, how="cross")
    full = grid.merge(d, on=[*GROUP, "compartment"], how="left", validate="one_to_one")
    full["observed"] = full.value.notna()
    full["fraction"] = full.value / (full.value + full.ceiling)
    full["contribution"] = full.weight * full.fraction
    full["observed_weight"] = full.weight.where(full.observed, 0.0)
    groups = [*GROUP, "category"]
    result = full.groupby(groups, as_index=False, sort=True).agg(
        lower=("contribution", "sum"),
        observed_weight=("observed_weight", "sum"),
        observed_compartments=("observed", "sum"),
        total_compartments=("compartment", "size"),
    )
    result["missing_weight"] = (1 - result.observed_weight).clip(0, 1)
    result["upper"] = result.lower + result.missing_weight
    result["sampling_uncertainty"] = np.where(
        result.source.eq("hpa_rna_consensus"),
        "Not estimated from consensus summary values",
        "See separate marginal donor intervals; no joint profile interval estimated",
    )
    return full, result


def donor_summaries(donors):
    """Median and q90 for independent donor IDs; do not invent matched donors.

    This helper requires one resolved value per gene/cohort/compartment/donor.
    Duplicate copies do not increase n; conflicting or technical replicate
    values must be resolved upstream. Median requires two independent donors;
    q90 requires ten. No cross-compartment joint interval is implied.
    """
    keys = ["gene_id", "source", "version", "taxonomy_id", "unit", "cohort", "compartment"]
    required = {*keys, "donor_id", "value"}
    if not required <= set(donors):
        raise ValueError("Explicit independent donor identities and cohort required")
    d = donors[list(required)].drop_duplicates()
    if d[[*keys, "donor_id"]].isna().any().any() or d[[*keys, "donor_id"]].eq("").any().any():
        raise ValueError("Missing donor identity or provenance")
    if d.duplicated([*keys, "donor_id"]).any():
        raise ValueError("Resolve technical replicates before counting independent donors")
    if not np.isfinite(d.value).all() or d.value.lt(0).any():
        raise ValueError("Invalid donor measurements")
    records = []
    for identity, group in d.groupby(keys, sort=True):
        n = len(group)
        for name, quantile, minimum in [("median", 0.5, 2), ("q90", 0.9, 10)]:
            low = high = np.nan
            if n >= minimum:
                # Marginal resampling: no invented matching across tissues.
                samples = np.random.default_rng(0).choice(group.value.to_numpy(), size=(2000, n))
                estimates = np.quantile(samples, quantile, axis=1)
                low, high = np.quantile(estimates, [0.025, 0.975])
            records.append(
                dict(zip(keys, identity))
                | {
                    "summary": name,
                    "donors": n,
                    "minimum_donors": minimum,
                    "value": float(group.value.quantile(quantile)) if n >= minimum else np.nan,
                    "status": "Measured" if n >= minimum else "Insufficient independent donors",
                    "bootstrap_lower": low,
                    "bootstrap_upper": high,
                    "bootstrap_resamples": 2000 if n >= minimum else 0,
                    "sampling_note": "Marginal percentile bootstrap; small selected cohorts have unstable population coverage; no joint tissue interval",
                }
            )
    return pd.DataFrame(records)

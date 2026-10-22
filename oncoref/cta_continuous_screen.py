# Licensed under the Apache License, Version 2.0.
"""Continuous CTA expression margins and explicitly supplied joint uncertainty.

This module does not fit a donor model or invent uncertainty from tissue means.
The point screen and statistical assessment are deliberately separate outputs.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from functools import lru_cache

import numpy as np
import pandas as pd

from .normal_tissues import (
    SCREENED_NORMAL_TISSUE_TIERS,
    normal_tissue_screening_tier,
    resolve_normal_tissue,
)

BASE_POLICY_ID = "oncoref.cta-ihc-adaptive-screen.2026-10-21.2"
POLICY_ID = "oncoref.cta-continuous-screen.2026-10-22.1"


@dataclass(frozen=True)
class ContinuousPolicy:
    """Provisional research parameters, not calibrated clinical limits.

    The additive allowance is fixed before examining a gene, independent of its
    measurement uncertainty. The inherited absolute ceilings always still apply.
    The donor-tail check uses a separate posterior for the 95th donor percentile;
    a tail concern causes Review, not a claim of demonstrated clinical toxicity.
    """

    allowance_ntpm: float = 0.1
    probability_cutoff: float = 0.95
    minimum_draws: int = 1000
    donor_quantile: float = 0.95

    def __post_init__(self):
        if not math.isfinite(self.allowance_ntpm) or self.allowance_ntpm < 0:
            raise ValueError("RNA allowance must be finite and nonnegative")
        if not math.isfinite(self.probability_cutoff) or not 0.5 < self.probability_cutoff < 1:
            raise ValueError("Probability cutoff must be strictly between 0.5 and 1")
        if isinstance(self.minimum_draws, bool) or not isinstance(self.minimum_draws, int):
            raise ValueError("Minimum draws must be an integer")
        if self.minimum_draws < 1000:
            raise ValueError("At least 1000 independent posterior draws are required")
        if self.donor_quantile != 0.95:
            raise ValueError("This policy specifies the 95th donor percentile")

    def manifest(self):
        return {
            "policy_id": POLICY_ID,
            "base_policy_id": BASE_POLICY_ID,
            **asdict(self),
            "rule": "N <= C AND N <= ((1-s)/s)*R + allowance_ntpm AND R > N",
            "units": "HPA v23 consensus nTPM",
            "status": "provisional research policy; not clinically validated",
        }


def _inclusive(values, limits):
    # One representable step accommodates algebraic rounding at inclusive bounds.
    # This is not a biological tolerance and does not weaken the strict R > N test.
    return values <= np.nextafter(limits, np.inf)


def _inclusive_margin(values, limits):
    margin = limits - values
    return margin.mask((margin < 0) & _inclusive(values, limits), 0.0)


@lru_cache(maxsize=512)
def _ihc_tier(tissue, cell):
    mapped = resolve_normal_tissue(
        tissue, cell_type=cell, species="human", source_id="hpa_normal_tissue", source_version="v23"
    )
    return normal_tissue_screening_tier(mapped, policy="ihc_adaptive_cta")


def _positive_ihc(audit):
    """Bulk-tissue positivity cannot erase an explicitly permitted cell compartment."""
    disallowed = pd.Series(False, index=audit.index)
    unresolved = pd.Series(False, index=audit.index)
    positive_levels = {"Low", "Medium", "High", "Ascending", "Descending"}
    for row in audit.loc[audit.ihc_evidence_status.eq("positive_ihc")].itertuples():
        try:
            contexts = json.loads(row.supporting_ihc)
            positives = [x for x in contexts if x["ihc_level"] in positive_levels]
            if not positives:
                raise ValueError("Positive IHC status without positive cell evidence")
            tiers = {_ihc_tier(x["tissue"], x["cell_type"]) for x in positives}
        except (AttributeError, KeyError, TypeError, json.JSONDecodeError) as exc:
            raise ValueError(
                "Positive IHC requires usable supporting cell-context evidence"
            ) from exc
        disallowed.loc[row.Index] = bool(tiers & set(SCREENED_NORMAL_TISSUE_TIERS))
        unresolved.loc[row.Index] = "unresolved" in tiers
    return disallowed, unresolved


def continuous_tissue_screen(audit, *, policy=ContinuousPolicy()):
    """Add continuous margins to the existing HPA tissue audit (one or many genes).

    R is the highest reproductive tissue mean; N is this normal tissue's mean;
    C and s are its existing local-IHC-dependent absolute and selectivity limits.
    Disposable expression never contributes to R. Missing values remain unknown.
    The caller still applies admission, coverage, nomination and external holds.
    """
    required = {
        "Ensembl_Gene_ID",
        "tissue",
        "source_name",
        "source_version",
        "screening_policy_id",
        "screening_tier",
        "rna_ntpm",
        "rna_ceiling_ntpm",
        "reproductive_rna_max_ntpm",
        "minimum_reproductive_selectivity",
        "selectivity_status",
        "ihc_evidence_status",
    }
    if missing := required - set(audit.columns):
        raise ValueError(f"Missing audit columns: {sorted(missing)}")
    if (
        not audit.source_name.eq("hpa_rna_consensus").all()
        or not audit.source_version.eq("v23").all()
        or not audit.screening_policy_id.eq(BASE_POLICY_ID).all()
    ):
        raise ValueError("Continuous policy requires the pinned HPA v23 base policy")
    if not audit.index.is_unique or audit.duplicated(["Ensembl_Gene_ID", "tissue"]).any():
        raise ValueError("Duplicate gene/tissue audit rows")
    result = audit.copy()
    for column in [
        "rna_ntpm",
        "rna_ceiling_ntpm",
        "reproductive_rna_max_ntpm",
        "minimum_reproductive_selectivity",
    ]:
        result[column] = pd.to_numeric(result[column], errors="coerce")
    screened = result.screening_tier.isin(SCREENED_NORMAL_TISSUE_TIERS)
    n = result.rna_ntpm
    r = result.reproductive_rna_max_ntpm
    s = result.minimum_reproductive_selectivity
    c = result.rna_ceiling_ntpm
    valid = (
        screened
        & np.isfinite(n)
        & n.ge(0)
        & np.isfinite(r)
        & r.gt(0)
        & np.isfinite(c)
        & c.ge(0)
        & s.gt(0)
        & s.le(1)
        & result.selectivity_status.eq("measured")
    )
    allowed = ((1 - s) / s * r + policy.allowance_ntpm).where(valid)
    absolute = _inclusive(n, c)
    selectivity = _inclusive(n, allowed)
    enrichment = r.gt(n)
    positive, unresolved_ihc = _positive_ihc(result)
    positive &= screened
    unresolved_ihc &= screened
    result["disallowed_positive_ihc"] = positive
    result["unresolved_positive_ihc"] = unresolved_ihc
    result["continuous_policy_id"] = POLICY_ID
    result["allowance_ntpm"] = policy.allowance_ntpm
    result["continuous_selectivity_ceiling_ntpm"] = allowed
    result["continuous_selectivity_margin_ntpm"] = _inclusive_margin(n, allowed)
    result["absolute_margin_ntpm"] = _inclusive_margin(n, c).where(screened)
    result["reproductive_margin_ntpm"] = (r - n).where(valid)
    result["point_result"] = "Unassessed"
    result.loc[result.screening_tier.eq("permitted"), "point_result"] = "Not screened"
    result.loc[valid, "point_result"] = "Pass"
    result.loc[unresolved_ihc, "point_result"] = "Unassessed"
    result.loc[valid & ~(absolute & selectivity & enrichment) | positive, "point_result"] = "Fail"
    reasons = pd.Series("", index=result.index, dtype=object)
    for mask, reason in [
        (positive, "positive_ihc"),
        (unresolved_ihc, "unresolved_positive_ihc"),
        (valid & ~absolute, "absolute_ceiling"),
        (valid & ~selectivity, "continuous_selectivity"),
        (valid & ~enrichment, "reproductive_not_higher"),
        (screened & ~valid, "incomplete_expression_evidence"),
        (~screened & ~result.screening_tier.eq("permitted"), "unresolved_anatomy"),
    ]:
        reasons.loc[mask] += reason + ";"
    result["point_reasons"] = reasons.str.rstrip(";")
    return result


@dataclass
class JointExpressionEvidence:
    """Inputs from a separately fitted, checked model, never a donor bootstrap stub.

    means: one independent posterior draw per row, tissue names as columns, of
    the HPA-v23-compatible consensus tissue means. Reproductive reference selection
    is repeated inside every draw. donor_upper: aligned draws of each normal
    tissue's 95th donor percentile, including biological variation. No multiplication
    of marginal probabilities or independent reshuffling is performed here.

    Provenance must identify the gene, model, input hash and normalization, and
    assert model_checks_passed and independent_draws. These are caller attestations,
    not validation performed by this module. Invalid/unassessed model checks cannot
    produce Included or an RNA-based Excluded decision.
    """

    means: pd.DataFrame
    donor_upper: pd.DataFrame | None
    provenance: dict


def _panel(names, label):
    names = tuple(names)
    if not names or any(not isinstance(x, str) or not x for x in names):
        raise ValueError(f"Supply an independent, nonempty {label} tissue panel")
    if len(set(names)) != len(names):
        raise ValueError(f"Duplicate {label} tissues")
    return names


@lru_cache(maxsize=16)
def _validate_panels(normal, reproductive):
    if set(normal) & set(reproductive):
        raise ValueError("Normal and reproductive panels must be disjoint")
    for tissue in normal:
        mapped = resolve_normal_tissue(
            tissue, species="human", source_id="hpa_rna_consensus", source_version="v23"
        )
        if (
            normal_tissue_screening_tier(mapped, policy="ihc_adaptive_cta")
            not in SCREENED_NORMAL_TISSUE_TIERS
        ):
            raise ValueError(f"Not a mapped screened normal tissue: {tissue}")
    for tissue in reproductive:
        mapped = resolve_normal_tissue(
            tissue, species="human", source_id="hpa_rna_consensus", source_version="v23"
        )
        if "reproductive" not in mapped.scopes:
            raise ValueError(f"Not a mapped reproductive tissue: {tissue}")


def assess_joint_expression(
    audit,
    *,
    normal_tissues,
    reproductive_tissues,
    evidence=None,
    hard_exclusions=(),
    review_reasons=(),
    other_gates_passed=False,
    policy=ContinuousPolicy(),
):
    """Return Included / Excluded / Review and a separate point-estimate result.

    Both tissue panels must come from an independent coverage inventory, not the
    gene's available rows. Missing posterior data is Review. Strong joint evidence
    of failing the mean rules is Excluded. Included also requires the donor-tail
    check; tail concerns and intermediate probabilities are Review. Hard exclusions
    (including positive IHC) always override RNA evidence. This is conditional on
    the supplied model and measured tissue panel, never a probability of safety.

    The caller must explicitly confirm other_gates_passed after checking the full
    IHC table (including tissues without RNA), reproductive admission, coverage,
    anatomy and nomination/family/external reviews. This RNA audit alone cannot
    establish those gates. The default therefore cannot yield Included.
    """
    normal = _panel(normal_tissues, "normal")
    reproductive = _panel(reproductive_tissues, "reproductive")
    _validate_panels(normal, reproductive)
    point = continuous_tissue_screen(audit, policy=policy)
    if point.Ensembl_Gene_ID.nunique() != 1:
        raise ValueError("Supply exactly one gene's audit")
    gene = point.Ensembl_Gene_ID.iloc[0]
    rows = point.loc[point.screening_tier.isin(SCREENED_NORMAL_TISSUE_TIERS)].set_index("tissue")
    if set(rows.index) - set(normal):
        raise ValueError("Independent panel omits screened audit tissues")
    hard = set(hard_exclusions)
    if rows.disallowed_positive_ihc.any():
        hard.add("positive_ihc")
    reasons = set(review_reasons)
    if other_gates_passed is not True:
        reasons.add("other_gene_gates_not_confirmed")
    if rows.unresolved_positive_ihc.any():
        reasons.add("unresolved_positive_ihc")
    if rows.ihc_evidence_status.eq("reviewed_unreliable").any():
        reasons.add("reviewed_unreliable_ihc")
    if set(normal) - set(rows.index):
        reasons.add("missing_normal_rna_contexts")
    if rows.point_result.eq("Unassessed").any():
        reasons.add("incomplete_expression_evidence")
    out = {
        "Ensembl_Gene_ID": gene,
        "point_result": "Fail"
        if rows.point_result.eq("Fail").any()
        else "Unassessed"
        if reasons
        else "Pass",
        "result": "Review",
        "reason": "",
        "probability_mean_rules": None,
        "probability_with_donor_tail": None,
        "posterior_draws": 0,
        "model_id": "",
    }
    if hard:
        return {**out, "result": "Excluded", "reason": ";".join(sorted(hard))}
    if evidence is None:
        reasons.add("donor_uncertainty_unavailable")
        return {**out, "reason": ";".join(sorted(reasons))}
    meta = evidence.provenance
    if (
        meta.get("gene_id") != gene
        or meta.get("normalization") != "HPA-v23-nTPM"
        or meta.get("estimand") != "posterior_tissue_consensus_mean"
        or not isinstance(meta.get("model_id"), str)
        or not meta["model_id"].strip()
        or not isinstance(meta.get("input_sha256"), str)
        or len(meta["input_sha256"]) != 64
        or any(x not in "0123456789abcdef" for x in meta["input_sha256"])
    ):
        raise ValueError(
            "Posterior provenance must identify gene, model, input hash and nTPM scale"
        )
    out["model_id"] = meta["model_id"]
    if meta.get("model_checks_passed") is not True:
        reasons.add("model_checks_not_passed")
    if meta.get("independent_draws") is not True:
        reasons.add("independent_posterior_draws_not_established")
    means = evidence.means
    if not means.index.is_unique or not means.columns.is_unique:
        raise ValueError("Posterior draw IDs and tissue columns must be unique")
    out["posterior_draws"] = len(means)
    if len(means) < policy.minimum_draws:
        reasons.add("insufficient_posterior_draws")
    required = [*normal, *reproductive]
    if set(required) - set(means.columns):
        reasons.add("missing_posterior_tissues")
    if reasons:
        return {**out, "reason": ";".join(sorted(reasons))}
    values = means.loc[:, required].to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Posterior expression must be finite and nonnegative")
    # Repeating consensus values, or bootstrapping all-zero donors, cannot provide
    # uncertainty or establish biological absence.
    if means.loc[:, required].nunique().le(1).any():
        return {**out, "reason": "degenerate_posterior_review"}
    r = means.loc[:, reproductive].max(axis=1).to_numpy()[:, None]
    n = means.loc[:, normal].to_numpy(dtype=float)
    c = rows.loc[list(normal), "rna_ceiling_ntpm"].to_numpy(dtype=float)
    s = rows.loc[list(normal), "minimum_reproductive_selectivity"].to_numpy(dtype=float)
    allowance = (1 - s) / s * r + policy.allowance_ntpm
    passed = (_inclusive(n, c) & _inclusive(n, allowance) & (r > n)).all(axis=1)
    probability = float(passed.mean())
    out["probability_mean_rules"] = probability
    if probability <= 1 - policy.probability_cutoff:
        return {**out, "result": "Excluded", "reason": "joint_mean_rules_failed"}
    if evidence.donor_upper is None:
        return {**out, "reason": "donor_tail_uncertainty_unavailable"}
    tail = evidence.donor_upper
    if meta.get("donor_quantile") != policy.donor_quantile:
        raise ValueError("Donor-tail posterior must describe the 95th donor percentile")
    if not tail.columns.is_unique or not tail.index.equals(means.index):
        raise ValueError("Donor-tail draws must retain exact joint posterior alignment")
    if set(normal) - set(tail.columns):
        return {**out, "reason": "missing_donor_tail_tissues"}
    upper = tail.loc[:, normal].to_numpy(dtype=float)
    if not np.isfinite(upper).all() or (upper < 0).any():
        raise ValueError("Donor-tail expression must be finite and nonnegative")
    if tail.loc[:, normal].nunique().le(1).any():
        return {**out, "reason": "degenerate_donor_tail_review"}
    joint = float((passed & _inclusive(upper, c).all(axis=1)).mean())
    out["probability_with_donor_tail"] = joint
    if joint >= policy.probability_cutoff:
        return {**out, "result": "Included", "reason": "joint_expression_rules_passed"}
    return {
        **out,
        "reason": "donor_tail_review"
        if probability >= policy.probability_cutoff
        else "expression_uncertainty_overlaps_limits",
    }

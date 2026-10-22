"""Continuous boundaries and uncertainty must not manufacture expression absence."""

import json

import numpy as np
import pandas as pd
import pytest

from oncoref.cta_continuous_screen import (
    BASE_POLICY_ID,
    ContinuousPolicy,
    JointExpressionEvidence,
    assess_joint_expression,
    continuous_tissue_screen,
)


def audit(n=0.2, r=1, s=0.75, ceiling=2):
    return pd.DataFrame(
        [
            {
                "Ensembl_Gene_ID": "GENE",
                "tissue": "heart muscle",
                "source_name": "hpa_rna_consensus",
                "source_version": "v23",
                "screening_policy_id": BASE_POLICY_ID,
                "screening_tier": "vital",
                "rna_ntpm": n,
                "reproductive_rna_max_ntpm": r,
                "minimum_reproductive_selectivity": s,
                "rna_ceiling_ntpm": ceiling,
                "selectivity_status": "measured",
                "ihc_evidence_status": "unavailable",
            }
        ]
    )


def posterior(normal=("heart muscle",), count=1000):
    rng = np.random.default_rng(32)
    means = pd.DataFrame({t: rng.uniform(0.1, 0.2, count) for t in normal})
    means["testis"] = rng.uniform(9, 11, count)
    return JointExpressionEvidence(
        means=means,
        donor_upper=means.loc[:, list(normal)] * 2,
        provenance={
            "gene_id": "GENE",
            "normalization": "HPA-v23-nTPM",
            "estimand": "posterior_tissue_consensus_mean",
            "model_id": "synthetic-test-only",
            "input_sha256": "a" * 64,
            "model_checks_passed": True,
            "independent_draws": True,
            "donor_quantile": 0.95,
        },
    )


def assess(rows=None, evidence=None, **kwargs):
    return assess_joint_expression(
        audit() if rows is None else rows,
        normal_tissues=("heart muscle",),
        reproductive_tissues=("testis",),
        evidence=evidence,
        other_gates_passed=True,
        **kwargs,
    )


def test_continuity_at_one_ntpm_and_effect_of_absolute_scale():
    low = continuous_tissue_screen(audit(n=0.08, r=0.1)).iloc[0]
    high = continuous_tissue_screen(audit(n=8, r=10, ceiling=20)).iloc[0]
    assert low.point_result == "Pass"
    assert high.point_result == "Fail"
    a = continuous_tissue_screen(audit(n=0.999999, r=2)).iloc[0]
    b = continuous_tissue_screen(audit(n=1.000001, r=2)).iloc[0]
    assert (
        abs(a.continuous_selectivity_margin_ntpm - b.continuous_selectivity_margin_ntpm) < 0.000003
    )


def test_equality_and_absolute_ceilings_cannot_be_rescued_by_allowance():
    for n, r, ceiling, reason in [
        (0.1, 0.1, 2, "reproductive_not_higher"),
        (2.1, 100, 2, "absolute_ceiling"),
        (0, 0, 2, "incomplete_expression_evidence"),
    ]:
        row = continuous_tissue_screen(audit(n=n, r=r, ceiling=ceiling)).iloc[0]
        assert row.point_result != "Pass"
        assert reason in row.point_reasons
    # Inclusive arithmetic boundary, without rounding a displayed percentage.
    row = continuous_tissue_screen(audit(n=0.6, r=1, s=2 / 3)).iloc[0]
    assert row.point_result == "Pass"


def test_zero_allowance_recovers_ratio_and_positive_ihc_always_vetoes():
    row = continuous_tissue_screen(audit(n=0.4), policy=ContinuousPolicy(allowance_ntpm=0)).iloc[0]
    assert row.point_result == "Fail"
    assert continuous_tissue_screen(audit(n=0.4)).iloc[0].point_result == "Pass"
    rows = audit()
    rows["ihc_evidence_status"] = "positive_ihc"
    rows["supporting_ihc"] = json.dumps(
        [{"tissue": "heart muscle", "cell_type": "cardiomyocytes", "ihc_level": "Low"}]
    )
    assert continuous_tissue_screen(rows).iloc[0].point_result == "Fail"
    assert assess(rows, posterior())["result"] == "Excluded"


def test_permitted_hair_follicle_staining_is_not_a_whole_skin_veto():
    rows = audit().assign(
        tissue="skin", screening_tier="manageable", ihc_evidence_status="positive_ihc"
    )
    rows["supporting_ihc"] = json.dumps(
        [{"tissue": "skin 1", "cell_type": "hair follicles", "ihc_level": "Low"}]
    )
    result = continuous_tissue_screen(rows).iloc[0]
    assert result.point_result == "Pass"
    assert not result.disallowed_positive_ihc
    rows["supporting_ihc"] = json.dumps(
        [{"tissue": "skin 1", "cell_type": "keratinocytes", "ihc_level": "Low"}]
    )
    assert continuous_tissue_screen(rows).iloc[0].point_result == "Fail"


def test_exact_75_percent_boundary_does_not_fail_due_to_float_rounding():
    result = continuous_tissue_screen(
        audit(n=0.1, r=0.3), policy=ContinuousPolicy(allowance_ntpm=0)
    )
    assert result.iloc[0].point_result == "Pass"


@pytest.mark.parametrize("bad", [np.nan, np.inf, -1])
def test_invalid_values_are_unknown_and_bad_parameters_rejected(bad):
    row = continuous_tissue_screen(audit(n=bad)).iloc[0]
    assert row.point_result == "Unassessed"
    with pytest.raises(ValueError):
        ContinuousPolicy(allowance_ntpm=bad)


def test_missing_uncertainty_is_review_even_for_a_passing_point_estimate():
    result = assess()
    assert result["point_result"] == "Pass"
    assert result["result"] == "Review"
    assert result["probability_mean_rules"] is None
    assert result["reason"] == "donor_uncertainty_unavailable"
    failed = assess(audit(n=20))
    assert failed["point_result"] == "Fail"
    assert failed["result"] == "Review"
    assert assess(hard_exclusions=("external_safety_hold",))["result"] == "Excluded"


def test_complete_strong_joint_evidence_and_mean_failure():
    evidence = posterior()
    assert assess(evidence=evidence)["result"] == "Included"
    evidence.means["heart muscle"] += 5
    result = assess(evidence=evidence)
    assert result["result"] == "Excluded"
    assert result["probability_mean_rules"] == 0


def test_mean_precision_does_not_hide_donor_tail():
    evidence = posterior()
    evidence.donor_upper["heart muscle"] += 3
    result = assess(evidence=evidence)
    assert result["probability_mean_rules"] == 1
    assert result["probability_with_donor_tail"] == 0
    assert result["result"] == "Review"
    assert result["reason"] == "donor_tail_review"
    evidence.donor_upper = None
    assert assess(evidence=evidence)["reason"] == "donor_tail_uncertainty_unavailable"


def test_joint_probability_is_not_product_or_minimum_of_tissue_probabilities():
    rows = pd.concat([audit(), audit().assign(tissue="lung")], ignore_index=True)
    evidence = posterior(normal=("heart muscle", "lung"))
    # Each tissue passes 97%, but disjoint failures mean only 94% pass jointly.
    evidence.means.loc[:29, "heart muscle"] = 5
    evidence.means.loc[30:59, "lung"] = 5
    result = assess_joint_expression(
        rows,
        normal_tissues=("heart muscle", "lung"),
        reproductive_tissues=("testis",),
        evidence=evidence,
        other_gates_passed=True,
    )
    assert result["probability_mean_rules"] == 0.94
    assert result["result"] == "Review"


def test_reproductive_reference_is_reselected_per_joint_draw():
    evidence = posterior()
    evidence.means["testis"] *= 0.001
    evidence.means["ovary"] = np.linspace(9, 11, 1000)
    result = assess_joint_expression(
        audit(),
        normal_tissues=("heart muscle",),
        reproductive_tissues=("testis", "ovary"),
        evidence=evidence,
        other_gates_passed=True,
    )
    assert result["result"] == "Included"
    with pytest.raises(ValueError, match="reproductive tissue"):
        assess_joint_expression(
            audit(), normal_tissues=("heart muscle",), reproductive_tissues=("tonsil",)
        )


@pytest.mark.parametrize("defect", ["missing", "constant", "zero", "few", "qc", "dependent"])
def test_inadequate_posterior_never_becomes_included(defect):
    evidence = posterior(count=500 if defect == "few" else 1000)
    if defect == "missing":
        evidence.means = evidence.means.drop(columns="testis")
    elif defect in {"zero", "constant"}:
        evidence.means["heart muscle"] = 0 if defect == "zero" else 0.1
    elif defect == "qc":
        evidence.provenance["model_checks_passed"] = False
    elif defect == "dependent":
        evidence.provenance["independent_draws"] = False
    assert assess(evidence=evidence)["result"] == "Review"


def test_source_units_gene_and_draw_alignment_are_enforced():
    evidence = posterior()
    evidence.provenance["normalization"] = "GTEx TPM"
    with pytest.raises(ValueError, match="provenance"):
        assess(evidence=evidence)
    evidence = posterior()
    evidence.provenance["gene_id"] = "WRONG"
    with pytest.raises(ValueError, match="provenance"):
        assess(evidence=evidence)
    evidence = posterior()
    evidence.donor_upper = evidence.donor_upper.iloc[::-1]
    with pytest.raises(ValueError, match="alignment"):
        assess(evidence=evidence)
    with pytest.raises(ValueError, match="pinned HPA"):
        continuous_tissue_screen(audit().assign(source_version="v24"))


def test_incomplete_panel_and_other_review_flags_prevent_inclusion():
    evidence = posterior(normal=("heart muscle", "lung"))
    result = assess_joint_expression(
        audit(),
        normal_tissues=("heart muscle", "lung"),
        reproductive_tissues=("testis",),
        evidence=evidence,
    )
    assert result["result"] == "Review"
    assert "missing_normal_rna_contexts" in result["reason"]
    assert assess(evidence=evidence, review_reasons=("unresolved_anatomy",))["result"] == "Review"


def test_rna_draws_alone_cannot_confirm_the_other_gene_gates():
    result = assess_joint_expression(
        audit(),
        normal_tissues=("heart muscle",),
        reproductive_tissues=("testis",),
        evidence=posterior(),
    )
    assert result["result"] == "Review"
    assert "other_gene_gates_not_confirmed" in result["reason"]

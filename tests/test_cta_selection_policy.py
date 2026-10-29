"""Single RNA boundary and positive, evidence-backed Neither definitions."""

import pandas as pd
import pytest

from oncoref import cta_selection_policy as policy
from oncoref.cta_funnel import EvidenceLayer


def context(tier="vital", value=4, grade="Missing"):
    return {
        "gene_id": "g",
        "source": "hpa_single_cell_tissue",
        "version": "v23",
        "taxonomy_id": 9606,
        "stage": "adult",
        "state": "normal",
        "tissue": "colon",
        "cell_type": "cell",
        "category": tier,
        "mapping_status": "mapped",
        "negative_ihc_grade": grade,
        "value": value,
        "unit": "nTPM",
    }


@pytest.mark.parametrize("tier,ceiling", [("vital", 4), ("protected", 8), ("manageable", 16)])
def test_one_boundary_equality_and_anatomical_risk(tier, ceiling):
    rows = pd.DataFrame([context(tier, ceiling), context(tier, ceiling + 0.01)])
    assessed = policy.assess_contexts(rows)
    assert assessed.candidate_exclusion.tolist() == [False, True]
    assert assessed.unsafe_evidence.tolist() == [False, tier != "manageable"]
    assert "strong_boundary" not in assessed
    assert (
        policy.assess_contexts(pd.DataFrame([context(tier, 0, "Enhanced")])).candidate_ceiling.iloc[
            0
        ]
        == ceiling * 8
    )


def test_compartment_match_does_not_transfer_to_other_tissue_state_or_species():
    crosswalk = pd.DataFrame(
        [
            {
                "source": "hpa_single_cell_tissue",
                "version": "v23",
                "tissue": "colon",
                "cell_type": "distal enterocytes",
                "ihc_cell_type": "glandular cells",
                "ihc_match_source": "citation",
            }
        ]
    )
    base = context() | {"cell_type": "distal enterocytes"}
    rows = pd.DataFrame(
        [
            base,
            base | {"tissue": "heart muscle"},
            base | {"state": "stimulated"},
            base | {"taxonomy_id": 9615},
            base | {"stage": "fetal"},
            base | {"cell_type": "T-cells"},
        ]
    )
    ihc = pd.DataFrame(
        [
            {
                "Gene": "g",
                "Tissue": "colon",
                "Cell type": "glandular cells",
                "Level": "Not detected",
                "Reliability": "Enhanced",
            }
        ]
    )
    result = policy.match_ihc(rows, ihc, crosswalk=crosswalk, positive_controls={"g"})
    assert result.negative_ihc_grade.tolist() == ["Enhanced", *["Missing/uncertain"] * 5]
    assert result.ihc_match_scope.tolist()[0] == "reviewed_compartment"
    assert not result.ihc_exact_match.any()
    assert result.ihc_match_scope.iloc[-1] == "same_tissue_unresolved_cell"
    conflict = policy.match_ihc(
        rows.iloc[:1], ihc, crosswalk=crosswalk, positive_controls={"g"}, conflicting_genes={"g"}
    )
    assert conflict.negative_ihc_grade.iloc[0] == "Missing/uncertain"


def profiles():
    return pd.DataFrame(
        [
            {
                "entity": str(i),
                "required_coverage": True,
                "unresolved_evidence": False,
                "tumor_measured": True,
                "tumor_profile_complete": True,
                "developmental_measured": True,
                "normal_maximum_floor_ratio": 0.5,
                "reproductive_ntpm": 0.5,
                "developmental_all_lineage_max_tpm": 0.5,
                "tumor_maximum_floor_ratio": 0.5,
                "any_detected_protein": False,
                "positive_context": False,
                "vital_protected_exclusion": False,
                "manageable_tissue_count": 0,
            }
            for i in range(5)
        ]
    )


def test_missingness_failed_gate_and_manageable_excess_are_not_neither():
    facts = profiles()
    facts.loc[1, "required_coverage"] = False
    facts.loc[2, "tumor_maximum_floor_ratio"] = 2  # Failed association alone is insufficient.
    facts.loc[3, ["positive_context", "normal_maximum_floor_ratio"]] = [True, 20]
    facts.loc[4, ["normal_maximum_floor_ratio", "manageable_tissue_count"]] = [20, 2]
    facts = policy.neither_profiles(facts)
    yes, no = pd.Series(True, index=facts.index), pd.Series(False, index=facts.index)
    result = policy.classify(
        facts, [EvidenceLayer("candidate gate", yes, no)], entity_column="entity"
    )
    assert result.classifications.category.fillna("Unassigned").tolist() == [
        "Neither",
        "Unassigned",
        "Unassigned",
        "Unassigned",
        "Neither",
    ]
    assert result.classifications.neither_profile.tolist() == [
        "low_expression",
        "",
        "",
        "",
        "diffuse_manageable",
    ]
    assert (
        result.funnel[["CTA_eligible_so_far", "Unsafe", "Neither", "Unassigned"]]
        .sum(axis=1)
        .eq(5)
        .all()
    )


def test_unsafe_wins_regardless_of_order_and_cannot_be_rescued_by_neither():
    facts = policy.neither_profiles(profiles())
    yes, no = pd.Series(True, index=facts.index), pd.Series(False, index=facts.index)
    layers = [EvidenceLayer("gate", yes, no), EvidenceLayer("vital", yes, yes)]
    first = policy.classify(facts, layers, entity_column="entity").classifications
    second = policy.classify(facts, list(reversed(layers)), entity_column="entity").classifications
    assert first.category.eq("Unsafe").all()
    pd.testing.assert_series_equal(first.category, second.category)


def test_tumor_measurement_is_distinct_from_complete_low_expression_profile():
    facts = profiles().iloc[:1].copy()
    facts["tumor_profile_complete"] = False
    result = policy.neither_profiles(facts)
    assert result.tumor_measured.all()
    assert not result.neither_low_expression.any()


@pytest.mark.parametrize("name", sorted(policy.CANDIDATE_ONLY_LAYERS))
def test_admission_failures_cannot_be_wired_as_unsafe(name):
    facts = policy.neither_profiles(profiles())
    yes, no = pd.Series(True, index=facts.index), pd.Series(False, index=facts.index)
    result = policy.classify(
        facts, [EvidenceLayer(name, yes, no)], entity_column="entity"
    ).classifications
    assert not result.category.eq("Unsafe").any()
    with pytest.raises(ValueError, match="admission requirement, not safety evidence"):
        policy.classify(facts, [EvidenceLayer(name, yes, yes)], entity_column="entity")


def test_normal_safety_is_evaluated_even_after_cancer_admission_fails():
    facts = policy.neither_profiles(profiles())
    yes, no = pd.Series(True, index=facts.index), pd.Series(False, index=facts.index)
    layers = [
        EvidenceLayer("Tumor sample ≥2 TPM/pTPM", yes, no),
        EvidenceLayer("Normal bulk", yes, yes),
    ]
    result = policy.classify(facts, layers, entity_column="entity").classifications
    assert result.category.eq("Unsafe").all()
    assert result.Unsafe_reasons.eq("Normal bulk").all()


def test_association_does_not_depend_on_normal_expression_or_coverage():
    facts = pd.DataFrame(
        {
            "reproductive_ntpm": [5, 5, 4, 0, 100, 0],
            "developmental_mean_tpm": [0, 0, 0, 5, 0, 5],
            "healthy_maximum_ntpm": [4, 5, 0, 100, float("nan"), 0],
            "healthy_maximum_tpm": [100, 0, 0, 4, 0, float("nan")],
        }
    )
    result = policy.association_evidence(facts)
    assert result.positive_context.tolist() == [True, True, False, True, True, True]
    assert result.reproductive_route.tolist() == [True, True, False, False, True, False]
    assert result.developmental_route.tolist() == [False, False, False, True, False, True]
    assert "adult_restriction" not in result


def test_association_failure_does_not_make_expression_unsafe():
    facts = policy.neither_profiles(profiles())
    facts["neither_low_expression"] = False
    yes, no = pd.Series(True, index=facts.index), pd.Series(False, index=facts.index)
    result = policy.classify(
        facts,
        [EvidenceLayer("Reproductive / embryonic association", yes, no)],
        entity_column="entity",
    ).classifications
    assert result.category.isna().all()
    assert result.Unsafe_reasons.eq("").all()


def test_diffuseness_is_route_specific_and_missingness_is_not_failure():
    facts = policy.association_evidence(
        pd.DataFrame(
            {
                "protein_equivalence_group_id": list("abcdef"),
                "reproductive_ntpm": [5, 5, 1, 5, 5, 5],
                "developmental_mean_tpm": [0, 10, 10, 0, 0, 0],
            }
        )
    )
    scores = pd.DataFrame(
        {
            "entity_id": list("abcdef"),
            "source": "hpa_rna_consensus",
            "version": "v23",
            "unit": "nTPM",
            "taxonomy_id": 9606,
            "summary": "consensus",
            "diffuseness": [0.2, 0.2, 0.9, float("nan"), 0.1, 0.2],
            "reference_value": [5, 5, 1, 5, 5, 5],
            "reference_complete": [True, True, True, False, True, True],
            "normal_coverage": [1, 1, 1, 0.5, 1, 0.9],
        }
    )
    result = policy.concentration_evidence(facts, scores)
    assert result.reproductive_concentration_hold.tolist() == [
        True,
        False,
        False,
        False,
        False,
        False,
    ]
    assert result.adult_profile_status.tolist() == [
        "diffuse",
        "diffuse",
        "no_qualifying_adult_reference",
        "unassessed",
        "concentrated",
        "unassessed",
    ]
    assert result.reproductive_route.tolist() == [False, False, False, True, True, True]
    with pytest.raises(ValueError, match="HPA v23"):
        policy.concentration_evidence(facts, scores.assign(source="gtex_tissue_rna"))


def test_diffuse_reproductive_neither_requires_clear_safety_and_complete_coverage():
    facts = policy.neither_profiles(profiles())
    facts[["neither_low_expression", "neither_diffuse_manageable"]] = False
    facts.loc[3, "required_coverage"] = False
    facts.loc[4, "unresolved_evidence"] = True
    yes, no = pd.Series(True, index=facts.index), pd.Series(False, index=facts.index)
    layers = [
        EvidenceLayer("Reproductive concentration", yes, no),
        EvidenceLayer(
            "Normal bulk",
            pd.Series([False, True, True, False, False]),
            pd.Series([False, False, True, False, False]),
        ),
    ]
    first = policy.classify(facts, layers, entity_column="entity").classifications
    second = policy.classify(facts, list(reversed(layers)), entity_column="entity").classifications
    assert first.category.fillna("Unassigned").tolist() == [
        "Neither",
        "Unassigned",
        "Unsafe",
        "Unassigned",
        "Unassigned",
    ]
    assert first.neither_profile.iloc[0] == "diffuse_reproductive"
    pd.testing.assert_series_equal(first.category, second.category)

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

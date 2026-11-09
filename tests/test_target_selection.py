"""Shared membership, comparison-scale and sample-independence contracts."""

import itertools
import json
from html.parser import HTMLParser
from pathlib import Path

import pandas as pd
import pytest

from oncoref import target_selection as policy


def facts(**changes):
    return {
        "entity_id": "protein-A",
        "unsafe": False,
        **dict.fromkeys(policy.SHARED_CHECKS, True),
        "cancer_evidence": "qualifies",
        "reproductive_evidence": "qualifies",
        **changes,
    }


def observations(tumor_values=(2, 3), normal_values=(1, 0.5)):
    common = {
        "entity_id": "protein-A",
        "comparison_id": "panel-v1",
        "unit": "TPM",
        "qc_pass": True,
        "identity_resolved": True,
    }
    tumors = pd.DataFrame(
        [
            dict(
                common,
                source_id="tumor-study",
                sample_id=f"sample-{i}",
                participant_id=f"person-{i}",
                value=v,
            )
            for i, v in enumerate(tumor_values)
        ]
    )
    normals = pd.DataFrame(
        [
            dict(
                common, source_id="normal-study", context_id=f"context-{i}", scope="normal", value=v
            )
            for i, v in enumerate(normal_values)
        ]
    )
    return tumors, normals


def assess(tumors, normals, **kwargs):
    return policy.assess_cancer_rna(
        tumors,
        normals,
        required_normal_contexts={"context-0", "context-1"},
        comparison_id="panel-v1",
        comparability_reference="reviewed-comparison-receipt.json",
        unit="TPM",
        **kwargs,
    )


@pytest.mark.parametrize(
    "unsafe,cancer,reproductive", list(itertools.product([False, True], repeat=3))
)
def test_membership_set_algebra_and_unsafe_precedence(unsafe, cancer, reproductive):
    frame = pd.DataFrame(
        [
            facts(
                unsafe=unsafe,
                cancer_evidence="qualifies" if cancer else "does_not_qualify",
                reproductive_evidence="qualifies" if reproductive else "does_not_qualify",
                published_nomination=True,
                paper_links=["https://doi.org/example"],
            )
        ]
    )
    row = policy.classify(frame).iloc[0]
    assert row.is_cancer_antigen == (cancer and not unsafe)
    assert row.is_reproductive_antigen == (reproductive and not unsafe)
    assert row.is_cta == (row.is_cancer_antigen and row.is_reproductive_antigen)
    assert row.is_target_candidate == (row.is_cancer_antigen or row.is_reproductive_antigen)
    assert (row.category == "Unsafe") == unsafe
    assert row.published_nomination and row.paper_links == ["https://doi.org/example"]
    assert "category" not in frame  # Inputs remain available for other policy replays.


@pytest.mark.parametrize("check", policy.SHARED_CHECKS)
def test_every_shared_check_blocks_both_sets_without_creating_unsafe(check):
    row = policy.classify(pd.DataFrame([facts(**{check: False})])).iloc[0]
    assert row.category == "Other"
    assert not row.is_target_candidate
    assert row.subtype == ("Sequence" if check == "sequence_pass" else "Review")
    row = policy.classify(pd.DataFrame([facts(unsafe=True, **{check: False})])).iloc[0]
    assert row.category == "Unsafe"


def test_unassessed_association_is_not_negative_or_a_global_safety_hold():
    result = policy.classify(
        pd.DataFrame(
            [
                facts(reproductive_evidence="unassessed"),
                facts(cancer_evidence="unassessed"),
                facts(cancer_evidence="unassessed", reproductive_evidence="does_not_qualify"),
                facts(cancer_evidence="does_not_qualify", reproductive_evidence="does_not_qualify"),
            ]
        )
    )
    assert result.classification.tolist() == [
        "Cancer Antigen only",
        "Reproductive Antigen only",
        "Other:Review",
        "Other:Expression",
    ]


@pytest.mark.parametrize("bad", [None, float("nan"), "False", 0, 1])
def test_unknown_or_coerced_boolean_cannot_clear_an_eligibility_check(bad):
    with pytest.raises(ValueError, match="explicit booleans"):
        policy.classify(pd.DataFrame([facts(normal_screen_pass=bad)]))


@pytest.mark.parametrize(
    "normal,values,status,routes",
    [
        (0, [0.02, 0.03], "does_not_qualify", []),
        (0, [1, 1], "qualifies", ["two_independent_tumors_above_normal"]),
        (0.4, [4.99, 0], "does_not_qualify", []),
        (0.4, [5, 0], "qualifies", ["one_tumor_fivefold_above_normal"]),
        (3, [3, 3], "does_not_qualify", []),
        (3, [3.01, 3.01], "qualifies", ["two_independent_tumors_above_normal"]),
        (3, [14.99, 0], "does_not_qualify", []),
        (3, [15, 0], "qualifies", ["one_tumor_fivefold_above_normal"]),
        (
            3,
            [15, 4],
            "qualifies",
            ["two_independent_tumors_above_normal", "one_tumor_fivefold_above_normal"],
        ),
    ],
)
def test_ratio_boundaries_and_low_expression_floor(normal, values, status, routes):
    result = assess(*observations(values, (normal, normal)))
    assert result["status"] == status
    assert result["routes"] == routes
    assert result["normal_maximum"] == normal


def test_replicates_and_duplicate_cohort_views_are_one_participant():
    tumors, normals = observations()
    tumors["participant_id"] = "same-person"
    tumors["source_id"] = ["original-study", "reused-study"]
    result = assess(tumors, normals)
    assert result["status"] == "does_not_qualify"
    assert result["independent_participants_above_normal"] == 1
    tumors["participant_id"] = None
    assert assess(tumors, normals)["status"] == "does_not_qualify"
    tumors.loc[0, "value"] = 5
    assert assess(tumors, normals)["routes"] == ["one_tumor_fivefold_above_normal"]


@pytest.mark.parametrize(
    "field,bad", [("value", None), ("qc_pass", False), ("identity_resolved", False)]
)
def test_missing_or_failed_required_normals_never_become_zero(field, bad):
    tumors, normals = observations((100, 200))
    normals.loc[0, field] = bad
    result = assess(tumors, normals)
    assert result["status"] == "unassessed"
    assert result["missing_normal_contexts"] == ["context-0"]


def test_compartment_allowances_exclude_only_the_named_normal_compartment():
    tumors, normals = observations()
    extra = pd.DataFrame(
        [
            dict(normals.iloc[0], context_id="hair", scope="disposable", value=100),
            dict(normals.iloc[0], context_id="germ-cells", scope="reproductive", value=1000),
        ]
    )
    assert assess(tumors, pd.concat([normals, extra]))["status"] == "qualifies"
    # Immune cells in a reproductive organ retain their normal-cell scope.
    extra.loc[1, ["context_id", "scope", "value"]] = ["testis-T-cells", "normal", 4]
    result = assess(tumors, pd.concat([normals, extra]))
    assert result["status"] == "does_not_qualify"
    assert result["normal_maximum"] == 4


@pytest.mark.parametrize("field,bad", [("unit", "nTPM"), ("comparison_id", "other-panel")])
def test_incompatible_scales_cannot_supply_a_fold_change(field, bad):
    tumors, normals = observations()
    normals.loc[0, field] = bad
    with pytest.raises(ValueError, match="different units or comparison panels"):
        assess(tumors, normals)


def test_rejected_tumor_qc_and_conflicting_sample_identity():
    tumors, normals = observations((100, 2))
    tumors.loc[0, "qc_pass"] = False
    assert assess(tumors, normals)["status"] == "does_not_qualify"
    tumors.loc[1, "sample_id"] = tumors.loc[0, "sample_id"]
    with pytest.raises(ValueError, match="Conflicting participant_id"):
        assess(tumors, normals)


def test_generated_definition_is_the_runtime_contract():
    path = Path(__file__).parents[1] / "docs/audits/target-selection-20261109/definition.json"
    assert json.loads(path.read_text()) == policy.definition()


def test_report_embeds_both_figures_without_conflicting_svg_ids():
    class Parser(HTMLParser):
        def __init__(self):
            super().__init__()
            self.ids = []
            self.figures = 0

        def handle_starttag(self, tag, attrs):
            if tag == "svg":
                self.figures += 1
            if "id" in dict(attrs):
                self.ids.append(dict(attrs)["id"])

    parser = Parser()
    path = Path(__file__).parents[1] / "docs/audits/target-selection-20261109/report.html"
    parser.feed(path.read_text())
    assert parser.figures == 2
    assert len(parser.ids) == len(set(parser.ids))

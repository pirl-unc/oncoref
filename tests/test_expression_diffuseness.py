"""Concentration is a profile property, not a safety or expression-presence gate."""

import numpy as np
import pandas as pd
import pytest

from oncoref.expression_diffuseness import assess, hpa_bulk_panel


def example(values, reference=10, weights=None):
    meta = {
        "entity_id": "p",
        "source": "assay",
        "version": "v1",
        "taxonomy_id": 9606,
        "unit": "nTPM",
        "cohort": "cohort",
        "summary": "summary",
    }
    tissues = [str(i) for i in range(len(values))]
    d = pd.DataFrame([meta | {"compartment": t, "value": n} for t, n in zip(tissues, values)])
    r = pd.DataFrame([meta | {"reference_value": reference, "reference_complete": True}])
    p = pd.DataFrame({"compartment": tissues, "weight": weights or [1 / len(values)] * len(values)})
    return d, r, p


def test_focal_outlier_is_distinct_from_diffuse_expression():
    _, focal = assess(*example([1000] + [0] * 29))
    _, diffuse = assess(*example([10] * 30))
    assert focal.diffuseness.iloc[0] < 1 / 30
    assert diffuse.diffuseness.iloc[0] == pytest.approx(0.5)
    assert focal.review.iloc[0] == "concentrated"
    assert diffuse.review.iloc[0] == "diffuse"


def test_low_expression_is_continuous_and_scale_invariant():
    _, low = assess(*example([1] * 30))
    _, same_ratio = assess(*example([100] * 30, reference=1000))
    assert low.diffuseness.iloc[0] == pytest.approx(1 / 11)
    assert same_ratio.diffuseness.iloc[0] == pytest.approx(low.diffuseness.iloc[0])
    _, tiny = assess(*example([0.2] * 30, reference=0.1))
    assert tiny.review.iloc[0] == "weak_reference"


def test_missing_data_are_bounds_not_zeros_or_renormalized_measurements():
    _, missing = assess(*example([10, np.nan]))
    assert missing.lower.iloc[0] == pytest.approx(0.25)
    assert missing.upper.iloc[0] == pytest.approx(0.75)
    assert pd.isna(missing.diffuseness.iloc[0])
    _, absent = assess(*example([np.nan, np.nan]))
    assert absent.lower.iloc[0] == 0 and absent.upper.iloc[0] == 1
    assert absent.review.iloc[0] == "incomplete_normal_coverage"


def test_reference_source_and_coverage_cannot_be_borrowed():
    d, r, p = example([10, 10])
    r["source"] = "different_assay"
    _, unmatched = assess(d, r, p)
    assert unmatched.review.iloc[0] == "no_reference_signal"
    d, r, p = example([10, 10])
    r["reference_complete"] = False
    _, partial = assess(d, r, p)
    assert partial.lower.iloc[0] == 0 and partial.upper.iloc[0] == 0.5
    assert pd.isna(partial.diffuseness.iloc[0])
    assert partial.review.iloc[0] == "incomplete_reference"


def test_tissue_subdivision_preserves_weight_and_disposable_counts():
    _, original = assess(*example([5, 0], weights=[0.5, 0.5]))
    _, split = assess(*example([5, 5, 0], weights=[0.25, 0.25, 0.5]))
    assert original.diffuseness.iloc[0] == pytest.approx(split.diffuseness.iloc[0])
    panel = hpa_bulk_panel()
    assert len(panel) == 38 and panel.anatomy.nunique() == 29
    assert {"appendix", "gallbladder", "tonsil"} <= set(panel.compartment)
    assert panel[panel.anatomy.eq("brain")].weight.sum() == pytest.approx(1 / 29)


def test_zero_reference_and_conflicting_rows_do_not_fake_specificity():
    _, empty = assess(*example([0, 0], reference=0))
    assert empty.review.iloc[0] == "no_reference_signal"
    assert pd.isna(empty.diffuseness.iloc[0])
    d, r, p = example([1, 2])
    with pytest.raises(ValueError, match="Conflicting duplicate"):
        assess(pd.concat([d, d.iloc[[0]].assign(value=3)]), r, p)


def test_completely_absent_protein_remains_in_coverage_report():
    d, r, p = example([1, 2])
    _, summary = assess(d.iloc[:0], r, p)
    assert summary.entity_id.tolist() == ["p"]
    assert summary.normal_coverage.item() == 0
    assert summary.lower.item() == 0 and summary.upper.item() == 1
    assert pd.isna(summary.diffuseness.item())


def test_thymus_is_permissive_not_disposable_or_a_reproductive_reference():
    from oncoref.expression_diffuseness import hpa_bulk_roles, profile_role
    from oncoref.normal_tissues import normal_tissue_category, resolve_normal_tissue

    for species in ("human", "dog", "mouse"):
        thymus = resolve_normal_tissue("thymus", species=species)
        assert profile_role(thymus) == "permissive"
        assert normal_tissue_category(thymus) == "manageable"
        assert (
            profile_role(resolve_normal_tissue("thymus", species=species, stage="fetal"))
            == "reproductive"
        )
        for tissue in ("tonsil", "spleen", "lymph node"):
            annotated = resolve_normal_tissue(tissue, species=species)
            if annotated.mapping_status == "mapped":
                assert profile_role(annotated) == "normal"
    roles = hpa_bulk_roles().set_index("compartment")
    assert len(roles) == 50 and roles.index.is_unique
    assert roles.loc["thymus", "role"] == "permissive"
    assert roles.loc["thymus", "weight"] == 0
    assert roles[roles.role.eq("normal")].weight.sum() == pytest.approx(1)
    assert roles[roles.role.eq("reproductive")].shape[0] == 11
    assert "thymus" not in set(hpa_bulk_panel().compartment)

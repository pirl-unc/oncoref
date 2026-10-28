import numpy as np
import pandas as pd
import pytest

from oncoref.normal_expression_breadth import anatomical_panel, burden, donor_summaries


def profile(values=(1.0, 3.0)):
    return pd.DataFrame(
        {
            "gene_id": "g",
            "source": "hpa_rna_consensus",
            "version": "v23",
            "taxonomy_id": 9606,
            "unit": "nTPM",
            "cohort": "consensus",
            "summary": "consensus",
            "compartment": ["a", "b"],
            "value": values,
            "ceiling": 2.0,
        }
    )


def panel():
    return pd.DataFrame({"compartment": ["a", "b"], "category": "vital", "weight": [0.4, 0.6]})


def test_known_burden_missing_range_and_no_favorable_renormalization():
    _, result = burden(profile(), panel())
    assert result.lower.iloc[0] == pytest.approx(0.4 / 3 + 0.6 * 3 / 5)
    _, missing = burden(profile((1.0, np.nan)), panel())
    assert missing.lower.iloc[0] == pytest.approx(0.4 / 3)
    assert missing.upper.iloc[0] == pytest.approx(0.4 / 3 + 0.6)
    _, absent = burden(profile((np.nan, np.nan)), panel())
    assert absent.lower.iloc[0] == 0 and absent.upper.iloc[0] == 1


def test_monotonic_and_reproductive_values_cannot_rescue():
    d = profile()
    _, before = burden(d, panel())
    d["reproductive_rna"] = 1e9
    _, unrelated = burden(d, panel())
    pd.testing.assert_frame_equal(before, unrelated)
    d.loc[0, "value"] = 20
    _, after = burden(d, panel())
    assert after.lower.iloc[0] > before.lower.iloc[0]


def test_duplicates_do_not_change_weights_or_observation_counts():
    d = profile()
    _, result = burden(d, panel())
    _, copied = burden(pd.concat([d, d]), panel())
    pd.testing.assert_frame_equal(result, copied)
    wrong = d.copy()
    wrong.loc[0, "value"] = 5
    with pytest.raises(ValueError, match="Conflicting duplicate"):
        burden(pd.concat([d, wrong]), panel())


def test_compartment_splitting_conserves_weight():
    d, p = profile(), panel()
    _, before = burden(d, p)
    extra = d.iloc[[0]].copy()
    extra.compartment = "a2"
    p.loc[0, "weight"] = 0.2
    p.loc[2] = ["a2", "vital", 0.2]
    _, after = burden(pd.concat([d, extra]), p)
    assert after.lower.iloc[0] == pytest.approx(before.lower.iloc[0])


@pytest.mark.parametrize(
    "field,value",
    [
        ("unit", "TPM"),
        ("unit", "CP10K"),
        ("taxonomy_id", 10090),
        ("source", "gtex"),
        ("version", "v24"),
        ("value", -1),
        ("ceiling", 0),
    ],
)
def test_invalid_units_species_sources_and_measurements(field, value):
    d = profile()
    d[field] = value
    with pytest.raises(ValueError):
        burden(d, panel())


def test_fixed_panel_brain_has_one_organ_weight_and_explicit_gaps():
    p = anatomical_panel()
    assert p.groupby("category").weight.sum().to_dict() == pytest.approx(
        {"vital": 1, "protected": 1, "manageable": 1}
    )
    assert p.loc[p.organ.eq("brain"), "weight"].sum() == pytest.approx(
        p.loc[p.organ.eq("kidney"), "weight"].sum()
    )
    assert {"inner ear", "non-retinal eye", "cartilage"} <= set(p.compartment)
    assert "testis" not in set(p.compartment)


def test_donor_counts_are_independent_and_sparse_upper_quantiles_unavailable():
    d = pd.concat([profile().iloc[[0]]] * 3, ignore_index=True).drop(columns=["summary", "ceiling"])
    d["donor_id"] = ["d1", "d2", "d3"]
    d["value"] = [1, 2, 3]
    summaries = donor_summaries(pd.concat([d, d])).set_index("summary")
    assert summaries.loc["median", "value"] == 2
    assert summaries.loc["median", "donors"] == 3
    assert pd.isna(summaries.loc["q90", "value"])
    d.loc[1, "donor_id"] = "d1"
    with pytest.raises(ValueError, match="technical replicates"):
        donor_summaries(d)


def test_cohorts_and_summary_types_are_separate():
    d = profile()
    other = d.copy()
    other.cohort = "independent"
    other.value = 20
    _, scores = burden(pd.concat([d, other]), panel())
    assert len(scores) == 2
    assert scores.lower.nunique() == 2

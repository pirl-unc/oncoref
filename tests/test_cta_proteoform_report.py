"""Scientific boundary and display semantics for the proteoform CTA report."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

from oncoref.proteoforms import collapse_to_proteoforms


def load_script(name):
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parents[1] / "scripts" / f"{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


report = load_script("cta_proteoform_report")
plots = load_script("plot_cta_proteoform_report")


def test_nyeso_sum_and_transcriptome_threshold_are_recomputed():
    genes = pd.DataFrame(
        {
            "Ensembl_Gene_ID": ["A", "B", "ENSG00000268651", "ENSG00000184033", "C", "D"],
            "Symbol": ["A", "B", "CTAG1A", "CTAG1B", "C", "D"],
            "patient": [0.0, 7.0, 4.0, 4.0, 10.0, 1.0],
        }
    )
    collapsed = collapse_to_proteoforms(genes, scope="genome", sample_cols=["patient"])
    nyeso = collapsed[collapsed.proteoform_key.eq("NY-ESO-1")][["patient"]].to_numpy()
    assert len(collapsed) == 5 and nyeso.item() == 8
    assert collapsed.patient.sum() == genes.patient.sum()
    frames, cutoff = report.panel_metrics(
        collapsed[["patient"]].to_numpy(), nyeso, np.array([True]), (70,)
    )
    assert np.isclose(cutoff.item(), 7.8)
    assert not np.isclose(cutoff.item(), np.quantile(genes.patient, 0.7))
    assert frames[0].n_expressing.item() == 1
    assert frames[0].mean_expression_expressing.item() == 8


def test_partial_member_sum_and_exact_prevalence_boundaries():
    frame = pd.DataFrame(
        {
            "complete_measurement": [True, True, False],
            "available_sum_complete_measurement": [True] * 3,
            "n_patients": [20] * 3,
            "n_expressing": [10, 11, 20],
        }
    )
    assert report.passing(frame, 50).tolist() == [False, True, False]
    assert report.passing(frame, 50, allow_partial=True).tolist() == [False, True, True]


def test_coverage_has_fixed_denominator_including_unavailable_cohorts():
    cohorts = pd.DataFrame({"cancer_code": ["P", "S", "B", "C"]})
    groups = cohorts.assign(
        overlap_group=["P", "P", "B", "C"], cancer_type_group=["P", "P", "B", "C"]
    )
    data = cohorts.assign(
        proteoform_key="K",
        n_patients=20,
        n_expressing=[2, 3, 3, 0],
        complete_measurement=[True, True, True, False],
        available_sum_complete_measurement=[True, True, True, False],
    )
    row = report.breadth_rows(data, cohorts, groups).iloc[0]
    assert row.n_cohort_views_gt10 == 2
    assert row.fraction_cohort_views_gt10 == 0.5
    assert row.fraction_cancer_type_groups_gt10 == 2 / 3
    assert row.fraction_measured_cohort_views_gt10 == 2 / 3


def test_molecular_subtypes_merge_but_distinct_histologies_remain():
    codes = ["LAML", "LAML_RISK", "SARC_A", "SARC_B"]
    groups = pd.DataFrame({"cancer_code": codes, "overlap_group": codes})
    registry = pd.DataFrame(
        {
            "code": codes,
            "parent_code": ["", "LAML", "SARC", "SARC"],
            "ontology_level": ["type", "molecular_subtype", "type", "type"],
        }
    )
    cohorts = pd.DataFrame({"cancer_code": codes, "n_patients": [100, 150, 60, 70]})
    result = report.cancer_type_groups(groups, registry, cohorts).set_index("cancer_code")
    assert result.loc["LAML", "cancer_type_group"] == result.loc["LAML_RISK", "cancer_type_group"]
    assert result.loc["SARC_A", "cancer_type_group"] != result.loc["SARC_B", "cancer_type_group"]


def test_heatmaps_keep_subthreshold_nonzero_and_zero_in_included_columns():
    data = pd.DataFrame(
        {
            "proteoform_key": ["K", "L"] * 3,
            "cancer_code": ["A", "A", "B", "B", "C", "C"],
            "n_patients": [100] * 6,
            "n_expressing": [80, 1, 11, 0, 10, 0],
            "available_sum_complete_measurement": [True] * 6,
        }
    )
    shown = plots.heatmap_cells(data, ["K", "L"], ["A", "B", "C"])
    assert set(shown.cancer_code) == {"A", "B"}
    assert shown[shown.proteoform_key.eq("L")].n_expressing.tolist() == [1, 0]

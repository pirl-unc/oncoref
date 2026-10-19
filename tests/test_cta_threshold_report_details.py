"""Cohort-size and overlap semantics for the CTA analysis extension."""

import importlib.util
from pathlib import Path

import pandas as pd

spec = importlib.util.spec_from_file_location(
    "cta_details", Path(__file__).parents[1] / "scripts/cta_threshold_report_details.py"
)
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def test_minimum_size_is_inclusive_and_only_applies_to_passing_cohorts():
    frame = pd.DataFrame(
        {
            "complete_measurement": [True] * 4,
            "n_patients": [9, 10, 20, 1000],
            "n_expressing": [9, 8, 16, 700],
            "linear_tpm_comparable": [True, False, True, True],
        }
    )
    assert report.passes(frame, 70, 10).tolist() == [False, True, True, False]
    assert report.passes(frame, 70, 20, True).tolist() == [False, False, True, False]


def test_transitive_subtype_overlap_and_physical_namespaces():
    rows = []
    for code, source, patients in [
        ("P", "TREEHOUSE_POLYA_25_01", ["A", "B"]),
        ("S1", "TREEHOUSE_POLYA_25_01_SUBTYPE", ["A"]),
        ("S2", "TREEHOUSE_POLYA_25_01_SUBTYPE", ["B"]),
        ("OTHER", "UNRELATED_STUDY", ["A"]),
    ]:
        for patient in patients:
            rows.append(
                {
                    "cancer_code": code,
                    "source_cohort": source,
                    "patient_id": patient,
                    "patient_id_basis": "source_sample_unverified_patient",
                    "included_in_report": True,
                }
            )
    c = pd.DataFrame({"cancer_code": ["P", "S1", "S2", "OTHER"], "n_patients": [2, 1, 1, 1]})
    groups, edges = report.overlap_groups(pd.DataFrame(rows), c)
    mapping = groups.set_index("cancer_code").overlap_group
    assert mapping["P"] == mapping["S1"] == mapping["S2"]
    assert mapping["OTHER"] != mapping["P"]
    assert len(edges) == 2


def test_largest_support_is_not_highest_prevalence_and_ties_preserved():
    selected = pd.DataFrame(
        {
            report.ID: ["G"] * 3,
            "Symbol": ["GENE"] * 3,
            "cancer_code": ["SMALL", "BIG1", "BIG2"],
            "cancer_name": ["small", "large one", "large two"],
            "n_patients": [2, 100, 100],
            "n_expressing": [2, 80, 90],
            "fraction_expressing": [1, 0.8, 0.9],
            "mean_expression_expressing": [100, 20, 30],
            "expression_unit": ["biological_clean_TPM"] * 3,
            "linear_tpm_comparable": [True] * 3,
        }
    )
    groups = pd.DataFrame(
        {"cancer_code": ["SMALL", "BIG1", "BIG2"], "overlap_group": ["A", "B", "B"]}
    )
    result = report.summarize_genes(selected, groups).iloc[0]
    assert result.largest_passing_cohort_size == 100
    assert result.largest_passing_cohorts == "BIG1;BIG2"
    assert result.n_passing_cohort_views == 3
    assert result.n_overlap_adjusted_groups == 2
    assert result.passes_min_cohort_size and result.passes_large_linear_filter

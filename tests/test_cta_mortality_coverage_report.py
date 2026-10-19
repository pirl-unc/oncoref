"""Patient-level panel coverage must preserve coexpression and missing-data meaning."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location(
    "cta_mortality_report", Path(__file__).parents[1] / "scripts/cta_mortality_coverage_report.py"
)
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def test_panel_or_counts_coexpressing_patients_once_and_excludes_cutoff_ties():
    values = np.array([[11, 11, 10, 0], [11, 0, 11, 0]])
    result = report.panel_patient_coverage(values, np.full(4, 10), [True, True])
    assert result["n_expressing_any"] == 3
    assert result["fraction_expressing_any"] == 0.75
    assert not result["coverage_is_lower_bound"]


def test_missing_panel_members_produce_lower_bounds_not_negative_calls():
    values = np.array([[11, 0], [100, 100]])
    result = report.panel_patient_coverage(values, [10, 10], [True, False])
    assert result["fraction_expressing_any"] == 0.5
    assert result["fraction_expressing_any_upper_bound"] == 1
    assert result["coverage_is_lower_bound"]
    none = report.panel_patient_coverage(values, [10, 10], [False, False])
    assert np.isnan(none["fraction_expressing_any"])


def test_ranking_excludes_residual_and_keeps_named_categories():
    categories = list(report.LABELS)
    rows = pd.DataFrame(
        {
            "burden_category": ["other_and_unknown_primary", *categories],
            "world_mortality_pct": [90, *range(10, 0, -1)],
            "derivation_basis": ["residual", *(["direct_source"] * 10)],
        }
    )
    top = report.top_mortality_categories(rows)
    assert top.burden_category.tolist() == categories
    assert top.mortality_rank.tolist() == list(range(1, 11))

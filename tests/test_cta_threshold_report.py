"""Numerical boundary checks for the one-off CTA report."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location(
    "cta_report", Path(__file__).parents[1] / "scripts/cta_threshold_report.py"
)
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def test_within_patient_quantile_positive_means_and_ties():
    # Medians are 2 and 20. Values equal to the median do not express.
    values = np.array([[0, 0], [1, 10], [2, 20], [4, 80], [9, 90]], dtype=float)
    frames, cutoffs = report.threshold_metrics(values, [0, 2, 3], percentiles=(50,))
    assert np.array_equal(cutoffs, [[2, 20]])
    result = frames[0]
    assert result.n_expressing.tolist() == [0, 0, 2]
    assert result.mean_expression_expressing.iloc[:2].isna().all()
    assert result.mean_expression_expressing.iloc[2] == 42  # linear arithmetic mean


def test_exact_prevalence_boundaries_are_excluded():
    frame = pd.DataFrame(
        {"complete_measurement": [True] * 4, "n_expressing": [1, 2, 3, 4], "n_patients": [4] * 4}
    )
    assert report.qualifies(frame, 25, minimum_patients=1).tolist() == [False, True, True, True]
    assert report.qualifies(frame, 50, minimum_patients=1).tolist() == [False, False, True, True]
    assert report.qualifies(frame, 75, minimum_patients=1).tolist() == [False, False, False, True]


def test_missing_is_not_zero_and_incomplete_measurement_cannot_select():
    values = np.array([[0, 0], [9, np.nan], [1, 1]], dtype=float)
    frames, _ = report.threshold_metrics(values, [1], percentiles=(50,))
    row = frames[0].iloc[0]
    assert row.n_measured_patients == 1
    assert row.n_expressing == 1
    assert row.mean_expression_expressing == 9
    assert not report.qualifies(frames[0], 25).iloc[0]


def test_zero_expression_never_positive_and_patient_repeats_average():
    values = np.array([[0, 0, 0], [0, 0, 0], [2, 6, 10]], dtype=float)
    grouped, patients, _ = report.group_patients(
        values, ["TH01_0001_S01", "TH01_0001_S02", "TH01_0002_S01"], "TH", {}
    )
    assert patients == ["TH01_0001", "TH01_0002"]
    assert np.array_equal(grouped, [[0, 0], [0, 0], [4, 10]])
    frames, _ = report.threshold_metrics(grouped, [0, 2], percentiles=(30,))
    assert frames[0].n_expressing.tolist() == [0, 2]
    assert frames[0].n_patients.tolist() == [2, 2]
    assert report.donor_identity("THR24_2086_S02", "TH", {}) == (
        "THR24_2086",
        "Treehouse_donor_prefix",
    )

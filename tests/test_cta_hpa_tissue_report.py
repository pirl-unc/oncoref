"""Preserve distinct meanings of zero, missing, RNA sums and ordinal protein evidence."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location(
    "cta_hpa_report", Path(__file__).parents[1] / "scripts/cta_hpa_tissue_report.py"
)
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def test_rna_sum_requires_all_members_but_keeps_recorded_zero():
    data = pd.DataFrame({"testis": [4, 5], "liver": [0, 0], "brain": [3, np.nan]})
    summed = report.sum_member_rna(data)
    assert summed.testis == 9
    assert summed.liver == 0
    assert np.isnan(summed.brain)


def test_ihc_takes_strongest_cell_type_and_keeps_unavailable_distinct():
    observations = pd.DataFrame(
        {
            "Gene": ["A", "A", "A"],
            "Tissue": ["testis", "testis", "liver"],
            "Level": ["Low", "High", "Not detected"],
        }
    )
    matrix = report.ihc_matrix(observations, ["A", "B"], ["testis", "liver", "brain"])
    assert matrix.loc["A", "testis"] == 3
    assert matrix.loc["A", "liver"] == 0
    assert matrix.loc["B"].isna().all()
    assert np.isnan(matrix.loc["A", "brain"])


def test_display_group_keeps_thymus_in_other_and_breast_explicit():
    reproductive, other = report.ordered_tissues(["thymus", "testis", "breast", "liver"])
    assert reproductive == ["testis", "breast"]
    assert other == ["liver", "thymus"]


def test_rare_study_tissue_does_not_count_as_routine_panel():
    rows = [(f"G{i}", "lung", "Not detected") for i in range(10)]
    rows += [("G0", "retina", "Not detected")]
    observations = pd.DataFrame(rows, columns=["Gene", "Tissue", "Level"])
    assert report.routine_ihc_tissues(observations) == {"lung"}

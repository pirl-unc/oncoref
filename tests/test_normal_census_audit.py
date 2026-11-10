"""Adversarial checks for donor-resolved atlas inventory and count aggregation."""

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("cellxgene_census")
pa = pytest.importorskip("pyarrow")
SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location(
    "audit_normal_census", SCRIPTS / "audit_normal_census.py"
)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


@pytest.mark.parametrize(
    "stage,expected",
    [
        ("18-year-old stage", "adult"),
        ("17-year-old stage", "child"),
        ("12th week post-fertilization stage", "prenatal"),
        ("15-month-old stage", "child"),
        ("15-19 year-old", "unresolved"),
        ("postnatal stage", "unresolved"),
        ("unknown", "unresolved"),
        ("Carnegie stage 20", "prenatal"),
        ("third decade stage", "adult"),
        ("adult stage", "unresolved"),
        ("prime adult stage", "unresolved"),
        ("young adult stage", "unresolved"),
    ],
)
def test_age_is_not_inferred_from_cell_name_or_ambiguous_ranges(stage, expected):
    assert audit.age_scope(stage) == expected


class FakeQuery:
    def __init__(self, table):
        self.table = table

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def X(self, layer):
        assert layer == "raw"
        return SimpleNamespace(tables=lambda: iter([self.table]))


def run_summary(rows, values, denominator=3):
    table = pa.table({"soma_dim_0": rows, "soma_dim_1": [50] * len(rows), "soma_data": values})
    experiment = SimpleNamespace(axis_query=lambda *args, **kwargs: FakeQuery(table))
    census = {"census_data": {"homo_sapiens": experiment}}
    # Actual IDs are non-contiguous; neither cells nor features are positional.
    mapping = np.full(11, -1)
    mapping[[2, 7, 10]] = 0
    features = pd.DataFrame({"soma_joinid": [50], "feature_id": ["GENE"]})
    groups = pd.DataFrame({"n_cells": [denominator]})
    return audit.summarize_block(census, features, mapping, groups)


def test_fractional_counts_are_retained_and_flagged_not_rounded_to_umis():
    out = run_summary([2, 7, 10], [1.0, 2.0, 0.4]).iloc[0]
    assert out.n_positive == 3
    assert out.n_ge2 == 1
    assert out.n_fractional == 1
    assert out.count_sum == pytest.approx(3.4)
    assert out.count_max == 2


def test_expression_outside_pinned_inventory_fails():
    with pytest.raises(ValueError, match="pinned inventory"):
        run_summary([3], [1.0])


def test_invalid_count_or_impossible_denominator_fails():
    with pytest.raises(ValueError, match="Invalid sparse counts"):
        run_summary([2], [-1.0])
    with pytest.raises(ValueError, match="exceeds population"):
        run_summary([2, 7], [1.0, 1.0], denominator=1)

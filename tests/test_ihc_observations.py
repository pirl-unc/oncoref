"""Biological result categories remain distinct from validation and missing data."""

import pytest

from oncoref.ihc_observations import parse_staining_statement


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("Negative.", "not_detected"),
        ("Distinct cytoplasmic staining of epithelial cells.", "detected"),
        ("Weak staining (cross-reactive).", "reported_cross_reactivity"),
        ("Strong staining of epithelial cells; lymphocytes are negative.", "mixed"),
        ("A positive control should show strong staining.", "not_assessed"),
        ("", "not_assessed"),
    ],
)
def test_explicit_categorical_result(text, expected):
    assert parse_staining_statement(text)["result"] == expected


def test_numeric_bounds_and_original_scale_are_preserved():
    row = parse_staining_statement(
        "Weak nuclear staining in >20% of cells; 10–30% in other samples."
    )
    assert row["intensity"] == "low"
    assert row["localization"] == "nuclear"
    assert row["reported_percentages"] == [
        {"operator": ">", "low": 20.0, "high": 20.0},
        {"operator": "range", "low": 10.0, "high": 30.0},
    ]

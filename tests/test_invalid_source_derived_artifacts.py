# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Derived artifacts of an invalid-scale source are withheld by every reader (#604).

The legacy UNC NUTM1 source is registered ``invalid_legacy_tpm``. ``per_sample_expression``
already rejects it; every artifact derived from it (percentiles, representatives,
within-sample fractions, cohort reference summaries) must be withheld too, and a named
cohort must raise rather than come back empty.
"""

import warnings

import pandas as pd
import pytest

import oncoref
from oncoref import expression

INVALID = "NUTM"


@pytest.fixture(autouse=True)
def fresh_caches():
    expression._invalid_selected_source_codes.cache_clear()
    yield
    expression._invalid_selected_source_codes.cache_clear()


def test_manifest_marks_the_nutm_selected_source_invalid():
    invalid = expression._invalid_selected_source_codes()
    assert INVALID in invalid
    assert "UNC_NUTM1" in invalid[INVALID]
    assert "oncoref#604" in invalid[INVALID]
    table = expression._reference_summary_source_table()
    flagged = table[
        table.selected
        & table.source_scale_class.astype("string").isin(expression.INVALID_SOURCE_SCALE_CLASSES)
    ]
    # A source the registry calls invalid is never also called linear-TPM comparable.
    assert not flagged.linear_tpm_comparable.astype(bool).any()
    assert set(flagged.cancer_code.astype(str)) == set(invalid)


READERS = {
    "percentiles": lambda: oncoref.cohort_gene_percentiles(INVALID, auto_fetch=False),
    "within-sample fractions": lambda: oncoref.gene_within_sample_top_fraction(INVALID),
    "proteoform within-sample fractions": lambda: oncoref.proteoform_within_sample_top_fraction(
        INVALID, auto_fetch=False
    ),
    "gene representatives": lambda: oncoref.gene_representative_samples(INVALID),
    "proteoform representatives": lambda: oncoref.proteoform_representative_samples(INVALID),
    "cohort reference summary": lambda: oncoref.cancer_reference_expression(
        INVALID, auto_fetch=False
    ),
    "cohort reference summary, on_missing=empty": lambda: oncoref.cancer_reference_expression(
        INVALID, auto_fetch=False, on_missing="empty"
    ),
    "selected summary rows": lambda: oncoref.cancer_reference_expression(
        INVALID, auto_fetch=False, reference_source="summary_rows", sample_qc="all"
    ),
}


@pytest.mark.parametrize("reader", sorted(READERS))
def test_every_derived_reader_rejects_a_named_invalid_cohort(reader):
    with pytest.raises(ValueError, match="invalid expression scale"):
        READERS[reader]()


def test_named_cohort_is_never_silently_dropped_from_a_multi_cohort_request():
    with pytest.raises(ValueError, match=r"NUTM.*invalid expression scale"):
        oncoref.cancer_reference_expression(["GBM", INVALID], auto_fetch=False)


def test_cohort_listings_omit_invalid_cohorts(tmp_path, monkeypatch):
    for code in ("GBM", INVALID):
        (tmp_path / f"{code}.parquet").write_bytes(b"")
    monkeypatch.setattr(expression, "_local_shard_dirs", lambda *a, **k: (tmp_path,))
    assert expression.available_percentile_cohorts() == ["GBM"]
    assert expression.available_representative_cohorts() == ["GBM"]
    assert expression.available_within_sample_cohorts() == ["GBM"]


def test_default_request_listing_omits_invalid_cohorts_with_a_warning(monkeypatch):
    monkeypatch.setattr(
        expression, "_reference_available_codes_for_modes", lambda *a, **k: ["GBM", INVALID]
    )
    kwargs = {"reference_source": "artifact", "sample_qc": "pass"}
    with pytest.warns(UserWarning, match="Omitting cohort"):
        rows = expression._reference_expression_requests(None, ["tpm_clean"], **kwargs)
    assert [row["cancer_code"] for row in rows] == ["GBM"]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        described = expression._reference_expression_requests(
            None, ["tpm_clean"], describe_invalid=True, **kwargs
        )
    assert [row["cancer_code"] for row in described] == ["GBM", INVALID]


def test_availability_still_describes_the_invalid_source_without_raising():
    table = oncoref.cancer_reference_expression_availability(INVALID)
    row = table.iloc[0]
    assert row.source_scale_class == "invalid_legacy_tpm"
    assert not bool(row.linear_tpm_comparable)
    assert not bool(row.available)
    # The reason says the artifact is withheld, not that it is missing.
    assert row.missing_reason == "invalid_source_scale"
    all_sources = oncoref.cancer_reference_expression_availability(
        INVALID, reference_source="summary_rows_all", sample_qc="all", all_sources=True
    )
    assert set(all_sources.source_cohort) >= {"UNC_NUTM1"}


def test_valid_cohorts_and_codes_are_not_affected():
    expression._reject_invalid_source_code("GBM", context="test")  # does not raise
    assert "GBM" not in expression._invalid_selected_source_codes()


@pytest.fixture
def two_source_table(monkeypatch):
    columns = list(expression._reference_summary_source_table().columns)
    rows = [
        {"cancer_code": INVALID, "source_cohort": "UNC_NUTM1", "selected": True,
         "source_scale_class": "invalid_legacy_tpm", "linear_tpm_comparable": False,
         "source_type": "unc-case-series"},
        {"cancer_code": INVALID, "source_cohort": "TREEHOUSE_POLYA_25_01", "selected": False,
         "source_scale_class": "linear_rnaseq_tpm", "linear_tpm_comparable": True,
         "source_type": "treehouse-compendium"},
    ]  # fmt: skip
    table = pd.DataFrame(rows).reindex(columns=columns)
    monkeypatch.setattr(expression, "_reference_summary_source_table", lambda: table)
    seen = {}

    def record(code, sources, **kwargs):
        seen["sources"] = sorted(sources.source_cohort.astype(str))
        return pd.DataFrame()

    monkeypatch.setattr(expression, "_reference_summary_rows_for_sources", record)
    expression._invalid_selected_source_codes.cache_clear()
    return seen


def test_source_union_omits_the_invalid_source_unless_explicitly_requested(two_source_table):
    with pytest.warns(UserWarning, match=r"Omitting NUTM source\(s\) \['UNC_NUTM1'\]"):
        expression._reference_summary_all_rows(INVALID)
    assert two_source_table["sources"] == ["TREEHOUSE_POLYA_25_01"]
    expression._reference_summary_all_rows(INVALID, source_cohorts={"TREEHOUSE_POLYA_25_01"})
    assert two_source_table["sources"] == ["TREEHOUSE_POLYA_25_01"]
    # Naming the invalid source passes it through to the loader, which raises on it.
    expression._reference_summary_all_rows(INVALID, source_cohorts={"UNC_NUTM1"})
    assert "UNC_NUTM1" in two_source_table["sources"]


def test_loader_rejects_any_invalid_source_it_is_handed():
    sources = pd.DataFrame(
        {"source_cohort": ["UNC_NUTM1"], "source_scale_class": ["invalid_legacy_tpm"]}
    )
    with pytest.raises(ValueError, match="invalid expression scale"):
        expression._reference_summary_rows_for_sources(INVALID, sources)

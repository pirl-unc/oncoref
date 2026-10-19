"""Scientific boundaries in the reproductive/disposable report's filtering."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

spec = importlib.util.spec_from_file_location(
    "cta_scope_report", Path(__file__).parents[1] / "scripts/cta_scope_report.py"
)
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)

GID = "ENSG00000999999"


def audit(rna_values, ihc_values=(), *, missing_rna=(), unreliable=False):
    rna = pd.DataFrame(
        [(GID, "EXAMPLE", tissue, value) for tissue, value in rna_values],
        columns=["Gene", "Gene name", "Tissue", "nTPM"],
    )
    ihc = pd.DataFrame(
        [(GID, "EXAMPLE", tissue, cell, level, "Enhanced") for tissue, cell, level in ihc_values],
        columns=["Gene", "Gene name", "Tissue", "Cell type", "Level", "Reliability"],
    )
    has_protein = (
        any(row[2] in report.PROTEIN_DETECTED_LEVELS for row in ihc_values) and not unreliable
    )
    table = pd.DataFrame(
        [
            {
                "Ensembl_Gene_ID": GID,
                "Symbol": "EXAMPLE",
                "biotype": "protein_coding",
                "source_databases": "Wang2016_CT",
                "protein_reliability": "Enhanced" if has_protein else "no data",
                "protein_reproductive": "False",
                "rna_reproductive_tissue_scope": "extended",
                "rna_deflated_reproductive_frac": 0.0,
                "passes_filters": False,
                "never_expressed": False,
                "restriction_confidence": "MODERATE",
                "rna_restriction_level": "STRICT",
                "rna_restriction": "TESTIS",
                "rna_testis_ntpm": 10.0,
                "rna_ovary_ntpm": 0.0,
                "rna_placenta_ntpm": 0.0,
            }
        ]
    )
    observations = {
        name: report.annotate_measurements(frame, name, {GID})
        for name, frame in (("bulk_rna", rna), ("ihc", ihc))
    }
    expected = {
        "bulk_rna": {(tissue, "") for tissue, _ in rna_values}
        | {(tissue, "") for tissue in missing_rna},
        "ihc": {(tissue, cell) for tissue, cell, _ in ihc_values},
    }
    return report.compare_candidates(
        table,
        observations,
        unreliable_ids={GID} if unreliable else set(),
        expected_contexts=expected,
    ).set_index("scope")


def test_cell_resolved_disposable_gain_does_not_admit_whole_skin_rna():
    cells = [("skin 1", "melanocytes", "High")]
    result = audit([("testis", 100), ("skin", 0)], cells)
    assert not result.loc["reproductive", "adaptive_retained"]
    assert result.loc["reproductive_or_disposable", "adaptive_retained"]
    assert not result.loc["reproductive", "strict_retained"]
    assert result.loc["reproductive_or_disposable", "strict_retained"]
    skin = audit([("testis", 1), ("skin", 100)], cells)
    assert not skin.adaptive_retained.any()
    assert not skin.strict_retained.any()
    assert skin.outside_rna_tissues.eq("skin").all()


def test_trace_outside_rna_separates_adaptive_and_strict_policies():
    result = audit([("testis", 100), ("heart muscle", 0.1)])
    assert result.adaptive_retained.all()
    assert not result.strict_retained.any()
    assert result.outside_rna_max_ntpm.eq(0.1).all()


def test_missing_rna_is_not_negative_coverage():
    result = audit([("testis", 100)], missing_rna=["heart muscle"])
    assert result.strict_status.eq("unresolved").all()
    assert result.rna_routine_missing.eq(1).all()
    assert not result.strict_retained.any()


def test_missing_ihc_and_measured_negative_ihc_remain_distinct():
    missing = audit([("testis", 100), ("heart muscle", 0)])
    negative = audit(
        [("testis", 100), ("heart muscle", 0)], [("heart muscle", "myocytes", "Not detected")]
    )
    assert missing.strict_retained.all()
    assert negative.strict_retained.all()
    assert missing.ihc_status.eq("unavailable").all()
    assert negative.ihc_status.eq("no_positive_observation").all()
    assert missing.ihc_coverage.eq("unavailable").all()
    assert negative.ihc_coverage.eq("routine_complete").all()
    assert missing.absence_evidence.eq("not_established").all()


def test_zero_only_gene_never_passes_strict_vacuously():
    result = audit([("testis", 0), ("heart muscle", 0)])
    assert result.strict_status.eq("no_positive_observation").all()
    assert not result.strict_retained.any()


def test_entirely_missing_rna_keeps_unknown_maximum_and_fails_both_policies():
    result = audit([], missing_rna=["testis", "heart muscle"])
    assert result.rna_status.eq("unavailable").all()
    assert result.outside_rna_max_ntpm.isna().all()
    assert result.strict_status.eq("unresolved").all()
    assert not result.strict_retained.any()
    assert not result.adaptive_retained.any()


@pytest.mark.parametrize("level", ["Ascending", "Descending"])
def test_gradient_ihc_blocks_strict_even_when_legacy_adaptive_ignores_it(level):
    result = audit([("testis", 100), ("heart muscle", 0)], [("heart muscle", "myocytes", level)])
    assert result.adaptive_retained.all()
    assert result.strict_status.eq("outside_scope").all()
    assert not result.strict_retained.any()


def test_reviewed_unreliable_ihc_is_not_silently_absent_in_strict_screen():
    result = audit(
        [("testis", 100)], [("testis", "cells in seminiferous ducts", "High")], unreliable=True
    )
    assert result.adaptive_retained.all()
    assert result.strict_status.eq("unresolved").all()
    assert not result.strict_retained.any()


def test_unmapped_positive_and_invalid_measurement_block_claim():
    result = audit([("testis", 100), ("unknown organ", 2)])
    assert result.rna_status.eq("unresolved").all()
    assert not result.strict_retained.any()
    result = audit([("testis", 100), ("heart muscle", np.nan)])
    assert result.rna_routine_missing.eq(1).all()
    assert not result.strict_retained.any()


def test_mixed_reproductive_organs_stay_inside_the_union():
    result = audit([("testis", 100), ("ovary", 10)])
    assert result.adaptive_retained.all()
    assert result.strict_retained.all()


def test_source_identity_retains_namespaces_but_collapses_mapped_aliases():
    base = {
        "Ensembl_Gene_ID": "",
        "source_gene_id": "123",
        "source_symbol": "",
        "source_id_type": "entrez",
    }
    assert report.source_identity(base) != report.source_identity(
        {**base, "source_id_type": "symbol"}
    )
    assert report.source_identity({**base, "Ensembl_Gene_ID": GID}) == GID


def test_invalid_or_duplicate_rna_cannot_be_interpreted_as_absence():
    columns = ["Gene", "Gene name", "Tissue", "nTPM"]
    for values, match in (
        ([(GID, "EXAMPLE", "testis", -1)], "Invalid"),
        ([(GID, "EXAMPLE", "testis", 1)] * 2, "Duplicate"),
    ):
        with pytest.raises(ValueError, match=match):
            report.annotate_measurements(pd.DataFrame(values, columns=columns), "bulk_rna", {GID})

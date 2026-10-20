"""Boundaries of the independent RNA and IHC gates in the permitted-tissue report."""

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPTS = Path(__file__).parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location(
    "cta_permitted_report", SCRIPTS / "cta_permitted_report.py"
)
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)
GID = "ENSG00000999999"


def evaluate(
    rna,
    ihc=(),
    *,
    vital_limit=1,
    protected_limit=5,
    scope=report.EXPANDED,
    missing=(),
    unreliable=False,
):
    rows = []
    for tissue, value in rna:
        rows.append(
            {
                "Ensembl_Gene_ID": GID,
                "Symbol": "EXAMPLE",
                "modality": "bulk_rna",
                "tissue": tissue,
                "cell_type": "",
                "rna_ntpm": value,
                "ihc_level": "",
                "ihc_reliability": "",
            }
        )
    for tissue, cell, level in ihc:
        rows.append(
            {
                "Ensembl_Gene_ID": GID,
                "Symbol": "EXAMPLE",
                "modality": "ihc",
                "tissue": tissue,
                "cell_type": cell,
                "rna_ntpm": "",
                "ihc_level": level,
                "ihc_reliability": "Enhanced",
            }
        )
    observations = report.annotate_archive(pd.DataFrame(rows))
    expected = {
        "bulk_rna": {(tissue, "") for tissue, _ in rna} | {(tissue, "") for tissue in missing},
        "ihc": {(tissue, cell) for tissue, cell, _ in ihc},
    }
    return report.evaluate_gene(
        observations,
        scope=scope,
        vital_limit=vital_limit,
        protected_limit=protected_limit,
        expected=expected,
        unreliable=unreliable,
    )


def test_raw_rna_tier_boundaries_are_inclusive_and_per_tissue():
    assert evaluate([("testis", 100), ("heart muscle", 1), ("retina", 5)])["evidence_pass"]
    assert not evaluate([("testis", 100000), ("heart muscle", 1.1)])["evidence_pass"]
    assert not evaluate([("testis", 100000), ("retina", 5.1)])["evidence_pass"]


@pytest.mark.parametrize(
    "tissue,cell",
    [
        ("liver", "hepatocytes"),
        ("pancreas", "exocrine glandular cells"),
        ("retina", "pigment epithelial cells"),
        ("skeletal muscle", "myocytes"),
        ("heart muscle", "cardiomyocytes"),
        ("lung", "alveolar cells"),
        ("cerebral cortex", "neuronal cells"),
    ],
)
def test_low_ihc_is_independent_of_rna_in_both_protected_tiers(tissue, cell):
    result = evaluate([("testis", 100), (tissue, 0)], [(tissue, cell, "Low")])
    assert not result["evidence_pass"]
    assert "positive_ihc_outside_permitted" in result["failure_reasons"]
    assert "Low" in result["ihc_exclusions"]


@pytest.mark.parametrize("level", ["Medium", "High", "Ascending", "Descending"])
def test_all_positive_ihc_levels_remain_exclusions(level):
    assert not evaluate([("testis", 100), ("liver", 0)], [("liver", "hepatocytes", level)])[
        "evidence_pass"
    ]


@pytest.mark.parametrize(
    "tissue,cell",
    [
        ("skin 1", "keratinocytes"),
        ("salivary gland", "glandular cells"),
        ("thymus", "cortical cells"),
    ],
)
def test_requested_whole_organs_are_permitted(tissue, cell):
    rna_tissue = "skin" if tissue == "skin 1" else tissue
    result = evaluate([("testis", 100), (rna_tissue, 50)], [(tissue, cell, "High")])
    assert result["evidence_pass"]
    assert not evaluate(
        [("testis", 100), (rna_tissue, 50)], [(tissue, cell, "High")], scope=report.COMBINED
    )["evidence_pass"]


def test_missing_ihc_is_explicit_and_not_negative_staining():
    missing = evaluate([("testis", 100), ("liver", 0)])
    negative = evaluate([("testis", 100), ("liver", 0)], [("liver", "hepatocytes", "Not detected")])
    assert missing["evidence_pass"] and negative["evidence_pass"]
    assert missing["ihc_coverage"] == "unavailable"
    assert missing["ihc_forbidden_contexts_observed"] == 0
    assert negative["ihc_forbidden_contexts_observed"] == 1
    assert missing["absence_evidence"] == negative["absence_evidence"] == "not_established"


def test_incomplete_and_unreliable_evidence_cannot_pass():
    assert not evaluate([("testis", 100)], missing=["heart muscle"])["evidence_pass"]
    assert not evaluate([("testis", 100)], unreliable=True)["evidence_pass"]
    assert not evaluate([("testis", 100), ("head", 0.1)])["evidence_pass"]
    assert not evaluate([("testis", 100)], [("liver", "hepatocytes", "Not representative")])[
        "evidence_pass"
    ]


def test_positive_evidence_inside_permitted_scope_is_required():
    assert not evaluate([("testis", 0), ("liver", 0)])["evidence_pass"]
    assert not evaluate([("testis", 0.5), ("retina", 3)])["evidence_pass"]
    assert evaluate(
        [("testis", 0.5), ("retina", 3)], [("testis", "cells in seminiferous ducts", "Low")]
    )["evidence_pass"]

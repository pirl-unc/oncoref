"""Boundaries of the independent RNA and IHC gates in the permitted-tissue report."""

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

from oncoref.cta_normal_screen import adaptive_normal_rna_limits

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
    protected_limit=2,
    scope=report.EXPANDED,
    missing=(),
    unreliable=False,
    adaptive=False,
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
        rna_limits=(
            adaptive_normal_rna_limits(
                observations,
                routine_ihc_contexts=expected["ihc"],
                reviewed_unreliable=unreliable,
            ).rna_ceiling_ntpm
            if adaptive
            else None
        ),
    )


def test_raw_rna_tier_boundaries_are_inclusive_and_per_tissue():
    assert evaluate([("testis", 100), ("heart muscle", 1), ("thyroid gland", 5)])["evidence_pass"]
    assert not evaluate([("testis", 100000), ("heart muscle", 1.1)])["evidence_pass"]
    assert not evaluate([("testis", 100000), ("thyroid gland", 5.1)])["evidence_pass"]
    assert evaluate([("testis", 100), ("retina", 1)])["evidence_pass"]
    assert not evaluate([("testis", 100), ("retina", 1.1)])["evidence_pass"]
    assert evaluate([("testis", 100), ("salivary gland", 2)])["evidence_pass"]
    assert not evaluate([("testis", 100), ("salivary gland", 2.1)])["evidence_pass"]


def test_adaptive_missing_rna_is_unresolved_even_with_negative_ihc():
    result = evaluate(
        [],
        [("heart muscle", "cardiomyocytes", "Not detected")],
        missing=["heart muscle", "testis"],
        adaptive=True,
    )
    assert not result["evidence_pass"]
    assert "incomplete_routine_rna" in result["failure_reasons"]


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
        ("lymph node", "germinal center cells"),
        ("thyroid gland", "glandular cells"),
    ],
)
def test_low_ihc_is_independent_of_rna_in_all_screened_tiers(tissue, cell):
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
        ("spleen", "cells in white pulp"),
    ],
)
def test_former_whole_organ_allowances_are_now_screened(tissue, cell):
    rna_tissue = "skin" if tissue == "skin 1" else tissue
    result = evaluate([("testis", 100), (rna_tissue, 50)], [(tissue, cell, "High")])
    assert not result["evidence_pass"]
    assert "positive_ihc_outside_permitted" in result["failure_reasons"]
    assert "rna_above_tissue_ceiling" in result["failure_reasons"]


def test_caution_lymphoid_expression_is_screened_and_remains_flagged():
    ihc = [("spleen", "cells in white pulp", "High")]
    result = evaluate([("testis", 100), ("spleen", 50)], ihc, adaptive=True)
    assert not result["evidence_pass"]
    assert result["tissue_review_status"] == "required"
    assert "asplenia_infection_risk_review" in result["tissue_review_flags"]
    for tissue, cell in (
        ("lymph node", "germinal center cells"),
        ("bone marrow", "hematopoietic cells"),
    ):
        result = evaluate(
            [("testis", 100), ("spleen", 50), (tissue, 0)],
            [*ihc, (tissue, cell, "Low")],
            adaptive=True,
        )
        assert not result["evidence_pass"]
        assert "positive_ihc_outside_permitted" in result["failure_reasons"]
    assert not evaluate([("testis", 100), ("spleen", 50), ("lymph node", 5.1)])["evidence_pass"]


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
    assert not evaluate([("testis", 0.5), ("thyroid gland", 3)])["evidence_pass"]
    assert evaluate(
        [("testis", 0.5), ("thyroid gland", 3)], [("testis", "cells in seminiferous ducts", "Low")]
    )["evidence_pass"]


def test_adaptive_absolute_boundary_and_positive_veto_remain_independent():
    negative = [("heart muscle", "cardiomyocytes", "Not detected")]
    assert evaluate([("testis", 100), ("heart muscle", 6)], negative, adaptive=True)[
        "evidence_pass"
    ]
    assert not evaluate([("testis", 100), ("heart muscle", 6.1)], negative, adaptive=True)[
        "evidence_pass"
    ]
    result = evaluate(
        [("testis", 100), ("heart muscle", 5), ("liver", 0)],
        [*negative, ("liver", "hepatocytes", "Low")],
        adaptive=True,
    )
    assert not result["evidence_pass"]
    assert "positive_ihc_outside_permitted" in result["failure_reasons"]


@pytest.mark.parametrize("tissue,ceiling", [("salivary gland", 12), ("thyroid gland", 30)])
def test_adaptive_protected_and_caution_ceilings_are_inclusive(tissue, ceiling):
    negative = [(tissue, "glandular cells", "Not detected")]
    assert evaluate([("testis", 100), (tissue, ceiling)], negative, adaptive=True)["evidence_pass"]
    assert not evaluate([("testis", 100), (tissue, ceiling + 0.1)], negative, adaptive=True)[
        "evidence_pass"
    ]


@pytest.mark.parametrize(
    "tissue,cell",
    [
        ("tonsil", "germinal center cells"),
        ("appendix", "glandular cells"),
        ("gallbladder", "glandular cells"),
    ],
)
def test_removable_organ_allowances_do_not_extend_to_other_immune_or_glandular_sites(tissue, cell):
    assert evaluate([("testis", 100), (tissue, 50)], [(tissue, cell, "High")])["evidence_pass"]
    assert not evaluate(
        [("testis", 100), (tissue, 50), ("salivary gland", 0)],
        [(tissue, cell, "High"), ("salivary gland", "glandular cells", "Low")],
    )["evidence_pass"]


def test_cutaneous_melanocyte_allowance_does_not_clear_other_skin_cells():
    assert evaluate([("testis", 100), ("skin", 5)], [("skin 1", "melanocytes", "High")])[
        "evidence_pass"
    ]
    assert not evaluate(
        [("testis", 100), ("skin", 5)],
        [("skin 1", "melanocytes", "High"), ("skin 1", "keratinocytes", "Low")],
    )["evidence_pass"]

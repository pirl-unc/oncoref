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
    vital_limit=2,
    protected_limit=6,
    scope=report.EXPANDED,
    missing=(),
    unreliable=False,
    adaptive=False,
    selectivity=False,
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
    audit = (
        adaptive_normal_rna_limits(
            observations, routine_ihc_contexts=expected["ihc"], reviewed_unreliable=unreliable
        )
        if adaptive or selectivity
        else None
    )
    return report.evaluate_gene(
        observations,
        scope=scope,
        vital_limit=vital_limit,
        protected_limit=protected_limit,
        expected=expected,
        unreliable=unreliable,
        rna_limits=(audit.rna_ceiling_ntpm if adaptive else None),
        selectivity_audit=(audit if selectivity else None),
    )


def test_raw_rna_tier_boundaries_are_inclusive_and_per_tissue():
    assert evaluate([("testis", 100), ("heart muscle", 2), ("thyroid gland", 10)])["evidence_pass"]
    assert not evaluate([("testis", 100000), ("heart muscle", 2.1)])["evidence_pass"]
    assert not evaluate([("testis", 100000), ("thyroid gland", 10.1)])["evidence_pass"]
    assert evaluate([("testis", 100), ("retina", 2)])["evidence_pass"]
    assert not evaluate([("testis", 100), ("retina", 2.1)])["evidence_pass"]
    assert evaluate([("testis", 100), ("salivary gland", 6)])["evidence_pass"]
    assert not evaluate([("testis", 100), ("salivary gland", 6.1)])["evidence_pass"]


def test_adaptive_missing_rna_is_unresolved_even_with_negative_ihc():
    result = evaluate(
        [],
        [("heart muscle", "cardiomyocytes", "Not detected")],
        missing=["heart muscle", "testis"],
        adaptive=True,
    )
    assert not result["evidence_pass"]
    assert "incomplete_routine_rna" in result["failure_reasons"]


def test_combined_screen_with_no_rna_is_unresolved_not_an_alignment_error():
    result = evaluate(
        [],
        [("testis", "cells in seminiferous ducts", "Low")],
        missing=["testis", "heart muscle"],
        adaptive=True,
        selectivity=True,
    )
    assert not result["evidence_pass"]
    assert "incomplete_routine_rna" in result["failure_reasons"]
    assert "unresolved_reproductive_selectivity" in result["failure_reasons"]


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


def test_manageable_lymphoid_expression_is_screened_and_remains_flagged():
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


def test_positive_reproductive_evidence_is_required():
    assert not evaluate([("testis", 0), ("liver", 0)])["evidence_pass"]
    assert evaluate([("testis", 0.5), ("thyroid gland", 3)])["evidence_pass"]
    assert evaluate(
        [("testis", 0), ("thyroid gland", 3)], [("testis", "cells in seminiferous ducts", "Low")]
    )["evidence_pass"]


@pytest.mark.parametrize("value", [1e-12, 0.1, 0.8, 1.0])
def test_any_valid_positive_reproductive_rna_supplies_admission(value):
    result = evaluate([("testis", value), ("heart muscle", 0)], selectivity=True)
    assert result["positive_evidence_reproductive"]
    assert result["evidence_pass"]


@pytest.mark.parametrize("tissue", ["tonsil", "appendix", "gallbladder"])
def test_disposable_expression_is_exempt_but_cannot_supply_reproductive_evidence(tissue):
    only_disposable = evaluate([("testis", 0), (tissue, 10000)])
    assert not only_disposable["evidence_pass"]
    assert not only_disposable["positive_evidence_reproductive"]
    assert "insufficient_positive_reproductive_evidence" in only_disposable["failure_reasons"]
    with_reproductive = evaluate([("testis", 1), (tissue, 10000)])
    assert with_reproductive["evidence_pass"]
    assert with_reproductive["positive_evidence_reproductive"]


def test_disposable_ihc_cannot_supply_reproductive_evidence():
    result = evaluate([("testis", 0), ("skin", 0)], [("skin 1", "melanocytes", "High")])
    assert not result["evidence_pass"]
    assert not result["positive_evidence_reproductive"]
    assert not result["ihc_exclusions"]


def test_adaptive_absolute_boundary_and_positive_veto_remain_independent():
    negative = [("heart muscle", "cardiomyocytes", "Not detected")]
    assert evaluate([("testis", 100), ("heart muscle", 12)], negative, adaptive=True)[
        "evidence_pass"
    ]
    assert not evaluate([("testis", 100), ("heart muscle", 12.1)], negative, adaptive=True)[
        "evidence_pass"
    ]
    result = evaluate(
        [("testis", 100), ("heart muscle", 5), ("liver", 0)],
        [*negative, ("liver", "hepatocytes", "Low")],
        adaptive=True,
    )
    assert not result["evidence_pass"]
    assert "positive_ihc_outside_permitted" in result["failure_reasons"]


@pytest.mark.parametrize("tissue,ceiling", [("salivary gland", 24), ("thyroid gland", 60)])
def test_adaptive_protected_and_manageable_ceilings_are_inclusive(tissue, ceiling):
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


@pytest.mark.parametrize(
    "tissue,cell,reference,value",
    [
        ("heart muscle", "cardiomyocytes", 3, 2),
        ("salivary gland", "glandular cells", 11, 9),
        ("thyroid gland", "glandular cells", 1, 1),
    ],
)
def test_combined_screen_selectivity_boundaries_are_inclusive(tissue, cell, reference, value):
    negative = [(tissue, cell, "Not detected")]
    assert evaluate(
        [("testis", reference), (tissue, value)], negative, adaptive=True, selectivity=True
    )["evidence_pass"]
    result = evaluate(
        [("testis", reference - 0.001), (tissue, value)], negative, adaptive=True, selectivity=True
    )
    assert not result["evidence_pass"]
    assert "reproductive_selectivity_below_threshold" in result["failure_reasons"]
    assert not result["rna_exclusions"]


def test_large_reproductive_or_disposable_signal_cannot_override_other_gates():
    assert evaluate([("testis", 99), ("heart muscle", 1), ("tonsil", 1e6)], selectivity=True)[
        "evidence_pass"
    ]
    result = evaluate([("testis", 1e6), ("heart muscle", 2.01)], selectivity=True)
    assert "rna_above_tissue_ceiling" in result["failure_reasons"]
    assert not result["selectivity_exclusions"]
    result = evaluate(
        [("testis", 1e6), ("heart muscle", 0)],
        [("heart muscle", "cardiomyocytes", "Low")],
        selectivity=True,
    )
    assert "positive_ihc_outside_permitted" in result["failure_reasons"]


def test_zero_reproductive_rna_with_positive_ihc_leaves_selectivity_unresolved():
    result = evaluate(
        [("testis", 0), ("heart muscle", 0)],
        [("testis", "cells in seminiferous ducts", "Low")],
        selectivity=True,
    )
    assert result["positive_evidence_reproductive"]
    assert not result["evidence_pass"]
    assert "unresolved_reproductive_selectivity" in result["failure_reasons"]


def test_other_organ_negative_ihc_cannot_relax_heart_selectivity():
    result = evaluate(
        [("testis", 2), ("heart muscle", 1), ("kidney", 0)],
        [("kidney", "cells in glomeruli", "Not detected")],
        adaptive=True,
        selectivity=True,
    )
    assert not result["evidence_pass"]
    assert "heart muscle" in result["selectivity_exclusions"]
    assert "75%" in result["selectivity_exclusions"]

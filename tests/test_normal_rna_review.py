"""High normal RNA cannot be hidden by enrichment, negative IHC or Variable."""

import pandas as pd
import pytest

from oncoref import normal_rna_review as review


def observation(value=30.0, category="vital", mapping_status="mapped", **updates):
    return pd.DataFrame(
        [
            {
                "gene_id": "ENSG1",
                "source": "hpa_immune_cell",
                "version": "v23",
                "taxonomy_id": 9606,
                "modality": "rna",
                "stage": "adult",
                "unit": "nTPM",
                "value": value,
                "category": category,
                "mapping_status": mapping_status,
                "read_count": 1,
            }
            | updates
        ]
    )


@pytest.mark.parametrize("unit", ["TPM", "nTPM"])
def test_exact_boundary_and_single_measurement(unit):
    assert not review.findings(observation(30, unit=unit)).iloc[0].exclude
    assert review.findings(observation(30.00001, unit=unit)).iloc[0].exclude


@pytest.mark.parametrize("category", ["reproductive", "disposable"])
def test_known_permitted_expression_is_not_forbidden(category):
    assert not review.findings(observation(50000, category)).iloc[0].exclude
    assert review.findings(observation(50000, category, "unmapped")).iloc[0].exclude


@pytest.mark.parametrize(
    "updates",
    [
        {"unit": "CP10K"},
        {"unit": "counts"},
        {"value": float("nan")},
        {"value": -1},
        {"value": float("inf")},
        {"taxonomy_id": 10090},
        {"stage": "fetal"},
        {"modality": "ihc"},
        {"source": ""},
        {"version": None},
    ],
)
def test_invalid_or_incommensurate_evidence_is_rejected(updates):
    with pytest.raises(ValueError):
        review.findings(observation(**updates))


def baseline(category="Neither", detail="Variable"):
    return pd.DataFrame(
        [
            {
                "Gene ID": "ENSG1",
                "Policy": review.BASE_POLICY_ID,
                "Category": category,
                "Detail": detail,
                "Reason": "baseline",
                "Previous category": "CTA",
                "Reproductive positive control": True,
                "Reproductive RNA (nTPM)": 1000,
            }
        ]
    )


def summary(high=True):
    return pd.DataFrame(
        [
            {
                "Gene ID": "ENSG1",
                "Normal RNA assessed": True,
                "High normal RNA": high,
                "High normal RNA findings": "source: 31 nTPM",
                "Normal RNA maximum (TPM)": float("nan"),
                "Normal RNA maximum (nTPM)": 31 if high else 1,
            }
        ]
    )


@pytest.mark.parametrize(
    "category,detail",
    [
        ("Neither", "Variable"),
        ("Neither", "Antibody"),
        ("Neither", "Unassessed"),
        ("CTA", "Pass"),
    ],
)
def test_high_normal_rna_overrides_review_and_reproductive_support(category, detail):
    row = review.classify(baseline(category, detail), summary()).iloc[0]
    assert row.Category == "Unsafe" and row.Detail == "High RNA" and row.Forbidden
    assert not row.Target
    assert row["Previous category"] == category


def test_existing_safety_reason_is_retained_and_no_gene_is_promoted():
    assert review.classify(baseline("Unsafe", "Safety"), summary()).iloc[0].Detail == "Safety"
    assert review.classify(baseline(), summary(False)).iloc[0].Detail == "Variable"


def test_incomplete_audit_cannot_silently_clear_neither():
    with pytest.raises(ValueError, match="Every CTA and Neither"):
        review.classify(baseline(), summary().iloc[:0])


def test_summary_flag_cannot_hide_high_measurement():
    with pytest.raises(ValueError, match="disagree"):
        review.classify(baseline(), summary().assign(**{"High normal RNA": False}))

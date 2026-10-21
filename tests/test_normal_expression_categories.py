"""Expression categories preserve non-CTA genes and conservative sequence evidence."""

import pandas as pd
import pytest
from cta_permitted_report import annotate_archive

from oncoref.normal_expression_categories import classify_normal_expression


def classify(rna, ihc=(), missing=()):
    records = [
        {
            "modality": "bulk_rna",
            "tissue": tissue,
            "cell_type": "",
            "rna_ntpm": value,
            "ihc_level": "",
            "ihc_reliability": "",
        }
        for tissue, value in rna
    ] + [
        {
            "modality": "ihc",
            "tissue": tissue,
            "cell_type": cell,
            "rna_ntpm": "",
            "ihc_level": level,
            "ihc_reliability": "Enhanced",
        }
        for tissue, cell, level in ihc
    ]
    frame = pd.DataFrame(records).assign(Ensembl_Gene_ID="ENSG00000999999", Symbol="EXAMPLE")
    rows = annotate_archive(frame)
    return classify_normal_expression(
        rows,
        routine_rna_contexts={(t, "") for t, _ in rna} | {(t, "") for t in missing},
        routine_ihc_contexts={(t, c) for t, c, _ in ihc},
    )


def test_reproductive_expression_does_not_require_a_publication_nomination():
    result = classify([("testis", 100), ("heart muscle", 0)])
    assert result["category"] == "Reproductive"
    assert not result["sequence_specificity_established"]


def test_disposable_without_reproductive_expression_is_not_dropped():
    result = classify([("testis", 0), ("tonsil", 10), ("heart muscle", 0.8)])
    assert result["category"] == "Other" and result["other_subgroup"] == "Disposable"
    assert result["normal_expression_observed"] and result["conservative_normal_reference"]
    assert result["category_meaning"] == "below-threshold expression"


def test_silent_is_below_threshold_not_literal_zero_and_ihc_remains_unassessed():
    result = classify([("testis", 0), ("tonsil", 0), ("heart muscle", 0.5)])
    assert result["category"] == "Other" and result["other_subgroup"] == "Silent"
    assert result["ihc_coverage"] == "unavailable"
    assert result["protein_absence_evidence"] == "not_established"


def test_missing_required_rna_is_unassessed_never_silent():
    result = classify([("testis", 0)], missing=["heart muscle"])
    assert result["category"] == "Unassessed" and result["other_subgroup"] == ""
    assert result["conservative_normal_reference"]


@pytest.mark.parametrize("reproductive", [0, 100])
def test_any_positive_non_disposable_ihc_places_gene_in_normal(reproductive):
    result = classify([("testis", reproductive), ("liver", 0)], [("liver", "hepatocytes", "Low")])
    assert result["category"] == "Normal"
    assert "positive_ihc" in result["failure_reasons"]


def test_other_uses_same_local_ihc_absolute_ladder_without_reproductive_ratio():
    rna = [("testis", 0), ("heart muscle", 5)]
    assert classify(rna)["category"] == "Normal"
    relaxed = classify(rna, [("heart muscle", "cardiomyocytes", "Not detected")])
    assert relaxed["category"] == "Other"
    assert not relaxed["failure_reasons"]


def test_protein_only_reproductive_positive_cannot_invent_selectivity():
    result = classify(
        [("testis", 0), ("heart muscle", 0)], [("testis", "cells in seminiferous ducts", "Low")]
    )
    assert result["category"] == "Unassessed"
    assert "unresolved_reproductive_selectivity" in result["quality_flags"]


def test_reproductive_screen_label_does_not_erase_normal_sequence_reference_evidence():
    result = classify([("testis", 100), ("heart muscle", 0.1)])
    assert result["category"] == "Reproductive"
    assert result["normal_expression_observed"] and result["conservative_normal_reference"]

"""Admission evidence must never override normal-cell counterevidence."""

import pandas as pd
import pytest

from oncoref import normal_cell_screen as screen


def base(category="Neither", detail="Low"):
    return pd.DataFrame(
        {
            "Gene": ["G"],
            "Gene ID": ["ENSG1"],
            "Policy": [screen.BASE_POLICY_ID],
            "Category": [category],
            "Detail": [detail],
            "Previous CTA screen": ["Retained"],
            "Reproductive RNA (nTPM)": [1.2],
            "Reproductive positive control": [False],
            "Reason": ["baseline"],
        }
    )


def summary(**flags):
    row = {"Gene ID": "ENSG1"}
    row.update(
        dict.fromkeys(
            [
                "Cell RNA review",
                "Immune RNA review",
                "Corneal RNA review",
                "Antibody review",
                "Additional reproductive protein",
                "Cell evidence assessed",
            ],
            False,
        )
    )
    row.update(
        dict.fromkeys(
            [
                "Cell RNA findings",
                "Immune RNA findings",
                "Corneal RNA findings",
                "Antibody findings",
            ],
            "",
        )
    )
    return pd.DataFrame([row | flags])


def papers(*findings):
    return pd.DataFrame(
        [
            {"Gene ID": "ENSG1", "Finding": f, "Source": "https://example.org/study"}
            for f in findings
        ],
        columns=["Gene ID", "Finding", "Source"],
    )


def test_good_reproductive_protein_bypasses_rna_floor():
    result = screen.classify(base(), summary(**{"Additional reproductive protein": True}), papers())
    assert result.iloc[0].Category == "CTA"


def test_reproductive_and_tumor_facts_support_low_rna_but_tumor_alone_does_not():
    assert (
        screen.classify(base(), summary(), papers("reproductive", "recognition")).iloc[0].Category
        == "CTA"
    )
    assert screen.classify(base(), summary(), papers("tumor")).iloc[0].Category == "Neither"


@pytest.mark.parametrize(
    "category,detail",
    [
        ("Unsafe", "Protein"),
        ("Unsafe", "RNA"),
        ("Unsafe", "Safety"),
        ("Neither", "Variable"),
        ("Neither", "Unassessed"),
    ],
)
def test_support_does_not_rescue_prior_normal_exclusions_or_missing_data(category, detail):
    result = screen.classify(
        base(category, detail), summary(), papers("reproductive", "recognition")
    )
    assert result.iloc[0].Category == category and result.iloc[0].Detail == detail


def test_cell_rna_or_uncertain_ihc_hold_is_not_forbidden():
    for flag in ["Immune RNA review", "Cell RNA review", "Corneal RNA review", "Antibody review"]:
        result = screen.classify(base(), summary(**{flag: True}), papers("reproductive", "tumor"))
        assert result.iloc[0].Category == "Neither" and not result.iloc[0].Forbidden


def test_normal_protein_dominates_support_and_review():
    result = screen.classify(
        base(),
        summary(**{"Immune RNA review": True}),
        papers("reproductive", "recognition", "normal_protein"),
    )
    assert result.iloc[0].Category == "Unsafe" and result.iloc[0].Forbidden


def test_mobile_immune_cells_do_not_inherit_reproductive_permission():
    frame = pd.DataFrame(
        {
            "source": ["hpa_single_cell_tissue"] * 2,
            "version": ["v23"] * 2,
            "taxonomy_id": [9606] * 2,
            "tissue": ["fallopian tube", "testis"],
            "cell_type": ["granulocytes", "spermatogonia"],
            "gene_id": ["ENSG1"] * 2,
            "unit": ["nTPM"] * 2,
            "value": [30.0, 20.0],
            "read_count": [100, 100],
        }
    )
    result = screen.cluster_findings(frame)
    assert result.iloc[0].category == "vital" and result.iloc[0].review
    assert result.iloc[0].reproductive_cell_max == 20
    with pytest.raises(ValueError, match="human-specific"):
        screen.cluster_findings(frame.assign(taxonomy_id=10090))


def test_dice_recurrence_and_native_units():
    frame = pd.DataFrame(
        {
            "source": ["dice_cd4_stim"] * 4,
            "gene_id": ["ENSG1"] * 4,
            "gene_name": ["G"] * 4,
            "cell_type": ["CD4_STIM"] * 4,
            "state": ["activated"] * 4,
            "unit": ["TPM"] * 4,
            "taxonomy_id": [9606] * 4,
            "value": [5.0, 5.0, 5.0, 0.0],
        }
    )
    assert screen.dice_findings(frame).iloc[0].review
    assert not screen.dice_findings(frame.iloc[1:]).iloc[0].review
    with pytest.raises(ValueError, match="original DICE"):
        screen.dice_findings(frame.assign(unit="nTPM"))
    # HPA's summary of DICE does not count as independent replication.
    result = screen.immune_findings(frame.assign(source="hpa_immune_schmiedel", value=50))
    assert not result.review.any()


def test_cornea_requires_supported_groups_and_keeps_development_separate():
    frame = pd.DataFrame(
        {
            "unit": ["CP10K"] * 2,
            "stage": ["adult"] * 2,
            "taxonomy_id": [9606] * 2,
            "gene_name": ["G"] * 2,
            "cell_type": ["limbal cell"] * 2,
            "value": [1.0, 1.0],
            "positive_cells": [5, 5],
            "read_count": [20, 20],
        }
    )
    assert screen.cornea_findings(frame).iloc[0].review
    assert not screen.cornea_findings(frame.iloc[:1]).iloc[0].review
    with pytest.raises(ValueError, match="Adult corneal"):
        screen.cornea_findings(frame.assign(stage="fetal"))

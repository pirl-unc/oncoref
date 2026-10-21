"""Adaptive absolute ceilings must depend on local, complete negative evidence."""

import pandas as pd
import pytest

from oncoref.cta_normal_screen import adaptive_normal_rna_limits
from oncoref.normal_tissues import normal_tissue_screening_policy


def observations(tissue="heart muscle", value=5, ihc=()):
    rows = [
        {
            "Ensembl_Gene_ID": "ENSG_EXAMPLE",
            "modality": "bulk_rna",
            "source_name": "hpa_rna_consensus",
            "source_version": "v23",
            "tissue": tissue,
            "cell_type": "",
            "rna_ntpm": value,
            "ihc_level": "",
            "ihc_reliability": "",
            "valid_measurement": True,
        }
    ]
    for site, cell, level, reliability in ihc:
        rows.append(
            {
                **rows[0],
                "modality": "ihc",
                "source_name": "hpa_normal_tissue",
                "tissue": site,
                "cell_type": cell,
                "rna_ntpm": float("nan"),
                "ihc_level": level,
                "ihc_reliability": reliability,
            }
        )
    return pd.DataFrame(rows)


@pytest.mark.parametrize(
    "reliability,vital,protected,manageable",
    [
        ("Enhanced", 6, 12, 30),
        ("Supported", 3, 6, 15),
        ("Approved", 2, 4, 10),
        ("Uncertain", 1, 2, 5),
    ],
)
def test_approved_ladder_uses_local_negative_reliability(reliability, vital, protected, manageable):
    for tissue, cell, expected in (
        ("heart muscle", "cardiomyocytes", vital),
        ("salivary gland", "glandular cells", protected),
        ("thyroid gland", "glandular cells", manageable),
    ):
        rows = observations(tissue, ihc=[(tissue, cell, "Not detected", reliability)])
        result = adaptive_normal_rna_limits(rows, routine_ihc_contexts=[(tissue, cell)])
        assert result.iloc[0].rna_ceiling_ntpm == expected
        assert result.iloc[0].ihc_reliability_used == reliability


def test_missing_uncertain_unknown_or_unreliable_evidence_keeps_baseline():
    expected = [("heart muscle", "cardiomyocytes")]
    rows = observations()
    assert (
        adaptive_normal_rna_limits(rows, routine_ihc_contexts=expected).iloc[0].rna_ceiling_ntpm
        == 1
    )
    for reliability in ["Uncertain", "", "not a validation tier"]:
        rows = observations(ihc=[("heart muscle", "cardiomyocytes", "Not detected", reliability)])
        assert (
            adaptive_normal_rna_limits(rows, routine_ihc_contexts=expected).iloc[0].rna_ceiling_ntpm
            == 1
        )
    rows = observations(ihc=[("heart muscle", "cardiomyocytes", "Not detected", "Enhanced")])
    result = adaptive_normal_rna_limits(
        rows, routine_ihc_contexts=expected, reviewed_unreliable=True
    )
    assert result.iloc[0].rna_ceiling_ntpm == 1
    rows.loc[1, "valid_measurement"] = False
    assert (
        adaptive_normal_rna_limits(rows, routine_ihc_contexts=expected).iloc[0].ihc_evidence_status
        == "invalid_observation"
    )


def test_every_routine_cell_context_is_required_and_weakest_reliability_wins():
    expected = [("kidney", "cells in glomeruli"), ("kidney", "cells in tubules")]
    rows = observations(
        "kidney", ihc=[("kidney", "cells in glomeruli", "Not detected", "Enhanced")]
    )
    result = adaptive_normal_rna_limits(rows, routine_ihc_contexts=expected)
    assert result.iloc[0].rna_ceiling_ntpm == 1
    assert result.iloc[0].ihc_evidence_status == "incomplete_ihc_contexts"
    rows = observations(
        "kidney",
        ihc=[
            ("kidney", "cells in glomeruli", "Not detected", "Enhanced"),
            ("kidney", "cells in tubules", "Not detected", "Supported"),
        ],
    )
    result = adaptive_normal_rna_limits(rows, routine_ihc_contexts=expected)
    assert result.iloc[0].rna_ceiling_ntpm == 3


@pytest.mark.parametrize("level", ["Low", "Medium", "High", "Ascending", "Descending"])
def test_any_local_positive_including_extra_cell_context_prevents_relaxation(level):
    rows = observations(
        ihc=[
            ("heart muscle", "cardiomyocytes", "Not detected", "Enhanced"),
            ("heart muscle", "other cells", level, "Enhanced"),
        ]
    )
    result = adaptive_normal_rna_limits(
        rows, routine_ihc_contexts=[("heart muscle", "cardiomyocytes")]
    )
    assert result.iloc[0].rna_ceiling_ntpm == 1
    assert result.iloc[0].ihc_evidence_status == "positive_ihc"


def test_other_organs_and_narrower_regions_cannot_clear_rna_tissue():
    for tissue, ihc_tissue in [("heart muscle", "kidney"), ("basal ganglia", "caudate")]:
        rows = observations(
            tissue, ihc=[(ihc_tissue, "neuronal cells", "Not detected", "Enhanced")]
        )
        result = adaptive_normal_rna_limits(
            rows, routine_ihc_contexts=[(ihc_tissue, "neuronal cells")]
        )
        assert result.iloc[0].rna_ceiling_ntpm == 1
        assert result.iloc[0].ihc_evidence_status == "no_routine_ihc_contexts"


def test_equivalent_tissue_labels_require_all_named_source_contexts():
    expected = [("stomach 1", "glandular cells"), ("stomach 2", "glandular cells")]
    rows = observations(
        "stomach",
        ihc=[
            ("stomach 1", "glandular cells", "Not detected", "Enhanced"),
            ("stomach 2", "glandular cells", "Not detected", "Enhanced"),
        ],
    )
    result = adaptive_normal_rna_limits(rows, routine_ihc_contexts=expected)
    assert result.iloc[0].rna_ceiling_ntpm == 12
    assert result.iloc[0].required_ihc_context_count == 2
    assert (
        adaptive_normal_rna_limits(rows.iloc[:2], routine_ihc_contexts=expected)
        .iloc[0]
        .rna_ceiling_ntpm
        == 2
    )


def test_permitted_and_unresolved_contexts_have_no_fabricated_ceiling():
    for tissue, status in [("testis", "permitted"), ("head", "unresolved")]:
        result = adaptive_normal_rna_limits(observations(tissue), routine_ihc_contexts=[])
        assert pd.isna(result.iloc[0].rna_ceiling_ntpm)
        assert result.iloc[0].ihc_evidence_status == status


def test_source_species_and_gene_boundaries_are_enforced():
    rows = observations()
    with pytest.raises(ValueError, match="human"):
        adaptive_normal_rna_limits(rows, routine_ihc_contexts=[], species="dog")
    rows.loc[0, "source_version"] = "v25"
    with pytest.raises(ValueError, match="source or version"):
        adaptive_normal_rna_limits(rows, routine_ihc_contexts=[])
    rows = pd.concat([observations(), observations()], ignore_index=True)
    rows.loc[1, "Ensembl_Gene_ID"] = "ENSG_OTHER"
    with pytest.raises(ValueError, match="one gene"):
        adaptive_normal_rna_limits(rows, routine_ihc_contexts=[])


def test_policy_inherits_anatomy_without_changing_fixed_limits():
    fixed = normal_tissue_screening_policy()
    adaptive = normal_tissue_screening_policy("ihc_adaptive_cta")
    assert fixed["vital"] == adaptive["vital"]
    assert fixed["protected"] == adaptive["protected"]
    assert fixed["manageable"] == adaptive["manageable"]
    assert "negative_ihc_rna_limits" not in fixed
    adaptive["negative_ihc_rna_limits"]["Enhanced"]["vital"] = 999
    assert (
        normal_tissue_screening_policy("ihc_adaptive_cta")["negative_ihc_rna_limits"]["Enhanced"][
            "vital"
        ]
        == 6
    )


def test_bulk_skin_cannot_borrow_negativity_from_selected_compartments():
    required = [
        ("skin 1", "keratinocytes"),
        ("skin 1", "melanocytes"),
        ("skin 2", "epidermal cells"),
    ]
    negatives = [(tissue, cell, "Not detected", "Enhanced") for tissue, cell in required]
    rows = observations("skin", ihc=negatives)
    complete = adaptive_normal_rna_limits(rows, routine_ihc_contexts=required)
    assert complete.iloc[0].screening_tier == "manageable"
    assert complete.iloc[0].rna_ceiling_ntpm == 30
    assert complete.iloc[0].screening_policy_id.endswith("2026-10-20.4")
    missing = adaptive_normal_rna_limits(rows.iloc[:-1], routine_ihc_contexts=required)
    assert missing.iloc[0].rna_ceiling_ntpm == 5
    # Allowed melanocyte staining cannot establish a negative bulk-skin context.
    rows.loc[2, "ihc_level"] = "Low"
    positive = adaptive_normal_rna_limits(rows, routine_ihc_contexts=required)
    assert positive.iloc[0].rna_ceiling_ntpm == 5
    assert positive.iloc[0].ihc_evidence_status == "positive_ihc"


def test_missing_rna_preserves_audit_schema_without_inventing_measurements():
    rows = observations(ihc=[("heart muscle", "cardiomyocytes", "Not detected", "Enhanced")])
    expected = [("heart muscle", "cardiomyocytes")]
    full = adaptive_normal_rna_limits(rows, routine_ihc_contexts=expected)
    for no_rna in (rows.iloc[1:], rows.iloc[:0]):
        empty = adaptive_normal_rna_limits(no_rna, routine_ihc_contexts=expected)
        assert empty.empty
        assert empty.columns.equals(full.columns)
        assert empty.relaxed.sum() == 0

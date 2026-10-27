"""Cell identity changes risk tiers without removing evidence or tissue limits."""

import numpy as np
import pandas as pd
import pytest

from oncoref.compartment_rna_policy import assess_contexts, assess_imported_protein_reviews


def context(**changes):
    return {
        "source": "hpa_single_cell_tissue",
        "taxonomy_id": 9606,
        "stage": "adult",
        "tissue": "heart muscle",
        "cell_type": "fibroblasts",
        "category": "vital",
        "mapping_status": "mapped",
        "negative_ihc_grade": "Missing/uncertain",
        "unit": "nTPM",
        "value": 12.0,
    } | changes


def test_specific_compartments_leave_critical_neighbors_strict():
    d = pd.DataFrame(
        [
            context(),
            context(cell_type="cardiomyocytes"),
            context(cell_type="endothelial cells"),
            context(
                tissue="colon", cell_type="undifferentiated cells", category="protected", value=20.3
            ),
            context(tissue="colon", cell_type="enteric neurons", category="vital", value=20.3),
            context(tissue="bronchus", cell_type="basal respiratory cells", value=24.9),
        ]
    )
    r = assess_contexts(d)
    assert r.candidate_ceiling.tolist() == [16, 4, 4, 32, 4, 4]
    assert r.candidate_exclusion.tolist() == [False, True, True, False, True, True]
    assert r.stem_progenitor_caution.tolist() == [False, False, False, True, False, False]
    pd.testing.assert_frame_equal(d, r[d.columns])


def test_risk_and_ihc_axes_and_native_units_remain_linked():
    d = pd.DataFrame(
        [
            context(
                source={
                    "nTPM": "hpa_immune_cell",
                    "TPM": "hpa_immune_monaco",
                    "CP10K": "cornea_collin_adult_counts",
                }[unit],
                category=tier,
                cell_type="immune",
                unit=unit,
                negative_ihc_grade=grade,
            )
            for tier in ["vital", "protected", "manageable"]
            for unit in ["nTPM", "TPM", "CP10K"]
            for grade in ["Missing/uncertain", "Approved", "Supported", "Enhanced"]
        ]
    )
    r = assess_contexts(d)
    for _, part in r.groupby(["category", "unit"], sort=False):
        assert part.candidate_ceiling.iloc[1] == 2 * part.candidate_ceiling.iloc[0]
        assert part.candidate_ceiling.iloc[3] == 8 * part.candidate_ceiling.iloc[0]
    assert r.strong_boundary.eq(2 * r.candidate_ceiling).all()


def test_unknown_contexts_missing_values_and_positive_protein_are_not_erased():
    d = pd.DataFrame(
        [
            context(mapping_status="unmapped"),
            context(value=np.nan),
            context(value=0, ihc_positive=True),
            context(value=16),
            context(value=32),
            context(category="reproductive", tissue="testis", value=100),
        ]
    )
    d["ihc_positive"] = d.ihc_positive.fillna(False).astype(bool)
    r = assess_contexts(d)
    assert r.candidate_exclusion.tolist() == [True, False, True, False, True, False]
    assert r.unsafe_evidence.tolist() == [False, False, True, False, True, False]
    assert pd.isna(r.value.iloc[1])
    with pytest.raises(ValueError, match="adult human"):
        assess_contexts(pd.DataFrame([context(stage="fetal")]))


def test_protein_intensity_rule_is_independent_of_gene_name_and_source_access():
    d = pd.DataFrame(
        {
            "Gene": ["a", "b", "c", "d"],
            "Review ID": ["a", "b", "c", "d"],
            "Finding": "Normal protein",
            "Status": "Accepted",
            "Review depth": [
                "Primary full text",
                "Primary full text and abstract",
                "Primary abstract",
                "Migrated prior cited assessment",
            ],
        }
    )
    intensities = pd.DataFrame(
        {
            "Review ID": ["a", "b", "c", "d"],
            "protein_level": ["Low", "High", "Low", "Detected_unquantified"],
            "intensity_source": "primary_source",
            "intensity_locator": "results",
        }
    )
    r = assess_imported_protein_reviews(d, intensities)
    assert r.retained_protein_exclusion.tolist() == [False, True, False, True]
    with pytest.raises(ValueError, match="exactly once"):
        assess_imported_protein_reviews(d, intensities.iloc[:3])


def test_low_protein_is_flagged_but_cannot_erase_rna_or_gain_negative_ihc_credit():
    d = pd.DataFrame(
        [
            context(value=0, ihc_positive=True, Level="Low"),
            context(value=100, ihc_positive=True, Level="Low"),
            context(value=0, ihc_positive=True, Level="Medium"),
        ]
    )
    r = assess_contexts(d)
    assert r.low_normal_protein_flag.tolist() == [True, True, False]
    assert r.candidate_exclusion.tolist() == [False, True, True]
    assert r.unsafe_evidence.tolist() == [False, True, True]
    assert r.negative_ihc_grade.eq("Missing/uncertain").all()

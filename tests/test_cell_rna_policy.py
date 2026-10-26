"""Two decisions, anatomical context and uncertainty must remain independent."""

import numpy as np
import pandas as pd
import pytest

from oncoref import cell_rna_policy as policy


def context(**changes):
    row = {
        "gene_id": "ENSG1",
        "source": "hpa_single_cell_tissue",
        "version": "v23",
        "taxonomy_id": 9606,
        "tissue": "heart muscle",
        "cell_type": "cardiomyocytes",
        "stage": "adult",
        "state": "unspecified",
        "unit": "nTPM",
        "value": 3.0,
        "category": "vital",
        "mapping_status": "mapped",
        "ihc_positive": False,
        "negative_ihc_grade": "Missing/uncertain",
    }
    return pd.DataFrame([row | changes])


def thresholds():
    return pd.DataFrame(
        [
            {
                "assay": "hpa_single_cell_tissue",
                "unit": "nTPM",
                "category": "vital",
                "cta_cutoff": 2.0,
                "reference_cutoff": 5.0,
            },
            {
                "assay": "hpa_single_cell_tissue",
                "unit": "nTPM",
                "category": "manageable",
                "cta_cutoff": 10.0,
                "reference_cutoff": 20.0,
            },
        ]
    )


FACTORS = {"Missing/uncertain": 1, "Approved": 1.25, "Supported": 1.5, "Enhanced": 2}


@pytest.mark.parametrize(
    "value,assessment,cta,reference",
    [
        (1.99, "Below CTA cutoff", False, False),
        (2, "CTA exclusion", True, False),
        (4.99, "CTA exclusion", True, False),
        (5, "RNA exclusion", True, True),
        (np.nan, "Unassessed", False, False),
    ],
)
def test_two_cutoffs_and_exact_boundaries(value, assessment, cta, reference):
    r = policy.assess(context(value=value), thresholds(), ihc_factors=FACTORS).iloc[0]
    assert (r.assessment, r.cta_excluded, r.reference_review) == (assessment, cta, reference)


def test_anatomical_tier_and_negative_ihc_change_thresholds():
    assert (
        policy.assess(context(category="manageable"), thresholds(), ihc_factors=FACTORS)
        .iloc[0]
        .assessment
        == "Below CTA cutoff"
    )
    r = policy.assess(
        context(negative_ihc_grade="Enhanced"), thresholds(), ihc_factors=FACTORS
    ).iloc[0]
    assert r.cta_cutoff == 4 and r.reference_cutoff == 10
    assert r.assessment == "Below CTA cutoff"


def test_reproductive_maximum_does_not_rescue_normal_signal():
    r = policy.assess(
        context(value=100, reproductive_cell_max=10000), thresholds(), ihc_factors=FACTORS
    )
    assert r.iloc[0].assessment == "RNA exclusion"


def test_missing_anatomy_and_positive_ihc_are_not_negative_evidence():
    r = policy.assess(
        context(mapping_status="unmapped", value=0), thresholds(), ihc_factors=FACTORS
    )
    assert r.iloc[0].assessment == "Unassessed"
    r = policy.assess(context(ihc_positive=True, value=np.nan), thresholds(), ihc_factors=FACTORS)
    assert r.iloc[0].assessment == "Protein exclusion"


@pytest.mark.parametrize(
    "changes",
    [
        {"unit": "TPM"},
        {"taxonomy_id": 10090},
        {"stage": "fetal"},
        {"source": "hpa_immune_schmiedel"},
        {"value": -1},
        {"value": np.inf},
    ],
)
def test_no_unit_species_stage_or_duplicate_study_transfer(changes):
    with pytest.raises(ValueError):
        policy.assess(context(**changes), thresholds(), ihc_factors=FACTORS)


def ihc(**changes):
    return pd.DataFrame(
        [
            {
                "Gene": "ENSG1",
                "Tissue": "heart muscle",
                "Cell type": "cardiomyocytes",
                "Level": "Not detected",
                "Reliability": "Enhanced",
            }
            | changes
        ]
    )


@pytest.mark.parametrize(
    "change,controls,conflicts",
    [
        ({"cell_type": "fibroblasts"}, {"ENSG1"}, ()),
        ({"state": "activated"}, {"ENSG1"}, ()),
        ({"source": "hpa_immune_cell"}, {"ENSG1"}, ()),
        ({}, set(), ()),
        ({}, {"ENSG1"}, ("ENSG1",)),
    ],
)
def test_unmatched_cells_states_missing_controls_and_conflicts_get_no_relaxation(
    change, controls, conflicts
):
    r = policy.match_negative_ihc(
        context(**change), ihc(), positive_controls=controls, conflicting_genes=conflicts
    )
    assert r.iloc[0].negative_ihc_grade == "Missing/uncertain"


def test_exact_matched_negative_requires_reliable_grade():
    r = policy.match_negative_ihc(context(), ihc(), positive_controls={"ENSG1"})
    assert r.iloc[0].negative_ihc_grade == "Enhanced"
    r = policy.match_negative_ihc(
        context(), ihc(Reliability="Uncertain"), positive_controls={"ENSG1"}
    )
    assert r.iloc[0].negative_ihc_grade == "Missing/uncertain"


def test_positive_ihc_is_not_relaxed():
    r = policy.match_negative_ihc(context(value=0), ihc(Level="Low"), positive_controls={"ENSG1"})
    result = policy.assess(r, thresholds(), ihc_factors=FACTORS).iloc[0]
    assert result.assessment == "Protein exclusion"


def test_require_complete_valid_explicit_configuration():
    for t in [
        thresholds().iloc[0:0],
        thresholds().assign(reference_cutoff=1),
        pd.concat([thresholds()] * 2),
    ]:
        with pytest.raises(ValueError):
            policy.assess(context(), t, ihc_factors=FACTORS)
    with pytest.raises(ValueError):
        policy.assess(context(), thresholds(), ihc_factors=FACTORS | {"Enhanced": 0.5})


def test_recurrence_counts_distinct_groups_and_missing_is_not_zero():
    from scripts.cta_cell_policy_report import summarize

    d = pd.concat(
        [context(value=v).assign(sample=s) for s, v in [(1, 10), (2, 2), (3, 0)]], ignore_index=True
    )
    assert summarize(d, support=d.value.ge(1), rank=2, unique="sample").value.iloc[0] == 2
    assert pd.isna(summarize(d, support=d.value.ge(1), rank=3, unique="sample").value.iloc[0])
    with pytest.raises(ValueError, match="inflate recurrence"):
        summarize(pd.concat([d, d.iloc[[0]]]), support=True, rank=2, unique="sample")

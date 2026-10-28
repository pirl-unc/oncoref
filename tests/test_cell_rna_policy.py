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
    "tissue", ["testis", "placenta", "ovary", "prostate", "breast", "endometrium", "fallopian tube"]
)
def test_reproductive_immune_neutrality_is_local_and_preserves_measurement(tissue):
    d = context(tissue=tissue, cell_type="macrophages", concept_id="blood", value=99)
    neutral = policy.reproductive_immune_contexts(d)
    assert neutral.neutral_reproductive_immune.iloc[0]
    assert neutral.prior_context_category.iloc[0] == "vital"
    assert neutral.value.iloc[0] == 99
    assessed = policy.assess(neutral, thresholds(), ihc_factors=FACTORS)
    assert not assessed.cta_excluded.iloc[0] and not assessed.reference_review.iloc[0]
    assert assessed.assessment.iloc[0].startswith("Neutral")
    assert d.category.iloc[0] == "vital"
    with pytest.raises(ValueError, match="already"):
        policy.reproductive_immune_contexts(neutral)


@pytest.mark.parametrize(
    "tissue,stage",
    [
        ("blood", "adult"),
        ("heart muscle", "adult"),
        ("lung", "juvenile"),
        ("unknown", "unspecified"),
    ],
)
def test_immune_lineage_elsewhere_does_not_inherit_reproductive_permission(tissue, stage):
    d = context(tissue=tissue, stage=stage, cell_type="macrophages", concept_id="blood", value=99)
    result = policy.reproductive_immune_contexts(d)
    assert not result.neutral_reproductive_immune.iloc[0]
    assert result.category.iloc[0] == "vital"


@pytest.mark.parametrize("stage", ["embryonic", "fetal"])
def test_explicit_developmental_immune_context_is_neutral_without_relabeling_age(stage):
    d = context(tissue="liver", stage=stage, cell_type="t-cells", concept_id="blood")
    result = policy.reproductive_immune_contexts(d)
    assert result.neutral_reproductive_immune.iloc[0]
    assert result.neutral_basis.iloc[0] == "explicit_developmental_stage"
    assert result.stage.iloc[0] == stage


def test_reproductive_site_alone_does_not_invent_immune_identity():
    d = context(tissue="testis", cell_type="unknown", concept_id="")
    assert not policy.reproductive_immune_contexts(d).neutral_reproductive_immune.iloc[0]
    d["source_cell_group"] = "Blood & immune cells"
    assert policy.reproductive_immune_contexts(d).neutral_reproductive_immune.iloc[0]


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


@pytest.mark.parametrize(
    "value,excluded,review", [(6, False, False), (6.01, True, False), (12, True, True)]
)
def test_linked_immune_boundary_is_shared_and_inclusive(value, excluded, review):
    d = context(source="dice_cd4_stim", unit="TPM", cell_type="CD4", value=value)
    d = policy.protected_immune_contexts(d)
    r = policy.assess(
        d,
        policy.linked_thresholds(),
        ihc_factors=dict.fromkeys(FACTORS, 1),
        inclusive_ceiling=True,
    ).iloc[0]
    assert r.category == "protected" and r.cta_cutoff == 6 and r.reference_cutoff == 12
    assert r.cta_excluded == excluded and r.reference_review == review


def test_linked_tissue_grade_values_share_registry_without_assay_specific_caps():
    from oncoref.normal_tissues import normal_tissue_screening_policy

    grid = policy.linked_thresholds()
    ladder = normal_tissue_screening_policy("ihc_adaptive_cta")["negative_ihc_rna_limits"]
    for r in grid.itertuples():
        grade = "Missing" if r.negative_ihc_grade == "Missing/uncertain" else r.negative_ihc_grade
        expected = ladder[grade][r.category]
        assert r.nominal_per_million_ceiling == expected
        assert r.cta_cutoff == expected * r.native_units_per_nominal_unit
        assert r.reference_cutoff == 2 * r.cta_cutoff
    with pytest.raises(ValueError, match="cannot also"):
        policy.assess(context(), grid, ihc_factors=FACTORS)


def test_linked_immune_tier_preserves_neutrality_nonimmune_anatomy_and_missing_mapping():
    d = pd.concat(
        [
            context(cell_type="fibroblasts", concept_id="heart"),
            context(cell_type="macrophages", concept_id="blood"),
            context(tissue="testis", cell_type="macrophages", concept_id="blood"),
            context(mapping_status="unmapped", concept_id="", category="unresolved"),
        ],
        ignore_index=True,
    )
    r = policy.protected_immune_contexts(policy.reproductive_immune_contexts(d))
    assert r.category.tolist() == [
        "vital",
        "protected",
        policy.NEUTRAL_REPRODUCTIVE_IMMUNE,
        "unresolved",
    ]
    with pytest.raises(ValueError, match="already"):
        policy.protected_immune_contexts(r)


def test_linked_cornea_rescales_counts_without_relabeling_tpm():
    d = context(source="cornea_collin_adult_counts", unit="CP10K", value=0.02)
    r = policy.assess(
        d,
        policy.linked_thresholds(),
        ihc_factors=dict.fromkeys(FACTORS, 1),
        inclusive_ceiling=True,
    ).iloc[0]
    assert r.unit == "CP10K" and r.cta_cutoff == 0.02 and not r.cta_excluded
    assert r.nominal_per_million_ceiling == 2


def test_two_axis_grid_doubles_both_axes_and_shares_bulk_cell_boundaries():
    ladder = policy.threshold_ladder("two_axis_doubling")
    grades = ["Missing", *policy.GRADES]
    tiers = ["vital", "protected", "manageable"]
    matrix = np.array([[ladder[g][t] for g in grades] for t in tiers])
    assert matrix.tolist() == [[2, 4, 8, 16], [4, 8, 16, 32], [8, 16, 32, 64]]
    assert np.array_equal(matrix[:, 1:], 2 * matrix[:, :-1])
    assert np.array_equal(matrix[1:, :], 2 * matrix[:-1, :])
    for grade in grades:
        bulk = pd.DataFrame(
            {
                "screening_tier": tiers,
                "ihc_reliability_used": grade,
                "profile_positive_control": True,
                "rna_ntpm": [ladder[grade][t] for t in tiers],
            }
        )
        assessed = policy.apply_bulk_ladder(bulk, mode="two_axis_doubling")
        assert not assessed.local_exclusion.any()
        cell_grid = policy.linked_thresholds(mode="two_axis_doubling")
        g = "Missing/uncertain" if grade == "Missing" else grade
        c = cell_grid[
            cell_grid.assay.eq("hpa_single_cell_tissue") & cell_grid.negative_ihc_grade.eq(g)
        ].set_index("category")
        assert assessed.effective_ceiling_ntpm.tolist() == c.loc[tiers, "cta_cutoff"].tolist()
    with pytest.raises(ValueError, match="Unknown"):
        policy.threshold_ladder("unreviewed")


def test_bulk_geometric_grade_requires_positive_profile_and_honors_cap():
    d = pd.DataFrame(
        {
            "screening_tier": ["protected", "protected", "unresolved"],
            "ihc_reliability_used": "Enhanced",
            "profile_positive_control": [True, False, True],
            "rna_ntpm": [9, 9, 0],
        }
    )
    r = policy.apply_bulk_ladder(d, mode="two_axis_doubling", max_ihc_grade="Approved")
    assert r.effective_ceiling_ntpm.iloc[:2].tolist() == [8, 4]
    assert r.local_exclusion.all()


def test_recurrence_counts_distinct_groups_and_missing_is_not_zero():
    from scripts.cta_cell_policy_report import summarize

    d = pd.concat(
        [context(value=v).assign(sample=s) for s, v in [(1, 10), (2, 2), (3, 0)]], ignore_index=True
    )
    assert summarize(d, support=d.value.ge(1), rank=2, unique="sample").value.iloc[0] == 2
    assert pd.isna(summarize(d, support=d.value.ge(1), rank=3, unique="sample").value.iloc[0])
    with pytest.raises(ValueError, match="inflate recurrence"):
        summarize(pd.concat([d, d.iloc[[0]]]), support=True, rank=2, unique="sample")

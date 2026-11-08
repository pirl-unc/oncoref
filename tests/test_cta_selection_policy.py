"""Single RNA boundary and positive, evidence-backed Other definitions."""

import pandas as pd
import pytest

from oncoref import cta_selection_policy as policy
from oncoref.cta_funnel import EvidenceLayer


def context(tier="vital", value=4, grade="Missing"):
    return {
        "gene_id": "g",
        "source": "hpa_single_cell_tissue",
        "version": "v23",
        "taxonomy_id": 9606,
        "stage": "adult",
        "state": "normal",
        "tissue": "colon",
        "cell_type": "cell",
        "category": tier,
        "mapping_status": "mapped",
        "negative_ihc_grade": grade,
        "value": value,
        "unit": "nTPM",
    }


@pytest.mark.parametrize("tier,ceiling", [("vital", 4), ("protected", 8), ("manageable", 16)])
def test_one_boundary_equality_and_anatomical_risk(tier, ceiling):
    rows = pd.DataFrame([context(tier, ceiling), context(tier, ceiling + 0.01)])
    tissue = {"vital": "heart", "protected": "colon", "manageable": "thyroid gland"}[tier]
    rows["tissue"] = tissue
    assessed = policy.assess_contexts(rows)
    assert assessed.candidate_exclusion.tolist() == [False, True]
    assert assessed.unsafe_evidence.tolist() == [False, tier != "manageable"]
    assert "strong_boundary" not in assessed
    assert (
        policy.assess_contexts(
            pd.DataFrame([{**context(tier, 0, "Enhanced"), "tissue": tissue}])
        ).candidate_ceiling.iloc[0]
        == ceiling * 8
    )


def test_resolved_cells_are_screened_even_in_a_disposable_organ():
    rows = pd.DataFrame(
        [
            context("disposable", 17) | {"tissue": "tonsil", "cell_type": "B cells"},
            context("disposable", 9) | {"tissue": "tonsil", "cell_type": "T cells"},
        ]
    )
    assessed = policy.assess_contexts(rows)
    assert assessed.selection_tier.tolist() == ["manageable", "protected"]
    assert assessed.candidate_exclusion.tolist() == [True, True]
    assert assessed.unsafe_evidence.tolist() == [False, True]


def test_compartment_match_does_not_transfer_to_other_tissue_state_or_species():
    crosswalk = pd.DataFrame(
        [
            {
                "source": "hpa_single_cell_tissue",
                "version": "v23",
                "tissue": "colon",
                "cell_type": "distal enterocytes",
                "ihc_cell_type": "glandular cells",
                "ihc_match_source": "citation",
            }
        ]
    )
    base = context() | {"cell_type": "distal enterocytes"}
    rows = pd.DataFrame(
        [
            base,
            base | {"tissue": "heart muscle"},
            base | {"state": "stimulated"},
            base | {"taxonomy_id": 9615},
            base | {"stage": "fetal"},
            base | {"cell_type": "T-cells"},
        ]
    )
    ihc = pd.DataFrame(
        [
            {
                "Gene": "g",
                "Tissue": "colon",
                "Cell type": "glandular cells",
                "Level": "Not detected",
                "Reliability": "Enhanced",
            }
        ]
    )
    result = policy.match_ihc(rows, ihc, crosswalk=crosswalk, positive_controls={"g"})
    assert result.negative_ihc_grade.tolist() == ["Enhanced", *["Missing/uncertain"] * 5]
    assert result.ihc_match_scope.tolist()[0] == "reviewed_compartment"
    assert not result.ihc_exact_match.any()
    assert result.ihc_match_scope.iloc[-1] == "same_tissue_unresolved_cell"
    conflict = policy.match_ihc(
        rows.iloc[:1], ihc, crosswalk=crosswalk, positive_controls={"g"}, conflicting_genes={"g"}
    )
    assert conflict.negative_ihc_grade.iloc[0] == "Missing/uncertain"


def profiles():
    return pd.DataFrame(
        [
            {
                "entity": str(i),
                "required_coverage": True,
                "unresolved_evidence": False,
                "tumor_measured": True,
                "tumor_profile_complete": True,
                "developmental_measured": True,
                "normal_maximum_floor_ratio": 0.5,
                "reproductive_ntpm": 0.5,
                "developmental_all_lineage_max_tpm": 0.5,
                "tumor_maximum_floor_ratio": 0.5,
                "any_detected_protein": False,
                "positive_context": False,
                "vital_protected_exclusion": False,
                "manageable_tissue_count": 0,
            }
            for i in range(5)
        ]
    )


def test_every_protein_resolves_and_normal_holds_only_block_cta():
    facts = profiles()
    facts.loc[1, "required_coverage"] = False
    facts.loc[2, "tumor_maximum_floor_ratio"] = 2
    facts.loc[3, ["positive_context", "normal_maximum_floor_ratio"]] = [True, 20]
    facts.loc[4, ["normal_maximum_floor_ratio", "manageable_tissue_count"]] = [20, 2]
    facts = policy.other_profiles(facts)
    yes, no = pd.Series(True, index=facts.index), pd.Series(False, index=facts.index)
    result = policy.classify(
        facts,
        [
            EvidenceLayer("Reproductive / embryonic association", yes, no),
            EvidenceLayer("Normal bulk", pd.Series([False, False, False, True, True]), no),
        ],
        entity_column="entity",
    )
    assert result.classifications.category.tolist() == ["Other"] * 5
    assert result.classifications.subtype.tolist() == ["Nonreproductive"] * 5
    assert result.classifications.assessment_status.tolist() == [
        "screened",
        "incomplete_evidence",
        "screened",
        "normal_evidence_hold",
        "normal_evidence_hold",
    ]
    assert (
        result.funnel[["CTA_eligible_so_far", "Unsafe", "Other", "Unassigned"]]
        .sum(axis=1)
        .eq(5)
        .all()
    )


def test_unsafe_wins_regardless_of_order_and_cannot_be_rescued_by_neither():
    facts = policy.other_profiles(profiles())
    yes, no = pd.Series(True, index=facts.index), pd.Series(False, index=facts.index)
    layers = [EvidenceLayer("gate", yes, no), EvidenceLayer("vital", yes, yes)]
    first = policy.classify(facts, layers, entity_column="entity").classifications
    second = policy.classify(facts, list(reversed(layers)), entity_column="entity").classifications
    assert first.category.eq("Unsafe").all()
    pd.testing.assert_series_equal(first.category, second.category)


def test_missing_tumor_or_developmental_evidence_does_not_revoke_normal_low_profile():
    facts = profiles().iloc[:1].copy()
    facts["tumor_profile_complete"] = False
    facts["developmental_measured"] = False
    facts["developmental_all_lineage_max_tpm"] = float("nan")
    result = policy.other_profiles(facts)
    assert result.tumor_measured.all()
    assert result.other_low_expression.all()


@pytest.mark.parametrize("name", sorted(policy.CANDIDATE_ONLY_LAYERS))
def test_admission_failures_cannot_be_wired_as_unsafe(name):
    facts = policy.other_profiles(profiles())
    yes, no = pd.Series(True, index=facts.index), pd.Series(False, index=facts.index)
    result = policy.classify(
        facts, [EvidenceLayer(name, yes, no)], entity_column="entity"
    ).classifications
    assert not result.category.eq("Unsafe").any()
    with pytest.raises(ValueError, match="admission requirement, not safety evidence"):
        policy.classify(facts, [EvidenceLayer(name, yes, yes)], entity_column="entity")


def test_normal_safety_is_evaluated_even_after_cancer_admission_fails():
    facts = policy.other_profiles(profiles())
    yes, no = pd.Series(True, index=facts.index), pd.Series(False, index=facts.index)
    layers = [
        EvidenceLayer("Tumor sample ≥2 TPM/pTPM", yes, no),
        EvidenceLayer("Normal bulk", yes, yes),
    ]
    result = policy.classify(facts, layers, entity_column="entity").classifications
    assert result.category.eq("Unsafe").all()
    assert result.Unsafe_reasons.eq("Normal bulk").all()


def test_association_does_not_depend_on_normal_expression_or_coverage():
    facts = pd.DataFrame(
        {
            "reproductive_ntpm": [5, 5, 4, 0, 100, 0],
            "developmental_mean_tpm": [0, 0, 0, 5, 0, 5],
            "healthy_maximum_ntpm": [4, 5, 0, 100, float("nan"), 0],
            "healthy_maximum_tpm": [100, 0, 0, 4, 0, float("nan")],
        }
    )
    result = policy.association_evidence(facts)
    assert result.positive_context.tolist() == [True, True, False, True, True, True]
    assert result.reproductive_route.tolist() == [True, True, False, False, True, False]
    assert result.developmental_route.tolist() == [False, False, False, True, False, True]
    assert "adult_restriction" not in result


def test_association_failure_does_not_make_expression_unsafe():
    facts = policy.other_profiles(profiles())
    facts["other_low_expression"] = False
    yes, no = pd.Series(True, index=facts.index), pd.Series(False, index=facts.index)
    result = policy.classify(
        facts,
        [EvidenceLayer("Reproductive / embryonic association", yes, no)],
        entity_column="entity",
    ).classifications
    assert result.category.eq("Other").all()
    assert result.Unsafe_reasons.eq("").all()


def test_reproductive_ihc_route_requires_quality_anatomy_identity_and_all_members():
    facts = pd.DataFrame(
        {
            "member_gene_ids": [
                '["a"]',
                '["b"]',
                '["c"]',
                '["d"]',
                '["e", "f"]',
                '["g", "h"]',
                '["i"]',
            ],
            "Antibody concern": [False, False, False, False, False, False, True],
            "reproductive_ntpm": 0.1,
            "developmental_mean_tpm": 0,
        }
    )
    base = {
        "source": "hpa_normal_tissue",
        "version": "v23",
        "taxonomy_id": 9606,
        "mapping_status": "mapped",
        "category": "reproductive",
        "Level": "Low",
        "Reliability": "Approved",
    }

    observations = pd.DataFrame(
        [
            base | {"gene_id": "a"},
            base | {"gene_id": "b", "Reliability": "Uncertain"},
            base | {"gene_id": "c", "category": "disposable"},
            base | {"gene_id": "d", "taxonomy_id": 10090},
            base | {"gene_id": "e"},
            base | {"gene_id": "g"},
            base | {"gene_id": "h"},
            base | {"gene_id": "i"},
        ]
    )
    result, evidence = policy.reproductive_protein_evidence(facts, observations)
    result = policy.association_evidence(result)
    assert result.positive_context.tolist() == [True, False, False, False, False, True, False]
    assert result.reproductive_rna_association.eq(False).all()
    assert not evidence.loc[evidence.gene_id.eq("b"), "qualifying_reproductive_positive"].any()


def test_tumor_restricted_support_has_no_reproductive_requirement_and_no_safety_override():
    facts = profiles()
    facts["tumor_samples_ge10"] = [2, 1, float("nan"), 2, 2]
    facts["tumor_identity_resolved"] = [True, True, False, True, True]
    facts["reliable_tumor_protein"] = [False, False, True, False, False]
    facts["normal_restriction_reviewed"] = True
    facts.loc[3, "normal_maximum_floor_ratio"] = 1  # Strictly below 2, not <=2.
    facts = policy.tumor_restricted_evidence(policy.other_profiles(facts))
    yes, no = pd.Series(True, index=facts.index), pd.Series(False, index=facts.index)
    normal = pd.Series([False, False, False, False, True])
    result = policy.classify(
        facts,
        [
            EvidenceLayer("Reproductive / embryonic association", yes, no),
            EvidenceLayer("Normal bulk", normal, normal),
        ],
        entity_column="entity",
    )
    assert result.classifications.category.tolist() == [
        "Other",
        "Other",
        "Other",
        "Other",
        "Unsafe",
    ]
    assert (
        result.funnel[["CTA_eligible_so_far", "Other", "Unsafe", "Unassigned"]]
        .sum(axis=1)
        .eq(5)
        .all()
    )


def test_numerical_tumor_restriction_without_normal_context_review_remains_neither():
    facts = profiles().iloc[:1].copy()
    facts["tumor_samples_ge10"] = 2
    facts["tumor_identity_resolved"] = True
    facts["normal_restriction_reviewed"] = False
    facts["reliable_tumor_protein"] = False
    facts = policy.tumor_restricted_evidence(policy.other_profiles(facts))
    assert facts.tumor_restricted_profile.all()
    assert facts.tumor_restricted_lead.all()
    assert not facts.tumor_restricted_identity.any()
    result = policy.classify(
        facts,
        [
            EvidenceLayer(
                "Reproductive / embryonic association", pd.Series([True]), pd.Series([False])
            )
        ],
        entity_column="entity",
    )
    assert result.classifications.category.tolist() == ["Other"]


def test_diffuseness_is_route_specific_and_missingness_is_not_failure():
    facts = policy.association_evidence(
        pd.DataFrame(
            {
                "protein_equivalence_group_id": list("abcdef"),
                "reproductive_ntpm": [5, 5, 1, 5, 5, 5],
                "developmental_mean_tpm": [0, 10, 10, 0, 0, 0],
            }
        )
    )
    scores = pd.DataFrame(
        {
            "entity_id": list("abcdef"),
            "source": "hpa_rna_consensus",
            "version": "v23",
            "unit": "nTPM",
            "taxonomy_id": 9606,
            "summary": "consensus",
            "diffuseness": [0.2, 0.2, 0.9, float("nan"), 0.1, 0.2],
            "reference_value": [5, 5, 1, 5, 5, 5],
            "reference_complete": [True, True, True, False, True, True],
            "normal_coverage": [1, 1, 1, 0.5, 1, 0.9],
        }
    )
    result = policy.concentration_evidence(facts, scores)
    assert result.reproductive_concentration_hold.tolist() == [
        True,
        False,
        False,
        False,
        False,
        False,
    ]
    assert result.adult_profile_status.tolist() == [
        "diffuse",
        "diffuse",
        "no_qualifying_adult_reference",
        "unassessed",
        "concentrated",
        "unassessed",
    ]
    assert result.reproductive_route.tolist() == [False, False, False, True, True, True]
    with pytest.raises(ValueError, match="HPA v23"):
        policy.concentration_evidence(facts, scores.assign(source="gtex_tissue_rna"))


def test_diffuse_other_discloses_holds_and_incomplete_coverage_without_unsafe():
    facts = policy.other_profiles(profiles())
    facts[["other_low_expression", "other_disposable"]] = False
    facts["positive_context"] = True
    facts.loc[3, "required_coverage"] = False
    facts.loc[4, "unresolved_evidence"] = True
    yes, no = pd.Series(True, index=facts.index), pd.Series(False, index=facts.index)
    layers = [
        EvidenceLayer("Reproductive concentration", yes, no),
        EvidenceLayer(
            "Normal bulk",
            pd.Series([False, True, True, False, False]),
            pd.Series([False, False, True, False, False]),
        ),
    ]
    first = policy.classify(facts, layers, entity_column="entity").classifications
    second = policy.classify(facts, list(reversed(layers)), entity_column="entity").classifications
    assert first.category.tolist() == ["Other", "Other", "Unsafe", "Other", "Other"]
    assert first.subtype.tolist() == ["Diffuse", "Diffuse", "", "Diffuse", "Diffuse"]
    assert first.assessment_status.tolist() == [
        "screened",
        "normal_evidence_hold",
        "screened",
        "incomplete_evidence",
        "incomplete_evidence",
    ]
    pd.testing.assert_series_equal(first.category, second.category)


def test_other_subtypes_are_exhaustive_and_normal_holds_never_become_unsafe():
    assert set(policy.OTHER_SUBTYPE_RULES) == {
        "Nonreproductive",
        "Reproductive",
        "Diffuse",
        "Held",
        "Contained",
        "Nested",
        "Translation",
    }
    names = ["Nonreproductive", "Reproductive", "Diffuse"]
    facts = pd.concat([profiles().iloc[:1]] * len(names), ignore_index=True)
    facts["entity"] = names
    facts["positive_context"] = [False, True, True]
    yes = pd.Series(True, index=facts.index)
    no = pd.Series(False, index=facts.index)
    layers = [
        EvidenceLayer("Reproductive / embryonic association", pd.Series([True, False, False]), no),
        EvidenceLayer("Cancer expression", pd.Series([True, True, False]), no),
        EvidenceLayer("Reproductive concentration", pd.Series([True, True, True]), no),
        EvidenceLayer("Normal", no, no),
    ]
    result = policy.classify(facts, layers, entity_column="entity").classifications
    assert result.subtype.tolist() == names
    assert result.classification.tolist() == ["Other:" + name for name in names]
    # Target gates pass but a normal-tissue hold blocks CTA: Other:Held, not Unsafe.
    passing = facts.iloc[:1].copy()
    passing["positive_context"] = True
    one = pd.Series([True])
    gates = [
        EvidenceLayer(
            "Reproductive / embryonic association", pd.Series([False]), pd.Series([False])
        ),
        EvidenceLayer("Cancer expression", pd.Series([False]), pd.Series([False])),
        EvidenceLayer("Reproductive concentration", pd.Series([False]), pd.Series([False])),
    ]
    held = policy.classify(
        passing, [*gates, EvidenceLayer("Normal", one, ~one)], entity_column="entity"
    ).classifications
    assert held.classification.tolist() == ["Other:Held"]
    assert held.assessment_status.tolist() == ["normal_evidence_hold"]
    clear = policy.classify(
        passing, [*gates, EvidenceLayer("Normal", ~one, ~one)], entity_column="entity"
    ).classifications
    assert clear.classification.tolist() == ["CTA"]
    excluded = policy.classify(
        facts, [*layers[:-1], EvidenceLayer("Normal", yes, yes)], entity_column="entity"
    ).classifications
    assert excluded.classification.eq("Unsafe").all()
    assert excluded.subtype.eq("").all()


def test_adult_disposable_pattern_is_independent_of_developmental_expression():
    facts = profiles()
    facts["disposable_maximum_floor_ratio"] = 10.0
    facts["outside_disposable_floor_ratio"] = 0.2
    facts["protein_outside_disposable"] = False
    facts["developmental_measured"] = [False, True, True, True, False]
    facts["developmental_all_lineage_max_tpm"] = [float("nan"), 0, 1000, 1, float("nan")]
    result = policy.other_profiles(facts)
    assert result.other_disposable.all()
    facts.loc[3, "outside_disposable_floor_ratio"] = 1.0
    facts.loc[4, "protein_outside_disposable"] = True
    assert policy.other_profiles(facts).other_disposable.tolist() == [
        True,
        True,
        True,
        False,
        False,
    ]


def _passing_facts(n):
    facts = pd.concat([profiles().iloc[:1]] * n, ignore_index=True)
    facts["entity"] = [f"e{i}" for i in range(n)]
    facts["positive_context"] = True
    return facts


def _gates(n):
    no = pd.Series([False] * n)
    return [
        EvidenceLayer("Reproductive / embryonic association", no, no),
        EvidenceLayer("Cancer expression", no, no),
        EvidenceLayer("Reproductive concentration", no, no),
        EvidenceLayer("Normal", no, pd.Series([False] * (n - 1) + [True])),
    ]


def test_cta_mostly_inside_unsafe_sequence_is_contained_not_unsafe():
    facts = policy.other_profiles(_passing_facts(5))
    # e4 is Unsafe; overlap is the fraction of each entity's 9-mers inside Unsafe proteins.
    fraction = pd.Series([0.0, 0.94, 0.95, float("nan"), 1.0], index=facts.index)
    seen = []

    def overlap(unsafe):
        seen.append(unsafe.tolist())
        return fraction

    result = policy.classify(
        facts, _gates(5), entity_column="entity", unsafe_overlap=overlap
    ).classifications
    assert result.classification.tolist() == [
        "CTA",
        "CTA",
        "Other:Contained",  # the 95% boundary is inclusive
        "CTA",  # undefined overlap (no 9-mer) is never contained
        "Unsafe",  # Unsafe is untouched by its own overlap
    ]
    assert seen[0] == [False, False, False, False, True]
    assert result.subtype.tolist() == ["", "", "Contained", "", ""]
    assert result.assessment_status.tolist() == ["screened"] * 5
    assert result.unsafe_sequence_overlap_pct.tolist()[:3] == [0.0, 94.0, 95.0]


def test_containment_adds_a_funnel_stage_and_an_independent_effect():
    facts = policy.other_profiles(_passing_facts(5))
    fraction = pd.Series([0.0, 0.96, 0.95, 0.1, 1.0], index=facts.index)
    result = policy.classify(
        facts, _gates(5), entity_column="entity", unsafe_overlap=lambda unsafe: fraction
    )
    stage = result.funnel[result.funnel.layer.eq("Unsafe sequence containment")].iloc[0]
    assert (stage.new_candidate_removals, stage.CTA_eligible_so_far) == (2, 2)
    effect = result.independent_effects[
        result.independent_effects.layer.eq("Unsafe sequence containment")
    ]
    assert effect.CTA_added_if_layer_removed.iloc[0] == 2
    # The final assignment row reconciles with the per-protein categories.
    final = result.funnel.iloc[-1]
    assert (final.CTA_eligible_so_far, final.Unsafe, final.Other) == (2, 1, 2)


def test_without_an_overlap_input_classification_is_unchanged():
    facts = policy.other_profiles(_passing_facts(5))
    plain = policy.classify(facts, _gates(5), entity_column="entity")
    assert plain.classifications.classification.tolist() == ["CTA"] * 4 + ["Unsafe"]
    assert plain.classifications.unsafe_sequence_overlap_pct.isna().all()
    assert "Unsafe sequence containment" not in set(plain.funnel.layer)


def test_overlap_must_align_with_the_evidence_index():
    facts = policy.other_profiles(_passing_facts(5))
    with pytest.raises(ValueError, match="align"):
        policy.classify(
            facts,
            _gates(5),
            entity_column="entity",
            unsafe_overlap=lambda unsafe: pd.Series([0.0] * 5, index=range(10, 15)),
        )


def test_ninemer_overlap_counts_only_unsafe_sequence_and_is_cached():
    import numpy as np

    from oncoref.cta_sequence_overlap import unsafe_overlap

    groups = [frozenset({"a", "b", "c", "d"}), frozenset({"a", "b", "c", "x"}), frozenset()]
    fn = unsafe_overlap(groups, pd.RangeIndex(3))
    fraction = fn(pd.Series([True, False, False]))
    assert fraction.iloc[1] == pytest.approx(0.75)
    assert np.isnan(fraction.iloc[2])
    # Mutating a returned series cannot corrupt later calls.
    fraction.iloc[1] = 0
    assert fn(pd.Series([True, False, False])).iloc[1] == pytest.approx(0.75)
    assert fn(pd.Series([False, False, False])).iloc[1] == 0


def test_nested_microprotein_is_other_nested_never_cta():
    facts = policy.other_profiles(_passing_facts(3))
    facts["nested_microprotein"] = [False, True, False]
    layers = _gates(3)
    layers[0] = EvidenceLayer(
        "Reproductive / embryonic association", facts.nested_microprotein, layers[0][2]
    )
    result = policy.classify(facts, layers, entity_column="entity").classifications
    assert result.classification.tolist()[:2] == ["CTA", "Other:Nested"]


def test_translation_failure_blocks_admission_without_erasing_unsafe():
    facts = pd.concat([profiles().iloc[:1]] * 2, ignore_index=True)
    facts["entity"] = ["held", "unsafe"]
    facts["positive_context"] = True
    no = pd.Series(False, index=facts.index)
    yes = ~no
    layers = [
        EvidenceLayer("Translation qualification", yes, no),
        EvidenceLayer("Normal", no, pd.Series([False, True])),
    ]
    result = policy.classify(facts, layers, entity_column="entity").classifications
    assert result.classification.tolist() == ["Other:Translation", "Unsafe"]
    assert result.assessment_status.tolist() == ["translation_review", "screened"]
    with pytest.raises(ValueError, match="admission requirement"):
        policy.classify(
            facts, [EvidenceLayer("Translation qualification", yes, yes)], entity_column="entity"
        )


def test_translation_admission_requires_support_and_exact_record_qualification():
    facts = pd.DataFrame(
        {
            "translation_support_eligible": [True, False, True, False],
            "translation_qualified": [True, True, False, False],
            "representative_selection_approved": [True, True, True, True],
        }
    )
    assert policy.translation_admission(facts).tolist() == [True, False, False, False]
    assert not policy.translation_admission(
        facts.assign(representative_selection_approved=False)
    ).any()
    with pytest.raises(ValueError, match="boolean translation evidence"):
        policy.translation_admission(facts.assign(translation_support_eligible=None))

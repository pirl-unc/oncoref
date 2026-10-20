# Licensed under the Apache License, Version 2.0.
from dataclasses import replace

import pytest

from oncoref.normal_tissue_gene_filters import (
    NormalTissueCoverageRequirement,
    NormalTissueEvidence,
    evaluate_normal_tissue_gene_filter,
    normal_tissue_gene_filter_profiles,
)
from oncoref.normal_tissues import resolve_normal_tissue

PREGNANCY = "male_reproductive_maternal_fetal_excluded"
MALE_PRESERVING = "female_developmental_male_reproductive_excluded"


def evidence(tissue, *, stage="adult", sex="male", result="detected", species="human", **kwargs):
    a = resolve_normal_tissue(tissue, species=species, stage=stage, sex=sex, **kwargs)
    return NormalTissueEvidence(
        "ENSEMBL:example",
        a,
        f"{tissue}:{stage}:{sex}",
        "RNA",
        "test-detection-v1",
        result,
        f"test:{tissue}:{stage}:{sex}",
    )


def panel(profile):
    # Synthetic minimal coverage exercises logic; not a biological coverage recommendation.
    specs = (
        [
            ("heart", "adult", "female", None),
            ("uterus", "adult", "female", "pregnant"),
            ("embryo", "embryonic", None, None),
            ("testis", "fetal", "male", None),
            ("ovary", "fetal", "female", None),
            ("placenta", None, None, None),
            ("amnion", None, None, None),
        ]
        if profile == PREGNANCY
        else [
            ("testis", "adult", "male", None),
            ("heart", "adult", "male", None),
            ("mammary gland", "adult", "male", None),
        ]
    )
    negatives = [
        evidence(t, stage=s, sex=sex, result="not_detected", physiological_state=state)
        for t, s, sex, state in specs
    ]
    required = [
        NormalTissueCoverageRequirement(x.annotation, x.context_id, x.assay, x.detection_policy_id)
        for x in negatives
    ]
    return negatives, required


@pytest.mark.parametrize(
    "profile,tissue,sex",
    [
        (PREGNANCY, "testis", "male"),
        (MALE_PRESERVING, "ovary", "female"),
    ],
)
def test_positive_only_cannot_pass_a_patient_context_filter(profile, tissue, sex):
    positive = evidence(tissue, sex=sex)
    result = evaluate_normal_tissue_gene_filter([positive], profile=profile)
    assert result["positive_pattern_compatible"]
    assert result["status"] == "insufficient_evidence"
    assert result["clinical_safety"] == "not_assessed"
    negatives, required = panel(profile)
    result = evaluate_normal_tissue_gene_filter(
        [positive, *negatives],
        profile=profile,
        required_contexts=required,
        coverage_panel_id="synthetic-test-panel-v1",
    )
    assert result["passes_declared_panel"]
    assert result["clinical_safety"] == "not_assessed"


@pytest.mark.parametrize(
    "tissue,stage,sex",
    [
        ("testis", "fetal", "male"),
        ("ovary", "fetal", "female"),
        ("placenta", None, None),
        ("uterus", "adult", "female"),
        ("hair follicle", "adult", "female"),
    ],
)
def test_pregnancy_profile_has_no_developmental_or_disposable_exception(tissue, stage, sex):
    rows = [evidence("testis"), evidence(tissue, stage=stage, sex=sex)]
    result = evaluate_normal_tissue_gene_filter(rows, profile=PREGNANCY)
    assert result["status"] == "excluded"
    assert result["excluded_positive_indices"] == (1,)


@pytest.mark.parametrize(
    "tissue,stage,sex",
    [
        ("testis", "adult", "male"),
        ("testis", "fetal", "male"),
        ("prostate", "adult", "male"),
        ("mammary gland", "adult", "male"),
    ],
)
def test_male_preservation_vetoes_male_reproductive_and_known_male_mammary(tissue, stage, sex):
    result = evaluate_normal_tissue_gene_filter(
        [evidence("ovary", sex="female"), evidence(tissue, stage=stage, sex=sex)],
        profile=MALE_PRESERVING,
    )
    assert result["status"] == "excluded"


def test_female_and_fetal_pattern_is_a_union_with_explicit_male_anatomy_veto():
    result = evaluate_normal_tissue_gene_filter(
        [evidence("ovary", sex="female"), evidence("liver", stage="fetal", sex=None)],
        profile=MALE_PRESERVING,
    )
    assert result["positive_pattern_compatible"]
    assert not result["passes_declared_panel"]


@pytest.mark.parametrize(
    "change",
    [
        {"context_id": "another-sample"},
        {"assay": "IHC"},
        {"detection_policy_id": "another-threshold"},
        {"result": "unassessed"},
    ],
)
def test_negative_calls_cannot_clear_another_sample_assay_or_rule(change):
    negatives, required = panel(PREGNANCY)
    negatives[0] = replace(negatives[0], **change)
    result = evaluate_normal_tissue_gene_filter(
        [evidence("testis"), *negatives],
        profile=PREGNANCY,
        required_contexts=required,
        coverage_panel_id="synthetic-test-panel-v1",
    )
    assert not result["passes_declared_panel"]
    assert 0 in result["missing_required_context_indices"]


def test_positive_in_any_source_overrides_a_complete_negative_panel():
    negatives, required = panel(PREGNANCY)
    result = evaluate_normal_tissue_gene_filter(
        [evidence("testis"), *negatives, evidence("heart", sex="female")],
        profile=PREGNANCY,
        required_contexts=required,
        coverage_panel_id="test",
    )
    assert result["status"] == "excluded"


def test_unknown_stage_or_sex_never_establishes_postnatal_sex_restriction():
    for row in (
        evidence("testis", stage=None),
        evidence("testis", sex=None),
        evidence("testis", stage="new unreviewed stage"),
    ):
        result = evaluate_normal_tissue_gene_filter([row], profile=PREGNANCY)
        assert not result["positive_pattern_compatible"]
        assert result["status"] == "insufficient_evidence"


def test_species_genes_and_taxonomy_are_not_pooled():
    with pytest.raises(ValueError, match="different species"):
        evaluate_normal_tissue_gene_filter(
            [evidence("testis"), evidence("testis", species="dog")],
            profile=PREGNANCY,
        )
    with pytest.raises(ValueError, match="one gene"):
        evaluate_normal_tissue_gene_filter(
            [evidence("testis"), replace(evidence("testis"), gene_id="ENSEMBL:other")],
            profile=PREGNANCY,
        )
    row = evidence("testis")
    with pytest.raises(ValueError, match="current taxonomy"):
        evaluate_normal_tissue_gene_filter(
            [replace(row, annotation=replace(row.annotation, taxonomy_id="old"))],
            profile=PREGNANCY,
        )
    profiles = normal_tissue_gene_filter_profiles()
    profiles[PREGNANCY]["label"] = "changed"
    assert normal_tissue_gene_filter_profiles()[PREGNANCY]["label"] != "changed"


def test_negative_source_version_and_sex_do_not_substitute_for_required_context():
    negatives, required = panel(PREGNANCY)
    for changed in (
        replace(negatives[0].annotation, source_version="another-version"),
        replace(negatives[0].annotation, sex="male"),
    ):
        rows = [replace(negatives[0], annotation=changed), *negatives[1:]]
        result = evaluate_normal_tissue_gene_filter(
            [evidence("testis"), *rows],
            profile=PREGNANCY,
            required_contexts=required,
            coverage_panel_id="test",
        )
        assert not result["passes_declared_panel"]
        assert 0 in result["missing_required_context_indices"]


def test_coverage_panel_fingerprint_changes_with_source_context_not_list_order():
    negatives, required = panel(PREGNANCY)

    def evaluate(contexts):
        return evaluate_normal_tissue_gene_filter(
            [evidence("testis"), *negatives],
            profile=PREGNANCY,
            required_contexts=contexts,
            coverage_panel_id="test",
        )["coverage_panel_sha256"]

    original = evaluate(required)
    assert evaluate(reversed(required)) == original
    altered = [replace(required[0], context_id="another-sample"), *required[1:]]
    assert evaluate(altered) != original

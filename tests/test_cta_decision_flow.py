"""One classification vocabulary without losing normal holds or member identity."""

import json
from pathlib import Path

import pandas as pd
import pytest

from oncoref import cta_decision_flow as flow
from oncoref import cta_selection

ROOT = Path(__file__).resolve().parents[1]


def test_uniform_view_preserves_all_categories_loci_and_evidence():
    released = cta_selection.results()
    result = cta_selection.decisions()
    assert len(result) == 19893
    assert sum(map(len, result.gene_ids)) == 20099
    assert result.category.value_counts().to_dict() == {"Unsafe": 18259, "Other": 1027, "CTA": 607}
    assert set(result[result.category.eq("Other")].subtype) == {"Review", "Expression", "Sequence"}
    for column in [
        "entity_id",
        "category",
        "gene_ids",
        "reasons",
        "assessment_status",
        "evidence_flags",
    ]:
        pd.testing.assert_series_equal(result[column], released[column])
    pd.testing.assert_series_equal(
        result.source_classification, released.classification, check_names=False
    )
    pd.testing.assert_frame_equal(flow.apply(result), result)
    assert result[result.category.eq("CTA")].subtype.eq("").all()


def test_review_takes_precedence_over_expression_and_sequence():
    rows = pd.DataFrame(
        [
            {
                "category": "Other",
                "subtype": label,
                "classification": "Other:" + label,
                "assessment_status": status,
            }
            for label, status in [
                ("Nonreproductive", "normal_evidence_hold"),
                ("Nested", "translation_review"),
                ("Diffuse", "incomplete_evidence"),
                ("Nested", "screened"),
                ("Reproductive", "screened"),
            ]
        ]
    )
    result = flow.apply(rows)
    assert result.subtype.tolist() == ["Review", "Review", "Review", "Sequence", "Expression"]
    assert rows.subtype.iloc[0] == "Nonreproductive"  # no mutation
    for column, value, message in [
        ("category", "Discovery", "assessed"),
        ("subtype", "Unknown", "subtype"),
        ("category", "CTA", "completed"),
    ]:
        invalid = rows.copy()
        invalid.loc[0, column] = value
        with pytest.raises(ValueError, match=message):
            flow.apply(invalid)


def test_api_filters_only_standard_other_reasons():
    review = cta_selection.decisions(subtype="Review")
    assert not review.empty and review.classification.eq("Other:Review").all()
    with pytest.raises(ValueError, match="subcategories"):
        cta_selection.decisions(category="CTA", subtype="Review")
    with pytest.raises(ValueError, match="subcategories"):
        cta_selection.decisions(subtype="Translation")


def test_ssu72_coding_rule_and_all_inventory_exports():
    audit = ROOT / "docs/audits/cta-candidate-evidence-20261108"
    proof = pd.read_csv(audit / "ssu72-admission.csv")
    assert set(proof.protein_group_label) == {f"SSU72L{i}" for i in range(1, 6)}
    assert proof.classification.eq("CTA").all() and proof.length.eq(194).all()
    assert proof.developmental_mean_tpm.ge(5).all()
    assert (
        proof[
            [
                "cancer_association",
                "required_coverage",
                "translation_support_eligible",
                "translation_qualified",
                "representative_selection_approved",
            ]
        ]
        .all()
        .all()
    )
    assert proof.exclusion_reasons.isna().all()
    proteins = pd.read_csv(audit / "uniform-protein-decisions.csv.gz")
    genes = pd.read_csv(audit / "uniform-gene-decisions.csv.gz")
    assert len(proteins) == 19893 and len(genes) == 20099 and genes.gene_id.is_unique
    assert proteins.category.value_counts().to_dict() == {
        "Unsafe": 18259,
        "Other": 1034,
        "CTA": 600,
    }
    assert proteins.classification.isin(
        ["CTA", "Unsafe", "Other:Review", "Other:Sequence", "Other:Expression"]
    ).all()
    assert json.loads((audit / "decision-flow.json").read_text()) == flow.definition()

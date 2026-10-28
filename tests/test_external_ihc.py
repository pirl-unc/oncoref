import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from oncoref.external_ihc import assess_evidence, reconcile_hpa


def observation(**updates):
    return {
        "evidence_id": "independent-normal-1",
        "gene_id": "ENSG00000164458",
        "tissue": "lung",
        "cell_type": "macrophages",
        "stage": "adult",
        "specimen_state": "normal",
        "taxonomy_id": 9606,
        "antibody_id": "fixture-clone",
        "reagent_lineage": "fixture-independent-lineage",
        "source_url": "https://example.org/verified-ihc",
        "source_locator": "Table 2",
        "source_kind": "peer_reviewed",
        "assay": "IHC",
        "result": "not_detected",
        "review_status": "reviewed",
        "positive_control": True,
        "specificity_validation": "genetic_ihc",
        "validation_locator": "Figure 1",
        "exact_cell_coverage": True,
        "normal_donors": 3,
        "independent_of_hpa": True,
        "notes": "Synthetic evidence for a policy test, not biological data.",
        **updates,
    }


def hpa(**updates):
    return pd.DataFrame(
        [
            {
                "gene_id": "ENSG00000164458",
                "tissue": "Lung",
                "cell_type": "Macrophages",
                "stage": "adult",
                "specimen_state": "normal",
                "taxonomy_id": 9606,
                "Antibody": "HPA003322",
                "Level": "Medium",
                "Reliability": "Uncertain",
                **updates,
            }
        ]
    )


def test_context_validated_negative_can_supersede_uncertain_without_erasing_source():
    result = reconcile_hpa(hpa(), pd.DataFrame([observation()]))
    assert result.external_decision.iloc[0] == "uncertain_superseded_by_negative"
    assert result.Level.iloc[0] == "Medium" and result.Reliability.iloc[0] == "Uncertain"
    assert "independent-normal-1" in result.external_evidence_ids.iloc[0]


@pytest.mark.parametrize(
    "changes",
    [
        {"cell_type": "unspecified", "exact_cell_coverage": False},
        {"specificity_validation": "genetic_western_blot"},
        {"normal_donors": None},
        {"normal_donors": 1},
        {"stage": "fetal"},
        {"specimen_state": "tumor_adjacent"},
        {"taxonomy_id": 10090},
        {"independent_of_hpa": False},
        {"positive_control": False},
        {"review_status": "pending"},
        {"validation_locator": ""},
        {"tissue": "liver"},
        {"assay": "WB"},
    ],
)
def test_incomplete_or_mismatched_evidence_cannot_clear_uncertain(changes):
    result = reconcile_hpa(hpa(), pd.DataFrame([observation(**changes)]))
    assert result.external_decision.iloc[0] == "unresolved"


def test_reliable_positives_and_validated_disagreement_are_retained():
    e = pd.DataFrame([observation()])
    assert (
        reconcile_hpa(hpa(Reliability="Supported"), e).external_decision.iloc[0]
        == "conflicting_validated_evidence"
    )
    positive = observation(
        evidence_id="second",
        antibody_id="other-clone",
        reagent_lineage="other-lineage",
        result="detected",
        normal_donors=1,
    )
    assert (
        reconcile_hpa(hpa(), pd.DataFrame([positive])).external_decision.iloc[0]
        == "normal_expression_supported"
    )
    assert (
        reconcile_hpa(hpa(), pd.DataFrame([observation(), positive])).external_decision.iloc[0]
        == "conflicting_validated_evidence"
    )


def test_duplicate_records_and_invented_identity_rejected():
    with pytest.raises(ValueError, match="Unique"):
        assess_evidence(pd.DataFrame([observation(), observation()]))
    with pytest.raises(ValueError, match="outside"):
        assess_evidence(pd.DataFrame([observation(gene_id="ENSG99999999999")]))


def test_same_reagent_is_not_independent_even_if_claimed_in_a_record():
    row = observation(antibody_id="HPA003322")
    assert reconcile_hpa(hpa(), pd.DataFrame([row])).external_decision.iloc[0] == "unresolved"


def test_supplier_parser_keeps_descriptive_and_crossreactive_staining_unadjudicated():
    path = Path(__file__).resolve().parents[1] / "scripts/collect_external_ihc.py"
    spec = importlib.util.spec_from_file_location("collect_external_ihc", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    html = """<h1>Example</h1><h3>Staining Pattern in Normal Tissues</h3>
      <table><tr><td>Organ</td><td>Lung</td><td>Negative.</td></tr>
      <tr><td>Heart</td><td>Weak staining (cross-reactive).</td></tr></table>
      <a href="https://www.proteinatlas.org/ENSG00000164458">Target</a>
      <h3>Staining Pattern in Relevant Tumor Types</h3>
      <a href="https://www.proteinatlas.org/ENSG00000141510">Related product</a>"""
    record, rows = module.parse_panel(html, "https://example.org")
    assert record["source_gene_id"] == "ENSG00000164458"
    assert [r["result"] for r in rows] == ["not_detected", "reported_cross_reactivity"]
    assert not any(r["exact_cell_coverage"] for r in rows)

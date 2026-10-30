"""Clinical causality, shared peptide scope and bounded RNA credit."""

import json
from copy import deepcopy

import pandas as pd
import pytest

from oncoref import reference_data, tcr_trials
from oncoref.cta_selection_policy import assess_contexts


def test_reviewed_products_and_shared_peptide_not_family_transfer():
    p = tcr_trials.products().set_index("product_id")
    assert p.loc["lete_cel_c259", "rna_credit"]
    assert not p.loc["magea3_a3a", "rna_credit"]
    assert not p.loc["magea3_a2_neuro", "rna_credit"]
    assert not p.loc["tbi1301", "rna_credit"]
    targets = p.loc["lete_cel_c259", "targets"]
    assert {t["symbol"] for t in targets} == {"CTAG1A", "CTAG1B", "CTAG2"}
    assert {t["peptide"] for t in targets} == {"SLLMWITQC"}
    assert all(t["sequence_matches"] for t in targets)
    for ids in [["ENSG00000126890"], ["ENSG00000184033", "ENSG00000268651"]]:
        result = tcr_trials.protein_credit(ids)
        assert result["rna_credit"]
        assert result["cohort_ids"] == ("lete_mrcls_20",)
    assert not tcr_trials.protein_credit(["ENSG00000126890", "unknown"])["rna_credit"]
    assert not tcr_trials.protein_credit(["ENSG00000126890"], taxon_id=10090)["rna_credit"]
    assert not tcr_trials.protein_credit(["ENSG00000221867"])["rna_credit"]


def test_registration_unknown_or_missing_denominator_cannot_supply_credit():
    d = tcr_trials.evidence()
    c = next(c for c in d["cohorts"] if c["product_id"] == "ima203")
    c.update(review_status="pending", grade45_specificity="unassessed", patients_treated=None)
    assert not tcr_trials.products(record=d).set_index("product_id").loc["ima203", "rna_credit"]
    c["grade45_specificity"] = "none_identified"
    with pytest.raises(ValueError, match="denominator"):
        tcr_trials.validate(d)


def test_later_specificity_signal_blocks_product_even_with_favorable_cohort():
    d = tcr_trials.evidence()
    e = deepcopy(next(e for e in d["events"] if e["product_id"] == "ima203"))
    e.update(event_id="later_event", cohort_id=None, specificity_attribution="suspected", grade=5)
    d["events"].append(e)
    assert not tcr_trials.products(record=d).set_index("product_id").loc["ima203", "rna_credit"]
    e["grade"] = None
    assert not tcr_trials.products(record=d).set_index("product_id").loc["ima203", "rna_credit"]
    e["specificity_attribution"] = "unassessed"
    assert not tcr_trials.products(record=d).set_index("product_id").loc["ima203", "rna_credit"]


def test_mutations_do_not_leak_and_events_preserve_uncertainty():
    d = tcr_trials.evidence()
    d["products"].clear()
    assert tcr_trials.evidence()["products"]
    events = tcr_trials.events().set_index("event_id")
    assert pd.isna(events.loc["lete_gbs_cases", "grade"])
    assert events.loc["lete_late_cardiac_death", "treatment_attribution"] == "potentially_related"
    assert events.loc["lete_late_cardiac_death", "specificity_attribution"] == "not_established"
    assert events.loc["magea3_titin_deaths", "recognized_protein"] == "TTN"


def test_source_pins_match_uniform_registry():
    for source in tcr_trials.evidence()["sources"].values():
        spec = reference_data.REFERENCE_SOURCES[source["source_id"]]
        assert spec["pins"][source["version"]]["sha256"] == source["sha256"]
        assert spec["urls"][source["version"]] == source["artifact_url"]


def context(grade, tier="vital", value=6, credit=True):
    return {
        "source": "hpa_rna_consensus",
        "taxonomy_id": 9606,
        "stage": "adult",
        "tissue": "heart muscle",
        "cell_type": "bulk",
        "category": tier,
        "mapping_status": "mapped",
        "negative_ihc_grade": grade,
        "value": value,
        "unit": "nTPM",
        "clinical_rna_credit": credit,
        "clinical_product_ids": "lete_cel_c259" if credit else "",
    }


@pytest.mark.parametrize("tier,base", [("vital", 4), ("protected", 8), ("manageable", 16)])
def test_single_step_cap_and_original_ihc_preserved(tier, base):
    grades = ["Missing", "Approved", "Supported", "Enhanced"]
    result = assess_contexts(pd.DataFrame([context(g, tier, value=base * 2) for g in grades]))
    assert result.baseline_ceiling.tolist() == [base, base * 2, base * 4, base * 8]
    assert result.candidate_ceiling.tolist() == [base * 2, base * 4, base * 8, base * 8]
    assert result.negative_ihc_grade.tolist() == grades
    assert result.clinical_credit_applied.tolist() == [True, True, True, False]
    assert not result.candidate_exclusion.any()  # inclusive ceiling
    uncredited = assess_contexts(pd.DataFrame([context("Missing", tier, base * 2, False)]))
    assert uncredited.candidate_exclusion.iloc[0]
    assert uncredited.unsafe_evidence.iloc[0] == (tier != "manageable")


def test_credit_cannot_fill_missing_rna_or_unresolved_anatomy():
    row = context("Missing", value=float("nan"))
    result = assess_contexts(pd.DataFrame([row]))
    assert pd.isna(result.value.iloc[0])
    row.update(mapping_status="unresolved")
    assert assess_contexts(pd.DataFrame([row])).candidate_exclusion.iloc[0]
    row = context("Missing")
    row.pop("clinical_product_ids")
    with pytest.raises(ValueError, match="provenance"):
        assess_contexts(pd.DataFrame([row]))
    row = context("Missing", credit="False")
    with pytest.raises(ValueError, match="boolean"):
        assess_contexts(pd.DataFrame([row]))


def test_protein_annotation_requires_all_member_loci():
    facts = pd.DataFrame(
        {
            "member_gene_ids": [
                json.dumps(["ENSG00000184033", "ENSG00000268651"]),
                json.dumps(["ENSG00000184033", "unknown"]),
            ]
        }
    )
    annotated = tcr_trials.annotate_proteins(facts)
    assert annotated.clinical_rna_credit.tolist() == [True, False]
    assert "clinical_rna_credit" not in facts

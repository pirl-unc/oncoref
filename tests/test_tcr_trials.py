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
        "tissue": {"vital": "heart muscle", "protected": "colon", "manageable": "thyroid gland"}[
            tier
        ],
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


def test_potency_proven_gold_standard_calibration_and_cli(capsys):
    from oncoref.cli import main as cli_main

    gs_all = tcr_trials.gold_standard(potency_proven_only=False)
    gs_potent = tcr_trials.gold_standard(potency_proven_only=True)
    assert len(gs_all) == 16
    assert len(gs_potent) == 14
    assert {"magea1_ima202", "ct83_kklc1_phase1"} == (
        set(gs_all["benchmark_id"]) - set(gs_potent["benchmark_id"])
    )

    ev = tcr_trials.evaluate_against_gold_standard().set_index("benchmark_id")
    assert ev.loc["magea1_ima202", "benchmark_alignment"] == "excluded_potency_unproven"
    assert ev.loc["ct83_kklc1_phase1", "benchmark_alignment"] == "excluded_potency_unproven"
    for bid in (
        "ny_eso_1_1g4_c259",
        "ctag2_lete_cel_c259",
        "magea4_afami_cel",
        "prame_ima203",
        "magea10_adp_a2m10",
    ):
        assert ev.loc[bid, "oncoref_category"] == "CTA"
        assert ev.loc[bid, "benchmark_alignment"] == "true_safe_cta"

    for bid in (
        "ceacam5_nci_tcr",
        "mlana_dmf5_tcr",
        "pmel_gp100_tcr_imm",
        "tyr_til1383i_tcr",
        "magea3_a2_neuro",
        "magea3_a1_titin",
    ):
        assert bool(ev.loc[bid, "effective_oncoref_unsafe"]) is True
        assert ev.loc[bid, "benchmark_alignment"] == "true_unsafe_toxic"

    for bid in ("magea3_dp4_6f9", "afp_adp_a2afp", "wt1_ttcr_c4"):
        assert ev.loc[bid, "oncoref_category"] == "Unsafe"
        assert ev.loc[bid, "benchmark_alignment"] == "conservative_unsafe_epitope_tolerated"

    summary = tcr_trials.gold_standard_summary()
    assert summary["potency_proven_entries"] == 14
    assert summary["severe_specificity_toxic_entries"] == 6
    assert summary["severe_specificity_toxic_primary_unsafe"] == 6
    assert summary["severe_specificity_toxic_flagged_unsafe"] == 6
    assert summary["false_safe_misses"] == 0
    assert summary["primary_target_sensitivity_to_toxicity"] == 1.0
    assert summary["sensitivity_to_clinical_toxicity"] == 1.0
    assert summary["potency_proven_safe_cta_entries"] == 5
    assert summary["potency_proven_safe_cta_retained"] == 5
    assert summary["recall_on_safe_ctas"] == 1.0
    assert set(summary["conservative_non_cta_targets"]) == {"AFP", "MAGEA3", "WT1"}

    # Guard against string-bool coercion or contradicting stratum/toxicity flags
    bad_bool = gs_all.copy()
    bad_bool["potency_proven"] = bad_bool["potency_proven"].astype(str)
    with pytest.raises(ValueError, match="boolean"):
        tcr_trials.validate_gold_standard(bad_bool)

    bad_stratum = gs_all.copy()
    bad_stratum.loc[0, "severe_specificity_toxicity"] = True
    with pytest.raises(ValueError, match="inconsistent"):
        tcr_trials.validate_gold_standard(bad_stratum)

    assert cli_main(["tcr-trials", "--table", "benchmark_summary"]) == 0
    cli_payload = json.loads(capsys.readouterr().out)
    assert cli_payload["false_safe_misses"] == 0
    assert cli_payload["sensitivity_to_clinical_toxicity"] == 1.0

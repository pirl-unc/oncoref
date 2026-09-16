"""Scientific curation regressions: gaps, distinct arms, and audit completeness."""

import csv
import json
from pathlib import Path
from statistics import median

import pandas as pd
import pytest
from scripts import audit_tmb_ici_apd1 as audit
from scripts.recompute_neuroblastoma_tmb import summarize
from scripts.recompute_normal_like_breast_tmb import summarize as summarize_normal_like

from oncoref import apd1, ici, tmb
from oncoref.load_dataset import get_data


@pytest.mark.parametrize("code", ["UCEC_POLE", "UCEC_CNL", "LUAD_STK11", "MBL", "DIPG"])
def test_rejected_response_anchor_blocks_numeric_fallback(code):
    assert ici.cancer_ici_response(code) is None
    assert apd1.cancer_apd1_response(code) is None
    assert code not in apd1.cancer_apd1_response(include_inherited=True)
    record = apd1.cancer_apd1_response_record(code)
    assert record["inheritance_kind"] == "direct_missing"
    assert record["missing_reason"]
    assert record["apd1_orr_pct"] is None
    assert ici.pooled_ici_response(code, include_alternates=True)["pooled_pct"] is None


@pytest.mark.parametrize("code", ["STAD_MSI", "RCC", "PITNET"])
def test_remaining_tmb_gap_stays_missing_even_with_numeric_ancestor(code):
    assert tmb.cancer_tmb(code) is None
    source = tmb.resolve_tmb_source(code)
    assert source["inheritance_kind"] == "direct_missing"
    assert source["source_review_status"] == "audited_gap"


@pytest.mark.parametrize(
    "code,value,statistic,estimate_type,status",
    [
        ("HL", 7.66, "unspecified", "reported_summary", "source_checked"),
        ("UCEC_CNL", 2.9, "unspecified", "reported_summary", "source_checked"),
        ("UCEC_CNH", 2.3, "unspecified", "reported_summary", "source_checked"),
        ("CTCL", 3.5, "median", "subtype_proxy", "source_checked_proxy"),
        ("LUAD_EGFR", 3.8, "median", "broader_cohort_proxy", "source_checked_proxy"),
        ("BRCA_Normal", 1.28, "median", "sample_recomputed_median", "source_checked"),
        ("HCL", 0.2, "approximate", "approximate_capture_normalized", "approximation_reviewed"),
    ],
)
def test_best_effort_estimates_preserve_statistic_and_population_limits(
    code, value, statistic, estimate_type, status
):
    assert tmb.cancer_tmb(code) == value
    row = tmb.cancer_tmb_record(code)
    assert row["inheritance_kind"] == "direct"
    assert row["tmb_statistic"] == statistic
    assert row["estimate_type"] == estimate_type
    assert row["source_review_status"] == status
    assert row["missing_reason"] is None
    assert row["source_locator"] and row["source_scope"] and row["tmb_assay"]
    if statistic in {"unspecified", "approximate"}:
        assert row["median_tmb_mut_mb"] is None
        assert row["mean_tmb_mut_mb"] is None
        assert row["estimate_tmb_mut_mb"] == value


def test_normal_like_calculation_and_hcl_approximation_remain_reproducible():
    audit_dir = Path(__file__).resolve().parents[1] / "docs/audits"
    with (audit_dir / "tmb-normal-like-sample-rates.csv").open(newline="") as handle:
        samples = list(csv.DictReader(handle))
    result = summarize_normal_like(samples)
    assert (
        result == json.loads((audit_dir / "tmb-normal-like-recomputed.json").read_text())["results"]
    )
    assert result["n"] == 36
    assert {r["subtype"] for r in samples} == {"BRCA_Normal"}
    row = tmb.cancer_tmb_record("BRCA_Normal")
    assert row["median_tmb_mut_mb"] == round(result["median"], 2) == 1.28
    assert row["mean_tmb_mut_mb"] == round(result["mean"], 2) == 1.93
    assert (
        tmb.tmb_evidence_fields("BRCA_Normal", 1.93, statistic="mean")["estimate_type"]
        == "sample_recomputed_mean"
    )
    assert row["confidence"] == "low"
    with pytest.raises(ValueError, match="Duplicate"):
        summarize_normal_like([*samples, samples[0]])
    with (audit_dir / "tmb-hcl-source-counts.csv").open(newline="") as handle:
        counts = list(csv.DictReader(handle))
    hcl = json.loads((audit_dir / "tmb-hcl-approximation.json").read_text())
    assert len(counts) == 5
    for sensitivity in hcl["sensitivity"]:
        count = median(int(r[sensitivity["count_definition"]]) for r in counts)
        assert count == sensitivity["median_count"]
        assert [count / d for d in hcl["nominal_capture_mb"]] == sensitivity[
            "median_count_per_nominal_mb"
        ]
    assert hcl["selected_rounded_estimate_mut_mb"] == tmb.cancer_tmb("HCL") == 0.2
    assert tmb.cancer_tmb_record("HCL")["confidence"] == "low"


def test_proxy_median_and_same_cohort_mean_stay_distinct():
    row = tmb.cancer_tmb_record("LUAD_EGFR")
    assert row["tmb_mut_mb"] == row["median_tmb_mut_mb"] == 3.8
    assert row["mean_tmb_mut_mb"] == 5.6
    assert row["n_samples"] == 383
    assert row["source_scope"] == "egfr_mutant_lung_proxy_for_luad"
    row = tmb.cancer_tmb_record("CTCL")
    assert row["n_samples"] == 67
    assert row["source_scope"] == "mycosis_fungoides_proxy_for_ctcl"


def test_typed_estimate_columns_are_consistent_and_inventory_flags_survive():
    raw = get_data("cancer-tmb")
    selected = raw[raw["estimate_tmb_mut_mb"].notna()]
    assert selected["estimate_statistic"].isin({"unspecified", "approximate"}).all()
    assert selected[["median_tmb_mut_mb", "mean_tmb_mut_mb"]].isna().all().all()
    assert raw.loc[raw["estimate_tmb_mut_mb"].isna(), "estimate_statistic"].isna().all()
    rows = {r["cancer_code"]: r for r in audit.audit() if r["dataset"] == "cancer-tmb"}
    assert "summary_statistic_unspecified" in rows["HL"]["findings"]
    assert rows["HL"]["metric"] == "TMB_UNSPECIFIED"
    assert "population_proxy" in rows["CTCL"]["findings"]
    assert "approximate_tmb_estimate" in rows["HCL"]["findings"]
    assert rows["HCL"]["metric"] == "TMB_APPROXIMATE"


def test_tmb_review_covers_exactly_the_curated_rows_and_gates_median_provenance():
    raw = get_data("cancer-tmb")
    review = get_data("cancer-tmb-source-audit")
    assert review["cancer_code"].is_unique
    assert set(review["cancer_code"]) == set(raw["cancer_code"])
    frame = tmb.cancer_tmb_df()
    checked = frame[frame["source_review_status"] == "source_checked"]
    assert checked["source_locator"].notna().all()
    assert checked["tmb_assay"].notna().all()
    assert set(checked["estimate_type"]) == {
        "published_median",
        "published_mean",
        "sample_recomputed_median",
        "reported_summary",
    }
    unpublished = frame[frame["source_review_status"] != "source_checked"]
    assert "published_median" not in set(unpublished["estimate_type"])
    assert "published_mean" not in set(unpublished["estimate_type"])
    assert frame.set_index("cancer_code").loc["FL", "median_tmb_mut_mb"] == 5.05
    assert frame.set_index("cancer_code").loc["BCC", "median_tmb_mut_mb"] == 47.3


def test_recovered_source_medians_preserve_genomic_and_response_distinction():
    rows = tmb.cancer_tmb_df().set_index("cancer_code")
    assert rows.loc["UVM", "median_tmb_mut_mb"] == 0.34
    assert rows.loc["UVM", "pmid_doi"] == "PMID:31930041"
    assert rows.loc["UVM", "n_samples"] == 80
    pole = tmb.resolve_tmb_source("UCEC_POLE")
    assert pole["median_tmb_mut_mb"] == 150.8
    assert pole["n_samples"] == 61
    assert pole["estimate_type"] == "published_median"
    assert pole["source_scope"] == "pathogenic_pole_including_multiple_classifiers"
    assert ici.cancer_ici_response("UCEC_POLE") is None
    assert apd1.cancer_apd1_response("UCEC_POLE") is None


def test_recovered_uveal_anchor_is_a_prospective_pd1_cohort():
    assert ici.cancer_ici_response("UVM") == 11.7
    assert apd1.cancer_apd1_response("UVM") == 11.7
    pool = ici.pooled_ici_response("UVM")
    assert pool["selected_regimen"] == "PD-1"
    assert pool["n_studies"] == 1
    row = ici.cancer_ici_response_estimates_df().set_index("estimate_id").loc["ICI-e23fe5ef46-01"]
    assert row["ref"] == "PMID:31175402"
    assert row["responders"] == 2
    assert row["metric_n"] == 17
    assert "first-line" in row["setting"]


def test_stk11_subgroup_denominators_and_regimens_are_not_whole_cohorts():
    rows = ici.cancer_ici_response_estimates_df().set_index("estimate_id")
    su2c = rows.loc["ICI-2e06c5f7d3-01"]
    assert su2c["metric_n"] == 54
    assert su2c["responders"] == 4
    assert su2c["source_n"] == 174
    assert su2c["regimen"] == "PD-1+/-CTLA-4"
    cm057 = rows.loc["ICI-1c401c554d-01"]
    assert cm057["metric_n"] == 6
    assert cm057["responders"] == 0
    assert cm057["source_n"] == 24
    assert pd.isna(cm057["ci_high"])
    assert cm057["ci_basis"] == "not_reported"
    assert cm057["trial_nct"] == "NCT01673867"


def test_withdrawn_claim_ledgers_cover_all_original_withdrawals():
    audit_dir = Path(__file__).resolve().parents[1] / "docs/audits"
    with (audit_dir / "tmb-source-recovery.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == len({r["cancer_code"] for r in rows}) == 20
    assert all(r["original_source_assessment"] and r["decision_notes"] for r in rows)
    with (audit_dir / "ici-withdrawn-anchor-verification.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {r["cancer_code"] for r in rows} == {
        "COAD",
        "READ",
        "UCEC",
        "UCEC_CNL",
        "UCEC_CNH",
        "UCEC_POLE",
        "LUAD_STK11",
        "DIPG",
        "MBL",
        "UVM",
    }
    assert all(r["claimed_source_result"] and r["source_locator"] for r in rows)


def test_neuroblastoma_reanalysis_reproduces_published_overall_and_separates_subgroups():
    audit_dir = Path(__file__).resolve().parents[1] / "docs/audits"
    with (audit_dir / "tmb-pugh2013-sample-rates.csv").open(newline="") as handle:
        samples = list(csv.DictReader(handle))
    assert len(samples) == len({r["case_id"] for r in samples}) == 240
    assert sum(r["mycn_amp_status"] == "9" for r in samples) == 5
    results = summarize(samples)
    provenance = json.loads((audit_dir / "tmb-pugh2013-recomputed.json").read_text())
    assert results == provenance["results"]
    by_group = {(r["group"], r["measure"]): r for r in results}
    assert round(by_group["all_high_risk", "total_exonic_mut_mb"]["median"], 2) == 0.60
    assert round(by_group["all_high_risk", "nonsilent_mut_mb"]["median"], 2) == 0.48
    amplified = by_group["MYCNamp", "nonsilent_mut_mb"]
    assert amplified["n"] == 77
    row = tmb.cancer_tmb_record("NBL_MYCNamp")
    assert row["tmb_mut_mb"] == round(amplified["median"], 2) == 0.43
    assert row["estimate_type"] == "sample_recomputed_median"
    assert row["source_scope"] == "high_risk_mycn_amplified_cohort"
    nonamp = tmb.cancer_tmb_record("NBL_MYCNnonamp")
    assert nonamp["tmb_mut_mb"] == 0.66
    assert nonamp["n_samples"] == 58
    assert nonamp["pmid_doi"] == "PMID:33172452"
    assert nonamp["estimate_type"] == "published_median"


def test_chalmers_claims_match_extracted_published_rows_and_specimen_counts():
    path = Path(__file__).resolve().parents[1] / "docs/audits/tmb-chalmers-source-rows.csv"
    rows = tmb.cancer_tmb_df().set_index("cancer_code")
    with path.open(newline="") as handle:
        for source in csv.DictReader(handle):
            row = rows.loc[source["cancer_code"]]
            assert row["median_tmb_mut_mb"] == float(source["median_tmb_mut_mb"])
            assert row["n_samples"] == int(source["n_specimens"])
            assert row["source_review_status"] == "source_checked"
            assert source["source_locator"] in row["source_locator"]
            assert row["tmb_assay"] == "targeted_panel_coding_including_synonymous"


def test_pool_defaults_to_one_primary_regimen_and_blocks_overlapping_trial_updates():
    pool = ici.pooled_ici_response("CRC_MSI")
    assert pool["selected_regimen"] == "PD-1"
    assert {s["role"] for s in pool["sources"]} == {"primary"}
    assert {s["regimen"] for s in pool["sources"]} == {"PD-1"}
    assert pool["n_studies"] == 1
    overlap = ici.pooled_ici_response("BLCA", regimen="PD-1", include_alternates=True)
    assert not overlap["poolable"]
    assert overlap["pooled_pct"] is None
    assert overlap["pooling_block_reason"] == "potentially_overlapping_cohorts"
    assert overlap["n_studies"] > 1


def test_controls_and_mixed_agent_cohorts_never_enter_ici_pool():
    for code, regimen in [("CRC_MSI", "non-ICI"), ("UVM", "PD-(L)1")]:
        pool = ici.pooled_ici_response(
            code, regimen=regimen, include_alternates=True, verified_only=False
        )
        assert pool["n_studies"] == 0
        assert pool["pooled_pct"] is None


def test_response_denominators_are_evaluable_populations_and_case_ci_is_absent():
    frame = ici.cancer_ici_response_estimates_df().set_index("estimate_id")
    mds = frame.loc["ICI-a5110c10af-01"]
    assert float(mds["source_n"]) == 28
    assert float(mds["metric_n"]) == 27
    assert float(mds["ci_high"]) == 12.5  # 95% Wilson for 0/27, not source 90% CI.
    assert float(frame.loc["ICI-e680dc5342-01", "metric_n"]) == 153
    assert float(frame.loc["ICI-ba043c8be9-01", "responders"]) == 90
    case = frame.loc["ICI-64a41d348a-01"]
    assert case["metric"] == "PR_COUNT"
    assert case["value_basis"] == "reported_context"
    assert pd.isna(case["ci_low"]) and pd.isna(case["ci_high"])


def test_full_inventory_has_no_unresolved_arithmetic_or_anchor_mismatch():
    rows = audit.audit()
    assert len(rows) == sum(
        len(audit.read_table(audit.ROOT / "oncoref/data", name)) for name in audit.TABLES
    )
    errors = {
        "invalid_numeric_value",
        "negative_value",
        "percentage_out_of_bounds",
        "invalid_denominator",
        "invalid_numerator",
        "reversed_ci",
        "value_outside_ci",
        "count_rate_mismatch",
        "anchor_reference_mismatch",
        "anchor_value_differs_from_endpoint",
        "missing_primary_endpoint",
        "drug_regimen_mismatch",
        "nonreported_population_anchor",
    }
    for row in rows:
        assert not errors.intersection(row["findings"].split(";")), row
    # A clean arithmetic result must not certify an unreviewed TMB source.
    legacy = next(r for r in rows if r["dataset"] == "cancer-tmb" and r["cancer_code"] == "CRC_MSI")
    assert "tmb_source_not_revalidated" in legacy["findings"]


def test_inventory_catches_wrong_denominator_and_drug_identity(tmp_path):
    for name in (*audit.TABLES, "cancer-tmb-source-audit"):
        rows = audit.read_table(audit.ROOT / "oncoref/data", name)
        if name == "cancer-ici-response-estimates":
            target = next(r for r in rows if r["estimate_id"] == "ICI-e045c0d841-01")
            target.update(metric_n="22", drug="nivolumab+ipilimumab", regimen="PD-1")
        with (tmp_path / f"{name}.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    result = next(r for r in audit.audit(tmp_path) if r["row_key"] == "ICI-e045c0d841-01")
    assert {"drug_regimen_mismatch", "count_rate_mismatch"} <= set(result["findings"].split(";"))


def test_checked_in_inventory_is_reproducible():
    path = Path(__file__).resolve().parents[1] / "docs/audits/tmb-ici-apd1.csv"
    with path.open(newline="") as handle:
        assert list(csv.DictReader(handle)) == audit.audit()

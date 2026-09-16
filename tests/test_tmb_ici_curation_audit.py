"""Scientific curation regressions: gaps, distinct arms, and audit completeness."""

import csv
from pathlib import Path

import pandas as pd
import pytest
from scripts import audit_tmb_ici_apd1 as audit

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


@pytest.mark.parametrize("code", ["BRCA_Normal", "UCEC_CNL", "RB", "HL", "CTCL", "LUAD_EGFR"])
def test_rejected_tmb_median_stays_missing_even_with_numeric_ancestor(code):
    assert tmb.cancer_tmb(code) is None
    source = tmb.resolve_tmb_source(code)
    assert source["inheritance_kind"] == "direct_missing"
    assert source["source_review_status"] == "rejected_population_median"


def test_tmb_review_covers_exactly_the_curated_rows_and_gates_median_provenance():
    raw = get_data("cancer-tmb")
    review = get_data("cancer-tmb-source-audit")
    assert review["cancer_code"].is_unique
    assert set(review["cancer_code"]) == set(raw["cancer_code"])
    frame = tmb.cancer_tmb_df()
    checked = frame[frame["source_review_status"] == "source_checked"]
    assert checked["source_locator"].notna().all()
    assert checked["tmb_assay"].notna().all()
    assert set(checked["estimate_type"]) == {"published_median"}
    unpublished = frame[frame["source_review_status"] != "source_checked"]
    assert "published_median" not in set(unpublished["estimate_type"])
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

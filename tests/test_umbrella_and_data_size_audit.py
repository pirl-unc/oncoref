# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Tests for umbrella expression coverage, proxy cohort audit, and wheel size audit (#737, #738)."""

from __future__ import annotations

from pathlib import Path

import pytest

import oncoref
from oncoref import data_manifest


def test_cancer_reference_expression_availability_include_aggregates():
    default_avail = oncoref.cancer_reference_expression_availability(
        reference_source="summary_rows", sample_qc="all"
    )
    assert "NSCLC" not in set(default_avail["requested_code"])
    assert "THYM_EPITHELIAL" not in set(default_avail["requested_code"])

    with_agg = oncoref.cancer_reference_expression_availability(
        reference_source="summary_rows", sample_qc="all", include_aggregates=True
    )
    agg_rows = with_agg[with_agg["request_kind"].eq("aggregate_member")]
    requested_umbrellas = set(agg_rows["requested_code"])
    expected_targets = {
        "THYM_EPITHELIAL",
        "NSCLC",
        "CRC",
        "CRC_MSI",
        "BTC",
        "RCC",
        "SARC",
        "SARC_ESS",
        "SARC_LPS",
        "SARC_RMS",
        "NEN",
        "NET",
        "NEC",
        "NEC_LUNG",
    }
    assert expected_targets <= requested_umbrellas

    with pytest.raises(ValueError, match="include_aggregates=True"):
        oncoref.cancer_reference_expression_availability(
            reference_source="summary_rows_all",
            sample_qc="all",
            all_sources=True,
            include_aggregates=True,
        )


def test_umbrella_expression_coverage_surfaces_thym_complete_and_proxy_members():
    cov = oncoref.umbrella_expression_coverage().set_index("cancer_code")
    assert "THYM_EPITHELIAL" in cov.index
    assert "THYMCA" not in cov.index
    assert "NSCLC" in cov.index
    assert "SARC_ESS" in cov.index
    assert "SARC_LPS" in cov.index

    thym_ep = cov.loc["THYM_EPITHELIAL"]
    assert thym_ep["available_member_codes"] == "THYM;THYMCA"
    assert thym_ep["unbacked_child_subtypes"] == ""
    assert int(thym_ep["total_reference_samples"]) == 119
    assert thym_ep["coverage_status"] == "complete_linear_tpm"

    assert cov.loc["NSCLC", "coverage_status"] == "complete_linear_tpm"
    assert cov.loc["SARC_ESS", "coverage_status"] == "includes_proxy_scale_members"
    assert cov.loc["SARC_LPS", "coverage_status"] == "partial_member_gap"
    assert cov.loc["SARC_LPS", "unbacked_child_subtypes"] == "SARC_MPLPS"
    assert cov.loc["SARC_LPS", "non_comparable_member_codes"] == "SARC_MYXLPS"


def test_non_comparable_expression_cohorts_inventories_all_9_non_linear_cohorts():
    all_nc = oncoref.non_comparable_expression_cohorts()
    assert set(all_nc["cancer_code"]) == {
        "CMN",
        "CTCL",
        "EPN",
        "HCL",
        "MTC",
        "NUTM",
        "SARC_ESS_HG",
        "SARC_ESS_LG",
        "SARC_MYXLPS",
    }
    assert "public_rnaseq_audit_status" in all_nc.columns
    assert all_nc["public_rnaseq_audit_status"].astype(str).str.strip().ne("").all()
    nutm = all_nc.set_index("cancer_code").loc["NUTM"]
    assert bool(nutm["is_invalid_withheld"]) is True
    assert nutm["allowed_uses"] == "none_withheld"

    proxy_only = oncoref.non_comparable_expression_cohorts(include_invalid=False)
    assert len(proxy_only) == 8
    assert "NUTM" not in set(proxy_only["cancer_code"])
    assert "bulk_rnaseq_cpm_proxy" not in set(proxy_only["scale_family"])
    assert set(proxy_only["allowed_uses"]) == {"within_sample_rank_and_pathway_only"}


def test_package_data_size_audit_accounts_for_every_shipped_file():
    data_root = Path(oncoref.__file__).resolve().parent / "data"
    disk_files = [p for p in data_root.rglob("*") if p.is_file()]
    expected_bytes = sum(p.stat().st_size for p in disk_files)

    audit = oncoref.package_data_size_audit()
    assert len(audit) == len(disk_files)
    assert int(audit["size_bytes"].sum()) == expected_bytes
    assert set(audit["retention_tier"]) == set(data_manifest.RETENTION_TIERS)

    by_path = audit.set_index("relative_path")
    assert by_path.loc["target-selection.csv.gz", "retention_tier"] == "core_runtime"
    assert bool(by_path.loc["target-selection.csv.gz", "shipped_in_wheel"]) is True
    assert (
        by_path.loc["translation-qualification.csv.gz", "retention_tier"]
        == "active_calibration_or_audit"
    )
    assert bool(by_path.loc["translation-qualification.csv.gz", "shipped_in_wheel"]) is False
    assert (
        by_path.loc["cta-selection-audit-human.json.gz", "retention_tier"]
        == "archival_compatibility_snapshot"
    )
    assert bool(by_path.loc["cta-selection-audit-human.json.gz", "shipped_in_wheel"]) is False

    wheel_only_bytes = int(audit.loc[audit["shipped_in_wheel"], "size_bytes"].sum())
    assert wheel_only_bytes < 33_500_000

    summary = oncoref.package_data_size_summary()
    assert int(summary["n_files"].sum()) == len(disk_files)
    assert int(summary["total_bytes"].sum()) == expected_bytes
    # Core runtime stays comfortably under 25 MB once archival & calibration tiers are separated.
    core_mb = float(summary.loc[summary["retention_tier"].eq("core_runtime"), "total_mb"].iloc[0])
    assert core_mb < 25.0

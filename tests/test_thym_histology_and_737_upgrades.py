"""Regression tests for TCGA-THYM WHO histology split (THYM/THYMCA) and #737 linear-TPM upgrades."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from scripts.build_tcga_thym_histology import (
    SOURCE_COHORT,
    histology_call,
    membership,
    verified_bytes,
)

import oncoref
from oncoref import cancer_types, samples, source_matrices


@pytest.mark.parametrize(
    ("pan_tt", "tcga_hd", "pan_icdo", "tcga_icdo", "expected_code", "expected_reason"),
    [
        (
            "Thymoma, Type C",
            "Thymoma; Type C",
            "8586/3",
            "8586/3",
            "THYMCA",
            "concordant_thymic_carcinoma",
        ),
        ("Thymoma, Type Ab", "Thymoma; Type AB", "8582/3", "8582/3", "THYM", "concordant_thymoma"),
        ("", "Thymoma; Type B2", "", "8584/3", "THYM", "concordant_thymoma"),
        (
            "Thymoma, Type C",
            "Thymoma; Type B3",
            "8586/3",
            "8585/3",
            "indeterminate",
            "conflicting_histology_fields",
        ),
        ("", "", "", "", "indeterminate", "missing_histology_metadata"),
        ("Neuroendocrine tumor", "", "", "", "indeterminate", "unrecognized_who_histology"),
        (
            "Thymoma, Type C",
            "Thymoma; Type C",
            "9999/3",
            "9999/3",
            "indeterminate",
            "unrecognized_icdo_morphology",
        ),
    ],
)
def test_histology_call_requires_concordance_and_rejects_conflicts(
    pan_tt, tcga_hd, pan_icdo, tcga_icdo, expected_code, expected_reason
):
    code, _subtype, reason = histology_call(pan_tt, tcga_hd, pan_icdo, tcga_icdo)
    assert code == expected_code
    assert reason == expected_reason


def test_membership_audits_primary_barcodes_and_rejects_duplicate_patients():
    pan_sample = pd.DataFrame(
        {
            "SAMPLE_ID": ["TCGA-AA-0001-01", "TCGA-AA-0002-01", "TCGA-AA-0003-06"],
            "TUMOR_TYPE": ["Thymoma, Type Ab", "Thymoma, Type C", "Thymoma, Type C"],
        }
    )
    pan_patient = pd.DataFrame(
        {
            "PATIENT_ID": ["TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003"],
            "ICD_O_3_HISTOLOGY": ["8582/3", "8586/3", "8586/3"],
        }
    )
    tcga_patient = pd.DataFrame(
        {
            "PATIENT_ID": ["TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003"],
            "HISTOLOGICAL_DIAGNOSIS": ["Thymoma; Type AB", "Thymoma; Type C", "Thymoma; Type C"],
            "ICD_O_3_HISTOLOGY": ["8582/3", "8586/3", "8586/3"],
        }
    )
    audit = membership(
        ["TCGA-AA-0001-01", "TCGA-AA-0002-01", "TCGA-AA-0003-06"],
        pan_sample,
        pan_patient,
        tcga_patient,
    ).set_index("sample_id")
    assert audit.loc["TCGA-AA-0001-01", "routed_cancer_code"] == "THYM"
    assert audit.loc["TCGA-AA-0002-01", "routed_cancer_code"] == "THYMCA"
    assert not bool(audit.loc["TCGA-AA-0003-06", "included"])
    assert audit.loc["TCGA-AA-0003-06", "selection_reason"] == "non_primary_sample"


def test_checked_in_thym_histology_audit_and_registry_reconcile():
    root = Path(__file__).resolve().parents[1] / "docs/audits"
    audit = pd.read_csv(root / "thym-histology-membership.csv", keep_default_na=False)
    build_json = json.loads((root / "thym-histology-build.json").read_text())

    assert len(audit) == 119
    assert audit["included"].astype(bool).all()
    assert audit["routed_cancer_code"].value_counts().to_dict() == {"THYM": 108, "THYMCA": 11}
    assert build_json["n_samples_by_code"] == {"THYM": 108, "THYMCA": 11}
    assert build_json["n_indeterminate"] == 0

    assert cancer_types.resolve_cancer_type("thymic carcinoma") == "THYMCA"
    assert cancer_types.resolve_cancer_type("thymoma") == "THYM"
    records = cancer_types.cancer_type_records(["THYM", "THYMCA", "THYM_EPITHELIAL"]).set_index(
        "code"
    )
    assert records.loc["THYM", "source_matrix_cohort"] == SOURCE_COHORT
    assert records.loc["THYMCA", "source_matrix_cohort"] == SOURCE_COHORT
    assert int(records.loc["THYM", "source_matrix_n_samples"]) == 108
    assert int(records.loc["THYMCA", "source_matrix_n_samples"]) == 11
    assert records.loc["THYMCA", "reference_source"] == "own_cohort"
    assert bool(records.loc["THYMCA", "is_classification_target"]) is True
    assert source_matrices.source_sample_namespace(SOURCE_COHORT) == "TREEHOUSE_POLYA_25_01"

    sm_thym = samples.samples_for_cancer_code("THYM")
    sm_thymca = samples.samples_for_cancer_code("THYMCA")
    assert len(sm_thym) == 108
    assert len(sm_thymca) == 11
    assert set(sm_thym["sample_id"]).isdisjoint(set(sm_thymca["sample_id"]))


def test_gse125285_bcc_and_cscc_upgraded_to_linear_rnaseq_tpm():
    for code, expected_n in (("BCC", 25), ("cSCC", 10)):
        meta = oncoref.cancer_reference_expression_source_metadata(code)
        assert meta["source_cohort"] == "GSE125285_BCC_CSCC"
        assert meta["source_scale_class"] == "linear_rnaseq_tpm"
        assert meta["linear_tpm_comparable"] is True
        assert meta["tpm_proxy"] is False
        assert int(meta["n_reference_samples"]) == expected_n
        manifest_rows = samples.samples_for_cancer_code(code)
        assert len(manifest_rows) == expected_n
        assert set(manifest_rows["raw_unit"]) == {"TPM"}


def test_verified_bytes_rejects_tampered_file(tmp_path):
    p = tmp_path / "tampered.tsv"
    p.write_text("bad")
    with pytest.raises(ValueError, match="SHA-256"):
        verified_bytes(p, "0" * 64)

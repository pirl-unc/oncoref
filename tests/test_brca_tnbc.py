"""Receptor-defined TNBC must never borrow the PAM50 basal-like identity."""

from pathlib import Path

import pandas as pd
import pytest
from scripts.build_brca_tnbc import (
    SOURCE_COHORT,
    membership,
    receptor_call,
    tmb_rows,
    verified_bytes,
)

from oncoref import cancer_types, ici, source_matrices


@pytest.mark.parametrize(
    "er,pr,ihc,fish,state,her2",
    [
        ("Negative", "Negative", "Negative", "[Not Available]", "TNBC", "Negative"),
        ("Negative", "Negative", "Equivocal", "Negative", "TNBC", "Negative"),
        ("Negative", "Negative", "[Not Available]", "Negative", "TNBC", "Negative"),
        ("Negative", "Negative", "Indeterminate", "Negative", "TNBC", "Negative"),
        ("Negative", "Negative", "Equivocal", "[Not Available]", "indeterminate", "Indeterminate"),
        ("Negative", "Negative", "Negative", "Positive", "indeterminate", "Indeterminate"),
        ("Negative", "Negative", "Positive", "Negative", "indeterminate", "Indeterminate"),
        ("[Not Available]", "Negative", "Negative", "Negative", "indeterminate", "Negative"),
        ("Negative", "Indeterminate", "Negative", "Negative", "indeterminate", "Negative"),
        ("Negative", "Negative", "Equivocal", "Positive", "non_TNBC", "Positive"),
        (
            "Positive",
            "[Not Available]",
            "Equivocal",
            "[Not Available]",
            "non_TNBC",
            "Indeterminate",
        ),
        ("Negative", "Positive", "Negative", "Negative", "non_TNBC", "Negative"),
    ],
)
def test_receptor_calls_preserve_unknowns_and_conflicts(er, pr, ihc, fish, state, her2):
    assert receptor_call(er, pr, ihc, fish)[:2] == (state, her2)


def _clinical():
    receptors = pd.DataFrame(
        {
            "PATIENT_ID": ["TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003"],
            "ER_STATUS_BY_IHC": ["Negative", "Positive", "Negative"],
            "PR_STATUS_BY_IHC": ["Negative"] * 3,
            "IHC_HER2": ["Negative", "Negative", "Equivocal"],
            "HER2_FISH_STATUS": ["[Not Available]"] * 3,
        }
    )
    pam50 = pd.DataFrame(
        {
            "PATIENT_ID": receptors["PATIENT_ID"],
            "SUBTYPE": ["BRCA_LumA", "BRCA_Basal", "BRCA_Basal"],
        }
    )
    return receptors, pam50


def test_membership_uses_receptors_not_pam50_and_preserves_every_column():
    receptors, pam50 = _clinical()
    ids = [
        "TCGA-AA-0001-01",
        "TCGA-AA-0001-06",
        "TCGA-AA-0002-01",
        "TCGA-AA-0003-01",
        "TCGA-AA-0004-01",
    ]
    audit = membership(ids, receptors, pam50).set_index("sample_id")
    assert audit.index.tolist() == ids
    assert audit.index[audit["included"]].tolist() == ["TCGA-AA-0001-01"]
    assert audit.loc[ids[1], "selection_reason"] == "non_primary_sample"
    assert audit.loc[ids[2], "receptor_state"] == "non_TNBC"
    assert audit.loc[ids[3], "receptor_state"] == "indeterminate"
    assert audit.loc[ids[4], "receptor_reason"] == "missing_receptor_metadata"
    assert audit.loc[ids[4], "pam50"] == "unavailable"
    # Neither changing nor removing PAM50 calls can alter receptor selection.
    pam50["SUBTYPE"] = pd.NA
    assert membership(ids, receptors, pam50)["included"].tolist() == audit["included"].tolist()


def test_ambiguous_identifiers_fail_instead_of_duplicating_patients():
    receptors, pam50 = _clinical()
    with pytest.raises(ValueError, match="duplicate receptor"):
        membership(["TCGA-AA-0001-01"], pd.concat([receptors, receptors]), pam50)
    with pytest.raises(ValueError, match="Duplicate expression"):
        membership(["TCGA-AA-0001-01"] * 2, receptors, pam50)
    with pytest.raises(ValueError, match="barcode"):
        membership(["TCGA-AA-0001"], receptors, pam50)


def test_tmb_uses_exact_selected_primary_samples_and_retains_missingness():
    audit = pd.DataFrame(
        {
            "sample_id": [f"TCGA-AA-000{i}-01" for i in range(1, 6)],
            "patient_id": [f"TCGA-AA-000{i}" for i in range(1, 6)],
            "included": [True, True, True, True, False],
        }
    )
    samples = audit.iloc[[0, 1, 2, 4]].rename(
        columns={
            "sample_id": "SAMPLE_ID",
            "patient_id": "PATIENT_ID",
        }
    )
    samples["TMB_NONSYNONYMOUS"] = [2.5, 10, float("inf"), 99]
    result = tmb_rows(audit, samples, {"TCGA-AA-0001-01", "TCGA-AA-0003-01", "TCGA-AA-0005-01"})
    assert result["tmb_included"].tolist() == [True, False, False, False]
    assert result["nonsynonymous_mut_mb"].dropna().tolist() == [2.5]
    assert result["tmb_reason"].tolist() == [
        "included",
        "not_mutation_profile_eligible",
        "missing_or_invalid_TMB",
        "missing_sample_metadata",
    ]
    samples.loc[samples.index[0], "PATIENT_ID"] = "TCGA-ZZ-9999"
    with pytest.raises(ValueError, match="identity mismatch"):
        tmb_rows(audit, samples, set())


def test_unreviewed_input_is_rejected(tmp_path):
    path = tmp_path / "metadata.tsv"
    path.write_text("changed input")
    with pytest.raises(ValueError, match="SHA-256"):
        verified_bytes(path, "0" * 64)


def test_registry_and_response_populations_stay_distinct():
    assert cancer_types.resolve_cancer_type("TNBC") == "BRCA_TNBC"
    records = cancer_types.cancer_type_records(["BRCA_TNBC", "BRCA_Basal"]).set_index("code")
    assert records.loc["BRCA_TNBC", "parent_code"] == "BRCA"
    assert records.loc["BRCA_TNBC", "source_matrix_cohort"] == SOURCE_COHORT
    assert records.loc["BRCA_TNBC", "source_matrix_n_samples"] == 157
    assert records.loc["BRCA_Basal", "source_matrix_n_samples"] == 172
    assert source_matrices.source_sample_namespace(SOURCE_COHORT) == (
        source_matrices.source_sample_namespace(records.loc["BRCA_Basal", "source_matrix_cohort"])
    )
    assert ici.cancer_ici_response("BRCA_TNBC") == 5.3
    assert ici.cancer_ici_response("BRCA_Basal") is None
    assert ici.resolve_ici_response_source("BRCA_Basal")["missing_reason"] == (
        "tnbc_trials_do_not_estimate_pam50_basal_response"
    )
    data = Path(__file__).resolve().parents[1] / "oncoref/data"
    for name in (
        "cancer-ici-response.csv",
        "cancer-apd1-response.csv",
        "cancer-ici-response-estimates.csv",
    ):
        frame = pd.read_csv(data / name)
        trials = frame[frame["trial_name"].str.contains("KEYNOTE-086|PCD4989g", na=False)]
        assert not trials.empty
        assert set(trials["cancer_code"]) == {"BRCA_TNBC"}
    endpoints = pd.read_csv(data / "cancer-ici-response-estimates.csv")
    assert len(endpoints[endpoints["cancer_code"].eq("BRCA_TNBC")]) == 19


def test_checked_in_membership_audit_reconciles_with_registry_and_tmb():
    root = Path(__file__).resolve().parents[1] / "docs/audits"
    audit = pd.read_csv(root / "brca-receptor-membership.csv")
    assert len(audit) == 1099
    primary = audit[audit["sample_type"].eq("primary")]
    assert primary["receptor_state"].value_counts().to_dict() == {
        "non_TNBC": 858,
        "TNBC": 157,
        "indeterminate": 77,
    }
    selected = audit[audit["included"]]
    assert not selected["patient_id"].duplicated().any()
    assert selected["pam50"].eq("BRCA_Basal").sum() == 115
    rates = pd.read_csv(root / "brca-tnbc-tmb.csv")
    assert set(rates["sample_id"]) == set(selected["sample_id"])
    measured = rates.loc[rates["tmb_included"], "nonsynonymous_mut_mb"]
    assert len(measured) == 156
    assert measured.median() == pytest.approx(2.1166666665)


@pytest.mark.parametrize("strict_pd1", [True, False])
def test_response_tmb_plot_joins_the_tnbc_population(monkeypatch, strict_pd1):
    from oncoref import plots

    monkeypatch.setattr(plots, "_family_scatter", lambda rows, **kwargs: rows)
    points = {code: (tmb, orr) for code, tmb, orr in plots.apd1_vs_tmb(strict_pd1=strict_pd1)}
    assert "BRCA_Basal" not in points
    assert "BRCA" not in points
    assert points["BRCA_TNBC"] == (2.12, 5.3)

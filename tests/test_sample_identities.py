"""Patient evidence must survive matrix joins without inflating replication."""

import pandas as pd
import pytest

from oncoref import expression, samples


def resolve(source, ids):
    return samples.sample_identities(pd.DataFrame({"source_cohort": source, "sample_id": ids}))


def test_real_case_mappings_and_treehouse_repeat_libraries():
    beat = resolve("BEATAML_OHSU_2022", ["BA2000R", "BA2003R"])
    assert beat.participant_id.tolist() == ["BEATAML_OHSU_2022:2304", "BEATAML_OHSU_2022:2491"]
    tree = resolve("TREEHOUSE_POLYA_25_01", ["TH06_0613_S01", "TH06_0613_S02"])
    assert tree.participant_id.nunique() == 1
    assert tree.sample_key.nunique() == 2
    assert tree.participant_id.iloc[0] == "TREEHOUSE:TH06_0613"
    alternate = resolve("TREEHOUSE_RIBOD_25_01", ["TH06_0613_S03"])
    assert alternate.participant_id.iloc[0] == tree.participant_id.iloc[0]


def test_unverified_geo_case_label_does_not_establish_patient():
    result = resolve("GSE270638_MENINGIOMA_2024", ["I55", "I60"])
    assert result.participant_id.eq("").all()
    assert result.patient_identity_reason.eq("case_field_not_verified_as_patient").all()
    assert resolve("unknown_source", ["TH06_0613_S01"]).participant_id.iloc[0] == ""


def test_cll_source_exclusions_do_not_become_independent_patients():
    result = resolve("CLLMAP_2022", ["CRC-0001", "CRC-0007", "GCLL-0136"])
    assert result.participant_id.tolist() == ["CLLMAP_2022:CRC-0001", "", ""]
    assert result.patient_identity_references.iloc[0] == '["https://www.cllmap.org/downloads.html"]'


def test_derived_views_share_physical_and_patient_keys():
    a = resolve("TREEHOUSE_POLYA_25_01", ["TCGA-AB-1234-01A"])
    b = resolve("TREEHOUSE_POLYA_25_01_TCGA_SAMPLES", ["TCGA-AB-1234-01A"])
    assert a.sample_key.tolist() == b.sample_key.tolist()
    assert a.participant_id.tolist() == b.participant_id.tolist() == ["TCGA:TCGA-AB-1234"]


@pytest.mark.parametrize("source", ["TCGA", "TARGET", "MMRF"])
def test_project_name_alone_does_not_establish_a_barcode(source):
    assert resolve(source, ["sample1"]).participant_id.iloc[0] == ""


def test_conflicts_and_namespaces_are_preserved(monkeypatch):
    monkeypatch.setattr(
        samples,
        "sample_manifest",
        lambda: pd.DataFrame(
            [
                {
                    "source_cohort": "BEATAML_OHSU_2022",
                    "sample_id": "x",
                    "case_id": "1",
                    "included": True,
                },
            ]
        ),
    )
    monkeypatch.setattr(
        samples,
        "sample_molecular_provenance",
        lambda: pd.DataFrame(
            [
                {"source_cohort": "BEATAML_OHSU_2022", "sample_id": "x", "donor_id": "2"},
                {"source_cohort": "A", "sample_id": "x", "donor_id": "1"},
                {"source_cohort": "B", "sample_id": "x", "donor_id": "1"},
            ]
        ),
    )
    conflict = resolve("BEATAML_OHSU_2022", ["x"])
    assert conflict.participant_id.iloc[0] == ""
    assert conflict.patient_identity_status.iloc[0] == "conflict"
    assert (
        conflict.patient_identity_candidates.iloc[0]
        == '["BEATAML_OHSU_2022:1", "BEATAML_OHSU_2022:2"]'
    )
    assert resolve("A", ["x"]).participant_id.iloc[0] != resolve("B", ["x"]).participant_id.iloc[0]


def test_matrix_metadata_uses_actual_columns_and_identical_qc(monkeypatch):
    calls = []

    def matrix(code, **kwargs):
        calls.append((code, kwargs))
        return pd.DataFrame(
            {
                "Ensembl_Gene_ID": ["ENSG1"],
                "Symbol": ["G"],
                "TH06_0613_S02": [1],
                "TH06_0613_S01": [2],
            }
        )

    monkeypatch.setattr(expression, "per_sample_expression", matrix)
    monkeypatch.setattr(
        expression.source_matrices,
        "cohort_info",
        lambda code: {"source_cohort": "TREEHOUSE_POLYA_25_01"},
    )
    result = expression.sample_metadata("aml", auto_fetch=False, sample_qc="pass")
    assert calls == [("LAML", {"normalize": "tpm_raw", "auto_fetch": False, "sample_qc": "pass"})]
    assert result.sample_id.tolist() == ["TH06_0613_S02", "TH06_0613_S01"]
    assert result.participant_id.nunique() == 1
    assert result.attrs["sample_qc"] == "pass"


def test_empty_input_and_invalid_identifiers():
    assert resolve("A", []).empty
    with pytest.raises(ValueError, match="nonempty"):
        resolve("A", [None])

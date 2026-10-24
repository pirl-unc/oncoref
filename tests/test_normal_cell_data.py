"""Native evidence fidelity and the shared cache/API contract."""

import gzip
import xml.etree.ElementTree as ET

import pandas as pd
import pytest

from oncoref import catalog, reference_data
from oncoref import normal_cell_data as data


def test_all_new_sources_have_immutable_pins_and_catalog_routes():
    for name, spec in reference_data._normal_manifest["artifacts"].items():
        assert len(spec["sha256"]) == len(spec["archive_sha256"]) == 64
        assert spec["bytes"] > 0
        assert catalog.dataset(name).kind == spec["kind"]
        assert reference_data.resolve_version(name) == spec["version"]


def test_primary_ihc_preserves_uncertainty_without_inventing_donor_scores():
    entry = ET.fromstring("""<entry url="http://v23.proteinatlas.org/g">
      <name>G</name><identifier db="Ensembl" id="ENSG1"/>
      <tissueExpression technology="IHC" assayType="tissue"><verification type="reliability" description="Conflicting antibodies">uncertain</verification></tissueExpression>
      <antibody id="AB1"><antigenSequence>ACDE</antigenSequence>
        <tissueExpression technology="IHC" assayType="tissue"><verification type="validation">approved</verification>
          <data><tissue>Testis</tissue><tissueCell><cellType>Spermatocytes</cellType><level type="staining">high</level></tissueCell>
            <patient><patientId>P1</patientId></patient><patient><patientId>P2</patientId></patient>
          </data></tissueExpression></antibody></entry>""")
    rows = data._entry_ihc(entry)
    assert len(rows) == 1
    assert rows[0]["Profile reliability"] == "Uncertain"
    assert rows[0]["Antibody validation"] == "Approved"
    assert rows[0]["Donor IDs"] == "P1;P2"
    assert rows[0]["Level"] == "High"


def test_summary_restores_uncertain_profile_missing_from_xml():
    frame = pd.DataFrame({"Gene": ["G1", "G2"], "Profile reliability": ["Missing", "Enhanced"]})
    result = data._ihc_with_summary(frame, pd.Series({"G1": "Uncertain", "G2": "Supported"}))
    assert result["Profile reliability"].tolist() == ["Uncertain", "Enhanced"]
    assert result["XML profile reliability"].tolist() == ["Missing", "Enhanced"]
    assert result["Profile reliability source"].tolist() == ["hpa_gene_summary", "hpa_antibody_xml"]


def test_dice_values_are_not_matched_donor_records(monkeypatch):
    native = pd.DataFrame(
        {
            "Feature_name": ["ENSG1.2", "ENSG2.1"],
            "Transcript_Length(bp)": [100, 200],
            "Additional_annotations": ["G1;protein_coding", "G2;protein_coding"],
            "0": [8.0, 20.0],
            "1": [0.0, 3.0],
        }
    )
    monkeypatch.setattr(data, "table", lambda *a, **k: native.copy())
    result = data.observations("dice_cd4_stim")
    assert result.unit.eq("TPM").all() and result.state.eq("activated").all()
    assert result.cell_label.str.contains("T cell").all()
    assert result.donor.eq("").all() and result["sample"].eq("").all()
    assert set(result.gene_id) == {"ENSG1", "ENSG2"}
    assert sorted(result.value) == [0, 3, 8, 20]
    assert len(data._gene_filter(native, {"ENSG1"})) == 1
    assert len(data._gene_filter(native, {"G2"})) == 1


def test_cornea_normalization_uses_all_genes_and_includes_zero_cells(tmp_path, monkeypatch):
    path = tmp_path / "counts.csv.gz"
    with gzip.open(path, "wt") as f:
        f.write('"","cell_1","cell_2"\n"A",2,0\n"B",8,20\n')
    monkeypatch.setattr(reference_data, "ensure", lambda *a, **k: path)
    monkeypatch.setattr(reference_data, "provenance", lambda *a, **k: {"sha256": "a" * 64})
    monkeypatch.setattr(
        data,
        "cornea_cells",
        lambda: pd.DataFrame(
            {
                "cell": ["cell_1", "cell_2"],
                "cell_type": ["limbal cell", "limbal cell"],
                "sample_group": ["1", "2"],
            }
        ),
    )
    result = data.cornea_summary(genes=["A"])
    assert result.value.tolist() == [2000.0, 0.0]
    assert result.positive_cells.tolist() == [1, 0]
    assert result.unit.eq("CP10K").all()
    assert result.donor.eq("").all()
    assert data.cornea_summary(genes=["absent"]).empty


def test_cornea_annotation_mismatch_is_not_silently_dropped(tmp_path, monkeypatch):
    path = tmp_path / "counts.csv.gz"
    with gzip.open(path, "wt") as f:
        f.write('"","cell_1"\n"A",2\n')
    monkeypatch.setattr(reference_data, "ensure", lambda *a, **k: path)
    monkeypatch.setattr(reference_data, "provenance", lambda *a, **k: {"sha256": "b" * 64})
    monkeypatch.setattr(
        data,
        "cornea_cells",
        lambda: pd.DataFrame(
            {"cell": ["other_1"], "cell_type": ["limbal cell"], "sample_group": ["1"]}
        ),
    )
    with pytest.raises(ValueError, match="do not match"):
        data.cornea_summary()

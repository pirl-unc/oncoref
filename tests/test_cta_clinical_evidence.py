"""Clinical nominations retain locus resolution and independent safety decisions."""

import gzip
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from cta_clinical_curation_report import clinical_disposition

from oncoref.gene_ids import canonical_gene_space
from oncoref.load_dataset import get_data


def test_clinical_and_safety_records_have_canonical_coding_identities():
    space = canonical_gene_space().set_index("ensembl_gene_id")
    for name in ["cta-trial-evidence", "cta-external-safety-reviews"]:
        rows = get_data(name)
        assert rows.Ensembl_Gene_ID.is_unique
        for row in rows.itertuples():
            assert space.loc[row.Ensembl_Gene_ID, "symbol"] == row.Symbol
            assert space.loc[row.Ensembl_Gene_ID, "biotype"] == "protein_coding"
            assert row.source_url.startswith("https://")
            assert row.limitations


def test_registry_intervention_evidence_is_pinned_and_distinguished_from_family_labels():
    rows = get_data("cta-trial-evidence").set_index("Symbol")
    path = Path(__file__).resolve().parents[1] / "docs/audits/cta-curation-20261021"
    records = json.loads(gzip.decompress((path / "trial-registry-snapshot.json.gz").read_bytes()))
    ids = {r["identificationModule"]["nctId"] for r in records}
    for source in rows.source_anchor:
        if source.startswith("NCT"):
            assert source in ids
    assert rows.loc["ROPN1", "evidence_kind"] == "family_label_clinical"
    assert rows.loc["ROPN1B", "evidence_kind"] == "family_label_clinical"
    assert "shared" in rows.loc["CTAG1A", "target_resolution"]
    assert "shared" in rows.loc["CTAG1B", "target_resolution"]
    assert rows.loc["CLDN6", "scope"] == "oncofetal"
    assert rows.loc["MAGEA12", "evidence_kind"] == "clinical_safety_context"


def test_normal_brain_hold_is_distinct_from_receptor_cross_reactivity_review():
    rows = get_data("cta-external-safety-reviews").set_index("Symbol")
    assert rows.loc["MAGEA12", "action"] == "hold"
    assert rows.loc["MAGEA3", "action"] == "product_review"


@pytest.mark.parametrize(
    "hold,retained,reasons,expected",
    [
        (True, False, "", "external_safety_hold"),
        (False, True, "", "retained"),
        (
            False,
            False,
            "positive_ihc_outside_permitted;rna_above_tissue_ceiling",
            "normal_protein_counterevidence",
        ),
        (False, False, "rna_above_tissue_ceiling", "rna_threshold_review"),
        (False, False, "reproductive_selectivity_below_threshold", "selectivity_reference_review"),
        (False, False, "incomplete_routine_rna", "coverage_or_curation_review"),
    ],
)
def test_clinical_triage_preserves_distinct_evidence_strengths(hold, retained, reasons, expected):
    row = SimpleNamespace(external_safety_hold=hold, retained=retained, failure_reasons=reasons)
    assert clinical_disposition(row) == expected

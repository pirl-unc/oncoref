"""Do not turn missing cells, cross-reactive probes or weak qPCR into absence."""

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from oncoref import cancer_only, gene_knowledge

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/audits/cancer-only-20261109"


def test_entire_measured_zero_reproductive_subset_is_reviewed_once():
    baseline = pd.read_csv(
        ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz",
        low_memory=False,
    )
    expected = baseline[baseline.tumor_restricted_profile & baseline.reproductive_ntpm.eq(0)]
    result = cancer_only.records()
    assert len(result) == 47 and result.gene_id.is_unique
    assert set(result.entity_id) == set(expected.protein_equivalence_group_id)
    assert result.reproductive_bulk_ntpm.eq(0).all()
    assert result.reproductive_cell_ntpm.isna().all()  # never filled with zero
    assert result.developmental_max_tpm.gt(0).all()
    assert not result.confirmed_cancer_only.any()
    assert set(result.classification) <= {"Other:Review", "Other:Expression", "Unsafe"}
    assert result.notes.map(bool).all() and result.sources.map(bool).all()
    for observations in result.tumor_observations:
        assert len(observations) == len({o["sample_key"] for o in observations}) == 2
        assert all(o["unit"] == "TPM" and o["value"] >= 10 for o in observations)
    assert cancer_only.audit()["counts"]["numerical_screen"] == 66


def test_copy_cutoff_is_not_cq_cutoff_and_assays_are_not_duplicated_donors():
    table = pd.read_csv(REPORT / "olfactory-assays.csv")
    r = table[table.protein.eq("OR5D18")].iloc[0]
    assert r.donors_cq_lt35 == 7 and r.donors_ge5_copies == 0
    assert r.cross_reactive_partner == "OR5D16"
    assert len(table[table.protein.eq("OR6V1")]) == 2
    assert table.donors_ge5_copies.dropna().between(0, 26).all()
    assert table[table.protein.eq("OR2B2")].assay_id.isna().all()
    assert table[table.protein.eq("OR52E8")].cross_reactive_partner.eq("OR52E6").all()


def test_normal_assays_create_review_without_inventing_protein_or_unsafe():
    table = cancer_only.records().set_index("protein")
    for symbol in [
        "OR1L4",
        "OR5C1",
        "OR52E8",
        "OR5D18",
        "TAS2R42",
        "USP17L8",
        "USP17L13",
        "OR11H1",
        "OR51B6",
    ]:
        assert table.loc[symbol, "classification"] == "Other:Review"
        records = gene_knowledge.findings(symbol, include_updates=True)
        new = records[records.id.str.startswith("GK-20261109-")]
        assert not new.empty and new.effect.eq("review").all()
        assert not new.tissue.eq("source-specified normal tissue").any()
    assert table.loc["IFNA4", "classification"] == table.loc["IFNW1", "classification"] == "Unsafe"
    assert table.hpa_ms.eq("Not available").sum() == 46
    assert table.loc["KRTAP21-1", "hpa_ms"] == "Evidence at protein level"
    assert table.loc["OR8S1", "reported_fetal_cells"] == 1
    placental = gene_knowledge.findings("CLDN25", include_updates=True)
    placental = placental[placental.id.eq("GK-20261109-CLDN25-HPA")].iloc[0]
    assert placental.stage == "adult and/or placental; see tissue"
    assert placental.effect == "review" and "Melanoma" not in placental.tissue
    assert all(
        v["cell"] != "Tumor cells" for entries in table.hpa_positive_staining for v in entries
    )
    assert all(a["reliability"] == "uncertain" for a in table.loc["IFNA4", "hpa_antibodies"])


def test_accessor_and_receipts_preserve_exact_gene_and_source_identity():
    result = cancer_only.records("OR5D18")
    pd.testing.assert_frame_equal(result, cancer_only.records(result.iloc[0].gene_id))
    assert cancer_only.records("not-a-gene").empty
    with pytest.raises(ValueError, match="nonempty"):
        cancer_only.records("")
    manifest = json.loads((REPORT / "manifest.json").read_text())
    for section, base in [("inputs", ROOT), ("outputs", REPORT)]:
        for name, digest in manifest[section].items():
            assert hashlib.sha256((base / name).read_bytes()).hexdigest() == digest, name

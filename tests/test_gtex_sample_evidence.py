import gzip
import sys
from pathlib import Path

import pandas as pd
import pytest

from oncoref import reference_data

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))
import gtex_sample_evidence as gtex


def test_versioned_rows_and_donors(tmp_path):
    matrix = tmp_path / "tissue.gz"
    with gzip.open(matrix, "wt") as handle:
        handle.write("#1.2\n3 2\nName\tDescription\tGTEX-A-1\tGTEX-A-2\n")
        handle.write("g.2\tgene\t2\t9\ng.1\tgene\t7\t4\nh.1\tother\t0\tnan\n")
    values, donors = gtex.sample_matrix(matrix, {b"g.2": "g", b"h.1": "h"})
    assert list(donors) == ["GTEX-A", "GTEX-A"]
    assert values.loc["g"].tolist() == [2, 9]
    assert pd.isna(values.loc["h", "GTEX-A-2"])
    with pytest.raises(ValueError, match="Duplicate GTEx locus"):
        gtex.sample_matrix(matrix, {b"g.2": "g", b"g.1": "g"})


def test_catalog_contains_all_primary_tissue_shards_and_the_pilot():
    pins = reference_data.REFERENCE_SOURCES[gtex.SOURCE]["pins"]
    assert len([v for v in pins if v.startswith("v11-")]) == 52
    assert not any("cells_" in v for v in pins)
    assert (
        pins["v11-testis"]["sha256"]
        == "d2f63d00c9cfd7151e108d6521b47142435eef6fc79371ec6a637e9ce23130fd"
    )
    assert all(pin["unit"] == "TPM" and len(pin["sha256"]) == 64 for pin in pins.values())


def test_live_gtex_witnesses_match_the_symmetric_ten_protein_audit():
    from oncoref import target_selection

    audits = Path(__file__).parents[1] / "docs/audits"
    rebuild = audits / "target-selection-20261110/rebuild"
    leads = pd.read_csv(audits / "gtex-reproductive-witnesses/symmetric-comparison.csv")
    catalog = target_selection.results().set_index("protein_group_label").loc[leads.protein]
    witnesses = pd.read_csv(
        rebuild / "reproductive-witnesses.csv.gz",
        usecols=[
            "entity_id",
            "tissue",
            "source_id",
            "participant_id",
            "value",
        ],
    )
    testis = witnesses[witnesses.source_id.eq(gtex.SOURCE_ID) & witnesses.tissue.eq("testis")]
    donors = testis.groupby("entity_id").participant_id.nunique()
    assert catalog.entity_id.map(donors).tolist() == leads.testis_donors_at_tenfold.tolist()
    bounds = catalog[["normal_reference_tpm", "disposable_max_tpm"]].max(axis=1).clip(lower=1)
    assert leads.threshold_tpm.tolist() == pytest.approx(bounds * target_selection.ENRICHMENT_FOLD)
    assert leads.classification.tolist() == catalog.classification.tolist()
    assert leads.classification.value_counts().to_dict() == {
        "CTA": 10,
    }

    qc = pd.read_csv(rebuild / "gtex-library-qc.csv.gz", keep_default_na=False)
    assert qc.sample_id.is_unique and len(qc) == 18_637
    assert qc.participant_id.nunique() == 946 and qc.source_tissue.nunique() == 52


def test_site_p99_uses_one_value_per_measured_donor():
    values = pd.DataFrame([[0, 2, 8], [None, None, None], [None, None, 3]], index=["a", "b", "c"])
    result = gtex.donor_p99(values, pd.Index(["first", "first", "second"]))
    assert result.to_dict() == pytest.approx({"a": 7.93, "c": 3.0})
    with pytest.raises(ValueError, match="donor"):
        gtex.donor_p99(values, pd.Index(["first", None, "second"]))


def test_live_normal_reference_matches_independent_p99_sensitivity():
    from oncoref import target_selection

    audit = Path(__file__).parents[1] / "docs/audits"
    expected = pd.read_csv(
        audit / "gtex-cross-tissue-signals/protein-sensitivity.csv.gz"
    ).set_index("entity_id")
    actual = target_selection.results().set_index("entity_id").loc[expected.index]
    assert actual[["normal_reference_tpm", "disposable_max_tpm"]].max(
        axis=1
    ).tolist() == pytest.approx(expected.ceiling_p99.tolist())
    frame = pd.read_csv(
        audit / "target-selection-20261110/rebuild/normal-comparison.csv.gz",
        usecols=["source_id", "scope", "sample_id", "participant_id", "source_context"],
    )
    p99 = frame[frame.source_id.eq("gtex_donor_p99:v11")]
    assert p99.source_context.nunique() == 43
    assert set(p99.scope) <= {"normal", "disposable"}
    assert p99[["sample_id", "participant_id"]].isna().all().all()
    assert set(frame.loc[frame.source_id.eq(gtex.SOURCE_ID), "scope"]) == {"reproductive"}

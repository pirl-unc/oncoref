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


def test_site_p95_uses_one_value_per_measured_donor():
    values = pd.DataFrame([[0, 2, 8], [None, None, None], [None, None, 3]], index=["a", "b", "c"])
    result = gtex.donor_percentile(values, pd.Index(["first", "first", "second"]))
    assert result.to_dict() == pytest.approx({"a": 7.65, "c": 3.0})
    with pytest.raises(ValueError, match="donor"):
        gtex.donor_percentile(values, pd.Index(["first", None, "second"]))


def test_live_ceiling_reproduces_accepted_normal_measurements():
    from oncoref import target_selection

    audit = Path(__file__).parents[1] / "docs/audits"
    frame = pd.read_csv(
        audit / "target-selection-20261110/rebuild/normal-comparison.csv.gz",
        usecols=[
            "entity_id",
            "value",
            "source_id",
            "scope",
            "sample_id",
            "participant_id",
            "source_context",
            "qc_pass",
            "identity_resolved",
        ],
    )
    accepted = frame[
        frame.qc_pass & frame.identity_resolved & frame.scope.isin(["normal", "disposable"])
    ]
    expected = accepted.groupby("entity_id").value.max()
    actual = target_selection.results().set_index("entity_id").loc[expected.index]
    assert actual[["normal_reference_tpm", "disposable_max_tpm"]].max(
        axis=1
    ).tolist() == pytest.approx(expected.tolist())
    p95 = frame[frame.source_id.eq("gtex_donor_p95:v11")]
    assert p95.source_context.nunique() == 43
    assert set(p95.scope) <= {"normal", "disposable"}
    assert p95[["sample_id", "participant_id"]].isna().all().all()
    assert set(frame.loc[frame.source_id.eq(gtex.SOURCE_ID), "scope"]) == {"reproductive"}


def test_published_markers_need_a_complete_coordinated_foreign_panel_and_peers():
    from oncoref import normal_library_qc as qc

    markers = [gene for _, symbols in qc.CARRYOVER_PANELS.values() for gene in symbols]
    values = pd.DataFrame(0.0, index=markers, columns=["a", "b", "c", "isolated", "coordinated"])
    values.loc["PRSS1", "isolated"] = 5
    values["coordinated"] = 10
    reference = pd.DataFrame({"code": ["P"] * 4 + ["E"] * 2, "scale": 1000.0}, index=markers)
    flags = qc.published_carryover(values, "liver", reference)
    assert flags.all(axis=1).tolist() == [False, False, False, False, True]
    assert not flags.loc["isolated"].any()
    assert (
        qc.published_carryover(values, "Esophagus_Muscularis", reference).carryover_E.isna().all()
    )
    assert (
        qc.published_carryover(values.drop(index="KRT4"), "liver", reference)
        .carryover_E.isna()
        .all()
    )
    assert qc.published_carryover(values.iloc[:, :2], "liver", reference).isna().all().all()
    assert not qc.published_carryover(values * 0.01, "liver", reference).any().any()


def test_reproductive_witness_anatomy_and_mammary_state_match_source_libraries():
    import json

    from oncoref import host_context

    rebuild = Path(__file__).parents[1] / "docs/audits/target-selection-20261110/rebuild"
    witnesses = pd.read_csv(rebuild / "reproductive-witnesses.csv.gz")
    individual = witnesses[witnesses.source_id.eq(gtex.SOURCE_ID)]
    uterus = individual[individual.source_context.eq("Uterus")]
    assert len(uterus) and set(uterus.tissue) == set(uterus.compartment) == {"uterus"}
    assert set(host_context.HostContext(protect=("uterus",)).protected) == {"uterus", "endometrium"}
    assert set(host_context.HostContext(protect=("endometrium",)).protected) == {"endometrium"}
    breast = individual[individual.tissue.eq("breast")]
    audit = pd.read_csv(rebuild / "gtex-library-qc.csv.gz").set_index("sample_id")
    assert len(breast) and breast.sample_id.map(audit.mammary_state_inference).equals(
        breast.mammary_state_inference
    )
    for record in breast.mammary_state_inference.unique():
        state = json.loads(record)
        assert state["provider_confirmed"] is False and len(state["marker_tpm"]) == 4


def test_cross_source_qc_counts_preserve_assessment_coverage():
    rebuild = Path(__file__).parents[1] / "docs/audits/target-selection-20261110/rebuild"
    report = pd.read_csv(rebuild / "normal-rna/carryover-summary.csv").set_index("source")
    for source, filename in [
        ("GTEx v11 tissue", "gtex-library-qc.csv.gz"),
        ("HPA v23 tissue", "normal-library-qc.csv"),
        ("HPA v23 sorted immune", "normal-rna/immune-library-qc.csv"),
    ]:
        data = pd.read_csv(rebuild / filename)
        assert data.sample_id.is_unique
        flags = data[["carryover_P", "carryover_E"]].astype("boolean")
        row = report.loc[source]
        assert (row.libraries, row.assessed, row.either) == (
            len(flags),
            flags.notna().any(axis=1).sum(),
            flags.any(axis=1).sum(),
        )
        for code in ("P", "E"):
            assert row[code] == flags["carryover_" + code].sum()
            assert row[code + "_assessed"] == flags["carryover_" + code].count()
    assert 0 < report.at["HPA v23 tissue", "assessed"] < report.at["HPA v23 tissue", "libraries"]

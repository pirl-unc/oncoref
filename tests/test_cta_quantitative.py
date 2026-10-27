"""Common measurements, not familiarity or publication counts, govern admission."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from oncoref import cta_quantitative as policy

AUDIT = Path(__file__).resolve().parents[1] / "docs/audits/cta-quantitative-20261026"


def facts(**changes):
    row = {
        "gene_id": "g1",
        "symbol": "UNFAMILIAR",
        "taxonomy_id": 9606,
        "reference": "ref",
        "reproductive_ntpm": 10,
        "adult_bulk_max_ntpm": 2,
        "developmental_mean_tpm": np.nan,
        "developmental_cells": np.nan,
        "developmental_positive_cells": np.nan,
        "developmental_fraction": np.nan,
        "normal_bulk_tissues": 50,
        "cell_records": 426,
        "cell_supported_records": 0,
        "cell_max_ntpm": np.nan,
        "hpa_immune_cell_records": 19,
        "hpa_immune_cell_max_nTPM": 0,
        "hpa_immune_monaco_records": 30,
        "hpa_immune_monaco_max_TPM": 0,
        "dice_sources": 15,
        "dice_third_max_tpm": 0,
        "cornea_records": 84,
        "cornea_second_max_cp10k": np.nan,
        "tumor_floor_ptpm": 1,
        "tumor_minimum_samples": 3,
        "tumor_minimum_prevalence": 0.05,
        "tumor_positive": True,
        "ihc_detected": False,
    }
    row.update(changes)
    return pd.DataFrame([row])


def selected(**changes):
    return bool(policy.classify(facts(**changes)).selected.iloc[0])


def test_identity_and_publication_fields_cannot_admit_failed_measurements():
    assert selected()
    assert not selected(symbol="PRAME", reproductive_ntpm=0, publications=999, reviewed=True)
    assert selected(symbol="UNFAMILIAR", publications=0, reviewed=False)
    assert not selected(symbol="TBXT", reproductive_ntpm=0, developmental_mean_tpm=100)


def test_development_requires_supported_lineage_but_not_invented_embryo_replication():
    d = facts(
        reproductive_ntpm=0,
        developmental_mean_tpm=5,
        developmental_cells=20,
        developmental_positive_cells=5,
        developmental_fraction=0.25,
        developmental_embryos=1,
    )
    r = policy.classify(d)
    assert r.selected.iloc[0] and r.developmental_embryos.iloc[0] == 1
    d["developmental_positive_cells"] = 4
    assert not policy.classify(d).selected.iloc[0]


@pytest.mark.parametrize(
    "changes",
    [
        {"normal_bulk_tissues": 49},
        {"cell_records": 425},
        {"hpa_immune_cell_records": 18},
        {"hpa_immune_monaco_records": 29},
        {"dice_sources": 14},
        {"cornea_records": np.nan},
        {"cell_supported_records": 1, "cell_max_ntpm": np.nan},
        {"hpa_immune_cell_max_nTPM": np.nan},
    ],
)
def test_missing_coverage_is_not_absence(changes):
    assert not selected(**changes)


def test_no_supported_signal_is_distinct_from_missing_source_records():
    assert selected(cell_supported_records=0, cell_max_ntpm=np.nan)
    assert not selected(cell_records=np.nan, cell_supported_records=0, cell_max_ntpm=np.nan)


@pytest.mark.parametrize(
    "changes",
    [
        {"adult_bulk_max_ntpm": 10.1},
        {"cell_max_ntpm": 30, "cell_supported_records": 1},
        {"hpa_immune_cell_max_nTPM": 5},
        {"hpa_immune_monaco_max_TPM": 5},
        {"dice_third_max_tpm": 5},
        {"cornea_second_max_cp10k": 1},
        {"ihc_detected": True},
        {"tumor_positive": False},
    ],
)
def test_each_adult_and_tumor_gate_applies_without_overrides(changes):
    assert not selected(**changes)


def test_reproductive_maximum_cannot_hide_absolute_normal_expression():
    assert not selected(reproductive_ntpm=1000000, adult_bulk_max_ntpm=11)
    assert not selected(reproductive_ntpm=1000000, cell_max_ntpm=30, cell_supported_records=1)
    assert not selected(reproductive_ntpm=5, adult_bulk_max_ntpm=5)
    assert selected(
        reproductive_ntpm=11, adult_bulk_max_ntpm=10, cell_max_ntpm=29.9, cell_supported_records=1
    )


def test_gene_wise_annotations_do_not_change_selection():
    d = policy.results()
    expected = policy.classify(d).selected
    d["symbol"] = "RENAMED"
    d["gene_id"] = ["other-" + str(i) for i in range(len(d))]
    d["Expression category"] = "Unsafe"
    d["Clinical qualification"] = True
    d["Cancer-germline papers"] = 0
    assert policy.classify(d).selected.equals(expected)


def test_frozen_reprocessing_and_controls():
    d = policy.results()
    assert len(d) == 20097 and d.selected.sum() == 447
    assert policy.classify(d).selected.equals(d.selected)
    s = d[d.selected]
    expected = {
        "PRAME",
        "MAGEA4",
        "CTAG1A",
        "CTAG1B",
        "CTAG2",
        "XAGE1A",
        "XAGE1B",
        "SSX2",
        "CT83",
        "TBXT",
    }
    assert expected <= set(s.symbol)
    assert not {"MAGEA3", "MAGEA12", "AFP", "CLDN6", "GPC3", "WT1", "BIRC5", "EGFR", "ERBB2"} & set(
        s.symbol
    )
    assert s["Expression category"].value_counts().to_dict() == {
        "Neither": 210,
        "CTA": 127,
        "Unsafe": 110,
    }
    assert s.ihc_missing.sum() == 271
    assert s.rna_coverage_complete.all()


def test_audit_artifact_hashes_and_published_result_match():
    manifest = json.loads((AUDIT / "run-manifest.json").read_text())
    for filename, expected in manifest["outputs"].items():
        assert hashlib.sha256((AUDIT / filename).read_bytes()).hexdigest() == expected
    published = Path(policy.__file__).parent / "data/cta-quantitative-gene-sets.csv.gz"
    assert published.read_bytes() == (AUDIT / "gene-results.csv.gz").read_bytes()


@pytest.mark.parametrize(
    "changes", [{"taxonomy_id": 10090}, {"tumor_floor_ptpm": 2}, {"cell_max_ntpm": -1}]
)
def test_bad_reference_units_or_summary_contract_rejected(changes):
    with pytest.raises(ValueError):
        policy.classify(facts(**changes))

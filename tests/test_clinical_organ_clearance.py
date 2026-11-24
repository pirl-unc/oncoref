# Licensed under the Apache License, Version 2.0.
"""Regression tests for #766 clinical pMHC-directed organ-safety tier & CTAG2 isoform contract."""

from __future__ import annotations

import json
from copy import deepcopy

import pytest

from oncoref import (
    antigen_sets,
    clinical_organ_clearance,
    human_screening_inputs,
    sequence_bundle,
    target_selection,
)


def test_ctag2_isoform_specificity_lage1a_vs_lage1b_and_tsarina_columns():
    doc = clinical_organ_clearance.ledger()
    tbl = clinical_organ_clearance.clearances_table().set_index("clearance_id")

    # CTAG2 lete-cel clearance covers ENST00000369585 (LAGE-1a, 180 aa, SLLMWITQC)
    # and explicitly excludes ENST00000247306 (LAGE-1b, 210 aa, lacks SLLMWITQC)
    ctag2_row = tbl.loc["COC-20261123-CTAG2-LETE-CEL"]
    assert ctag2_row["status"] == "cleared"
    assert ctag2_row["clearance_scope"] == "epitope"
    assert json.loads(ctag2_row["clinical_cleared_epitopes"]) == ["SLLMWITQC"]
    assert json.loads(ctag2_row["clinical_cleared_transcripts"]) == ["ENST00000369585"]
    assert json.loads(ctag2_row["uncleared_coding_transcripts"]) == ["ENST00000247306"]
    assert json.loads(ctag2_row["clinical_cleared_organs"]) == ["lung", "urinary bladder"]
    assert "HLA-A*02:01" in json.loads(ctag2_row["clinical_hla_alleles"])

    # Structured tsarina columns exposed directly on target_selection.results() and sequence_bundle
    ts_ctag2 = target_selection.results("CTAG2").iloc[0]
    for col in clinical_organ_clearance.CLINICAL_CLEARANCE_COLUMNS:
        assert col in ts_ctag2.index
    assert ts_ctag2["clinical_clearance_scope"] == "epitope"
    assert json.loads(ts_ctag2["clinical_cleared_transcripts"]) == ["ENST00000369585"]
    assert json.loads(ts_ctag2["clinical_cleared_organs"]) == ["lung", "urinary bladder"]

    disp = sequence_bundle.target_selection_dispositions().set_index("protein_group_label")
    assert disp.loc["CTAG2", "clinical_clearance_scope"] == "epitope"
    assert json.loads(disp.loc["CTAG2", "clinical_cleared_epitopes"]) == ["SLLMWITQC"]

    cta_prots = antigen_sets.proteins("CTA").set_index("protein_group_label")
    assert json.loads(cta_prots.loc["NY-ESO-1", "clinical_cleared_epitopes"]) == ["SLLMWITQC"]
    assert json.loads(cta_prots.loc["PRAME", "clinical_cleared_epitopes"]) == ["SLLQHLIGL"]
    assert json.loads(cta_prots.loc["MAGEA4", "clinical_cleared_epitopes"]) == ["GVYDGREHTV"]

    # Xq28 #711 attribution recorded alongside clinical clearance
    xq28 = clinical_organ_clearance.ctag2_xq28_attribution()
    assert xq28["selected_transcript"] == "ENST00000369585"
    assert xq28["uncleared_isoform"] == "ENST00000247306"
    assert xq28["external_unsafe_overlap_pct"] < 1.0
    assert doc["schema"] == "oncoref.clinical-organ-clearances"


def test_bbb_barrier_vaccine_modality_and_toxicity_negative_controls_fail_closed():
    base = clinical_organ_clearance.ledger()

    # 1. Attempting to clear a BBB-protected organ (brain) via systemic IV fails closed
    bad_bbb = deepcopy(base)
    bad_bbb["entries"][0]["cleared_organs"] = ["lung", "brain"]
    with pytest.raises(ValueError, match="barrier-protected"):
        clinical_organ_clearance.validate(bad_bbb)

    # 2. Attempting to claim blanket gene-level clearance fails closed
    bad_scope = deepcopy(base)
    bad_scope["entries"][0]["clearance_scope"] = "gene"
    with pytest.raises(ValueError, match="epitope"):
        clinical_organ_clearance.validate(bad_scope)

    # 3. Attempting to use a peptide vaccine modality fails closed
    bad_mod = deepcopy(base)
    bad_mod["entries"][0]["modality"] = "peptide_vaccine"
    with pytest.raises(ValueError, match=r"modality|TCR-T"):
        clinical_organ_clearance.validate(bad_mod)

    # 4. Attempting to claim heart clearance on lete_cel_c259 (recorded day-161 cardiac event) fails closed
    bad_heart = deepcopy(base)
    bad_heart["entries"][0]["cleared_organs"] = ["lung", "heart"]
    with pytest.raises(ValueError, match="heart"):
        clinical_organ_clearance.validate(bad_heart)

    # 5. Negative toxicity controls (MAGEA3/titin, MAGEA12 neuro, TBI-1301) cannot list cleared_organs
    bad_tox = deepcopy(base)
    neg = next(e for e in bad_tox["entries"] if e["product_id"] == "magea3_a2_neuro")
    neg["cleared_organs"] = ["brain"]
    with pytest.raises(ValueError, match="Blocked entry"):
        clinical_organ_clearance.validate(bad_tox)

    # 6. Unreviewed HLA allele or versioned/malformed ENST ID fails closed
    bad_hla = deepcopy(base)
    bad_hla["entries"][0]["hla_alleles"] = ["HLA-B*07:02"]
    with pytest.raises(ValueError, match="hla_alleles"):
        clinical_organ_clearance.validate(bad_hla)

    bad_tx = deepcopy(base)
    bad_tx["entries"][0]["cleared_transcripts"] = ["ENST00000369585.1"]
    with pytest.raises(ValueError, match="ENST"):
        clinical_organ_clearance.validate(bad_tx)


def test_ctag2_advances_to_cta_under_clinical_organ_clearance():
    res = human_screening_inputs.evaluate_human_catalog(
        include_nuclear_spillover_qc=True,
        include_clinical_clearances=True,
    )
    cat = res.catalog.set_index("protein_group_label")
    assert cat.loc["CTAG2", "classification"] == "CTA"
    assert float(cat.loc["CTAG2", "normal_reference_tpm"]) == pytest.approx(6.92457, rel=1e-4)
    assert float(cat.loc["CTAG2", "tumor_max_tpm"]) == pytest.approx(604.332, rel=1e-4)
    assert cat.loc["CTAG2", "clinical_clearance_scope"] == "epitope"
    assert json.loads(cat.loc["CTAG2", "clinical_cleared_transcripts"]) == ["ENST00000369585"]
    # Negative controls remain strictly Unsafe
    assert cat.loc["MAGEA3", "classification"] == "Unsafe"
    assert cat.loc["MAGEA12", "classification"] == "Unsafe"

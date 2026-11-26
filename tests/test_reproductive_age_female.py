"""Tests for reproductive-age female GTEx stratification & GSE107746 adult oocyte/follicle RNA."""

from __future__ import annotations

import pytest

from oncoref import (
    antigen_sets,
    gtex,
    human,
    human_screening_inputs,
    reproductive_age_female,
    target_selection,
)


def test_gtex_donor_age_stratifies_reproductive_vs_postmenopausal_xx_donors():
    all_ages = gtex.donor_age()
    assert len(all_ages) >= 940
    assert {
        "donor_id",
        "reported_sex",
        "chromosomal_sex",
        "age_bracket",
        "age_midpoint",
        "is_reproductive_age",
        "hardy_scale",
    } <= set(all_ages.columns)
    assert all_ages["donor_id"].is_unique

    xx_repro = gtex.donor_age(reproductive_age_only=True, sex="XX")
    assert 100 <= len(xx_repro) <= 130
    assert (xx_repro["chromosomal_sex"] == "XX").all()
    assert xx_repro["is_reproductive_age"].all()
    assert (xx_repro["age_midpoint"] <= 49).all()
    assert set(xx_repro["donor_id"]) == reproductive_age_female.reproductive_age_female_donor_ids()

    xx_all = gtex.donor_age(sex="XX")
    assert len(xx_all) > len(xx_repro)
    assert int((~xx_all["is_reproductive_age"]).sum()) >= 190

    with pytest.raises(ValueError):
        gtex.donor_age(sex="INVALID")


def test_gse107746_adult_oocyte_follicle_panel_sha256_and_paralog_guards():
    receipt = reproductive_age_female.adult_oocyte_follicle_receipt()
    assert receipt["panel_id"] == reproductive_age_female.ADULT_FOLLICLE_PANEL_ID
    assert receipt["study_accession"] == "GSE107746"
    assert receipt["library_qc"]["libraries_total"] == 151
    assert receipt["library_qc"]["libraries_qc_pass"] == 128
    assert receipt["identity_mapping"]["paralog_kmer_ceiling"] == 0.05

    df = reproductive_age_female.adult_oocyte_follicle_rna()
    assert len(df) == receipt["identity_mapping"]["reported_rows_ge_0_5_tpm"]
    assert df["entity_id"].is_unique
    by_sym = df.set_index("protein_group_label")

    # Genuine adult human oocyte transcripts quantified at high stage-replicated TPM
    for ooc_gene in (
        "H1-8",
        "PRR23E",
        "MOS",
        "CPSF4L",
        "BMP15",
        "IZUMO1R",
        "ACCSL",
        "WEE2",
        "ZP4",
        "ZAR1L",
        "KHDC3L",
    ):
        assert ooc_gene in by_sym.index
        assert float(by_sym.loc[ooc_gene, "oocyte_supported_tpm"]) >= 25.0

    # Paralog-ambiguous / host-overlapped genes (FUT5, OCM2, NLRP8, DPPA3) are excluded from Cufflinks FPKM
    for excluded_paralog in ("FUT5", "OCM2", "NLRP8", "DPPA3"):
        assert excluded_paralog not in by_sym.index


def test_reproductive_age_female_promotions_and_single_rule_ablation():
    cat = target_selection.results().set_index("protein_group_label")
    sex_df = human.sex_restricted_antigens().set_index("protein_group_label")

    # MOS promoted CA only -> CTA-MF (734 TPM adult oocytes + 10.8 nTPM spermatocytes);
    # CPSF4L promoted CA only -> CTA-F (266 TPM adult oocytes, male <= 3.7 TPM)
    assert cat.loc["MOS", "classification"] == "CTA"
    assert sex_df.loc["MOS", "sex_restriction_code"] == "CTA-MF"
    assert cat.loc["CPSF4L", "classification"] == "CTA"
    assert sex_df.loc["CPSF4L", "sex_restriction_code"] == "CTA-F"

    # Adult oocyte expression expands CTA-F from 8 to 12 (with OR1N1 held in RA-F-only by the
    # monogenic chemosensory GPCR guard) and moves BMP15 + OR52D1 + OR1N1 into RA-F-only
    cta_f = set(antigen_sets.proteins("CTA-F")["protein_group_label"])
    assert len(cta_f) == 12
    assert {
        "CSN2",
        "LALBA",
        "DPPA3",
        "KHDC3L",
        "ZAR1L",
        "ZP4",
        "VRTN",
        "PRR23E",
        "CPSF4L",
        "IZUMO1R",
        "ACCSL",
        "TBXT",
    } == cta_f
    assert sex_df.loc["OR1N1", "sex_restriction_code"] == "RA-F-only"
    assert sex_df.loc["BMP15", "sex_restriction_code"] == "RA-F-only"
    assert sex_df.loc["OR52D1", "sex_restriction_code"] == "RA-F-only"
    assert sex_df.loc["CIMIP2A", "sex_restriction_code"] == "CTA-MF"
    assert sex_df.loc["H1-8", "sex_restriction_code"] == "RA-MF-only"

    # Orthogonal single-rule ablation turning off include_reproductive_age_female (under active GPCR guard)
    ablated = human_screening_inputs.evaluate_human_catalog(
        include_reproductive_age_female=False
    ).catalog.set_index("protein_group_label")
    assert ablated.loc["MOS", "classification"] == "Cancer Antigen only"
    assert ablated.loc["CPSF4L", "classification"] == "Cancer Antigen only"
    assert ablated.loc["DCAF4L2", "sex_restriction_code"] == "CTA-M"
    assert ablated.loc["SPATA3", "sex_restriction_code"] == "RA-M-only"
    assert ablated.loc["DEFB123", "sex_restriction_code"] == "CTA-MF"
    assert ablated.loc["TBXT", "sex_restriction_code"] == "CTA-P"
    assert int(ablated["category"].eq("CTA").sum()) == 320

    # Turning off both reproductive-age female and chemosensory GPCR guard restores exact 1.51.0 counts (323 CTA)
    ablated_151 = human_screening_inputs.evaluate_human_catalog(
        include_reproductive_age_female=False,
        include_chemosensory_gpcr_guard=False,
    ).catalog.set_index("protein_group_label")
    assert int(ablated_151["category"].eq("CTA").sum()) == 323
    assert ablated_151.loc[
        ablated_151["sex_restriction_code"].ne(""), "sex_restriction_code"
    ].value_counts().to_dict() == {
        "CTA-M": 236,
        "RA-M-only": 123,
        "CTA-P": 52,
        "CTA-MF": 27,
        "RA-P-only": 13,
        "CTA-F": 8,
        "RA-MF-only": 3,
        "RA-F-only": 3,
    }

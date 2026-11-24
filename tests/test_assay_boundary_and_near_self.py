"""Tests for assay-boundary adjudications (#768) and 1-mismatch Hamming / tryptic-MS audit (#769)."""

from __future__ import annotations

import pandas as pd
import pytest

from oncoref import (
    antigen_sets,
    assay_boundary_adjudications,
    human,
    peptide_subtraction,
    sequence_bundle,
    target_selection,
)


@pytest.fixture(scope="module")
def integrated_cat() -> pd.DataFrame:
    return human.evaluate_integrated_catalog(
        include_home_panels=True,
        include_adjudications=True,
        include_737_cohorts=True,
        include_multiplex_bleed_qc=True,
        include_nuclear_spillover_qc=True,
        include_clinical_clearances=True,
        include_repeat_family_equivalence=True,
        include_assay_boundary_fixes=True,
    ).catalog.set_index("protein_group_label")


def test_713_non_admixed_tumor_witness_moves_11_spermatid_genes_while_keeping_prm1_and_spata19_cta(
    integrated_cat: pd.DataFrame,
):
    # Without #737 cSCC, SPATA19 lacks a >=10x non-TGCT witness and dynamically moves to RA only
    cat_no_737 = human.evaluate_integrated_catalog(
        include_multiplex_bleed_qc=True,
        include_assay_boundary_fixes=True,
        include_737_cohorts=False,
    ).catalog.set_index("protein_group_label")
    assert cat_no_737.loc["SPATA19", "classification"] == "Reproductive Antigen only"
    assert (
        cat_no_737.loc["SPATA19", "assay_boundary_768_disposition"]
        == "reclassified_ra_only_admixed_tgct_witness"
    )

    # With full integration (#737 cSCC GSM3567709 = 13.28 TPM >= 10.0 TPM), both PRM1 and SPATA19 stay CTA
    assert integrated_cat.loc["PRM1", "classification"] == "CTA"
    assert integrated_cat.loc["SPATA19", "classification"] == "CTA"
    assert integrated_cat.loc["SPATA19", "tumor_top_cancer_code"] == "cSCC"
    assert (
        integrated_cat.loc["SPATA19", "assay_boundary_768_disposition"]
        == "retained_cta_independent_somatic_witness"
    )

    # PRM2, CETN1, PGK2 and the remaining 8 AAFL sole-witness spermatid markers (11 total) move to RA only
    for ra_lbl in (
        "PRM2",
        "CETN1",
        "PGK2",
        "BOD1L2",
        "CCDC70",
        "TSSK2",
        "PRR30",
        "ZNRF4",
        "IQCF5",
        "BPIFA3",
        "CYLC2",
    ):
        assert integrated_cat.loc[ra_lbl, "classification"] == "Reproductive Antigen only"
        assert (
            integrated_cat.loc[ra_lbl, "assay_boundary_768_disposition"]
            == "reclassified_ra_only_admixed_tgct_witness"
        )


def test_mageb6_ihc_dkkl1_cornea_and_foxr2_spata31a3_germline_tf_rules(
    integrated_cat: pd.DataFrame,
):
    # MAGEB6 (RNA ~0, complete RNA coverage, unreliable single-antibody IHC) -> CTA
    assert integrated_cat.loc["MAGEB6", "classification"] == "CTA"
    assert integrated_cat.loc["MAGEB6", "assay_boundary_768_disposition"] == "ihc_rna_zero_cleared"
    # NLRP9 (somatic RNA 2.32 TPM >= 0.5) and XAGE2 (other Unsafe gates) remain Unsafe
    assert integrated_cat.loc["NLRP9", "classification"] == "Unsafe"
    assert integrated_cat.loc["XAGE2", "classification"] == "Unsafe"

    # DKKL1 (only 8 multi-UMI cells across 21,343 adult cornea cells) -> CTA
    assert integrated_cat.loc["DKKL1", "classification"] == "CTA"
    assert (
        integrated_cat.loc["DKKL1", "assay_boundary_768_disposition"]
        == "adult_cornea_sparse_umi_cleared"
    )
    # Multi-UMI corneal genes (>= 24 cells with >=2 UMIs) stay Unsafe
    for unsafe_cornea in ("INSL4", "GPHB5", "QRFP", "C20orf202"):
        assert integrated_cat.loc[unsafe_cornea, "classification"] == "Unsafe"

    # Protamine-diluted germline TFs with >=2 GTEx Testis donors and germline scRNA nTPM >= 1.0 -> CTA
    for tf_lbl in ("FOXR2", "SPATA31A3", "PPP4R3C"):
        assert integrated_cat.loc[tf_lbl, "classification"] == "CTA"
        assert integrated_cat.loc[tf_lbl, "sex_restriction_code"] == "CTA-M"
        assert (
            integrated_cat.loc[tf_lbl, "assay_boundary_768_disposition"]
            == "germline_tf_protamine_diluted_qualified"
        )
    # Single-donor DEFB125 (1 donor < 2) fails the tsarina donor-floor guard and stays CA only
    assert integrated_cat.loc["DEFB125", "classification"] == "Cancer Antigen only"


def test_near_self_1mm_hamming_policy_and_ctage1_mia2_tryptic_ms_audit(
    integrated_cat: pd.DataFrame,
):
    # Verify PeptideIndex + hamming_1mm_policy() gives pure 1-substitution Hamming without indels
    records = pd.DataFrame(
        [
            {"gene_id": "ENSG_UNSAFE", "sequence": "MKVLWAALLVTFLAGCQAKVEQAVETEPEPELR"},
        ]
    )
    idx_1mm = peptide_subtraction.PeptideIndex(
        records,
        policy=assay_boundary_adjudications.hamming_1mm_policy(),
        lengths=(9,),
    )
    assert len(idx_1mm.matches("KVLWAALLY")) == 1
    assert idx_1mm.matches("KVLWAALYY").empty
    assert idx_1mm.matches("KVLWALLVT").empty

    # Verify tsarina 1-mm columns across target_selection, antigen_sets, and sequence_bundle
    ts = target_selection.results().set_index("protein_group_label")
    disp = sequence_bundle.target_selection_dispositions().set_index("protein_group_label")
    ctas = antigen_sets.proteins("CTA").set_index("protein_group_label")

    for col in assay_boundary_adjudications.NEAR_SELF_1MM_COLUMNS:
        assert col in ts.columns
        assert col in disp.columns
        assert col in ctas.columns

    # Promoted former-Unsafe/Review CTAs (PLAC1, ZAN, MAGEB6, DKKL1, CTAG2, DUX4, FOXR2)
    # now have positive exact-unique & >=2-mm safe 9-mers directly on target_selection.results(),
    # antigen_sets.proteins("CTA"), and sequence_bundle.target_selection_dispositions()!
    for cleared_cta in ("PLAC1", "ZAN", "MAGEB6", "DKKL1", "CTAG2", "DUX4", "FOXR2"):
        assert int(ts.loc[cleared_cta, "exact_unique_9mers"]) > 0
        assert int(ts.loc[cleared_cta, "safe_2mm_9mers"]) > 0
        assert int(ctas.loc[cleared_cta, "exact_unique_9mers"]) == int(
            integrated_cat.loc[cleared_cta, "exact_unique_9mers"]
        )
        assert int(disp.loc[cleared_cta, "safe_2mm_9mers"]) == int(
            integrated_cat.loc[cleared_cta, "safe_2mm_9mers"]
        )

    # Newly Unsafe #735/#711 proteins (SP7, DAZL) report 0 unique 9-mers on both ts and integrated_cat
    assert int(ts.loc["SP7", "exact_unique_9mers"]) == 0
    assert int(ts.loc["DAZL", "exact_unique_9mers"]) == 0
    assert int(integrated_cat.loc["SP7", "exact_unique_9mers"]) == 0
    assert int(integrated_cat.loc["DAZL", "exact_unique_9mers"]) == 0

    # MAGEA6: 86 exact-unique 9-mers, of which >80% are 1-mm from Unsafe MAGEA3/MAGEA12 (16 safe >=2-mm)
    assert int(ts.loc["MAGEA6", "exact_unique_9mers"]) == 86
    assert int(ts.loc["MAGEA6", "safe_2mm_9mers"]) == 16
    assert float(ts.loc["MAGEA6", "unsafe_near_self_1mm_pct"]) > 80.0
    assert float(ts.loc["MAGEA2/B", "unsafe_near_self_1mm_pct"]) > 60.0

    # All 12 audited proteins (including PLAC1) are present in tryptic_ms_paralog_audit()
    ms_df = assay_boundary_adjudications.tryptic_ms_paralog_audit().set_index("protein_group_label")
    assert len(ms_df) == 12
    assert "PLAC1" in ms_df.index
    assert int(ms_df.loc["CTAGE1", "unsafe_shared_tryptic_peptides"]) >= 30
    assert "MIA2" in str(ms_df.loc["CTAGE1", "top_unsafe_ms_paralogs"])
    assert "TGM2" in str(ms_df.loc["TGM4", "top_unsafe_ms_paralogs"])
    assert "SPATA31C" in str(ms_df.loc["SPATA31A3", "top_unsafe_ms_paralogs"])

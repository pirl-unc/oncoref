# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Tests for the TCGA-TGCT somatic testicular admixture & AAFL witness audit (#713)."""

from __future__ import annotations

import oncoref
from oncoref import target_selection, tgct_localization


def test_tgct_reproductive_localization_audit_covers_all_154_libraries():
    full = oncoref.tgct_reproductive_localization_audit()
    assert len(full) == 154
    assert full["sample_id"].is_unique

    r_only = oncoref.tgct_reproductive_localization_audit(r_only=True)
    assert len(r_only) == 7
    assert r_only.iloc[0]["sample_id"] == "TCGA-2G-AAFL-01"
    aafl = r_only.iloc[0]
    assert aafl["testis_tissue_content"] > 0.04
    assert aafl["spermatid_gtex_log10_r"] > 0.98
    assert aafl["INSL3_tpm"] > 100.0
    assert aafl["CYP17A1_tpm"] > 20.0
    assert aafl["CLU_tpm"] > 1500.0
    assert aafl["INHA_tpm"] > 40.0
    assert aafl["inferred_malignant_program"] == "embryonal_carcinoma_dominant"
    assert aafl["localization_call"] == "admixed_seminiferous_and_interstitial_parenchyma"


def test_tgct_aafl_protein_support_partitions_all_34_aafl_peaking_catalog_proteins():
    ts = target_selection.results()
    expected_ids = set(
        ts.loc[
            ts["tumor_top_sample"].astype(str).str.contains("TCGA-2G-AAFL-01", na=False),
            "entity_id",
        ].astype(str)
    )
    support = oncoref.tgct_aafl_protein_support()
    assert len(support) == 34
    assert set(support["entity_id"]) == expected_ids
    assert (support["reference_testis_tpm"] > 0).all()
    assert support["aafl_to_reference_testis_ratio"].notna().all()

    sole = oncoref.tgct_admixture_sole_witness_entity_ids()
    restricted = oncoref.tgct_admixture_restricted_entity_ids()
    assert len(sole) == 11
    assert len(restricted) == 13
    assert sole < restricted

    by_label = support.set_index("protein_group_label")
    assert by_label.loc["TEX38", "support_tier"] == "confirmed_outside_r_tgct"
    assert by_label.loc["TEX38", "clean_cancer_codes_at_threshold"] == "NEC_MERKEL"
    assert by_label.loc["PRM1", "support_tier"] == "r_tgct_multi_library_only"
    assert by_label.loc["PRM2", "support_tier"] == "r_tgct_multi_library_only"
    for sym in (
        "BOD1L2",
        "CETN1",
        "CCDC70",
        "TSSK2",
        "PGK2",
        "PRR30",
        "ZNRF4",
        "SPATA19",
        "IQCF5",
        "BPIFA3",
        "CYLC2",
    ):
        assert by_label.loc[sym, "support_tier"] == tgct_localization.AAFL_ADMIXED_SOLE_WITNESS_TIER
        assert int(by_label.loc[sym, "non_aafl_witnesses_at_threshold"]) == 0

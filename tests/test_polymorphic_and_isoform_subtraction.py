# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Tests for germline coding variant subtraction, isoform selectivity, and ncORF/TE contracts (#736)."""

from __future__ import annotations

import pandas as pd
import pytest

from oncoref import antigen_sets, sequence_bundle
from oncoref.peptide_subtraction import (
    CodingVariant,
    PeptideIndex,
    apply_coding_variants,
    near_self_contract,
)


def test_coding_variant_augments_background_and_subtracts_common_polymorphic_self_peptides():
    bg = pd.DataFrame(
        [
            {
                "gene_id": "ENSG_SELF",
                "transcript_id": "ENST_SELF.3",
                "sequence": "MKLPVRGSACDEFGHIKLMNPQRSTVWYA",
            }
        ]
    )
    base_index = PeptideIndex(bg, lengths=[9])
    # Reference 9mer around 1-based pos 13 ('F' -> 'Y'): CDEFGHIKL vs CDEYGHIKL
    assert not base_index.matches("CDEFGHIKL").empty
    assert base_index.matches("CDEYGHIKL").empty

    common_var = CodingVariant(
        target_id="ENST_SELF",
        position=13,
        ref_aa="F",
        alt_aa="Y",
        af=0.025,
        source_id="gnomAD:v4.1",
    )
    rare_var = CodingVariant(
        target_id="ENST_SELF",
        position=14,
        ref_aa="G",
        alt_aa="A",
        af=0.0001,
        source_id="gnomAD:v4.1",
    )

    poly_index = base_index.with_coding_variants([common_var, rare_var], min_af=0.001)
    hits = poly_index.matches("CDEYGHIKL")
    assert len(hits) == 1
    assert bool(hits.iloc[0]["is_coding_variant"]) is True
    assert hits.iloc[0]["variant_id"] == "ENST_SELF:p.F13Y"
    assert float(hits.iloc[0]["variant_af"]) == pytest.approx(0.025)
    assert hits.iloc[0]["variant_source_id"] == "gnomAD:v4.1"
    # Rare variant below min_af=0.001 is not added to the default background:
    assert poly_index.matches("CDEFAAIKL").empty
    # Candidate query carrying the polymorphic self 9mer is now subtracted:
    assert poly_index.subtract("CDEYGHIKL", lengths=[9]).empty
    assert "population_polymorphism_support" in near_self_contract()


def test_coding_variant_fails_closed_on_ref_mismatch_and_unknown_target():
    bg = pd.DataFrame([{"transcript_id": "ENST1", "sequence": "ACDEFGHIKLMN"}])
    wrong_ref = CodingVariant(
        target_id="ENST1",
        position=3,
        ref_aa="W",
        alt_aa="Y",
        af=0.05,
        source_id="gnomAD:v4.1",
    )
    with pytest.raises(ValueError, match="reference mismatch"):
        apply_coding_variants(bg, [wrong_ref])

    missing_target = CodingVariant(
        target_id="ENST_UNKNOWN",
        position=1,
        ref_aa="A",
        alt_aa="V",
        af=0.05,
        source_id="gnomAD:v4.1",
    )
    with pytest.raises(ValueError, match="matched no background"):
        apply_coding_variants(bg, [missing_target])


def test_evaluate_isoform_selectivity_preserves_locus_unsafe_precedence():
    # Canonical Unsafe isoform shares prefix ACDEFGHIKLMN; tumor splice isoform has a novel junction WWWWYYYYY
    unsafe_bg = pd.DataFrame(
        [{"gene_id": "ENSG_MIXED", "transcript_id": "ENST_CANON", "sequence": "ACDEFGHIKLMNPQRST"}]
    )
    unsafe_idx = PeptideIndex(unsafe_bg, lengths=[9])

    candidates = pd.DataFrame(
        [
            {
                "gene_id": "ENSG_MIXED",
                "transcript_id": "ENST_TUMOR_ISO.2",
                "sequence": "ACDEFGHIKWWWWYYYYY",
            },
            {
                "gene_id": "ENSG_MIXED",
                "transcript_id": "ENST_SUBSET_ISO.1",
                "sequence": "ACDEFGHIKLMN",
            },
        ]
    )

    # 1. Without transcript-level normal quantification -> held fail-closed!
    held = antigen_sets.evaluate_isoform_selectivity(
        candidates, unsafe_index=unsafe_idx, lengths=[9]
    )
    assert (~held["qualifies_isoform_selective"]).all()
    assert set(held["isoform_status"]) == {"held_missing_normal_isoform_quantification"}

    # 2. With normal isoform TPM <= 1.0 and tumor isoform TPM >= 10.0:
    # ENST_TUMOR_ISO qualifies ONLY on its junction peptides, whereas ENST_SUBSET_ISO has 0 unique 9mers!
    normal_tpm = pd.DataFrame(
        {
            "transcript_id": ["ENST_TUMOR_ISO", "ENST_SUBSET_ISO"],
            "heart": [0.0, 0.0],
            "lung": [0.2, 0.1],
        }
    )
    tumor_tpm = pd.DataFrame(
        {
            "transcript_id": ["ENST_TUMOR_ISO", "ENST_SUBSET_ISO"],
            "max_tumor_tpm": [25.0, 30.0],
        }
    )
    evaluated = antigen_sets.evaluate_isoform_selectivity(
        candidates,
        unsafe_index=unsafe_idx,
        normal_isoform_tpm=normal_tpm,
        tumor_isoform_tpm=tumor_tpm,
        lengths=[9],
    ).set_index("transcript_id")

    assert bool(evaluated.loc["ENST_TUMOR_ISO.2", "qualifies_isoform_selective"]) is True
    assert (
        evaluated.loc["ENST_TUMOR_ISO.2", "isoform_status"]
        == "qualifies_isoform_selective_peptides_only"
    )
    assert "WWWWYYYYY" in evaluated.loc["ENST_TUMOR_ISO.2", "retained_unique_peptides"]
    assert "ACDEFGHIK" not in evaluated.loc["ENST_TUMOR_ISO.2", "retained_unique_peptides"]

    assert bool(evaluated.loc["ENST_SUBSET_ISO.1", "qualifies_isoform_selective"]) is False
    assert (
        evaluated.loc["ENST_SUBSET_ISO.1", "isoform_status"] == "excluded_no_unsafe_unique_peptides"
    )
    assert antigen_sets.isoform_selectivity_contract()["contract_id"] == (
        antigen_sets.ISOFORM_SELECTIVITY_CONTRACT_ID
    )


def test_ncorf_te_screening_contract_exposes_all_four_mandatory_gates():
    contract = sequence_bundle.ncorf_te_screening_contract()
    assert contract["contract_id"] == sequence_bundle.NCORF_TE_CONTRACT_ID
    assert set(contract["mandatory_gates"]) == {
        "translation_evidence",
        "multicopy_mappability",
        "normal_atlas_and_mtec_screen",
        "self_and_polymorphic_subtraction",
    }

# Licensed under the Apache License, Version 2.0.
"""Unit and empirical regression tests for #764 stoichiometric multiplex-bleed QC."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oncoref import human_screening_inputs
from oncoref import normal_library_qc as qc


def _synthetic_sentinel_ref() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "gene_id": "G_PRM1",
                "symbol": "PRM1",
                "source_median_tpm": 20000.0,
                "source_key": "testis",
                "source_organ": "Testis",
            },
            {
                "gene_id": "G_PRM2",
                "symbol": "PRM2",
                "source_median_tpm": 20000.0,
                "source_key": "testis",
                "source_organ": "Testis",
            },
            {
                "gene_id": "G_TNP1",
                "symbol": "TNP1",
                "source_median_tpm": 8000.0,
                "source_key": "testis",
                "source_organ": "Testis",
            },
            {
                "gene_id": "G_SMCP",
                "symbol": "SMCP",
                "source_median_tpm": 2000.0,
                "source_key": "testis",
                "source_organ": "Testis",
            },
        ]
    )


def test_single_gene_or_disproportionate_expression_never_flagged():
    ref = _synthetic_sentinel_ref()
    # s_clean: all zero
    # s_single: genuine single-gene ectopic expression of PRM1 at 12 TPM without PRM2/TNP1
    # s_disprop: PRM1=400 TPM (2e-2), PRM2=0.6 TPM (3e-5), TNP1=0.6 TPM (7.5e-5) -> fails bleed band / log10 IQR
    # s_bleed: proportional ~2e-4 testis bleed across PRM1 (4.0), PRM2 (4.2), TNP1 (1.6), SMCP (0.5)
    mat = pd.DataFrame(
        {
            "s_clean": [0.0, 0.0, 0.0, 0.0],
            "s_single": [12.0, 0.0, 0.0, 0.0],
            "s_disprop": [400.0, 0.6, 0.6, 0.0],
            "s_bleed": [4.0, 4.2, 1.6, 0.5],
        },
        index=["G_PRM1", "G_PRM2", "G_TNP1", "G_SMCP"],
    )
    flagged = qc.detect_stoichiometric_multiplex_bleed(mat, site="Thyroid", sentinel_reference=ref)
    assert list(flagged["sample_id"]) == ["s_bleed"]
    assert flagged.iloc[0]["sentinels_detected"] == 4
    assert flagged.iloc[0]["log10_iqr"] <= 0.1

    # Self-organ comparison (Testis vs testis) is skipped
    assert qc.detect_stoichiometric_multiplex_bleed(
        mat, site="Testis", sentinel_reference=ref
    ).empty

    # min_sentinels < 3 is rejected to forbid single-gene masking
    with pytest.raises(ValueError, match="min_sentinels must be >= 3"):
        qc.detect_stoichiometric_multiplex_bleed(
            mat, site="Thyroid", sentinel_reference=ref, min_sentinels=2
        )


def test_orthogonal_assay_veto_preserves_samples_even_when_flagged():
    samples = [f"GTEX-D{i:03d}-01" for i in range(20)]
    donors = [f"GTEX-D{i:03d}" for i in range(20)]
    # 18 clean donors at 0.2 TPM, 2 bleed-contaminated donors at 10.0 TPM -> raw P95 > 9 TPM
    vals = pd.Series([0.2] * 18 + [10.0, 10.0], index=samples)
    flagged = {samples[18], samples[19]}

    # When orthogonal scRNA/HPA bulk is negative, contaminated samples are masked -> P95 == 0.2
    clean_p95 = qc.compute_bleed_cleaned_donor_summary(
        vals, donors, flagged_samples=flagged, orthogonal_supported=False
    )
    assert np.isclose(clean_p95, 0.2)

    # When orthogonal scRNA/HPA bulk supports target-organ expression, NO samples are masked!
    vetoed_p95 = qc.compute_bleed_cleaned_donor_summary(
        vals, donors, flagged_samples=flagged, orthogonal_supported=True
    )
    assert vetoed_p95 > 9.0


def test_packaged_gtex_multiplex_bleed_tables_and_prm1_somatic_witness():
    t_sum = qc.multiplex_bleed_tissue_summary()
    deltas = qc.multiplex_bleed_gene_deltas()
    receipt = qc.multiplex_bleed_receipt()

    assert not t_sum.empty
    assert not deltas.empty
    assert receipt["bleed_qc_applied_count"] == int(deltas["bleed_qc_applied"].sum())
    assert receipt["bleed_qc_applied_count"] >= 2

    # PRM1 and PRM2 (>50,000 TPM in testis) GTEx Thyroid bleed is cleaned while preserving HPA somatic floors
    by_label = deltas.set_index("protein_group_label")
    assert "PRM1" in by_label.index
    prm1 = by_label.loc["PRM1"]
    assert float(prm1["source_peak_tpm"]) >= 10000.0
    assert bool(prm1["bleed_qc_applied"]) is True
    assert float(prm1["baseline_normal_reference_tpm"]) > 7.5
    assert float(prm1["clean_normal_reference_tpm"]) < 6.0
    # Female BRCA somatic witness (58.97 TPM) now clears 10x the bleed-cleaned normal ceiling!
    assert float(prm1["max_non_tgct_tumor_tpm"]) > 50.0
    assert float(prm1["max_non_tgct_tumor_tpm"]) >= 10.0 * max(
        float(prm1["clean_normal_reference_tpm"]), 1.0
    )

    # Integrated human catalog evaluation runs cleanly with include_multiplex_bleed_qc=True
    res = human_screening_inputs.evaluate_human_catalog(
        include_home_panels=True,
        include_adjudications=True,
        include_737_cohorts=True,
        include_multiplex_bleed_qc=True,
    )
    prm1_eval = res.catalog.set_index("protein_group_label").loc["PRM1"]
    assert prm1_eval["classification"] == "CTA"
    assert prm1_eval["multiplex_bleed_764_disposition"] == "multiplex_bleed_cleaned"
    assert float(prm1_eval["normal_reference_tpm"]) < 6.0

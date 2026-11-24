from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oncoref import gtex_deconvolution, normal_library_qc


def test_parse_gtex_sample_barcode_extracts_donor_and_sm_plate():
    bc = gtex_deconvolution.parse_gtex_sample_barcode("GTEX-1HSEH-0826-SM-B2LXE")
    assert bc.donor_id == "GTEX-1HSEH"
    assert bc.site_code == "0826"
    assert bc.sm_plate_prefix == "B2LX"
    assert bc.sm_well_id == "B2LXE"


def test_infer_chromosomal_sex_is_splash_robust_at_donor_level():
    y_genes = list(gtex_deconvolution.Y_CHROMOSOME_SENTINEL_GENES)
    # Female donor GTEX-FEM1 has 5 clean XX shards and 1 thyroid shard splashed with 8 TPM Y-RNA
    mat = pd.DataFrame(
        {
            "GTEX-FEM1-0126-SM-AAA1": [0.0, 0.0, 0.0, 0.0],
            "GTEX-FEM1-0226-SM-AAA2": [0.1, 0.0, 0.0, 0.0],
            "GTEX-FEM1-0326-SM-AAA3": [0.0, 0.1, 0.0, 0.0],
            "GTEX-FEM1-0826-SM-B2LXE": [4.0, 2.0, 1.0, 1.0],  # 8 TPM Y splash -> ambiguous sample
            "GTEX-MAL1-0126-SM-BBB1": [35.0, 20.0, 12.0, 15.0],
        },
        index=y_genes,
    )
    s_sex, d_sex = gtex_deconvolution.infer_chromosomal_sex(mat)
    assert d_sex["GTEX-FEM1"] == "XX"
    assert d_sex["GTEX-MAL1"] == "XY"
    assert s_sex["GTEX-FEM1-0826-SM-B2LXE"] == "ambiguous"


def test_lineage_disjoint_anchor_selection_excludes_shared_corpus_luteum_steroidogenesis():
    # STAR and CYP17A1 have low median in Ovary (post-menopausal majority) but high P95
    # (pre-menopausal corpus luteum), whereas CYP11B1/CYP11B2/MC2R are strictly adrenal!
    med = pd.DataFrame(
        {
            "adrenal_gland": [5000.0, 4800.0, 1200.0, 3600.0, 4900.0, 1.0, 0.5, 0.2],
            "ovary": [0.5, 0.2, 0.1, 3.0, 2.5, 0.1, 0.1, 0.1],
            "thyroid": [0.1, 0.1, 0.0, 0.2, 0.1, 4000.0, 2500.0, 1800.0],
            "lung": [0.1, 0.1, 0.0, 45.0, 30.0, 0.2, 0.1, 0.1],
        },
        index=["CYP11B1", "CYP11B2", "MC2R", "STAR", "CYP17A1", "TG", "TPO", "TSHR"],
    )
    p95 = pd.DataFrame(
        {
            "adrenal_gland": [9000.0, 8500.0, 2400.0, 7000.0, 8800.0, 2.0, 1.0, 0.5],
            "ovary": [0.8, 0.4, 0.2, 1700.0, 510.0, 0.2, 0.2, 0.2],
            "thyroid": [0.2, 0.2, 0.1, 0.5, 0.4, 7500.0, 4600.0, 3100.0],
            "lung": [0.2, 0.2, 0.1, 80.0, 65.0, 0.4, 0.2, 0.2],
        },
        index=med.index,
    )

    anchors = gtex_deconvolution.select_lineage_disjoint_anchors(med, p95)
    adrenal_genes = set(anchors.loc[anchors["source_site"].eq("adrenal_gland"), "gene_id"])
    assert adrenal_genes == {"CYP11B1", "CYP11B2", "MC2R"}
    assert "STAR" not in adrenal_genes
    assert "CYP17A1" not in adrenal_genes


def test_nnls_archetype_mixture_recovers_true_weights_and_native_retention():
    basis = pd.DataFrame(
        {
            "ovary": [120.0, 95.0, 80.0, 0.0, 0.0, 0.0],
            "lung": [0.0, 0.0, 0.0, 3400.0, 1800.0, 1850.0],
        },
        index=["FOXL2", "INHA", "WNT2B", "SFTPB", "SFTPA1", "SFTPA2"],
    )
    samples = pd.DataFrame(
        {
            "clean_ovary": [118.0, 96.0, 79.0, 0.0, 0.0, 0.0],
            "lung_in_ovary_swap": [1.0, 0.5, 0.0, 3600.0, 1900.0, 1950.0],
            "mixed_90_10": [108.0, 85.5, 72.0, 340.0, 180.0, 185.0],
        },
        index=basis.index,
    )
    w = gtex_deconvolution.solve_nonnegative_archetype_mixture(basis, samples, l1_penalty=1e-4)
    assert w.at["ovary", "clean_ovary"] > 0.95
    assert w.at["lung", "clean_ovary"] < 0.01
    assert w.at["ovary", "lung_in_ovary_swap"] < 0.02
    assert w.at["lung", "lung_in_ovary_swap"] > 1.0
    assert pytest.approx(0.90, abs=0.02) == w.at["ovary", "mixed_90_10"]
    assert pytest.approx(0.10, abs=0.02) == w.at["lung", "mixed_90_10"]


def test_deconvolve_distinguishes_swap_plate_splash_same_donor_co_plated_boundary_and_ectopic():
    anchors = pd.DataFrame(
        [
            {
                "gene_id": "G_PRM1",
                "symbol": "PRM1",
                "source_site": "testis",
                "source_family": "testis",
                "source_median_tpm": 25000.0,
                "source_p95_tpm": 55000.0,
                "max_off_lineage_p95_tpm": 0.5,
                "excluded_target_sites": "",
            },
            {
                "gene_id": "G_PRM2",
                "symbol": "PRM2",
                "source_site": "testis",
                "source_family": "testis",
                "source_median_tpm": 24000.0,
                "source_p95_tpm": 50000.0,
                "max_off_lineage_p95_tpm": 0.5,
                "excluded_target_sites": "",
            },
            {
                "gene_id": "G_TNP1",
                "symbol": "TNP1",
                "source_site": "testis",
                "source_family": "testis",
                "source_median_tpm": 8000.0,
                "source_p95_tpm": 18000.0,
                "max_off_lineage_p95_tpm": 0.3,
                "excluded_target_sites": "",
            },
            {
                "gene_id": "G_SFTPB",
                "symbol": "SFTPB",
                "source_site": "lung",
                "source_family": "lung",
                "source_median_tpm": 3400.0,
                "source_p95_tpm": 6500.0,
                "max_off_lineage_p95_tpm": 0.8,
                "excluded_target_sites": "",
            },
            {
                "gene_id": "G_SFTPA1",
                "symbol": "SFTPA1",
                "source_site": "lung",
                "source_family": "lung",
                "source_median_tpm": 1800.0,
                "source_p95_tpm": 4000.0,
                "max_off_lineage_p95_tpm": 0.5,
                "excluded_target_sites": "",
            },
            {
                "gene_id": "G_SFTPA2",
                "symbol": "SFTPA2",
                "source_site": "lung",
                "source_family": "lung",
                "source_median_tpm": 1850.0,
                "source_p95_tpm": 4100.0,
                "max_off_lineage_p95_tpm": 0.5,
                "excluded_target_sites": "",
            },
            {
                "gene_id": "G_KLK3",
                "symbol": "KLK3",
                "source_site": "prostate",
                "source_family": "prostate",
                "source_median_tpm": 4000.0,
                "source_p95_tpm": 9000.0,
                "max_off_lineage_p95_tpm": 0.4,
                "excluded_target_sites": "",
            },
            {
                "gene_id": "G_MSMB",
                "symbol": "MSMB",
                "source_site": "prostate",
                "source_family": "prostate",
                "source_median_tpm": 3500.0,
                "source_p95_tpm": 8000.0,
                "max_off_lineage_p95_tpm": 0.4,
                "excluded_target_sites": "",
            },
            {
                "gene_id": "G_ACP3",
                "symbol": "ACP3",
                "source_site": "prostate",
                "source_family": "prostate",
                "source_median_tpm": 1600.0,
                "source_p95_tpm": 3800.0,
                "max_off_lineage_p95_tpm": 0.3,
                "excluded_target_sites": "",
            },
        ]
    )
    samples = {
        # Source testis well from DON1 on plate B2LX
        "GTEX-DON1-0126-SM-B2LX1": ("testis", [25000, 24000, 8000, 0, 0, 0, 0, 0, 0]),
        # DIFFERENT donor DON2 thyroid well on same plate B2LX with 5% Testis splash -> plate_pipetting_splash!
        "GTEX-DON2-0226-SM-B2LX2": ("thyroid", [1250, 1200, 400, 0, 0, 0, 0, 0, 0]),
        # SAME donor DON1 thyroid well on same plate B2LX with 2% Testis -> same_donor_carryover (NOT plate splash!)
        "GTEX-DON1-0226-SM-B2LX3": ("thyroid", [500, 480, 160, 0, 0, 0, 0, 0, 0]),
        # Clean thyroid wells + one 0.03% index-hopped thyroid well + one single-gene ectopic PRM1 well
        "GTEX-DON3-0226-SM-CLN01": ("thyroid", [7.5, 7.2, 2.4, 0, 0, 0, 0, 0, 0]),
        "GTEX-DON4-0226-SM-CLN02": ("thyroid", [50.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0]),
        **{
            f"GTEX-D{i:03d}-0226-SM-CLN{i:02d}": ("thyroid", [0, 0, 0, 0, 0, 0, 0, 0, 0])
            for i in range(10, 25)
        },
        # Ovary sample that is actually a 100% Lung tube swap!
        "GTEX-SWAP-2426-SM-9JGGK": ("ovary", [0, 0, 0, 3500, 1900, 1950, 0, 0, 0]),
        **{
            f"GTEX-OV{i:02d}-2426-SM-OVP{i:02d}": ("ovary", [0, 0, 0, 0, 0, 0, 0, 0, 0])
            for i in range(10)
        },
        # Bladder sample carrying anatomical boundary prostate signal (must NOT be subtracted!)
        "GTEX-BLD1-1526-SM-BLDP1": ("bladder", [0, 0, 0, 0, 0, 0, 80.0, 70.0, 32.0]),
        **{
            f"GTEX-BL{i:02d}-1526-SM-BLP{i:02d}": ("bladder", [0, 0, 0, 0, 0, 0, 0, 0, 0])
            for i in range(10)
        },
    }
    mat = pd.DataFrame(
        {sid: vals for sid, (_, vals) in samples.items()},
        index=anchors["gene_id"],
        dtype=float,
    )
    sites = pd.Series({sid: st for sid, (st, _) in samples.items()})

    events, summary = gtex_deconvolution.deconvolve_gtex_samples(mat, sites, anchors)
    ev_by_id = {(r.sample_id, r.source_site): r for r in events.itertuples(index=False)}

    assert ev_by_id[("GTEX-DON2-0226-SM-B2LX2", "testis")].regime == "plate_pipetting_splash"
    assert ev_by_id[("GTEX-DON2-0226-SM-B2LX2", "testis")].co_plated_on_sm_batch is True
    assert ev_by_id[("GTEX-DON1-0226-SM-B2LX3", "testis")].regime == "same_donor_carryover"
    assert ev_by_id[("GTEX-DON3-0226-SM-CLN01", "testis")].regime == "flowcell_index_hopping"
    assert ("GTEX-DON4-0226-SM-CLN02", "testis") not in ev_by_id  # single-gene ectopic preserved!
    assert ev_by_id[("GTEX-SWAP-2426-SM-9JGGK", "lung")].regime == "tube_swap"
    assert bool(summary.set_index("sample_id").loc["GTEX-SWAP-2426-SM-9JGGK", "quarantined_swap"])
    assert ev_by_id[("GTEX-BLD1-1526-SM-BLDP1", "prostate")].regime == "anatomical_boundary"

    # Verify summarize_deconvolved_donor_tpm drops splash, subtracts hopping, preserves ectopic & boundary!
    thy_samps = [s for s, (st, _) in samples.items() if st == "thyroid"]
    prm1_thy = mat.loc["G_PRM1", thy_samps].copy()
    prm1_thy["GTEX-DON4-0226-SM-CLN02"] = 0.0
    res_clean = gtex_deconvolution.summarize_deconvolved_donor_tpm(
        prm1_thy,
        site="thyroid",
        events=events,
        source_medians={"testis": 25000.0},
        source_p95s={"testis": 55000.0},
    )
    assert res_clean["raw_max"] == 1250.0
    assert res_clean["value"] < 0.5

    res_orth = gtex_deconvolution.summarize_deconvolved_donor_tpm(
        prm1_thy,
        site="thyroid",
        events=events,
        source_medians={"testis": 25000.0},
        source_p95s={"testis": 55000.0},
        orthogonal_supported=True,
    )
    assert res_orth["value"] > 10.0

    bld_samps = [s for s, (st, _) in samples.items() if st == "bladder"]
    klk3_bld = mat.loc["G_KLK3", bld_samps]
    res_bld = gtex_deconvolution.summarize_deconvolved_donor_tpm(
        klk3_bld,
        site="bladder",
        events=events,
        source_medians={"prostate": 4000.0},
        source_p95s={"prostate": 9000.0},
    )
    assert res_bld["dropped_samples"] == 0
    assert res_bld["value"] > 5.0


def test_harrell_davis_quantile_monotonicity_and_affine_equivariance():
    rng = np.random.default_rng(42)
    x = rng.lognormal(mean=0.5, sigma=0.8, size=18)
    w = rng.uniform(0.6, 1.0, size=18)
    q90 = gtex_deconvolution.purity_weighted_harrell_davis_quantile(x, w, 0.90)
    q95 = gtex_deconvolution.purity_weighted_harrell_davis_quantile(x, w, 0.95)
    q95_scaled = gtex_deconvolution.purity_weighted_harrell_davis_quantile(3.0 * x + 2.5, w, 0.95)
    assert q90 <= q95 <= float(x.max())
    assert pytest.approx(3.0 * q95 + 2.5, rel=1e-9) == q95_scaled


def test_packaged_gtex_deconvolution_tables_and_lung_in_ovary_swap():
    receipt = gtex_deconvolution.gtex_deconvolution_receipt()
    assert receipt["shards_evaluated"] == 52
    assert receipt["samples_evaluated"] > 18000
    assert receipt["anchor_source_archetypes"] >= 18
    for must_have in (
        "lung",
        "prostate",
        "thyroid",
        "adrenal_gland",
        "testis",
        "liver",
        "pancreas",
    ):
        assert must_have in set(receipt["anchor_source_sites"])

    swaps = gtex_deconvolution.gtex_quarantined_swaps()
    assert "GTEX-1GZHY-2426-SM-9JGGK" in set(swaps["sample_id"])
    ov_swap = swaps[swaps["sample_id"].eq("GTEX-1GZHY-2426-SM-9JGGK")].iloc[0]
    assert ov_swap["target_site"] == "ovary"
    assert ov_swap["source_site"] == "lung"
    assert ov_swap["donor_sex"] == "XX"
    assert ov_swap["sample_sex"] == "XY"
    assert ov_swap["alpha_hat"] > 0.80

    splashes = gtex_deconvolution.gtex_sample_splash_events(include_boundary=False)
    # GTEX-1HSEH had his own testis co-plated on SM-B2LX -> classified as same_donor_carryover!
    thy_b2lx = splashes[
        splashes["sample_id"].eq("GTEX-1HSEH-0826-SM-B2LXE") & splashes["source_site"].eq("testis")
    ]
    assert len(thy_b2lx) == 1
    assert thy_b2lx.iloc[0]["regime"] == "same_donor_carryover"

    # GTEX-1477Z-0226-SM-5TDCI (BP-48852 / SM-5TDC) is a verified cross-donor testis -> thyroid plate splash!
    thy_5tdc = splashes[
        splashes["sample_id"].eq("GTEX-1477Z-0226-SM-5TDCI") & splashes["source_site"].eq("testis")
    ]
    assert len(thy_5tdc) == 1
    assert thy_5tdc.iloc[0]["regime"] == "plate_pipetting_splash"
    assert bool(thy_5tdc.iloc[0]["co_plated_on_sm_batch"]) is True
    assert bool(thy_5tdc.iloc[0]["cross_donor_batch_matched"]) is True

    # Verify GTEX-1I1GR-1426-SM-CY8H8 (BP-72061 cross-donor extraction splash) is caught as plate_pipetting_splash!
    eso_cy8h = splashes[
        splashes["sample_id"].eq("GTEX-1I1GR-1426-SM-CY8H8") & splashes["source_site"].eq("testis")
    ]
    assert len(eso_cy8h) == 1
    assert eso_cy8h.iloc[0]["regime"] == "plate_pipetting_splash"
    assert bool(eso_cy8h.iloc[0]["cross_donor_batch_matched"]) is True

    # Verify same-donor surgical co-biopsies are NEVER mislabeled as plate_pipetting_splash!
    only_plate_splash = splashes[splashes["regime"].eq("plate_pipetting_splash")]
    assert not (
        only_plate_splash["sample_id"].str.startswith("GTEX-P4QT-")
        & only_plate_splash["target_site"].eq("adipose_subcutaneous")
    ).any()
    assert not (
        only_plate_splash["sample_id"].str.startswith("GTEX-11P82-")
        & only_plate_splash["target_site"].eq("nerve_tibial")
    ).any()

    summary = gtex_deconvolution.gtex_deconvolution_tissue_summary()
    assert len(summary) == 52
    assert (summary["median_native_retention"] > 0.4).all()
    # Testis should have <5% plate_splash_samples (no false-positive whole_blood or cerebellum splash!)
    testis_row = summary.set_index("site").loc["testis"]
    assert int(testis_row["plate_splash_samples"]) < 20
    assert int(splashes["source_site"].eq("whole_blood").sum()) == 0
    assert set(normal_library_qc.quarantined_gtex_swap_sample_ids()) == set(swaps["sample_id"])


def test_review_edge_cases_versioned_sex_nan_nnls_and_smooth_betainc_hd():
    # 1. Versioned Ensembl Y-chromosome IDs & all-NaN columns
    versioned_y = [f"{g}.12" for g in gtex_deconvolution.Y_CHROMOSOME_SENTINEL_GENES]
    mat = pd.DataFrame(
        {
            "GTEX-FEM2-0126-SM-A1": [0.0, 0.0, 0.0, 0.0],
            "GTEX-MAL2-0126-SM-B1": [40.0, 25.0, 15.0, 10.0],
            "GTEX-MISS-0126-SM-C1": [np.nan, np.nan, np.nan, np.nan],
        },
        index=versioned_y,
    )
    s_sex, d_sex = gtex_deconvolution.infer_chromosomal_sex(mat)
    assert s_sex["GTEX-FEM2-0126-SM-A1"] == "XX"
    assert s_sex["GTEX-MAL2-0126-SM-B1"] == "XY"
    assert s_sex["GTEX-MISS-0126-SM-C1"] == "unassessed"
    assert d_sex["GTEX-MISS"] == "unassessed"

    # 2. NaN in basis or sample never poisons NNLS mixture
    basis = pd.DataFrame(
        {"t1": [100.0, np.nan, 80.0], "t2": [0.0, 90.0, 0.0]}, index=["g1", "g2", "g3"]
    )
    samps = pd.DataFrame({"s1": [100.0, 0.0, np.nan]}, index=["g1", "g2", "g3"])
    w = gtex_deconvolution.solve_nonnegative_archetype_mixture(basis, samps, l1_penalty=1e-4)
    assert np.isfinite(w.to_numpy()).all()
    assert w.at["t1", "s1"] > 0.3

    # 3. Smooth Lentz betainc Harrell-Davis quantile across n=79 vs n=80 (no step discontinuity!)
    v79 = np.array([0.0] * 78 + [100.0])
    v80 = np.array([0.0] * 79 + [100.0])
    q79 = gtex_deconvolution.purity_weighted_harrell_davis_quantile(v79, quantile=0.95)
    q80 = gtex_deconvolution.purity_weighted_harrell_davis_quantile(v80, quantile=0.95)
    assert q79 > 1.0 and q80 > 1.0
    assert abs(q79 - q80) < 0.5

    # 4. summarize_deconvolved_donor_tpm never drops zeros or high endogenous housekeeping TPM
    evs = pd.DataFrame(
        [
            {
                "sample_id": "GTEX-D1-0126-SM-P1",
                "target_site": "testis",
                "source_site": "pancreas",
                "regime": "plate_pipetting_splash",
                "alpha_hat": 0.01,
                "alpha_lower_q25": 0.008,
            }
        ]
    )
    vals = pd.Series(
        {
            "GTEX-D1-0126-SM-P1": 1500.0,  # endogenous ACTB (1% of 500 = 5 TPM < 20% of 1500)
            **{f"GTEX-D{i}-0126-SM-P{i}": 1400.0 for i in range(2, 10)},
        }
    )
    res_actb = gtex_deconvolution.summarize_deconvolved_donor_tpm(
        vals, site="testis", events=evs, source_medians={"pancreas": 500.0}
    )
    assert res_actb["dropped_samples"] == 0

    zero_vals = pd.Series({f"GTEX-D{i}-0126-SM-P{i}": 0.0 for i in range(1, 10)})
    res_zero = gtex_deconvolution.summarize_deconvolved_donor_tpm(
        zero_vals, site="testis", events=evs, source_medians={"pancreas": 25000.0}
    )
    assert res_zero["dropped_samples"] == 0
    assert res_zero["value"] == 0.0

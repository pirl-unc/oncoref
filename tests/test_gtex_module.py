from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from oncoref import bulk_qc, gtex, gtex_deconvolution, human, normal_tissues, reference_data


def test_gtex_qc_cleanup_and_normalize_are_three_orthogonal_axes():
    raw = gtex.resolve_spec(qc=None, cleanup=None, normalize=None)
    assert raw.is_raw_default is True

    # Pure QC rejection (drop cell lines/swaps & mask splashes, ZERO model subtraction!)
    qc_only = gtex.resolve_spec(qc="all", cleanup=None)
    assert qc_only.qc == frozenset(gtex.ALL_QC)
    assert qc_only.cleanup == frozenset()
    assert qc_only.subtract_index_hopping is False
    assert qc_only.nnls_deconvolution is False

    # Explicit subsets on each axis independently
    sub = gtex.resolve_spec(qc=["swaps"], cleanup=["index_hopping"])
    assert sub.quarantine_swaps is True
    assert sub.reject_splash is False
    assert sub.subtract_index_hopping is True
    assert sub.nnls_deconvolution is False

    # Maximal QC + maximal cleanup + separate clean_tpm normalization
    full = gtex.resolve_spec(qc="all", cleanup="all", normalize="clean_tpm", sex="XX")
    assert full.qc == frozenset(gtex.ALL_QC)
    assert full.cleanup == frozenset(gtex.ALL_CLEANUP)
    assert full.clean_tpm is True
    assert full.sex == "XX"

    # "all" is a pure alias for the full list/set of operations (identical GtexSpec, hash, & cache key)
    explicit_list = gtex.resolve_spec(
        qc=["cell_lines", "swaps", "splash"],
        cleanup=["index_hopping", "deconvolution"],
        normalize="clean_tpm",
        sex="XX",
    )
    explicit_set = gtex.GtexSpec(
        qc={"splash", "cell_lines", "swaps"},
        cleanup={"deconvolution", "index_hopping"},
        normalize="clean_tpm",
        sex="XX",
    )
    direct_alias = gtex.GtexSpec(qc="all", cleanup=["all"], normalize="clean_tpm", sex="XX")
    assert full == explicit_list == explicit_set == direct_alias
    assert hash(full) == hash(explicit_list) == hash(explicit_set) == hash(direct_alias)
    assert gtex.QC_ALIASES["all"] == frozenset(gtex.ALL_QC)
    assert gtex.CLEANUP_ALIASES["all"] == frozenset(gtex.ALL_CLEANUP)
    assert gtex.expand_qc_ops("all") == frozenset({"cell_lines", "swaps", "splash"})
    assert gtex.expand_cleanup_ops({"all"}) == frozenset({"index_hopping", "deconvolution"})

    # Sentinel override behavior on existing GtexSpec: preserves unset fields, overrides explicit None/bool
    base_spec = gtex.GtexSpec(
        qc="all",
        cleanup="all",
        normalize="clean_tpm",
        sex="XX",
        include_subcompartments=False,
        hpa_unique_tissues_only=True,
    )
    unchanged = gtex.resolve_spec(qc=base_spec)
    assert unchanged == base_spec
    cleared = gtex.resolve_spec(
        qc=base_spec,
        cleanup=None,
        normalize=None,
        sex=None,
        include_subcompartments=True,
        hpa_unique_tissues_only=False,
    )
    assert cleared.cleanup == frozenset()
    assert cleared.normalize is None
    assert cleared.sex is None
    assert cleared.include_subcompartments is True
    assert cleared.hpa_unique_tissues_only is False

    k_alias, _ = gtex.spec_cache_key(
        provider="direct",
        version="v11",
        granularity="samples",
        stat="p95",
        value_column="TPM",
        spec=full,
        upstream_sha256={"s": "abc"},
    )
    k_list, _ = gtex.spec_cache_key(
        provider="direct",
        version="v11",
        granularity="samples",
        stat="p95",
        value_column="TPM",
        spec=explicit_list,
        upstream_sha256={"s": "abc"},
    )
    assert k_alias == k_list

    with pytest.raises(ValueError, match="Unknown GTEx qc"):
        gtex.resolve_spec(qc=["index_hopping"])  # index_hopping belongs to cleanup, not qc!
    with pytest.raises(ValueError, match="Unknown GTEx cleanup"):
        gtex.resolve_spec(cleanup=["swaps"])  # swaps belongs to qc, not cleanup!


def test_gtex_versions_and_tissues_across_direct_and_hpa():
    assert gtex.versions("direct", granularity="medians") == ("v11", "v10", "v8")
    assert gtex.versions("direct", granularity="samples") == ("v11",)
    assert gtex.versions("hpa") == ("v23",)

    t_v11 = gtex.tissues("v11", include_subcompartments=True, include_cell_lines=False)
    assert int(t_v11["kind"].eq("primary_shard").sum()) == 52
    assert int(t_v11["kind"].eq("v11_subcompartment_median").sum()) == 14
    assert "Pancreas_Islets" in set(t_v11["tissue"])

    t_hpa = gtex.tissues(provider="hpa")
    assert "retina" in set(t_hpa["tissue"])
    assert t_hpa.set_index("tissue").loc["retina", "kind"] == "hpa_unique_eyegex"


def test_apply_summary_spec_cleans_hpa_splash_preserves_nan_and_rejects_sex_filter(
    monkeypatch, tmp_path
):
    monkeypatch.setenv("CANCERDATA_DATA_DIR", str(tmp_path))
    anchors = gtex.gtex_deconvolution.gtex_deconvolution_anchors()
    testis_a = anchors[anchors["source_site"].eq("testis")].head(5)
    prost_a = anchors[anchors["source_site"].eq("prostate")].head(4)

    rows = []
    for r in testis_a.itertuples(index=False):
        m = float(r.source_median_tpm)
        rows.append(
            {
                "Name": r.gene_id,
                "Description": r.symbol,
                "testis": m,
                "prostate": 0.0,
                "breast": 0.0006 * m,
                "skeletal muscle": 0.0002 * m,
                "urinary bladder": 0.0,
                "retina": np.nan,
            }
        )
    for r in prost_a.itertuples(index=False):
        m = float(r.source_median_tpm)
        rows.append(
            {
                "Name": r.gene_id,
                "Description": r.symbol,
                "testis": 0.0,
                "prostate": m,
                "breast": 0.0,
                "skeletal muscle": 0.0,
                "urinary bladder": 0.02 * m,
                "retina": 0.0,
            }
        )
    syn_hpa = pd.DataFrame(rows)

    hop_cleaned = gtex.apply_summary_spec(
        syn_hpa,
        spec=gtex.resolve_spec(cleanup=["index_hopping"]),
        provider="hpa",
        value_column="nTPM",
    )
    prm_sym = testis_a.iloc[0]["symbol"]
    klk_sym = prost_a.iloc[0]["symbol"]
    by_sym = hop_cleaned.set_index("Description")
    assert by_sym.at[prm_sym, "breast"] < 0.5
    assert by_sym.at[prm_sym, "skeletal muscle"] < 0.5
    assert np.isnan(by_sym.at[prm_sym, "retina"])  # Pre-existing NaN preserved!
    assert by_sym.at[klk_sym, "urinary bladder"] == pytest.approx(
        syn_hpa.set_index("Description").at[klk_sym, "urinary bladder"]
    )

    retina_only = gtex.apply_summary_spec(
        syn_hpa,
        spec=gtex.resolve_spec(qc="all", cleanup="all", hpa_unique_tissues_only=True),
        provider="hpa",
        value_column="nTPM",
    )
    assert list(retina_only.columns) == ["Name", "Description", "retina"]

    with pytest.raises(ValueError, match="Donor sex filtering"):
        gtex.apply_summary_spec(syn_hpa, spec=gtex.GtexSpec(sex="XX"), provider="hpa")


def test_offline_synthetic_gct_shard_qc_cleanup_clean_tpm_and_cache(monkeypatch, tmp_path: Path):
    """Offline CI-safe end-to-end test for sample_tpm, tissue_summary, clean_tpm, and caching."""
    deriv_dir = tmp_path / "derived_gtex"
    monkeypatch.setattr(gtex, "derived_cache_dir", lambda: deriv_dir)

    anchors = gtex.gtex_deconvolution.gtex_deconvolution_anchors()
    testis_a = anchors[anchors["source_site"].eq("testis")].head(4)
    swap_sid = "GTEX-1GZHY-2426-SM-9JGGK"
    # Pick a real splashes-flagged thyroid sample from packaged events
    ev_all = gtex.gtex_deconvolution.gtex_sample_splash_events(include_boundary=False)
    splash_row = ev_all[
        ev_all["target_site"].eq("thyroid")
        & ev_all["source_site"].eq("testis")
        & ev_all["regime"].eq("plate_pipetting_splash")
    ].iloc[0]
    splash_sid = str(splash_row["sample_id"])

    dsex = gtex.donor_sex()
    xx_donors = dsex[dsex["chromosomal_sex"].eq("XX")]["donor_id"].head(6).tolist()
    xy_donors = dsex[dsex["chromosomal_sex"].eq("XY")]["donor_id"].head(3).tolist()
    xx_samps = [f"{d}-0226-SM-SYN{i:02d}" for i, d in enumerate(xx_donors)]
    xy_samps = [f"{d}-0226-SM-SYM{i:02d}" for i, d in enumerate(xy_donors)]
    all_samps = [swap_sid, splash_sid, *xx_samps, *xy_samps]

    # Build synthetic median GCT and synthetic thyroid shard GCT.gz (with ACTB, RPLP0, FOXL2, & testis anchors)
    gene_rows = [
        ("ENSG00000075624.17", "ACTB", 800.0, 900.0),  # biological housekeeping
        ("ENSG00000089157.16", "RPLP0", 2500.0, 2200.0),  # censored ribosomal gene for clean_tpm!
        ("ENSG00000183770.6", "FOXL2", 0.0, 50.0),  # zero in testis, 50 in thyroid
        ("ENSG00000124107.6", "SLPI", 0.0, 0.0),  # true zero in both!
    ]
    for r in testis_a.itertuples(index=False):
        gene_rows.append((f"{r.gene_id}.4", str(r.symbol), float(r.source_median_tpm), 0.0))

    med_lines = [
        "#1.2",
        f"{len(gene_rows)}\t4",
        "Name\tDescription\tTestis\tThyroid\tPancreas_Islets\tCells_Cultured_fibroblasts",
    ]
    for gid_v, sym, t_med, thy_med in gene_rows:
        med_lines.append(f"{gid_v}\t{sym}\t{t_med}\t{thy_med}\t10.0\t5.0")
    med_path = tmp_path / "syn_medians.gct"
    med_path.write_text("\n".join(med_lines) + "\n", encoding="utf-8")

    shard_lines = [
        "#1.2",
        f"{len(gene_rows)}\t{len(all_samps)}",
        "\t".join(["Name", "Description", *all_samps]),
    ]
    for gid_v, sym, t_med, thy_med in gene_rows:
        vals = []
        for s in all_samps:
            if s == swap_sid:
                vals.append(9999.0)
            elif s == splash_sid:
                # Splash sample has endogenous ACTB=1200 (must NOT be masked!), true zero SLPI=0.0 (must stay 0.0!),
                # and splashed testis anchors at 5% of testis median (MUST be masked to NaN!)
                if sym == "ACTB":
                    vals.append(1200.0)
                elif sym == "SLPI":
                    vals.append(0.0)
                elif t_med > 100.0 and thy_med == 0.0:
                    vals.append(0.05 * t_med)
                else:
                    vals.append(thy_med if thy_med > 0 else 100.0)
            elif s == xx_samps[0] and t_med > 100.0 and thy_med == 0.0:
                # 0.04% index hopping in an XX donor (<0.5% band)
                vals.append(0.0004 * t_med)
            else:
                vals.append(thy_med)
        shard_lines.append("\t".join([gid_v, sym, *(f"{v:.4f}" for v in vals)]))

    shard_path = tmp_path / "syn_thyroid.gct.gz"
    with gzip.open(shard_path, "wt", encoding="utf-8") as fh:
        fh.write("\n".join(shard_lines) + "\n")

    def fake_ensure(source_name: str, ver: str | None = None) -> Path:
        if source_name == gtex.SOURCE_ID:
            return med_path
        if source_name == gtex.SHARD_SOURCE_ID:
            return shard_path
        raise AssertionError(f"Unexpected ensure({source_name!r}, {ver!r})")

    monkeypatch.setattr(reference_data, "ensure", fake_ensure)

    # 1. qc="all", cleanup=None: drops swap, masks splashed testis anchor to NaN, preserves ACTB=1200 & SLPI=0.0!
    qc_only = gtex.sample_tpm("thyroid", qc="all", cleanup=None)
    assert swap_sid not in qc_only.columns
    by_sym = qc_only.set_index("Description")
    assert by_sym.at["ACTB", splash_sid] == pytest.approx(1200.0)
    assert by_sym.at["SLPI", splash_sid] == 0.0
    assert np.isnan(by_sym.at[testis_a.iloc[0]["symbol"], splash_sid])

    # 2. normalize="clean_tpm": exact parity between genes=None and genes=["RPLP0", "ACTB"]!
    clean_all = gtex.sample_tpm("thyroid", qc="all", cleanup="all", normalize="clean_tpm")
    clean_sub = gtex.sample_tpm(
        "thyroid", genes=["RPLP0", "ACTB"], qc="all", cleanup="all", normalize="clean_tpm"
    )
    pd.testing.assert_frame_equal(
        clean_all[clean_all["Description"].isin(["ACTB", "RPLP0"])].reset_index(drop=True),
        clean_sub.reset_index(drop=True),
    )
    assert np.isnan(clean_all.set_index("Description").at[testis_a.iloc[0]["symbol"], splash_sid])

    # 3. tissue_summary preserves cleanup=frozenset() when a GtexSpec is passed via qc=!
    custom_spec = gtex.GtexSpec(qc=["swaps"], cleanup=None, sex="XX")
    sum_df = gtex.tissue_summary("thyroid", qc=custom_spec)
    assert sum_df.attrs["spec"]["cleanup"] == []
    assert sum_df.attrs["spec"]["sex"] == "XX"

    # 4. Empty tissues_select=[] returns empty metadata frame cleanly; unknown tissue raises ValueError
    empty_sum = gtex.tissue_summary([], version="v11")
    assert list(empty_sum.columns) == ["Name", "Description"]
    with pytest.raises(ValueError, match="Unknown GTEx"):
        gtex.tissue_summary("not_a_real_organ", version="v10", stat="median")
    with pytest.raises(ValueError, match="Donor sex filtering"):
        gtex.tissue_rna("v11", qc=gtex.GtexSpec(sex="XX"))

    # 5. Derived Parquet cache write, reuse, and clear
    gtex.clear_derived_cache()
    m1 = gtex.tissue_rna("v11", qc=["cell_lines"], cleanup=None, include_subcompartments=False)
    assert list(m1.columns) == ["Name", "Description", "Testis", "Thyroid"]
    assert len(gtex.cached_specs()) == 1
    assert gtex.clear_derived_cache() >= 2
    assert gtex.cached_specs().empty


def test_real_gtex_v11_sample_tpm_qc_only_vs_cleanup_and_derived_cache(monkeypatch, tmp_path):
    required = [
        (gtex.SOURCE_ID, "v11"),
        (gtex.SHARD_SOURCE_ID, "v11-ovary"),
        (gtex.SHARD_SOURCE_ID, "v11-breast_mammary_tissue"),
        (gtex.SHARD_SOURCE_ID, "v11-testis"),
    ]
    if not all(reference_data.is_cached(src, ver) for src, ver in required):
        pytest.skip("Real GTEx v11 shards not cached locally")

    deriv_dir = tmp_path / "derived_gtex_real"
    monkeypatch.setattr(gtex, "derived_cache_dir", lambda: deriv_dir)

    # Verify testis qc="all" preserves 100% of endogenous ACTB, GAPDH, SOX9, & PRM1 samples!
    testis_df = gtex.sample_tpm(
        "testis", genes=["ACTB", "GAPDH", "SOX9", "PRM1"], qc="all", cleanup=None
    )
    t_samps = [c for c in testis_df.columns if c not in ("Name", "Description")]
    assert len(t_samps) == 414
    assert int(testis_df[t_samps].notna().all(axis=1).sum()) == 4

    raw_ov = gtex.sample_tpm("ovary", genes=["SFTPB", "FOXL2"], qc=None, cleanup=None)
    assert "GTEX-1GZHY-2426-SM-9JGGK" in raw_ov.columns
    qc_ov = gtex.sample_tpm("ovary", genes=["SFTPB", "FOXL2"], qc=["swaps"], cleanup=None)
    assert "GTEX-1GZHY-2426-SM-9JGGK" not in qc_ov.columns
    surviving = [c for c in qc_ov.columns if c not in ("Name", "Description")]
    pd.testing.assert_series_equal(
        qc_ov.set_index("Description").loc["FOXL2", surviving],
        raw_ov.set_index("Description").loc["FOXL2", surviving],
    )

    xx_br = gtex.tissue_summary(
        "breast_mammary_tissue",
        genes=["PRM1", "TGM4", "BTN1A1"],
        stat="host_exposure",
        qc="all",
        cleanup="all",
        sex="XX",
    ).set_index("Description")
    raw_br = gtex.tissue_summary(
        "breast_mammary_tissue",
        genes=["PRM1", "TGM4", "BTN1A1"],
        stat="max",
        qc=None,
        cleanup=None,
    ).set_index("Description")
    assert raw_br.at["PRM1", "breast_mammary_tissue"] > 50.0
    assert xx_br.at["PRM1", "breast_mammary_tissue"] < 2.0
    assert xx_br.at["TGM4", "breast_mammary_tissue"] < 1.0


def test_comprehensive_normal_tissue_source_grid():
    grid = normal_tissues.normal_tissue_source_grid()
    assert len(grid) >= 40
    by_cid = grid.set_index("concept_id")
    assert by_cid.at["retina", "gtex_v11_shards"] == ""
    assert "retina" in by_cid.at["retina", "gtex_via_hpa_v8"]
    assert "EyeGEx" in by_cid.at["retina", "notes"]
    # Heart, substantia_nigra, brain, and endometrium all resolve GTEx v11 shards without bogus EyeGEx notes!
    for present_cid in ("heart", "substantia_nigra", "brain", "endometrium"):
        assert by_cid.at[present_cid, "gtex_v11_shards"] != "", present_cid
        assert "EyeGEx" not in by_cid.at[present_cid, "notes"], present_cid
    assert by_cid.at["placenta", "gtex_v11_shards"] == ""
    assert "placenta" in by_cid.at["placenta", "hpa_surgical_biopsies"]
    assert "Pancreas_Islets" in by_cid.at["pancreas", "gtex_v11_subcompartments"]
    assert "olfactory_epithelium" in by_cid.index
    assert by_cid.at["olfactory_epithelium", "specialized_and_developmental"] != ""
    missing = human.gtex.coverage_grid(missing_from_gtex=True)
    assert "heart" not in set(missing["concept_id"])
    assert len(missing) >= 10


def test_real_cached_hpa_gtex_v23_cleans_prm1_splash_and_preserves_retina(monkeypatch, tmp_path):
    if not (
        reference_data.is_cached("rna_tissue_gtex", "v23")
        and reference_data.is_cached("rna_brain_gtex", "v23")
    ):
        pytest.skip("HPA GTEx v23 tables not cached locally")

    monkeypatch.setattr(gtex, "derived_cache_dir", lambda: tmp_path / "derived_hpa")
    raw_hpa = gtex.tissue_rna(
        provider="hpa", version="v23", qc=None, cleanup=None, value_column="nTPM"
    )
    clean_hpa = gtex.tissue_rna(
        provider="hpa", version="v23", qc="all", cleanup="all", value_column="nTPM"
    )
    prm_raw = raw_hpa[raw_hpa["Description"].eq("PRM1")].iloc[0]
    prm_clean = clean_hpa[clean_hpa["Description"].eq("PRM1")].iloc[0]

    assert float(prm_raw["breast"]) > 8.0
    assert float(prm_raw["skeletal muscle"]) > 4.0
    assert float(prm_clean["breast"]) < 2.5
    assert float(prm_clean["skeletal muscle"]) < 2.5
    assert float(prm_clean["testis"]) > 10000.0
    assert "retina" in clean_hpa.columns


def test_gtex_subpackage_shared_bulk_qc_and_publication_figure_export(tmp_path: Path):
    # Subpackage & shared bulk_qc wiring
    assert gtex.deconvolution.gtex_quarantined_swaps is gtex_deconvolution.gtex_quarantined_swaps
    assert (
        gtex.deconvolution.solve_nonnegative_archetype_mixture
        is bulk_qc.solve_nonnegative_archetype_mixture
    )
    assert (
        gtex.deconvolution.purity_weighted_harrell_davis_quantile
        is bulk_qc.purity_weighted_harrell_davis_quantile
    )
    assert bulk_qc.expand_ops(
        "all", allowed=("a", "b"), aliases={"all": frozenset({"a", "b"})}, param_name="qc"
    ) == frozenset({"a", "b"})

    # Publication figure export writes both 300-dpi PNG and vector SVG
    png, svg = gtex.plots.plot_normal_source_coverage_grid(tmp_path / "cov_grid")
    assert png.is_file() and png.stat().st_size > 10_000
    assert svg.is_file() and svg.stat().st_size > 5_000

    # Defensive edge cases: Mapping rejected in expand_ops, NaN preserved in betainc,
    # read-only cache fallback, and stat / value_column validation across providers
    with pytest.raises(TypeError, match="dict/mapping"):
        gtex.resolve_spec(qc={"swaps": False})  # type: ignore[arg-type]
    beta_out = bulk_qc.regularized_betainc(2.0, 2.0, np.array([0.0, 0.5, np.nan, 1.0]))
    assert np.isnan(beta_out[2]) and beta_out[1] == pytest.approx(0.5)
    with pytest.raises(ValueError, match="arithmetic means only"):
        gtex.tissue_summary(provider="hpa", stat="p95")
    with pytest.raises(ValueError, match="'TPM' only"):
        gtex.tissue_rna("v11", provider="direct", value_column="nTPM")

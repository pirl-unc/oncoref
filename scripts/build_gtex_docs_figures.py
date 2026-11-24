#!/usr/bin/env python3
"""Generate publication-ready GTEx QC/deconvolution figures & evaluation tables (#778)."""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import gtex, reference_data  # noqa: E402

AUDIT_DIR = ROOT / "docs/audits/gtex-holistic-deconvolution-20261124"


def _extract_stoichiometry_and_hopping() -> tuple[list[dict], pd.DataFrame, pd.DataFrame]:
    """Extract real multi-gene stoichiometry traces, hopping sample points, and shard sex counts."""
    anchors = gtex.deconvolution.gtex_deconvolution_anchors()
    dsex_map = gtex.donor_sex().set_index("donor_id")["chromosomal_sex"].to_dict()

    # Showcase 4 real libraries across 4 orders of magnitude
    showcase_specs = [
        (
            "ovary",
            "lung",
            "GTEX-1GZHY-2426-SM-9JGGK",
            "XY-Lung swap in XX-Ovary (133%)",
            1.335,
        ),
        (
            "skin_sun_exposed_lower_leg",
            "pituitary",
            "GTEX-145LU-0126-SM-5TDCB",
            "Pituitary splash in Skin (23.8%, SM-5TDC)",
            0.2378,
        ),
        (
            "esophagus_muscularis",
            "testis",
            "GTEX-1I1GR-1426-SM-CY8H8",
            "Testis splash in Esophagus (10.7%, BP-72061)",
            0.1072,
        ),
        (
            "thyroid",
            "testis",
            None,  # pick a clean ~0.03% index-hopped thyroid sample
            "Testis index hopping in Thyroid (0.03%)",
            0.0003,
        ),
    ]

    traces: list[dict] = []
    hop_rows: list[dict] = []

    for tgt_site, src_site, target_sid, label, fallback_alpha in showcase_specs:
        sub_a = anchors[anchors["source_site"].eq(src_site)].copy()
        excl_ok = [
            tgt_site not in {x for x in str(ex).split(";") if x}
            for ex in sub_a["excluded_target_sites"]
        ]
        sub_a = sub_a.loc[excl_ok]
        df_s = gtex.sample_tpm(tgt_site, genes=sub_a["gene_id"].tolist(), qc=None, cleanup=None)
        df_s["_gid"] = df_s["Name"].astype(str).str.split(".").str[0]
        df_s = df_s.drop_duplicates("_gid").set_index("_gid")

        if target_sid is None:
            s_cols = [c for c in df_s.columns if c not in ("Name", "Description")]
            meds = sub_a.set_index("gene_id")["source_median_tpm"]
            common = [g for g in meds.index if g in df_s.index]
            mat = df_s.loc[common, s_cols].to_numpy(dtype=float)
            m_vec = meds.loc[common].to_numpy(dtype=float)[:, None]
            det = mat >= 0.20
            for j, sc_id in enumerate(s_cols):
                if int(det[:, j].sum()) >= 4:
                    rv = mat[det[:, j], j] / m_vec[det[:, j], 0]
                    ah = float(np.median(rv))
                    lr = np.log10(np.maximum(rv, 1e-12))
                    iqr = float(np.percentile(lr, 75.0) - np.percentile(lr, 25.0))
                    if 2e-5 <= ah < 5e-3 and iqr <= 0.55:
                        hop_rows.append({"alpha_hat": ah, "log10_iqr": iqr})
                        if (
                            target_sid is None
                            and 2.2e-4 <= ah <= 4.5e-4
                            and int(det[:, j].sum()) >= 5
                        ):
                            target_sid = sc_id
                            fallback_alpha = ah

        if target_sid and target_sid in df_s.columns:
            sx, sy = [], []
            for r in sub_a.itertuples(index=False):
                if r.gene_id in df_s.index:
                    val = float(df_s.at[r.gene_id, target_sid])
                    if val >= 0.08:
                        sx.append(float(r.source_median_tpm))
                        sy.append(val)
            if sx:
                traces.append(
                    {
                        "label": label,
                        "source_tpm": sx,
                        "recipient_tpm": sy,
                        "alpha_hat": fallback_alpha,
                    }
                )

    # Also sample index-hopping points from liver, adrenal_gland, colon_transverse, ovary, breast
    for extra_site in (
        "liver",
        "adrenal_gland",
        "colon_transverse",
        "ovary",
        "breast_mammary_tissue",
    ):
        for src_s in ("testis", "pancreas", "muscle_skeletal", "lung"):
            if gtex.deconvolution.is_anatomical_boundary(src_s, extra_site):
                continue
            sub_a = anchors[anchors["source_site"].eq(src_s)]
            excl_ok = [
                extra_site not in {x for x in str(ex).split(";") if x}
                for ex in sub_a["excluded_target_sites"]
            ]
            sub_a = sub_a.loc[excl_ok]
            if len(sub_a) < 3:
                continue
            df_s = gtex.sample_tpm(
                extra_site, genes=sub_a["gene_id"].tolist(), qc=None, cleanup=None
            )
            df_s["_gid"] = df_s["Name"].astype(str).str.split(".").str[0]
            df_s = df_s.drop_duplicates("_gid").set_index("_gid")
            s_cols = [c for c in df_s.columns if c not in ("Name", "Description")]
            meds = sub_a.set_index("gene_id")["source_median_tpm"]
            common = [g for g in meds.index if g in df_s.index]
            mat = df_s.loc[common, s_cols].to_numpy(dtype=float)
            m_vec = meds.loc[common].to_numpy(dtype=float)[:, None]
            det = mat >= 0.25
            for j in np.where(det.sum(axis=0) >= 3)[0]:
                rv = mat[det[:, j], j] / m_vec[det[:, j], 0]
                ah = float(np.median(rv))
                if 2e-5 <= ah < 5e-3:
                    lr = np.log10(np.maximum(rv, 1e-12))
                    iqr = float(np.percentile(lr, 75.0) - np.percentile(lr, 25.0))
                    if iqr <= 0.60:
                        hop_rows.append({"alpha_hat": ah, "log10_iqr": iqr})

    # Compute exact XX vs XY biopsy counts across 12 representative shards from header lines only
    rep_sites = [
        "ovary",
        "uterus",
        "vagina",
        "breast_mammary_tissue",
        "thyroid",
        "lung",
        "liver",
        "heart_left_ventricle",
        "muscle_skeletal",
        "whole_blood",
        "prostate",
        "testis",
    ]
    sex_rows = []
    swaps_set = set(gtex.deconvolution.gtex_quarantined_swaps()["sample_id"])
    for st in rep_sites:
        p = reference_data.pinned_path("gtex_individual_tissue_tpm", f"v11-{st}")
        with gzip.open(p, "rb") as fh:
            fh.readline()
            fh.readline()
            samps = [
                s
                for s in fh.readline().decode("utf-8").rstrip("\n").split("\t")[2:]
                if s not in swaps_set
            ]
        xx_n = sum(
            1
            for s in samps
            if dsex_map.get(gtex.deconvolution.parse_gtex_sample_barcode(s).donor_id) == "XX"
        )
        xy_n = sum(
            1
            for s in samps
            if dsex_map.get(gtex.deconvolution.parse_gtex_sample_barcode(s).donor_id) == "XY"
        )
        sex_rows.append({"site": st, "xx": xx_n, "xy": xy_n})

    return traces, pd.DataFrame(hop_rows), pd.DataFrame(sex_rows)


def _compute_evaluation_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute real before/after evaluation tables on GTEx v11 and HPA v23."""
    # 1. Panel A: Direct GTEx v11 contaminated gene x shard benchmarks
    bench = [
        ("ovary", "SFTPB", "max", None, "SFTPB in Ovary (max, XY-Lung swap)"),
        ("breast_mammary_tissue", "PRM1", "host_exposure", "XX", "PRM1 in XX Breast (host P99)"),
        ("breast_mammary_tissue", "TGM4", "host_exposure", "XX", "TGM4 in XX Breast (host P99)"),
        ("ovary", "PRM1", "host_exposure", "XX", "PRM1 in Ovary (host P99)"),
        ("thyroid", "PRM1", "p95", None, "PRM1 in Thyroid (somatic P95)"),
        ("thyroid", "PRM2", "p95", None, "PRM2 in Thyroid (somatic P95)"),
        ("esophagus_muscularis", "PRM1", "p99", None, "PRM1 in Esoph. Musc. (P99)"),
    ]
    eval_rows = []
    for site, sym, stat_m, sex_f, lbl in bench:
        r_raw = gtex.tissue_summary(
            site,
            genes=[sym],
            stat="max" if stat_m == "host_exposure" else stat_m,
            qc=None,
            cleanup=None,
            use_cache=False,
        )
        r_qc = gtex.tissue_summary(
            site, genes=[sym], stat=stat_m, qc="all", cleanup=None, sex=sex_f, use_cache=False
        )
        r_full = gtex.tissue_summary(
            site, genes=[sym], stat=stat_m, qc="all", cleanup="all", sex=sex_f, use_cache=False
        )
        eval_rows.append(
            {
                "site": site,
                "symbol": sym,
                "statistic": stat_m,
                "sex": sex_f or "all",
                "label": lbl,
                "raw_tpm": round(float(r_raw.set_index("Description").at[sym, site]), 4),
                "qc_only_tpm": round(float(r_qc.set_index("Description").at[sym, site]), 4),
                "full_clean_tpm": round(float(r_full.set_index("Description").at[sym, site]), 4),
            }
        )
    eval_df = pd.DataFrame(eval_rows)

    # 2. Panel B: HPA v23 arithmetic-mean GTEx v8 PRM1 before vs after cleanup
    raw_hpa = gtex.tissue_rna(
        provider="hpa", version="v23", qc=None, cleanup=None, value_column="nTPM", use_cache=False
    )
    clean_hpa = gtex.tissue_rna(
        provider="hpa", version="v23", qc="all", cleanup="all", value_column="nTPM", use_cache=False
    )
    prm_raw = raw_hpa[raw_hpa["Description"].eq("PRM1")].iloc[0]
    prm_clean = clean_hpa[clean_hpa["Description"].eq("PRM1")].iloc[0]
    hpa_tissues = [
        "breast",
        "skeletal muscle",
        "thyroid gland",
        "colon",
        "adipose tissue",
        "retina",
    ]
    hpa_df = pd.DataFrame(
        [
            {
                "tissue": t,
                "raw_ntpm": round(float(prm_raw[t]), 4),
                "clean_ntpm": round(float(prm_clean[t]), 4),
            }
            for t in hpa_tissues
            if t in raw_hpa.columns
        ]
    )

    # 3. Panel C: Endogenous housekeeping, native lineage TFs, and surgical boundary negative controls
    ctrl_specs = [
        ("testis", "ACTB", "p95", "ACTB in Testis (housekeeping)"),
        ("testis", "GAPDH", "p95", "GAPDH in Testis (housekeeping)"),
        ("testis", "SOX9", "p95", "SOX9 in Testis (Sertoli TF)"),
        ("ovary", "FOXL2", "p95", "FOXL2 in Ovary (granulosa TF)"),
        ("bladder", "KLK3", "p95", "KLK3 in Bladder (urogenital boundary)"),
        ("bladder", "MSMB", "p95", "MSMB in Bladder (urogenital boundary)"),
    ]
    ctrl_rows = []
    for site, sym, st_m, lbl in ctrl_specs:
        raw_s = gtex.sample_tpm(site, genes=[sym], qc=None, cleanup=None)
        cln_s = gtex.sample_tpm(site, genes=[sym], qc="all", cleanup="all")
        raw_cols = [c for c in raw_s.columns if c not in ("Name", "Description")]
        cln_cols = [c for c in cln_s.columns if c not in ("Name", "Description")]
        tot_d = len({gtex.deconvolution.parse_gtex_sample_barcode(s).donor_id for s in raw_cols})
        fin_samps = [
            s for s in cln_cols if np.isfinite(float(cln_s.set_index("Description").at[sym, s]))
        ]
        kept_d = len({gtex.deconvolution.parse_gtex_sample_barcode(s).donor_id for s in fin_samps})

        v_raw = float(
            gtex.tissue_summary(
                site, genes=[sym], stat=st_m, qc=None, cleanup=None, use_cache=False
            )
            .set_index("Description")
            .at[sym, site]
        )
        v_cln = float(
            gtex.tissue_summary(
                site, genes=[sym], stat=st_m, qc="all", cleanup="all", use_cache=False
            )
            .set_index("Description")
            .at[sym, site]
        )
        ctrl_rows.append(
            {
                "site": site,
                "symbol": sym,
                "label": lbl,
                "raw_val": round(v_raw, 3),
                "clean_val": round(v_cln, 3),
                "kept_donors": kept_d,
                "total_donors": tot_d,
            }
        )
    controls_df = pd.DataFrame(ctrl_rows)

    # 4. Panel D: Cross-tissue housekeeping CV across all 52 GTEx v11 primary shards before vs after clean_tpm
    raw_m = gtex.tissue_rna(
        "v11",
        qc=["cell_lines"],
        cleanup=None,
        normalize=None,
        include_subcompartments=False,
        use_cache=False,
    )
    cln_m = gtex.tissue_rna(
        "v11",
        qc=["cell_lines"],
        cleanup="all",
        normalize="clean_tpm",
        include_subcompartments=False,
        use_cache=False,
    )
    shard_cols = [c for c in raw_m.columns if c not in ("Name", "Description")]
    hk_genes = ["ACTB", "B2M", "GAPDH", "PPIA", "EEF1A1", "PGK1", "YWHAZ", "TBP"]
    hk_rows = []
    raw_by_sym = raw_m.drop_duplicates("Description").set_index("Description")
    cln_by_sym = cln_m.drop_duplicates("Description").set_index("Description")
    for sym in hk_genes:
        if sym not in raw_by_sym.index:
            continue
        rv = raw_by_sym.loc[sym, shard_cols].to_numpy(dtype=float)
        cv_raw = float(np.std(rv) / np.maximum(np.mean(rv), 1e-9))
        cv_v = cln_by_sym.loc[sym, shard_cols].to_numpy(dtype=float)
        cv_cln = float(np.std(cv_v) / np.maximum(np.mean(cv_v), 1e-9))
        hk_rows.append(
            {
                "symbol": sym,
                "raw_cv": round(cv_raw, 4),
                "clean_tpm_cv": round(cv_cln, 4),
                "cv_reduction_pct": round(100.0 * (cv_raw - cv_cln) / cv_raw, 2),
            }
        )
    hk_cv_df = pd.DataFrame(hk_rows)
    return eval_df, hpa_df, controls_df, hk_cv_df


def main() -> None:
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    traces, hop_df, sex_df = _extract_stoichiometry_and_hopping()
    if "--reuse-tables" in sys.argv and (AUDIT_DIR / "evaluation-before-after.csv").is_file():
        eval_df = pd.read_csv(AUDIT_DIR / "evaluation-before-after.csv")
        hpa_df = pd.read_csv(AUDIT_DIR / "hpa-gtex-v8-prm1-cleanup.csv")
        controls_df = pd.read_csv(AUDIT_DIR / "negative-controls-retention.csv")
        hk_cv_df = pd.read_csv(AUDIT_DIR / "clean-tpm-housekeeping-cv.csv")
    else:
        eval_df, hpa_df, controls_df, hk_cv_df = _compute_evaluation_tables()
        eval_df.to_csv(AUDIT_DIR / "evaluation-before-after.csv", index=False, lineterminator="\n")
        hpa_df.to_csv(AUDIT_DIR / "hpa-gtex-v8-prm1-cleanup.csv", index=False, lineterminator="\n")
        controls_df.to_csv(
            AUDIT_DIR / "negative-controls-retention.csv", index=False, lineterminator="\n"
        )
        hk_cv_df.to_csv(
            AUDIT_DIR / "clean-tpm-housekeeping-cv.csv", index=False, lineterminator="\n"
        )

    f1 = gtex.plots.plot_contamination_regimes(
        AUDIT_DIR / "gtex-contamination-regimes", hopping_sample=hop_df
    )
    f2 = gtex.plots.plot_splash_matrix_and_stoichiometry(
        AUDIT_DIR / "gtex-splash-matrix-and-stoichiometry", stoichiometry_traces=traces
    )
    f3 = gtex.plots.plot_donor_sex_and_mammary_bias(
        AUDIT_DIR / "gtex-donor-sex-and-mammary-bias", shard_sex_counts=sex_df
    )
    f4 = gtex.plots.plot_cleanup_evaluation(
        AUDIT_DIR / "gtex-cleanup-evaluation",
        eval_df=eval_df,
        hpa_df=hpa_df,
        controls_df=controls_df,
        hk_cv_df=hk_cv_df,
    )
    f5 = gtex.plots.plot_normal_source_coverage_grid(AUDIT_DIR / "gtex-normal-source-coverage-grid")

    print(
        json.dumps(
            {
                "figures": [p.name for pair in (f1, f2, f3, f4, f5) for p in pair],
                "eval_rows": len(eval_df),
                "hpa_rows": len(hpa_df),
                "controls_rows": len(controls_df),
                "hk_cv_rows": len(hk_cv_df),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

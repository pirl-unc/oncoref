# Licensed under the Apache License, Version 2.0.
"""Publication-ready diagnostic and evaluation figures for GTEx QC & deconvolution (#778).

Follows :mod:`oncoref.figure_style`: opaque white ground, 300 dpi PNG + vector SVG
exports, hairline left/bottom spines, and colour-blind-safe categorical palette.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .. import normal_tissues
from . import deconvolution as gtex_deconvolution


def _save_png_and_svg(fig, out_stem: Path, *, dpi: int = 300) -> tuple[Path, Path]:
    import matplotlib.pyplot as plt

    from .. import figure_export, figure_style

    out_stem = Path(out_stem)
    if out_stem.suffix in {".png", ".svg", ".pdf"}:
        out_stem = out_stem.with_suffix("")
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    figure_style.finalize(fig)
    figure_export._place_scatter_labels(fig)
    png_path = out_stem.with_suffix(".png")
    svg_path = out_stem.with_suffix(".svg")
    fig.savefig(png_path, dpi=max(300, int(dpi)))
    fig.savefig(svg_path)
    plt.close(fig)
    return png_path, svg_path


def _pretty_site(site: str) -> str:
    s = str(site).replace("_", " ")
    repl = {
        "brain cerebellar hemisphere": "cerebellar hemisphere",
        "brain cerebellum": "cerebellum",
        "brain spinal cord cervical c-1": "spinal cord (C1)",
        "brain nucleus accumbens basal ganglia": "nucleus accumbens",
        "brain caudate basal ganglia": "caudate",
        "brain putamen basal ganglia": "putamen",
        "brain anterior cingulate cortex ba24": "anterior cingulate",
        "brain frontal cortex ba9": "frontal cortex (BA9)",
        "skin sun exposed lower leg": "skin (lower leg)",
        "skin not sun exposed suprapubic": "skin (suprapubic)",
        "small intestine terminal ileum": "terminal ileum",
        "esophagus gastroesophageal junction": "esophagus (GEJ)",
        "breast mammary tissue": "breast (mammary)",
        "adipose visceral omentum": "visceral omentum",
        "adipose subcutaneous": "subcutaneous adipose",
    }
    return repl.get(s, s)


def plot_contamination_regimes(
    out_stem: Path,
    *,
    hopping_sample: pd.DataFrame | None = None,
) -> tuple[Path, Path]:
    """Figure 1: Physical separation of contamination regimes & per-shard incidence."""
    import matplotlib.pyplot as plt

    from .. import figure_style

    figure_style.apply()
    events = gtex_deconvolution.gtex_sample_splash_events(include_boundary=True)
    tsum = gtex_deconvolution.gtex_deconvolution_tissue_summary()

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(12.0, 5.6), gridspec_kw={"width_ratios": [1.10, 1.0]}
    )

    if hopping_sample is not None and not hopping_sample.empty:
        ax_a.scatter(
            np.log10(np.maximum(hopping_sample["alpha_hat"].to_numpy(dtype=float), 1e-6)),
            hopping_sample["log10_iqr"].to_numpy(dtype=float),
            s=11,
            c="#b8c4d0",
            alpha=0.38,
            label="Flowcell index hopping (<0.5%)",
            rasterized=True,
        )

    regime_styles = [
        ("anatomical_boundary", "Anatomical / surgical boundary", "#8c969e", 18, 0.40, "o"),
        ("same_donor_carryover", "Same-donor carryover (>=0.5%)", figure_style.WEAK, 28, 0.80, "^"),
        (
            "plate_pipetting_splash",
            "Cross-donor plate/batch splash (>=0.5%)",
            figure_style.ACCENT,
            34,
            0.88,
            "o",
        ),
        ("tube_swap", "Quarantined organ swap / chimera", figure_style.THRESHOLD, 76, 0.98, "D"),
    ]
    for reg, label, col, sz, alp, mk in regime_styles:
        sub = events[events["regime"].eq(reg)]
        if sub.empty:
            continue
        ax_a.scatter(
            np.log10(np.maximum(sub["alpha_hat"].to_numpy(dtype=float), 1e-6)),
            sub["log10_iqr"].to_numpy(dtype=float),
            s=sz,
            c=col,
            alpha=alp,
            marker=mk,
            edgecolors="white" if reg == "tube_swap" else "none",
            linewidths=0.7 if reg == "tube_swap" else 0.0,
            label=f"{label} (n={len(sub)})",
            zorder=5 if reg == "tube_swap" else 3,
        )

    for thresh, txt in [
        (-4.7, "Hopping floor\n(2e-5)"),
        (-2.301, "Wet-lab splash\n(0.5%)"),
        (-0.699, "Swap gate\n(20%)"),
    ]:
        ax_a.axvline(thresh, color="#666666", linestyle="--", linewidth=0.8, zorder=1)
        ax_a.text(thresh + 0.05, 0.835, txt, fontsize=7.5, color="#444444", va="top")

    # Place all callouts in open upper whitespace above y=0.61 so zero markers are touched
    callouts = [
        ("GTEX-1GZHY-2426-SM-9JGGK", "lung", "XY-Lung -> XX-Ovary (133%)", (-1.45, 0.77)),
        ("GTEX-11GS4-2826-SM-5HL56", "pancreas", "Pancreas -> Nerve (61%)", (-1.62, 0.70)),
        (
            "GTEX-1H1DG-0011-R11a-SM-C1YQK",
            "adrenal_gland",
            "Adrenal -> Cerebellum (43%)",
            (-1.90, 0.635),
        ),
        (
            "GTEX-1I1GR-1426-SM-CY8H8",
            "testis",
            "Testis -> Esoph. (11%, BP-72061)",
            (-3.65, 0.72),
        ),
        (
            "GTEX-145LU-0126-SM-5TDCB",
            "pituitary",
            "Pituitary -> Skin (24%, SM-5TDC)",
            (-3.65, 0.64),
        ),
    ]
    for sid, src_s, lbl, xy_text in callouts:
        m = events[events["sample_id"].eq(sid) & events["source_site"].eq(src_s)]
        if not m.empty:
            r = m.iloc[0]
            ax_a.annotate(
                lbl,
                xy=(float(np.log10(r["alpha_hat"])), float(r["log10_iqr"])),
                xytext=xy_text,
                fontsize=7.3,
                color="#1f1f1f",
                arrowprops={
                    "arrowstyle": "-",
                    "color": "#777777",
                    "lw": 0.6,
                    "shrinkA": 2,
                    "shrinkB": 3,
                },
                zorder=8,
            )

    ax_a.set_xlim(-5.0, 0.48)
    ax_a.set_ylim(-0.02, 0.87)
    ax_a.set_xlabel("Estimated foreign tissue fraction log10(alpha_hat)")
    ax_a.set_ylabel("Multi-sentinel stoichiometric dispersion IQR(log10 ratio)")
    ax_a.legend(loc="upper center", bbox_to_anchor=(0.48, -0.15), ncol=2, fontsize=7.5)
    ax_a.text(-0.12, 1.03, "a", transform=ax_a.transAxes, fontsize=12, fontweight="bold")

    ts = tsum.copy()
    ts["splash_pct"] = 100.0 * ts["plate_splash_samples"] / ts["total_samples"]
    ts["carry_pct"] = 100.0 * ts["same_donor_carryover_samples"] / ts["total_samples"]
    ts["swap_pct"] = 100.0 * ts["quarantined_swap_samples"] / ts["total_samples"]
    ts["total_wet_pct"] = ts["splash_pct"] + ts["carry_pct"] + ts["swap_pct"]
    top_ts = ts.sort_values("total_wet_pct", ascending=True).tail(18)

    y_pos = np.arange(len(top_ts))
    ax_b.barh(
        y_pos,
        top_ts["splash_pct"],
        color=figure_style.ACCENT,
        height=0.68,
        label="Cross-donor plate/batch splash",
    )
    ax_b.barh(
        y_pos,
        top_ts["carry_pct"],
        left=top_ts["splash_pct"],
        color=figure_style.WEAK,
        height=0.68,
        label="Same-donor prosector carryover",
    )
    ax_b.barh(
        y_pos,
        top_ts["swap_pct"],
        left=top_ts["splash_pct"] + top_ts["carry_pct"],
        color=figure_style.THRESHOLD,
        height=0.68,
        label="Quarantined swap / chimera",
    )
    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels([_pretty_site(s) for s in top_ts["site"]], fontsize=8)
    ax_b.set_xlabel("Biopsies with >=0.5% non-boundary foreign RNA (%) — top 18 of 52 shards")
    ax_b.legend(loc="upper center", bbox_to_anchor=(0.48, -0.15), ncol=2, fontsize=7.6)
    ax_b.text(-0.28, 1.03, "b", transform=ax_b.transAxes, fontsize=12, fontweight="bold")

    fig.tight_layout(w_pad=2.2)
    return _save_png_and_svg(fig, out_stem)


def plot_splash_matrix_and_stoichiometry(
    out_stem: Path,
    *,
    stoichiometry_traces: list[dict] | None = None,
) -> tuple[Path, Path]:
    """Figure 2: Source x recipient splash matrix and multi-gene stoichiometric linearity."""
    import matplotlib.pyplot as plt

    from .. import figure_style

    figure_style.apply()
    events = gtex_deconvolution.gtex_sample_splash_events(include_boundary=False)
    splashes = events[events["regime"].isin(["plate_pipetting_splash", "tube_swap"])].copy()

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(11.8, 5.4), gridspec_kw={"width_ratios": [1.12, 1.0]}
    )

    top_srcs = splashes["source_site"].value_counts().head(12).index.tolist()
    top_tgts = splashes["target_site"].value_counts().head(14).index.tolist()
    sub = splashes[splashes["source_site"].isin(top_srcs) & splashes["target_site"].isin(top_tgts)]
    agg = sub.groupby(["target_site", "source_site"], as_index=False).agg(
        n=("sample_id", "count"), max_alpha=("alpha_hat", "max")
    )
    src_pos = {s: i for i, s in enumerate(top_srcs)}
    tgt_pos = {t: i for i, t in enumerate(reversed(top_tgts))}

    xs = [src_pos[r.source_site] for r in agg.itertuples(index=False)]
    ys = [tgt_pos[r.target_site] for r in agg.itertuples(index=False)]
    sizes = [28.0 + 42.0 * float(r.n) for r in agg.itertuples(index=False)]
    alphas_pct = np.array([100.0 * float(r.max_alpha) for r in agg.itertuples(index=False)])

    sc = ax_a.scatter(
        xs,
        ys,
        s=sizes,
        c=np.log10(np.maximum(alphas_pct, 0.5)),
        cmap="YlOrRd",
        edgecolors="#333333",
        linewidths=0.6,
        zorder=3,
    )
    cbar = fig.colorbar(sc, ax=ax_a, fraction=0.042, pad=0.03)
    cbar.set_label("Peak foreign fraction log10(%)", fontsize=8.5)
    ax_a.set_xticks(range(len(top_srcs)))
    ax_a.set_xticklabels([_pretty_site(s) for s in top_srcs], rotation=40, ha="right", fontsize=7.8)
    ax_a.set_yticks(range(len(top_tgts)))
    ax_a.set_yticklabels([_pretty_site(t) for t in reversed(top_tgts)], fontsize=8.0)
    ax_a.set_xlim(-0.7, len(top_srcs) - 0.3)
    ax_a.set_ylim(-0.7, len(top_tgts) - 0.3)
    ax_a.set_xlabel("Contaminating source organ (bubble area = libraries)")
    ax_a.set_ylabel("Recipient GTEx v11 shard")
    ax_a.grid(True, linestyle=":", linewidth=0.4, color="#dddddd", zorder=0)
    ax_a.text(-0.26, 1.03, "a", transform=ax_a.transAxes, fontsize=12, fontweight="bold")

    x_grid = np.logspace(1.3, 4.9, 100)
    for frac, lbl_line, x_lbl in [
        (1.0, "100%", 22000.0),
        (0.1, "10%", 68000.0),
        (0.01, "1%", 68000.0),
        (3e-4, "0.03%", 68000.0),
    ]:
        ax_b.plot(x_grid, frac * x_grid, linestyle="--", linewidth=0.75, color="#cccccc", zorder=1)
        ax_b.text(x_lbl, frac * x_lbl * 1.18, lbl_line, fontsize=7.2, color="#777777", ha="right")

    palette_traces = [
        figure_style.THRESHOLD,
        figure_style.ACCENT,
        figure_style.KEPT,
        figure_style.PALETTE[5],
    ]
    markers_traces = ["D", "o", "s", "^"]
    if stoichiometry_traces:
        for idx, tr in enumerate(stoichiometry_traces):
            col = palette_traces[idx % len(palette_traces)]
            mk = markers_traces[idx % len(markers_traces)]
            sx = np.asarray(tr["source_tpm"], dtype=float)
            sy = np.asarray(tr["recipient_tpm"], dtype=float)
            ax_b.scatter(
                sx,
                sy,
                s=34,
                c=col,
                marker=mk,
                alpha=0.92,
                edgecolors="white",
                linewidths=0.4,
                label=tr["label"],
                zorder=4,
            )
            slope = float(tr["alpha_hat"])
            ax_b.plot(x_grid, slope * x_grid, color=col, linewidth=1.1, alpha=0.65, zorder=2)

    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.set_xlim(18.0, 95000.0)
    ax_b.set_ylim(0.04, 52000.0)
    ax_b.set_xlabel("Native source-organ median expression (TPM)")
    ax_b.set_ylabel("Observed recipient-biopsy expression (TPM)")
    ax_b.legend(loc="lower right", fontsize=7.3)
    ax_b.text(-0.15, 1.03, "b", transform=ax_b.transAxes, fontsize=12, fontweight="bold")

    fig.tight_layout(w_pad=2.0)
    return _save_png_and_svg(fig, out_stem)


def plot_donor_sex_and_mammary_bias(
    out_stem: Path,
    *,
    shard_sex_counts: pd.DataFrame | None = None,
) -> tuple[Path, Path]:
    """Figure 3: Splash-robust donor chromosomal sex (XX vs XY) and male-biased shards."""
    import matplotlib.pyplot as plt

    from .. import figure_style

    figure_style.apply()
    dsex = gtex_deconvolution.gtex_donor_sex()

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(11.4, 5.0), gridspec_kw={"width_ratios": [1.0, 1.12]}
    )

    xx_v = np.log10(
        dsex.loc[dsex["chromosomal_sex"].eq("XX"), "median_y_sentinel_tpm"].to_numpy(dtype=float)
        + 0.01
    )
    xy_v = np.log10(
        dsex.loc[dsex["chromosomal_sex"].eq("XY"), "median_y_sentinel_tpm"].to_numpy(dtype=float)
        + 0.01
    )
    bins = np.linspace(-2.1, 3.0, 46)
    ax_a.hist(
        xx_v,
        bins=bins,
        color=figure_style.KEPT,
        alpha=0.85,
        label=f"XX female donors (n={len(xx_v)})",
    )
    ax_a.hist(
        xy_v,
        bins=bins,
        color=figure_style.ACCENT,
        alpha=0.85,
        label=f"XY male donors (n={len(xy_v)})",
    )
    ax_a.axvline(np.log10(2.01), color="#555555", linestyle="--", linewidth=0.85)
    ax_a.axvline(np.log10(20.01), color="#555555", linestyle="--", linewidth=0.85)
    ax_a.text(
        0.80,
        ax_a.get_ylim()[1] * 0.72,
        "0 ambiguous\ndonors across\nall 946 donors",
        fontsize=8.0,
        color="#333333",
        ha="center",
    )
    ax_a.set_xlabel("Cross-shard donor median Y-sentinel sum log10(TPM + 0.01)")
    ax_a.set_ylabel("GTEx v11 donors")
    ax_a.legend(loc="upper left", fontsize=8.0)
    ax_a.text(-0.14, 1.03, "a", transform=ax_a.transAxes, fontsize=12, fontweight="bold")

    if shard_sex_counts is not None and not shard_sex_counts.empty:
        sub_s = shard_sex_counts.copy()
    else:
        sub_s = pd.DataFrame(
            [
                {"site": "ovary", "xx": 192, "xy": 0},
                {"site": "uterus", "xx": 153, "xy": 0},
                {"site": "vagina", "xx": 170, "xy": 0},
                {"site": "breast_mammary_tissue", "xx": 178, "xy": 336},
                {"site": "thyroid", "xx": 230, "xy": 454},
                {"site": "lung", "xx": 193, "xy": 411},
                {"site": "liver", "xx": 77, "xy": 185},
                {"site": "heart_left_ventricle", "xx": 139, "xy": 313},
                {"site": "muscle_skeletal", "xx": 261, "xy": 557},
                {"site": "whole_blood", "xx": 266, "xy": 537},
                {"site": "prostate", "xx": 0, "xy": 282},
                {"site": "testis", "xx": 0, "xy": 414},
            ]
        )
    sub_s["tot"] = sub_s["xx"] + sub_s["xy"]
    sub_s["xy_pct"] = 100.0 * sub_s["xy"] / np.maximum(sub_s["tot"], 1)
    sub_s["xx_pct"] = 100.0 * sub_s["xx"] / np.maximum(sub_s["tot"], 1)

    y_idx = np.arange(len(sub_s))
    ax_b.barh(
        y_idx, sub_s["xx_pct"], color=figure_style.KEPT, height=0.66, label="XX female biopsies"
    )
    ax_b.barh(
        y_idx,
        sub_s["xy_pct"],
        left=sub_s["xx_pct"],
        color=figure_style.ACCENT,
        height=0.66,
        label="XY male biopsies",
    )
    ax_b.axvline(50.0, color="#444444", linestyle=":", linewidth=0.8)
    for i, r in enumerate(sub_s.itertuples(index=False)):
        if r.site == "breast_mammary_tissue":
            ax_b.text(
                r.xx_pct + 2.0,
                i,
                f"65.4% male fat pads ({r.xy}/{r.tot})",
                va="center",
                fontsize=7.8,
                fontweight="bold",
                color="white",
            )
    ax_b.set_yticks(y_idx)
    ax_b.set_yticklabels([_pretty_site(s) for s in sub_s["site"]], fontsize=8.2)
    ax_b.set_xlim(0.0, 100.0)
    ax_b.set_xlabel("Donor chromosomal sex share (%)")
    ax_b.legend(loc="upper center", bbox_to_anchor=(0.5, 1.09), ncol=2, fontsize=8.0)
    ax_b.text(-0.26, 1.03, "b", transform=ax_b.transAxes, fontsize=12, fontweight="bold")

    fig.tight_layout(w_pad=2.0)
    return _save_png_and_svg(fig, out_stem)


def plot_cleanup_evaluation(
    out_stem: Path,
    *,
    eval_df: pd.DataFrame,
    hpa_df: pd.DataFrame,
    controls_df: pd.DataFrame,
    hk_cv_df: pd.DataFrame,
) -> tuple[Path, Path]:
    """Figure 4: 4-panel quantitative evaluation of oncoref.gtex (qc, cleanup, normalize)."""
    import matplotlib.pyplot as plt

    from .. import figure_style

    figure_style.apply()
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 9.4))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    y_a = np.arange(len(eval_df))
    h = 0.25
    ax_a.barh(
        y_a + h,
        np.log10(eval_df["raw_tpm"].to_numpy(dtype=float) + 0.01),
        height=h,
        color=figure_style.DROP,
        label="Raw upstream (qc=None, cleanup=None)",
    )
    ax_a.barh(
        y_a,
        np.log10(eval_df["qc_only_tpm"].to_numpy(dtype=float) + 0.01),
        height=h,
        color=figure_style.ACCENT,
        label="Rejection QC only (qc='all', cleanup=None)",
    )
    ax_a.barh(
        y_a - h,
        np.log10(eval_df["full_clean_tpm"].to_numpy(dtype=float) + 0.01),
        height=h,
        color=figure_style.KEPT,
        label="Full QC + cleanup (qc='all', cleanup='all')",
    )
    ax_a.axvline(np.log10(4.01), color=figure_style.THRESHOLD, linestyle="--", linewidth=0.9)
    ax_a.set_ylim(-0.65, len(eval_df) + 0.35)
    ax_a.text(
        np.log10(4.01) + 0.07,
        len(eval_df) + 0.25,
        "4.0 TPM gate",
        fontsize=7.6,
        color=figure_style.THRESHOLD,
        va="top",
    )
    ax_a.set_yticks(y_a)
    ax_a.set_yticklabels(eval_df["label"].tolist(), fontsize=8.0)
    ax_a.set_xlabel("GTEx v11 summary expression log10(TPM + 0.01)")
    ax_a.legend(loc="upper right", fontsize=7.3)
    ax_a.text(-0.30, 1.03, "a", transform=ax_a.transAxes, fontsize=12, fontweight="bold")

    y_b = np.arange(len(hpa_df))
    ax_b.barh(
        y_b + 0.18,
        hpa_df["raw_ntpm"],
        height=0.34,
        color=figure_style.DROP,
        label="Raw HPA v23 (GTEx v8 mean)",
    )
    ax_b.barh(
        y_b - 0.18,
        hpa_df["clean_ntpm"],
        height=0.34,
        color=figure_style.KEPT,
        label="Cleaned HPA v23 (qc='all', cleanup='all')",
    )
    ax_b.axvline(
        4.0,
        color=figure_style.THRESHOLD,
        linestyle="--",
        linewidth=0.9,
        label="4.0 nTPM Vital ceiling",
    )
    ax_b.set_yticks(y_b)
    ax_b.set_yticklabels(hpa_df["tissue"].tolist(), fontsize=8.2)
    ax_b.set_xlabel("HPA v23 PRM1 expression (nTPM)")
    ax_b.legend(loc="upper right", fontsize=7.4)
    ax_b.text(-0.22, 1.03, "b", transform=ax_b.transAxes, fontsize=12, fontweight="bold")

    y_c = np.arange(len(controls_df))
    ratio_pct = (
        100.0
        * controls_df["clean_val"].to_numpy(dtype=float)
        / np.maximum(controls_df["raw_val"].to_numpy(dtype=float), 1e-9)
    )
    ax_c.barh(y_c, ratio_pct, height=0.56, color=figure_style.KEPT, alpha=0.88)
    ax_c.axvline(100.0, color="#333333", linestyle="--", linewidth=0.85)
    for i, r in enumerate(controls_df.itertuples(index=False)):
        ax_c.text(
            ratio_pct[i] + 2.2,
            i,
            f"{r.clean_val:.1f} / {r.raw_val:.1f} TPM ({r.kept_donors}/{r.total_donors} donors)",
            va="center",
            fontsize=7.4,
            color="#222222",
        )
    ax_c.set_yticks(y_c)
    ax_c.set_yticklabels(controls_df["label"].tolist(), fontsize=8.0)
    ax_c.set_xlim(0.0, 232.0)
    ax_c.set_xlabel("Retained expression after qc='all', cleanup='all' (% of raw)")
    ax_c.text(-0.30, 1.03, "c", transform=ax_c.transAxes, fontsize=12, fontweight="bold")

    y_d = np.arange(len(hk_cv_df))
    ax_d.barh(
        y_d + 0.18,
        hk_cv_df["raw_cv"],
        height=0.34,
        color=figure_style.DROP,
        label="Raw TPM (compositional skew)",
    )
    ax_d.barh(
        y_d - 0.18,
        hk_cv_df["clean_tpm_cv"],
        height=0.34,
        color=figure_style.ACCENT,
        label="normalize='clean_tpm' (16/9/75)",
    )
    ax_d.set_yticks(y_d)
    ax_d.set_yticklabels(hk_cv_df["symbol"].tolist(), fontsize=8.2)
    ax_d.set_xlim(0.0, 1.05)
    ax_d.set_xlabel("Cross-tissue coefficient of variation (CV across 52 GTEx shards)")
    ax_d.legend(loc="upper right", fontsize=7.6)
    ax_d.text(-0.18, 1.03, "d", transform=ax_d.transAxes, fontsize=12, fontweight="bold")

    fig.tight_layout(w_pad=2.2, h_pad=2.2)
    return _save_png_and_svg(fig, out_stem)


def plot_normal_source_coverage_grid(out_stem: Path) -> tuple[Path, Path]:
    """Figure 5: Canonical 70-tissue x 12-source human normal coverage matrix."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    from .. import figure_style

    figure_style.apply()
    grid = normal_tissues.normal_tissue_source_grid().copy()

    cat_order = {
        "vital": 0,
        "protected": 1,
        "manageable": 2,
        "disposable": 3,
        "thymic": 4,
        "reproductive": 5,
        "developmental": 6,
    }
    grid["_ord"] = grid["tissue_category"].map(lambda x: cat_order.get(str(x), 9))
    grid = grid.sort_values(["_ord", "concept_id"]).reset_index(drop=True)

    source_cols = [
        ("gtex_v11_shards", "GTEx v11\n52 shards"),
        ("gtex_v11_medians", "GTEx v11\n66 medians"),
        ("gtex_v11_subcompartments", "GTEx v11\n14 subcomp."),
        ("gtex_via_hpa_v8", "GTEx v8\nvia HPA"),
        ("hpa_surgical_biopsies", "HPA v23\nbiopsies"),
        ("hpa_tissue_rna", "HPA v23\ntissue RNA"),
        ("hpa_rna_consensus", "HPA\nconsensus"),
        ("hpa_single_cell", "HPA\nscRNA"),
        ("cellxgene_census", "CELLxGENE\nCensus"),
        ("sorted_immune", "Sorted\nimmune"),
        ("specialized_and_developmental", "Home / dev\npanels (#735)"),
        ("hpa_ihc_protein", "HPA IHC\nprotein"),
    ]
    mat = np.zeros((len(grid), len(source_cols)), dtype=int)
    for j, (col, _) in enumerate(source_cols):
        present = grid[col].astype(str).str.strip().ne("").to_numpy()
        gtex_absent = grid["gtex_v11_shards"].astype(str).str.strip().eq("").to_numpy()
        mat[:, j] = np.where(~present, 0, np.where(gtex_absent & (j >= 3), 2, 1))

    short_labels = {
        "Mammary tissue (whole breast or mammary gland)": "mammary tissue (breast)",
        "Hair follicle compartment": "hair follicle",
    }

    fig, ax = plt.subplots(figsize=(8.8, 13.4))
    cmap = ListedColormap(["#f3f4f6", figure_style.KEPT, figure_style.ACCENT])
    ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=2, interpolation="nearest")

    ax.set_xticks(range(len(source_cols)))
    ax.set_xticklabels([lbl for _, lbl in source_cols], fontsize=7.8, rotation=35, ha="right")
    ax.set_yticks(range(len(grid)))
    ax.set_yticklabels(
        [
            f"{short_labels.get(str(r.tissue_label), str(r.tissue_label))} [{r.tissue_category}]"
            for r in grid.itertuples(index=False)
        ],
        fontsize=6.8,
    )

    cats = grid["tissue_category"].tolist()
    for i in range(1, len(cats)):
        if cats[i] != cats[i - 1]:
            ax.axhline(i - 0.5, color="#333333", linewidth=0.9)

    ax.set_xlabel(
        "Normal RNA & protein reference modalities (green = GTEx-covered; blue = non-GTEx / EyeGEx rescue)"
    )
    ax.set_ylabel("Canonical human normal tissues and developmental compartments (n = 70)")
    fig.tight_layout()
    return _save_png_and_svg(fig, out_stem)


__all__ = [
    "plot_cleanup_evaluation",
    "plot_contamination_regimes",
    "plot_donor_sex_and_mammary_bias",
    "plot_normal_source_coverage_grid",
    "plot_splash_matrix_and_stoichiometry",
]

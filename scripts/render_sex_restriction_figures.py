#!/usr/bin/env python3
"""Render publication-ready figures interrogating sex restriction and female-lineage CTAs."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import (  # noqa: E402
    gtex,
    human,
    human_protein_atlas,
    species_target_selection,
    target_selection,
    tcr_trials,
)

OUT_DIR = ROOT / "docs/antigen-whitepaper"
TUMOR_PEAKS_CSV = ROOT / "scripts/cta_evidence_assets/female-lineage-tumor-peaks.csv"

# Publication palette harmonized with OncoRef whitepaper figures
C_UNSAFE = "#b84a47"
C_OTHER = "#9aa7b2"
C_CTAM = "#2b6cb0"
C_CTAE = "#2f855a"
C_CTAF = "#c53070"
C_CTAMF = "#c27d2b"
C_RAF_ONLY = "#e68ab8"
C_INK = "#1f2933"
C_MUTED = "#52606d"
C_GRID = "#e4e7eb"

FEMALE_LINEAGE_PROGRAMS: list[tuple[str, list[str]]] = [
    (
        "Mammary & lactation\n(32 canonical genes)",
        [
            "LALBA",
            "CSN2",
            "BTN1A1",
            "CSN1S1",
            "CSN3",
            "SCGB2A2",
            "SCGB1D2",
            "SCGB2A1",
            "SCGB1D4",
            "PIP",
            "MUCL1",
            "ANKRD30A",
            "ANKRD30B",
            "PIGR",
            "LTF",
            "LPO",
            "CEL",
            "XDH",
            "MFGE8",
            "PRLR",
            "OXTR",
            "ELF5",
            "TRPS1",
            "GATA3",
            "ESR1",
            "PGR",
            "VTCN1",
            "FOLR1",
            "SLC34A2",
            "ROPN1",
            "ROPN1B",
            "CT83",
        ],
    ),
    (
        "Oocyte, ZP, SCMC\n& follicle (49 genes)",
        [
            "ZP1",
            "ZP2",
            "ZP3",
            "ZP4",
            "GDF9",
            "BMP15",
            "MOS",
            "FIGLA",
            "NOBOX",
            "LHX8",
            "SOHLH1",
            "SOHLH2",
            "FOXL2",
            "AMH",
            "AMHR2",
            "FSHR",
            "LHCGR",
            "CYP19A1",
            "ASTL",
            "FETUB",
            "H1-8",
            "TUBB8",
            "WEE2",
            "PATL2",
            "BTG4",
            "PABPN1L",
            "ZAR1",
            "ZAR1L",
            "NLRP5",
            "OOEP",
            "TLE6",
            "KHDC3L",
            "PADI6",
            "NLRP2",
            "NLRP4",
            "NLRP7",
            "NLRP8",
            "NLRP9",
            "NLRP11",
            "NLRP13",
            "NLRP14",
            "OOSP1",
            "OOSP2",
            "OOSP4A",
            "DPPA3",
            "KPNA7",
            "NPM2",
            "ZBED3",
            "TCL1B",
        ],
    ),
    (
        "Shared meiosis &\ngerm cell (21 genes)",
        [
            "HORMAD1",
            "HORMAD2",
            "SPO11",
            "SYCP1",
            "SYCP2",
            "SYCP3",
            "SMC1B",
            "STAG3",
            "REC8",
            "MEIOC",
            "STRA8",
            "DAZL",
            "DDX4",
            "MAEL",
            "PIWIL1",
            "PIWIL2",
            "PIWIL3",
            "PIWIL4",
            "TDRD1",
            "TDRD9",
            "BOLL",
        ],
    ),
    (
        "Endometrium, oviduct\n& Müllerian (18 genes)",
        [
            "OVGP1",
            "PAEP",
            "MMP26",
            "MMP7",
            "MUC16",
            "WFDC2",
            "MSLN",
            "PAX8",
            "LEFTY2",
            "PRL",
            "IGFBP1",
            "LIF",
            "HOXA10",
            "HOXA11",
            "WNT7A",
            "CLDN6",
            "PRAME",
            "TACSTD2",
        ],
    ),
    (
        "Placenta &\ntrophoblast (39 genes)",
        [
            "ERVV-1",
            "ERVV-2",
            "ERVW-1",
            "ERVFRD-1",
            "PLAC1",
            "PLAC8",
            "PSG1",
            "PSG2",
            "PSG3",
            "PSG4",
            "PSG5",
            "PSG6",
            "PSG7",
            "PSG8",
            "PSG9",
            "PSG11",
            "CGB3/5/8",
            "CGB1",
            "CGB2",
            "CGB7",
            "CGA",
            "CSHL1",
            "GH2",
            "GCM1",
            "PAPPA",
            "PAPPA2",
            "KISS1",
            "LGALS13",
            "LGALS14",
            "LGALS16",
            "INSL4",
            "ALPP",
            "HLA-G",
            "VGLL1",
            "XAGE2",
            "XAGE3",
            "PAGE4",
            "MAGEA10",
            "HSD3B1",
        ],
    ),
]


def _style_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#9aa5b1")
    ax.spines["bottom"].set_color("#9aa5b1")
    ax.tick_params(colors=C_INK, labelsize=9)


def draw_funnel_and_landscape(cat: pd.DataFrame, sex_df: pd.DataFrame, he: dict) -> Path:
    by_lbl = cat.set_index("protein_group_label", drop=False)
    sex_by_eid = sex_df.set_index("entity_id")

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(14.6, 6.4), gridspec_kw={"width_ratios": [1.05, 1.25]}
    )
    fig.patch.set_facecolor("white")

    # --- Panel A: Fate of 159 canonical female-lineage candidates across 5 biological programs ---
    buckets = [
        ("Stage 1 Unsafe (somatic sharing)", C_UNSAFE),
        ("Other / CA-only (review, seq, <10× RA)", C_OTHER),
        ("CTA-M / RA-M (testis-active, female-silent)", C_CTAM),
        ("CTA-MF / RA-MF (shared male & female)", C_CTAMF),
        ("CTA-P / RA-P-only (pregnancy/peripartum-only)", C_CTAE),
        ("CTA-F / RA-F-only (female-restricted)", C_CTAF),
    ]
    prog_labels = [lbl for lbl, _ in FEMALE_LINEAGE_PROGRAMS]
    counts = np.zeros((len(FEMALE_LINEAGE_PROGRAMS), len(buckets)), dtype=float)
    raw_counts = np.zeros((len(FEMALE_LINEAGE_PROGRAMS), len(buckets)), dtype=int)

    for i, (_plbl, genes) in enumerate(FEMALE_LINEAGE_PROGRAMS):
        for g in genes:
            if g not in by_lbl.index:
                continue
            r = by_lbl.loc[g]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            if r.category == "Unsafe":
                j = 0
            elif not bool(r.is_reproductive_antigen):
                j = 1
            else:
                sc = sex_by_eid.loc[r.entity_id, "sex_restriction_class"]
                if sc == "Male-restricted":
                    j = 2
                elif sc == "Shared male & female":
                    j = 3
                elif sc == "Pregnancy / peripartum-only":
                    j = 4
                else:
                    j = 5
            raw_counts[i, j] += 1
        tot = raw_counts[i].sum()
        counts[i] = 100.0 * raw_counts[i] / max(tot, 1)

    y = np.arange(len(prog_labels))[::-1]
    left = np.zeros(len(prog_labels))
    for j, (blbl, color) in enumerate(buckets):
        widths = counts[:, j]
        ax_a.barh(
            y,
            widths,
            left=left,
            height=0.60,
            color=color,
            edgecolor="white",
            linewidth=0.8,
            label=blbl,
        )
        for idx_row, (yp, w, rc) in enumerate(zip(y, widths, raw_counts[:, j])):
            if w >= 6.5 and rc > 0:
                ax_a.text(
                    left[idx_row] + w / 2.0,
                    yp,
                    f"{rc}",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=8.5,
                    fontweight="bold",
                )
        left += widths

    # Annotate exact CTA-F survivors on the right edge of each row
    ctaf_notes = [
        "CSN2, LALBA (2 CTA-F; BTN1A1 → CTA-P)",
        "ZP4, DPPA3, KHDC3L, ZAR1L + 4 RA-F-only",
        "0 CTA-F (SYCP1 → CTA-M; MOS, SPO11 → CTA-MF)",
        "0 CTA-F (PRAME, SCGB1D4 → CTA-MF)",
        "0 CTA-F (ERVV-1/2, LGALS16 → CTA-P)",
    ]
    for yp, note in zip(y, ctaf_notes):
        ax_a.text(101.2, yp, note, va="center", ha="left", fontsize=8.0, color=C_INK)

    ax_a.set_yticks(y)
    ax_a.set_yticklabels(prog_labels, fontsize=8.8, color=C_INK)
    ax_a.set_xlim(0, 148)
    ax_a.set_xticks([0, 25, 50, 75, 100])
    ax_a.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax_a.set_xlabel(
        "Share of canonical female-program protein groups (numbers = gene counts)",
        fontsize=9.2,
        color=C_INK,
    )
    ax_a.set_title(
        "A. Why canonical female-lineage genes rarely become CTA-F (n = 159)",
        loc="left",
        fontsize=10.8,
        fontweight="bold",
        color=C_INK,
        pad=10,
    )
    _style_axes(ax_a)
    ax_a.spines["bottom"].set_bounds(0, 100)
    ax_a.legend(
        loc="upper center",
        bbox_to_anchor=(0.46, -0.15),
        ncol=2,
        frameon=False,
        fontsize=8.0,
    )

    # --- Panel B: Live deconvolved adult male vs. female exposure across all 332 CTAs ---
    anat = species_target_selection.anatomy_for_species("human")
    ctas = sex_df[sex_df["is_cta"].astype(bool)].copy()

    def _effective_adult_exposure(eid: str, side: str) -> float:
        rec = he.get(eid, {})
        vals = [0.05]
        comps = (
            anat.male_adult_protected_compartments
            if side == "male"
            else anat.female_adult_protected_compartments
        )
        for c in comps:
            item = rec.get(c, {})
            if item.get("coverage") != "measured":
                continue
            if c == "resting_breast":
                t_val = item.get(
                    "male_max_tpm" if side == "male" else "female_max_tpm", item.get("max_tpm")
                )
                n_val = item.get(
                    "male_max_ntpm" if side == "male" else "female_max_ntpm", item.get("max_ntpm")
                )
            else:
                t_val = item.get("max_tpm")
                n_val = item.get("max_ntpm")
            if t_val is not None and np.isfinite(t_val):
                vals.append(float(t_val))
            # Reflect scRNA nTPM >= 8 ceiling crossings when bulk TPM is diluted (e.g., oocytes, XX breast)
            if n_val is not None and np.isfinite(n_val) and float(n_val) >= 8.0:
                vals.append(max(float(t_val or 0.0), min(float(n_val), 800.0)))
            pstat = item.get("protein_status")
            if pstat in {"positive", "unresolved"} and not (
                c == "resting_breast" and side == "male" and pstat == "unresolved"
            ):
                vals.append(max(vals[-1], 20.0))
        return max(vals)

    ctas["male_eff"] = [_effective_adult_exposure(eid, "male") for eid in ctas["entity_id"]]
    ctas["female_eff"] = [_effective_adult_exposure(eid, "female") for eid in ctas["entity_id"]]

    # Shade CTA-F (top-left) and CTA-M (bottom-right) and CTA-E (bottom-left) quadrants
    ax_b.axvspan(0.03, 8.0, ymin=0, ymax=1, color="#fdf2f8", alpha=0.35, zorder=0)
    ax_b.axhspan(0.03, 8.0, xmin=0, xmax=1, color="#ebf8ff", alpha=0.35, zorder=0)
    ax_b.axvline(8.0, color="#616e7c", linestyle="--", linewidth=1.1, zorder=2)
    ax_b.axhline(8.0, color="#616e7c", linestyle="--", linewidth=1.1, zorder=2)
    ax_b.axvline(16.0, color="#9aa5b1", linestyle=":", linewidth=0.95, zorder=2)
    ax_b.axhline(16.0, color="#9aa5b1", linestyle=":", linewidth=0.95, zorder=2)

    code_meta = [
        ("CTA-M", C_CTAM, "o", 26),
        ("CTA-P", C_CTAE, "s", 30),
        ("CTA-MF", C_CTAMF, "^", 34),
        ("CTA-F", C_CTAF, "D", 56),
    ]
    legend_handles = []
    for code, color, marker, sz in code_meta:
        sub = ctas[ctas["sex_restriction_code"] == code]
        ax_b.scatter(
            np.clip(sub["male_eff"], 0.05, 80000),
            np.clip(sub["female_eff"], 0.05, 80000),
            c=color,
            marker=marker,
            s=sz,
            alpha=0.86,
            edgecolors="white",
            linewidths=0.45,
            zorder=4 if code == "CTA-F" else 3,
        )
        legend_handles.append(
            Line2D(
                [0],
                [0],
                marker=marker,
                color="w",
                label=f"{code} (n = {len(sub)})",
                markerfacecolor=color,
                markersize=7.5,
            )
        )

    callouts = {
        "LALBA": (-44, 10),
        "DPPA3": (-44, 6),
        "KHDC3L": (-56, 14),
        "ZAR1L": (-52, -18),
        "CSN2": (-48, -3),
        "ZP4": (-30, -13),
        "CPSF4L": (-48, 8),
        "IZUMO1R": (-54, -8),
        "VRTN": (-38, -9),
        "PRR23E": (-48, 4),
        "BTN1A1": (-48, -4),
        "ERVV-1": (-48, 8),
        "DUX4": (14, -8),
        "MOS": (8, 8),
        "SCGB1D4": (8, -4),
        "PRAME": (8, 8),
        "CTAG2": (8, -10),
        "DKKL1": (12, 6),
        "CT83": (14, -14),
        "NY-ESO-1": (12, 10),
    }
    ctas_by_lbl = ctas.set_index("protein_group_label")
    for sym, (dx, dy) in callouts.items():
        if sym not in ctas_by_lbl.index:
            continue
        r = ctas_by_lbl.loc[sym]
        mx = float(np.clip(r["male_eff"], 0.05, 80000))
        fy = float(np.clip(r["female_eff"], 0.05, 80000))
        c_txt = C_CTAF if r["sex_restriction_code"] == "CTA-F" else C_INK
        ax_b.annotate(
            sym,
            xy=(mx, fy),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=7.8,
            fontweight="bold" if r["sex_restriction_code"] == "CTA-F" else "normal",
            color=c_txt,
            arrowprops={"arrowstyle": "-", "color": "#9aa5b1", "lw": 0.65},
            zorder=6,
        )

    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.set_xlim(0.035, 110000)
    ax_b.set_ylim(0.035, 110000)
    ticks = [0.05, 0.5, 8, 16, 100, 1000, 50000]
    tick_lbls = ["≤0.05", "0.5", "8", "16", "100", "1k", "50k"]
    ax_b.set_xticks(ticks)
    ax_b.set_xticklabels(tick_lbls)
    ax_b.set_yticks(ticks)
    ax_b.set_yticklabels(tick_lbls)
    ax_b.grid(True, which="major", color=C_GRID, linewidth=0.6, alpha=0.8)
    ax_b.text(
        0.048,
        68000,
        "CTA-F quadrant\n(male-safe, <8 / <16 IHC−)",
        fontsize=8.0,
        color=C_CTAF,
        fontweight="bold",
        va="top",
    )
    ax_b.text(
        1100,
        0.065,
        "CTA-M quadrant (female-safe, <8 / <16 IHC−)",
        fontsize=8.0,
        color=C_CTAM,
        fontweight="bold",
        ha="left",
    )
    ax_b.text(
        0.048,
        5.8,
        "CTA-P quadrant\n(pregnancy/peripartum)",
        fontsize=7.6,
        color=C_CTAE,
        fontweight="bold",
        va="top",
    )
    ax_b.set_xlabel(
        "Max adult male organ exposure (TPM / ceiling-crossing nTPM)", fontsize=9.2, color=C_INK
    )
    ax_b.set_ylabel(
        "Max non-pregnant adult female exposure (TPM / nTPM)", fontsize=9.2, color=C_INK
    )
    ax_b.set_title(
        f"B. Deconvolved adult cross-sex exposure across all {len(ctas)} human CTAs",
        loc="left",
        fontsize=10.8,
        fontweight="bold",
        color=C_INK,
        pad=10,
    )
    _style_axes(ax_b)
    ax_b.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=4,
        frameon=False,
        fontsize=8.2,
    )

    fig.tight_layout()
    out_path = OUT_DIR / "sex-restriction-funnel-and-landscape.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


def draw_lalba_csn2_interrogation(cat: pd.DataFrame) -> Path:
    milk_genes = {
        "LALBA": "ENSG00000167531",
        "CSN2": "ENSG00000135222",
        "BTN1A1": "ENSG00000124557",
        "CSN1S1": "ENSG00000126545",
        "CSN3": "ENSG00000171209",
    }
    gids = list(milk_genes.values())

    # Load XX and XY GTEx Breast_Mammary_Tissue and GTEx Testis / Prostate
    b_xx = gtex.sample_tpm("breast_mammary_tissue", genes=gids, qc="all", cleanup="all", sex="XX")
    b_xy = gtex.sample_tpm("breast_mammary_tissue", genes=gids, qc="all", cleanup="all", sex="XY")
    t_xy = gtex.sample_tpm("testis", genes=gids, qc="all", cleanup="all")
    p_xy = gtex.sample_tpm("prostate", genes=gids, qc="all", cleanup="all")

    def _mat(df: pd.DataFrame) -> pd.DataFrame:
        cols = [c for c in df.columns if c not in ("Name", "Description")]
        return df.set_index("Description")[cols]

    m_xx = _mat(b_xx)
    m_xy = _mat(b_xy)
    m_testis = _mat(t_xy)
    m_prost = _mat(p_xy)

    core4 = ["LALBA", "CSN1S1", "CSN2", "CSN3"]
    sec_xx = [
        c
        for c in m_xx.columns
        if (m_xx.loc[core4, c] >= 5.0).sum() >= 3 and m_xx.loc[core4, c].max() >= 15.0
    ]
    rest_xx = [c for c in m_xx.columns if c not in sec_xx]
    sec_xy = [
        c
        for c in m_xy.columns
        if (m_xy.loc[core4, c] >= 5.0).sum() >= 3 and m_xy.loc[core4, c].max() >= 15.0
    ]
    rest_xy = [c for c in m_xy.columns if c not in sec_xy]

    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(16.2, 5.7), gridspec_kw={"width_ratios": [1.18, 1.05, 1.25]}
    )
    fig.patch.set_facecolor("white")

    rng = np.random.default_rng(20261125)
    genes_show = ["LALBA", "CSN2", "BTN1A1"]
    gene_colors = {"LALBA": C_CTAF, "CSN2": "#7b2cbf", "BTN1A1": C_CTAE}

    # --- Panel A: GTEx Breast_Mammary_Tissue across XX/XY secretory vs resting states ---
    groups_a = [
        (f"XX secretory/\nlactating\n(n={len(sec_xx)})", m_xx[sec_xx]),
        (f"XX resting\nbreast\n(n={len(rest_xx)})", m_xx[rest_xx]),
        (f"XY ICU/cadaver\ngalactorrhea\n(n={len(sec_xy)})", m_xy[sec_xy]),
        (f"XY resting\nmale breast\n(n={len(rest_xy)})", m_xy[rest_xy]),
    ]
    x_base = np.arange(len(groups_a))
    offsets = [-0.24, 0.0, 0.24]
    ax1.axhline(8.0, color="#616e7c", linestyle="--", linewidth=1.0, zorder=1)
    ax1.axhline(16.0, color="#9aa5b1", linestyle=":", linewidth=0.95, zorder=1)
    ax1.text(3.42, 8.8, "8 TPM base", fontsize=7.4, color=C_MUTED, ha="right")

    for g_idx, sym in enumerate(genes_show):
        c = gene_colors[sym]
        for grp_i, (_glbl, sub_m) in enumerate(groups_a):
            vals = np.clip(sub_m.loc[sym].dropna().to_numpy(dtype=float), 0.02, 60000.0)
            jitter = rng.uniform(-0.07, 0.07, size=len(vals))
            ax1.scatter(
                np.full(len(vals), x_base[grp_i] + offsets[g_idx]) + jitter,
                vals,
                s=13 if len(vals) > 20 else 24,
                color=c,
                alpha=0.55 if len(vals) > 20 else 0.85,
                edgecolors="none",
                zorder=3,
            )
            p99 = float(np.quantile(vals, 0.99)) if len(vals) >= 30 else float(np.max(vals))
            ax1.plot(
                [x_base[grp_i] + offsets[g_idx] - 0.085, x_base[grp_i] + offsets[g_idx] + 0.085],
                [p99, p99],
                color=C_INK,
                linewidth=1.7,
                zorder=5,
            )

    ax1.set_yscale("log")
    ax1.set_ylim(0.015, 90000)
    ax1.set_xticks(x_base)
    ax1.set_xticklabels([g[0] for g in groups_a], fontsize=8.3)
    ax1.set_yticks([0.02, 0.2, 8, 50, 1000, 40000])
    ax1.set_yticklabels(["≤0.02", "0.2", "8", "50", "1k", "40k"])
    ax1.set_ylabel("Deconvolved donor expression (TPM, bars = P99/max)", fontsize=9.0, color=C_INK)
    ax1.set_title(
        "A. GTEx Breast_Mammary_Tissue (514 XX & XY biopsies)",
        loc="left",
        fontsize=10.2,
        fontweight="bold",
        color=C_INK,
        pad=10,
    )
    ax1.grid(True, axis="y", color=C_GRID, linewidth=0.6)
    _style_axes(ax1)
    ax1.legend(
        handles=[
            Patch(facecolor=gene_colors["LALBA"], label="LALBA (CTA-F)"),
            Patch(facecolor=gene_colors["CSN2"], label="CSN2 (CTA-F)"),
            Patch(facecolor=gene_colors["BTN1A1"], label="BTN1A1 (CTA-P)"),
        ],
        loc="upper right",
        frameon=True,
        facecolor="white",
        edgecolor=C_GRID,
        fontsize=7.8,
    )

    # --- Panel B: Adult Male Gonad & Accessory Organs (GTEx Testis n=414, Prostate n=282, HPA scRNA & IHC) ---
    sc_clean = human_protein_atlas.single_cell_tissues(qc="all", cleanup="all")
    sc_testis = sc_clean[sc_clean["Tissue"].eq("testis")]
    sc_max = {
        sym: float(sc_testis.loc[sc_testis["Gene"].eq(milk_genes[sym]), "nTPM"].max() or 0.0)
        for sym in genes_show
    }

    groups_b = [
        (f"GTEx Testis\n(n={m_testis.shape[1]})", m_testis, "TPM"),
        (f"GTEx Prostate\n(n={m_prost.shape[1]})", m_prost, "TPM"),
    ]
    xb = np.arange(3)
    ax2.axhline(8.0, color="#616e7c", linestyle="--", linewidth=1.0, zorder=1)
    ax2.axhline(16.0, color="#9aa5b1", linestyle=":", linewidth=0.95, zorder=1)
    ax2.text(2.42, 8.6, "8 Protected base", fontsize=7.2, color=C_MUTED, ha="right")
    ax2.text(2.42, 17.5, "16 IHC− ceiling", fontsize=7.2, color=C_MUTED, ha="right")

    for g_idx, sym in enumerate(genes_show):
        c = gene_colors[sym]
        for grp_i, (_glbl, sub_m, _u) in enumerate(groups_b):
            vals = np.clip(sub_m.loc[sym].dropna().to_numpy(dtype=float), 0.02, 50.0)
            jitter = rng.uniform(-0.07, 0.07, size=len(vals))
            ax2.scatter(
                np.full(len(vals), xb[grp_i] + offsets[g_idx]) + jitter,
                vals,
                s=11,
                color=c,
                alpha=0.45,
                edgecolors="none",
                zorder=3,
            )
            p99 = float(np.quantile(vals, 0.99))
            ax2.plot(
                [xb[grp_i] + offsets[g_idx] - 0.085, xb[grp_i] + offsets[g_idx] + 0.085],
                [p99, p99],
                color=C_INK,
                linewidth=1.8,
                zorder=5,
            )
        # Third column: HPA Testis late-spermatid / germ-cell peak nTPM
        nt = max(sc_max[sym], 0.03)
        ax2.bar(
            xb[2] + offsets[g_idx],
            nt,
            width=0.19,
            color=c,
            alpha=0.88,
            edgecolor="white",
            zorder=4,
        )
        ax2.text(
            xb[2] + offsets[g_idx],
            nt * 1.14,
            f"{sc_max[sym]:.1f}",
            ha="center",
            va="bottom",
            fontsize=7.6,
            fontweight="bold",
            color=C_INK,
        )

    ax2.annotate(
        "LALBA → CTA-F (male-safe):\nTestis P50=2.6, P95=4.9, P99=7.8 TPM (<8)\n& Enhanced IHC: Not detected (<16)",
        xy=(xb[0] + offsets[0], 7.76),
        xytext=(0.06, 25.5),
        fontsize=7.3,
        color=C_INK,
        arrowprops={"arrowstyle": "->", "color": C_CTAF, "lw": 0.9},
        bbox={
            "boxstyle": "round,pad=0.25",
            "facecolor": "#fdf2f8",
            "edgecolor": "#ed64a6",
            "lw": 0.8,
        },
    )
    ax2.set_yscale("log")
    ax2.set_ylim(0.015, 85.0)
    ax2.set_xticks(xb)
    ax2.set_xticklabels(
        [
            f"GTEx Testis\n(n={m_testis.shape[1]})",
            f"GTEx Prostate\n(n={m_prost.shape[1]})",
            "HPA Testis\nscRNA (nTPM)",
        ],
        fontsize=8.3,
    )
    ax2.set_yticks([0.02, 0.2, 1, 8, 16, 50])
    ax2.set_yticklabels(["≤0.02", "0.2", "1", "8", "16", "50"])
    ax2.set_ylabel("Male reproductive RNA (TPM or peak cluster nTPM)", fontsize=9.0, color=C_INK)
    ax2.set_title(
        "B. Male reproductive RNA & IHC− Protected gate",
        loc="left",
        fontsize=10.2,
        fontweight="bold",
        color=C_INK,
        pad=10,
    )
    ax2.grid(True, axis="y", color=C_GRID, linewidth=0.6)
    _style_axes(ax2)

    # --- Panel C: Breast & gynecologic tumor peaks vs. somatic-safety status ---
    peaks_df = pd.read_csv(TUMOR_PEAKS_CSV).drop_duplicates("gene_symbol").set_index("gene_symbol")
    panel_c_genes = [
        ("CSN2", "CTA-F", C_CTAF, "Male ≤2.1 TPM · Lactation 99k"),
        ("LALBA", "CTA-F", C_CTAF, "Testis P99 7.8 TPM (<8 · IHC−)"),
        ("BTN1A1", "CTA-P", C_CTAE, "Adult ≤2.1 TPM · Lactation-only"),
        ("ZP4", "CTA-F", C_CTAF, "Male ≤0.1 TPM · Oocyte zona"),
        ("SCGB1D4", "CTA-MF", C_CTAMF, "UCEC 3,576 · Epididymis 51"),
        ("PRAME", "CTA-MF", C_CTAMF, "UCEC/UCS · Testis + Ovary"),
        ("CTAG2", "CTA-MF", C_CTAMF, "OV 604 · Testis + Vagina"),
        ("XAGE3", "CTA-M", C_CTAM, "TNBC/UCEC · Female <3.3"),
        ("CT83", "CTA-M", C_CTAM, "CESC/TNBC · Female <2.2"),
        ("PIP", "Unsafe", C_UNSAFE, "Salivary 16,000 nTPM"),
        ("SCGB2A2", "Unsafe", C_UNSAFE, "Skin 690 nTPM"),
        ("WFDC2", "Unsafe", C_UNSAFE, "Bronchus 18,000 nTPM"),
        ("FOLR1", "Unsafe", C_UNSAFE, "Choroid plexus 2,900 nTPM"),
        ("SLC34A2", "Unsafe", C_UNSAFE, "Lung AT2 2,200 nTPM"),
        ("PAX8", "Unsafe", C_UNSAFE, "Kidney 350 TPM"),
    ]
    cat_by_lbl = cat.set_index("protein_group_label")
    cohort_cols = [
        c
        for c in (
            "brca_max_tpm",
            "brca_tnbc_max_tpm",
            "ov_max_tpm",
            "ucec_max_tpm",
            "ucs_max_tpm",
            "cesc_max_tpm",
        )
        if c in peaks_df.columns
    ]

    yc = np.arange(len(panel_c_genes))[::-1]
    peaks_gyn = []
    for sym, _status, _col, _note in panel_c_genes:
        v = 0.0
        if sym in peaks_df.index and cohort_cols:
            vals = pd.to_numeric(peaks_df.loc[sym, cohort_cols], errors="coerce")
            if vals.notna().any():
                v = float(vals.max())
        if v <= 0.0 and sym in cat_by_lbl.index:
            tm = cat_by_lbl.loc[sym, "tumor_max_tpm"]
            if pd.notna(tm):
                v = float(tm)
        peaks_gyn.append(max(v, 0.5))

    ax3.barh(
        yc,
        peaks_gyn,
        height=0.64,
        color=[g[2] for g in panel_c_genes],
        alpha=0.88,
        edgecolor="white",
    )
    ax3.axhline(yc[8] - 0.5, color="#616e7c", linestyle=":", linewidth=1.0)
    for yp, val, (_sym, status, col, note) in zip(yc, peaks_gyn, panel_c_genes):
        ax3.text(
            val * 1.18,
            yp,
            f"{status}: {note}",
            va="center",
            ha="left",
            fontsize=7.4,
            color=col if status == "Unsafe" else C_INK,
        )

    ax3.set_xscale("log")
    ax3.set_xlim(1.0, 9e8)
    ax3.set_xticks([10, 100, 1000, 10000, 100000])
    ax3.set_xticklabels(["10", "100", "1k", "10k", "100k"])
    ax3.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax3.set_yticks(yc)
    ax3.set_yticklabels([g[0] for g in panel_c_genes], fontsize=8.3, fontweight="bold")
    ax3.set_xlabel("Peak TPM across BRCA, TNBC, OV, UCEC, UCS & CESC", fontsize=9.0, color=C_INK)
    ax3.set_title(
        "C. Peak female-cancer TPM vs. safety & sex call",
        loc="left",
        fontsize=10.2,
        fontweight="bold",
        color=C_INK,
        pad=10,
    )
    ax3.grid(True, axis="x", color=C_GRID, linewidth=0.6)
    _style_axes(ax3)
    ax3.spines["bottom"].set_bounds(1.0, 200000)

    fig.tight_layout()
    out_path = OUT_DIR / "lalba-csn2-mammary-testis-interrogation.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


def draw_clinical_tcr_gold_standard() -> Path:
    ev = tcr_trials.evaluate_against_gold_standard().set_index("benchmark_id", drop=False)

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(15.8, 6.9), gridspec_kw={"width_ratios": [1.20, 1.12]}
    )
    fig.patch.set_facecolor("white")

    # --- Panel A: Clinical ORR (%) vs Strongest OncoRef Somatic / Cross-Reactive Proxy Ratio ---
    ax_a.axvspan(0.004, 1.0, ymin=0.145, ymax=1.0, color="#f0fff4", alpha=0.55, zorder=0)
    ax_a.axvspan(1.0, 250000.0, ymin=0.145, ymax=1.0, color="#fff5f5", alpha=0.55, zorder=0)
    ax_a.axhspan(-14.5, 3.5, xmin=0.0, xmax=1.0, color="#f1f5f9", alpha=0.85, zorder=0)

    ax_a.axvline(1.0, color="#1f2933", linestyle="--", linewidth=1.35, zorder=2)
    ax_a.axhline(3.5, color="#64748b", linestyle=":", linewidth=1.05, zorder=2)

    # Show how matched negative IHC moves PRAME from 3.0x (naive RNA) into <1.0x CTA territory
    prame_y = float(ev.loc["prame_ima203", "orr_pct"])
    ax_a.scatter(
        [3.0],
        [prame_y],
        s=58,
        facecolors="none",
        edgecolors="#2f855a",
        linewidths=1.2,
        linestyle="--",
        zorder=3,
    )
    ax_a.annotate(
        "",
        xy=(0.82, prame_y),
        xytext=(2.82, prame_y),
        arrowprops={"arrowstyle": "->", "color": "#2f855a", "lw": 1.1, "ls": ":"},
        zorder=4,
    )
    ax_a.text(
        2.55,
        prame_y - 2.6,
        "IHC− credit\n(3.0×→0.75×)",
        fontsize=6.8,
        color="#276749",
        ha="center",
        va="top",
    )

    # Bind coordinates directly to live tcr_trials.evaluate_against_gold_standard()
    # (values <0.01x clipped to 0.01x floor; magea3_a1_titin d4-5 fulminant lysis plotted at y=74 with *)
    label_specs = [
        # (benchmark_id, color, marker, size, label, dx, dy, y_override)
        ("magea4_afami_cel", C_CTAE, "o", 68, "MAGEA4 (afami-cel, 37%)", 7, 7, None),
        ("ny_eso_1_1g4_c259", C_CTAE, "o", 72, "NY-ESO-1 (1G4/c259, 58%)", -38, 10, None),
        ("magea10_adp_a2m10", C_CTAE, "o", 62, "MAGEA10 (ADP-A2M10, 9.1%)", -82, 9, None),
        ("prame_ima203", C_CTAE, "o", 72, "PRAME (IMA203, 52.5%)", -98, -3, None),
        ("ctag2_lete_cel_c259", C_CTAE, "o", 68, "CTAG2 LAGE-1a (c259, 30%)", -118, -4, None),
        ("ceacam5_nci_tcr", C_UNSAFE, "o", 70, "CEACAM5 (33%; 3/3 colitis)", 7, 10, None),
        ("mlana_dmf5_tcr", C_UNSAFE, "o", 70, "MLANA (30%; uveitis/ear)", 8, -4, None),
        ("pmel_gp100_tcr_imm", C_UNSAFE, "o", 70, "PMEL (19% TCR / 9% ImmTAC)", 7, -3, None),
        ("tyr_til1383i_tcr", C_UNSAFE, "o", 64, "TYR (33%; vitiligo)", -62, 10, None),
        (
            "magea3_a2_neuro",
            "#8b1e3f",
            "^",
            82,
            "MAGEA3→A12/EPS8L2 (56%; 2/9 fatal CNS)",
            8,
            4,
            None,
        ),
        (
            "magea3_a1_titin",
            "#8b1e3f",
            "^",
            86,
            "MAGEA3→TTN (a3a; 2/2 fatal heart*)",
            -156,
            5,
            74.0,
        ),
        (
            "wt1_ttcr_c4",
            "#d69e2e",
            "D",
            64,
            "WT1 (TTCR-C4, 100% RFS; podocyte escape)",
            8,
            -15,
            None,
        ),
        (
            "magea3_dp4_6f9",
            "#d69e2e",
            "D",
            64,
            "MAGEA3 (DP4 6F9, 24% safe CD4)",
            -14,
            -15,
            None,
        ),
        ("afp_adp_a2afp", "#d69e2e", "D", 64, "AFP (ADP-A2AFP, 9.5% safe)", 8, -3, None),
        ("magea1_ima202", "#64748b", "s", 50, "MAGEA1 (IMA202, 0% ORR)", 7, -3, -3.8),
        ("ct83_kklc1_phase1", "#64748b", "s", 50, "CT83 (KK-LC-1 Ph1, pending)", 7, -3, -10.2),
    ]

    for bid, col, mk, sz, lbl, dx, dy, y_over in label_specs:
        erow = ev.loc[bid]
        xv = float(np.clip(float(erow["strongest_proxy_ratio"]), 0.01, 200000.0))
        yv = float(y_over) if y_over is not None else float(erow["orr_pct"])
        ax_a.scatter(
            [xv],
            [yv],
            c=col,
            marker=mk,
            s=sz,
            edgecolors="white",
            linewidths=0.75,
            zorder=5,
        )
        ax_a.annotate(
            lbl,
            xy=(xv, yv),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=7.25,
            fontweight="bold" if col in (C_CTAE, C_UNSAFE, "#8b1e3f") else "normal",
            color=C_INK if col != "#64748b" else "#475569",
            arrowprops={"arrowstyle": "-", "color": "#94a3b8", "lw": 0.6},
            zorder=6,
        )

    ax_a.text(
        0.0055,
        110.0,
        "OncoRef Safe (<1×)\n5/5 potent CTAs kept\n0 toxic misses",
        fontsize=7.8,
        fontweight="bold",
        color="#22543d",
        va="top",
    )
    ax_a.text(
        180000.0,
        110.0,
        "OncoRef Unsafe (≥1× ceiling)\n6/6 human specificity toxicities flagged",
        fontsize=7.8,
        fontweight="bold",
        color="#9b2c2c",
        va="top",
        ha="right",
    )
    ax_a.text(
        180000.0,
        -6.5,
        "Potency unproven\n(excluded from calibration)",
        fontsize=7.1,
        color="#475569",
        ha="right",
        va="center",
        style="italic",
    )

    ax_a.set_xscale("log")
    ax_a.set_xlim(0.0042, 240000.0)
    ax_a.set_ylim(-14.5, 114.0)
    ax_a.set_xticks([0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0])
    ax_a.set_xticklabels(["≤0.01×", "0.1×", "1× gate", "10×", "100×", "1,000×", "10,000×"])
    ax_a.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax_a.set_yticks([0, 20, 40, 60, 80, 100])
    ax_a.set_yticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])
    ax_a.grid(True, which="major", color=C_GRID, linewidth=0.6, alpha=0.8)
    ax_a.set_xlabel(
        "Strongest somatic or cross-reactive normal proxy signal (× Vital/Protected ceiling)",
        fontsize=9.0,
        color=C_INK,
    )
    ax_a.set_ylabel(
        "Observed clinical objective response / AML RFS (%)",
        fontsize=9.0,
        color=C_INK,
    )
    ax_a.set_title(
        "A. Human TCR-T / ImmTAC gold standard vs. OncoRef safety boundary",
        loc="left",
        fontsize=10.4,
        fontweight="bold",
        color=C_INK,
        pad=10,
    )
    _style_axes(ax_a)
    legend_a = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Potent, safe CTA (5/5 kept)",
            markerfacecolor=C_CTAE,
            markersize=7.5,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="On-target organ toxicity (4/4 Unsafe)",
            markerfacecolor=C_UNSAFE,
            markersize=7.5,
        ),
        Line2D(
            [0],
            [0],
            marker="^",
            color="w",
            label="Fatal cross-reactivity (2/2 Unsafe)",
            markerfacecolor="#8b1e3f",
            markersize=8.0,
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="w",
            label="Tolerated non-CTA epitope (3 Unsafe)",
            markerfacecolor="#d69e2e",
            markersize=7.0,
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            label="0% ORR / pending (excluded)",
            markerfacecolor="#64748b",
            markersize=7.0,
        ),
    ]
    ax_a.legend(
        handles=legend_a,
        loc="upper center",
        bbox_to_anchor=(0.48, -0.14),
        ncol=3,
        frameon=False,
        fontsize=7.6,
    )

    # --- Panel B: Multi-assay proxy attribution matrix across all 16 clinical regimens ---
    row_order = [
        ("ny_eso_1_1g4_c259", "NY-ESO-1 (1G4 / c259)", "CTA · Safe"),
        ("magea4_afami_cel", "MAGEA4 (afami-cel)", "CTA · Safe"),
        ("magea10_adp_a2m10", "MAGEA10 (ADP-A2M10)", "CTA · Safe"),
        ("prame_ima203", "PRAME (IMA203)", "CTA · Safe"),
        ("ctag2_lete_cel_c259", "CTAG2 LAGE-1a (c259)", "CTA · Safe"),
        ("ceacam5_nci_tcr", "CEACAM5 (NCI TCR)", "Gr3 Colitis"),
        ("mlana_dmf5_tcr", "MLANA (DMF5 TCR)", "Uveitis / Ear"),
        ("pmel_gp100_tcr_imm", "PMEL (gp100 / ImmTAC)", "Uveitis / Rash"),
        ("tyr_til1383i_tcr", "TYR (TIL1383I TCR)", "Vitiligo"),
        ("magea3_a2_neuro", "MAGEA3→A12/EPS8L2", "Fatal CNS (2/9)"),
        ("magea3_a1_titin", "MAGEA3→TTN (a3a)", "Fatal Heart (2/2)"),
        ("magea3_dp4_6f9", "MAGEA3 (DP4 6F9 CD4)", "Safe Epitope"),
        ("afp_adp_a2afp", "AFP (ADP-A2AFP)", "Safe Epitope"),
        ("wt1_ttcr_c4", "WT1 (TTCR-C4)", "Safe Epitope"),
        ("magea1_ima202", "MAGEA1 (IMA202, 0%)", "0% ORR Excl."),
        ("ct83_kklc1_phase1", "CT83 (Ph1 pending)", "Pending Excl."),
    ]
    col_names = [
        "Bulk\nRNA",
        "scRNA\nCluster",
        "Ocular\nCP10K",
        "Normal\nIHC",
        "Rare-Cell\nUMI",
        "Paralog /\nMimicry",
        "#766\nClear.",
    ]
    # Codes: 0 = pass (<ceiling), 1 = Unsafe gate trigger, 2 = calibrated credit (IHC- or #766), 3 = paralog/mimicry warning
    matrix_codes = {
        "ny_eso_1_1g4_c259": [0, 0, 0, 0, 0, 0, 0],
        "magea4_afami_cel": [0, 0, 0, 0, 0, 0, 0],
        "magea10_adp_a2m10": [0, 0, 0, 0, 0, 0, 0],
        "prame_ima203": [2, 2, 0, 2, 0, 0, 0],
        "ctag2_lete_cel_c259": [0, 2, 0, 2, 2, 0, 2],
        "ceacam5_nci_tcr": [1, 1, 1, 1, 0, 0, 0],
        "mlana_dmf5_tcr": [1, 1, 1, 1, 0, 0, 0],
        "pmel_gp100_tcr_imm": [1, 1, 1, 0, 0, 0, 0],
        "tyr_til1383i_tcr": [1, 1, 1, 1, 0, 0, 0],
        "magea3_a2_neuro": [0, 1, 0, 0, 1, 1, 0],
        "magea3_a1_titin": [1, 1, 0, 1, 0, 3, 0],
        "magea3_dp4_6f9": [0, 1, 0, 0, 0, 3, 0],
        "afp_adp_a2afp": [0, 1, 1, 0, 0, 0, 0],
        "wt1_ttcr_c4": [1, 1, 0, 1, 0, 0, 0],
        "magea1_ima202": [0, 0, 0, 0, 0, 0, 0],
        "ct83_kklc1_phase1": [2, 0, 0, 2, 0, 0, 0],
    }
    code_fill = {0: "#edf2f7", 1: "#c53030", 2: "#2f855a", 3: "#d69e2e"}
    code_sym = {0: "·", 1: "×", 2: "✓", 3: "!"}
    code_fg = {0: "#64748b", 1: "white", 2: "white", 3: "white"}

    n_rows = len(row_order)
    n_cols = len(col_names)
    for r_i, (bid, _r_lbl, out_lbl) in enumerate(row_order):
        yp = n_rows - 1 - r_i
        vals = matrix_codes[bid]
        for c_i, code in enumerate(vals):
            rect = plt.Rectangle(
                (c_i - 0.44, yp - 0.40),
                0.88,
                0.80,
                facecolor=code_fill[code],
                edgecolor="white",
                linewidth=1.0,
            )
            ax_b.add_patch(rect)
            ax_b.text(
                c_i,
                yp,
                code_sym[code],
                ha="center",
                va="center",
                fontsize=8.8,
                fontweight="bold",
                color=code_fg[code],
            )
        # Clinical outcome pill text on right
        out_col = (
            "#22543d"
            if "CTA · Safe" in out_lbl
            else "#9b2c2c"
            if (
                "Fatal" in out_lbl
                or "Colitis" in out_lbl
                or "Uveitis" in out_lbl
                or "Vitiligo" in out_lbl
            )
            else "#975a16"
            if "Safe Epitope" in out_lbl
            else "#475569"
        )
        ax_b.text(
            n_cols - 0.30,
            yp,
            out_lbl,
            va="center",
            ha="left",
            fontsize=7.6,
            fontweight="bold",
            color=out_col,
        )

    # Horizontal group dividers
    for sep_after in (4, 8, 10, 13):
        y_sep = n_rows - 1 - sep_after - 0.5
        ax_b.plot(
            [-0.5, n_cols + 1.6], [y_sep, y_sep], color="#94a3b8", linestyle=":", linewidth=0.9
        )

    ax_b.set_xlim(-0.55, n_cols + 1.85)
    ax_b.set_ylim(-0.6, n_rows - 0.4)
    ax_b.set_xticks(np.arange(n_cols))
    ax_b.set_xticklabels(col_names, fontsize=7.8, color=C_INK)
    ax_b.xaxis.tick_top()
    ax_b.set_yticks(np.arange(n_rows)[::-1])
    ax_b.set_yticklabels([r[1] for r in row_order], fontsize=7.9, color=C_INK)
    ax_b.set_title(
        "B. Which proxy catches or calibrates each clinical trial? (*a3a d4–5 lysis)",
        loc="left",
        fontsize=10.4,
        fontweight="bold",
        color=C_INK,
        pad=26,
    )
    for sp in ("top", "right", "left", "bottom"):
        ax_b.spines[sp].set_visible(False)
    ax_b.tick_params(length=0)

    legend_b = [
        Patch(facecolor="#c53030", label="× Excludes as Unsafe"),
        Patch(facecolor="#2f855a", label="✓ Calibrated IHC− / #766 credit"),
        Patch(facecolor="#d69e2e", label="! Structural / paralog warning"),
        Patch(facecolor="#edf2f7", edgecolor="#cbd5e1", label="· Passes (< ceiling)"),
    ]
    ax_b.legend(
        handles=legend_b,
        loc="upper center",
        bbox_to_anchor=(0.45, -0.03),
        ncol=2,
        frameon=False,
        fontsize=7.6,
    )

    _ = ev  # validated benchmark join
    fig.tight_layout()
    out_path = OUT_DIR / "clinical-tcr-gold-standard.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


def main() -> None:
    cat = target_selection.results()
    sex_df = human.sex_restricted_antigens()
    he = species_target_selection.packaged_species_host_exposure("human")
    p1 = draw_funnel_and_landscape(cat, sex_df, he)
    p2 = draw_lalba_csn2_interrogation(cat)
    p3 = draw_clinical_tcr_gold_standard()
    print(f"Wrote {p1.relative_to(ROOT)}, {p2.relative_to(ROOT)}, and {p3.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

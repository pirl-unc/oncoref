#!/usr/bin/env python3
"""Render deterministic process, purity, and cross-species comparative whitepaper figures."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Patch

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from oncoref import canine, human, mouse  # noqa: E402

OUT_DIR = ROOT / "docs/antigen-whitepaper"

C_INK = "#1f2933"
C_MUTED = "#52606d"
C_GRID = "#e4e7eb"
C_UNSAFE = "#b84a47"
C_REVIEW = "#9aa7b2"
C_SEQ = "#64748b"
C_EXPR = "#cbd5e1"
C_RA = "#7b2cbf"
C_CA = "#d97706"
C_CTA = "#1d4ed8"
C_ECTOPIC = "#1e6091"
C_BOTH = "#2a9d8f"
C_LINEAGE = "#e9c46a"
C_CTAM = "#2b6cb0"
C_CTAP = "#2f855a"
C_CTAF = "#c53070"
C_CTAMF = "#c27d2b"


def _style_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#9aa5b1")
    ax.spines["bottom"].set_color("#9aa5b1")
    ax.tick_params(colors=C_INK, labelsize=8.6)


def draw_conceptual_pipeline_architecture(
    h_cat: pd.DataFrame, d_cat: pd.DataFrame, m_cat: pd.DataFrame
) -> Path:
    """Render Panel A (3-stage conceptual architecture) + Panel B (3-species proteome attrition)."""
    fig = plt.figure(figsize=(16.0, 8.0), facecolor="white")
    gs = fig.add_gridspec(2, 1, height_ratios=[1.12, 1.0], hspace=0.30)
    ax_top = fig.add_subplot(gs[0])
    ax_bot = fig.add_subplot(gs[1])

    # --- Panel A: Conceptual Process Diagram ---
    ax_top.set_xlim(0, 100)
    ax_top.set_ylim(0, 36)
    ax_top.axis("off")
    ax_top.text(
        0.0,
        34.8,
        "A. Unified Three-Stage Antigen Screening & Patient-Context Architecture (Human, Dog, Mouse)",
        fontsize=11.2,
        fontweight="bold",
        color=C_INK,
        va="top",
    )

    boxes = [
        (
            0.5,
            8.0,
            21.5,
            23.5,
            "#f8fafc",
            "#475569",
            "0. Reference Proteome & PEGs",
            [
                "• SHA-256 ReferenceIdentity",
                "  GRCh38 · CanFam3.1 · GRCm39",
                "• Collapse identical proteins",
                "  into PEG:<sha256[:20]>",
                "• Zero cross-species label",
                "  or unit transfer",
            ],
        ),
        (
            24.5,
            8.0,
            23.5,
            23.5,
            "#fff5f5",
            C_UNSAFE,
            "1. Multi-Modal Somatic Safety",
            [
                "• Somatic bulk P95 & sorted",
                "  activated T / NK / Mac RNA",
                "• Donor scRNA UMI rare-cell",
                "• Tissue IHC / MS proteomics",
                "• Vital ≤4 · Protected ≤8 TPM",
                "→ Breach = Unsafe (86–92%)",
            ],
        ),
        (
            50.5,
            8.0,
            23.5,
            23.5,
            "#f1f5f9",
            "#334155",
            "2. Eligibility, Purity & Homes",
            [
                "• ≥95% Unsafe 9-mer overlap",
                "  or readthrough → Other:Seq",
                "• Unscreened niche & OR/TAS",
                "  GPCRs → Unscreened tissue",
                "• Tumor gDNA QC (≤2.5% ORs)",
                "• Cell lines → purity sidecar",
            ],
        ),
        (
            76.5,
            8.0,
            23.0,
            23.5,
            "#eff6ff",
            C_CTA,
            "3. Symmetric 10× H₀ & Host",
            [
                "• H₀ = max(N, D, 1 TPM)",
                "• In-vivo tumor ≥10× H₀ → CA",
                "• Repro / embryo ≥10× H₀ → RA",
                "• Intersection CTA = CA ∩ RA",
                "• Tags ectopic vs lineage",
                "• HostContext → {M, F, P, MF}",
            ],
        ),
    ]

    for x, y, w, h, fc, ec, title, bullets in boxes:
        patch = FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.4,rounding_size=1.0",
            facecolor=fc,
            edgecolor=ec,
            linewidth=1.5,
        )
        ax_top.add_patch(patch)
        ax_top.text(
            x + 1.0,
            y + h - 1.8,
            title,
            fontsize=9.1,
            fontweight="bold",
            color=ec if ec != "#334155" else C_INK,
            va="top",
        )
        for bi, line in enumerate(bullets):
            ax_top.text(
                x + 1.0,
                y + h - 5.4 - bi * 3.0,
                line,
                fontsize=7.9,
                color=C_INK,
                va="top",
            )

    for x_start, x_end in ((22.3, 24.2), (48.3, 50.2), (74.3, 76.2)):
        ax_top.annotate(
            "",
            xy=(x_end, 19.8),
            xytext=(x_start, 19.8),
            arrowprops={"arrowstyle": "-|>", "color": "#475569", "lw": 1.6, "mutation_scale": 13},
        )

    # Bottom callout banner under Stage 1-3 showing exits
    exits = [
        (36.2, 2.8, "Exit: Unsafe & Other:Review", C_UNSAFE),
        (62.2, 2.8, "Exit: Other:Sequence & Unscreened tissue", "#475569"),
        (88.0, 2.8, "Output: CTA · CA only · RA only → Peptides", C_CTA),
    ]
    for cx, cy, txt, col in exits:
        ax_top.text(
            cx,
            cy,
            txt,
            ha="center",
            va="center",
            fontsize=7.7,
            fontweight="bold",
            color=col,
            bbox={
                "boxstyle": "round,pad=0.3",
                "facecolor": "white",
                "edgecolor": col,
                "lw": 0.9,
                "alpha": 0.95,
            },
        )

    # --- Panel B: Quantitative Three-Species Attrition & Candidate Yield ---
    species_specs = [
        (f"Human (GRCh38)\n{len(h_cat):,} PEGs · 8/8 pillars", h_cat),
        (f"Domestic dog (CanFam3.1)\n{len(d_cat):,} PEGs · 8/8 pillars", d_cat),
        (f"Laboratory mouse (GRCm39)\n{len(m_cat):,} PEGs · 8/8 pillars", m_cat),
    ]
    segments = [
        ("Stage 1 Unsafe", lambda df: df["category"].eq("Unsafe"), C_UNSAFE),
        ("Other:Review", lambda df: df["classification"].eq("Other:Review"), C_REVIEW),
        (
            "Other:Seq & Unscreened",
            lambda df: df["classification"].isin(["Other:Sequence", "Other:Unscreened tissue"]),
            C_SEQ,
        ),
        ("Other:Expression (<10×)", lambda df: df["classification"].eq("Other:Expression"), C_EXPR),
        (
            "RA only",
            lambda df: df["classification"].eq("Reproductive Antigen only"),
            C_RA,
        ),
        (
            "CA only",
            lambda df: df["classification"].eq("Cancer Antigen only"),
            C_CA,
        ),
        ("CTA (CA ∩ RA)", lambda df: df["classification"].eq("CTA"), C_CTA),
    ]

    y_pos = np.arange(len(species_specs))[::-1]
    left = np.zeros(len(species_specs))
    pct_mat = np.zeros((len(species_specs), len(segments)))
    cnt_mat = np.zeros((len(species_specs), len(segments)), dtype=int)

    for si, (_lbl, df) in enumerate(species_specs):
        tot = max(len(df), 1)
        for j, (_slbl, fn, _col) in enumerate(segments):
            c = int(fn(df).sum())
            cnt_mat[si, j] = c
            pct_mat[si, j] = 100.0 * c / tot

    for j, (slbl, _fn, col) in enumerate(segments):
        widths = pct_mat[:, j]
        ax_bot.barh(
            y_pos,
            widths,
            left=left,
            height=0.54,
            color=col,
            edgecolor="white",
            linewidth=0.8,
            label=slbl,
        )
        for si, (yp, w, cnt) in enumerate(zip(y_pos, widths, cnt_mat[:, j])):
            if w >= 4.2 and cnt > 0:
                fg = C_INK if col == C_EXPR else "white"
                ax_bot.text(
                    left[si] + w / 2.0,
                    yp,
                    f"{cnt:,}\n({w:.1f}%)",
                    ha="center",
                    va="center",
                    fontsize=7.6,
                    fontweight="bold",
                    color=fg,
                )
        left += widths

    for si, yp in enumerate(y_pos):
        cta_c = cnt_mat[si, 6]
        ca_c = cnt_mat[si, 5]
        ra_c = cnt_mat[si, 4]
        ax_bot.text(
            100.8,
            yp,
            f"Candidates: {cta_c:,} CTA  |  {ca_c:,} CA only  |  {ra_c:,} RA only",
            va="center",
            ha="left",
            fontsize=8.3,
            fontweight="bold",
            color=C_INK,
        )

    ax_bot.set_yticks(y_pos)
    ax_bot.set_yticklabels([s[0] for s in species_specs], fontsize=8.8, fontweight="bold")
    ax_bot.set_xlim(0, 136)
    ax_bot.set_xticks([0, 20, 40, 60, 80, 100])
    ax_bot.set_xticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])
    ax_bot.set_xlabel(
        "Share of reference protein-equivalence groups (PEGs)", fontsize=9.0, color=C_INK
    )
    ax_bot.set_title(
        "B. Empirical Proteome-Wide Attrition & Candidate Yield Across Species",
        loc="left",
        fontsize=11.0,
        fontweight="bold",
        color=C_INK,
        pad=8,
    )
    _style_axes(ax_bot)
    ax_bot.spines["bottom"].set_bounds(0, 100)
    ax_bot.legend(
        loc="upper center",
        bbox_to_anchor=(0.45, -0.22),
        ncol=7,
        frameon=False,
        fontsize=7.8,
    )

    out_path = OUT_DIR / "conceptual-pipeline-architecture.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


def draw_ectopic_vs_lineage_and_cell_lines(
    h_cat: pd.DataFrame, d_cat: pd.DataFrame, m_cat: pd.DataFrame
) -> Path:
    """Render Panel A (ectopic vs lineage), Panel B (in-vivo vs cell-line purity), Panel C (GPCRs)."""
    fig, (ax_a, ax_b, ax_c) = plt.subplots(
        1, 3, figsize=(16.5, 5.6), gridspec_kw={"width_ratios": [1.05, 1.22, 1.05]}
    )
    fig.patch.set_facecolor("white")

    # --- Panel A: Ectopic vs. Lineage-Retained CTA Activation Across Species ---
    sp_rows = [
        ("Human CTAs", h_cat[h_cat["is_cta"].astype(bool)]),
        ("Dog CTAs", d_cat[d_cat["is_cta"].astype(bool)]),
        ("Mouse CTAs", m_cat[m_cat["is_cta"].astype(bool)]),
        ("Human CA-only", h_cat[h_cat["classification"].eq("Cancer Antigen only")]),
        ("Dog CA-only", d_cat[d_cat["classification"].eq("Cancer Antigen only")]),
        ("Mouse CA-only", m_cat[m_cat["classification"].eq("Cancer Antigen only")]),
    ]
    modes = [
        ("Pure ectopic (non-homologous)", "ectopic", C_ECTOPIC),
        ("Ectopic + homologous", "ectopic_and_lineage_retained", C_BOTH),
        ("Homologous lineage-retained only", "lineage_retained", C_LINEAGE),
    ]
    ya = np.arange(len(sp_rows))[::-1]
    left_a = np.zeros(len(sp_rows))
    for mlbl, mkey, col in modes:
        pcts = []
        cnts = []
        for _rlbl, sub in sp_rows:
            tot = max(len(sub), 1)
            c = (
                int(sub["tumor_activation_mode"].eq(mkey).sum())
                if "tumor_activation_mode" in sub
                else 0
            )
            cnts.append(c)
            pcts.append(100.0 * c / tot)
        ax_a.barh(
            ya,
            pcts,
            left=left_a,
            height=0.58,
            color=col,
            edgecolor="white",
            linewidth=0.8,
            label=mlbl,
        )
        for idx_r, (yp, w, c) in enumerate(zip(ya, pcts, cnts)):
            if w >= 10.0 and c > 0:
                fg = C_INK if mkey == "lineage_retained" else "white"
                ax_a.text(
                    left_a[idx_r] + w / 2.0,
                    yp,
                    f"{c}",
                    ha="center",
                    va="center",
                    fontsize=7.8,
                    fontweight="bold",
                    color=fg,
                )
        left_a += np.array(pcts)

    for idx_r, (_rlbl, sub) in enumerate(sp_rows):
        ect_tot = int(
            sub["tumor_activation_mode"].isin(["ectopic", "ectopic_and_lineage_retained"]).sum()
        )
        ax_a.text(
            101.8,
            ya[idx_r],
            f"{ect_tot}/{len(sub)} ectopic",
            va="center",
            ha="left",
            fontsize=7.6,
            color=C_INK,
        )

    ax_a.axhline(2.5, color="#94a3b8", linestyle=":", linewidth=0.9)
    ax_a.set_yticks(ya)
    ax_a.set_yticklabels([r[0] for r in sp_rows], fontsize=8.4, fontweight="bold")
    ax_a.set_xlim(0, 134)
    ax_a.set_xticks([0, 25, 50, 75, 100])
    ax_a.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax_a.set_xlabel("Share of qualified antigens (numbers = PEGs)", fontsize=8.6, color=C_INK)
    ax_a.set_title(
        "A. Ectopic vs. homologous lineage retention",
        loc="left",
        fontsize=10.2,
        fontweight="bold",
        color=C_INK,
        pad=10,
    )
    _style_axes(ax_a)
    ax_a.spines["bottom"].set_bounds(0, 100)
    ax_a.legend(
        loc="upper center",
        bbox_to_anchor=(0.45, -0.18),
        ncol=1,
        frameon=False,
        fontsize=7.6,
    )

    # --- Panel B: In-Vivo Graft vs Pure In-Vitro Cell Line Purity (Murine CTAs vs Host TIL Granzymes) ---
    ax_b.axhspan(10.0, 8000.0, xmin=0.0, xmax=1.0, color="#f0fdf4", alpha=0.45, zorder=0)
    ax_b.axhspan(0.03, 3.0, xmin=0.28, xmax=1.0, color="#fff5f5", alpha=0.55, zorder=0)
    ax_b.axhline(10.0, color="#64748b", linestyle="--", linewidth=1.0, zorder=2)
    ax_b.axvline(10.0, color="#64748b", linestyle="--", linewidth=1.0, zorder=2)

    m_ctas = m_cat[m_cat["is_cta"].astype(bool)].copy()
    m_ctas["invivo_x"] = np.clip(m_ctas["tumor_max_tpm"].fillna(0.05), 0.05, 6000.0)
    m_ctas["cell_y"] = np.clip(m_ctas["cell_line_max_tpm"].fillna(0.05), 0.05, 6000.0)

    ect_mask = m_ctas["tumor_activation_mode"].isin(["ectopic", "ectopic_and_lineage_retained"])
    ax_b.scatter(
        m_ctas.loc[ect_mask, "invivo_x"],
        m_ctas.loc[ect_mask, "cell_y"],
        c=C_CTA,
        marker="o",
        s=52,
        alpha=0.88,
        edgecolors="white",
        linewidths=0.6,
        zorder=4,
        label=f"Ectopic murine CTAs (n={int(ect_mask.sum())})",
    )
    if (~ect_mask).any():
        ax_b.scatter(
            m_ctas.loc[~ect_mask, "invivo_x"],
            m_ctas.loc[~ect_mask, "cell_y"],
            c=C_LINEAGE,
            marker="s",
            s=52,
            alpha=0.92,
            edgecolors=C_INK,
            linewidths=0.6,
            zorder=4,
            label=f"Lineage-retained murine CTAs (n={int((~ect_mask).sum())})",
        )

    # Plot verified GSE334005 / pure-line benchmarks for excluded host LAK granzymes + Trap1a (P1A)
    host_til_pts = [
        ("Gzme (TIL NK)", 2138.0, 0.08),
        ("Gzmd (TIL NK)", 1126.0, 0.06),
        ("Gzmf (TIL NK)", 885.0, 2.26),
        ("Gzmg (TIL NK)", 443.0, 0.05),
    ]
    ax_b.scatter(
        [p[1] for p in host_til_pts],
        [p[2] for p in host_til_pts],
        c=C_UNSAFE,
        marker="^",
        s=70,
        edgecolors="white",
        linewidths=0.7,
        zorder=5,
        label="Excluded host LAK/NK granzymes (Unsafe)",
    )
    til_offsets = {
        "Gzmf (TIL NK)": (-56, 10),
        "Gzme (TIL NK)": (-26, 20),
        "Gzmd (TIL NK)": (-56, 6),
        "Gzmg (TIL NK)": (-84, -2),
    }
    for lbl_g, gx, gy in host_til_pts:
        dx_off, dy_off = til_offsets.get(lbl_g, (-48, 10))
        ax_b.annotate(
            lbl_g,
            xy=(gx, gy),
            xytext=(dx_off, dy_off),
            textcoords="offset points",
            fontsize=7.3,
            fontweight="bold",
            color=C_UNSAFE,
            arrowprops={"arrowstyle": "-", "color": "#cbd5e1", "lw": 0.6},
        )

    # Annotate all ectopic murine CTAs dynamically
    top_annot = m_ctas[ect_mask].sort_values("cell_y", ascending=False).head(6)
    used_lbls: set[str] = set()
    ect_offsets_by_lbl = {
        "Tex19.1": (-44, 8),
        "Hmga2": (8, -10),
        "Wfdc13": (-48, -11),
        "Pramel13": (8, 8),
        "Plac1": (-42, 10),
        "Serpinb9e": (8, -4),
    }
    for idx_a, r_m in enumerate(top_annot.itertuples(index=False)):
        lbl_m = str(r_m.protein_group_label)
        used_lbls.add(lbl_m)
        offs = ect_offsets_by_lbl.get(
            lbl_m, [(-44, 8), (8, -10), (-48, -11), (8, 8), (-42, 10), (8, -4)][idx_a % 6]
        )
        ax_b.annotate(
            lbl_m,
            xy=(float(r_m.invivo_x), float(r_m.cell_y)),
            xytext=offs,
            textcoords="offset points",
            fontsize=7.4,
            fontweight="bold",
            color=C_INK,
            arrowprops={"arrowstyle": "-", "color": "#94a3b8", "lw": 0.6},
        )
    lin_offsets_by_lbl = {"Pbsn": (-42, 16), "Sbp": (10, 8)}
    for lin_i, lin_r in enumerate(m_ctas[~ect_mask].head(2).itertuples(index=False)):
        lbl_lin = str(lin_r.protein_group_label)
        if lbl_lin not in used_lbls:
            ax_b.annotate(
                lbl_lin,
                xy=(float(lin_r.invivo_x), float(lin_r.cell_y)),
                xytext=lin_offsets_by_lbl.get(lbl_lin, (10, 8 if lin_i else -10)),
                textcoords="offset points",
                fontsize=7.3,
                color=C_INK,
                arrowprops={"arrowstyle": "-", "color": "#94a3b8", "lw": 0.6},
            )

    ax_b.text(
        0.06,
        4500.0,
        "Cancer-cell intrinsic\n(≥10 TPM in pure lines)",
        fontsize=7.6,
        fontweight="bold",
        color="#166534",
        va="top",
    )
    ax_b.text(
        1.15,
        2.35,
        "Host TIL / in-vivo-only\n(0–2 TPM in pure lines)",
        fontsize=7.3,
        fontweight="bold",
        color=C_UNSAFE,
        ha="left",
        va="top",
    )

    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.set_xlim(0.035, 7500.0)
    ax_b.set_ylim(0.035, 7500.0)
    ticks_b = [0.05, 1, 10, 100, 1000]
    lbls_b = ["≤0.05", "1", "10", "100", "1k"]
    ax_b.set_xticks(ticks_b)
    ax_b.set_xticklabels(lbls_b)
    ax_b.set_yticks(ticks_b)
    ax_b.set_yticklabels(lbls_b)
    ax_b.grid(True, which="major", color=C_GRID, linewidth=0.6)
    ax_b.set_xlabel("Peak in-vivo syngeneic / GEMM / carcinogen TPM", fontsize=8.6, color=C_INK)
    ax_b.set_ylabel("Peak pure in-vitro cancer cell-line TPM (9 lines)", fontsize=8.6, color=C_INK)
    ax_b.set_title(
        "B. In-vivo vs. pure cell-line expression separates TILs",
        loc="left",
        fontsize=10.2,
        fontweight="bold",
        color=C_INK,
        pad=10,
    )
    _style_axes(ax_b)
    ax_b.legend(
        loc="upper center",
        bbox_to_anchor=(0.48, -0.18),
        ncol=1,
        frameon=False,
        fontsize=7.5,
    )

    # --- Panel C: Monogenic Chemosensory GPCR Harmonization Across All 3 Species ---
    raw_gpcr_specs = [
        (
            "Human",
            h_cat,
            re.compile(r"^(?:OR\d+[A-Z]+\d+P?|TAAR\d+P?|TAS1R\d+|TAS2R\d+|VN1R\d+)$"),
        ),
        (
            "Dog",
            d_cat,
            re.compile(
                r"^(?:COR\d+|OR\d+[A-Z]*\d*|TAAR\d+|TAS1R\d+|TAS2R\d+|VN1R\d+)$", re.IGNORECASE
            ),
        ),
        (
            "Mouse",
            m_cat,
            re.compile(
                r"^(?:Olfr\d+|Or\d+[a-z]+\d+|Taar\d+[a-z]*|Tas1r\d+|Tas2r\d+|Vmn1r\d+|Vmn2r\d+)$",
                re.IGNORECASE,
            ),
        ),
    ]
    gpcr_specs = []
    for sp_name, df_sp, rx_sp in raw_gpcr_specs:
        m_sp = (
            df_sp["protein_group_label"]
            .astype(str)
            .map(lambda x, r=rx_sp: bool(r.match(x.split("/")[0])))
        )
        gpcr_specs.append((f"{sp_name}\n({int(m_sp.sum()):,} GPCRs)", df_sp.loc[m_sp], rx_sp))
    gpcr_cats = [
        ("Unsafe (somatic)", lambda s: s.eq("Unsafe"), C_UNSAFE),
        ("Other:Review / Seq", lambda s: s.isin(["Other:Review", "Other:Sequence"]), C_REVIEW),
        ("Other:Expression (<10×)", lambda s: s.eq("Other:Expression"), C_EXPR),
        ("Other:Unscreened tissue", lambda s: s.eq("Other:Unscreened tissue"), "#0f766e"),
        ("RA only (germline ≥10×)", lambda s: s.eq("Reproductive Antigen only"), C_RA),
        ("CA / CTA (admitted)", lambda s: s.isin(["Cancer Antigen only", "CTA"]), C_CTA),
    ]
    yc = np.arange(len(gpcr_specs))[::-1]
    left_c = np.zeros(len(gpcr_specs))
    for glbl, gfn, col in gpcr_cats:
        pcts_c = []
        cnts_c = []
        for _slbl, sub_df, _rx in gpcr_specs:
            sub = sub_df["classification"]
            tot = max(len(sub), 1)
            c = int(gfn(sub).sum())
            cnts_c.append(c)
            pcts_c.append(100.0 * c / tot)
        ax_c.barh(
            yc,
            pcts_c,
            left=left_c,
            height=0.52,
            color=col,
            edgecolor="white",
            linewidth=0.8,
            label=glbl,
        )
        for idx_r, (yp, w, c) in enumerate(zip(yc, pcts_c, cnts_c)):
            if w >= 9.0 and c > 0:
                fg = C_INK if col == C_EXPR else "white"
                ax_c.text(
                    left_c[idx_r] + w / 2.0,
                    yp,
                    f"{c}",
                    ha="center",
                    va="center",
                    fontsize=7.8,
                    fontweight="bold",
                    color=fg,
                )
        left_c += np.array(pcts_c)

    for _idx_r, yp in enumerate(yc):
        ax_c.text(
            101.8,
            yp,
            "0 CA / 0 CTA",
            va="center",
            ha="left",
            fontsize=8.0,
            fontweight="bold",
            color="#0f766e",
        )

    ax_c.set_yticks(yc)
    ax_c.set_yticklabels([g[0] for g in gpcr_specs], fontsize=8.5, fontweight="bold")
    ax_c.set_xlim(0, 132)
    ax_c.set_xticks([0, 25, 50, 75, 100])
    ax_c.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax_c.set_xlabel("Share of chemosensory GPCR PEGs", fontsize=8.6, color=C_INK)
    ax_c.set_title(
        "C. Monogenic chemosensory GPCR guard",
        loc="left",
        fontsize=10.2,
        fontweight="bold",
        color=C_INK,
        pad=10,
    )
    _style_axes(ax_c)
    ax_c.spines["bottom"].set_bounds(0, 100)
    ax_c.legend(
        loc="upper center",
        bbox_to_anchor=(0.45, -0.18),
        ncol=2,
        frameon=False,
        fontsize=7.3,
    )

    fig.tight_layout()
    out_path = OUT_DIR / "ectopic-vs-lineage-and-cell-lines.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


def draw_comparative_species_funnel_and_cohorts(
    h_cat: pd.DataFrame, d_cat: pd.DataFrame, m_cat: pd.DataFrame
) -> Path:
    """Render Panel A (dog cohorts), Panel B (mouse cohorts), Panel C (mammalian orthology & sex classes)."""
    fig, (ax_a, ax_b, ax_c) = plt.subplots(
        1, 3, figsize=(16.6, 6.4), gridspec_kw={"width_ratios": [1.18, 1.15, 1.05]}
    )
    fig.patch.set_facecolor("white")

    # --- Panel A: Spontaneous Companion-Dog Patient Cohorts ---
    d_coh = canine.tumor_cohorts().copy()
    d_ctas = d_cat[d_cat["is_cta"].astype(bool)]
    peak_cnt_dog = d_ctas["tumor_top_cancer_code"].value_counts().to_dict()
    d_coh["peak_ctas"] = d_coh["cancer_code"].map(lambda c: int(peak_cnt_dog.get(c, 0)))
    d_coh = d_coh.sort_values(["n_samples", "cancer_code"], ascending=[True, False]).reset_index(
        drop=True
    )

    ya = np.arange(len(d_coh))
    ax_a.barh(
        ya,
        d_coh["n_samples"],
        height=0.66,
        color="#2b6cb0",
        alpha=0.86,
        edgecolor="white",
    )
    max_dn = max(int(d_coh["n_samples"].max()), 10)
    for yi, r in enumerate(d_coh.itertuples(index=False)):
        pk_txt = (
            f" · {r.peak_ctas} peak CTA" + ("s" if r.peak_ctas != 1 else "")
            if r.peak_ctas > 0
            else ""
        )
        ax_a.text(
            r.n_samples + max_dn * 0.02,
            yi,
            f"n={r.n_samples}{pk_txt}",
            va="center",
            ha="left",
            fontsize=7.2,
            color=C_INK,
        )

    ax_a.set_yticks(ya)
    ax_a.set_yticklabels(d_coh["cancer_code"], fontsize=7.8, fontweight="bold")
    ax_a.set_xlim(0, max_dn * 1.48)
    ax_a.set_xlabel(
        f"QC-passing spontaneous canine patient tumors (total n = {int(d_coh['n_samples'].sum()):,})",
        fontsize=8.5,
        color=C_INK,
    )
    ax_a.set_title(
        f"A. Canine spontaneous cohorts ({len(d_coh)} cohorts + {len(canine.cell_lines())} lines)",
        loc="left",
        fontsize=10.0,
        fontweight="bold",
        color=C_INK,
        pad=10,
    )
    ax_a.grid(True, axis="x", color=C_GRID, linewidth=0.6)
    _style_axes(ax_a)

    # --- Panel B: Murine In-Vivo Syngeneic / GEMM / Carcinogen Cohorts ---
    m_coh = mouse.tumor_cohorts().copy()
    m_ctas = m_cat[m_cat["is_cta"].astype(bool)]
    peak_cnt_mouse = m_ctas["tumor_top_cancer_code"].value_counts().to_dict()
    m_coh["peak_ctas"] = m_coh["cancer_code"].map(lambda c: int(peak_cnt_mouse.get(c, 0)))
    m_coh = m_coh.sort_values(["n_samples", "cancer_code"], ascending=[True, False]).reset_index(
        drop=True
    )
    kind_colors = {
        "syngeneic_graft": "#0f766e",
        "primary_gemm": "#7b2cbf",
        "carcinogen_primary": "#d97706",
    }
    yb = np.arange(len(m_coh))
    ax_b.barh(
        yb,
        m_coh["n_samples"],
        height=0.64,
        color=[kind_colors.get(str(k), "#0f766e") for k in m_coh["model_kind"]],
        alpha=0.88,
        edgecolor="white",
    )
    max_mn = max(int(m_coh["n_samples"].max()), 5)
    for yi, r in enumerate(m_coh.itertuples(index=False)):
        pk_txt = f" · {r.peak_ctas} peak" if r.peak_ctas > 0 else ""
        ax_b.text(
            r.n_samples + max_mn * 0.025,
            yi,
            f"n={r.n_samples} ({r.mhc_haplotype}){pk_txt}",
            va="center",
            ha="left",
            fontsize=7.3,
            color=C_INK,
        )

    ax_b.set_yticks(yb)
    ax_b.set_yticklabels(m_coh["cancer_code"], fontsize=7.9, fontweight="bold")
    ax_b.set_xlim(0, max_mn * 1.65)
    ax_b.set_xlabel(
        f"QC-passing in-vivo murine tumors (total n = {int(m_coh['n_samples'].sum())}) + {len(mouse.cell_lines())} cell lines",
        fontsize=8.5,
        color=C_INK,
    )
    ax_b.set_title(
        f"B. Murine in-vivo models ({len(m_coh)} cohorts, 3 MHC backgrounds)",
        loc="left",
        fontsize=10.0,
        fontweight="bold",
        color=C_INK,
        pad=10,
    )
    ax_b.grid(True, axis="x", color=C_GRID, linewidth=0.6)
    _style_axes(ax_b)
    ax_b.legend(
        handles=[
            Patch(facecolor="#0f766e", label="Syngeneic / orthotopic graft"),
            Patch(facecolor="#7b2cbf", label="Primary autochthonous GEMM"),
            Patch(facecolor="#d97706", label="Carcinogen primary (DEN)"),
        ],
        loc="lower right",
        frameon=True,
        facecolor="white",
        edgecolor=C_GRID,
        fontsize=7.2,
    )

    # --- Panel C: 4-Way Physiological Sex-Restriction Partition Across Species ---
    sex_specs = [
        (f"Human\n({int(h_cat['is_cta'].sum())} CTAs)", h_cat[h_cat["is_cta"].astype(bool)]),
        (f"Dog\n({int(d_cat['is_cta'].sum())} CTAs)", d_cat[d_cat["is_cta"].astype(bool)]),
        (f"Mouse\n({int(m_cat['is_cta'].sum())} CTAs)", m_cat[m_cat["is_cta"].astype(bool)]),
    ]
    sex_classes = [
        ("CTA-M (male-restricted)", "CTA-M", C_CTAM),
        ("CTA-P (pregnancy/peripartum)", "CTA-P", C_CTAP),
        ("CTA-F (female-restricted)", "CTA-F", C_CTAF),
        ("CTA-MF (shared M & F)", "CTA-MF", C_CTAMF),
    ]
    yc = np.arange(len(sex_specs))[::-1]
    left_c = np.zeros(len(sex_specs))
    for slbl, scode, col in sex_classes:
        pcts = []
        cnts = []
        for _sp_lbl, sub in sex_specs:
            tot = max(len(sub), 1)
            c = int(sub["sex_restriction_code"].eq(scode).sum())
            cnts.append(c)
            pcts.append(100.0 * c / tot)
        ax_c.barh(
            yc,
            pcts,
            left=left_c,
            height=0.50,
            color=col,
            edgecolor="white",
            linewidth=0.8,
            label=slbl,
        )
        for idx_r, (yp, w, c) in enumerate(zip(yc, pcts, cnts)):
            if w >= 8.0 and c > 0:
                ax_c.text(
                    left_c[idx_r] + w / 2.0,
                    yp,
                    f"{c}\n({w:.0f}%)",
                    ha="center",
                    va="center",
                    fontsize=7.5,
                    fontweight="bold",
                    color="white",
                )
        left_c += np.array(pcts)

    # Also annotate conserved multi-species CTA count from cross-species-orthologs if present
    orth_csv = OUT_DIR / "cross-species-orthologs.csv"
    if orth_csv.is_file():
        odf = pd.read_csv(orth_csv)
        n_cons = int(odf["concordance_regime"].eq("Conserved multi-species CTA").sum())
        n_ra_cons = int(odf["concordance_regime"].eq("CTA + conserved RA-only ortholog").sum())
        ax_c.text(
            50.0,
            -0.68,
            f"Orthology across CTA union: {n_cons} multi-species CTAs +\n{n_ra_cons} CTA + conserved RA-only families",
            ha="center",
            va="top",
            fontsize=7.8,
            fontweight="bold",
            color=C_INK,
            bbox={
                "boxstyle": "round,pad=0.3",
                "facecolor": "#f8fafc",
                "edgecolor": "#cbd5e1",
                "lw": 0.9,
            },
        )

    ax_c.set_yticks(yc)
    ax_c.set_yticklabels([s[0] for s in sex_specs], fontsize=8.6, fontweight="bold")
    ax_c.set_xlim(0, 102)
    ax_c.set_ylim(-1.15, 2.6)
    ax_c.set_xticks([0, 25, 50, 75, 100])
    ax_c.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax_c.set_xlabel("Share of species CTAs by HostContext physiology", fontsize=8.5, color=C_INK)
    ax_c.set_title(
        "C. 4-way sex/peripartum partition across species",
        loc="left",
        fontsize=10.0,
        fontweight="bold",
        color=C_INK,
        pad=10,
    )
    _style_axes(ax_c)
    ax_c.spines["bottom"].set_bounds(0, 100)
    ax_c.legend(
        loc="upper center",
        bbox_to_anchor=(0.48, -0.18),
        ncol=2,
        frameon=False,
        fontsize=7.4,
    )

    fig.tight_layout()
    out_path = OUT_DIR / "comparative-species-funnel-and-cohorts.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


def main() -> None:
    h_cat = human.target_selection()
    d_cat = canine.target_selection()
    m_cat = mouse.target_selection()
    p1 = draw_conceptual_pipeline_architecture(h_cat, d_cat, m_cat)
    p2 = draw_ectopic_vs_lineage_and_cell_lines(h_cat, d_cat, m_cat)
    p3 = draw_comparative_species_funnel_and_cohorts(h_cat, d_cat, m_cat)
    print(f"Wrote:\n  {p1}\n  {p2}\n  {p3}")


if __name__ == "__main__":
    main()

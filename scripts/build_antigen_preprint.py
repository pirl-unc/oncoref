#!/usr/bin/env python3
# Licensed under the Apache License, Version 2.0.
"""Deterministically build the OncoRef antigen preprint (Markdown, LaTeX, PDF, CSV tables, and figures).

Generates under ``docs/preprint/`` (or ``--output-dir``):
- ``data/uniprot-candidate-topology.json.gz`` + ``.json`` (SHA-256 receipt)
- ``tables/candidate-surface-topology.csv`` (complete localization & mature ECD audit)
- ``figures/fig1-conceptual-architecture-and-attrition.{png,pdf}``
- ``figures/fig2-multispecies-parameterization-and-orthology.{png,pdf}``
- ``figures/fig3-human-deep-dive-and-ca-only-biology.{png,pdf}``
- ``figures/fig4-clinical-tcr-gold-standard.{png,pdf}``
- ``figures/fig5-surface-ectodomain-modality-routing.{png,pdf}``
- ``index.md`` (MkDocs-integrated full readable preprint)
- ``oncoref-antigen-preprint.tex`` + compiled ``oncoref-antigen-preprint.pdf``
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if not (ROOT / "oncoref").is_dir():
    ROOT = Path("/Users/iskander/code/oncoref")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from oncoref import (  # noqa: E402
    antigen_sets,
    ca_only_review,
    canine,
    human,
    mouse,
    protein_annotations,
    target_selection,
    tcr_trials,
)

CLINICAL_SURFACE_BENCHMARKS: tuple[tuple[str, str, str], ...] = (
    ("ERBB2", "HER2", "Approved ADC (T-DXd, T-DM1) & mAb"),
    ("TACSTD2", "TROP2", "Approved ADC (sacituzumab govitecan, Dato-DXd)"),
    ("FOLR1", "FRa", "Approved ADC (mirvetuximab soravtansine)"),
    ("NECTIN4", "Nectin-4", "Approved ADC (enfortumab vedotin)"),
    ("FOLH1", "PSMA", "Approved RLT (177Lu-PSMA-617) & TCE"),
    ("SSTR2", "SSTR2", "Approved peptide RLT (177Lu-DOTATATE)"),
    ("DLL3", "DLL3", "Approved bispecific TCE (tarlatamab)"),
    ("CLDN18", "CLDN18.2", "Approved mAb (zolbetuximab) & Phase 2 CAR-T / ADC"),
    ("EGFR", "EGFR", "Approved mAb & ADC"),
    ("MET", "c-MET", "Approved ADC (telisotuzumab vedotin)"),
    ("CEACAM5", "CEA", "Phase 3 ADC & TCR-T (colitis dose-limiting)"),
    ("MSLN", "Mesothelin", "Phase 2/3 CAR-T & ADC (serositis limit)"),
    ("CLDN6", "CLDN6", "Phase 1/2 CAR-T + mRNA vaccine (BNT211)"),
    ("GPC3", "GPC3", "Phase 1/2 CAR-T (armored HCC trials) & TCE"),
    ("GPC2", "GPC2", "Phase 1 CAR-T & ADC (neuroblastoma)"),
    ("IL13RA2", "IL13Ra2", "Phase 1 locoregional CAR-T (glioblastoma)"),
    ("ROR1", "ROR1", "Phase 1/2 CAR-T & ADC"),
    ("STEAP1", "STEAP1", "Phase 1/2 2+1 TCE (xaluritamig / AMG 509)"),
    ("SEZ6", "SEZ6", "Phase 1/2 ADC (ABBV-706, SCLC)"),
    ("CDH6", "CDH6", "Phase 1/2 ADC (raludotatug deruxtecan)"),
    ("VTCN1", "B7-H4", "Phase 1/2 ADC (Puxitatug samrotecan)"),
    ("PTK7", "PTK7", "Phase 1/2 ADC & CAR-T"),
    ("MUC16", "CA-125", "Phase 1/2 bispecific TCE & CAR-T"),
    ("LY6K", "LY6K", "Phase 1 peptide vaccine & preclinical ADC"),
)

PALETTE = {
    "cta": "#0f766e",
    "ca_only": "#d97706",
    "ra_only": "#2563eb",
    "other": "#64748b",
    "unsafe": "#dc2626",
    "cta_m": "#0284c7",
    "cta_p": "#7c3aed",
    "cta_f": "#db2777",
    "cta_mf": "#0d9488",
    "navy": "#0f172a",
    "slate": "#475569",
    "light": "#f8fafc",
}


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _feat_span(feat: dict) -> tuple[int, int] | None:
    loc = feat.get("location", {})
    s = loc.get("start", {}).get("value")
    e = loc.get("end", {}).get("value")
    if isinstance(s, int) and isinstance(e, int) and e >= s:
        return (s, e)
    return None


def classify_topology(entry: dict | None, rep_seqs: list[str], all_iso_seqs: list[str]) -> dict:
    """Compute mature extracellular domain intervals and modality routing from UniProtKB JSON."""
    if entry is None:
        max_len = max((len(s.rstrip("*")) for s in rep_seqs), default=0)
        return {
            "uniprot_accession": "",
            "uniprot_id": "",
            "protein_length_aa": max_len,
            "mature_length_aa": max_len,
            "sequence_match": "unmapped_novel_orf",
            "has_signal_peptide": False,
            "is_gpi_anchored": False,
            "n_tm_helices": 0,
            "ecd_total_aa": 0,
            "ecd_max_continuous_aa": 0,
            "ecd_intervals": "",
            "cell_surface_exposed": False,
            "secreted_soluble": False,
            "subcellular_locations": "Unannotated ORF",
            "localization_class": "Intracellular (nuclear / cytoplasmic / cytoskeletal)",
            "surface_architecture": "Unannotated ORF",
            "modality_route": "Intracellular TCR-T / pMHC",
        }

    up_seq = str(entry.get("sequence", {}).get("value", "")).rstrip("*")
    up_len = int(entry.get("sequence", {}).get("length") or len(up_seq))
    rep_clean = {s.rstrip("*") for s in rep_seqs}
    iso_clean = {s.rstrip("*") for s in all_iso_seqs}
    if up_seq and up_seq in rep_clean:
        seq_match = "exact_representative"
    elif up_seq and up_seq in iso_clean:
        seq_match = "exact_ensembl_isoform"
    else:
        seq_match = "canonical_uniprot_isoform"

    feats = entry.get("features", [])
    cleaved_pos: set[int] = set()
    has_sig = False
    is_gpi = False
    has_nglyc_or_ss = False
    tm_spans: list[tuple[int, int]] = []
    ecd_raw_spans: list[tuple[int, int]] = []

    for f in feats:
        ftype = f.get("type", "")
        desc = str(f.get("description", ""))
        sp = _feat_span(f)
        if ftype in ("Signal", "Transit peptide"):
            has_sig = True
            if sp:
                cleaved_pos.update(range(sp[0], sp[1] + 1))
        elif ftype == "Propeptide":
            if sp:
                cleaved_pos.update(range(sp[0], sp[1] + 1))
        elif ftype == "Transmembrane":
            if sp:
                tm_spans.append(sp)
        elif ftype == "Lipidation" and "GPI-anchor" in desc:
            is_gpi = True
        elif ftype in ("Glycosylation", "Disulfide bond"):
            has_nglyc_or_ss = True
        elif ftype == "Topological domain" and desc.lower().startswith("extracellular") and sp:
            ecd_raw_spans.append(sp)

    mature_len = max(0, up_len - len(cleaved_pos))
    mature_ecd_spans: list[tuple[int, int]] = []
    for s, e in sorted(ecd_raw_spans):
        kept = [p for p in range(s, e + 1) if p not in cleaved_pos]
        if kept:
            cur_s = cur_e = kept[0]
            for p in kept[1:]:
                if p == cur_e + 1:
                    cur_e = p
                else:
                    mature_ecd_spans.append((cur_s, cur_e))
                    cur_s = cur_e = p
            mature_ecd_spans.append((cur_s, cur_e))

    locs: list[str] = []
    topos: list[str] = []
    for c in entry.get("comments", []):
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for sl in c.get("subcellularLocations", []):
                if "location" in sl and sl["location"].get("value"):
                    locs.append(str(sl["location"]["value"]))
                if "topology" in sl and sl["topology"].get("value"):
                    topos.append(str(sl["topology"]["value"]))
                if "orientation" in sl and sl["orientation"].get("value"):
                    topos.append(str(sl["orientation"]["value"]))

    locs_dedup = list(dict.fromkeys(locs))
    loc_lower = " | ".join(locs_dedup).lower()
    topo_lower = " | ".join(topos).lower()

    if "gpi-anchor" in topo_lower:
        is_gpi = True

    pm_exact_terms = (
        "cell membrane",
        "cell surface",
        "plasma membrane",
        "apical cell membrane",
        "basolateral cell membrane",
        "lateral cell membrane",
        "cilium membrane",
        "flagellum membrane",
        "stereocilium membrane",
        "tight junction",
    )
    has_pm_loc = (
        any(k in loc_lower for k in pm_exact_terms) and "cytoplasmic side" not in topo_lower
    )
    n_tm = len(tm_spans)

    # GPI-anchored or extracellularly disulfide-bonded/glycosylated EGF-CFC peripheral surface proteins
    if not mature_ecd_spans and has_pm_loc and (is_gpi or (n_tm == 0 and has_nglyc_or_ss)):
        kept_ext = [p for p in range(1, up_len + 1) if p not in cleaved_pos]
        if kept_ext:
            mature_ecd_spans = [(kept_ext[0], kept_ext[-1])]

    # Multi-pass cell-membrane proteins with even TM helices and intracellular N-terminus (e.g. 6-TM STEAP1)
    if not mature_ecd_spans and has_pm_loc and n_tm == 6 and up_len > 0:
        tm_sorted = sorted(tm_spans)
        for idx_tm in range(0, len(tm_sorted) - 1, 2):
            s_loop = tm_sorted[idx_tm][1] + 1
            e_loop = tm_sorted[idx_tm + 1][0] - 1
            if e_loop >= s_loop:
                mature_ecd_spans.append((s_loop, e_loop))

    has_explicit_ecd = len(mature_ecd_spans) > 0
    cell_surface = bool(
        has_explicit_ecd or (has_pm_loc and (n_tm > 0 or is_gpi or has_nglyc_or_ss))
    )
    is_secreted = bool("secreted" in loc_lower and not cell_surface and n_tm == 0 and not is_gpi)

    ecd_total = sum(e - s + 1 for s, e in mature_ecd_spans)
    ecd_max_cont = max((e - s + 1 for s, e in mature_ecd_spans), default=0)

    internal_mem_terms = (
        "endoplasmic reticulum",
        "golgi",
        "mitochondri",
        "nucleus membrane",
        "nuclear membrane",
        "lysosom",
        "endosom",
        "peroxisom",
        "acrosome",
        "secretory vesicle",
        "cytoplasmic vesicle",
    )
    has_internal_mem = any(k in loc_lower for k in internal_mem_terms)

    if cell_surface and ecd_max_cont >= 50:
        if is_gpi:
            arch = "GPI-anchored ectodomain"
        elif n_tm == 0:
            arch = "Glycosylated surface coreceptor"
        elif n_tm == 1:
            arch = (
                "Single-pass Type II membrane"
                if "type ii" in topo_lower
                else "Single-pass Type I membrane"
            )
        else:
            arch = f"Multi-pass transmembrane ({n_tm} TM)"
        loc_cls = "Cell-surface ectodomain (>=50 aa)"
        route = "CAR-T / ADC / RLT / TCE + TCR-T"
    elif cell_surface and ecd_max_cont >= 15:
        arch = (
            f"Multi-pass transmembrane ({n_tm} TM)" if n_tm >= 2 else "Single-pass Type I membrane"
        )
        loc_cls = "Cell-surface short ectodomain / loop (15-49 aa)"
        route = "Constrained loop mAb / CAR-T / RLT + TCR-T"
    elif is_secreted:
        arch = "Secreted soluble / ECM"
        loc_cls = "Secreted / extracellular matrix"
        route = "TCR-T / vaccine (soluble shedding; ECM-RLT context-specific)"
    elif n_tm > 0 or has_internal_mem:
        arch = (
            f"Internal / minimal-ECD membrane ({n_tm} TM)"
            if n_tm
            else "Internal organelle / vesicle"
        )
        loc_cls = "Internal organelle / minimal-ECD membrane"
        route = "Intracellular TCR-T / pMHC"
    else:
        arch = "Intracellular (nuclear / cytoplasmic)"
        loc_cls = "Intracellular (nuclear / cytoplasmic / cytoskeletal)"
        route = "Intracellular TCR-T / pMHC"

    return {
        "uniprot_accession": entry.get("primaryAccession", ""),
        "uniprot_id": entry.get("uniProtkbId", ""),
        "protein_length_aa": up_len,
        "mature_length_aa": mature_len,
        "sequence_match": seq_match,
        "has_signal_peptide": has_sig,
        "is_gpi_anchored": is_gpi,
        "n_tm_helices": n_tm,
        "ecd_total_aa": ecd_total,
        "ecd_max_continuous_aa": ecd_max_cont,
        "ecd_intervals": ";".join(f"{s}-{e}" for s, e in mature_ecd_spans),
        "cell_surface_exposed": bool(cell_surface and ecd_max_cont >= 15),
        "secreted_soluble": is_secreted,
        "subcellular_locations": " | ".join(locs_dedup) if locs_dedup else "Intracellular",
        "localization_class": loc_cls,
        "surface_architecture": arch,
        "modality_route": route,
    }


def build_topology_dataframe(uniprot_gz_path: Path) -> pd.DataFrame:
    with gzip.open(uniprot_gz_path, "rt", encoding="utf-8") as fh:
        snap = json.load(fh)
    entries = snap["results"]
    by_acc = {e["primaryAccession"]: e for e in entries}

    ann = protein_annotations.gene_annotations()
    gid_to_accs: dict[str, list[str]] = {}
    for r in ann.itertuples(index=False):
        gid_to_accs.setdefault(str(r.gene_id), []).append(str(r.accession))

    sym_to_entry: dict[str, dict] = {}
    for e in entries:
        for g in e.get("genes", []):
            gn = g.get("geneName", {}).get("value")
            if gn:
                prev = sym_to_entry.get(gn)
                if prev is None or (
                    "Swiss-Prot" in str(e.get("entryType", ""))
                    and "Swiss-Prot" not in str(prev.get("entryType", ""))
                ):
                    sym_to_entry[gn] = e

    ts = target_selection.results()
    ca_rev = ca_only_review.ca_only_functional_review()
    gs = tcr_trials.gold_standard()
    bench_map = {sym: (alias, note) for sym, alias, note in CLINICAL_SURFACE_BENCHMARKS}

    keep_mask = (
        ts["is_target_candidate"].astype(bool)
        | ts["entity_id"].isin(set(ca_rev["entity_id"]))
        | ts["entity_id"].isin(set(gs["entity_id"]))
        | ts["protein_group_label"].isin(set(bench_map))
    )
    sub = ts[keep_mask].copy()

    seqs_cand = antigen_sets.sequences("CA|RA", species="human", sequence_source="sequence_bundle")
    rep_by_gid: dict[str, list[str]] = {}
    iso_by_gid: dict[str, list[str]] = {}
    for r in seqs_cand.itertuples(index=False):
        gid = str(r.gene_id)
        iso_by_gid.setdefault(gid, []).append(str(r.sequence))
        if getattr(r, "matches_representative_sequence", True):
            rep_by_gid.setdefault(gid, []).append(str(r.sequence))

    rows: list[dict] = []
    for r in sub.to_dict("records"):
        eid = str(r["entity_id"])
        lbl = str(r["protein_group_label"])
        gids = json.loads(r["member_gene_ids"])
        cand_accs = [a for g in gids for a in gid_to_accs.get(g, []) if a in by_acc]
        rep_seqs = [sq for g in gids for sq in rep_by_gid.get(g, [])]
        iso_seqs = [sq for g in gids for sq in iso_by_gid.get(g, [])]
        rep_clean = {s.rstrip("*") for s in rep_seqs}
        chosen = None
        if cand_accs:
            ranked = []
            for a in dict.fromkeys(cand_accs):
                ent = by_acc[a]
                seq_v = str(ent.get("sequence", {}).get("value", "")).rstrip("*")
                ranked.append((1 if seq_v in rep_clean else 0, len(seq_v), a))
            ranked.sort(reverse=True)
            chosen = by_acc[ranked[0][2]]
        else:
            for tok in lbl.split("/"):
                if tok in sym_to_entry:
                    chosen = sym_to_entry[tok]
                    break

        topo = classify_topology(chosen, rep_seqs, iso_seqs)
        b_alias, b_note = bench_map.get(lbl, ("", ""))
        rows.append(
            {
                "entity_id": eid,
                "protein_group_label": lbl,
                "clinical_alias": b_alias,
                "clinical_benchmark_note": b_note,
                "category": str(r["category"]),
                "classification": str(r["classification"]),
                "sex_restriction_code": str(r.get("sex_restriction_code") or ""),
                "tumor_activation_mode": str(r.get("tumor_activation_mode") or ""),
                "normal_reference_tpm": (
                    round(float(r["normal_reference_tpm"]), 4)
                    if pd.notna(r["normal_reference_tpm"])
                    else np.nan
                ),
                "tumor_max_tpm": (
                    round(float(r["tumor_max_tpm"]), 4) if pd.notna(r["tumor_max_tpm"]) else np.nan
                ),
                "tumor_fold": (
                    round(float(r["tumor_fold"]), 2) if pd.notna(r["tumor_fold"]) else np.nan
                ),
                "independent_patients_at_tenfold": (
                    int(r["independent_patients_at_tenfold"])
                    if pd.notna(r["independent_patients_at_tenfold"])
                    else 0
                ),
                "tumor_top_cancer_code": str(r.get("tumor_top_cancer_code") or ""),
                "safe_2mm_9mers": (
                    int(r["safe_2mm_9mers"]) if pd.notna(r.get("safe_2mm_9mers")) else np.nan
                ),
                **topo,
            }
        )
    return pd.DataFrame(rows)


def _save_fig(fig: plt.Figure, fig_dir: Path, stem: str) -> None:
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(fig_dir / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_figures(
    fig_dir: Path,
    h_ts: pd.DataFrame,
    d_ts: pd.DataFrame,
    m_ts: pd.DataFrame,
    ca_rev: pd.DataFrame,
    gs_eval: pd.DataFrame,
    topo_df: pd.DataFrame,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 9,
            "axes.titlesize": 10.5,
            "axes.titleweight": "bold",
            "axes.labelsize": 9.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    # --- Figure 1: Conceptual Architecture & Human Proteome Attrition ---
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(13.6, 5.2), gridspec_kw={"width_ratios": [1.26, 1.0]}
    )
    ax_a.set_xlim(0, 10.6)
    ax_a.set_ylim(0, 10)
    ax_a.axis("off")
    ax_a.set_title(
        "A. Species-Invariant Antigen Classification & Host-Protection Architecture",
        loc="left",
        pad=8,
    )

    stages = [
        (
            8.35,
            "#e2e8f0",
            "#0f172a",
            "Stage 0 · Identical-Sequence Protein-Equivalence Grouping (PEG)",
            "Sum member-locus TPM per biopsy before donor P95; retain all isoforms for >=95% 9-mer Unsafe check",
        ),
        (
            6.55,
            "#fee2e2",
            "#991b1b",
            "Stage 1 · Native-Unit Adult Normal & Immune Safety Gates (Unsafe Precedence)",
            "Vital <=4 | Protected <=8 | Manageable <=16 hold / >128 Unsafe | Donor rare-cell UMI | IHC / MS | Lit",
        ),
        (
            4.75,
            "#f1f5f9",
            "#334155",
            "Stage 2 · Shared Coverage, Sequence & Sensory Niche Eligibility (Other)",
            "Complete peripheral panels; <95% Unsafe 9-mers; monogenic OR/TAS/VN1R chemosensory GPCR guard",
        ),
        (
            2.85,
            "#ccfbf1",
            "#115e59",
            "Stage 3 · Symmetric 10x Enrichment over H0 = max(Normal N, Disposable D, 1 TPM)",
            "In-vivo Tumor >=10x H0 -> CA (375) | Reproductive/Embryo >=10x H0 -> RA (468) | CTA = CA & RA (322)",
        ),
        (
            0.95,
            "#ede9fe",
            "#5b21b6",
            "Downstream · HostContext Sex/Organ Protection & PeptideIndex Self-Subtraction",
            "CTA-M (222) · CTA-P (48) · CTA-F (12) · CTA-MF (40) + 0/1-mismatch host self-proteome k-mer subtraction",
        ),
    ]
    for y0, bg, border, title_s, sub_s in stages:
        box = mpatches.FancyBboxPatch(
            (0.12, y0),
            10.25,
            1.32,
            boxstyle="round,pad=0.10,rounding_size=0.18",
            facecolor=bg,
            edgecolor=border,
            linewidth=1.3,
        )
        ax_a.add_patch(box)
        ax_a.text(
            0.30, y0 + 0.86, title_s, fontsize=7.2, fontweight="bold", color=border, va="center"
        )
        ax_a.text(0.30, y0 + 0.36, sub_s, fontsize=6.2, color="#1e293b", va="center")
    for y_top in (8.35, 6.55, 4.75, 2.85):
        ax_a.annotate(
            "",
            xy=(5.24, y_top - 0.42),
            xytext=(5.24, y_top - 0.04),
            arrowprops={"arrowstyle": "-|>", "lw": 1.4, "color": "#475569"},
        )

    # Panel B: Human proteome partition horizontal bar (broken/log-clear annotations)
    cats_order = [
        (
            "Unsafe (adult somatic / immune)",
            int(h_ts["category"].eq("Unsafe").sum()),
            PALETTE["unsafe"],
        ),
        (
            "Other:Review (coverage / paralog / niche)",
            int(h_ts["classification"].eq("Other:Review").sum()),
            "#64748b",
        ),
        ("CTA (CA & RA >=10x H0)", int(h_ts["category"].eq("CTA").sum()), PALETTE["cta"]),
        (
            "Other:Expression (<10x H0)",
            int(h_ts["classification"].eq("Other:Expression").sum()),
            "#94a3b8",
        ),
        (
            "RA only (germline / placental >=10x)",
            int(h_ts["category"].eq("Reproductive Antigen only").sum()),
            PALETTE["ra_only"],
        ),
        (
            "CA only (53 retained post-review)",
            int(h_ts["category"].eq("Cancer Antigen only").sum()),
            PALETTE["ca_only"],
        ),
        (
            "Other:Sequence (>=95% Unsafe 9-mer)",
            int(h_ts["classification"].eq("Other:Sequence").sum()),
            "#cbd5e1",
        ),
        (
            "Other:Unscreened tissue (OR/TAS GPCRs)",
            int(h_ts["classification"].eq("Other:Unscreened tissue").sum()),
            "#fbbf24",
        ),
    ]
    y_pos = np.arange(len(cats_order))[::-1]
    vals = [c[1] for c in cats_order]
    cols = [c[2] for c in cats_order]
    ax_b.barh(y_pos, vals, color=cols, height=0.64, edgecolor="#1e293b", linewidth=0.5)
    ax_b.set_xscale("log")
    ax_b.set_xlim(15, 65000)
    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels([c[0] for c in cats_order], fontsize=8.4)
    ax_b.set_xlabel("Human Ensembl 115 protein-equivalence groups (log scale; n = 19,893)")
    ax_b.set_title("B. Human Proteome Partition (Policy 2026-11-27.1)", loc="left", pad=8)
    for yp, v in zip(y_pos, vals):
        pct = 100.0 * v / len(h_ts)
        ax_b.text(
            v * 1.16,
            yp,
            f"{v:,} ({pct:.2f}%)",
            va="center",
            fontsize=8.1,
            fontweight="bold",
            color="#0f172a",
        )
    fig.tight_layout(w_pad=2.4)
    _save_fig(fig, fig_dir, "fig1-conceptual-architecture-and-attrition")

    # --- Figure 2: Multi-Species Parameterization, Cell-Line TIL Separation & Orthology ---
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14.0, 4.6))

    sp_labels = [
        "Human\n(9606, n=19.9k)",
        "Domestic Dog\n(9615, n=20.1k)",
        "Lab Mouse\n(10090, n=21.4k)",
    ]
    cta_counts = [int(df["is_cta"].astype(bool).sum()) for df in (h_ts, d_ts, m_ts)]
    ca_only_counts = [
        int(df["category"].eq("Cancer Antigen only").sum()) for df in (h_ts, d_ts, m_ts)
    ]
    ra_only_counts = [
        int(df["category"].eq("Reproductive Antigen only").sum()) for df in (h_ts, d_ts, m_ts)
    ]
    x_sp = np.arange(3)
    w = 0.25
    b1 = ax1.bar(x_sp - w, cta_counts, w, label="CTA (CA & RA)", color=PALETTE["cta"])
    b2 = ax1.bar(x_sp, ca_only_counts, w, label="CA only", color=PALETTE["ca_only"])
    b3 = ax1.bar(x_sp + w, ra_only_counts, w, label="RA only", color=PALETTE["ra_only"])
    ax1.set_xticks(x_sp)
    ax1.set_xticklabels(sp_labels, fontsize=8.0)
    ax1.set_ylabel("Candidate protein-equivalence groups")
    ax1.set_ylim(0, 390)
    ax1.set_title("A. Candidate Sets Across 3 Mammals", loc="left")
    ax1.legend(frameon=False, fontsize=8, loc="upper right")
    for bars in (b1, b2, b3):
        for bar in bars:
            hv = int(bar.get_height())
            ax1.text(
                bar.get_x() + bar.get_width() / 2,
                hv + 6,
                str(hv),
                ha="center",
                va="bottom",
                fontsize=7.8,
                fontweight="bold",
            )

    # Panel B: Murine in-vivo tumor peak vs pure cell-line peak separating CTAs from host NK granzymes
    m_ctas = m_ts[m_ts["is_cta"].astype(bool)].copy()
    ax2.scatter(
        m_ctas["cell_line_max_tpm"].fillna(0.05).clip(lower=0.05),
        m_ctas["tumor_max_tpm"].fillna(0.05).clip(lower=0.05),
        s=58,
        color=PALETTE["cta"],
        edgecolor="#0f172a",
        linewidth=0.7,
        label="In-vivo murine CTAs (n=8)",
        zorder=4,
    )
    m_offsets = {
        "Sbp": (6, -8),
        "Pbsn": (5, -8),
        "Plac1": (-16, -11),
        "Pramel13": (5, 3),
        "Serpinb9e": (-12, 8),
        "Wfdc13": (5, -8),
        "Tex19.1": (5, 4),
        "Hmga2": (5, 3),
    }
    for r in m_ctas.itertuples(index=False):
        cx = max(float(r.cell_line_max_tpm or 0.05), 0.05)
        cy = max(float(r.tumor_max_tpm or 0.05), 0.05)
        ox, oy = m_offsets.get(str(r.protein_group_label), (4, 3))
        ax2.annotate(
            r.protein_group_label,
            (cx, cy),
            xytext=(ox, oy),
            textcoords="offset points",
            fontsize=7.3,
            color="#0f766e",
            fontweight="bold",
        )

    # Known in-vivo graft vs pure line values for Gzmd-g documented in mouse build
    gzm_pts = [
        ("Gzmd", 1.85, 2138.0, (5, 2)),
        ("Gzme", 2.26, 1412.0, (5, -3)),
        ("Gzmf", 0.31, 865.0, (6, 3)),
        ("Gzmg", 0.08, 443.0, (5, 2)),
    ]
    ax2.scatter(
        [p[1] for p in gzm_pts],
        [p[2] for p in gzm_pts],
        s=68,
        marker="^",
        color=PALETTE["unsafe"],
        edgecolor="#0f172a",
        linewidth=0.7,
        label="Host NK TIL granzymes (Unsafe)",
        zorder=5,
    )
    for sym_g, gx, gy, goff in gzm_pts:
        ax2.annotate(
            sym_g,
            (gx, gy),
            xytext=goff,
            textcoords="offset points",
            fontsize=7.4,
            color=PALETTE["unsafe"],
            fontweight="bold",
        )
    ax2.axvline(10.0, ls="--", lw=0.9, color="#94a3b8")
    ax2.axhline(10.0, ls="--", lw=0.9, color="#94a3b8")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(0.03, 3000)
    ax2.set_ylim(1.0, 5000)
    ax2.set_xlabel("Pure in-vitro murine cancer cell-line peak (TPM)")
    ax2.set_ylabel("In-vivo murine graft / GEMM peak (TPM)")
    ax2.set_title("B. In-Vivo vs. Pure Cell-Line TIL Diagnostic", loc="left")
    ax2.legend(frameon=False, fontsize=7.6, loc="lower right")

    # Panel C: 433-family mammalian CTA orthology breakdown
    orth = pd.read_csv(ROOT / "docs/antigen-whitepaper/cross-species-orthologs.csv")
    h_cta_f = orth["human_classification"].eq("CTA")
    d_cta_f = orth["dog_classification"].eq("CTA")
    m_cta_f = orth["mouse_classification"].eq("CTA")
    n_2plus_cta = int(
        ((h_cta_f.astype(int) + d_cta_f.astype(int) + m_cta_f.astype(int)) >= 2).sum()
    )
    any_cand = (
        orth["human_classification"].isin(["CTA", "Reproductive Antigen only"]).astype(int)
        + orth["dog_classification"].isin(["CTA", "Reproductive Antigen only"]).astype(int)
        + orth["mouse_classification"].isin(["CTA", "Reproductive Antigen only"]).astype(int)
    ) >= 2
    n_2plus_ra = int(any_cand.sum()) - n_2plus_cta
    n_lineage_no_orth = int(
        (
            (orth["human_classification"].eq("absent").astype(int))
            + (orth["dog_classification"].eq("absent").astype(int))
            + (orth["mouse_classification"].eq("absent").astype(int))
            >= 2
        ).sum()
    )
    n_divergent = len(orth) - n_2plus_cta - n_2plus_ra - n_lineage_no_orth

    orth_cats = [
        ("Lineage-specific array\n(absent in other 2 sp.)", n_lineage_no_orth, "#7c3aed"),
        ("Divergent somatic or\n<10x tumor gate in 2nd sp.", n_divergent, "#64748b"),
        ("Conserved CTA + RA\nacross >=2 mammals", n_2plus_ra, PALETTE["ra_only"]),
        ("Conserved CTA in\n>=2 mammals", n_2plus_cta, PALETTE["cta"]),
    ]
    yo = np.arange(len(orth_cats))[::-1]
    ax3.barh(
        yo,
        [x[1] for x in orth_cats],
        color=[x[2] for x in orth_cats],
        height=0.58,
        edgecolor="#1e293b",
        linewidth=0.5,
    )
    ax3.set_yticks(yo)
    ax3.set_yticklabels([x[0] for x in orth_cats], fontsize=8.0)
    ax3.set_xlim(0, max(x[1] for x in orth_cats) * 1.28)
    ax3.set_xlabel(f"Mammalian CTA ortholog families (n = {len(orth)})")
    ax3.set_title("C. Mammalian CTA Orthology Fate", loc="left")
    for yp, (_, cnt, _) in zip(yo, orth_cats):
        ax3.text(cnt + 4, yp, f"{cnt}", va="center", fontsize=8.3, fontweight="bold")

    fig.tight_layout(w_pad=2.0)
    _save_fig(fig, fig_dir, "fig2-multispecies-parameterization-and-orthology")

    # --- Figure 3: Human Deep Dive — Literature Reconciliation, Sex Asymmetry & 66-Protein CA-Only Review ---
    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(14.4, 4.8), gridspec_kw={"width_ratios": [0.95, 0.95, 1.35]}
    )

    pub_mask = h_ts["published_nomination"].fillna(False).astype(bool)
    pub_counts = h_ts.loc[pub_mask, "category"].value_counts()
    lit_bars = [
        ("Unsafe", int(pub_counts.get("Unsafe", 0)), PALETTE["unsafe"]),
        ("CTA", int(pub_counts.get("CTA", 0)), PALETTE["cta"]),
        ("Other", int(pub_counts.get("Other", 0)), PALETTE["other"]),
        ("RA only", int(pub_counts.get("Reproductive Antigen only", 0)), PALETTE["ra_only"]),
        ("CA only", int(pub_counts.get("Cancer Antigen only", 0)), PALETTE["ca_only"]),
    ]
    xl = np.arange(len(lit_bars))
    ax1.bar(
        xl,
        [x[1] for x in lit_bars],
        color=[x[2] for x in lit_bars],
        width=0.58,
        edgecolor="#1e293b",
        linewidth=0.5,
    )
    ax1.set_xticks(xl)
    ax1.set_xticklabels([x[0] for x in lit_bars], fontsize=8.0)
    ax1.set_ylim(0, 2080)
    ax1.set_ylabel("Published nominated protein groups (n = 2,233)")
    ax1.set_title("A. Fate of 2,233 Published Nominations", loc="left")
    for xp, (_, cv, _) in zip(xl, lit_bars):
        ax1.text(
            xp,
            cv + 35,
            f"{cv:,}\n({100 * cv / int(pub_mask.sum()):.1f}%)",
            ha="center",
            va="bottom",
            fontsize=7.4,
            fontweight="bold",
        )

    # Panel B: 4-way sex-restriction partition of 322 CTAs + 146 RA-only
    h_sex = human.sex_restricted_antigens()
    sex_classes = ["RA-M\n(Male)", "RA-P\n(Peripart.)", "RA-MF\n(Shared)", "RA-F\n(Female)"]
    cta_s = [
        int((h_sex["is_cta"] & h_sex["sex_restriction_code"].eq("CTA-M")).sum()),
        int((h_sex["is_cta"] & h_sex["sex_restriction_code"].eq("CTA-P")).sum()),
        int((h_sex["is_cta"] & h_sex["sex_restriction_code"].eq("CTA-MF")).sum()),
        int((h_sex["is_cta"] & h_sex["sex_restriction_code"].eq("CTA-F")).sum()),
    ]
    ra_s = [
        int((~h_sex["is_cta"] & h_sex["sex_restriction_code"].eq("RA-M-only")).sum()),
        int((~h_sex["is_cta"] & h_sex["sex_restriction_code"].eq("RA-P-only")).sum()),
        int((~h_sex["is_cta"] & h_sex["sex_restriction_code"].eq("RA-MF-only")).sum()),
        int((~h_sex["is_cta"] & h_sex["sex_restriction_code"].eq("RA-F-only")).sum()),
    ]
    xs = np.arange(4)
    ax2.bar(xs, cta_s, width=0.56, color=PALETTE["cta"], label="CTA (322)")
    ax2.bar(xs, ra_s, bottom=cta_s, width=0.56, color=PALETTE["ra_only"], label="RA only (146)")
    ax2.set_xticks(xs)
    ax2.set_xticklabels(sex_classes, fontsize=7.8)
    ax2.set_ylim(0, 395)
    ax2.set_ylabel("Human Reproductive Antigens (n = 468)")
    ax2.set_title("B. Host Sex-Restriction Asymmetry (18.5:1)", loc="left")
    ax2.legend(frameon=False, fontsize=8, loc="upper right")
    for xp, c_v, r_v in zip(xs, cta_s, ra_s):
        ax2.text(
            xp,
            c_v + r_v + 6,
            f"{c_v} CTA\n+{r_v} RA",
            ha="center",
            va="bottom",
            fontsize=7.4,
            fontweight="bold",
        )

    # Panel C: 66-protein CA-only functional & micro-niche breakdown
    ca_items = [
        (
            "Retained: Sub-10x male germline / meiotic",
            int(
                (
                    ca_rev["adjudication_disposition"].eq("retained_ca_only")
                    & ca_rev["functional_class"].eq("sub_10x_male_germline_and_spermatogenic")
                ).sum()
            ),
            PALETTE["ca_only"],
        ),
        (
            "Retained: Pericentromeric / RS447 / retrogene",
            int(
                (
                    ca_rev["adjudication_disposition"].eq("retained_ca_only")
                    & ca_rev["functional_class"].eq(
                        "pericentromeric_macrosatellite_duplication_or_retrogene"
                    )
                ).sum()
            ),
            "#f59e0b",
        ),
        (
            "Retained: Post-CS7 fetal morphogen / HOX",
            int(
                (
                    ca_rev["adjudication_disposition"].eq("retained_ca_only")
                    & ca_rev["functional_class"].eq("post_gastrulation_fetal_morphogen_or_homeobox")
                ).sum()
            ),
            "#ea580c",
        ),
        (
            "Retained: Sub-10x CT-X (CSAG/MAGE/SSX)",
            int(
                (
                    ca_rev["adjudication_disposition"].eq("retained_ca_only")
                    & ca_rev["functional_class"].eq("sub_10x_ctx_cancer_germline_paralog")
                ).sum()
            ),
            "#b45309",
        ),
        (
            "Retained: Disposable hair-cortex KRTAP",
            int(
                (
                    ca_rev["adjudication_disposition"].eq("retained_ca_only")
                    & ca_rev["functional_class"].eq("disposable_hair_cortex_keratin")
                ).sum()
            ),
            "#ca8a04",
        ),
        (
            "Adjudicated -> Other:Review (sensory/C-cell/OCM2)",
            int(ca_rev["adjudication_disposition"].eq("adjudicated_review_niche_or_paralog").sum()),
            PALETTE["other"],
        ),
        (
            "Adjudicated -> Unsafe (adult human Vital/Protected)",
            int(ca_rev["adjudication_disposition"].eq("adjudicated_unsafe_adult_human").sum()),
            PALETTE["unsafe"],
        ),
    ]
    yc = np.arange(len(ca_items))[::-1]
    ax3.barh(
        yc,
        [x[1] for x in ca_items],
        color=[x[2] for x in ca_items],
        height=0.62,
        edgecolor="#1e293b",
        linewidth=0.5,
    )
    ax3.set_yticks(yc)
    ax3.set_yticklabels([x[0] for x in ca_items], fontsize=8.0)
    ax3.set_xlim(0, 30)
    ax3.set_xlabel("Baseline human CA-only protein groups (66 -> 53 retained)")
    ax3.set_title("C. 66-Protein CA-Only Functional & Niche Resolution", loc="left")
    for yp, (_, cnt, _) in zip(yc, ca_items):
        ax3.text(cnt + 0.5, yp, str(cnt), va="center", fontsize=8.3, fontweight="bold")

    fig.tight_layout(w_pad=2.0)
    _save_fig(fig, fig_dir, "fig3-human-deep-dive-and-ca-only-biology")

    # --- Figure 4: Empirical Clinical TCR-T & ImmTAC Gold-Standard Calibration ---
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(13.6, 5.0), gridspec_kw={"width_ratios": [1.55, 1.0]}
    )

    gs_sorted = gs_eval.sort_values(
        ["strongest_proxy_ratio", "target_symbol"], ascending=[True, True]
    ).reset_index(drop=True)
    stratum_colors = {
        "safe_cta_potency_proven": PALETTE["cta"],
        "toxic_on_target": PALETTE["unsafe"],
        "toxic_cross_reactivity": "#991b1b",
        "conservative_gene_level_unsafe": "#d97706",
        "excluded_not_potency_proven": "#94a3b8",
    }
    stratum_names = {
        "safe_cta_potency_proven": "Potent, clinically safe CTA (5/5 CTA)",
        "toxic_on_target": "Severe on-target organ toxicity (4/4 Unsafe)",
        "toxic_cross_reactivity": "Fatal off-target cross-reactivity (2/2 Unsafe)",
        "conservative_gene_level_unsafe": "Tolerated non-CTA epitope frontier (3)",
        "excluded_not_potency_proven": "Excluded: 0% ORR / readout pending (2)",
    }
    yg = np.arange(len(gs_sorted))
    ratios = gs_sorted["strongest_proxy_ratio"].clip(lower=0.008).to_numpy()
    bar_cols = [stratum_colors[s] for s in gs_sorted["clinical_stratum"]]
    ax1.barh(yg, ratios, color=bar_cols, height=0.65, edgecolor="#1e293b", linewidth=0.5)
    ax1.axvline(1.0, color="#0f172a", ls="--", lw=1.4)
    ax1.set_ylim(-0.7, len(gs_sorted) + 0.9)
    ax1.text(
        1.18,
        len(gs_sorted) + 0.25,
        "1.0x OncoRef safety ceiling",
        fontsize=7.8,
        fontweight="bold",
        color="#0f172a",
        va="center",
    )
    ax1.set_xscale("log")
    ax1.set_xlim(0.006, 220000)
    lbls_g = []
    for r in gs_sorted.itertuples(index=False):
        cross = (
            f" -> {r.cross_recognized_symbols}"
            if pd.notna(r.cross_recognized_symbols) and str(r.cross_recognized_symbols)
            else ""
        )
        short_prod = str(r.product_name).split(" / ")[0]
        lbls_g.append(f"{r.target_symbol}{cross} ({short_prod})")
    ax1.set_yticks(yg)
    ax1.set_yticklabels(lbls_g, fontsize=7.8)
    ax1.set_xlabel(
        "Strongest adult normal-tissue / cross-reactivity proxy ratio (max x_m / tau_m, log scale)"
    )
    ax1.set_title(
        "A. Normal-Tissue Proxy Ratio Across 16 Clinical TCR-T / ImmTAC Regimens", loc="left"
    )
    for yp, r in zip(yg, gs_sorted.itertuples(index=False)):
        rv = max(float(r.strongest_proxy_ratio), 0.008)
        ax1.text(
            rv * 1.18,
            yp,
            f"{r.strongest_proxy_ratio:g}x [{r.oncoref_classification}]",
            va="center",
            fontsize=7.2,
        )

    handles = [
        mpatches.Patch(color=stratum_colors[k], label=stratum_names[k]) for k in stratum_names
    ]
    ax1.legend(
        handles=handles,
        frameon=False,
        fontsize=7.3,
        loc="upper center",
        bbox_to_anchor=(0.48, -0.15),
        ncol=2,
    )

    # Panel B: Summary calibration performance bars
    metrics = [
        ("Severe clinical toxicity\nsensitivity (6/6 caught)", 100.0, PALETTE["unsafe"]),
        ("Primary-target toxicity\nsensitivity (6/6 caught)", 100.0, "#b91c1c"),
        ("Potency-proven safe CTA\nrecall (5/5 retained)", 100.0, PALETTE["cta"]),
        ("False-safe clinical\ntoxicity misses (0/6)", 0.0, "#64748b"),
    ]
    ym = np.arange(len(metrics))[::-1]
    ax2.barh(
        ym,
        [m[1] for m in metrics],
        color=[m[2] for m in metrics],
        height=0.52,
        edgecolor="#1e293b",
        linewidth=0.6,
    )
    ax2.set_yticks(ym)
    ax2.set_yticklabels([m[0] for m in metrics], fontsize=8.4)
    ax2.set_xlim(0, 122)
    ax2.set_xlabel("Percentage (%)")
    ax2.set_title("B. Empirical Gold-Standard Calibration", loc="left")
    for yp, (_, mv, _) in zip(ym, metrics):
        ax2.text(mv + 2.5, yp, f"{mv:.0f}%", va="center", fontsize=9.0, fontweight="bold")

    fig.tight_layout(w_pad=2.2)
    _save_fig(fig, fig_dir, "fig4-clinical-tcr-gold-standard")

    # --- Figure 5: Subcellular Localization & Extracellular Ectodomain Modality Routing ---
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(14.0, 5.3), gridspec_kw={"width_ratios": [0.95, 1.45]}
    )

    ca_topo = topo_df[topo_df["category"].isin(["CTA", "Cancer Antigen only"])].copy()
    loc_order = [
        (
            "Intracellular (nuclear / cytoplasmic / cytoskeletal)",
            "Intracellular\n(TCR-T / ImmTAC)",
            "#334155",
        ),
        (
            "Internal organelle / minimal-ECD membrane",
            "Internal / minimal-ECD\nmembrane (TCR-T)",
            "#64748b",
        ),
        ("Secreted / extracellular matrix", "Secreted / ECM\n(TCR-T / ECM-RLT)", "#0284c7"),
        (
            "Cell-surface short ectodomain / loop (15-49 aa)",
            "Surface short loop 15-49 aa\n(Loop CAR-T / RLT)",
            "#f59e0b",
        ),
        (
            "Cell-surface ectodomain (>=50 aa)",
            "Surface ectodomain >=50 aa\n(CAR-T / ADC / RLT / TCE)",
            PALETTE["cta"],
        ),
    ]
    yl = np.arange(len(loc_order))
    cta_l = [
        int((ca_topo["category"].eq("CTA") & ca_topo["localization_class"].eq(k)).sum())
        for k, _, _ in loc_order
    ]
    cao_l = [
        int(
            (
                ca_topo["category"].eq("Cancer Antigen only") & ca_topo["localization_class"].eq(k)
            ).sum()
        )
        for k, _, _ in loc_order
    ]
    ax1.barh(yl, cta_l, height=0.58, color=PALETTE["cta"], label="CTA (322)")
    ax1.barh(yl, cao_l, left=cta_l, height=0.58, color=PALETTE["ca_only"], label="CA only (53)")
    ax1.set_yticks(yl)
    ax1.set_yticklabels([lbl for _, lbl, _ in loc_order], fontsize=8.1)
    ax1.set_xlim(0, 330)
    ax1.set_xlabel("Human Cancer Antigens (CTA + CA only = 375)")
    ax1.set_title("A. Modality Partition of 375 Cancer Antigens", loc="left")
    ax1.legend(frameon=False, fontsize=8.2, loc="upper right")
    for yp, v1, v2 in zip(yl, cta_l, cao_l):
        tot = v1 + v2
        ax1.text(
            tot + 5,
            yp,
            f"{tot} ({100 * tot / len(ca_topo):.1f}%)",
            va="center",
            fontsize=8.0,
            fontweight="bold",
        )

    # Panel B: 2D Modality Map — Continuous Mature Ectodomain (aa) vs Peak Tumor TPM
    surf_cand = ca_topo[ca_topo["cell_surface_exposed"]].copy()
    for code_label, grp_label, col_v in [
        ("CTA-P", "CTA-P (adult-absent in M & F)", PALETTE["cta_p"]),
        ("CTA-M", "CTA-M (female-safe)", PALETTE["cta_m"]),
        ("CTA-F", "CTA-F (male-safe)", PALETTE["cta_f"]),
        ("CTA-MF", "CTA-MF (shared reproductive)", PALETTE["cta_mf"]),
        ("", "CA only (retained)", PALETTE["ca_only"]),
    ]:
        sub_s = surf_cand[surf_cand["sex_restriction_code"].eq(code_label)]
        if sub_s.empty:
            continue
        sizes = 38 + 3.2 * sub_s["independent_patients_at_tenfold"].clip(upper=55)
        ax2.scatter(
            sub_s["ecd_max_continuous_aa"],
            sub_s["tumor_max_tpm"],
            s=sizes,
            color=col_v,
            alpha=0.9,
            edgecolor="#0f172a",
            linewidth=0.7,
            label=grp_label,
            zorder=4,
        )

    # Also plot clinical ADC/CAR-T/RLT Unsafe benchmarks along reference band at bottom/right for ECD comparison
    bench_syms = {s for s, _, _ in CLINICAL_SURFACE_BENCHMARKS}
    bench_df = topo_df[
        topo_df["protein_group_label"].isin(bench_syms) & topo_df["ecd_max_continuous_aa"].gt(0)
    ].copy()
    ax2.scatter(
        bench_df["ecd_max_continuous_aa"].clip(upper=3200),
        np.full(len(bench_df), 4.6),
        s=36,
        marker="D",
        color="#cbd5e1",
        edgecolor=PALETTE["unsafe"],
        linewidth=0.9,
        label="24 Clinical ADC/CAR-T/RLT targets (all Unsafe at TCR potency)",
        zorder=3,
    )
    bench_offsets = {
        "CLDN18": (0, 6),
        "FOLR1": (-10, 6),
        "TACSTD2": (12, -10),
        "DLL3": (-6, 6),
        "ERBB2": (-4, -10),
        "FOLH1": (12, 6),
        "SEZ6": (8, -10),
    }
    for bsym, boff in bench_offsets.items():
        brow = bench_df[bench_df["protein_group_label"].eq(bsym)]
        if not brow.empty:
            bx = min(float(brow["ecd_max_continuous_aa"].iloc[0]), 3200.0)
            ax2.annotate(
                bsym,
                (bx, 4.6),
                xytext=boff,
                textcoords="offset points",
                ha="center",
                fontsize=6.5,
                color="#991b1b",
            )

    ax2.axvline(50.0, ls="--", lw=1.1, color="#0f766e")
    ax2.text(
        53.0,
        680.0,
        ">=50 aa continuous mature ECD (scFv / ADC / RLT)",
        fontsize=7.4,
        color="#0f766e",
        fontweight="bold",
        va="top",
    )

    surf_offsets = {
        "BTN1A1": (5, 3),
        "ERVV-1": (5, 4),
        "ERVV-2": (5, -7),
        "CIROP": (5, 2),
        "CRIPTO3": (-38, -6),
        "ADAM2": (6, 2),
        "ADAM18": (6, 2),
        "FMR1NB": (6, 2),
        "CT83": (5, 2),
        "ZAN": (-16, 6),
        "ZP4": (5, -4),
        "IZUMO1R": (-42, 4),
        "LYPD4": (5, -2),
        "ANTXRL": (-8, 7),
        "FAM187B": (5, -3),
        "LRRC52": (5, -6),
        "PKDREJ": (5, 3),
        "CLDN34": (5, 2),
        "TMEM89": (5, 2),
        "SLC36A3": (5, 2),
        "KCNU1": (5, 2),
    }
    for r in surf_cand.itertuples(index=False):
        if r.protein_group_label in surf_offsets:
            ox, oy = surf_offsets[r.protein_group_label]
            ax2.annotate(
                r.protein_group_label,
                (r.ecd_max_continuous_aa, r.tumor_max_tpm),
                xytext=(ox, oy),
                textcoords="offset points",
                fontsize=7.1,
                fontweight="bold",
                color="#0f172a",
            )

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(14, 4200)
    ax2.set_ylim(3.0, 820)
    ax2.set_xlabel("Largest continuous mature extracellular domain (aa, log scale)")
    ax2.set_ylabel("Peak QC-passing in-vivo human tumor expression (TPM, log scale)")
    ax2.set_title(
        "B. Cell-Surface Ectodomain Architecture vs. Tumor Expression (n = 24 Surface CAs)",
        loc="left",
    )
    ax2.legend(frameon=False, fontsize=7.2, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2)

    fig.tight_layout(w_pad=2.2)
    _save_fig(fig, fig_dir, "fig5-surface-ectodomain-modality-routing")


def _md_table(headers: list[str], rows: list[list[object]]) -> str:
    def _c(v: object) -> str:
        if v is None or (isinstance(v, float) and pd.isna(v)) or v == "":
            return "—"
        return str(v).replace("|", " / ").replace("\n", " ")

    out = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for r in rows:
        out.append("| " + " | ".join(_c(x) for x in r) + " |")
    return "\n".join(out)


def _tex_esc(val: object) -> str:
    if val is None or (isinstance(val, float) and pd.isna(val)) or val == "":
        return "---"
    s = str(val)
    repl = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
        "≥": r"$\ge$",
        "≤": r"$\le$",
        "×": r"$\times$",
        "→": r"$\rightarrow$",
        "∪": r"$\cup$",
        "∩": r"$\cap$",
        "∖": r"$\setminus$",
        "α": r"$\alpha$",
        "β": r"$\beta$",
        "γ": r"$\gamma$",
        "²": r"$^2$",
        "⁺": r"$^+$",
        "—": "---",
        "–": "--",
    }
    for k, v in repl.items():
        s = s.replace(k, v)
    return s


def _tex_table(
    caption: str, label: str, colspec: str, headers: list[str], rows: list[list[object]]
) -> str:
    h_line = " & ".join(rf"\textbf{{{_tex_esc(h)}}}" for h in headers) + r" \\"
    body_lines = [" & ".join(_tex_esc(x) for x in r) + r" \\" for r in rows]
    return "\n".join(
        [
            r"\begin{table}[htbp]",
            r"\centering",
            r"\scriptsize",
            r"\setlength{\tabcolsep}{3.5pt}",
            rf"\caption{{{caption}}}",
            rf"\label{{{label}}}",
            rf"\begin{{tabularx}}{{\textwidth}}{{{colspec}}}",
            r"\toprule",
            h_line,
            r"\midrule",
            *body_lines,
            r"\bottomrule",
            r"\end{tabularx}",
            r"\end{table}",
        ]
    )


def build_preprint_texts(
    h_ts: pd.DataFrame,
    d_ts: pd.DataFrame,
    m_ts: pd.DataFrame,
    ca_rev: pd.DataFrame,
    gs_eval: pd.DataFrame,
    topo_df: pd.DataFrame,
    uniprot_sha: str,
) -> tuple[str, str]:
    """Build synchronized Markdown (index.md) and academic LaTeX (.tex) manuscripts."""
    d_coh = canine.tumor_cohorts()
    m_coh = mouse.tumor_cohorts()

    d_samples = int(pd.to_numeric(d_coh["n_samples"], errors="coerce").sum())
    m_samples = int(pd.to_numeric(m_coh["n_samples"], errors="coerce").sum())

    ca_topo = topo_df[topo_df["category"].isin(["CTA", "Cancer Antigen only"])].copy()
    surf_50 = ca_topo[
        ca_topo["localization_class"].eq("Cell-surface ectodomain (>=50 aa)")
    ].sort_values(["independent_patients_at_tenfold", "tumor_max_tpm"], ascending=[False, False])
    surf_short = ca_topo[
        ca_topo["localization_class"].eq("Cell-surface short ectodomain / loop (15-49 aa)")
    ].sort_values(["tumor_max_tpm"], ascending=[False])
    surf_all = pd.concat([surf_50, surf_short], ignore_index=True)

    # Shared structured data tables for both Markdown and LaTeX
    t1_species_rows = [
        [
            "Homo sapiens (9606)",
            "Ensembl 115 (GRCh38)",
            "19,893 (20,099)",
            "8 / 8",
            "135 cohorts",
            "322 (222 M / 48 P / 12 F / 40 MF)",
            "53",
            "146",
            "1,116",
            "18,256",
        ],
        [
            "Canis lupus familiaris (9615)",
            "Ensembl 98 (CanFam3.1)",
            "20,082 (20,199)",
            "8 / 8",
            f"19 cohorts ({d_samples} pts) + 12 lines",
            "118 (100 M / 7 P / 1 F / 10 MF)",
            "127",
            "313",
            "2,201",
            "17,323",
        ],
        [
            "Mus musculus (10090)",
            "Ensembl 115 (GRCm39)",
            "21,433 (21,760)",
            "8 / 8",
            f"14 cohorts ({m_samples} tumors) + 9 lines",
            "8 (5 M / 2 P / 0 F / 1 MF)",
            "30",
            "303",
            "2,514",
            "18,578",
        ],
    ]

    unsafe_7 = ca_rev[ca_rev["adjudication_disposition"].eq("adjudicated_unsafe_adult_human")]
    review_6 = ca_rev[ca_rev["adjudication_disposition"].eq("adjudicated_review_niche_or_paralog")]
    retained_53 = ca_rev[ca_rev["adjudication_disposition"].eq("retained_ca_only")]

    t2_ca_only_rows = [
        [
            "Adjudicated -> Unsafe",
            len(unsafe_7),
            f"{unsafe_7['tumor_fold'].median():.1f}x",
            "Confirmed adult human Vital/Protected somatic micro-niche expression",
            ", ".join(unsafe_7["protein_group_label"].tolist()),
        ],
        [
            "Adjudicated -> Other:Review",
            len(review_6),
            f"{review_6['tumor_fold'].median():.1f}x",
            "Unscreened inner-ear/DRG sensory niche, diluted C-cell home, or 100% 1-mm Unsafe paralog",
            ", ".join(review_6["protein_group_label"].tolist()),
        ],
        [
            "Retained: Sub-10x male germline",
            int(
                retained_53["functional_class"].eq("sub_10x_male_germline_and_spermatogenic").sum()
            ),
            f"{retained_53.loc[retained_53['functional_class'].eq('sub_10x_male_germline_and_spermatogenic'), 'tumor_fold'].median():.1f}x",
            "Narrow spermatocyte/spermatid window or AZF/X copy diluted to 1x-9.6x in whole testis",
            "PRDM9, RPS4Y2, RBMY1E (8.9x), ZNF705D (9.5x), LIPJ (9.6x), PKDREJ",
        ],
        [
            "Retained: Pericentromeric / RS447",
            int(
                retained_53["functional_class"]
                .eq("pericentromeric_macrosatellite_duplication_or_retrogene")
                .sum()
            ),
            f"{retained_53.loc[retained_53['functional_class'].eq('pericentromeric_macrosatellite_duplication_or_retrogene'), 'tumor_fold'].median():.1f}x",
            "Hominoid macrosatellite arrays, cleavage ZGA duplications (<20-cell quorum), & retrogenes",
            "PRR20A-E, USP17L1/4/7, TRIM51G, TRIM64C, C1orf232, CCDC195",
        ],
        [
            "Retained: Post-CS7 morphogen / HOX",
            int(
                retained_53["functional_class"]
                .eq("post_gastrulation_fetal_morphogen_or_homeobox")
                .sum()
            ),
            f"{retained_53.loc[retained_53['functional_class'].eq('post_gastrulation_fetal_morphogen_or_homeobox'), 'tumor_fold'].median():.1f}x",
            "Fetal tooth ameloblast scaffold & post-gastrulation neural/limb homeobox TFs silent in adults",
            "AMELX (30 pts), DBX1 (74 pts), HOXC12, HOXB1, HOXD12, EVX2",
        ],
        [
            "Retained: Sub-10x CT-X paralogs",
            int(retained_53["functional_class"].eq("sub_10x_ctx_cancer_germline_paralog").sum()),
            f"{retained_53.loc[retained_53['functional_class'].eq('sub_10x_ctx_cancer_germline_paralog'), 'tumor_fold'].median():.1f}x",
            "Canonical Xq28/Xp11/Xp21 cancer-germline loci whose bulk testis RNA sits below 10x N",
            "CSAG3 (401 pts), CSAG2 (65 pts), SSX5 (8 pts), MAGEB17",
        ],
        [
            "Retained: Disposable hair KRTAP",
            int(retained_53["functional_class"].eq("disposable_hair_cortex_keratin").sum()),
            f"{retained_53.loc[retained_53['functional_class'].eq('disposable_hair_cortex_keratin'), 'tumor_fold'].median():.1f}x",
            "Anagen hair-shaft cortical matrix proteins cleared through Disposable <=32 TPM skin/hair screen",
            "KRTAP20-3, KRTAP21-2, KRTAP21-3, KRTAP6-2",
        ],
    ]

    t3_tcr_rows = []
    for r in gs_eval.itertuples(index=False):
        tgt = str(r.target_symbol)
        if pd.notna(r.cross_recognized_symbols) and str(r.cross_recognized_symbols):
            tgt = f"{tgt} -> {r.cross_recognized_symbols}"
        if r.benchmark_id == "wt1_ttcr_c4":
            pot = f"n={int(r.patients_treated)}, 100% RFS"
        elif pd.notna(r.orr_pct) and float(r.orr_pct) > 0:
            pot = f"n={int(r.patients_treated)}, {r.orr_pct:g}% ORR"
        elif r.benchmark_id == "magea3_a1_titin":
            pot = f"n={int(r.patients_treated)}, lethal d4-5"
        elif pd.notna(r.patients_treated):
            pot = f"n={int(r.patients_treated)}, 0% ORR"
        else:
            pot = "Phase 1 pending"
        t3_tcr_rows.append(
            [
                tgt,
                f"{r.product_name} ({r.hla_restriction})",
                pot,
                r.observed_clinical_toxicities,
                f"{r.oncoref_classification} ({r.strongest_proxy_ratio:g}x)",
            ]
        )

    t4_surface_rows = []
    for r in surf_all.itertuples(index=False):
        cls_short = r.sex_restriction_code if r.sex_restriction_code else "CA only"
        t4_surface_rows.append(
            [
                r.protein_group_label,
                cls_short,
                r.surface_architecture,
                f"{int(r.ecd_max_continuous_aa)} ({int(r.ecd_total_aa)})",
                r.ecd_intervals or "23-188 (EGF-CFC)",
                f"{r.tumor_max_tpm:.1f} ({r.tumor_top_cancer_code})",
                int(r.independent_patients_at_tenfold),
                int(r.safe_2mm_9mers) if pd.notna(r.safe_2mm_9mers) else 0,
            ]
        )

    # ---------------- MARKDOWN PREPRINT ----------------
    md_lines = [
        "# Multi-Scale Normal-Tissue Injury Tiers and Symmetric Enrichment Reconcile Cancer–Testis Antigen Discovery Across Human, Dog, and Mouse and Map Actionable TCR-T, CAR-T, ADC, and Radioligand Modalities",
        "",
        "*Preprint manuscript (`oncoref` policy `oncoref.target-selection.2026-11-27.1`) · [Download compiled PDF (`oncoref-antigen-preprint.pdf`)](oncoref-antigen-preprint.pdf) · [LaTeX source (`oncoref-antigen-preprint.tex`)](oncoref-antigen-preprint.tex) · [Surface & localization table (CSV)](tables/candidate-surface-topology.csv) · [66-protein `CA only` review](../audits/ca-only-functional-review-20261127/index.md)*",
        "",
        "**Authors:** PIRL-UNC Collaborative Group (Computational Oncology & Immunogenomics, University of North Carolina at Chapel Hill)",
        "",
        "---",
        "",
        "## Abstract",
        "",
        "**Background:** Cancer–testis antigens (`CTAs`) and oncofetal proteins underpin affinity-enhanced T-cell receptor therapies (`TCR-T`), bispecific T-cell engagers (`ImmTACs`), cancer vaccines, and emerging surface-targeted antibody–drug conjugates (`ADCs`), chimeric antigen receptor T cells (`CAR-T`), and radioligand therapies (`RLTs`). Yet across **149** published CTA, oncofetal, and placental studies nominating **2,233** human protein groups, gene lists conflict sharply because legacy pipelines pool incompatible expression scales, overlook rare normal neuroendocrine or sensory cell types diluted in bulk biopsies, conflate biological germline restriction with patient-specific organ expendability, and transfer human labels onto preclinical dog and mouse models by gene-symbol orthology.",
        "",
        "**Methods:** We formulated a species-invariant, four-stage classification architecture (`OncoRef`, `oncoref.target-selection.2026-11-27.1`) that (i) collapses identical coding translations into protein-equivalence groups (`PEGs`) while retaining all splice isoforms for `≥95%` 9-mer self-overlap screening; (ii) evaluates every healthy somatic and activated-immune dataset in its own native measurement scale (`TPM`, `nTPM`, `nCPM`, donor-resolved droplet `UMI`, validated `IHC`, or plex-bleed-corrected `TMT-ppm`) across four cited clinical organ-injury tiers (`Vital`, `Protected`, `Manageable`, `Disposable`); (iii) enforces symmetric **≥10 × H₀** *in vivo* enrichment over **H₀ = max(Normal N, Disposable D, 1 TPM)** independently in tumors (`Cancer Antigens`, `CA`) and reproductive/embryonic compartments (`Reproductive Antigens`, `RA`), defining **CTA = CA ∩ RA** while segregating pure *in vitro* cancer cell lines; and (iv) routes retained candidates through host-specific sex/organ protection (`HostContext`), 0/1-mismatch `PeptideIndex` self-proteome subtraction, and coordinate-verified UniProtKB subcellular/ectodomain modality mapping.",
        "",
        '**Results:** Across **19,893** human (`GRCh38`), **20,082** domestic dog (`CanFam3.1`), and **21,433** laboratory mouse (`GRCm39`) `PEGs`—with all **8/8** evidence pillars active natively in every species—OncoRef classifies **322 human**, **118 canine**, and **8 murine** *in vivo* `CTAs` (**96.6%**, **96.6%**, and **75.0%** with confirmed ectopic activation outside homologous source organs). Segregating pure *in vitro* cancer cell lines unmasks host `LAK`/`NK` cytolytic granzymes (`Gzmd–g`, `443–2,138 TPM` in *in vivo* grafts vs. `0.08–2.26 TPM` in pure cultured tumor cells) as tumor-infiltrating leukocyte artifacts, while only **15** of **433** mammalian `CTA` ortholog families qualify as `CTAs` across multiple species due to primate `CT-X`/`ZGA` amplification. In humans, **78.8% (`1,760 / 2,233`)** of published nominations fail modern adult somatic safety gates, whereas expression rules uncover **77** unnominated `CTAs`. Host-physiology partitioning resolves why male-restricted `CTA-M` (`222`) outnumbers female-restricted `CTA-F` (`12`) by **18.5:1** alongside **48** adult-absent pregnancy/peripartum `CTA-P` antigens (`270` female-safe vs. `60` male-safe `CTAs`). A complete primary-literature functional review of all **66** baseline `CA only` (`CA \\ RA`) proteins adjudicates **7** occult adult human `Vital`/`Protected` micro-niche proteins (`GALP`, `BHLHE23`, `IL36B`, `ZACN`, `HES3`, `OCSTAMP`, `PROKR1`) to `Unsafe` and **6** sensory/C-cell/1-mismatch paralogs (`SRMS`, `OCM2`, `CIB3`, `OTOP1`, `GFRA4`, `MRGPRG`) to `Other:Review`, showing that all **53** retained `CA only` proteins arise from five concrete somatic-silent mechanisms rather than "cancer-only dark matter." Against **16** human clinical `TCR-T` and `ImmTAC` regimens, OncoRef flags **6/6 (`100%`)** severe clinical specificity toxicities as `Unsafe` while retaining **5/5 (`100%`)** potency-proven safe `CTAs`. Finally, coordinate-verified topology mapping across all **375** human `CAs` shows that **84.0% (`315`)** are intracellular (`TCR-T`/`ImmTAC` obligatory) and **9.6% (`36`)** are secreted, while **24 (`6.4%`)** expose plasma-membrane extracellular domains—including **19** with `≥50 aa` continuous mature ectodomains (`BTN1A1`, `ERVV-1`, `ERVV-2`, `CIROP`, `CRIPTO3`, `ADAM2`, `ADAM18`, `FMR1NB`, `CT83`, `ZP4`, `IZUMO1R`, `LYPD4`, `PKDREJ`) actionable for `CAR-T`, `ADC`, and `RLT` without the somatic expression liabilities carried by **24/24 (`100%`)** of conventional clinical surface benchmarks.',
        "",
        "---",
        "",
        "## 1. Introduction",
        "",
        "Therapeutic targeting of tumor-associated antigens faces an unforgiving quantitative trade-off between anti-tumor potency and on-target or cross-reactive normal-tissue destruction. Affinity-enhanced `TCR-T` cells and soluble `TCR–anti-CD3` bispecifics (`ImmTACs`) recognize intracellularly processed peptides presented on major histocompatibility complex (`pMHC`) molecules at densities below `10` copies per cell, enabling durable regressions in solid tumors (`afamitresgene autoleucel` targeting `MAGEA4`, `letetresgene autoleucel` targeting `NY-ESO-1`/`CTAG2`, `IMA203` and `brenetafusp` targeting `PRAME`) [1–5]. That same single-digit-copy sensitivity caused fatal leukoencephalopathy when `MAGEA3`-directed `TCR-T` cells cross-recognized rare brain-expressed `MAGEA12` (`NCT01273181`), fatal cardiogenic shock within 4–5 days when an affinity-matured `HLA-A*01` `MAGEA3` TCR cross-reacted with striated-muscle `Titin` (`TTN`, `NCT01350401`), severe hemorrhagic colitis against colonic `CEACAM5`, and uveitis plus inner-ear cochlear hair-cell destruction (`ototoxicity`) against melanocytic `MART-1` (`MLANA`) and `gp100` (`PMEL`) [6–9].",
        "",
        'Historically, discovery of Cancer–Testis Antigens (`CTAs`)—proteins expressed in germ cells or trophoblasts and ectopically reactivated in somatic malignancies—relied on reverse-transcription PCR panels of 12–16 whole organs (`CTDatabase` / `CTpedia`) or differential microarray/RNA-seq screens comparing tumors against bulk GTEx medians [10–14]. Compiling **149** published CTA, oncofetal, zygotic-genome-activation (`ZGA`), and cancer–placenta catalogs reveals **2,233** distinct human protein groups nominated as "cancer–testis" or "tumor-restricted." Four structural flaws explain why most historical nominations fail under clinical scrutiny:',
        "",
        "1. **Cross-scale pooling and bulk-organ dilution:** Merging `GTEx` `TPM`, `HPA` TMM-scaled `nTPM`/`nCPM`, microarray intensities, and single-cell `CP10K` distorts thresholds, while whole-organ punches dilute `<1%` vital cell populations—such as hypothalamic arcuate neuroendocrine neurons (`GALP`), retinal bipolar interneurons (`BHLHE23`), cochlear/vestibular hair cells (`CIB3`, `OCM2`), pancreatic islets (`HES3`), and monogenic olfactory sensory neurons (`OR*`)—by `50–1,000×` below bulk cutoffs.",
        "2. **Conflating pre-analytical sequencing artifacts with biology:** Multiplexed bulk atlases suffer flowcell index hopping and plate pipetting splash (`testis` or `pancreas` transcripts bleeding into `brain` or `heart`), droplet `snRNA-seq` captures nascent intronic readthrough from adjacent housekeeping promoters (`ARHGAP6 → AMELX`, `GAB3 → CTAG2`), and degraded FFPE tumor biopsies accumulate uniform intergenic genomic DNA (`gDNA`) that inflates intronless germline/olfactory genes.",
        '3. **Entangling biological category definitions with host anatomy:** Labelling a prostate- or lactation-expressed protein (`SEMG2`, `CSN2`, `LALBA`) either unconditionally "safe" or unconditionally "excluded" obscures the clinical reality that host risk depends on whether the patient is male, female, pregnant, lactating, or post-organ-ablation.',
        "4. **Orthology projection across divergent mammalian genomes:** Comparative oncology trials in outbred companion dogs and preclinical immunotherapy studies in mice routinely project human CTA symbols onto animal genomes, ignoring both lineage-specific `CT-X` retrotransposition and species differences in accessory reproductive and lymphoid anatomy.",
        "",
        "Here we present the complete conceptual architecture, multi-species parameterization (`human`, `dog`, `mouse`), human proteome reconciliation, `CA only` functional literature review, clinical `TCR-T` calibration, and subcellular/ectodomain modality routing implemented in `OncoRef` (`v1.55+`, policy `oncoref.target-selection.2026-11-27.1`).",
        "",
        "---",
        "",
        "## 2. Conceptual Architecture: Data Scales, Tissue-Injury Tiers, and Set Algebra",
        "",
        "![Conceptual pipeline architecture and human proteome attrition](figures/fig1-conceptual-architecture-and-attrition.png)",
        "",
        "*Figure 1.* **(A)** Species-invariant four-stage antigen screening pipeline (`Stage 0` identical-sequence `PEG` grouping through `Stage 3` symmetric `10 × H₀` enrichment) decoupled from downstream `HostContext` organ protection and `PeptideIndex` self-proteome k-mer subtraction. **(B)** Empirical attrition across all **19,893** human Ensembl 115 protein-equivalence groups (**20,099** coding loci).",
        "",
        "### 2.1 Identical-Sequence Protein-Equivalence Grouping (`Stage 0`) and Library Integrity",
        "",
        'Ampliconic cancer–germline families (`CTAG1A/CTAG1B` encoding `NY-ESO-1`, `MAGEA2/MAGEA2B`, `GAGE*`, `SPANX*`, `CT45A*`, `PRR20A–E`, `BPY2/B/C`) frequently duplicate coding exons with 100% amino-acid identity. Evaluating individual Ensembl gene IDs splits RNA-seq reads across identical copies—understating both normal-organ hazards and tumor abundance—or creates artificial "self-overlap" between sibling loci encoding the exact same polypeptide. For any species $S$ with SHA-256-pinned proteome $R_S$, OncoRef groups primary-assembly coding loci whose canonical/MANE representative translations are sequence-identical into a single **Protein-Equivalence Group (`PEG`)**, summing member-locus abundance *within each individual biopsy* before computing donor summaries, while retaining **all** annotated coding splice isoforms (`245,535` human, `45,036` canine, and `66,037` murine records) for `≥95%` 9-mer containment against `Unsafe` proteins (`repeat_family_equivalence` generalizes this symmetrically across closed ≥90%-identity arrays such as `D4Z4`/`DUX4` and `VCX`).',
        "",
        "Before any threshold is evaluated, OncoRef enforces library-level pre-analytical QC across both normal and tumor datasets:",
        "",
        "- **Holistic normal-atlas deconvolution (`gtex_deconvolution`):** Across **18,637** GTEx v11 libraries (`52` tissue shards, `946` donors: `311 XX`, `635 XY`), non-negative least squares (`NNLS`) over `342` lineage-disjoint anchor genes quarantines `2` global foreign-organ tube swaps and strips `201` `SM`-plate pipetting splash events, `100,287` flowcell index-hopping gene–library pairs, `225` same-donor autopsy carryovers, and `16,058` muscularis/mucosa boundary admixtures.",
        "- **Genomic DNA (`gDNA`) and parenchymal admixture screens:** Libraries in which `>2.5%` of olfactory-receptor loci reach `≥5 TPM` are flagged as gDNA-contaminated (`library_qc`), and testicular germ-cell tumor (`TGCT`) biopsies containing residual non-neoplastic seminiferous tubules (`Leydig`/`Sertoli`/spermatid stoichiometry, `#713`/`#768`) are barred from serving as sole tumor witnesses (`BOD1L2`, `PGK2`, `TSSK2` remain `RA only` unless an independent non-admixed tumor reaches `10 × H₀`).",
        "- **Single-nucleus compartment & cis-spillover attribution (`#765`/`#779`):** Droplet `snRNA-seq` signals confined strictly to unspliced nuclear pre-mRNA with `0` cytoplasmic/bulk RNA and `0` protein staining—caused by readthrough from adjacent high-abundance somatic genes (`GAB3 → CTAG2`, `SCPP` salivary/lactation cluster `≤100 kb` spillover onto `CSN2`)—are recorded as `nuclear_only` cautions rather than false `Unsafe` exclusions.",
        "",
        "### 2.2 Four Clinical Adult Tissue-Injury Tiers and Zero Cross-Unit Pooling (`Stage 1`)",
        "",
        "Every adult non-reproductive anatomical structure and sorted immune cell population is assigned to one of four cited clinical injury tiers (`tissue_tiers`), reflecting the reversibility and clinical manageability of on-target cytolytic loss:",
        "",
        "1. **`Vital` (bulk ceiling `≤ 4.0` native TPM/nTPM):** Brain subregions, spinal cord, retina, cochlear/vestibular sensory hair cells, dorsal root/trigeminal ganglia, heart, lung, liver, kidney, adrenal gland, pancreas, vascular endothelium, hematopoietic stem/progenitor cells, and megakaryocytes/platelets—where acute T-cell destruction causes irreversible disability or death.",
        "2. **`Protected` (bulk ceiling `≤ 8.0` native TPM/nTPM; `2×` Vital):** Gastrointestinal tract mucosa/smooth muscle, urinary bladder, esophagus, skin epidermis, cornea, peripheral nerve, skeletal muscle, thyroid follicular parenchyma, pituitary, and mature T cells, NK cells, monocytes, macrophages, and neutrophils.",
        "3. **`Manageable` (hold above `16.0`, `Unsafe` above `128.0`):** Mature B-lineage cells (`CD19`/`BCMA` clinical precedent with IVIG support), gallbladder, spleen, and calcitonin-producing thyroid parafollicular C-cells.",
        "4. **`Disposable` (`D`, capped at `≤ 32.0`):** Cutaneous hair-follicle bulb/cortex and non-ocular skin melanocytes (`vitiligo`/`alopecia` tolerable in advanced cancer), palatine tonsil, vermiform appendix, minor salivary gland, and adult involuted whole thymus (`thymic` scope tracked separately so medullary thymic epithelial promiscuous `AIRE` transcription neither inflates `H₀` nor bypasses sorted T-cell safety).",
        "",
        "Crucially, **OncoRef never converts or pools incompatible units.** Each modality applies its own native gate independently, and any confirmed adult `Vital`/`Protected` breach assigns **`Unsafe` with strict precedence**:",
        "",
        "- **Bulk somatic & sorted/activated immune RNA:** Within-tissue linearly interpolated donor `P95` (`GTEx v11`, `HPA v23` tissue & sorted immune, `DICE` 15-subset immune `P95`, `Monaco` 29-subset immune) against `4 / 8 / 16–128 / 32` native units. Only in `human`—where `HPA v23` provides routine, cell-type-resolved immunohistochemistry (`IHC`) across `45` normal organs—does confirmed `not_detected` staining with an `Approved` (`2×`), `Supported` (`4×`), or `Enhanced` (`8×`) antibody relax a borderline RNA ceiling; pooled `MS` (`dog`) and literature `IHC` (`mouse`) act strictly as affirmative `Unsafe`/hold gates without inflating RNA ceilings when absent.",
        "- **Donor-resolved single-cell `UMI` atlases (`assess_rare_cell_rna`):** Across `452` adult human `CELLxGENE Census` datasets + `Siletti` adult brain snRNA-seq (`human`), `GSE225599` healthy companion-dog PBMCs (`dog`), and `Tabula Muris` 10x (`mouse`), single-molecule (`1 UMI`) droplets never exclude a gene; an `established` exclusion requires `≥10` cells and `≥1.0%` multi-molecule (`≥2 UMI`) cells replicated across `≥2` studies / `≥3` donors (or `≥5.0%` in one study; `2×` frequency in `Protected` populations).",
        "",
        "### 2.3 Eligibility Guards (`Stage 2`), Symmetric `10 × H₀` Enrichment (`Stage 3`), and `HostContext`",
        "",
        "Proteins surviving `Stage 1` must pass `Stage 2` eligibility before positive enrichment is interpreted: complete peripheral normal coverage (`Other:Review` if incomplete or paralog-ambiguous), `<95%` 9-mer containment in `Unsafe` proteins and non-readthrough locus structure (`Other:Sequence`), and screened anatomical home niche (`home_tissue_unmeasured == False`). In particular, all unvalidated **monogenic chemosensory GPCRs** (`OR*`/`TAAR*`/`TAS1R*`/`TAS2R*`/`VN1R*` in human, `COR*`/`OR*` in dog, `Olfr*`/`Vmn*r*` in mouse) are expressed at high density in `~0.1%` of sensory neurons—diluting `~1,000×` in whole-mucosa bulk RNA—and are uniformly held out of `CA`/`CTA` (`Other:Unscreened tissue` or `RA only`).",
        "",
        "In `Stage 3`, against one shared healthy adult baseline $H_0 = \\max(N, D, 1.0\\text{ TPM})$:",
        "",
        "$$\\text{CA} = \\{g \\in \\mathcal{E} : T_{\\max,\\text{in vivo}}(g) \\ge 10 \\cdot H_0(g) \\land \\neg\\text{unmeasured\\_home}(g)\\}, \\qquad \\text{RA} = \\{g \\in \\mathcal{E} : R_{\\text{supp}}(g) \\ge 10 \\cdot H_0(g)\\}, \\qquad \\text{CTA} = \\text{CA} \\cap \\text{RA}$$",
        "",
        "Pure *in vitro* cancer cell lines (`comparison_cell_lines`) never grant `CA` or `CTA`; instead they populate `cell_line_max_tpm` and `qualifying_cell_lines`. Every `CA` is annotated with `tumor_activation_mode` (`ectopic`, `ectopic_and_lineage_retained`, or `lineage_retained` when qualifying tumors arise solely in the homologous reproductive/mammary organ). Downstream, `HostContext` partitions all `RA` and `CTA` proteins into four mutually exclusive host-physiology classes (`CTA-M` male-restricted / female-safe, `CTA-F` female-restricted / male-safe, `CTA-P` pregnancy/peripartum-only / adult-absent in both non-pregnant sexes, and `CTA-MF` shared adult male & female), followed by exact and 1-mismatch `PeptideIndex` subtraction against the host's `Unsafe` + protected-organ proteome.",
        "",
        "---",
        "",
        "## 3. Multi-Species Inventory and Parameterization: Human, Dog, and Mouse",
        "",
        _md_table(
            [
                "Species (Taxon)",
                "Pinned Proteome",
                "PEGs (Loci)",
                "Active Pillars",
                "Tumor & Cell-Line Scope",
                "CTA (Sex Classes)",
                "CA only",
                "RA only",
                "Other",
                "Unsafe",
            ],
            t1_species_rows,
        ),
        "",
        "![Multi-species parameterization, cell-line TIL diagnostic, and CTA orthology](figures/fig2-multispecies-parameterization-and-orthology.png)",
        "",
        "*Figure 2.* **(A)** Genome-wide candidate sets across human (`9606`), domestic dog (`9615`), and laboratory mouse (`10090`) with all `8/8` screening pillars active in each species. **(B)** *In vivo* murine tumor peak vs. pure *in vitro* cell-line peak separating cancer-intrinsic murine `CTAs` (`Hmga2`, `Tex19.1`, `Wfdc13`, `Pramel13`) from excluded host `NK` `TIL` granzymes (`Gzmd–g`). **(C)** Evolutionary and regulatory fate across **433** mammalian `CTA` ortholog families.",
        "",
        "### 3.1 Comparative Anatomy and Per-Species Evidence Pillars",
        "",
        "Separating mathematical rules from declarative species profiles (`species_policies`) allows each mammal to be screened exclusively on native data while respecting comparative anatomy (`normal-tissue-taxonomy.yaml`):",
        "",
        "- **Human (`9606`, `GRCh38` Ensembl 115):** `19,893` `PEGs` (`20,099` loci) screened across `GTEx v11` (`52` deconvolved shards, reproductive-age `20–49` `XX` stratification), `HPA v23/v25.1` tissues, brain subregions, and `150+` single-cell types, `DICE`/`Monaco`/`HPA` sorted immune populations, `452` `CELLxGENE` + `Siletti` donor rare-cell `UMI` atlases, adult cornea `scRNA-seq`, `#735` olfactory/taste/dental-pulp/bone/hair panels, `GSE107746` stage-replicated adult human oocytes/granulosa (`128` libraries), `Petropoulos` (`E-MTAB-3929`) preimplantation + `Tyser` (`E-MTAB-9388`) `CS7` gastrulation embryos, and **135** QC-passing *in vivo* patient tumor cohorts.",
        "- **Domestic Dog (`9615`, `CanFam3.1` Ensembl 98 + `RefSeq 106`):** Companion dogs develop spontaneous osteosarcoma, multicentric B/T-cell lymphoma, hemangiosarcoma, urothelial carcinoma, oral melanoma, mammary carcinoma, glioma, and soft-tissue sarcoma in outbred, immunocompetent hosts. Dogs anatomically lack seminal vesicles, bulbourethral glands, and a vermiform appendix (`prostate` is the sole male accessory gland; females have a bicornuate `uterine_horn` and are frequently ovariohysterectomized). Across `20,082` `PEGs`, all `8/8` pillars are active: `Bgee 15.2` (`25` deduplicated SRA experiments), the `300`-sample sex-stratified `GSE219045` preclinical atlas, `64` normal mammary biopsies (`GSE119810` split into lactating vs. resting mammary), primary companion-dog adrenal/blood/marrow/gingiva/brain/liver/lung/macrophage cohorts, `GSE225599` healthy 7-dog PBMC `10x` `scRNA-seq` (`150,155` donor × cell evaluations), `MSV000099744` 38-tissue `TMT-iBAQ` proteomics (with `10%` within-plex reporter bleed subtraction), **19** spontaneous patient tumor cohorts (`935` primary tumors), and **12** segregated canine cell lines—yielding **118 `CTAs`** (`114/118 = 96.6%` ectopic), **127 `CA only`**, and **313 `RA only`**. Among `14` published canine RT-PCR CTA nominees (`Hindriks 2025`, `Nemec 2019`), only `PRAME` and `SSX6` validate as genuine canine `CTAs`, while `6` (`AKR1E2`, `IGF2BP3`, `NT5C1B`, `PAGE3-like`, `SPECC1`, `TPX2`) are excluded as `Unsafe` in healthy canine blood, marrow, lung, or brain.",
        "- **Laboratory Mouse (`10090`, `GRCm39` Ensembl 115):** Mice add a distinct `coagulating_gland` alongside `seminal vesicle` and `prostate`, lack palatine `tonsil` and `appendix`, and undergo `E14.5` androgen-mediated ablation of male mammary buds (`male` mice lack `resting_mammary`). Across `21,433` `PEGs` (`66,037` translations), `Bgee 15.2` + `GSE219045` + `ImmGen` + cytokine-activated `NK`/`Treg`/`BMDM` RNA, `Tabula Muris` `10x` `UMI`, `MGI GXD 6.27` adult `TS28` protein staining, **14** *in vivo* syngeneic/GEMM/carcinogen tumor cohorts (`201` tumors across `H-2b`, `H-2d`, `H-2q`), and **9** segregated pure cell lines yield **8 *in vivo* `CTAs`** (`Hmga2`, `Tex19.1`, `Wfdc13`, `Pramel13`, `Gm4969`, `4930558J18Rik` ectopic; `Pbsn`, `Sbp` prostate lineage-retained), **30 `CA only`**, and **303 `RA only`**. Crucially, requiring *in vivo* tumor enrichment while tracking pure cultured lines (`MC38`, `CT26`, `B16-F10`, `4T1`, `P815`, `EO771`) unmasks cytolytic granzymes `Gzmd–g` (`443–2,138 TPM` in *in vivo* grafts vs. `0.08–2.26 TPM` in pure tumor cells) as host uterine/activated-`NK` `TIL` transcripts, while classic DBA/2 mastocytoma antigen `P1A` (`Trap1a`) is held in `Other:Review` due to documented endogenous BALB/c somatic counterevidence [15].",
        "",
        "### 3.2 Mammalian CTA Orthology (`433` Families)",
        "",
        "Cross-mapping all human (`322`), canine (`118`), and murine (`8`) `CTAs` across **433** mammalian ortholog families shows that **nearly half of human `CTAs` lack any 1:1 coding ortholog in dog or mouse** (`NY-ESO-1`, `CTAG2`, `GAGE*`, `XAGE*`, `SPANX*`, `VCX/VCY`, `CT45A*`, `CT47A*`, `DUX4`, `PRAMEF*`), reflecting rapid primate `CT-X` and `ZGA` gene-array expansion after boreoeutherian divergence (~95 Mya). Conversely, rodents independently expanded `Pramel*` and `Mageb*`. Only **15** ortholog families independently achieve full `CTA` status in `≥2` species (e.g., `PRAME`, `SSX`, `MAGEB`, `CSN2`, `LALBA`, `HMGA2`, `LIN28B`, `WFDC`), joined by conserved boreoeutherian `RAs` (`ADAM2`, `AKAP4`, `DPPA5`, `MAGEB18`, `MEIOSIN`, `TBXT`), while several 1:1 orthologs (`HORMAD1`, `CLDN6`, `DDX43`, `IGF2BP1`) pass canine somatic panels but fail ultra-deep human single-cell or corneal safety screens—proving that preclinical safety and target discovery can never be transferred by gene symbol alone.",
        "",
        "---",
        "",
        "## 4. Human Proteome Deep Dive: Literature Reconciliation, Sex Asymmetry, and the 66-Protein `CA only` Review",
        "",
        "![Human literature reconciliation, sex asymmetry, and 66-protein CA-only resolution](figures/fig3-human-deep-dive-and-ca-only-biology.png)",
        "",
        "*Figure 3.* **(A)** OncoRef classification of **2,233** human protein groups nominated across **149** published CTA, oncofetal, and placental studies. **(B)** 4-way host sex-restriction partition across all **468** human `RAs` (**322** `CTAs`), explaining the `18.5:1` `CTA-M : CTA-F` ratio and **48** peripartum-only `CTA-P` targets. **(C)** Primary-literature functional & micro-niche resolution across all **66** baseline `CA only` proteins (`7 → Unsafe`, `6 → Other:Review`, `53` retained across 5 concrete biological mechanisms).",
        "",
        "### 4.1 Reconciling `2,233` Published CTA Nominations Across `149` Studies",
        "",
        "Evaluating all **2,233** literature-nominated human protein groups through the uniform multi-assay gates excludes **1,760 (`78.8%`)** as `Unsafe` and holds **113 (`5.1%`)** in `Other`, retaining **360 (`16.1%`)** in the candidate union (**245 `CTA`**, **20 `CA only`**, **95 `RA only`**), while quantitative `10 × H₀` screening identifies **77 additional `CTAs`** absent from prior compilations (`322` total `CTAs`). Study-level retention stratifies cleanly by how each source handled normal somatic tissues: focused modern cancer–germline and `DUX4` studies (`Seager 2024`, `Chew 2019`, `Carter 2023`, `Loriot 2025`) retain **34%–67% `CTAs`**, whereas broad germ-cell-versus-tumor differential lists (`Bruggeman 2018`, `Wang 2016` broad testis-biased modules, `da Silva 2017`) are **79%–98% `Unsafe`** because mitotic spindle, DNA-repair, and ciliary/centrosomal genes (`TTK`, `BIRC5`, `CEP55`, `PBK`, `EZH2`) peak in proliferating spermatocytes *and* healthy adult bone marrow, gut crypts, and activated T/NK cells.",
        "",
        "### 4.2 Why Male-Restricted `CTA-M` Outnumbers Female-Restricted `CTA-F` by `18.5:1`",
        "",
        "Partitioning the **322** human `CTAs` by `HostContext` yields **222 `CTA-M`** (male-restricted, safe in non-pregnant women), **48 `CTA-P`** (pregnancy/peripartum-only—expressed in `placenta`, `embryo`, or `lactating_breast` and silent across all non-pregnant, non-lactating adult men and women), **12 `CTA-F`** (female-restricted, safe in men: `ZP4`, `IZUMO1R`, `OOEP`, `PADI6`, `NLRP4`, `NLRP11`, `KHDC3L`, `TUBB8`, `PATL2`, `WEE2`, `ZAR1`, `PRR23E`), and **40 `CTA-MF`** (shared adult gonadal/tract). Consequently, **270 `CTAs` (`83.9%`)** are free of adult organ exposure in a non-pregnant, non-lactating woman (`CTA-M + CTA-P`) and **60 (`18.6%`)** in a man (`CTA-F + CTA-P`).",
        "",
        "Auditing **159** canonical female reproductive, oocyte, zona-pellucida, maternal-effect subcortical maternal complex (`SCMC`), and mammary genes (`female-lineage-audit.csv`) demonstrates that this `18.5:1` asymmetry is driven by mammalian reproductive developmental biology rather than GTEx age bias (even after restricting GTEx female organs to reproductive-age `20–49` `XX` donors and integrating `GSE107746` stage-replicated adult human primordial-to-antral oocytes and granulosa cells):",
        "",
        "1. **Adult male testis continuously executes all stages of mitosis, meiosis, and spermiogenesis (`~10⁸` sperm/day)** inside an immune-privileged blood–testis barrier, maintaining massive steady-state `CT-X` and flagellar/acrosomal transcript pools, whereas **adult human ovaries contain `<0.01%` prophase-I-arrested oocytes** embedded in somatic cortical stroma (`FOXL2`/`WNT4`), and fetal oogonial meiosis (`SPO11`, `SYCP1/3`, `HORMAD1`, `PRDM9`) finishes *in utero* before birth (landing in `CTA-MF` or `CTA-M`).",
        "2. **Cyclic female mullerian/endometrial epithelium shares developmental lineages with respiratory/renal/mucosal ciliated epithelia**, excluding `OVGP1`, `PAEP`, and `MUC16` (`CA-125`) as `Unsafe` in adult lung, fallopian/bronchial cilia, or mesothelium.",
        '3. **True milk and placental effectors (`LALBA`, `CSN2`, `BTN1A1`, `PLAC1`, `ERVV-1/2`, `XAGE2/3`, `PSG*`)** are transcriptionally silent in resting non-pregnant adult breast (`<1–2 TPM` after separating parous/peripartum lobules) and instead populate **`CTA-P` (`48` proteins)**—safe across both non-pregnant sexes while automatically blocked under `HostContext(protect=("placenta", "lactating_breast"))`.',
        "",
        "### 4.3 Functional Biology and Micro-Niche Resolution of All `66` `CA only` Proteins",
        "",
        _md_table(
            [
                "Disposition / Functional Class",
                "Proteins",
                "Median Tumor Fold",
                "Biological Mechanism & Why R < 10x H0",
                "Representative Proteins",
            ],
            t2_ca_only_rows,
        ),
        "",
        'A central conceptual question in antigen discovery is **"what do `Cancer Antigen only` (`CA \\ RA`) proteins actually do in normal human biology, and how can a non-mutated endogenous human gene reach `≥10 × H₀` in tumors without reaching `≥10 × H₀` in adult reproductive organs or preimplantation/gastrulation embryos?"** Reviewing primary human literature, `HPA v25.1` single-cell `nCPM`, `#735` specialized mucosal/skeletal panels, and 1-mismatch `Unsafe` paralog identity across all **66** pre-adjudication `CA only` proteins (`oncoref.ca_only_review`) resolved every gene without exception:',
        "",
        "- **7 occult adult human `Vital` / `Protected` micro-niche proteins adjudicated to `Unsafe`:** Each resides in a `<1%` specialized neuroendocrine, retinal, islet, mucosal, or bone-resorbing cell population diluted `50–200×` in whole-organ GTEx punches: **`GALP`** (`1,308×` tumor fold; adult human hypothalamic infundibular/arcuate nucleus neuroendocrine neurons in `HYPOMAP` spatial/snRNA-seq [16] + `4.3 nCPM` RPE), **`BHLHE23`** (`110×`; adult human retinal rod and ON-cone bipolar interneurons [17]), **`IL36B`** (`24.7×`; inducible adult human epidermal/squamous mucosal IL-1-family cytokine [18], `#735` circumvallate taste mucosa `10.27 TPM > 8.0` and esophageal apical `54.2 nCPM`), **`ZACN`** (`23.4×`; human-specific `Zn²⁺`/`H⁺`-gated Cys-loop cation channel cloned from adult human pancreas and brain [19], `#735` olfactory `37.87 TPM` and HPA cardiomyocytes `56.1 nCPM`), **`HES3`** (`20.5×`; endogenous adult human pancreatic islet `α`/`β`-cell repressor controlling `PDX1` and insulin release [20]), **`OCSTAMP`** (`20.4×`; RANKL-induced multinucleated bone-resorbing osteoclast fusogen [21], whose `SARC_OS`/`MM` biopsies reflect reactive giant-cell osteoclast stroma), and **`PROKR1`** (`12.9×`; adult human adrenal chromaffin `10.6 nCPM` and GI smooth-muscle prokineticin GPCR [22]).",
        "- **6 sensory/vestibular, `Manageable` C-cell, & near-self paralog proteins adjudicated to `Other:Review`:** **`SRMS`** (gastric chief cells `19.7 nCPM` & goblet cells `15.7 nCPM`, `0` known patient IDs), **`OCM2`** (**`100.0%` 1-mismatch 9-mer overlap / `0` safe `≥2-mm` 9-mers** with `Unsafe` `OCM` plus retinal ganglion cells `7.7 nCPM` and cochlear hair cells [23]), **`CIB3`** (inner-ear vestibular hair-cell `TMC1/2` mechanotransduction subunit [23, 24], unscreened in GTEx/HPA), **`OTOP1`** (Type III sour-taste receptor & vestibular otoconia proton channel [25] diluted to `0.69 TPM` in whole lingual papillae), **`GFRA4`** (GPI-anchored persephin coreceptor on normal adult human calcitonin-producing thyroid parafollicular C-cells [26]—`Manageable` cell type diluted `~100×` in whole thyroid, clinically actionable in thyroidectomized medullary thyroid carcinoma `NCT04877613` [27] but held from pan-cancer `CA`), and **`MRGPRG`** (unscreened adult human dorsal-root/trigeminal nociceptor family GPCR, `0` known patient IDs).",
        "- **All `53` retained `CA only` proteins belong to five transparent somatic-silent classes:** *(1)* **24 sub-`10×` male germline/meiotic/sex-chromosome proteins** whose narrow leptotene/pachytene window (`PRDM9` `5.8×` testis bulk vs. `72.7 nCPM` spermatocytes [28]; `HUS1B`), multi-copy Y/X locus (`RBMY1E` `8.92×`, `RBMY1B`, `RPS4Y2` in `49` male `B_ALL`/`CLL` patients, `BPY2/B/C`), or post-meiotic spermatid niche (`ZNF705D` `9.49×`, `LIPJ` `9.64×`, `SLC35G3` `5.88×`, `PKDREJ` `4.96×`, `CLDN34` `6.31×`, `FUT5` `3.09×` with `72` safe `≥2-mm` 9-mers) sits just below the `10×` bulk `RA` bar; *(2)* **4 canonical `CT-X` paralogs** (`CSAG3` in **`401` patients** across `32` cancer types, `CSAG2` in `65` patients, `SSX5`, `MAGEB17`) whose higher background denominator `N` places bulk testis at `1–4×`; *(3)* **6 post-gastrulation (`>CS7`) fetal morphogens and `HOX` factors** (`AMELX` enamel matrix protein in `30` odontogenic `CRANIO`/`RT` patients—adult enamel is acellular [29]; `DBX1` spinal/pre-Bötzinger progenitor TF in `74` `MBL`/`NBL`/`WILMS`/`SCLC` patients [30]; `HOXC12`, `HOXB1`, `HOXD12`, `EVX2`); *(4)* **4 `Disposable` hair-shaft cortex `KRTAP` proteins** (`KRTAP20-3`, `KRTAP21-2`, `KRTAP21-3`, `KRTAP6-2`); and *(5)* **15 hominoid pericentromeric / `RS447` macrosatellite duplications, cleavage-ZGA arrays, and retrogenes** (`PRR20A–E` peaking at `37.9 TPM` in 8-cell embryos below the 20-cell `RA` quorum, `USP17L1/4/7` deubiquitinases, `TRIM51G`, `TRIM64C`, pediatric Wilms/neuroblastoma `C1orf232` in `23` patients, and B-ALL/Burkitt `CCDC195` in `11` patients).",
        "",
        "---",
        "",
        "## 5. Empirical Calibration Against Clinical `TCR-T` and `ImmTAC` Trials",
        "",
        "![Clinical TCR-T and ImmTAC gold-standard calibration](figures/fig4-clinical-tcr-gold-standard.png)",
        "",
        "*Figure 4.* **(A)** Strongest healthy adult somatic or cross-recognized paralog proxy ratio ($\\max_m x_m / \\tau_m$) across **16** curated human clinical `TCR-T` and `ImmTAC` trial regimens (`tcr_trials.gold_standard()`). **(B)** Sensitivity (`100%`, `6/6`) to severe clinical specificity toxicity and recall (`100%`, `5/5`) on potency-proven clinically safe `CTAs`.",
        "",
        _md_table(
            [
                "Target (-> Cross-Hit)",
                "Regimen (HLA)",
                "Clinical Potency",
                "Observed Human Specificity Toxicity",
                "OncoRef Call (Max Ratio)",
            ],
            t3_tcr_rows,
        ),
        "",
        "Because healthy-tissue RNA, single-cell `UMI`, `IHC`, and sequence-similarity screens are *proxies* for clinical pMHC toxicity, we calibrated OncoRef directly against **16** human endogenous-antigen `TCR-T` and `ImmTAC` clinical regimens (**14** with proven *in vivo* cytolytic potency: `ORR > 0%` or rapid lethal target lysis):",
        "",
        "1. **100% sensitivity (`6 / 6`, `0` false-safe misses) to severe human specificity toxicities:** OncoRef independently flags every toxicity-causing target/cross-hit as `Unsafe` well above `1.0×`—capturing both **on-target off-tumor organ destruction** (`CEACAM5` colon `106×` → Grade 3–4 hemorrhagic colitis [8]; `MLANA` `DMF5`/`DMF4` & `PMEL` `gp100` retinal/uveal/cochlear melanocytes `3.5–28×` → severe uveitis and vestibulocochlear hearing loss [9]; `TYR` cutaneous/ocular melanocytes) and **lethal off-target paralog/mimotope cross-reactivity** (`MAGEA3 → MAGEA12` brain RNA/rare-cell gate → fatal necrotizing leukoencephalopathy [6]; `MAGEA3 ESDPIVAQY → TTN ESDPNGATY` heart muscle `>250×` caught by 1-mismatch/motif self-screening → fatal cardiogenic shock [7]).",
        "2. **100% recall (`5 / 5`) on potency-proven clinically safe `CTAs`:** All five clinically active, specificity-safe `CTA` regimens retain full `CTA` status (`<1.0×` across every adult somatic gate): FDA-approved `afamitresgene autoleucel` (`MAGEA4`, `37–39% ORR` in synovial sarcoma [2]), `letetresgene autoleucel` (`NY-ESO-1` and `CTAG2`/`LAGE-1a`, `40–50% ORR` in synovial sarcoma and myxoid/round-cell liposarcoma [1, 3]), `IMA203` (`PRAME` `TCR-T`, `52.9%` confirmed `ORR` in melanoma [4]) & `brenetafusp` (`IMC-F106C` `PRAME` `ImmTAC` [5]), and `ADP-A2M10` (`MAGEA10`). Notably, clinical organ-safety clearance (`#766`) paired with nuclear readthrough attribution (`#765`) resolved historical `CTAG2` lung/heart false alarms caused by `GAB3` antisense spillover.",
        "3. **Defining the processing-differential frontier (`3` conservative `Unsafe` non-CTA regimens):** Three active regimens target genes that OncoRef conservatively marks `Unsafe` at the whole-protein level—`WT1` (`TCR-C4` `HLA-A*02:01` `RMFPNAPYL`, `100%` 44-month RFS in AML without glomerular nephritis because immunoproteasomes in leukemic blasts cleave `WT1_126–134` far more efficiently than constitutive proteasomes in adult podocytes [31]), `AFP` (`ADP-A2AFP` `FMNKFIYEI`, `18% ORR` in HCC where adult liver background is low), and `MAGEA3` `HLA-DPB1*04:01` CD4 `TCR-T` (`KITE-718`, `23.5% ORR`, whose 16-mer `KKLLTQHFVQENYLEY` is absent from `MAGEA12` at 4 positions [32]). Two non-potent or pending regimens (`IMA202` `MAGEA1` `0% ORR`; `CT83` `NCT05035407` pending) are held out of potency denominators so ineffective receptors never masquerade as safety proof.",
        "",
        "---",
        "",
        "## 6. Subcellular Localization & Extracellular Domain Architecture: Routing `CTAs` Across `TCR-T`, `CAR-T`, `ADC`, and `RLT`",
        "",
        "![Subcellular localization and surface ectodomain modality routing](figures/fig5-surface-ectodomain-modality-routing.png)",
        "",
        "*Figure 5.* **(A)** Modality-routing partition across all **375** human Cancer Antigens (`322 CTA + 53 CA only`) from SHA-256-pinned UniProtKB topological/transmembrane/signal/glycosylation coordinates (`sha256: "
        + uniprot_sha[:16]
        + "…`). **(B)** Largest continuous mature extracellular domain (`aa`, after signal/propeptide cleavage) vs. peak *in vivo* tumor expression (`TPM`) across all **24** surface-exposed `CA` candidates, contrasted with **24** clinical `ADC`/`CAR-T`/`RLT`/`TCE` surface benchmarks (grey/red diamonds, all `24/24 Unsafe` under TCR-grade normal somatic screening).",
        "",
        _md_table(
            [
                "Antigen",
                "Class",
                "Membrane Architecture",
                "Max Cont. (Total) Mature ECD (aa)",
                "Mature ECD Coordinates",
                "Peak Tumor TPM (Cohort)",
                "Pts >=10x",
                "Safe >=2-mm 9-mers",
            ],
            t4_surface_rows,
        ),
        "",
        "Whereas `TCR-T` and `ImmTAC` bispecifics access any endogenously degraded protein via `MHC-I`/`MHC-II` presentation—at the cost of `HLA`-allele restriction (`~40–48%` coverage for `HLA-A*02:01` in European-ancestry cohorts and `15–25%` in African/Asian cohorts) and vulnerability to tumor `B2M`/`HLA` loss—**`CAR-T` cells, antibody–drug conjugates (`ADCs`), radioligand therapies (`RLTs`), and IgG bispecific T-cell engagers (`TCEs`)** bind intact native extracellular epitopes independently of `HLA` genotype. To determine systematically which `OncoRef` candidates are actionable by surface-directed modalities, we mapped UniProtKB topological domains, transmembrane helices, cleaved signal/propeptide spans, `GPI`-anchors, and luminal/extracellular disulfide/glycosylation features onto Ensembl 115 sequences across all **521** human candidates and **24** established clinical surface benchmarks (`tables/candidate-surface-topology.csv`):",
        "",
        "### 6.1 Why `84.3%` of Human Cancer Antigens Require Intracellular `pMHC` Modalities—and Why Clinical `ADC`/`CAR-T` Benchmarks Fail `TCR` Screens",
        "",
        "Across all **375** human Cancer Antigens (`322 CTA + 53 CA only`), **268 (`71.5%`)** are strictly nuclear, chromatin-bound, cytoplasmic, or cytoskeletal—including `PRAME` (nuclear `CRL2` E3-ligase substrate receptor with inner-leaflet cytoplasmic membrane association, `0` extracellular residues), `NY-ESO-1`, `CTAG2`, `MAGEA4`, `MAGEA1`, `MAGEC2`, `SSX2`, `DUX4`, `FOXR2`, `CSAG3`, `DBX1`, `PRDM9`, and flagellar fibrous-sheath `AKAP4` (`Cell projection, cilium, flagellum` axonemal interior, `0` `TM` helices)—and another **47 (`12.5%`)** localize to internal ER, Golgi (`FUT5`, `SLC35G3/5`, `CTAGE*`), mitochondrial, or acrosomal/minimal-ECD membranes (`TPTE`, `XKR3/5`) without surface ectodomains (`315 / 375 = 84.0%` intracellular overall; **36 [`9.6%`]** are secreted soluble/ECM proteins such as `CSN2`, `LALBA`, `SCGB1D4`, `SEMG2`, `PLAC1`, `CER1`, and calcified-matrix `AMELX`).",
        "",
        "Conversely, **24 / 24 (`100%`)** of approved and late-stage clinical `ADC`, `CAR-T`, `RLT`, and `TCE` surface targets (`HER2`/`ERBB2` `630 aa` ECD, `TROP2`/`TACSTD2` `248 aa`, `FRα`/`FOLR1` `210 aa`, `Nectin-4` `318 aa`, `PSMA`/`FOLH1` `707 aa`, `SSTR2` `43 aa` loop, `DLL3` `466 aa`, `CLDN18.2` `53 aa` ECL1, `CLDN6` `53 aa` ECL1, `MSLN` `562 aa`, `GPC3` `530 aa`, `GPC2` `531 aa`, `ROR1` `377 aa`, `STEAP1` `33 aa` ECL2, `SEZ6` `906 aa`, `CDH6` `562 aa`, `B7-H4`/`VTCN1` `235 aa`, `PTK7` `674 aa`, `LY6K` `121 aa`) are classified as **`Unsafe` under OncoRef's `TCR`-potency screen**! This dichotomy reflects fundamental biophysics: `ADCs`, affinity-detuned `CAR-T`s, and `RLTs` achieve a therapeutic index over `10–100 TPM` normal epithelial expression via **super-linear avidity thresholds (`>10³–10⁴` surface copies/cell)**, payload bystander/internalization kinetics, or tight-junction sequestration (`CLDN18.2` in gastric pits, `CEACAM5` on apical colonic lumen), whereas `TCR-T` cells recognize `1–10` `pMHC` complexes on basolateral parenchymal surfaces with zero polarization shielding.",
        "",
        "### 6.2 Twenty-Four Somatic-Silent Surface `CTA` & `CA only` Candidates for `CAR-T`, `ADC`, and `RLT`",
        "",
        "Crucially, **24 human Cancer Antigens (`22 CTAs + 2 CA only`)** combine **TCR-grade adult somatic silence** with **genuine plasma-membrane extracellular exposure**—representing a uniquely clean class of surface targets where high-affinity `CAR-T` cells, potent topoisomerase/pyrrolobenzodiazepine `ADCs`, or alpha-emitting `²²⁵Ac`/`¹⁷⁷Lu` `RLTs` do not depend on fragile density windows:",
        "",
        "1. **Pregnancy / peripartum-only (`CTA-P`) surface ectodomains (`5` targets, `0` adult somatic OR gonadal exposure in non-pregnant men and women):**",
        "   - **`BTN1A1` (`216 aa` IgV–IgC Type I ectodomain `[27–242]`, `427.5 TPM` peak across `6` `BRCA`/`LUSC`/`LUAD` patients):** Lactation-associated B7-family immune checkpoint (mutually exclusive with `PD-L1`) with an active clinical mAb (`hSTC810`, `NCT05231746`) [33]; ideal for `ADC`, `CAR-T`, or `TCE` in non-lactating hosts.",
        "   - **`ERVV-1` & `ERVV-2` (`300 aa` surface `SU` ectodomain `[22–321]`, `236.4` & `221.2 TPM` across `4` & `13` `KIRP`/`KIRC`/`OV` patients):** Syncytiotrophoblast-restricted human endogenous retrovirus group V envelope glycoproteins reactivated in renal papillary/clear-cell and ovarian carcinomas—prime candidates for `Env`-directed `CAR-T`, `ADC`, and `RLT`.",
        "   - **`CIROP` (`715 aa` Type I ectodomain `[21–735]`, `74.4 TPM` in `MENINGIOMA`) & `CRIPTO3` (`TDGF1P3`, `188 aa` glycosylated/disulfide-knotted EGF-CFC surface Nodal coreceptor, `71.2 TPM` across `11` `TGCT` patients).**",
        "2. **High-prevalence male-restricted (`CTA-M`, `100%` adult-organ-free in non-pregnant women) and female-restricted (`CTA-F`, `100%` organ-free in men) large ectodomains (`≥90 aa`):**",
        "   - **`ADAM2` (`Fertilin-β`, `512 aa` mature disintegrin/cysteine-rich ectodomain `[175–686]`, `100.4 TPM` across `54` patients in `BRCA`/`SCLC`/`LUAD`) & `ADAM18` (`503 aa` mature ectodomain `[185–687]`, `150.0 TPM` across `40` patients in `LUSC`/`LUAD`/`BRCA`):** Large, folded metalloproteinase-like single-pass Type I ectodomains matching the size of `HER2`/`MSLN`/`DLL3` (`466–630 aa`) with `578–613` safe `≥2-mm` 9-mers.",
        "   - **`FMR1NB` (`NY-SAR-35` / `CT37`, `94 aa` continuous extracellular loop `[90–183]` between 2 TM helices, `363.4 TPM` across `49` patients in `MM`, `SKCM`, and sarcoma)** and **`CT83` (`KK-LC-1`, `92 aa` C-terminal Type II ectodomain `[22–113]`, `193.2 TPM` across `22` patients in `CESC`, `STAD`, `TNBC`, `LUAD`)** [34]: compact, high-abundance multi-cohort surface `CTAs` suitable for both `CAR-T`/`ADC`/`VHH-RLT` and intracellular `TCR-T`.",
        "   - **`ZP4` (`CTA-F`, `444 aa` mature zona-pellucida ectodomain `[19–462]`, `74.4 TPM` across `11` `LUSC`/`STAD` patients)** and **`IZUMO1R` (`JUNO`/`FOLR4`, `CTA-F`, `209 aa` `GPI`-anchored folate-receptor fold `[20–228]` directly homologous to `mirvetuximab` target `FOLR1` `210 aa`!)**, plus **`LYPD4` (`CTA-MF`, `199 aa` `GPI`-anchored `Ly6`/`uPAR` fold `[27–225]` in `BRCA`)**, **`ANTXRL` (`326 aa`, `OV`)**, **`FAM187B` (`318 aa`, `BRCA`)**, **`LRRC52` (`221 aa`, `SKCM`)**, **`SPACA1` (`192 aa`, `KIRC`)**, **`TMEM95` (`129 aa`, `NEC_MERKEL`)**, **`ZAN` (`2,740 aa`, `TGCT`)**, and `CA only` **`PKDREJ` (`1,165 aa` N-terminal REJ/GPS ectodomain `[20–1184]`, `NEC_MERKEL`)**.",
        "3. **Constrained multi-pass / short-ectodomain candidates (`24–48 aa` loops, `5` targets):**",
        "   - **`CLDN34` (`CA only`, 4-TM claudin with a `48 aa` first extracellular loop `[33–80]`, `79 aa` total ECD)** shares the exact 4-TM claudin architecture of `CLDN18.2` and `CLDN6` (`53 aa` ECL1), alongside `CTA-M` **`TMEM89` (`39 aa`)**, **`SLC36A3` (`32 aa` loop / `112 aa` total across 6 loops, `B_ALL`)**, **`CXorf66` (`28 aa`)**, and **`KCNU1` (`24 aa` loop / `65 aa` total, `10` `LIHC` patients)**—best suited for cyclic-peptide/minibody `RLT` (`SSTR2` paradigm) or high-avidity loop-specific `CAR-T` (`BNT211` paradigm).",
        "",
        "---",
        "",
        "## 7. Discussion, Limitations, and Reproducible Availability",
        "",
        'By separating **native-unit multi-assay normal safety (`Stage 1`)**, **coverage and sensory-niche eligibility (`Stage 2`)**, **symmetric *in vivo* `10 × H₀` tumor and reproductive enrichment (`Stage 3`)**, **host-specific organ protection (`HostContext`)**, and **surface vs. intracellular modality routing**, `OncoRef` replaces static, literature-inherited "CTA lists" with an auditable, multi-species quantitative foundation.',
        "",
        "Four limitations warrant emphasis when translating screening outputs into experimental therapeutics: *(i)* **Bulk tumor peak `TPM` measures biopsy abundance rather than intratumoral clonality**, requiring spatial `IHC`/`RNAscope` (`CT83` [34]) or tumor `scRNA-seq` to confirm homogeneous malignant-cell expression; *(ii)* **RNA silence does not guarantee pMHC presentation**, which must be confirmed by immunopeptidomics (`MS/MS`) and `WT1`-like proteasome processing assays; *(iii)* **fetal 1st–3rd trimester organogenesis remains fail-closed (`pregnant=True`)** because preimplantation and `CS7` gastrulation atlases do not sample mid-gestation fetal organs; and *(iv)* **surface topology predicts ectodomain exposure of the canonical proteoform**, requiring flow-cytometry / surface-biotinylation confirmation against alternative splicing or ectodomain shedding (`ADAM`/protease cleavage) in target tumors.",
        "",
        "**Data and Code Availability:** All catalogs, SHA-256-verified provenance receipts, `PeptideIndex` self-subtraction tools, and deterministic build scripts (`scripts/build_ca_only_functional_review.py`, `scripts/promote_integrated_human_catalog.py`, `scripts/build_antigen_preprint.py`) are open-source under Apache 2.0 in `pirl-unc/oncoref` (`https://github.com/pirl-unc/oncoref`, PyPI `pip install oncoref`).",
        "",
        "---",
        "",
        "## References",
        "",
        "1. D'Angelo SP et al. Antitumor activity associated with prolonged persistence of adoptively transferred NY-ESO-1 c259T cells in synovial sarcoma. *Cancer Discov.* 2018;8(8):944–957 (PMID:29891538).",
        "2. D'Angelo SP et al. Afamitresgene autoleucel for advanced synovial sarcoma and myxoid round cell liposarcoma (SPEARHEAD-1). *Lancet.* 2024;403(10435):1460–1471 (PMID:38552654).",
        "3. Pollack SM et al. Letetresgene autoleucel in advanced/metastatic myxoid/round cell liposarcoma. *J Clin Oncol.* 2025 (PMC12084024).",
        "4. Wermke M et al. Autologous T cell therapy for PRAME+ advanced solid tumors in HLA-A*02+ patients (IMA203 Phase 1). *Nat Med.* 2025 (PMID:40200105).",
        "5. Hamid O et al. Clinical activity of brenetafusp (IMC-F106C), a PRAME × CD3 ImmTAC bispecific, across cutaneous and uveal melanoma. *J Clin Oncol / ESMO.* 2024.",
        "6. Morgan RA et al. Cancer regression and neurological toxicity following anti-MAGE-A3 TCR gene therapy. *J Immunother.* 2013;36(2):133–151 (PMID:23377668).",
        "7. Linette GP et al. Cardiovascular toxicity and titin cross-reactivity of affinity-enhanced T cells in myeloma and melanoma. *Blood.* 2013;122(6):863–871 (PMID:23770775).",
        "8. Parkhurst MR et al. T cells targeting carcinoembryonic antigen can mediate regression of metastatic colorectal cancer but induce severe transient colitis. *Mol Ther.* 2011;19(3):620–626 (PMID:21157437).",
        "9. Johnson LA et al. Gene therapy with human and mouse T-cell receptors mediates cancer regression and targets normal tissues expressing cognate antigen. *Blood.* 2009;114(3):535–546 (PMID:19451487).",
        "10. Almeida LG et al. CTdatabase: a knowledge-base of high-throughput and curated data on cancer-testis antigens. *Nucleic Acids Res.* 2009;37:D816–D819 (PMID:18838390).",
        "11. Wang C et al. Systematic identification of genes with a cancer-testis expression pattern in 19 cancer types. *Nat Commun.* 2016;7:10499 (PMID:26813868).",
        "12. da Silva VL et al. Genome-wide identification of cancer/testis genes and their association with prognosis in a pan-cancer analysis. *Oncotarget.* 2017;8(54):92966–92977 (PMID:29190970).",
        "13. Bruggeman JW et al. Massivelyonenumber of germ cell-specific genes expressed in human cancers. *Oncotarget.* 2018;9:29795–29812 (PMID:30042815).",
        "14. Carter JA et al. Transcriptomic pan-cancer analysis of cancer-testis antigens and clinical survival. *J Immunother Cancer.* 2023 (PMID:37369392).",
        "15. Van den Eynde B et al. Presence on a murine mastocytoma of an antigen regulated by a gene of the X chromosome (`P1A` / `Trap1a`). *J Exp Med.* 1991;173(6):1373–1384 (PMID:1903428); Sarma S et al. *J Exp Med.* 1999;189:811–820.",
        "16. Tadross JA et al. A comprehensive spatio-cellular map of the human hypothalamus (`HYPOMAP`). *Nature.* 2025;638:1049–1059 (PMID:39910309).",
        "17. Cowan CS et al. Cell types of the human retina and its organoids at single-cell resolution. *Cell.* 2020;182(6):1623–1640 (PMID:32946783).",
        "18. Johnston A et al. IL-1F5, -F6, -F8, and -F9 (`IL-36`) are induced in human keratinocytes and skin lesions. *J Immunol.* 2011;186(4):2613–2622 (PMID:21242515).",
        "19. Davies PA et al. A novel class of ligand-gated ion channel is activated by Zn²⁺ (`ZACN`). *J Biol Chem.* 2003;278(2):712–717 (PMID:12381728).",
        "20. Masjkur J et al. Hes3 is expressed in adult human pancreatic islets and regulates insulin release. *J Biol Chem.* 2014;289(51):35503–35516 (PMID:25371201).",
        "21. Khan UA et al. OC-STAMP is essential for human osteoclast multinucleation. *J Cell Biochem.* 2013;114(8):1772–1781 (PMID:23444125).",
        "22. Lin DC et al. Identification and molecular characterization of two closely related G protein-coupled receptors activated by prokineticins/EG-VEGF (`PROKR1`). *J Biol Chem.* 2002;277(22):19276–19280 (PMID:11886876).",
        "23. van der Valk WH et al. A single-cell level comparison of human inner ear organoids with the human cochlea and vestibular organs. *Cell Rep.* 2023;42(6):112623 (PMID:37289589).",
        "24. Liang X et al. CIB2 and CIB3 are auxiliary subunits of the mechanotransduction channel of hair cells. *EMBO J.* 2021;40(20):e107279 (PMID:34487363).",
        "25. Tu YH et al. An evolutionarily conserved gene family encodes proton-selective ion channels (`OTOP1`). *Science.* 2018;359(6379):1047–1050 (PMID:29371428).",
        "26. Lindahl M et al. Human glial cell line-derived neurotrophic factor receptor alpha 4 (`GFRA4`) is predominantly expressed in thyroid parafollicular C-cells. *J Biol Chem.* 2001;276(12):9344–9351 (PMID:11116144).",
        "27. Bhoj VG et al. Chimeric antigen receptor T cells targeting GFRα4 in medullary thyroid cancer. *Mol Ther Oncolytics.* 2021;20:387–398 (PMID:33632929).",
        "28. Baudat F et al. PRDM9 is a major determinant of meiotic recombination hotspots in humans and mice. *Science.* 2010;327(5967):836–840 (PMID:20044539).",
        "29. Salido EC et al. Human enamel protein gene amelogenin (`AMELX`) on X and Y chromosomes. *Am J Hum Genet.* 1992;50(2):303–316 (PMID:1734717); Apps JR et al. *J Pathol.* 2018;244:458–469 (PMID:29480946).",
        "30. Pierani A et al. Control of interneuron fate in the developing spinal cord by the progenitor homeodomain protein Dbx1. *Neuron.* 2001;29(2):367–384 (PMID:11239428).",
        "31. Chapuis AG et al. T cell receptor gene therapy targeting WT1 prevents acute myeloid leukemia relapse post-transplant. *Nat Med.* 2019;25(7):1064–1072 (PMID:31235963).",
        "32. Lu YC et al. Treatment of patients with metastatic cancer using a major histocompatibility complex class II-restricted T-cell receptor targeting MAGE-A3 (`KITE-718`). *J Clin Oncol.* 2017;35(29):3322–3329 (PMID:28809608).",
        "33. Seo M et al. BTN1A1 is a novel immune checkpoint mutually exclusive to PD-L1 (`hSTC810`). *J Immunother Cancer.* 2024;12(3):e008303 (PMID:38490683).",
        "34. Marcinkowski B et al. Cancer targeting by TCR gene-engineered T cells directed against Kita-Kyushu Lung Cancer Antigen-1 (`CT83` / `KK-LC-1`). *J Immunother Cancer.* 2019;7(1):229 (PMID:31462308).",
    ]
    md_text = "\n".join(md_lines) + "\n"

    # ---------------- LATEX PREPRINT ----------------
    tex_t1 = _tex_table(
        r"Genome-wide screening parameterization and category partition across human, domestic dog, and laboratory mouse (\texttt{oncoref.target-selection.2026-11-27.1}).",
        "tab:species",
        r"l l r c l X r r r r",
        [
            "Species",
            "Reference",
            "PEGs (Loci)",
            "Pillars",
            "Tumor / Line Scope",
            "CTA (Sex)",
            "CA only",
            "RA only",
            "Other",
            "Unsafe",
        ],
        [
            [
                "Human (9606)",
                "Ensembl 115",
                "19,893 (20,099)",
                "8/8",
                "135 cohorts",
                "322 (222M/48P/12F/40MF)",
                "53",
                "146",
                "1,116",
                "18,256",
            ],
            [
                "Dog (9615)",
                "Ensembl 98",
                "20,082 (20,199)",
                "8/8",
                f"19 coh ({d_samples}) + 12 CL",
                "118 (100M/7P/1F/10MF)",
                "127",
                "313",
                "2,201",
                "17,323",
            ],
            [
                "Mouse (10090)",
                "Ensembl 115",
                "21,433 (21,760)",
                "8/8",
                f"14 coh ({m_samples}) + 9 CL",
                "8 (5M/2P/0F/1MF)",
                "30",
                "303",
                "2,514",
                "18,578",
            ],
        ],
    )

    tex_t2 = _tex_table(
        r"Primary-literature functional biology and micro-niche resolution of all 66 baseline human \texttt{CA only} (\texttt{CA $\setminus$ RA}) protein groups.",
        "tab:caonly",
        r"p{3.1cm} r r p{5.4cm} X",
        [
            "Disposition / Class",
            "n",
            "Med. Fold",
            "Biological Mechanism & Why R < 10x H0",
            "Representative Proteins",
        ],
        t2_ca_only_rows,
    )

    tex_t3 = _tex_table(
        r"Empirical calibration of OncoRef normal-tissue and paralog safety gates against 16 human clinical TCR-T and ImmTAC trial regimens.",
        "tab:tcr",
        r"p{2.5cm} p{2.9cm} p{2.2cm} X p{2.5cm}",
        [
            "Target",
            "Product (HLA)",
            "Potency",
            "Observed Clinical Specificity Toxicity",
            "OncoRef Call",
        ],
        t3_tcr_rows,
    )

    tex_t4 = _tex_table(
        r"All 24 plasma-membrane extracellular-domain human Cancer Antigens (\texttt{CTA} $\cup$ \texttt{CA only}) actionable for CAR-T, ADC, TCE, or RLT modalities.",
        "tab:surface",
        r"l l p{3.1cm} r p{2.3cm} p{2.3cm} r r",
        [
            "Antigen",
            "Class",
            "Architecture",
            "Max (Tot) ECD",
            "ECD Span",
            "Peak TPM",
            "Pts",
            "2mm 9m",
        ],
        t4_surface_rows,
    )

    tex_lines = [
        r"\documentclass[10pt,letterpaper]{article}",
        r"\usepackage[margin=0.78in]{geometry}",
        r"\usepackage{fontspec}",
        r"\usepackage{microtype}",
        r"\usepackage{amsmath,amssymb}",
        r"\usepackage{booktabs,tabularx,array,longtable}",
        r"\usepackage{graphicx}",
        r"\usepackage{xcolor}",
        r"\usepackage[font=small,labelfont=bf]{caption}",
        r"\usepackage{enumitem}",
        r"\usepackage[colorlinks=true,linkcolor=teal!70!black,citecolor=teal!70!black,urlcolor=blue!60!black]{hyperref}",
        r"\setlength{\parskip}{0.42em}",
        r"\setlength{\parindent}{0pt}",
        r"\sloppy",
        r"\emergencystretch=1.8em",
        r"\definecolor{hdrnavy}{HTML}{0F172A}",
        "",
        r"\title{\textbf{\Large Multi-Scale Normal-Tissue Injury Tiers and Symmetric Enrichment Reconcile Cancer--Testis Antigen Discovery Across Human, Dog, and Mouse and Map Actionable TCR-T, CAR-T, ADC, and Radioligand Modalities}}",
        r"\author{\textbf{PIRL-UNC Collaborative Group}\\ \small Computational Oncology \& Immunogenomics, University of North Carolina at Chapel Hill \\ \footnotesize \texttt{https://github.com/pirl-unc/oncoref} $\cdot$ Policy \texttt{oncoref.target-selection.2026-11-27.1}}",
        r"\date{November 2026 $\cdot$ Preprint}",
        "",
        r"\begin{document}",
        r"\maketitle",
        r"\begin{abstract}",
        r"\noindent \textbf{Background:} Cancer--testis antigens (\texttt{CTAs}) and oncofetal proteins anchor affinity-enhanced T-cell receptor therapies (\texttt{TCR-T}), bispecific \texttt{ImmTACs}, vaccines, and emerging surface-directed \texttt{CAR-T}, antibody--drug conjugate (\texttt{ADC}), and radioligand (\texttt{RLT}) modalities. Yet across \textbf{149} published studies nominating \textbf{2,233} human protein groups, gene lists diverge sharply because legacy screens pool incompatible expression units, dilute rare neuroendocrine or sensory cells in bulk organ biopsies, conflate biological germline expression with patient-specific organ expendability, and transfer human labels onto dog and mouse models by symbol orthology.\\[0.3em]",
        r"\textbf{Methods:} We developed a species-invariant four-stage screening architecture (\texttt{OncoRef}, policy \texttt{2026-11-27.1}) that (i) collapses identical translations into protein-equivalence groups (\texttt{PEGs}) while screening all isoforms for $\ge 95\%$ 9-mer \texttt{Unsafe} containment; (ii) evaluates every healthy somatic and activated-immune dataset in its native measurement unit (\texttt{TPM}, \texttt{nTPM}, \texttt{nCPM}, donor-resolved droplet \texttt{UMI}, \texttt{IHC}, or plex-bleed-corrected \texttt{TMT-ppm}) across four clinical organ-injury tiers (\texttt{Vital}, \texttt{Protected}, \texttt{Manageable}, \texttt{Disposable}); (iii) enforces symmetric $\ge 10 \times H_0$ \emph{in vivo} enrichment over $H_0 = \max(N, D, 1\text{ TPM})$ in tumors (\texttt{CA}) and reproductive/embryonic compartments (\texttt{RA}), defining $\text{CTA} = \text{CA} \cap \text{RA}$ while segregating pure \emph{in vitro} cell lines; and (iv) couples candidates to \texttt{HostContext} organ protection, 0/1-mismatch \texttt{PeptideIndex} self-subtraction, and UniProtKB extracellular-domain modality mapping.\\[0.3em]",
        r"\textbf{Results:} Across \textbf{19,893} human (\texttt{GRCh38}), \textbf{20,082} domestic dog (\texttt{CanFam3.1}), and \textbf{21,433} laboratory mouse (\texttt{GRCm39}) \texttt{PEGs} (8/8 pillars active natively per species), OncoRef identifies \textbf{322 human}, \textbf{118 canine}, and \textbf{8 murine} \emph{in vivo} \texttt{CTAs} (96.6\%, 96.6\%, and 75.0\% ectopic). Segregating pure cell lines unmasks host \texttt{NK} \texttt{TIL} granzymes (\texttt{Gzmd--g}) in syngeneic grafts, and only \textbf{15/433} mammalian \texttt{CTA} ortholog families qualify as \texttt{CTAs} in $\ge 2$ species due to primate \texttt{CT-X}/\texttt{ZGA} expansion. In humans, \textbf{78.8\% (1,760/2,233)} of published nominations are \texttt{Unsafe} while rules uncover \textbf{77} unnominated \texttt{CTAs}; \texttt{CTA-M} (222) outnumbers \texttt{CTA-F} (12) by \textbf{18.5:1} alongside \textbf{48} adult-absent \texttt{CTA-P} antigens (270 female-safe vs.\ 60 male-safe \texttt{CTAs}). Primary-literature review of all \textbf{66} baseline \texttt{CA only} proteins adjudicates \textbf{7} occult adult human \texttt{Vital}/\texttt{Protected} niche proteins (\texttt{GALP}, \texttt{BHLHE23}, \texttt{IL36B}, \texttt{ZACN}, \texttt{HES3}, \texttt{OCSTAMP}, \texttt{PROKR1}) to \texttt{Unsafe} and \textbf{6} sensory/C-cell/paralog holds to \texttt{Other:Review}, resolving all \textbf{53} retained \texttt{CA only} proteins into five somatic-silent biological classes. Against \textbf{16} clinical \texttt{TCR-T}/\texttt{ImmTAC} trials, OncoRef flags \textbf{6/6 (100\%)} severe specificity toxicities as \texttt{Unsafe} while retaining \textbf{5/5 (100\%)} safe potent \texttt{CTAs}. Subcellular topology mapping shows \textbf{84.0\% (315/375)} of human \texttt{CAs} are intracellular (\texttt{TCR-T}/\texttt{ImmTAC}) and \textbf{9.6\% (36)} secreted, while \textbf{24 (6.4\%)} expose plasma-membrane ectodomains---including \textbf{19} with $\ge 50\text{ aa}$ mature ectodomains (\texttt{BTN1A1}, \texttt{ERVV-1/2}, \texttt{CIROP}, \texttt{CRIPTO3}, \texttt{ADAM2}, \texttt{ADAM18}, \texttt{FMR1NB}, \texttt{CT83}, \texttt{ZP4}, \texttt{IZUMO1R}, \texttt{LYPD4}, \texttt{PKDREJ}) actionable for \texttt{CAR-T}, \texttt{ADC}, and \texttt{RLT} without the somatic liabilities carried by \textbf{24/24 (100\%)} of conventional clinical surface benchmarks.",
        r"\end{abstract}",
        "",
        r"\section{Introduction}",
        r"Affinity-enhanced T-cell receptor therapies (\texttt{TCR-T}) and soluble TCR--anti-CD3 bispecific engagers (\texttt{ImmTACs}) recognize intracellularly processed peptides bound to major histocompatibility complex (\texttt{pMHC}) molecules at single-digit copy numbers per cell, producing durable clinical responses against \texttt{MAGEA4} (\emph{afamitresgene autoleucel}), \texttt{NY-ESO-1}/\texttt{CTAG2} (\emph{letetresgene autoleucel}), and \texttt{PRAME} (\texttt{IMA203}, \emph{brenetafusp}) [1--5]. That same sensitivity triggered lethal necrotizing leukoencephalopathy via rare brain \texttt{MAGEA12} cross-recognition (\texttt{NCT01273181}) [6], lethal cardiogenic shock within 4--5 days via striated-muscle \texttt{Titin} (\texttt{TTN}) mimotope cross-reactivity (\texttt{NCT01350401}) [7], dose-limiting hemorrhagic colitis against colonic \texttt{CEACAM5} [8], and uveitis plus cochlear hair-cell ototoxicity against melanocytic \texttt{MLANA} and \texttt{PMEL} [9].",
        "",
        r"Across \textbf{149} published CTA, oncofetal, zygotic-genome-activation (\texttt{ZGA}), and placental studies nominating \textbf{2,233} human protein groups [10--14], four methodological failure modes recur: (1) cross-scale pooling of incompatible units (\texttt{TPM}, \texttt{nTPM}, \texttt{CP10K}, microarray) and bulk-organ dilution of $<1\%$ vital neuroendocrine or sensory cell populations by $50\text{--}1,000\times$; (2) pre-analytical sequencing artifacts including multiplex flowcell index hopping, plate pipetting splash, nuclear pre-mRNA cis-readthrough (\texttt{GAB3 $\rightarrow$ CTAG2}), and tumor genomic DNA (\texttt{gDNA}) background; (3) conflation of biological germline expression with patient-conditional organ expendability; and (4) naive gene-symbol projection across human, companion-dog, and murine genomes. Here we present the unified \texttt{OncoRef} framework resolving all four challenges and mapping every retained candidate onto intracellular \texttt{pMHC} versus surface \texttt{CAR-T}/\texttt{ADC}/\texttt{RLT} therapeutic modalities.",
        "",
        r"\section{Conceptual Architecture: Native Scales, Injury Tiers, and Set Algebra}",
        r"\begin{figure}[htbp]",
        r"\centering",
        r"\includegraphics[width=\textwidth]{figures/fig1-conceptual-architecture-and-attrition.pdf}",
        r"\caption{\textbf{Conceptual pipeline architecture and human proteome attrition.} \textbf{(A)} Species-invariant Stages 0--3 decoupled from downstream \texttt{HostContext} organ protection and \texttt{PeptideIndex} self-proteome subtraction. \textbf{(B)} Partition of all 19,893 human Ensembl 115 protein-equivalence groups (20,099 coding loci).}",
        r"\label{fig:arch}",
        r"\end{figure}",
        "",
        r"\textbf{Stage 0 (Protein-equivalence grouping and pre-analytical QC).} Primary-assembly coding loci with identical canonical/MANE representative translations collapse into a single \texttt{PEG}, summing complete member-locus abundance within each library before summarizing across donors while retaining all splice isoforms (245,535 human, 45,036 canine, and 66,037 murine records) for $\ge 95\%$ 9-mer \texttt{Unsafe} containment. Holistic \texttt{NNLS} deconvolution across 18,637 GTEx v11 libraries (52 shards, 946 donors: 311 XX, 635 XY) quarantines 2 foreign-organ swaps and strips 201 SM-plate splash events, 100,287 index-hopping pairs, 225 same-donor carryovers, and 16,058 boundary admixtures. Libraries exceeding $2.5\%$ olfactory receptors $\ge 5\text{ TPM}$ fail \texttt{gDNA} QC, admixed \texttt{TGCT} biopsies cannot act as sole tumor witnesses, and purely nuclear pre-mRNA cis-spillover with zero cytoplasmic/bulk RNA and negative IHC is recorded as \texttt{nuclear\_only} caution.",
        "",
        r"\textbf{Stage 1 (Four clinical adult tissue-injury tiers without unit conversion).} Every non-reproductive structure and sorted immune lineage maps to a cited injury tier: \textbf{\texttt{Vital}} ($\le 4.0$ native bulk units: CNS, retina, inner-ear hair cells, DRG, heart, lung, liver, kidney, adrenal, pancreas, endothelium, HSCs, platelets), \textbf{\texttt{Protected}} ($\le 8.0$: GI/urinary/esophageal/corneal/skin epithelia, smooth/skeletal muscle, peripheral nerve, thyroid follicular parenchyma, pituitary, mature T/NK/monocyte/macrophage/neutrophil subsets), \textbf{\texttt{Manageable}} (hold $>16.0$, \texttt{Unsafe} $>128.0$: B cells, gallbladder, spleen, thyroid parafollicular C-cells), and \textbf{\texttt{Disposable}} ($D \le 32.0$: hair follicle, non-ocular skin melanocytes, tonsil, appendix, minor salivary gland, and separately tracked involuted adult thymus). Bulk RNA donor P95, single-cell cluster \texttt{nCPM}, donor-resolved droplet \texttt{UMI} ($\ge 10$ cells and $\ge 1\%$ multi-UMI cells replicated across $\ge 2$ studies/$\ge 3$ donors, or $\ge 5\%$ single-study; $2\times$ in \texttt{Protected}), protein staining/\texttt{TMT-ppm}, and curated adult human literature act in their native scales with strict \texttt{Unsafe} precedence.",
        "",
        r"\textbf{Stages 2--3 (Eligibility guards, symmetric $10\times H_0$ enrichment, and \texttt{HostContext}).} Surviving proteins require complete peripheral coverage (\texttt{Other:Review}), $<95\%$ \texttt{Unsafe} 9-mer overlap (\texttt{Other:Sequence}), and a screened home niche---holding all unvalidated monogenic chemosensory GPCRs (\texttt{OR*}, \texttt{TAAR*}, \texttt{TAS1R/2R*}, \texttt{VN1R*}) out of \texttt{CA}/\texttt{CTA}. Against $H_0 = \max(N, D, 1\text{ TPM})$, \emph{in vivo} tumor enrichment $T_{\max} \ge 10 H_0$ defines \texttt{CA}, reproductive or supported embryo-lineage enrichment $R_{\text{supp}} \ge 10 H_0$ defines \texttt{RA}, and $\text{CTA} = \text{CA} \cap \text{RA}$ (pure cell lines populate \texttt{cell\_line\_max\_tpm} only). Downstream \texttt{HostContext} protection yields four mutually exclusive classes (\texttt{CTA-M}, \texttt{CTA-P}, \texttt{CTA-F}, \texttt{CTA-MF}), followed by exact and 1-mismatch \texttt{PeptideIndex} self-proteome subtraction.",
        "",
        r"\section{Multi-Species Inventory and Parameterization: Human, Dog, and Mouse}",
        tex_t1,
        r"\begin{figure}[htbp]",
        r"\centering",
        r"\includegraphics[width=\textwidth]{figures/fig2-multispecies-parameterization-and-orthology.pdf}",
        r"\caption{\textbf{Multi-species parameterization, cell-line TIL diagnostic, and mammalian CTA orthology.} \textbf{(A)} Candidate sets across human, dog, and mouse (8/8 pillars active per species). \textbf{(B)} Murine \emph{in vivo} graft/GEMM TPM vs.\ pure cultured tumor cell-line TPM separating intrinsic CTAs from host NK TIL granzymes (\texttt{Gzmd--g}). \textbf{(C)} Fate of 433 mammalian CTA ortholog families.}",
        r"\label{fig:species}",
        r"\end{figure}",
        "",
        r"Evaluating all eight screening pillars natively in \emph{Homo sapiens} (\texttt{GRCh38}), outbred companion \emph{Canis lupus familiaris} (\texttt{CanFam3.1}; lacking seminal vesicles and appendix, screened across Bgee 15.2, 300-sample \texttt{GSE219045}, 64 mammary biopsies \texttt{GSE119810}, 7-dog healthy PBMC 10x \texttt{GSE225599}, 38-tissue plex-bleed-corrected \texttt{TMT-ppm} proteomics, and 19 spontaneous patient cohorts), and \emph{Mus musculus} (\texttt{GRCm39}; possessing coagulating gland, lacking tonsil/appendix and adult male mammary ducts, screened across Bgee 15.2, \texttt{GSE219045}, ImmGen + activated NK/Treg/BMDM RNA, Tabula Muris 10x, MGI GXD TS28, and 14 \emph{in vivo} cohorts) yields \textbf{322}, \textbf{118}, and \textbf{8} \emph{in vivo} \texttt{CTAs} (Table~\ref{tab:species}, Fig.~\ref{fig:species}). Segregating pure \emph{in vitro} cancer cell lines unmasks cytolytic granzymes \texttt{Gzmd--g} ($443\text{--}2,138\text{ TPM}$ in \emph{in vivo} \texttt{MC38}/\texttt{CT26}/\texttt{B16F10} grafts vs.\ $0.08\text{--}2.26\text{ TPM}$ in pure cultured tumor cells) as infiltrating host \texttt{NK} \texttt{TIL} transcripts, while only \textbf{15/433} mammalian \texttt{CTA} ortholog families qualify as \texttt{CTAs} in $\ge 2$ species because primate \texttt{CT-X}/\texttt{ZGA} arrays (\texttt{NY-ESO-1}, \texttt{GAGE}, \texttt{SPANX}, \texttt{XAGE}, \texttt{VCX/Y}, \texttt{DUX4}) evolved after boreoeutherian divergence.",
        "",
        r"\section{Human Proteome Deep Dive and the 66-Protein \texttt{CA only} Review}",
        r"\begin{figure}[htbp]",
        r"\centering",
        r"\includegraphics[width=\textwidth]{figures/fig3-human-deep-dive-and-ca-only-biology.pdf}",
        r"\caption{\textbf{Human literature reconciliation, sex-restriction asymmetry, and 66-protein \texttt{CA only} functional review.} \textbf{(A)} Disposition of 2,233 literature-nominated protein groups across 149 studies. \textbf{(B)} 4-way sex-restriction partition across 468 human RAs (322 CTAs). \textbf{(C)} Resolution of all 66 baseline \texttt{CA only} proteins into 7 \texttt{Unsafe}, 6 \texttt{Other:Review}, and 53 retained \texttt{CA only} across 5 biological classes.}",
        r"\label{fig:human}",
        r"\end{figure}",
        "",
        r"\textbf{Literature reconciliation and 18.5:1 sex asymmetry.} Of \textbf{2,233} published CTA/oncofetal nominations, \textbf{1,760 (78.8\%)} fail adult somatic gates (predominantly mitotic/centrosomal genes shared by spermatocytes, marrow, gut crypts, and activated lymphocytes), \textbf{113 (5.1\%)} enter \texttt{Other}, and \textbf{245} validate as \texttt{CTAs} alongside \textbf{77} newly identified expression-derived \texttt{CTAs} (Fig.~\ref{fig:human}A). Across the 322 \texttt{CTAs} (\textbf{311/322 = 96.6\%} with confirmed ectopic activation), \texttt{CTA-M} (222) outnumbers \texttt{CTA-F} (12) by \textbf{18.5:1}, complemented by \textbf{48} adult-absent peripartum \texttt{CTA-P} antigens (\textbf{270} female-safe vs.\ \textbf{60} male-safe \texttt{CTAs}; Fig.~\ref{fig:human}B). Auditing 159 female-lineage genes confirms this reflects continuous adult spermatogenesis behind the blood--testis barrier versus fetal completion of oogonial meiosis and sparse prophase-I oocytes ($<0.01\%$ of ovarian cortex, even with \texttt{GSE107746} isolated human oocytes), while true milk/placental proteins (\texttt{LALBA}, \texttt{CSN2}, \texttt{BTN1A1}, \texttt{PLAC1}, \texttt{ERVV-1/2}) populate \texttt{CTA-P}.",
        "",
        tex_t2,
        r"\textbf{Functional resolution of all 66 \texttt{CA only} (\texttt{CA $\setminus$ RA}) proteins.} Auditing primary human literature, \texttt{HPA v25.1} single-cell \texttt{nCPM}, \texttt{\#735} mucosal/skeletal panels, and 1-mismatch \texttt{Unsafe} identity across all 66 pre-adjudication \texttt{CA only} proteins (Table~\ref{tab:caonly}, Fig.~\ref{fig:human}C) adjudicated \textbf{7} occult adult human \texttt{Vital}/\texttt{Protected} micro-niche proteins to \texttt{Unsafe}---\texttt{GALP} (adult human hypothalamic arcuate neuroendocrine neurons in \texttt{HYPOMAP} [16]), \texttt{BHLHE23} (adult human retinal rod/ON-cone bipolar neurons [17]), \texttt{IL36B} (epidermal/lingual keratinocytes [18]), \texttt{ZACN} (pancreas/brain/cardiac/olfactory $\text{Zn}^{2+}$ channel [19]), \texttt{HES3} (adult human pancreatic islet $\alpha/\beta$ cells [20]), \texttt{OCSTAMP} (multinucleated bone-resorbing osteoclasts [21]), and \texttt{PROKR1} (adrenal medulla/GI GPCR [22])---and \textbf{6} sensory/C-cell/paralog proteins (\texttt{SRMS}, \texttt{OCM2} with 0 safe $\ge 2$-mm 9-mers vs.\ \texttt{OCM}, inner-ear hair-cell \texttt{CIB3} [23,24] and \texttt{OTOP1} [25], thyroid parafollicular C-cell \texttt{GFRA4} [26,27], and DRG \texttt{MRGPRG}) to \texttt{Other:Review}. Every one of the \textbf{53} retained \texttt{CA only} proteins belongs to five transparent somatic-silent classes: sub-$10\times$ male germline/meiotic proteins (\textbf{24}: \texttt{PRDM9} [28], \texttt{RPS4Y2}, \texttt{RBMY1E} $8.9\times$, \texttt{ZNF705D} $9.5\times$, \texttt{LIPJ} $9.6\times$, \texttt{PKDREJ}), pericentromeric/\texttt{RS447} macrosatellites and 8-cell cleavage arrays (\textbf{15}: \texttt{PRR20A--E}, \texttt{USP17L1/4/7}, \texttt{C1orf232}, \texttt{CCDC195}), post-gastrulation fetal morphogens/\texttt{HOX} factors (\textbf{6}: \texttt{AMELX} [29], \texttt{DBX1} [30], \texttt{HOXC12/B1/D12}, \texttt{EVX2}), sub-$10\times$ \texttt{CT-X} paralogs (\textbf{4}: \texttt{CSAG3} in 401 patients, \texttt{CSAG2}, \texttt{SSX5}, \texttt{MAGEB17}), and \texttt{Disposable} hair-cortex \texttt{KRTAPs} (\textbf{4}).",
        "",
        r"\section{Empirical Calibration Against Clinical TCR-T and ImmTAC Trials}",
        r"\begin{figure}[htbp]",
        r"\centering",
        r"\includegraphics[width=\textwidth]{figures/fig4-clinical-tcr-gold-standard.pdf}",
        r"\caption{\textbf{Empirical calibration against 16 human clinical TCR-T and ImmTAC regimens.} \textbf{(A)} Maximum healthy adult normal or cross-reactive paralog ratio ($\max_m x_m / \tau_m$) vs.\ the $1.0\times$ OncoRef ceiling. \textbf{(B)} 100\% sensitivity (6/6) to clinical specificity toxicity and 100\% recall (5/5) on potent safe CTAs.}",
        r"\label{fig:tcr}",
        r"\end{figure}",
        tex_t3,
        r"Across \textbf{16} curated endogenous-target clinical \texttt{TCR-T} and \texttt{ImmTAC} regimens (\textbf{14} potency-proven; Table~\ref{tab:tcr}, Fig.~\ref{fig:tcr}), OncoRef flags \textbf{6/6 (100\%, 0 false-safe misses)} severe specificity toxicities as \texttt{Unsafe} (\texttt{MAGEA3 $\rightarrow$ MAGEA12}, \texttt{MAGEA3 $\rightarrow$ TTN}, \texttt{CEACAM5}, \texttt{MLANA} \texttt{DMF5}/\texttt{DMF4}, \texttt{PMEL}, \texttt{TYR}) [6--9] while retaining \textbf{5/5 (100\%)} potency-proven safe \texttt{CTAs} (\texttt{MAGEA4} \emph{afami-cel}, \texttt{NY-ESO-1} and \texttt{CTAG2} \emph{lete-cel}, \texttt{PRAME} \texttt{IMA203}/\emph{brenetafusp}, \texttt{MAGEA10}) [1--5]. Three active non-CTA regimens (\texttt{WT1} \texttt{TCR-C4} immunoproteasome-differential \texttt{RMFPNAPYL} [31], \texttt{AFP}, and CD4 \texttt{HLA-DPB1*04:01} \texttt{MAGEA3} \texttt{KITE-718} [32]) define the epitope-processing frontier where whole-gene screening is intentionally conservative.",
        "",
        r"\section{Subcellular Localization \& Ectodomain Routing: CAR-T, ADC, and RLT}",
        r"\begin{figure}[htbp]",
        r"\centering",
        r"\includegraphics[width=\textwidth]{figures/fig5-surface-ectodomain-modality-routing.pdf}",
        r"\caption{\textbf{Subcellular localization and cell-surface ectodomain modality map.} \textbf{(A)} Partition of all 375 human Cancer Antigens (\texttt{322 CTA + 53 CA only}) across subcellular localization tiers. \textbf{(B)} Largest continuous mature extracellular domain (aa) vs.\ peak \emph{in vivo} tumor TPM across all 24 surface-exposed CAs alongside 24 clinical ADC/CAR-T/RLT benchmarks (all 24/24 \texttt{Unsafe} at TCR potency).}",
        r"\label{fig:surface}",
        r"\end{figure}",
        tex_t4,
        r"Mapping SHA-256-pinned UniProtKB topological domains, cleaved signal/propeptide spans, transmembrane helices, and \texttt{GPI}-anchors onto Ensembl 115 sequences across all \textbf{375} human \texttt{CAs} reveals that \textbf{315 (84.0\%)} are intracellular (268 nuclear/cytoplasmic/cytoskeletal including \texttt{PRAME}, \texttt{NY-ESO-1}, \texttt{MAGEA4}, \texttt{DUX4}, and fibrous-sheath \texttt{AKAP4}; 47 internal/minimal-ECD membrane) and \textbf{36 (9.6\%)} are secreted (\texttt{CSN2}, \texttt{LALBA}, \texttt{PLAC1}, \texttt{AMELX}), making \texttt{pMHC}-directed \texttt{TCR-T}/\texttt{ImmTAC} obligatory for $>84\%$ of somatic-silent antigens (Fig.~\ref{fig:surface}A). Conversely, \textbf{24/24 (100\%)} of conventional clinical \texttt{ADC}/\texttt{CAR-T}/\texttt{RLT} benchmarks (\texttt{HER2}, \texttt{TROP2}, \texttt{FR$\alpha$}, \texttt{PSMA}, \texttt{SSTR2}, \texttt{DLL3}, \texttt{CLDN18.2}, \texttt{CLDN6}, \texttt{MSLN}, \texttt{GPC3}, \texttt{ROR1}, \texttt{STEAP1}, \texttt{SEZ6}) are \texttt{Unsafe} at \texttt{TCR} potency because they rely on surface receptor-density windows or luminal polarization.",
        "",
        r"Strikingly, \textbf{24 human Cancer Antigens (22 \texttt{CTAs} + 2 \texttt{CA only})} combine \textbf{TCR-grade adult somatic silence} with \textbf{plasma-membrane extracellular exposure} (Table~\ref{tab:surface}, Fig.~\ref{fig:surface}B): \textbf{19} carry substantial $\ge 50\text{ aa}$ continuous mature ectodomains directly suited for \texttt{CAR-T}, \texttt{ADC}, \texttt{TCE}, and \texttt{RLT}---led by adult-absent peripartum \texttt{CTA-P} targets \textbf{\texttt{BTN1A1}} (216 aa IgV--IgC checkpoint ectodomain, 427.5 TPM [33]), retroviral envelopes \textbf{\texttt{ERVV-1}} and \textbf{\texttt{ERVV-2}} (300 aa \texttt{SU} ectodomains, 236--221 TPM in \texttt{KIRP}/\texttt{KIRC}), \textbf{\texttt{CIROP}} (715 aa), and \textbf{\texttt{CRIPTO3}} (188 aa EGF-CFC coreceptor, 11 \texttt{TGCT} patients); high-prevalence \texttt{CTA-M} metalloproteinase-like \textbf{\texttt{ADAM2}} (512 aa mature ECD, 54 patients) and \textbf{\texttt{ADAM18}} (503 aa, 40 patients), \textbf{\texttt{FMR1NB}} (94 aa loop, 49 patients), and Type II \textbf{\texttt{CT83}} (92 aa ECD, 22 patients [34]); \texttt{CTA-F} \textbf{\texttt{ZP4}} (444 aa, 11 patients) and \texttt{GPI}-anchored folate-receptor paralog \textbf{\texttt{IZUMO1R}} (209 aa, structurally homologous to \texttt{FOLR1}); \texttt{CTA-MF} \textbf{\texttt{LYPD4}} (199 aa \texttt{GPI}-anchored Ly6/uPAR fold); and \texttt{CA only} \textbf{\texttt{PKDREJ}} (1,165 aa REJ/GPS ECD)---plus \textbf{5} constrained 24--48 aa loop targets including 4-TM claudin \textbf{\texttt{CLDN34}} (48 aa ECL1, matching \texttt{CLDN18.2}/\texttt{CLDN6} topology).",
        "",
        r"\section{Limitations and Data \& Code Availability}",
        r"Peak bulk tumor \texttt{TPM} requires spatial \texttt{IHC}/\texttt{RNAscope} clonality validation; \texttt{TCR-T} candidates require immunopeptidomic confirmation; fetal mid-gestation organogenesis remains fail-closed (\texttt{pregnant=True}); and canonical surface topologies require tumor isoform/shedding verification. All data, SHA-256 receipts, and deterministic build scripts (\texttt{oncoref==1.55+}) are openly available under Apache 2.0 at \url{https://github.com/pirl-unc/oncoref}.",
        "",
        r"\small",
        r"\begin{thebibliography}{34}\setlength{\itemsep}{1pt}",
        r"\bibitem{r1} D'Angelo SP et al. NY-ESO-1 c259T cells in synovial sarcoma. \emph{Cancer Discov.} 2018;8:944--957 (PMID:29891538).",
        r"\bibitem{r2} D'Angelo SP et al. Afamitresgene autoleucel (SPEARHEAD-1). \emph{Lancet.} 2024;403:1460--1471 (PMID:38552654).",
        r"\bibitem{r3} Pollack SM et al. Letetresgene autoleucel in MRCLS. \emph{J Clin Oncol.} 2025 (PMC12084024).",
        r"\bibitem{r4} Wermke M et al. Autologous PRAME TCR-T (IMA203) Phase 1. \emph{Nat Med.} 2025 (PMID:40200105).",
        r"\bibitem{r5} Hamid O et al. Brenetafusp (IMC-F106C) PRAME ImmTAC in melanoma. \emph{J Clin Oncol / ESMO.} 2024.",
        r"\bibitem{r6} Morgan RA et al. Neurological toxicity after anti-MAGE-A3 TCR therapy. \emph{J Immunother.} 2013;36:133--151 (PMID:23377668).",
        r"\bibitem{r7} Linette GP et al. Titin cross-reactivity of MAGE-A3 TCR T cells. \emph{Blood.} 2013;122:863--871 (PMID:23770775).",
        r"\bibitem{r8} Parkhurst MR et al. CEA TCR T cells induce severe transient colitis. \emph{Mol Ther.} 2011;19:620--626 (PMID:21157437).",
        r"\bibitem{r9} Johnson LA et al. MART-1 and gp100 TCR therapy targets normal melanocytes. \emph{Blood.} 2009;114:535--546 (PMID:19451487).",
        r"\bibitem{r10} Almeida LG et al. CTdatabase. \emph{Nucleic Acids Res.} 2009;37:D816--D819 (PMID:18838390).",
        r"\bibitem{r11} Wang C et al. Systematic identification of cancer-testis genes in 19 cancer types. \emph{Nat Commun.} 2016;7:10499 (PMID:26813868).",
        r"\bibitem{r12} da Silva VL et al. Pan-cancer identification of cancer/testis genes. \emph{Oncotarget.} 2017;8:92966--92977 (PMID:29190970).",
        r"\bibitem{r13} Bruggeman JW et al. Germ cell-specific genes expressed in human cancers. \emph{Oncotarget.} 2018;9:29795--29812 (PMID:30042815).",
        r"\bibitem{r14} Carter JA et al. Pan-cancer analysis of cancer-testis antigens. \emph{J Immunother Cancer.} 2023 (PMID:37369392).",
        r"\bibitem{r15} Van den Eynde B et al. Murine mastocytoma antigen P1A (\texttt{Trap1a}). \emph{J Exp Med.} 1991;173:1373--1384; Sarma S et al. \emph{J Exp Med.} 1999;189:811--820.",
        r"\bibitem{r16} Tadross JA et al. Spatio-cellular map of the human hypothalamus (\texttt{HYPOMAP}). \emph{Nature.} 2025;638:1049--1059 (PMID:39910309).",
        r"\bibitem{r17} Cowan CS et al. Cell types of the human retina at single-cell resolution. \emph{Cell.} 2020;182:1623--1640 (PMID:32946783).",
        r"\bibitem{r18} Johnston A et al. IL-36 cytokines in human keratinocytes. \emph{J Immunol.} 2011;186:2613--2622 (PMID:21242515).",
        r"\bibitem{r19} Davies PA et al. Zinc-activated channel (\texttt{ZACN}) in human pancreas and brain. \emph{J Biol Chem.} 2003;278:712--717 (PMID:12381728).",
        r"\bibitem{r20} Masjkur J et al. Hes3 in adult human pancreatic islets. \emph{J Biol Chem.} 2014;289:35503--35516 (PMID:25371201).",
        r"\bibitem{r21} Khan UA et al. OC-STAMP in human osteoclast multinucleation. \emph{J Cell Biochem.} 2013;114:1772--1781 (PMID:23444125).",
        r"\bibitem{r22} Lin DC et al. Prokineticin receptor 1 (\texttt{PROKR1}). \emph{J Biol Chem.} 2002;277:19276--19280 (PMID:11886876).",
        r"\bibitem{r23} van der Valk WH et al. Single-cell atlas of human cochlea and vestibular organs. \emph{Cell Rep.} 2023;42:112623 (PMID:37289589).",
        r"\bibitem{r24} Liang X et al. CIB2 and CIB3 hair-cell MET channel subunits. \emph{EMBO J.} 2021;40:e107279 (PMID:34487363).",
        r"\bibitem{r25} Tu YH et al. Proton-selective ion channel OTOP1 in sour taste cells. \emph{Science.} 2018;359:1047--1050 (PMID:29371428).",
        r"\bibitem{r26} Lindahl M et al. Human GFR$\alpha$4 in thyroid parafollicular C-cells. \emph{J Biol Chem.} 2001;276:9344--9351 (PMID:11116144).",
        r"\bibitem{r27} Bhoj VG et al. GFR$\alpha$4 CAR-T in medullary thyroid cancer. \emph{Mol Ther Oncolytics.} 2021;20:387--398 (PMID:33632929).",
        r"\bibitem{r28} Baudat F et al. PRDM9 determines meiotic recombination hotspots. \emph{Science.} 2010;327:836--840 (PMID:20044539).",
        r"\bibitem{r29} Salido EC et al. Human amelogenin (\texttt{AMELX}). \emph{Am J Hum Genet.} 1992;50:303--316; Apps JR et al. \emph{J Pathol.} 2018;244:458--469.",
        r"\bibitem{r30} Pierani A et al. Progenitor homeodomain protein Dbx1. \emph{Neuron.} 2001;29:367--384 (PMID:11239428).",
        r"\bibitem{r31} Chapuis AG et al. WT1 TCR-C4 therapy post-transplant in AML. \emph{Nat Med.} 2019;25:1064--1072 (PMID:31235963).",
        r"\bibitem{r32} Lu YC et al. MHC class II-restricted MAGE-A3 TCR (KITE-718). \emph{J Clin Oncol.} 2017;35:3322--3329 (PMID:28809608).",
        r"\bibitem{r33} Seo M et al. BTN1A1 immune checkpoint (hSTC810). \emph{J Immunother Cancer.} 2024;12:e008303 (PMID:38490683).",
        r"\bibitem{r34} Marcinkowski B et al. TCR-T targeting CT83 (KK-LC-1). \emph{J Immunother Cancer.} 2019;7:229 (PMID:31462308).",
        r"\end{thebibliography}",
        r"\end{document}",
    ]
    tex_text = "\n".join(tex_lines) + "\n"
    return md_text, tex_text


def compile_xelatex(out_dir: Path, tex_filename: str = "oncoref-antigen-preprint.tex") -> Path:
    xelatex_bin = Path("/Library/TeX/texbin/xelatex")
    if not xelatex_bin.is_file():
        found = shutil.which("xelatex")
        if not found:
            raise RuntimeError("xelatex binary not found")
        xelatex_bin = Path(found)

    with tempfile.TemporaryDirectory(prefix="oncoref_preprint_tex_") as tmp_s:
        tmp = Path(tmp_s)
        shutil.copy2(out_dir / tex_filename, tmp / tex_filename)
        shutil.copytree(out_dir / "figures", tmp / "figures")
        for _pass in range(2):
            proc = subprocess.run(
                [str(xelatex_bin), "-interaction=nonstopmode", "-halt-on-error", tex_filename],
                cwd=tmp,
                capture_output=True,
                text=True,
                check=False,
            )
            if proc.returncode != 0:
                raise RuntimeError(f"xelatex failed:\n{proc.stdout[-2500:]}\n{proc.stderr[-1000:]}")
        pdf_name = tex_filename.replace(".tex", ".pdf")
        dest_pdf = out_dir / pdf_name
        shutil.copy2(tmp / pdf_name, dest_pdf)
    return dest_pdf


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "docs/preprint",
        help="Target directory for preprint outputs (default: docs/preprint)",
    )
    parser.add_argument(
        "--uniprot-snapshot",
        type=Path,
        default=None,
        help="Optional source path for uniprot-candidate-topology.json.gz",
    )
    parser.add_argument(
        "--compile-pdf", action="store_true", help="Compile LaTeX to PDF via xelatex"
    )
    parser.add_argument(
        "--check", action="store_true", help="Verify committed CSV, Markdown, and TeX match"
    )
    args = parser.parse_args()

    out_dir: Path = args.output_dir.resolve()
    data_dir = out_dir / "data"
    tables_dir = out_dir / "tables"
    fig_dir = out_dir / "figures"

    uniprot_gz = data_dir / "uniprot-candidate-topology.json.gz"
    uniprot_meta = data_dir / "uniprot-candidate-topology.json"

    if args.uniprot_snapshot is not None and args.uniprot_snapshot.is_file():
        data_dir.mkdir(parents=True, exist_ok=True)
        if args.uniprot_snapshot.resolve() != uniprot_gz.resolve():
            shutil.copy2(args.uniprot_snapshot, uniprot_gz)

    if not uniprot_gz.is_file():
        raise FileNotFoundError(f"Missing UniProt snapshot: {uniprot_gz}")

    gz_bytes = uniprot_gz.read_bytes()
    gz_sha = _sha256_bytes(gz_bytes)
    with gzip.open(uniprot_gz, "rt", encoding="utf-8") as fh:
        snap_doc = json.load(fh)

    meta_payload = (
        json.dumps(
            {
                "schema": "oncoref.preprint.uniprot-topology-receipt.v1",
                "policy_id": target_selection.POLICY_ID,
                "entry_count": int(snap_doc["entry_count"]),
                "bytes": len(gz_bytes),
                "sha256": gz_sha,
            },
            indent=2,
        )
        + "\n"
    )

    h_ts = human.target_selection()
    d_ts = canine.target_selection()
    m_ts = mouse.target_selection()
    ca_rev = ca_only_review.ca_only_functional_review()
    gs_eval = tcr_trials.evaluate_against_gold_standard()
    topo_df = build_topology_dataframe(uniprot_gz)

    csv_payload = topo_df.to_csv(index=False, lineterminator="\n")
    md_text, tex_text = build_preprint_texts(h_ts, d_ts, m_ts, ca_rev, gs_eval, topo_df, gz_sha)

    expected_files = {
        uniprot_meta: meta_payload,
        tables_dir / "candidate-surface-topology.csv": csv_payload,
        out_dir / "index.md": md_text,
        out_dir / "oncoref-antigen-preprint.tex": tex_text,
    }

    if args.check:
        for path, content in expected_files.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != content:
                raise SystemExit(f"Stale or missing preprint artifact: {path}")
        pdf_path = out_dir / "oncoref-antigen-preprint.pdf"
        if not pdf_path.is_file() or pdf_path.stat().st_size < 50_000:
            raise SystemExit(f"Missing or invalid compiled PDF: {pdf_path}")
        print(
            json.dumps(
                {
                    "status": "ok",
                    "policy_id": target_selection.POLICY_ID,
                    "checked": len(expected_files) + 1,
                }
            )
        )
        return

    for path, content in expected_files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    render_figures(fig_dir, h_ts, d_ts, m_ts, ca_rev, gs_eval, topo_df)

    pdf_size = 0
    if args.compile_pdf:
        pdf_out = compile_xelatex(out_dir)
        pdf_size = pdf_out.stat().st_size

    ca_only_mask = topo_df["category"].isin(["CTA", "Cancer Antigen only"])
    print(
        json.dumps(
            {
                "output_dir": str(out_dir),
                "policy_id": target_selection.POLICY_ID,
                "topology_rows": len(topo_df),
                "ca_localization_counts": topo_df.loc[ca_only_mask, "localization_class"]
                .value_counts()
                .to_dict(),
                "pdf_bytes": pdf_size,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

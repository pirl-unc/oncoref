#!/usr/bin/env python3
# Licensed under the Apache License, Version 2.0.
"""Render the live multi-species interactive antigen explorer (`docs/antigen-explorer.html`)."""

from __future__ import annotations

import argparse
import html
import json
import sys
from collections.abc import Iterable
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from oncoref import canine, human, mouse, target_selection, tcr_trials  # noqa: E402
from oncoref.species_target_selection import (  # noqa: E402
    classify_tumor_support,
    packaged_species_results,
)

OUT_HTML = ROOT / "docs/antigen-explorer.html"

HUMAN_ANCHORS = (
    "MAGEA3",
    "MAGEA12",
    "MLANA",
    "PMEL",
    "TYR",
    "CEACAM5",
    "WT1",
    "AFP",
    "DAZL",
    "GFY",
    "SP7",
    "PANX3",
    "CAGE1",
    "CTCFL",
    "DSPP",
    "TBATA",
    "OR8K5",
    "TAS1R2",
    "OR7C2",
    "KRT9",
    "PSMB11",
    "POTEC",
    "PSG8",
    "INSL4",
)

MOUSE_ANCHORS = (
    "Trap1a",
    "Gzmd",
    "Gzme",
    "Gzmf",
    "Gzmg",
    "Pmel",
    "Mlana",
    "Tyr",
    "Magea4",
    "Sp7",
    "Gfy",
)


def _num_or_none(val: object, digits: int = 2) -> float | None:
    if val is None or pd.isna(val) or val == "":
        return None
    return round(float(val), digits)


def _collect_nomination_symbols(candidates: Iterable[object]) -> set[str]:
    out: set[str] = set()
    for cand in candidates:
        if isinstance(cand, dict):
            for k in ("id", "label", "symbol", "reported_gene_symbol"):
                if cand.get(k):
                    out.add(str(cand[k]).strip())
            for tgt in cand.get("mapping", {}).get("targets", []) or []:
                if isinstance(tgt, dict) and tgt.get("symbol"):
                    out.add(str(tgt["symbol"]).strip())
        else:
            for k in ("symbol", "reported_gene_symbol", "id", "label"):
                v = getattr(cand, k, None)
                if v:
                    out.add(str(v).strip())
    return out


def _int_val(val: object) -> int:
    if val is None or pd.isna(val) or val == "":
        return 0
    return int(float(val))


def _species_payload(species_key: str, display_title: str) -> dict[str, object]:
    cat = packaged_species_results(species_key)
    tcr_by_gene: dict[str, str] = {}
    if species_key == "human":
        gs = tcr_trials.gold_standard()
        for r in gs.itertuples(index=False):
            tgt = str(r.target_symbol)
            note = (
                f"{r.product_name} ({r.hla_restriction}): "
                f"{r.clinical_stratum} — {r.observed_clinical_toxicities}"
            )
            tcr_by_gene[tgt] = f"{tcr_by_gene[tgt]}; {note}" if tgt in tcr_by_gene else note
        anchor_set = set(HUMAN_ANCHORS) | set(tcr_by_gene)
        n_cohorts = len(human.tumor_cohorts())
        n_lines = len(human.cell_lines())
    elif species_key == "dog":
        anchor_set = _collect_nomination_symbols(canine.cta_candidates())
        n_cohorts = len(canine.tumor_cohorts())
        n_lines = len(canine.cell_lines())
    else:
        anchor_set = set(MOUSE_ANCHORS) | _collect_nomination_symbols(mouse.cta_candidates())
        n_cohorts = len(mouse.tumor_cohorts())
        n_lines = len(mouse.cell_lines())

    anchor_eids: set[str] = set()
    nom_alias_by_eid: dict[str, list[str]] = {}
    for sym in sorted(anchor_set):
        matched = packaged_species_results(species_key, sym)
        for mrow in matched.itertuples(index=False):
            meid = str(mrow.entity_id)
            anchor_eids.add(meid)
            if str(mrow.protein_group_label) != sym:
                nom_alias_by_eid.setdefault(meid, []).append(sym)

    keep_mask = (
        cat["is_target_candidate"].astype(bool)
        | cat["entity_id"].isin(anchor_eids)
        | cat["protein_group_label"].isin(anchor_set)
        | (cat["published_nomination"].fillna(False).astype(bool) & (species_key != "human"))
    )
    sub = cat[keep_mask].copy()
    sub["_sort_cta"] = sub["is_cta"].astype(int)
    sub["_sort_cand"] = sub["is_target_candidate"].astype(int)
    sub["_sort_tum"] = pd.to_numeric(sub["tumor_max_tpm"], errors="coerce").fillna(-1.0)
    sub["_sort_rep"] = pd.to_numeric(sub["reproductive_supported_tpm"], errors="coerce").fillna(
        -1.0
    )
    sub = sub.sort_values(
        ["_sort_cta", "_sort_cand", "_sort_tum", "_sort_rep", "protein_group_label"],
        ascending=[False, False, False, False, True],
    )

    records: list[dict[str, object]] = []
    for r in sub.itertuples(index=False):
        raw_lbl = str(r.protein_group_label)
        aliases = nom_alias_by_eid.get(str(r.entity_id), [])
        lbl = f"{raw_lbl} ({', '.join(aliases)})" if aliases else raw_lbl
        is_ca = bool(r.is_cancer_antigen)
        s_10x = _int_val(getattr(r, "samples_at_tenfold", 0))
        p_10x = _int_val(getattr(r, "independent_patients_at_tenfold", 0))
        support = str(getattr(r, "tumor_support", "") or "") or classify_tumor_support(
            is_cancer_antigen=is_ca,
            samples_at_tenfold=s_10x,
            independent_patients_at_tenfold=p_10x,
        )
        mems = json.loads(str(r.member_gene_ids))
        reason = str(
            getattr(r, "decision_reason", "")
            or getattr(r, "Unsafe_reasons", "")
            or (
                "Passes normal-tissue safety and >=10x tumor & reproductive enrichment"
                if bool(r.is_cta)
                else str(r.classification)
            )
        )
        records.append(
            {
                "gene": lbl,
                "eid": str(r.entity_id),
                "genes": mems,
                "cat": str(r.category),
                "cls": str(r.classification),
                "sex": str(getattr(r, "sex_restriction_code", "") or ""),
                "mode": str(getattr(r, "tumor_activation_mode", "") or ""),
                "sup": support,
                "norm": _num_or_none(getattr(r, "normal_reference_tpm", None), 3),
                "disp": _num_or_none(getattr(r, "disposable_max_tpm", None), 3),
                "tum": _num_or_none(getattr(r, "tumor_max_tpm", None), 2),
                "tfold": _num_or_none(getattr(r, "tumor_fold", None), 1),
                "cohort": str(getattr(r, "tumor_top_cancer_code", "") or ""),
                "pat": p_10x,
                "samp": s_10x,
                "ctypes": _int_val(getattr(r, "cancer_types_at_tenfold", 0)),
                "rep": _num_or_none(getattr(r, "reproductive_supported_tpm", None), 2),
                "rfold": _num_or_none(getattr(r, "reproductive_fold", None), 1),
                "rcomp": str(getattr(r, "reproductive_qualifying_compartments", "") or ""),
                "cl": _num_or_none(getattr(r, "cell_line_max_tpm", None), 2),
                "cl_models": str(getattr(r, "qualifying_cell_lines", "") or ""),
                "home": str(getattr(r, "home_tissue", "") or ""),
                "u9": _int_val(getattr(r, "exact_unique_9mers", 0))
                if pd.notna(getattr(r, "exact_unique_9mers", None))
                else None,
                "s2mm": _int_val(getattr(r, "safe_2mm_9mers", 0))
                if pd.notna(getattr(r, "safe_2mm_9mers", None))
                else None,
                "mm1_pct": _num_or_none(getattr(r, "unsafe_near_self_1mm_pct", None), 1),
                "tcr": tcr_by_gene.get(raw_lbl, ""),
                "reason": reason,
            }
        )

    ctas = cat[cat["is_cta"].astype(bool)]
    return {
        "key": species_key,
        "title": display_title,
        "taxon_id": int(cat.attrs.get("taxon_id", 9606)),
        "policy_id": str(cat.attrs.get("policy_id", target_selection.POLICY_ID)),
        "total_proteins": len(cat),
        "cta_count": len(ctas),
        "ca_only_count": int(cat["classification"].eq("Cancer Antigen only").sum()),
        "ra_only_count": int(cat["classification"].eq("Reproductive Antigen only").sum()),
        "other_count": int(cat["category"].eq("Other").sum()),
        "unsafe_count": int(cat["category"].eq("Unsafe").sum()),
        "cta_m": int(ctas["sex_restriction_code"].eq("CTA-M").sum()),
        "cta_p": int(ctas["sex_restriction_code"].eq("CTA-P").sum()),
        "cta_f": int(ctas["sex_restriction_code"].eq("CTA-F").sum()),
        "cta_mf": int(ctas["sex_restriction_code"].eq("CTA-MF").sum()),
        "cohorts": n_cohorts,
        "cell_lines": n_lines,
        "rows": records,
    }


def render_html() -> str:
    data_bundle = {
        "human": _species_payload("human", "Human (Homo sapiens)"),
        "dog": _species_payload("dog", "Dog (Canis lupus familiaris)"),
        "mouse": _species_payload("mouse", "Mouse (Mus musculus)"),
    }
    embedded_json = json.dumps(data_bundle, separators=(",", ":"), ensure_ascii=True).replace(
        "</", "<\\/"
    )
    pol_esc = html.escape(target_selection.POLICY_ID)

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>OncoRef Interactive Antigen Explorer — Human, Dog &amp; Mouse</title>
<style>
  :root {{
    --bg: #f8fafc;
    --card: #ffffff;
    --ink: #0f172a;
    --muted: #475569;
    --line: #e2e8f0;
    --accent: #0f766e;
    --cta-bg: #dcfce7; --cta-fg: #14532d;
    --ca-bg: #dbeafe; --ca-fg: #1e3a8a;
    --ra-bg: #ede9fe; --ra-fg: #4c1d95;
    --oth-bg: #fef3c7; --oth-fg: #78350f;
    --uns-bg: #fee2e2; --uns-fg: #7f1d1d;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    background: var(--bg);
    color: var(--ink);
    line-height: 1.45;
  }}
  header {{
    background: #0f172a;
    color: #f8fafc;
    padding: 1.15rem 1.5rem;
    border-bottom: 3px solid #0d9488;
  }}
  header .top-row {{
    max-width: 1420px;
    margin: 0 auto;
    display: flex;
    flex-wrap: wrap;
    justify-content: space-between;
    align-items: baseline;
    gap: 0.75rem;
  }}
  header h1 {{
    margin: 0;
    font-size: 1.28rem;
    font-weight: 700;
  }}
  header .nav-links a {{
    color: #99f6e4;
    text-decoration: none;
    margin-left: 1rem;
    font-size: 0.88rem;
  }}
  header .nav-links a:hover {{ text-decoration: underline; }}
  main {{
    max-width: 1420px;
    margin: 1rem auto 2.5rem;
    padding: 0 1.25rem;
  }}
  .species-tabs {{
    display: flex;
    flex-wrap: wrap;
    gap: 0.5rem;
    margin-bottom: 0.85rem;
  }}
  .species-tab {{
    padding: 0.55rem 1rem;
    border: 1px solid var(--line);
    border-radius: 8px;
    background: var(--card);
    color: var(--ink);
    font-weight: 600;
    cursor: pointer;
    font-size: 0.92rem;
  }}
  .species-tab.active {{
    background: #0f766e;
    color: #ffffff;
    border-color: #0f766e;
  }}
  .kpi-strip {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
    gap: 0.6rem;
    margin-bottom: 0.9rem;
  }}
  .kpi {{
    background: var(--card);
    border: 1px solid var(--line);
    border-radius: 8px;
    padding: 0.6rem 0.8rem;
  }}
  .kpi .num {{ font-size: 1.22rem; font-weight: 700; }}
  .kpi .lbl {{ font-size: 0.78rem; color: var(--muted); }}
  .controls {{
    background: var(--card);
    border: 1px solid var(--line);
    border-radius: 10px;
    padding: 0.85rem 1rem;
    margin-bottom: 0.85rem;
  }}
  .pill-row {{
    display: flex;
    flex-wrap: wrap;
    gap: 0.4rem;
    margin-bottom: 0.65rem;
  }}
  .pill {{
    padding: 0.34rem 0.72rem;
    border-radius: 999px;
    border: 1px solid #cbd5e1;
    background: #f8fafc;
    color: #1e293b;
    font-size: 0.82rem;
    font-weight: 600;
    cursor: pointer;
  }}
  .pill.active {{
    background: #134e4a;
    color: #ffffff;
    border-color: #134e4a;
  }}
  .filter-bar {{
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 0.9rem;
  }}
  .search-box {{
    flex: 1 1 260px;
    padding: 0.45rem 0.7rem;
    border: 1px solid #cbd5e1;
    border-radius: 6px;
    font-size: 0.88rem;
  }}
  .toggle {{
    display: inline-flex;
    align-items: center;
    gap: 0.35rem;
    font-size: 0.84rem;
    color: #1e293b;
    cursor: pointer;
    user-select: none;
  }}
  .btn-export {{
    padding: 0.42rem 0.8rem;
    border-radius: 6px;
    border: 1px solid #0f766e;
    background: #ffffff;
    color: #0f766e;
    font-weight: 600;
    font-size: 0.82rem;
    cursor: pointer;
  }}
  .btn-export:hover {{ background: #f0fdfa; }}
  .explain-banner {{
    font-size: 0.84rem;
    color: var(--muted);
    margin: 0.45rem 0 0;
  }}
  .table-wrap {{
    background: var(--card);
    border: 1px solid var(--line);
    border-radius: 10px;
    overflow-x: auto;
  }}
  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.84rem;
  }}
  thead th {{
    background: #f1f5f9;
    text-align: left;
    padding: 0.58rem 0.65rem;
    border-bottom: 1px solid var(--line);
    cursor: pointer;
    white-space: nowrap;
    user-select: none;
    color: #334155;
  }}
  thead th:hover {{ background: #e2e8f0; }}
  tbody tr.data-row {{
    border-bottom: 1px solid #f1f5f9;
    cursor: pointer;
  }}
  tbody tr.data-row:hover {{ background: #f8fafc; }}
  tbody td {{
    padding: 0.5rem 0.65rem;
    vertical-align: middle;
  }}
  .num-col {{ text-align: right; font-variant-numeric: tabular-nums; }}
  .badge {{
    display: inline-block;
    padding: 0.14rem 0.48rem;
    border-radius: 5px;
    font-size: 0.75rem;
    font-weight: 700;
    white-space: nowrap;
  }}
  .b-cta {{ background: var(--cta-bg); color: var(--cta-fg); }}
  .b-ca  {{ background: var(--ca-bg);  color: var(--ca-fg); }}
  .b-ra  {{ background: var(--ra-bg);  color: var(--ra-fg); }}
  .b-oth {{ background: var(--oth-bg); color: var(--oth-fg); }}
  .b-uns {{ background: var(--uns-bg); color: var(--uns-fg); }}
  .sex-tag {{
    display: inline-block;
    padding: 0.1rem 0.42rem;
    border-radius: 4px;
    font-size: 0.74rem;
    font-weight: 600;
    background: #f1f5f9;
    color: #0f172a;
    border: 1px solid #cbd5e1;
  }}
  .mode-ect {{ color: #0f766e; font-weight: 600; }}
  .mode-lin {{ color: #64748b; }}
  tr.detail-row td {{
    background: #f8fafc;
    padding: 0.85rem 1.1rem;
    border-bottom: 2px solid #cbd5e1;
  }}
  .detail-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
    gap: 0.75rem;
    font-size: 0.82rem;
  }}
  .detail-box {{
    background: #ffffff;
    border: 1px solid var(--line);
    border-radius: 6px;
    padding: 0.6rem 0.75rem;
  }}
  .detail-box h4 {{
    margin: 0 0 0.35rem;
    font-size: 0.78rem;
    text-transform: uppercase;
    letter-spacing: 0.03em;
    color: var(--muted);
  }}
  code {{
    background: #f1f5f9;
    padding: 0.1rem 0.3rem;
    border-radius: 4px;
    font-size: 0.8rem;
  }}
</style>
</head>
<body>
<header>
  <div class="top-row">
    <div>
      <h1>OncoRef Interactive Antigen Explorer</h1>
      <div style="font-size:0.82rem;color:#cbd5e1;">
        Policy <code>{pol_esc}</code> · Symmetric <code>&ge;10 &times; max(Normal N, Disposable D, 1 TPM)</code> across Human, Dog &amp; Mouse
      </div>
    </div>
    <nav class="nav-links">
      <a href="index.html">Home</a>
      <a href="antigen-whitepaper/index.html">Antigen Whitepaper</a>
      <a href="species-antigen-whitepaper/index.html">Multi-Species Whitepaper</a>
      <a href="api/index.html">Python &amp; CLI API</a>
    </nav>
  </div>
</header>

<main>
  <div class="species-tabs" id="speciesTabs"></div>
  <div class="kpi-strip" id="kpiStrip"></div>

  <div class="controls">
    <div class="pill-row" id="presetPills"></div>
    <div class="filter-bar">
      <input type="search" id="searchInput" class="search-box"
             placeholder="Search gene symbol (MAGEA4, PRAME, DUX4, LALBA, Hmga2), Ensembl ID, cohort (OV, SKCM, OSA), or organ...">
      <label class="toggle"><input type="checkbox" id="chkEctopic"> Ectopic tumor activation only</label>
      <label class="toggle"><input type="checkbox" id="chkReplicated"> Replicated (&ge;2 patients/samples)</label>
      <label class="toggle"><input type="checkbox" id="chkCellLine"> Pure cell line &ge;10&times;</label>
      <button type="button" class="btn-export" id="btnExport">Export CSV</button>
    </div>
    <p class="explain-banner" id="presetExplain"></p>
  </div>

  <div style="font-size:0.82rem;color:#475569;margin-bottom:0.4rem;" id="resultCount"></div>

  <div class="table-wrap">
    <table>
      <thead>
        <tr>
          <th data-col="gene">Gene / Protein</th>
          <th data-col="cls">Category</th>
          <th data-col="sex">Physiology</th>
          <th data-col="mode">Activation</th>
          <th data-col="norm" class="num-col">Normal H&#8320;</th>
          <th data-col="tum" class="num-col">Tumor Peak TPM</th>
          <th data-col="tfold" class="num-col">Tumor Fold</th>
          <th data-col="cohort">Top Cohort</th>
          <th data-col="pat" class="num-col">Patients &ge;10&times;</th>
          <th data-col="rcomp">Reproductive / Home Niche</th>
          <th data-col="rep" class="num-col">Repro Peak TPM</th>
        </tr>
      </thead>
      <tbody id="tableBody"></tbody>
    </table>
  </div>
</main>

<script id="explorerData" type="application/json">{embedded_json}</script>
<script>
(function() {{
  const DATA = JSON.parse(document.getElementById("explorerData").textContent);
  const PRESETS = [
    {{ id: "CTA", label: "CTAs (CA ∩ RA)", desc: "Proteins passing somatic safety with >=10x enrichment in both in-vivo tumors and reproductive/embryonic tissues." }},
    {{ id: "FEMALE_SAFE", label: "Female-safe CTAs (CTA-M + CTA-P)", desc: "Physiological female-safe partition (CTA-M + CTA-P): CTAs silent (<8 TPM/nTPM, <16 IHC-) across non-pregnant, non-lactating adult female organs (antigens.table('CTA', sex='female') further enforces >=10x effective fold over protected female exposure)." }},
    {{ id: "MALE_SAFE", label: "Male-safe CTAs (CTA-F + CTA-P)", desc: "Physiological male-safe partition (CTA-F + CTA-P): CTAs silent (<8 TPM/nTPM, <16 IHC-) across adult male reproductive and resting mammary organs (antigens.table('CTA', sex='male') further enforces >=10x effective fold over protected male exposure)." }},
    {{ id: "CTA_P", label: "Adult-absent (CTA-P)", desc: "Pregnancy/peripartum-restricted CTAs (embryo, placenta, or lactating mammary only) — silent across all non-pregnant, non-lactating adults." }},
    {{ id: "CTA_F", label: "Female-restricted (CTA-F)", desc: "Oocyte, ovarian, uterine, or resting female mammary CTAs absent from adult male organs." }},
    {{ id: "CA_ONLY", label: "Cancer Antigen only", desc: "Somatic-safe proteins >=10x in in-vivo tumors without >=10x reproductive enrichment." }},
    {{ id: "RA_ONLY", label: "Reproductive Antigen only", desc: "Somatic-safe proteins >=10x in reproductive/embryonic organs without a qualifying in-vivo tumor sample yet." }},
    {{ id: "CANDIDATES", label: "All candidates (CA ∪ RA)", desc: "Complete union of CTA, Cancer Antigen only, and Reproductive Antigen only." }},
    {{ id: "ANCHORS", label: "Clinical & literature exclusions", desc: "Notable clinical TCR-T / ImmTAC benchmarks and historical literature nominations excluded as Unsafe or held in Other." }}
  ];

  let curSpecies = "human";
  let curPreset = "CTA";
  let sortCol = "tum";
  let sortDesc = true;
  let expandedKey = null;

  const tabsEl = document.getElementById("speciesTabs");
  const kpiEl = document.getElementById("kpiStrip");
  const pillsEl = document.getElementById("presetPills");
  const explainEl = document.getElementById("presetExplain");
  const countEl = document.getElementById("resultCount");
  const tbodyEl = document.getElementById("tableBody");
  const searchEl = document.getElementById("searchInput");
  const chkEct = document.getElementById("chkEctopic");
  const chkRep = document.getElementById("chkReplicated");
  const chkCl = document.getElementById("chkCellLine");

  function esc(s) {{
    return String(s ?? "").replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
  }}

  function fmtNum(v, d) {{
    if (v === null || v === undefined || Number.isNaN(v)) return "—";
    return Number(v).toLocaleString(undefined, {{ minimumFractionDigits: 0, maximumFractionDigits: d }});
  }}

  function badgeClass(cat) {{
    if (cat === "CTA") return "b-cta";
    if (cat === "Cancer Antigen only") return "b-ca";
    if (cat === "Reproductive Antigen only") return "b-ra";
    if (cat === "Unsafe") return "b-uns";
    return "b-oth";
  }}

  function matchesPreset(r) {{
    if (curPreset === "CTA") return r.cat === "CTA";
    if (curPreset === "FEMALE_SAFE") return r.sex === "CTA-M" || r.sex === "CTA-P";
    if (curPreset === "MALE_SAFE") return r.sex === "CTA-F" || r.sex === "CTA-P";
    if (curPreset === "CTA_P") return r.sex === "CTA-P";
    if (curPreset === "CTA_F") return r.sex === "CTA-F";
    if (curPreset === "CA_ONLY") return r.cls === "Cancer Antigen only";
    if (curPreset === "RA_ONLY") return r.cls === "Reproductive Antigen only";
    if (curPreset === "CANDIDATES") return r.cat === "CTA" || r.cls === "Cancer Antigen only" || r.cls === "Reproductive Antigen only";
    if (curPreset === "ANCHORS") return r.cat === "Unsafe" || r.cat === "Other";
    return true;
  }}

  function getFilteredRows() {{
    const sp = DATA[curSpecies];
    const q = searchEl.value.trim().toLowerCase();
    return sp.rows.filter(r => {{
      if (!matchesPreset(r)) return false;
      if (chkEct.checked && !(r.mode === "ectopic" || r.mode === "ectopic_and_lineage_retained")) return false;
      if (chkRep.checked && r.sup !== "replicated") return false;
      if (chkCl.checked && !r.cl_models) return false;
      if (q) {{
        const hay = (r.gene + " " + r.genes.join(" ") + " " + r.cohort + " " + r.rcomp + " " + r.sex + " " + r.cls + " " + r.reason + " " + r.tcr).toLowerCase();
        if (!hay.includes(q)) return false;
      }}
      return true;
    }}).sort((a, b) => {{
      let va = sortCol === "norm" ? (a.norm !== null ? Math.max(a.norm ?? 0, a.disp ?? 0, 1.0) : null) : a[sortCol];
      let vb = sortCol === "norm" ? (b.norm !== null ? Math.max(b.norm ?? 0, b.disp ?? 0, 1.0) : null) : b[sortCol];
      if (va === null || va === undefined || va === "") va = sortDesc ? -Infinity : Infinity;
      if (vb === null || vb === undefined || vb === "") vb = sortDesc ? -Infinity : Infinity;
      if (typeof va === "number" && typeof vb === "number") return sortDesc ? (vb - va) : (va - vb);
      return sortDesc ? String(vb).localeCompare(String(va)) : String(va).localeCompare(String(vb));
    }});
  }}

  function renderTabsAndKpis() {{
    tabsEl.innerHTML = Object.values(DATA).map(sp =>
      `<button type="button" class="species-tab ${{sp.key === curSpecies ? "active" : ""}}" data-sp="${{sp.key}}">` +
      `${{esc(sp.title)}} · ${{sp.cta_count}} CTAs</button>`
    ).join("");
    tabsEl.querySelectorAll("button").forEach(btn => {{
      btn.onclick = () => {{ curSpecies = btn.dataset.sp; expandedKey = null; renderAll(); }};
    }});

    const sp = DATA[curSpecies];
    const femSafe = sp.cta_m + sp.cta_p;
    const maleSafe = sp.cta_f + sp.cta_p;
    kpiEl.innerHTML = [
      [sp.total_proteins.toLocaleString(), "Coding protein groups"],
      [sp.cta_count.toLocaleString(), `CTAs (CA ∩ RA)`],
      [femSafe.toLocaleString(), `Female-safe (${{sp.cta_m}} M + ${{sp.cta_p}} P)`],
      [maleSafe.toLocaleString(), `Male-safe (${{sp.cta_f}} F + ${{sp.cta_p}} P)`],
      [sp.ca_only_count.toLocaleString(), "CA only"],
      [sp.ra_only_count.toLocaleString(), "RA only"],
      [`${{sp.cohorts}} / ${{sp.cell_lines}}`, "In-vivo cohorts / cell lines"]
    ].map(([n, l]) => `<div class="kpi"><div class="num">${{esc(n)}}</div><div class="lbl">${{esc(l)}}</div></div>`).join("");

    pillsEl.innerHTML = PRESETS.map(p =>
      `<button type="button" class="pill ${{p.id === curPreset ? "active" : ""}}" data-p="${{p.id}}">${{esc(p.label)}}</button>`
    ).join("");
    pillsEl.querySelectorAll("button").forEach(btn => {{
      btn.onclick = () => {{ curPreset = btn.dataset.p; expandedKey = null; renderAll(); }};
    }});
    const pObj = PRESETS.find(p => p.id === curPreset);
    explainEl.textContent = pObj ? pObj.desc : "";
  }}

  function renderTable() {{
    const rows = getFilteredRows();
    countEl.textContent = `Showing ${{rows.length.toLocaleString()}} proteins in ${{DATA[curSpecies].title}} (click any row for full evidence & Python snippet)`;
    const htmlRows = [];
    for (const r of rows) {{
      const h0 = Math.max(r.norm ?? 0, r.disp ?? 0, 1.0);
      const modeLbl = r.mode === "ectopic" ? "ectopic" : (r.mode === "ectopic_and_lineage_retained" ? "ectopic + lineage" : (r.mode === "lineage_retained" ? "lineage-retained" : "—"));
      const modeCls = (r.mode === "ectopic" || r.mode === "ectopic_and_lineage_retained") ? "mode-ect" : "mode-lin";
      const niche = r.rcomp || (r.home ? `home: ${{r.home}}` : "—");
      const patTxt = r.pat > 0 ? `${{r.pat}} (${{r.samp}}s)` : (r.samp > 0 ? `${{r.samp}}s` : "—");
      const rowKey = curSpecies + ":" + r.eid;
      htmlRows.push(
        `<tr class="data-row" data-rk="${{esc(rowKey)}}">` +
        `<td><strong>${{esc(r.gene)}}</strong></td>` +
        `<td><span class="badge ${{badgeClass(r.cat)}}">${{esc(r.cls)}}</span></td>` +
        `<td>${{r.sex ? `<span class="sex-tag">${{esc(r.sex)}}</span>` : "—"}}</td>` +
        `<td class="${{modeCls}}">${{esc(modeLbl)}}</td>` +
        `<td class="num-col">${{r.norm !== null ? fmtNum(h0, 2) : "—"}}</td>` +
        `<td class="num-col">${{fmtNum(r.tum, 1)}}</td>` +
        `<td class="num-col">${{r.tfold !== null ? fmtNum(r.tfold, 1) + "×" : "—"}}</td>` +
        `<td>${{esc(r.cohort || "—")}}</td>` +
        `<td class="num-col">${{esc(patTxt)}}</td>` +
        `<td>${{esc(niche)}}</td>` +
        `<td class="num-col">${{fmtNum(r.rep, 1)}}</td>` +
        `</tr>`
      );
      if (expandedKey === rowKey) {{
        const seqTxt = r.u9 !== null
          ? `Exact unique 9-mers: <strong>${{r.u9}}</strong> · &ge;2-mismatch safe 9-mers: <strong>${{r.s2mm}}</strong> · 1-mismatch near-self: <strong>${{fmtNum(r.mm1_pct, 1)}}%</strong>`
          : "Use <code>antigens.subtract(...)</code> for host self-proteome k-mer subtraction.";
        const clTxt = r.cl !== null
          ? `Peak cell-line TPM: <strong>${{fmtNum(r.cl, 1)}}</strong> (${{esc(r.cl_models || "below 10x")}})`
          : (DATA[curSpecies].cell_lines === 0 ? "Cell-line panel not evaluated in packaged human catalog." : "No in-vitro cell line >=10x recorded.");
        htmlRows.push(
          `<tr class="detail-row"><td colspan="11"><div class="detail-grid">` +
          `<div class="detail-box"><h4>Identity &amp; Decision</h4>` +
          `<div>Ensembl loci: <code>${{esc(r.genes.join(", "))}}</code></div>` +
          `<div style="margin-top:0.25rem;">Reason: ${{esc(r.reason)}}</div>` +
          `${{r.tcr ? `<div style="margin-top:0.3rem;color:#0f766e;"><strong>Clinical TCR benchmark:</strong> ${{esc(r.tcr)}}</div>` : ""}}</div>` +
          `<div class="detail-box"><h4>Baseline, Tumor &amp; Cell Lines</h4>` +
          `<div>Normal P95 N: <strong>${{fmtNum(r.norm, 3)}}</strong> · Disposable max D: <strong>${{fmtNum(r.disp, 3)}}</strong></div>` +
          `<div>Qualifying cancer types: <strong>${{r.ctypes}}</strong> · Replication: <strong>${{esc(r.sup || "none")}}</strong></div>` +
          `<div>${{clTxt}}</div></div>` +
          `<div class="detail-box"><h4>Sequence Safety &amp; Python API</h4>` +
          `<div>${{seqTxt}}</div>` +
          `<div style="margin-top:0.35rem;"><code>from oncoref import antigens; antigens.explain("${{esc(r.gene)}}", species="${{curSpecies}}")</code></div>` +
          `</div></div></td></tr>`
        );
      }}
    }}
    tbodyEl.innerHTML = htmlRows.join("");
    tbodyEl.querySelectorAll("tr.data-row").forEach(tr => {{
      tr.onclick = () => {{
        expandedKey = expandedKey === tr.dataset.rk ? null : tr.dataset.rk;
        renderTable();
      }};
    }});
  }}

  function renderAll() {{
    renderTabsAndKpis();
    renderTable();
  }}

  document.querySelectorAll("thead th[data-col]").forEach(th => {{
    th.onclick = () => {{
      const c = th.dataset.col;
      if (sortCol === c) sortDesc = !sortDesc;
      else {{ sortCol = c; sortDesc = (c !== "gene" && c !== "cls"); }}
      renderTable();
    }};
  }});
  searchEl.oninput = renderTable;
  chkEct.onchange = renderTable;
  chkRep.onchange = renderTable;
  chkCl.onchange = renderTable;

  document.getElementById("btnExport").onclick = () => {{
    const rows = getFilteredRows();
    const cols = ["gene", "cls", "sex", "mode", "sup", "norm", "disp", "tum", "tfold", "cohort", "pat", "samp", "rep", "rfold", "rcomp", "cl", "cl_models"];
    const lines = [cols.join(",")];
    for (const r of rows) {{
      lines.push(cols.map(k => `"${{String(r[k] ?? "").replace(/"/g, '""')}}"`).join(","));
    }}
    const blob = new Blob([lines.join("\\n") + "\\n"], {{ type: "text/csv;charset=utf-8" }});
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `oncoref-${{curSpecies}}-${{curPreset.toLowerCase()}}.csv`;
    a.click();
  }};

  renderAll();
}})();
</script>
</body>
</html>
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Verify docs/antigen-explorer.html is up to date without writing",
    )
    args = parser.parse_args()
    rendered = render_html()
    if args.check:
        if not OUT_HTML.exists() or OUT_HTML.read_text(encoding="utf-8") != rendered:
            print(
                "docs/antigen-explorer.html is out of date; run scripts/render_antigen_explorer.py"
            )
            return 1
        return 0
    OUT_HTML.write_text(rendered, encoding="utf-8")
    print(f"Wrote {OUT_HTML} ({len(rendered):,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

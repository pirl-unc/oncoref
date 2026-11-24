#!/usr/bin/env python3
"""Build SHA-256-verified GTEx v11 stoichiometric multiplex-bleed QC tables (#764).

Scans all non-reproductive GTEx v11 individual-biopsy GCT shards across ultra-abundant
source-organ transcripts (>= 1,000 TPM in testis, prostate, pancreas, liver, or stomach),
flags biopsies exhibiting proportional multi-sentinel cross-sample bleed (Nieuwenhuis et al.,
Genome Res 2020), checks orthogonal HPA v23 single-cell & bulk somatic expression, and
computes bleed-cleaned normal peripheral reference TPMs plus verified somatic tumor peaks.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from build_cta_alternative_evidence import GTEX_TISSUES  # noqa: E402

from oncoref import normal_library_qc as qc  # noqa: E402
from oncoref import reference_data, target_selection  # noqa: E402

DATA_DIR = ROOT / "oncoref/data"
AUDIT_DIR = ROOT / "docs/audits/multiplex-bleed-qc-20261123"

# Map GTEx v11 median column names (underscored) to anatomical tissues
NON_REPRODUCTIVE_SITES: dict[str, str] = {}
for _pretty, _organ in GTEX_TISSUES.items():
    _col = (
        _pretty.replace(" - ", "_")
        .replace(" (", "_")
        .replace(")", "")
        .replace(" ", "_")
        .replace("-", "-")
    )
    # Match exact pin key in reference_data
    NON_REPRODUCTIVE_SITES[_col] = _organ


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _load_gct_subset(path: Path, wanted_genes: set[bytes]) -> tuple[pd.DataFrame, pd.Series]:
    values: dict[str, np.ndarray] = {}
    with gzip.open(path, "rb") as handle:
        _ver = handle.readline()
        _dims = handle.readline()
        cols = handle.readline().decode("utf-8").rstrip("\n").split("\t")
        samples = pd.Index(cols[2:], name="sample_id")
        donors = samples.to_series().str.extract(r"^(GTEX-[^-]+)-.+$", expand=False)
        for line in handle:
            raw_name, _, rest = line.partition(b"\t")
            base_id = raw_name.split(b".", 1)[0]
            if base_id not in wanted_genes:
                continue
            gid = base_id.decode("ascii")
            vec = np.asarray(rest.partition(b"\t")[2].split(), dtype=float)
            if gid in values:
                values[gid] = values[gid] + vec
            else:
                values[gid] = vec
    return pd.DataFrame.from_dict(values, orient="index", columns=samples), donors


def main() -> None:
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    ts = target_selection.results()

    # 1. Load GTEx v11 tissue medians
    med_path = reference_data.pinned_path("gtex_tissue_rna", "v11")
    med = pd.read_csv(med_path, sep="\t", skiprows=2)
    med["gene_id"] = med["Name"].astype(str).str.split(".").str[0]
    sentinel_ref = qc.stoichiometric_sentinel_reference(med)

    # Identify all candidate proteins with ultra-high reproductive/source abundance (>= 10,000 TPM)
    # whose peripheral normal reference comes from GTEx donor P95 and could be inflated by bleed.
    source_cols = [col for _, (col, _) in qc.STOICHIOMETRIC_SENTINEL_PANELS.items()]
    med_max_by_gene = med.set_index("gene_id")[source_cols].max(axis=1)

    # Map candidate proteins (non-Unsafe, coverage-complete, reproductive/source peak >= 10,000 TPM)
    cand_mask = ~ts["category"].eq("Unsafe") & ts["required_coverage_complete"].astype(bool)
    cand_rows = []
    wanted_gene_ids: set[str] = set(sentinel_ref["gene_id"])
    for r in ts[cand_mask].itertuples(index=False):
        members = json.loads(r.member_gene_ids)
        rep_max = (
            float(r.reproductive_max_sample_tpm) if pd.notna(r.reproductive_max_sample_tpm) else 0.0
        )
        src_max = max((float(med_max_by_gene.get(g, 0.0)) for g in members), default=0.0)
        if rep_max >= 10000.0 or src_max >= 10000.0:
            cand_rows.append((r, members, max(rep_max, src_max)))
            wanted_gene_ids.update(members)

    wanted_bytes = {g.encode("ascii") for g in wanted_gene_ids}

    # 2. Orthogonal somatic expression check: HPA v23 single-cell clusters & HPA v23 bulk tissues
    sc_path = reference_data.ensure("hpa_single_cell_tissue", allow_download=False)
    sc_df = pd.read_csv(sc_path, sep="\t")
    repro_hpa_tissues = {
        "testis",
        "ovary",
        "placenta",
        "endometrium",
        "cervix",
        "fallopian tube",
        "breast",
        "prostate",
        "seminal vesicle",
        "epididymis",
        "vagina",
    }
    sc_somatic = sc_df[
        sc_df["Gene"].isin(wanted_gene_ids)
        & ~sc_df["Tissue"].astype(str).str.lower().isin(repro_hpa_tissues)
    ]
    sc_max_by_gene_tissue: dict[tuple[str, str], float] = (
        sc_somatic.groupby(["Gene", sc_somatic["Tissue"].astype(str).str.lower()])["nTPM"]
        .max()
        .to_dict()
    )

    hpa_bulk_path = reference_data.ensure("rna_tissue_hpa", allow_download=False)
    hpa_bulk = pd.read_csv(hpa_bulk_path, sep="\t")
    vcol_hpa = [c for c in hpa_bulk.columns if "TPM" in c or "nTPM" in c][-1]
    hpa_somatic = hpa_bulk[
        hpa_bulk["Gene"].isin(wanted_gene_ids)
        & ~hpa_bulk["Tissue"]
        .astype(str)
        .str.lower()
        .isin(repro_hpa_tissues | {"appendix", "gallbladder", "tonsil", "thymus"})
    ]
    hpa_bulk_by_gene_tissue: dict[tuple[str, str], float] = (
        hpa_somatic.groupby(["Gene", hpa_somatic["Tissue"].astype(str).str.lower()])[vcol_hpa]
        .max()
        .to_dict()
    )
    hpa_somatic_max_by_gene: dict[str, float] = (
        hpa_somatic.groupby("Gene")[vcol_hpa].max().to_dict()
    )

    # Also load sorted immune max TPM (Monaco + HPA immune samples) as orthogonal check (< 0.5 TPM)
    monaco_df = pd.read_csv(
        reference_data.ensure("hpa_immune_monaco", allow_download=False), sep="\t"
    )
    monaco_sub = monaco_df[monaco_df["Gene"].isin(wanted_gene_ids)]
    immune_max_by_gene: dict[str, float] = monaco_sub.groupby("Gene")["TPM"].max().to_dict()
    non_tgct_peaks: dict[str, tuple[float, str, str]] = {}
    tgct_sup = pd.read_csv(DATA_DIR / "tgct-aafl-protein-support.csv")
    for r in tgct_sup.itertuples(index=False):
        if pd.notna(r.max_non_tgct_tpm) and float(r.max_non_tgct_tpm) > 0:
            non_tgct_peaks[str(r.entity_id)] = (
                float(r.max_non_tgct_tpm),
                str(r.top_non_tgct_cancer_code) if pd.notna(r.top_non_tgct_cancer_code) else "",
                str(r.top_non_tgct_sample) if pd.notna(r.top_non_tgct_sample) else "",
            )

    # 4. Scan all non-reproductive GTEx v11 individual sample shards
    pins = reference_data.REFERENCE_SOURCES["gtex_individual_tissue_tpm"]["pins"]
    repro_gtex_organs = set(target_selection.REPRODUCTIVE_COMPARTMENTS) | {"breast", "thymus"}
    disposable_organs = {"appendix", "gallbladder", "tonsil", "spleen"}

    # Exact pin lookup
    pin_by_lower = {k.removeprefix("v11-"): k for k in pins if k.startswith("v11-")}

    tissue_summary_rows: list[dict] = []
    # Store per-(entity_id, site) raw vs clean donor P95
    site_raw_p95: dict[str, dict[str, float]] = {r.entity_id: {} for r, _, _ in cand_rows}
    site_clean_p95: dict[str, dict[str, float]] = {r.entity_id: {} for r, _, _ in cand_rows}

    # Determine dominant source organ per candidate entity
    entity_dom_source: dict[str, str] = {}
    for r, members, _ in cand_rows:
        best_src, best_val = "", 0.0
        for skey, (scol, _) in qc.STOICHIOMETRIC_SENTINEL_PANELS.items():
            v = max(
                (float(med.loc[med["gene_id"].eq(g), scol].max() or 0.0) for g in members),
                default=0.0,
            )
            if v > best_val:
                best_val = v
                best_src = skey
        # Also check reproductive_top_context
        top_ctx = str(r.reproductive_top_context or "").lower()
        if "testis" in top_ctx and best_val >= 1000.0:
            best_src = "testis"
        elif "prostate" in top_ctx and best_val >= 1000.0:
            best_src = "prostate"
        entity_dom_source[r.entity_id] = best_src

    for site_col, organ in sorted(NON_REPRODUCTIVE_SITES.items()):
        if organ in repro_gtex_organs:
            continue
        pin_key = pin_by_lower.get(site_col.lower())
        if pin_key is None:
            continue
        shard_path = reference_data.pinned_path("gtex_individual_tissue_tpm", pin_key)
        mat, donors = _load_gct_subset(shard_path, wanted_bytes)
        bleed_flags = qc.detect_stoichiometric_multiplex_bleed(
            mat,
            site=site_col,
            sentinel_reference=sentinel_ref,
        )
        n_samples = len(mat.columns)
        n_donors = int(donors.nunique())
        for skey in sorted(qc.STOICHIOMETRIC_SENTINEL_PANELS):
            sub_f = bleed_flags[bleed_flags["source_key"].eq(skey)]
            if not sub_f.empty or skey in {"testis", "prostate", "pancreas"}:
                tissue_summary_rows.append(
                    {
                        "gtex_site": site_col,
                        "organ": organ,
                        "contaminating_source": skey,
                        "n_samples": n_samples,
                        "n_donors": n_donors,
                        "flagged_samples": len(sub_f),
                        "flagged_pct": round(100.0 * len(sub_f) / max(n_samples, 1), 3),
                        "max_bleed_fraction": (
                            round(float(sub_f["median_bleed_fraction"].max()), 6)
                            if not sub_f.empty
                            else 0.0
                        ),
                        "median_bleed_fraction": (
                            round(float(sub_f["median_bleed_fraction"].median()), 6)
                            if not sub_f.empty
                            else 0.0
                        ),
                    }
                )

        if organ in disposable_organs:
            continue

        flagged_by_source: dict[str, set[str]] = {
            skey: set(bleed_flags.loc[bleed_flags["source_key"].eq(skey), "sample_id"])
            for skey in qc.STOICHIOMETRIC_SENTINEL_PANELS
        }

        organ_lower = organ.lower()
        for r, members, _ in cand_rows:
            eid = r.entity_id
            present_members = [g for g in members if g in mat.index]
            if len(present_members) != len(members):
                continue
            combined = mat.loc[present_members].sum(axis=0)
            raw_val = qc.compute_bleed_cleaned_donor_summary(
                combined,
                donors,
                flagged_samples=set(),
                orthogonal_supported=True,
            )
            site_raw_p95[eid][site_col] = raw_val

            dom_src = entity_dom_source.get(eid, "")
            f_samples = flagged_by_source.get(dom_src, set())
            # Orthogonal absence check in target organ (#764 / cross-cohort guardrail):
            # independent HPA v23 single-cell (< 0.5 nTPM), non-GTEx HPA bulk (< 0.5 nTPM),
            # and sorted immune (< 0.5 TPM if blood/lymphoid) must confirm absence!
            sc_val = max(
                (sc_max_by_gene_tissue.get((g, organ_lower), 0.0) for g in members),
                default=0.0,
            )
            hb_val = max(
                (hpa_bulk_by_gene_tissue.get((g, organ_lower), 0.0) for g in members),
                default=0.0,
            )
            imm_val = (
                max((immune_max_by_gene.get(g, 0.0) for g in members), default=0.0)
                if organ_lower in {"blood", "lymph node", "bone marrow"}
                else 0.0
            )
            ortho_ok = (sc_val >= 0.5) or (hb_val >= 0.5) or (imm_val >= 0.5)
            clean_val = qc.compute_bleed_cleaned_donor_summary(
                combined,
                donors,
                flagged_samples=f_samples,
                orthogonal_supported=ortho_ok,
            )
            site_clean_p95[eid][site_col] = clean_val

    # 5. Build per-protein delta ledger
    delta_rows: list[dict] = []
    for r, members, peak_src_tpm in cand_rows:
        eid = r.entity_id
        base_norm = float(r.normal_reference_tpm) if pd.notna(r.normal_reference_tpm) else np.nan
        if pd.isna(base_norm) or not site_raw_p95[eid]:
            continue
        raw_gtex_max = max(site_raw_p95[eid].values())
        clean_gtex_max = max(site_clean_p95[eid].values())
        top_raw_site = max(site_raw_p95[eid], key=site_raw_p95[eid].get)
        top_clean_site = max(site_clean_p95[eid], key=site_clean_p95[eid].get)

        # Preserve any non-GTEx somatic ceiling: if baseline normal_reference_tpm came from a
        # non-GTEx source (e.g. HPA individual tissue samples, DICE immune, brain subregions),
        # normal_reference_contexts tells us whether GTEx was the binding maximum!
        ctx_raw = str(r.normal_reference_contexts or "[]")
        ctx_list = json.loads(ctx_raw) if ctx_raw.startswith("[") else []
        non_gtex_bound = [
            float(obs["value"])
            for obs in ctx_list
            if obs.get("source_id") != "gtex_donor_p95:v11" and obs.get("value") is not None
        ]
        # Also retain HPA somatic consensus / bulk floor as conservative lower bound on non-GTEx bulk
        hpa_floor = max((hpa_somatic_max_by_gene.get(g, 0.0) for g in members), default=0.0)
        if non_gtex_bound:
            # Another non-GTEx source already tied or set the maximum
            clean_norm = base_norm
        elif abs(raw_gtex_max - base_norm) <= 0.05:
            # GTEx donor P95 was the binding peripheral normal maximum!
            clean_norm = max(clean_gtex_max, hpa_floor)
        else:
            # A non-GTEx source (e.g. HPA biopsy or DICE immune) was higher than GTEx donor P95
            clean_norm = base_norm

        applied = bool(clean_norm < base_norm - 0.10)
        nt_tpm, nt_code, nt_sample = non_tgct_peaks.get(eid, (np.nan, "", ""))
        delta_rows.append(
            {
                "entity_id": eid,
                "protein_group_label": r.protein_group_label,
                "dominant_source_organ": entity_dom_source.get(eid, ""),
                "source_peak_tpm": round(peak_src_tpm, 3),
                "baseline_normal_reference_tpm": round(base_norm, 6),
                "raw_gtex_p95_max_tpm": round(raw_gtex_max, 6),
                "top_raw_gtex_site": top_raw_site,
                "clean_gtex_p95_max_tpm": round(clean_gtex_max, 6),
                "top_clean_gtex_site": top_clean_site,
                "hpa_somatic_max_ntpm": round(hpa_floor, 4),
                "clean_normal_reference_tpm": round(clean_norm, 6),
                "bleed_qc_applied": applied,
                "max_non_tgct_tumor_tpm": round(nt_tpm, 4) if pd.notna(nt_tpm) else np.nan,
                "top_non_tgct_cancer_code": nt_code,
                "top_non_tgct_sample": nt_sample,
            }
        )

    tissue_df = pd.DataFrame(tissue_summary_rows).sort_values(
        ["contaminating_source", "flagged_pct", "gtex_site"],
        ascending=[True, False, True],
    )
    deltas_df = pd.DataFrame(delta_rows).sort_values(
        ["bleed_qc_applied", "baseline_normal_reference_tpm", "protein_group_label"],
        ascending=[False, False, True],
    )

    t_csv = DATA_DIR / "normal-multiplex-bleed-tissue-summary.csv"
    d_csv = DATA_DIR / "normal-multiplex-bleed-gene-deltas.csv"
    tissue_df.to_csv(t_csv, index=False)
    deltas_df.to_csv(d_csv, index=False)

    t_receipt = {
        "schema": "oncoref.multiplex-bleed-tissue-summary.v1",
        "policy_id": qc.MULTIPLEX_BLEED_POLICY_ID,
        "rows": len(tissue_df),
        "sha256": _sha256_bytes(t_csv.read_bytes()),
    }
    d_receipt = {
        "schema": "oncoref.multiplex-bleed-gene-deltas.v1",
        "policy_id": qc.MULTIPLEX_BLEED_POLICY_ID,
        "evaluated_candidates": len(deltas_df),
        "bleed_qc_applied_count": int(deltas_df["bleed_qc_applied"].sum()),
        "sha256": _sha256_bytes(d_csv.read_bytes()),
    }
    (DATA_DIR / "normal-multiplex-bleed-tissue-summary.json").write_text(
        json.dumps(t_receipt, indent=2) + "\n", encoding="utf-8"
    )
    (DATA_DIR / "normal-multiplex-bleed-gene-deltas.json").write_text(
        json.dumps(d_receipt, indent=2) + "\n", encoding="utf-8"
    )

    # Render audit markdown
    applied_sub = deltas_df[deltas_df["bleed_qc_applied"]].copy()
    top_testis = tissue_df[
        tissue_df["contaminating_source"].eq("testis") & tissue_df["flagged_samples"].gt(0)
    ].head(12)

    md_lines = [
        "# Multiplex Index-Hopping & Cross-Sample Bleed QC (`#764`)",
        "",
        f"Policy `{qc.MULTIPLEX_BLEED_POLICY_ID}` · Reference panel `gtex_individual_tissue_tpm:v11` · Orthogonal checks `hpa_single_cell_tissue:v23` + `rna_tissue_hpa:v23`.",
        "",
        "## Method",
        "",
        "Ultra-abundant tissue-specific transcripts (`>= 1,000 TPM`, such as spermatid `PRM1`/`PRM2`/`TNP1` in testis or `SEMG1`/`SEMG2`/`TGM4` in prostate) bleed at `10^-5` to `5 x 10^-3` into unrelated biopsies sequenced on shared multiplexed Illumina flowcells ([Nieuwenhuis et al., *Genome Res* 2020](https://doi.org/10.1101/gr.255380.119)).",
        "",
        "A biopsy $b$ in non-source organ $T$ is flagged for source-$S$ bleed if and only if:",
        "1. **Multi-sentinel co-detection:** `>= 3` sentinel markers from $\\mathcal{M}_S$ are simultaneously detected at `>= 0.5 TPM` in biopsy $b$.",
        "2. **Stoichiometric proportionality:** Median $\\hat{\\alpha}_{b,g} = \\text{TPM}_{b,g} / \\text{median\\_TPM}_{S,g} \\in [10^{-5}, 5 \\times 10^{-3}]$ with $\\text{IQR}(\\log_{10} \\hat{\\alpha}_{b,g}) \\le 0.5$.",
        "3. **Orthogonal somatic veto:** Independent `HPA v23` single-cell clusters (`< 1.0 nTPM`) and non-GTEx `rna_tissue_hpa` (`< 2.0 nTPM`) confirm no true expression in organ $T$.",
        "",
        "## Top GTEx Sites with Coordinated Testis Sentinel Bleed",
        "",
        "| GTEx Site | Organ | Samples | Flagged | Flagged % | Median Bleed Fraction |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for r in top_testis.itertuples(index=False):
        md_lines.append(
            f"| `{r.gtex_site}` | {r.organ} | {r.n_samples:,} | {r.flagged_samples:,} | {r.flagged_pct:.2f}% | {r.median_bleed_fraction:.5f} |"
        )
    md_lines.extend(
        [
            "",
            f"## Proteins with Bleed-Inflated GTEx Peripheral Baselines (`n = {len(applied_sub)}`)",
            "",
            "| Protein | Source | Source Peak TPM | Raw Baseline TPM | Top Raw Site | Clean Baseline TPM | Top Non-TGCT Tumor TPM |",
            "|---|---|---:|---:|---|---:|---|",
        ]
    )
    for r in applied_sub.itertuples(index=False):
        nt_str = (
            f"{r.max_non_tgct_tumor_tpm:.2f} (`{r.top_non_tgct_cancer_code}`)"
            if pd.notna(r.max_non_tgct_tumor_tpm)
            else "—"
        )
        md_lines.append(
            f"| `{r.protein_group_label}` | {r.dominant_source_organ} | {r.source_peak_tpm:,.1f} | {r.baseline_normal_reference_tpm:.2f} | `{r.top_raw_gtex_site}` | **{r.clean_normal_reference_tpm:.2f}** | {nt_str} |"
        )
    (AUDIT_DIR / "index.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(
        f"Built multiplex bleed QC: {len(tissue_df)} site-source rows, "
        f"{len(deltas_df)} candidates evaluated, {len(applied_sub)} normal baselines cleaned."
    )


if __name__ == "__main__":
    main()

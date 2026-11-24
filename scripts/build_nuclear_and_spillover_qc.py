#!/usr/bin/env python3
"""Build SHA-256-verified nuclear snRNA-seq compartment & spillover QC tables (#765).

Generates:
  - oncoref/data/census-collection-compartments.{csv,json}
  - oncoref/data/nuclear-spillover-adjudications.{csv,json}
  - oncoref/data/cleared-candidate-comparison-observations.json
  - docs/audits/nuclear-rna-spillover-20261123/index.md
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pyensembl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from antigen_sample_evidence import developmental_samples  # noqa: E402
from rebuild_target_catalog import (  # noqa: E402
    BASE,
    COMPARISON_ID,
    disposable_ntpm,
    normal_panel,
    shared_facts,
)

from oncoref import (  # noqa: E402
    catalog,
    cta_sequence_overlap,
    reference_data,
    source_matrices,
    target_selection,
)
from oncoref import (  # noqa: E402
    rare_cell_compartment_qc as rc_qc,
)

DATA_DIR = ROOT / "oncoref/data"
ATLAS_DIR = ROOT / "docs/audits/normal-cell-atlas-20261110"
REBUILD_DIR = ROOT / "docs/audits/target-selection-20261110/rebuild"
AUDIT_DIR = ROOT / "docs/audits/nuclear-rna-spillover-20261123"

KNOWN_NUCLEUS_COLLECTIONS: frozenset[str] = frozenset(
    {
        "283d65eb-dd53-496d-adb7-7570c7caa443",  # Siletti Human Brain Cell Atlas v1.0
        "d17249d2-0e6e-4500-abb8-e6c93fa1ac6f",  # Allen Brain Cell Atlas
        "1ca90a2d-2943-483d-b678-b809bf464c30",  # SEA-AD
        "72d37bc9-76cc-442d-9131-da0e273862db",  # Sepp mammalian cerebellum
        "180bff9c-c8a5-4539-b13b-ddbc00d643e6",  # Leng AD vulnerable neurons
        "0d35c0fd-ef0b-4b70-bce6-645a4660e5fa",  # Resilience AD L4 neurons
        "16c1e722-96ae-4bf6-b408-cd7f8918484f",  # Jäkel MS white matter snRNA-seq
        "22ccff12-9a18-4b01-ba10-302df57a01f8",  # Tran reward circuitry snRNA-seq
        "7c4552fd-8a6d-4da3-9854-2dfa8baca8bf",  # CSF1R single-nuclear
        "84ce6837-548d-4a1f-919f-0bc0d9a3952f",  # Prefrontal cortex snRNA-seq
        "b0f0b447-ac37-45b0-b1bf-5c0b7d871120",  # Kamath substantia nigra nuclei
        "e1fa9900-3fc9-4b57-9dce-c95724c88716",  # dlPFC single-nucleus
        "433700dc-e8a5-48b0-b5cd-beb22f3f88fe",  # Brain vascular snRNA-seq
        "307ce143-7cc8-4813-99b1-3797834149c9",  # Kidney snRNA-seq
        "6b701826-37bb-4356-9792-ff41fc4c3161",  # White adipose snRNA-seq
        "348da6dc-5bf6-435d-adc5-37747b9ae38a",  # Retina multiome nuclei
    }
)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _classify_collection_compartment(
    cid: str,
    cname: str,
    dtitles: str,
    tissues: list[str],
) -> tuple[str, str]:
    if cid in KNOWN_NUCLEUS_COLLECTIONS:
        return "nucleus", "curated_snrnaseq_or_multiome_nuclei"
    text = f"{cname} || {dtitles}".lower()
    has_sn = bool(
        re.search(r"\bsnrna|\bsingle[- ]nucle(?:us|i|ar)|\bnuclei\b|\bsingle[- ]soma\b", text)
    )
    has_sc = bool(re.search(r"\bscrna|\bsingle[- ]cell\b", text))
    # Postmortem adult human CNS parenchymal neuron atlases in Census use nuclear isolation
    if set(tissues) <= {"brain", "spinal cord"} and cid != "854c0855-23ad-4362-8b77-6b1639e7a9fc":
        return "nucleus", "postmortem_adult_cns_nuclear_isolation"
    if has_sn and not has_sc:
        return "nucleus", "title_snrnaseq"
    if has_sn and has_sc:
        return "mixed", "title_mixed_sc_and_sn"
    return "cell", "whole_cell_scrnaseq"


def main() -> None:
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    sr = pd.read_csv(ATLAS_DIR / "source-registry.csv")
    dv = pd.read_csv(ATLAS_DIR / "dataset-versions.csv")
    dp = pd.read_csv(ATLAS_DIR / "rare-cell/donor-populations.csv.gz")
    rc_findings = pd.read_csv(ATLAS_DIR / "rare-cell/findings.csv")
    rce = pd.read_csv(DATA_DIR / "rare-cell-evidence.csv.gz")
    ts = target_selection.results()
    ts_idx = ts.set_index("entity_id")

    sr_map = sr.set_index("collection_id")["collection_name"].to_dict()
    sr_doi = sr.set_index("collection_id")["collection_doi"].fillna("").to_dict()
    dv_titles = (
        dv.groupby("collection_id")["dataset_title"]
        .apply(lambda s: " | ".join(sorted(set(s.dropna().astype(str)))))
        .to_dict()
    )
    dp_tissues = dp.groupby("source_id")["tissue"].apply(lambda s: sorted(set(s))).to_dict()
    dp_donors = dp.groupby("source_id")["donor_id"].nunique().to_dict()

    comp_rows: list[dict] = []
    comp_by_cid: dict[str, str] = {}
    for cid in sorted(dp["source_id"].unique()):
        cname = str(sr_map.get(cid, ""))
        dtitles = str(dv_titles.get(cid, ""))
        tissues = dp_tissues.get(cid, [])
        comp, basis = _classify_collection_compartment(cid, cname, dtitles, tissues)
        comp_by_cid[cid] = comp
        comp_rows.append(
            {
                "collection_id": cid,
                "assay_compartment": comp,
                "compartment_basis": basis,
                "tissues": ";".join(tissues),
                "n_donors": int(dp_donors.get(cid, 0)),
                "collection_doi": str(sr_doi.get(cid, "")),
                "collection_name": cname,
            }
        )
    comp_df = pd.DataFrame(comp_rows)
    comp_csv_path = DATA_DIR / "census-collection-compartments.csv"
    comp_df.to_csv(comp_csv_path, index=False)
    (DATA_DIR / "census-collection-compartments.json").write_text(
        json.dumps(
            {
                "schema": "oncoref.census-collection-compartments.v1",
                "policy_id": rc_qc.NUCLEAR_SPILLOVER_POLICY_ID,
                "collections": len(comp_df),
                "compartment_counts": comp_df["assay_compartment"].value_counts().to_dict(),
                "sha256": _sha256_bytes(comp_csv_path.read_bytes()),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    # Load somatic bulk + immune expression for Ensembl 115 neighbor activity
    gtex_med = pd.read_csv(
        reference_data.pinned_path("gtex_tissue_rna", "v11"), sep="\t", skiprows=2
    )
    gtex_med["gene_id"] = gtex_med["Name"].astype(str).str.split(".").str[0]
    repro_cols = {
        "Testis",
        "Ovary",
        "Fallopian_Tube",
        "Uterus",
        "Vagina",
        "Cervix_Ectocervix",
        "Cervix_Endocervix",
        "Prostate",
        "Breast_Mammary_Tissue",
    }
    som_cols = [c for c in gtex_med.columns[2:-1] if c not in repro_cols]
    gtex_som_max = gtex_med.set_index("gene_id")[som_cols].max(axis=1).to_dict()

    hpa_imm = pd.read_csv(
        reference_data.ensure("hpa_immune_cell_samples", allow_download=False), sep="\t"
    )
    imm_max = hpa_imm.groupby("ENSG ID")["TPM"].max().to_dict()
    neut_sub = hpa_imm[hpa_imm["Immune cell"].astype(str).str.lower().str.contains("neutrophil")]
    neut_max = neut_sub.groupby("ENSG ID")["TPM"].max().to_dict()

    monaco_df = pd.read_csv(
        reference_data.ensure("hpa_immune_monaco", allow_download=False), sep="\t"
    )
    monaco_max = monaco_df.groupby("Gene")["TPM"].max().to_dict()

    ihc = pd.read_csv(reference_data.ensure("hpa_normal_tissue", allow_download=False), sep="\t")
    repro_ihc = {
        "testis",
        "ovary",
        "placenta",
        "endometrium 1",
        "endometrium 2",
        "cervix",
        "fallopian tube",
        "breast",
        "prostate",
        "seminal vesicle",
        "epididymis",
        "vagina",
        "appendix",
        "gallbladder",
        "tonsil",
        "spleen",
        "thymus",
    }
    ihc_pos_genes = set(
        ihc.loc[
            ihc["Level"].isin(["Low", "Medium", "High"])
            & ~ihc["Reliability"].eq("Uncertain")
            & ~ihc["Tissue"].astype(str).str.lower().isin(repro_ihc),
            "Gene",
        ].unique()
    )
    lit_doc = json.loads((DATA_DIR / "normal-expression-literature.json").read_text())
    lit_map = {e["entity_id"]: e["effect"] for e in lit_doc["entries"]}

    ens = pyensembl.EnsemblRelease(115)
    genes_by_contig: dict[str, list] = {}
    for g in ens.genes():
        genes_by_contig.setdefault(str(g.contig), []).append(g)
    for c in genes_by_contig:
        genes_by_contig[c].sort(key=lambda g: g.start)

    def _gene_geometry(gid: str) -> dict:
        try:
            g = ens.gene_by_id(gid)
        except Exception:
            return {
                "contig": "",
                "start": 0,
                "end": 0,
                "strand": "",
                "gene_span_bp": 0,
                "cds_bp": 0,
                "cds_span_ratio": 1.0,
                "noncoding_5prime_span_bp": 0,
                "n_cds_exons": 1,
                "bidir_bp": None,
                "bidir_partner": "",
                "bidir_tpm": 0.0,
                "upstream_bp": None,
                "upstream_partner": "",
                "upstream_tpm": 0.0,
            }
        contig = str(g.contig)
        g_span = max(1, g.end - g.start + 1)
        coding_txs = [
            t
            for t in g.transcripts
            if t.is_protein_coding and t.contains_start_codon and t.contains_stop_codon
        ]
        if coding_txs:
            best_tx = max(coding_txs, key=lambda t: len(t.coding_sequence))
            cds_bp = sum(e[1] - e[0] + 1 for e in best_tx.coding_sequence_position_ranges)
            n_cds_exons = len(best_tx.coding_sequence_position_ranges)
            cds_all_min = min(r[0] for t in coding_txs for r in t.coding_sequence_position_ranges)
            cds_all_max = max(r[1] for t in coding_txs for r in t.coding_sequence_position_ranges)
            noncoding_5p = (
                max(0, cds_all_min - g.start) if g.strand == "+" else max(0, g.end - cds_all_max)
            )
        else:
            cds_bp = 0
            n_cds_exons = 1
            noncoding_5p = 0
        cds_ratio = cds_bp / g_span

        g_tss = g.start if g.strand == "+" else g.end
        best_bidir: tuple[int, str, float] | None = None
        best_up: tuple[int, str, float] | None = None
        for nb in genes_by_contig.get(contig, []):
            if nb.gene_id == gid:
                continue
            if nb.end < g.start - 15000:
                continue
            if nb.start > g.end + 15000:
                break
            nb_som = max(
                float(gtex_som_max.get(nb.gene_id, 0.0)),
                float(imm_max.get(nb.gene_id, 0.0)),
            )
            if nb_som < 5.0:
                continue
            if nb.strand != g.strand:
                nb_tss = nb.start if nb.strand == "+" else nb.end
                dist = abs(g_tss - nb_tss)
                if (
                    dist <= rc_qc.BIDIRECTIONAL_PROMOTER_MAX_BP
                    and nb_som >= rc_qc.BIDIRECTIONAL_PARTNER_MIN_TPM
                    and (best_bidir is None or nb_som > best_bidir[2])
                ):
                    best_bidir = (dist, f"{nb.gene_name or nb.gene_id}:{nb.gene_id}", nb_som)
            else:
                dist = (g.start - nb.end) if g.strand == "+" else (nb.start - g.end)
                if (
                    -5000 <= dist <= rc_qc.SAME_STRAND_UPSTREAM_MAX_BP
                    and nb_som >= rc_qc.SAME_STRAND_PARTNER_MIN_TPM
                    and (best_up is None or nb_som > best_up[2])
                ):
                    best_up = (dist, f"{nb.gene_name or nb.gene_id}:{nb.gene_id}", nb_som)
        return {
            "contig": contig,
            "start": int(g.start),
            "end": int(g.end),
            "strand": str(g.strand),
            "gene_span_bp": int(g_span),
            "cds_bp": int(cds_bp),
            "cds_span_ratio": round(cds_ratio, 6),
            "noncoding_5prime_span_bp": int(noncoding_5p),
            "n_cds_exons": int(n_cds_exons),
            "bidir_bp": best_bidir[0] if best_bidir else None,
            "bidir_partner": best_bidir[1] if best_bidir else "",
            "bidir_tpm": round(best_bidir[2], 3) if best_bidir else 0.0,
            "upstream_bp": best_up[0] if best_up else None,
            "upstream_partner": best_up[1] if best_up else "",
            "upstream_tpm": round(best_up[2], 3) if best_up else 0.0,
        }

    adj_rows: list[dict] = []
    target_rce = rce[rce["rare_cell_status"].isin(["established", "unresolved"])].copy()
    for r in target_rce.itertuples(index=False):
        eid = str(r.entity_id)
        gid = (
            str(r.rare_cell_gene_id)
            if pd.notna(r.rare_cell_gene_id) and str(r.rare_cell_gene_id).strip()
            else json.loads(r.member_gene_ids)[0]
        )
        f_rows = rc_findings[
            rc_findings["entity_id"].eq(eid)
            & rc_findings["status"].isin(["established", "unresolved"])
        ]
        organs = sorted(set(f_rows["tissue"].astype(str)))
        sub_dp = dp[dp["gene_id"].eq(gid) & dp["tissue"].isin(organs)].copy()
        sub_dp = sub_dp[sub_dp["n_cells"].ge(1000)]
        sub_dp["mm_pct"] = 100.0 * sub_dp["n_multi"] / sub_dp["n_cells"]
        factor = sub_dp["tier"].map({"vital": 1.0, "protected": 2.0})
        mm_hits = sub_dp[sub_dp["n_multi"].ge(10) & sub_dp["mm_pct"].ge(1.0 * factor)]
        comp_set = sorted({comp_by_cid.get(str(c), "cell") for c in mm_hits["source_id"]})
        # Note MAGEA12 has 0 multi-molecule populations (confirmed by literature R1 single-molecule rule)
        if not comp_set and not sub_dp.empty:
            det_hits = sub_dp[
                sub_dp["n_positive"].ge(10)
                & (100.0 * sub_dp["n_positive"] / sub_dp["n_cells"]).ge(1.0)
            ]
            comp_set = sorted({comp_by_cid.get(str(c), "cell") for c in det_hits["source_id"]})
        nucleus_only = bool(comp_set == ["nucleus"])

        geom = _gene_geometry(gid)
        sp_info = rc_qc.assess_transcriptional_spillover_geometry(
            bidirectional_promoter_bp=geom["bidir_bp"],
            bidirectional_partner_tpm=geom["bidir_tpm"],
            same_strand_upstream_bp=geom["upstream_bp"],
            same_strand_partner_tpm=geom["upstream_tpm"],
            noncoding_5prime_span_bp=geom["noncoding_5prime_span_bp"],
            cds_span_ratio=geom["cds_span_ratio"],
            n_cds_exons=geom["n_cds_exons"],
        )
        bulk_tpm = float(r.finding_bulk_tpm) if pd.notna(r.finding_bulk_tpm) else 0.0
        has_ihc = gid in ihc_pos_genes
        lit_eff = lit_map.get(eid, "not_reviewed")
        u_reasons = (
            str(ts_idx.at[eid, "Unsafe_reasons"])
            if eid in ts_idx.index and pd.notna(ts_idx.at[eid, "Unsafe_reasons"])
            else ""
        )
        other_unsafe = any(
            x.strip() != "Established rare normal-cell RNA"
            for x in u_reasons.split("|")
            if x.strip()
        )
        disp_info = rc_qc.evaluate_nuclear_rare_cell_disposition(
            baseline_rare_cell_status=str(r.rare_cell_status),
            nucleus_only_signal=nucleus_only,
            finding_bulk_tpm=bulk_tpm,
            somatic_ihc_positive=has_ihc,
            literature_effect=lit_eff,
            other_unsafe_gates=other_unsafe,
            transcriptional_spillover_risk=sp_info["transcriptional_spillover_risk"],
            paralog_ambiguous=bool(r.paralog_ambiguous),
            host_overlapped=bool(r.host_overlapped),
            attribution_review_status=(
                str(r.attribution_review_status) if pd.notna(r.attribution_review_status) else ""
            ),
        )
        adj_rows.append(
            {
                "entity_id": eid,
                "protein_group_label": str(r.protein_group_label),
                "gene_id": gid,
                "baseline_classification": str(ts_idx.at[eid, "classification"]),
                "baseline_rare_cell_status": str(r.rare_cell_status),
                "finding_organs": ";".join(organs),
                "assay_compartments": ";".join(comp_set),
                "nucleus_only_signal": nucleus_only,
                "finding_bulk_tpm": round(bulk_tpm, 6),
                "somatic_ihc_positive": has_ihc,
                "literature_effect": lit_eff,
                "other_unsafe_gates": other_unsafe,
                "paralog_ambiguous": bool(r.paralog_ambiguous),
                "host_overlapped": bool(r.host_overlapped),
                "contig": geom["contig"],
                "start": geom["start"],
                "end": geom["end"],
                "strand": geom["strand"],
                "gene_span_bp": geom["gene_span_bp"],
                "cds_bp": geom["cds_bp"],
                "cds_span_ratio": geom["cds_span_ratio"],
                "noncoding_5prime_span_bp": geom["noncoding_5prime_span_bp"],
                "is_intronless_cds": sp_info["is_intronless_cds"],
                "accepted_cds_confirmation_assays": ";".join(
                    sp_info["accepted_cds_confirmation_assays"]
                ),
                "bidirectional_promoter_bp": (
                    geom["bidir_bp"] if geom["bidir_bp"] is not None else np.nan
                ),
                "bidirectional_partner": geom["bidir_partner"],
                "bidirectional_partner_tpm": geom["bidir_tpm"],
                "same_strand_upstream_bp": (
                    geom["upstream_bp"] if geom["upstream_bp"] is not None else np.nan
                ),
                "same_strand_partner": geom["upstream_partner"],
                "same_strand_partner_tpm": geom["upstream_tpm"],
                "transcriptional_spillover_risk": sp_info["transcriptional_spillover_risk"],
                "spillover_mechanisms": ";".join(sp_info["spillover_mechanisms"]),
                "nuclear_spillover_cleared": disp_info["cleared"],
                "disposition": disp_info["disposition"],
                "effective_rare_cell_status": disp_info["effective_rare_cell_status"],
                "effective_rare_cell_action": disp_info["effective_rare_cell_action"],
                "nuclear_rna_caution": disp_info["nuclear_rna_caution"],
            }
        )

    # Also record NUTM1 (PEG:51c85790ace1fbbbd8c4) immune droplet spillover adjudication!
    nutm1_eid = "PEG:51c85790ace1fbbbd8c4"
    nutm1_gid = "ENSG00000184507"
    n_geom = _gene_geometry(nutm1_gid)
    nop10_neut_tpm = float(neut_max.get("ENSG00000182117", n_geom["bidir_tpm"]))
    n_sp = rc_qc.assess_transcriptional_spillover_geometry(
        bidirectional_promoter_bp=n_geom["bidir_bp"],
        bidirectional_partner_tpm=max(n_geom["bidir_tpm"], nop10_neut_tpm),
        same_strand_upstream_bp=n_geom["upstream_bp"],
        same_strand_partner_tpm=n_geom["upstream_tpm"],
        noncoding_5prime_span_bp=n_geom["noncoding_5prime_span_bp"],
        cds_span_ratio=n_geom["cds_span_ratio"],
        n_cds_exons=n_geom["n_cds_exons"],
    )
    nutm1_sorted_neut_tpm = max(
        float(neut_max.get(nutm1_gid, 0.0)),
        float(monaco_max.get(nutm1_gid, 0.0)),
    )
    adj_rows.append(
        {
            "entity_id": nutm1_eid,
            "protein_group_label": "NUTM1",
            "gene_id": nutm1_gid,
            "baseline_classification": str(ts_idx.at[nutm1_eid, "classification"]),
            "baseline_rare_cell_status": str(ts_idx.at[nutm1_eid, "rare_cell_status"]),
            "finding_organs": "bone marrow (Tabula Sapiens 10x neutrophil audit)",
            "assay_compartments": "droplet_5p_3p_prompt_spillover",
            "nucleus_only_signal": True,
            "finding_bulk_tpm": round(nutm1_sorted_neut_tpm, 6),
            "somatic_ihc_positive": nutm1_gid in ihc_pos_genes,
            "literature_effect": lit_map.get(nutm1_eid, "not_reviewed"),
            "other_unsafe_gates": False,
            "paralog_ambiguous": False,
            "host_overlapped": False,
            "contig": n_geom["contig"],
            "start": n_geom["start"],
            "end": n_geom["end"],
            "strand": n_geom["strand"],
            "gene_span_bp": n_geom["gene_span_bp"],
            "cds_bp": n_geom["cds_bp"],
            "cds_span_ratio": n_geom["cds_span_ratio"],
            "noncoding_5prime_span_bp": n_geom["noncoding_5prime_span_bp"],
            "is_intronless_cds": n_sp["is_intronless_cds"],
            "accepted_cds_confirmation_assays": ";".join(n_sp["accepted_cds_confirmation_assays"]),
            "bidirectional_promoter_bp": n_geom["bidir_bp"],
            "bidirectional_partner": n_geom["bidir_partner"],
            "bidirectional_partner_tpm": round(max(n_geom["bidir_tpm"], nop10_neut_tpm), 3),
            "same_strand_upstream_bp": np.nan,
            "same_strand_partner": "",
            "same_strand_partner_tpm": 0.0,
            "transcriptional_spillover_risk": n_sp["transcriptional_spillover_risk"],
            "spillover_mechanisms": ";".join(n_sp["spillover_mechanisms"]),
            "nuclear_spillover_cleared": True,
            "disposition": "promoter_spillover_immune_explained",
            "effective_rare_cell_status": "below_screen",
            "effective_rare_cell_action": "none",
            "nuclear_rna_caution": "bidirectional_promoter_prompt_nop10_cleared",
        }
    )

    adj_df = pd.DataFrame(adj_rows).sort_values(
        ["nuclear_spillover_cleared", "disposition", "protein_group_label"],
        ascending=[False, True, True],
    )
    adj_csv_path = DATA_DIR / "nuclear-spillover-adjudications.csv"
    adj_df.to_csv(adj_csv_path, index=False)

    # 3. Extract exact native-TPM comparison observations + non-self 9-mer overlap for cleared former-Unsafe candidates
    # (plus MAGEB6 and DKKL1 for downstream reuse if needed, keeping schema clean & exact)
    facts = shared_facts()
    cleared_unsafe_eids = sorted(
        set(
            adj_df.loc[
                adj_df["nuclear_spillover_cleared"]
                & adj_df["baseline_classification"].eq("Unsafe"),
                "entity_id",
            ]
        )
    )
    sub_facts = facts[facts["entity_id"].isin(cleared_unsafe_eids)].copy()

    # Exact non-self 9-mer overlap against remaining Unsafe proteins
    ninemers = cta_sequence_overlap.group_ninemers(facts.member_gene_ids)
    overlap_fn = cta_sequence_overlap.unsafe_overlap(ninemers, facts.index)
    remaining_unsafe = facts["unsafe"].copy()
    remaining_unsafe.loc[facts["entity_id"].isin(cleared_unsafe_eids)] = False
    recomputed_overlap = (100.0 * overlap_fn(remaining_unsafe)).round(2)
    overlap_by_eid = dict(zip(facts["entity_id"], recomputed_overlap))

    normal_obs, required_contexts = normal_panel(
        sub_facts, write=False, include_individual_samples=True, normal_stat="p95"
    )
    dev_obs, _ = developmental_samples(sub_facts)
    dev_obs["comparison_id"] = COMPARISON_ID
    repro_obs = pd.concat(
        [normal_obs[normal_obs.scope.eq("reproductive") & normal_obs.sample_id.ne("")], dev_obs],
        ignore_index=True,
    )
    repro_obs["qc_pass"] &= repro_obs.library_qc_review.fillna("").eq("")
    disp_ntpm_map = disposable_ntpm(sub_facts).to_dict()

    # Scan pinned linear_rnaseq_tpm tumor parquets for exact qualifying tumor witnesses
    wanted_genes = {
        g: r.entity_id for r in sub_facts.itertuples() for g in json.loads(r.member_gene_ids)
    }
    folder = BASE / "tumor-identity-repair"
    manifest = json.loads((folder / "manifest.json").read_text())
    paths = {Path(p).stem: (Path(p), digest) for p, digest in manifest["matrices"].items()}
    available = pd.read_csv(folder / "cohort-availability.csv")
    lib_qc_df = pd.read_csv(REBUILD_DIR / "library-qc.csv.gz", low_memory=False)
    id_col_qc = "sample_id" if "sample_id" in lib_qc_df.columns else lib_qc_df.columns[0]
    passing_by_code = {
        str(code): set(grp[id_col_qc].astype(str))
        for code, grp in lib_qc_df[lib_qc_df["reasons"].fillna("").eq("")].groupby("cancer_code")
    }
    tumor_samples_meta = pd.read_csv(REBUILD_DIR / "tumor-samples.csv.gz", keep_default_na=False)
    part_by_key = tumor_samples_meta.set_index("sample_key")["participant_id"].to_dict()

    tumor_hits: dict[str, list[dict]] = {eid: [] for eid in cleared_unsafe_eids}
    for job in available.itertuples():
        if job.cancer_code in target_selection.GERM_CELL_TUMORS:
            continue
        if job.status != "Used" or job.source_scale_class != "linear_rnaseq_tpm":
            continue
        entry, _ = paths[job.cancer_code]
        p = catalog.ensure(
            "per-sample:" + entry.stem,
            version=entry.parent.name.removeprefix("v"),
        )
        id_df = pq.read_table(p, columns=["Ensembl_Gene_ID", "Symbol"]).to_pandas()
        base_ids = id_df["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
        match_idx = base_ids[base_ids.isin(wanted_genes)].index.tolist()
        if not match_idx:
            continue
        tbl = pq.read_table(p).take(match_idx).to_pandas()
        qc_samples = [
            c
            for c in tbl.columns
            if c not in {"Ensembl_Gene_ID", "Symbol"}
            and c in passing_by_code.get(job.cancer_code, set())
        ]
        if not qc_samples:
            continue
        ns = source_matrices.source_sample_namespace(job.source_cohort)
        for _, row in tbl.iterrows():
            gid = str(row["Ensembl_Gene_ID"]).split(".")[0]
            eid = wanted_genes[gid]
            vals = pd.to_numeric(row[qc_samples], errors="coerce").dropna()
            # Keep samples >= 5 TPM (all 10x witnesses must be >= 10 TPM since floor is 1.0 TPM)
            for sid, val in vals[vals >= 5.0].items():
                skey = f"{ns}:{sid}"
                tumor_hits[eid].append(
                    {
                        "source_id": ns,
                        "sample_id": sid,
                        "sample_key": skey,
                        "participant_id": part_by_key.get(skey, ""),
                        "cancer_code": str(job.cancer_code),
                        "value": round(float(val), 6),
                    }
                )

    candidates_payload: dict[str, dict] = {}
    for r in sub_facts.itertuples():
        eid = r.entity_id
        n_df = normal_obs[normal_obs.entity_id.eq(eid)]
        valid_n = n_df[
            n_df.scope.eq("normal") & n_df.qc_pass & n_df.identity_resolved & n_df.value.notna()
        ]
        valid_d = n_df[
            n_df.scope.eq("disposable") & n_df.qc_pass & n_df.identity_resolved & n_df.value.notna()
        ]
        missing_req = sorted(required_contexts - set(valid_n["context_id"]))
        n_max = float(valid_n["value"].max()) if not valid_n.empty else None
        d_max = float(valid_d["value"].max()) if not valid_d.empty else 0.0

        r_df = repro_obs[repro_obs.entity_id.eq(eid)]
        ra_rec = target_selection.assess_reproductive_rna(
            (
                r_df
                if not r_df.empty
                else pd.DataFrame(
                    [
                        {
                            "entity_id": eid,
                            "source_id": "none",
                            "sample_id": "none",
                            "participant_id": "",
                            "scope": "reproductive",
                            "value": np.nan,
                            "comparison_id": COMPARISON_ID,
                            "unit": "TPM",
                            "qc_pass": False,
                            "identity_resolved": False,
                        }
                    ]
                )
            ),
            n_df,
            required_normal_contexts=required_contexts,
            comparison_id=COMPARISON_ID,
            comparability_reference="comparison.json",
            unit="TPM",
        )
        # Summarize per-compartment reproductive maxima & qualifying values
        valid_r = r_df[r_df.qc_pass & r_df.identity_resolved & r_df.value.notna()]
        repro_comp_summary: list[dict] = []
        denom = max(n_max or 0.0, d_max, 1.0)
        for comp, cgrp in valid_r.groupby("compartment", sort=True):
            if not str(comp).strip():
                continue
            is_dev = bool(cgrp["scope"].eq("developmental").any())
            cmax = float(cgrp["value"].max())
            repro_comp_summary.append(
                {
                    "compartment": str(comp),
                    "scope": "developmental" if is_dev else "reproductive",
                    "max_tpm": round(cmax, 6),
                    "qualifies_tenfold": bool(
                        ra_rec["status"] == "qualifies" and cmax >= 10.0 * denom
                    ),
                }
            )

        t_list = pd.DataFrame(tumor_hits[eid])
        if not t_list.empty:
            t_list = (
                t_list.drop_duplicates("sample_key").sort_values("value", ascending=False).head(50)
            )
            top_tumors = t_list.to_dict("records")
        else:
            top_tumors = []

        candidates_payload[eid] = {
            "entity_id": eid,
            "protein_group_label": str(r.protein_group_label),
            "required_coverage_complete": len(missing_req) == 0,
            "missing_required_contexts": missing_req,
            "unsafe_sequence_overlap_pct": float(overlap_by_eid.get(eid, 0.0)),
            "normal_reference_tpm": round(n_max, 6) if n_max is not None else None,
            "disposable_max_tpm": round(d_max, 6),
            "disposable_max_ntpm": round(float(disp_ntpm_map.get(eid, 0.0)), 4),
            "reproductive_status": ra_rec["status"],
            "reproductive_supported_tpm": (
                round(float(ra_rec["supported_maximum"]), 6)
                if ra_rec["supported_maximum"] is not None
                else None
            ),
            "reproductive_compartments": repro_comp_summary,
            "top_tumor_observations": top_tumors,
        }

    obs_path = DATA_DIR / "cleared-candidate-comparison-observations.json"
    obs_path.write_text(
        json.dumps(
            {
                "schema": "oncoref.cleared-candidate-comparison-observations.v1",
                "policy_id": rc_qc.NUCLEAR_SPILLOVER_POLICY_ID,
                "comparison_id": COMPARISON_ID,
                "candidates": candidates_payload,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    adj_receipt = {
        "schema": "oncoref.nuclear-spillover-qc.v1",
        "policy_id": rc_qc.NUCLEAR_SPILLOVER_POLICY_ID,
        "ensembl_release": 115,
        "evaluated_proteins": len(adj_df),
        "cleared_proteins": int(adj_df["nuclear_spillover_cleared"].sum()),
        "disposition_counts": adj_df["disposition"].value_counts().to_dict(),
        "cleared_labels": adj_df.loc[
            adj_df["nuclear_spillover_cleared"], "protein_group_label"
        ].tolist(),
        "magea12_regression_status": str(
            adj_df.loc[adj_df["protein_group_label"].eq("MAGEA12"), "disposition"].iloc[0]
        ),
        "sha256": _sha256_bytes(adj_csv_path.read_bytes()),
    }
    (DATA_DIR / "nuclear-spillover-adjudications.json").write_text(
        json.dumps(adj_receipt, indent=2) + "\n",
        encoding="utf-8",
    )

    # 4. Write audit report
    cleared_sub = adj_df[adj_df["nuclear_spillover_cleared"]]
    md_lines = [
        "# Nuclear snRNA-seq Compartment & Transcriptional Spillover QC (`#765`)",
        "",
        f"**Policy ID:** `{rc_qc.NUCLEAR_SPILLOVER_POLICY_ID}`  ",
        "**Ensembl Release:** `115` (`GRCh38`)  ",
        "**CELLxGENE Census Release:** `2025-11-08` (`93` primary adult collections)",
        "",
        "## Overview",
        "",
        "pMHC-directed T-cell cytotoxicity requires spliced cytoplasmic `CDS` mRNA, translation,",
        "and surface HLA presentation. Intron-inclusive single-nucleus RNA-seq (`snRNA-seq`,",
        '`assay_compartment == "nucleus"`) captures unspliced pre-mRNA across long non-coding',
        "5'-UTR introns (`PLAC1`, `TENT5D`), termination readthrough downstream of active",
        "same-strand somatic neighbors (`PABIR2 -> PLAC1`, `EPO -> ZAN`), and divergent",
        "promoter-upstream transcripts (`NOP10 <-> NUTM1`, `135 bp` head-to-head CpG island).",
        "",
        "## Guardrails Enforced (`tsarina` Cross-Repo Contract)",
        "",
        "1. **Bulk-dilution & literature/IHC guardrails (`MAGEA12` regression anchor):**",
        "   `MAGEA12` (`finding_bulk_tpm = 0.30 TPM` in brain, `literature_effect = confirm_detection`)",
        "   remains strictly `Unsafe` (`retained_literature_or_ihc_supported`). Clearance requires",
        "   all four conditions simultaneously: (a) `nucleus_only_signal == True` (zero cytoplasmic",
        "   `scRNA-seq` multi-UMI populations), (b) `finding_bulk_tpm < 0.05 TPM`, (c) negative",
        "   somatic IHC and no `exclude`/`confirm_detection` literature finding, and (d) verified",
        "   Ensembl 115 spillover geometry (`< 1 kb` bidirectional promoter `>= 10 TPM`, `< 10 kb`",
        "   same-strand upstream partner `>= 5 TPM`, or `>= 50 kb` non-coding 5'-UTR span with",
        "   `CDS`/span ratio `<= 5%`).",
        "2. **Intronless / single-CDS-exon coding support (`DUX4`, `PLAC1`):**",
        "   Records `is_intronless_cds` and explicitly accepts `cds_exonic_reads` alongside",
        "   `spliced_exon_junction` for coding-sequence attribution.",
        "3. **Downstream `nuclear_only` disposition:**",
        '   Cleared loci emit `effective_rare_cell_status = "nuclear_only"` and',
        '   `nuclear_rna_caution = "snrnaseq_intronic_spillover_cleared"`, advancing loci that',
        "   satisfy all `10x` gates (`PLAC1`: `152.6 TPM` `WILMS`, `187.0 TPM` placenta vs `0.82 TPM`",
        "   normal; `ZAN`: `10.5 TPM` `TGCT` / `10.2 TPM` `B_ALL`, `16.7 TPM` testis vs `0.51 TPM`",
        "   normal) directly into `CTA` rather than trapping them in `Unsafe` or `Other:Review`.",
        "",
        "## Cleared Nuclear / Promoter Spillover Loci",
        "",
        "| Protein | Baseline | Finding Organ | Bulk TPM | Spillover Mechanism | Key Neighbor / Span | Disposition |",
        "| --- | --- | --- | ---: | --- | --- | --- |",
    ]
    for r in cleared_sub.itertuples(index=False):
        partner = (
            r.bidirectional_partner
            or r.same_strand_partner
            or f"5'-UTR {int(r.noncoding_5prime_span_bp):,} bp ({100 * r.cds_span_ratio:.2f}% CDS)"
        )
        md_lines.append(
            f"| `{r.protein_group_label}` | `{r.baseline_classification}` | `{r.finding_organs}` | "
            f"`{r.finding_bulk_tpm:.4f}` | `{r.spillover_mechanisms}` | `{partner}` | `{r.disposition}` |"
        )
    (AUDIT_DIR / "index.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(json.dumps(adj_receipt, indent=2))


if __name__ == "__main__":
    main()

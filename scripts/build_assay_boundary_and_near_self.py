#!/usr/bin/env python3
"""Build SHA-256-verified assay-boundary (#768) & 1-mm Hamming / tryptic-MS (#769) tables.

Generates:
  - oncoref/data/assay-boundary-adjudications.{csv,json}
  - oncoref/data/assay-boundary-comparison-observations.json
  - oncoref/data/candidate-near-self-1mm.{csv,json}
  - oncoref/data/tryptic-ms-paralog-audit.{csv,json}
  - docs/audits/assay-boundary-near-self-20261123/index.md
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from rebuild_target_catalog import (  # noqa: E402
    BASE,
    COMPARISON_ID,
    LEGACY_RARE,
    disposable_ntpm,
    normal_panel,
    shared_facts,
)

from oncoref import (  # noqa: E402
    assay_boundary_adjudications as ab,
)
from oncoref import (  # noqa: E402
    catalog,
    cta_sequence_overlap,
    curation,
    normal_adjudications,
    normal_cell_data,
    normal_library_qc,
    rare_cell_compartment_qc,
    reference_data,
    sequence_bundle,
    source_matrices,
    target_selection,
    tgct_localization,
    unmeasured_home_panels,
)
from oncoref.gene_matching import match_reference_genes, resolve_measurement_rows  # noqa: E402
from oncoref.load_dataset import get_data  # noqa: E402

DATA_DIR = ROOT / "oncoref/data"
REBUILD_DIR = ROOT / "docs/audits/target-selection-20261110/rebuild"
AUDIT_DIR = ROOT / "docs/audits/assay-boundary-near-self-20261123"


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_csv_with_receipt(df: pd.DataFrame, stem: str, extra_meta: dict) -> dict:
    csv_path = DATA_DIR / f"{stem}.csv"
    json_path = DATA_DIR / f"{stem}.json"
    df.to_csv(csv_path, index=False)
    digest = _sha256_bytes(csv_path.read_bytes())
    receipt = {
        "schema": f"oncoref.{stem}.v1",
        "policy_id": ab.ASSAY_BOUNDARY_POLICY_ID,
        "rows": len(df),
        "sha256": digest,
        **extra_meta,
    }
    json_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def _developmental_samples_readonly(sub_facts: pd.DataFrame) -> pd.DataFrame:
    """Compute Petropoulos/Tyser developmental samples without mutating frozen audit files."""
    sources = [
        ("petropoulos_E-MTAB-3929", "petropoulos_tpm", "petropoulos_design"),
        ("tyser_E-MTAB-9388", "tyser_tpm", "tyser_design"),
    ]
    rows: list[dict] = []
    for source, archive_role, design_role in sources:
        archive_path = curation.reference_input(archive_role)
        design_path = curation.reference_input(design_role)
        with zipfile.ZipFile(archive_path) as archive:
            name = next(n for n in archive.namelist() if n.endswith(".mtx"))
            genes = pd.read_csv(
                archive.open(name + "_rows"), sep="\t", header=None, names=["gene_id", "label"]
            )
            cells = pd.read_csv(archive.open(name + "_cols"), sep="\t", header=None)[0]
            matrix = mmread(archive.open(name)).tocsr()
        identity = match_reference_genes(
            genes,
            id_column="gene_id",
            symbol_column="label",
            source=source,
            source_version="Ensembl106_filtered_TPM",
        )
        decisions = resolve_measurement_rows(
            identity, equivalent=lambda a, b, values=matrix: (values[a] != values[b]).nnz == 0
        )
        chosen = decisions[decisions.row_position.notna()].reset_index(drop=True)
        matrix = matrix[chosen.row_position.astype(int).to_numpy()]
        meta = pd.read_csv(design_path, sep="\t").set_index("Assay").loc[cells]
        if source.startswith("petropoulos"):
            meta["embryo"] = meta["Sample Characteristic[individual]"]
            meta["stage"] = meta["Sample Characteristic[age]"]
            meta["lineage"] = (
                meta["Sample Characteristic[age]"]
                + " / "
                + meta["Factor Value[inferred lineage]"].fillna("unspecified lineage")
            )
        else:
            meta["embryo"] = "single-embryo"
            meta["stage"] = "Carnegie stage 7"
            meta["lineage"] = (
                "Carnegie stage 7 / "
                + meta["Sample Characteristic[inferred cell type - authors labels]"]
            )
        positions = dict(zip(chosen.gene_id, range(len(chosen))))
        groups = list(meta.groupby(["embryo", "lineage"]).indices.items())
        for protein in sub_facts.itertuples():
            members = json.loads(protein.member_gene_ids)
            complete = all(g in positions for g in members)
            values = (
                np.asarray(matrix[[positions[g] for g in members]].sum(axis=0)).ravel()
                if complete
                else None
            )
            for (embryo, lineage), indices in groups:
                stage = str(meta.iloc[indices].stage.iloc[0])
                n_positive = int((values[indices] >= 1).sum()) if complete else 0
                donors_required = target_selection.developmental_donor_requirement(
                    len(indices), n_positive
                )
                supported = complete and donors_required > 0
                rows.append(
                    {
                        "entity_id": protein.entity_id,
                        "source_id": source,
                        "sample_id": f"{embryo} / {lineage}",
                        "participant_id": f"{source}:{embryo}",
                        "tissue": "embryo",
                        "compartment": "embryo",
                        "source_context": lineage,
                        "cell_type": lineage,
                        "scope": "developmental",
                        "statistic": "mean native TPM within one embryo and lineage",
                        "observation_kind": "embryo_lineage_mean",
                        "unit": "TPM",
                        "value": float(values[indices].mean()) if complete else np.nan,
                        "qc_pass": bool(supported),
                        "identity_resolved": complete,
                        "cells": len(indices),
                        "positive_cells": n_positive,
                        "developmental_stage": stage,
                        "required_embryos": donors_required,
                        "comparison_id": COMPARISON_ID,
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    frozen_ts_path = (
        ROOT / "tests/fixtures/release-1.50-inputs/oncoref/data/target-selection.csv.gz"
    )
    if frozen_ts_path.is_file():
        ts = pd.read_csv(frozen_ts_path, compression="gzip", low_memory=False)
    else:
        ts = target_selection._baseline_2026_11_20_frame()
    ts_by_eid = ts.set_index("entity_id")
    facts = shared_facts()
    facts_by_eid = facts.set_index("entity_id")
    kept_unsafe_cols = [
        c for c in facts.columns if c.startswith("unsafe:") and c[7:] not in LEGACY_RARE
    ]

    adj_rows: list[dict] = []

    # ------------------------------------------------------------------
    # Rule A (#713 / #768): Non-admixed tumor witness gating across R-TGCT
    # ------------------------------------------------------------------
    bleed_deltas = {
        str(eid): 10.0 * max(float(val), 1.0)
        for eid, val in normal_library_qc.multiplex_bleed_gene_deltas(only_applied=True)
        .set_index("entity_id")["clean_normal_reference_tpm"]
        .items()
    }
    obs_737 = get_data("human-737-tumor-observations")
    peak_737_by_eid: dict[str, tuple[float, str, str]] = {}
    for _, grp737 in obs_737.groupby("entity_id"):
        top737 = grp737.sort_values("value", ascending=False).iloc[0]
        peak_737_by_eid[str(top737["entity_id"])] = (
            float(top737["value"]),
            str(top737["cancer_code"]),
            str(top737["sample_id"]),
        )
    tgct_restr = tgct_localization.tgct_admixture_restricted_entity_ids(include_multi_r_tgct=True)
    tgct_sup = tgct_localization.tgct_aafl_protein_support()
    for r in tgct_sup.itertuples(index=False):
        eid = str(r.entity_id)
        if eid not in tgct_restr:
            continue
        eff_thresh = float(bleed_deltas.get(eid, r.single_tumor_threshold_tpm))
        nt_tpm = float(r.max_non_tgct_tpm) if pd.notna(r.max_non_tgct_tpm) else 0.0
        nt_code = str(r.top_non_tgct_cancer_code) if pd.notna(r.top_non_tgct_cancer_code) else ""
        nt_samp = str(r.top_non_tgct_sample) if pd.notna(r.top_non_tgct_sample) else ""
        if eid in peak_737_by_eid and peak_737_by_eid[eid][0] > nt_tpm:
            nt_tpm, nt_code, nt_samp = peak_737_by_eid[eid]
        ev = ab.evaluate_non_admixed_tumor_witness(
            entity_id=eid,
            max_non_tgct_tpm=nt_tpm,
            single_tumor_threshold_tpm=eff_thresh,
            admixture_restricted_ids=tgct_restr,
        )
        adj_rows.append(
            {
                "entity_id": eid,
                "protein_group_label": str(r.protein_group_label),
                "rule_family": "tgct_admixture_witness_713",
                "baseline_classification": str(ts_by_eid.at[eid, "classification"]),
                "metric_primary": round(nt_tpm, 6),
                "metric_secondary": round(eff_thresh, 6),
                "evidence_detail": (
                    f"top_non_tgct={nt_code}:{nt_samp} "
                    f"({nt_tpm:.4f} vs thresh {eff_thresh:.4f} TPM); support_tier={r.support_tier}"
                ),
                "action_applied": bool(ev["admixture_sole_witness_excluded"]),
                "disposition": ev["disposition"],
            }
        )

    # ------------------------------------------------------------------
    # Rule B (#768): RNA-zero single-antibody HPA IHC reconciliation
    # ------------------------------------------------------------------
    ihc_unrel = get_data("cta-ihc-unreliable")
    ihc_cand_eids: list[str] = []
    for r in ihc_unrel.itertuples(index=False):
        m = ts[ts["representative_gene_id"].eq(r.Ensembl_Gene_ID)]
        if m.empty:
            continue
        mr = m.iloc[0]
        eid = str(mr["entity_id"])
        if mr["category"] != "Unsafe":
            continue
        sfr = facts_by_eid.loc[eid]
        fired = [c for c in kept_unsafe_cols if bool(sfr[c])]
        sole_ihc = (
            fired == ["unsafe:HPA normal protein"] and str(mr["rare_cell_status"]) != "established"
        )
        ihc_cand_eids.append(eid)
        adj_rows.append(
            {
                "entity_id": eid,
                "protein_group_label": str(mr["protein_group_label"]),
                "rule_family": "ihc_rna_zero_reconciliation",
                "baseline_classification": str(mr["classification"]),
                "metric_primary": float(mr["healthy_maximum_tpm"] or 0.0),
                "metric_secondary": float(mr["healthy_maximum_ntpm"] or 0.0),
                "evidence_detail": f"fired={fired}; reason={r.reason}",
                "action_applied": False,  # updated below after normal_panel verifies required_coverage_complete!
                "disposition": "pending_coverage_check"
                if sole_ihc
                else "retained_unsafe_ihc_or_rna_discordant_guard_failed",
            }
        )

    # ------------------------------------------------------------------
    # Rule C (#768): Adult cornea multi-UMI floor (>=10 cells with >=2 UMIs)
    # ------------------------------------------------------------------
    cornea_cands: list[tuple[str, str]] = []
    for eid, sfr in facts_by_eid.iterrows():
        mr = ts_by_eid.loc[eid]
        if mr["category"] != "Unsafe" or str(mr["rare_cell_status"]) == "established":
            continue
        fired = [c for c in kept_unsafe_cols if bool(sfr[c])]
        if fired == ["unsafe:Adult cornea"]:
            cornea_cands.append((str(eid), str(mr["protein_group_label"])))

    cornea_p = reference_data.ensure("cornea_collin_adult_counts", allow_download=False)
    want_cornea_syms = {lbl for _, lbl in cornea_cands}
    cornea_counts: dict[str, tuple[int, int, float]] = {}
    rows_iter = normal_cell_data._count_rows(cornea_p)
    cells_hdr = next(rows_iter)
    n_cornea_cells = len(cells_hdr)
    for gname, cnts in rows_iter:
        if gname in want_cornea_syms:
            cornea_counts[gname] = (
                int((cnts >= 1).sum()),
                int((cnts >= 2).sum()),
                float(cnts.max()),
            )

    cornea_cleared_eids: list[str] = []
    for eid, lbl in cornea_cands:
        n_pos, n_ge2, max_umi = cornea_counts.get(lbl, (0, 0, 0.0))
        ev = ab.evaluate_adult_cornea_multi_umi(
            sole_cornea_unsafe_gate=True,
            rare_cell_established=False,
            n_cells_ge2_umi=n_ge2,
        )
        if ev["cleared"]:
            cornea_cleared_eids.append(eid)
        adj_rows.append(
            {
                "entity_id": eid,
                "protein_group_label": lbl,
                "rule_family": "adult_cornea_multi_umi_floor",
                "baseline_classification": "Unsafe",
                "metric_primary": float(n_ge2),
                "metric_secondary": float(n_pos),
                "evidence_detail": f"GSE155683 adult cornea ({n_cornea_cells} cells): n_ge1={n_pos}, n_ge2={n_ge2}, max_umi={max_umi:.0f}",
                "action_applied": bool(ev["cleared"]),
                "disposition": ev["disposition"],
            }
        )

    # ------------------------------------------------------------------
    # Rule E (#779): Droplet ambient milk soup & <=100 kb SCPP cis-spillover
    # ------------------------------------------------------------------
    from pyensembl import EnsemblRelease

    from oncoref import (
        cta_selection_policy,
        gtex,
        human_protein_atlas,
        normal_tissues,
        tissue_tiers,
    )

    ens115 = EnsemblRelease(115)
    sc_raw_full = human_protein_atlas.single_cell_tissues(qc=None, cleanup=None)
    sc_clean = human_protein_atlas.single_cell_tissues(qc="all", cleanup="all")
    hpa_bulk_clean = human_protein_atlas.rna_tissue_hpa(qc="splash", cleanup="all")
    hpa_ihc_all = human_protein_atlas.normal_tissue()

    repro_sc_only = {"testis", "ovary", "placenta", "endometrium"}

    def _filter_somatic_sc(df_sc: pd.DataFrame) -> pd.DataFrame:
        return df_sc[
            ~df_sc["Tissue"].str.lower().isin(repro_sc_only)
            & ~(
                df_sc["Tissue"].str.lower().eq("breast")
                & df_sc["Cell type"].str.contains("breast", case=False, na=False)
            )
        ]

    somatic_sc_raw = _filter_somatic_sc(sc_raw_full)
    somatic_sc_all = _filter_somatic_sc(sc_clean)
    dom_sc_rows = somatic_sc_all[somatic_sc_all["nTPM"] >= ab.SC_CIS_SPILLOVER_MIN_NEIGHBOR_NTPM]
    dom_by_tc: dict[tuple[str, str], list[tuple[str, int, int, float, str]]] = {}
    for d_tis, d_cl, d_gid, d_sym, d_ntpm in zip(
        dom_sc_rows["Tissue"],
        dom_sc_rows["Cluster"],
        dom_sc_rows["Gene"],
        dom_sc_rows["Gene name"],
        dom_sc_rows["nTPM"],
    ):
        try:
            dg = ens115.gene_by_id(str(d_gid))
            dom_by_tc.setdefault((str(d_tis), str(d_cl)), []).append(
                (str(dg.contig), int(dg.start), int(dg.end), float(d_ntpm), str(d_sym))
            )
        except Exception:
            pass

    sc_raw_by_gene = {str(g): grp for g, grp in somatic_sc_raw.groupby("Gene", sort=False)}
    sc_clean_by_gene = {str(g): grp for g, grp in somatic_sc_all.groupby("Gene", sort=False)}
    hpa_bulk_by_gene = {
        (str(r.Gene), str(r.Tissue).lower()): float(r.TPM)
        for r in hpa_bulk_clean.itertuples(index=False)
        if pd.notna(r.TPM)
    }
    pos_ihc_gt = set()
    pos_ihc_mask = hpa_ihc_all["Level"].isin(["Low", "Medium", "High"]) & hpa_ihc_all[
        "Reliability"
    ].isin(["Approved", "Supported", "Enhanced"])
    for ir in hpa_ihc_all.loc[pos_ihc_mask].itertuples(index=False):
        pos_ihc_gt.add((str(ir.Gene), str(ir.Tissue).lower()))

    tier_cache: dict[tuple[str, str], float] = {}

    def _cluster_base_ceiling(tissue_lbl: str, cell_lbl: str) -> float:
        key = (tissue_lbl.lower(), cell_lbl.lower())
        if key not in tier_cache:
            ann = normal_tissues.resolve_normal_tissue(
                tissue_lbl, species="human", source_id="hpa_single_cell_tissue", context="normal"
            )
            organ = (ann.concept_id or tissue_lbl.lower()).replace("_", " ")
            tier = tissue_tiers.consequences(
                organ, cell_lbl.lower(), stage="adult", resolution="cell"
            )["tier"]
            tier_cache[key] = float(
                cta_selection_policy.TIER_BASES.get(
                    tier, cta_selection_policy.MANAGEABLE_UNSAFE_FLOOR
                )
            )
        return tier_cache[key]

    sc_cleared_eids: list[str] = []
    audit_sc_controls = {"CSN2", "CSN1S1", "CSN3", "ANKRD30A", "MAGEA3"}
    for eid, sfr in facts_by_eid.iterrows():
        mr = ts_by_eid.loc[eid]
        if mr["category"] != "Unsafe":
            continue
        lbl = str(mr["protein_group_label"])
        fired = [c for c in kept_unsafe_cols if bool(sfr[c])]
        sole_sc = fired == ["unsafe:HPA tissue clusters"]
        rc_est = str(mr["rare_cell_status"]) == "established"
        if not sole_sc and lbl not in audit_sc_controls:
            continue

        mems = json.loads(str(mr["member_gene_ids"]))
        raw_parts = [sc_raw_by_gene[g] for g in mems if g in sc_raw_by_gene]
        if raw_parts:
            grp_raw = (
                pd.concat(raw_parts, ignore_index=True)
                .groupby(["Tissue", "Cluster", "Cell type"], as_index=False)["nTPM"]
                .sum()
            )
        else:
            grp_raw = pd.DataFrame(columns=["Tissue", "Cluster", "Cell type", "nTPM"])

        sub_parts = [sc_clean_by_gene[g] for g in mems if g in sc_clean_by_gene]
        if sub_parts:
            grp_sc = (
                pd.concat(sub_parts, ignore_index=True)
                .groupby(["Tissue", "Cluster", "Cell type"], as_index=False)["nTPM"]
                .sum()
            )
        else:
            grp_sc = pd.DataFrame(columns=["Tissue", "Cluster", "Cell type", "nTPM"])

        raw_above_keys: set[tuple[str, str]] = set()
        for r_tis, r_cl, r_ct, r_ntpm in zip(
            grp_raw["Tissue"], grp_raw["Cluster"], grp_raw["Cell type"], grp_raw["nTPM"]
        ):
            if float(r_ntpm) >= _cluster_base_ceiling(str(r_tis), str(r_ct)):
                raw_above_keys.add((str(r_tis), str(r_cl)))

        clean_above_keys: set[tuple[str, str]] = set()
        unexplained = 0
        max_ratio = 0.0
        max_bulk_tpm = 0.0
        any_ihc_pos = False
        flagged_organs: set[str] = set()
        mem_coords = []
        for g in mems:
            try:
                mg = ens115.gene_by_id(g)
                mem_coords.append((str(mg.contig), int(mg.start), int(mg.end)))
            except Exception:
                pass

        for s_tis, s_cl, s_ct, s_ntpm in zip(
            grp_sc["Tissue"], grp_sc["Cluster"], grp_sc["Cell type"], grp_sc["nTPM"]
        ):
            ceil = _cluster_base_ceiling(str(s_tis), str(s_ct))
            if float(s_ntpm) < ceil:
                continue
            clean_above_keys.add((str(s_tis), str(s_cl)))
            t_low = str(s_tis).lower()
            flagged_organs.add(t_low)
            dom_list = dom_by_tc.get((str(s_tis), str(s_cl)), [])
            cis_sum = 0.0
            for d_chr, d_st, d_en, d_nt, _ in dom_list:
                if any(
                    c_chr == d_chr
                    and max(0, max(c_st - d_en, d_st - c_en)) <= ab.SC_CIS_SPILLOVER_MAX_DIST_BP
                    for c_chr, c_st, c_en in mem_coords
                ):
                    cis_sum += d_nt
            if cis_sum >= ab.SC_CIS_SPILLOVER_MIN_NEIGHBOR_NTPM:
                rat = float(s_ntpm) / cis_sum
                max_ratio = max(max_ratio, rat)
                if rat > ab.SC_CIS_SPILLOVER_MAX_RATIO:
                    unexplained += 1
            else:
                unexplained += 1

        soup_cleared = len(raw_above_keys - clean_above_keys)
        flagged_organs |= {t.lower() for t, _ in raw_above_keys if t.lower() != "breast"}
        for org in flagged_organs:
            b_sum = sum(hpa_bulk_by_gene.get((g, org), 0.0) for g in mems)
            max_bulk_tpm = max(max_bulk_tpm, b_sum)
            if any((g, org) in pos_ihc_gt for g in mems):
                any_ihc_pos = True

        ev_sc = ab.evaluate_sc_cluster_soup_and_spillover(
            sole_sc_cluster_unsafe_gate=sole_sc,
            rare_cell_established=rc_est,
            unexplained_somatic_clusters_above_ceiling=unexplained,
            max_cis_spillover_ratio=max_ratio,
            flagged_organ_hpa_bulk_tpm=max_bulk_tpm,
            flagged_organ_ihc_positive=any_ihc_pos,
            raw_somatic_clusters_above_ceiling=len(raw_above_keys),
            soup_cleared_clusters=soup_cleared,
        )
        if ev_sc["cleared"]:
            sc_cleared_eids.append(str(eid))
        if ev_sc["cleared"] or lbl in audit_sc_controls:
            adj_rows.append(
                {
                    "entity_id": str(eid),
                    "protein_group_label": lbl,
                    "rule_family": "sc_cluster_soup_and_cis_spillover",
                    "baseline_classification": "Unsafe",
                    "metric_primary": float(unexplained),
                    "metric_secondary": round(float(max_ratio), 8),
                    "evidence_detail": (
                        f"fired={fired}; flagged_organs={sorted(flagged_organs)}; "
                        f"unexplained_clusters={unexplained}; max_cis_ratio={max_ratio:.6f}; "
                        f"hpa_bulk_tpm={max_bulk_tpm:.4f}; ihc_pos={any_ihc_pos}"
                    ),
                    "action_applied": bool(ev_sc["cleared"]),
                    "disposition": ev_sc["disposition"],
                }
            )

    # ------------------------------------------------------------------
    # Extract exact native-TPM comparison observations for former-Unsafe candidates
    # ------------------------------------------------------------------
    check_unsafe_eids = sorted(set(ihc_cand_eids) | set(cornea_cleared_eids) | set(sc_cleared_eids))
    sub_facts = facts[facts["entity_id"].isin(check_unsafe_eids)].copy()

    ninemers = cta_sequence_overlap.group_ninemers(facts.member_gene_ids)
    overlap_fn = cta_sequence_overlap.unsafe_overlap(ninemers, facts.index)
    remaining_unsafe = facts["unsafe"].copy()
    remaining_unsafe.loc[facts["entity_id"].isin(check_unsafe_eids)] = False
    recomputed_overlap = (100.0 * overlap_fn(remaining_unsafe)).round(2)
    overlap_by_eid = dict(zip(facts["entity_id"], recomputed_overlap))

    normal_obs, required_contexts = normal_panel(
        sub_facts, write=False, include_individual_samples=True, normal_stat="p95"
    )
    if sc_cleared_eids:
        sc_facts = sub_facts[sub_facts["entity_id"].isin(sc_cleared_eids)]
        sc_gids = sorted({g for s in sc_facts["member_gene_ids"] for g in json.loads(s)})
        clean_gtex_p95 = gtex.tissue_summary(
            genes=sc_gids, stat="p95", qc="all", cleanup="all"
        ).copy()
        clean_gtex_p95["gid"] = clean_gtex_p95["Name"].astype(str).str.split(".").str[0]
        clean_gtex_by_gid = clean_gtex_p95.set_index("gid")
        site_cols_p95 = [
            c for c in clean_gtex_p95.columns if c not in ("Name", "Description", "gid")
        ]
        for sc_r in sc_facts.itertuples(index=False):
            sc_eid = str(sc_r.entity_id)
            sc_mems = json.loads(str(sc_r.member_gene_ids))
            if not all(g in clean_gtex_by_gid.index for g in sc_mems):
                continue
            sum_by_site = clean_gtex_by_gid.loc[sc_mems, site_cols_p95].sum(axis=0)
            mask_p95 = (
                normal_obs["entity_id"].eq(sc_eid)
                & normal_obs["source_id"].eq("gtex_donor_p95:v11")
                & normal_obs["scope"].isin(["normal", "disposable"])
            )
            for idx_n in normal_obs.index[mask_p95]:
                st_key = str(normal_obs.at[idx_n, "source_context"]).lower()
                if st_key in sum_by_site.index and pd.notna(sum_by_site[st_key]):
                    normal_obs.at[idx_n, "value"] = float(sum_by_site[st_key])
    dev_obs = _developmental_samples_readonly(sub_facts)
    repro_obs = pd.concat(
        [normal_obs[normal_obs.scope.eq("reproductive") & normal_obs.sample_id.ne("")], dev_obs],
        ignore_index=True,
    )
    if "library_qc_review" in repro_obs.columns:
        repro_obs["qc_pass"] &= repro_obs.library_qc_review.fillna("").eq("")
    disp_ntpm_map = disposable_ntpm(sub_facts).to_dict()

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

    tumor_hits: dict[str, list[dict]] = {eid: [] for eid in check_unsafe_eids}
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
        eid = str(r.entity_id)
        n_df = normal_obs[normal_obs.entity_id.eq(eid)]
        valid_n = n_df[
            n_df.scope.eq("normal") & n_df.qc_pass & n_df.identity_resolved & n_df.value.notna()
        ]
        valid_d = n_df[
            n_df.scope.eq("disposable") & n_df.qc_pass & n_df.identity_resolved & n_df.value.notna()
        ]
        missing_req = sorted(required_contexts - set(valid_n["context_id"]))
        cov_ok = len(missing_req) == 0
        n_max = float(valid_n["value"].max()) if not valid_n.empty else None
        d_max = float(valid_d["value"].max()) if not valid_d.empty else 0.0

        # Finalize IHC row with exact required_coverage_complete & normal_reference_tpm
        for arow in adj_rows:
            if arow["entity_id"] == eid and arow["rule_family"] == "ihc_rna_zero_reconciliation":
                sole_ok = arow["disposition"] == "pending_coverage_check"
                som_tpm = max(float(n_max or 0.0), float(arow["metric_primary"]))
                som_ntpm = float(arow["metric_secondary"])
                ev_ihc = ab.evaluate_ihc_rna_reconciliation(
                    unreliable_ihc_curated=True,
                    sole_ihc_unsafe_gate=sole_ok,
                    required_coverage_complete=cov_ok,
                    somatic_rna_max_tpm=som_tpm,
                    somatic_rna_max_ntpm=som_ntpm,
                )
                arow["metric_primary"] = round(som_tpm, 6)
                arow["action_applied"] = bool(ev_ihc["cleared"])
                arow["disposition"] = ev_ihc["disposition"]

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
        top_tumors = (
            t_list.drop_duplicates("sample_key")
            .sort_values("value", ascending=False)
            .head(50)
            .to_dict("records")
            if not t_list.empty
            else []
        )
        candidates_payload[eid] = {
            "entity_id": eid,
            "protein_group_label": str(r.protein_group_label),
            "required_coverage_complete": cov_ok,
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

    # ------------------------------------------------------------------
    # Rule D (#768): Protamine-diluted germline TF qualification
    # ------------------------------------------------------------------
    ca_only = ts[
        ts["classification"].eq("Cancer Antigen only")
        & (
            pd.to_numeric(ts["normal_reference_tpm"], errors="coerce")
            <= ab.GERMLINE_TF_MAX_SOMATIC_TPM
        )
        & (
            pd.to_numeric(ts["reproductive_max_sample_tpm"], errors="coerce")
            >= ab.GERMLINE_TF_MIN_BULK_TPM
        )
    ].copy()
    hsc = pd.read_csv(
        reference_data.ensure("hpa_single_cell_tissue", allow_download=False), sep="\t"
    )
    gtex_testis_p = reference_data.pinned_path("gtex_individual_tissue_tpm", "v11-testis")
    want_ca_gids = {g for s in ca_only["member_gene_ids"] for g in json.loads(s)}
    testis_donors_ge3: dict[str, int] = {}
    with gzip.open(gtex_testis_p, "rt") as fh:
        fh.readline()
        fh.readline()
        samples = pd.Index(fh.readline().rstrip("\n").split("\t")[2:])
        donors = samples.str.extract(r"^(GTEX-[^-]+)-.+$", expand=False)
        for line in fh:
            gid = line.partition("\t")[0].split(".")[0]
            if gid in want_ca_gids:
                vec = np.asarray(line.rstrip("\n").split("\t")[2:], dtype=float)
                dmed = pd.Series(vec, index=samples).groupby(donors, sort=False).median()
                testis_donors_ge3[gid] = int((dmed >= 3.0).sum())

    repro_sc_tissues = {"testis", "ovary", "placenta", "endometrium"}
    for r in ca_only.itertuples(index=False):
        eid = str(r.entity_id)
        gid = str(r.representative_gene_id)
        sub_sc = hsc[hsc["Gene"].eq(gid)]
        germ_sc = sub_sc[sub_sc["Tissue"].str.lower().isin(repro_sc_tissues)]
        som_sc = sub_sc[~sub_sc["Tissue"].str.lower().isin(repro_sc_tissues)]
        germ_ntpm = float(germ_sc["nTPM"].max()) if not germ_sc.empty else 0.0
        som_ntpm = float(som_sc["nTPM"].max()) if not som_sc.empty else 0.0
        n_don = testis_donors_ge3.get(gid, 0)
        ev_tf = ab.evaluate_germline_tf_qualification(
            baseline_classification=str(r.classification),
            normal_reference_tpm=float(r.normal_reference_tpm),
            germline_max_tpm=float(r.reproductive_max_sample_tpm),
            germline_donors_ge3_tpm=n_don,
            germline_sc_ntpm=germ_ntpm,
            somatic_sc_ntpm=som_ntpm,
        )
        adj_rows.append(
            {
                "entity_id": eid,
                "protein_group_label": str(r.protein_group_label),
                "rule_family": "germline_tf_protamine_dilution",
                "baseline_classification": str(r.classification),
                "metric_primary": round(float(r.reproductive_max_sample_tpm), 6),
                "metric_secondary": float(n_don),
                "evidence_detail": (
                    f"N={r.normal_reference_tpm:.4f} TPM, germline_max={r.reproductive_max_sample_tpm:.3f} TPM, "
                    f"gtex_testis_donors_ge3={n_don}, hpa_germline_sc_ntpm={germ_ntpm:.2f}, somatic_sc_ntpm={som_ntpm:.2f}"
                ),
                "action_applied": bool(ev_tf["qualifies"]),
                "disposition": ev_tf["disposition"],
            }
        )

    adj_df = (
        pd.DataFrame(adj_rows)
        .sort_values(
            ["rule_family", "action_applied", "protein_group_label"], ascending=[True, False, True]
        )
        .reset_index(drop=True)
    )

    # ------------------------------------------------------------------
    # #769: 1-mm Hamming near-self & I==L tryptic-MS paralog audit
    # Compute BOTH baseline (2026-11-20.2) and integrated Unsafe subtraction!
    # ------------------------------------------------------------------
    seq_tbl = sequence_bundle.sequences_table()
    nuc_cleared_unsafe_eids = set(
        rare_cell_compartment_qc.nuclear_spillover_adjudications(only_cleared=True)
        .loc[lambda df: df["baseline_classification"].eq("Unsafe"), "entity_id"]
        .astype(str)
    )
    ab_cleared_unsafe_eids = set(
        adj_df.loc[
            adj_df["action_applied"] & adj_df["baseline_classification"].eq("Unsafe"),
            "entity_id",
        ].astype(str)
    )
    cleared_all_unsafe_eids = nuc_cleared_unsafe_eids | ab_cleared_unsafe_eids

    newly_unsafe_735_711_eids = set(
        unmeasured_home_panels.unmeasured_home_protein_adjudications(only_ca_cta_gated=True)
        .loc[
            lambda df: df["home_panel_disposition"].eq("home_expressed_above_tier_ceiling"),
            "entity_id",
        ]
        .astype(str)
    ) | set(
        normal_adjudications.rare_cell_candidate_adjudications()
        .loc[
            lambda df: df["attribution_status"].eq("confounder_rejected"),
            "entity_id",
        ]
        .astype(str)
    )

    baseline_unsafe_eids = set(ts.loc[ts["category"].eq("Unsafe"), "entity_id"].astype(str))
    integrated_unsafe_eids = (
        baseline_unsafe_eids - cleared_all_unsafe_eids
    ) | newly_unsafe_735_711_eids

    gid_to_label = {
        g: str(r.protein_group_label)
        for r in ts.itertuples(index=False)
        for g in json.loads(r.member_gene_ids)
    }

    def _index_unsafe_set(u_eids: set[str]) -> tuple[set[str], dict[str, set[str]]]:
        u_gids = {
            g
            for s in ts.loc[ts["entity_id"].isin(u_eids), "member_gene_ids"]
            for g in json.loads(s)
        }
        u_rows = seq_tbl[seq_tbl["gene_id"].isin(u_gids)]
        u_9mers: set[str] = set()
        u_tr_il: dict[str, set[str]] = {}
        for g_id, seq in zip(u_rows["gene_id"], u_rows["sequence"]):
            lbl = gid_to_label.get(str(g_id), str(g_id))
            for m in re.finditer(r"[ACDEFGHIKLMNPQRSTVWY]+", str(seq)):
                seg = m.group()
                for i in range(len(seg) - 8):
                    u_9mers.add(seg[i : i + 9])
            for tp in ab.tryptic_peptides(str(seq)):
                u_tr_il.setdefault(tp.replace("I", "L"), set()).add(lbl)
        return u_9mers, u_tr_il

    base_unsafe_9mers, base_unsafe_tryptic_il = _index_unsafe_set(baseline_unsafe_eids)
    int_unsafe_9mers, int_unsafe_tryptic_il = _index_unsafe_set(integrated_unsafe_eids)

    # Candidate proteins: union of baseline non-Unsafe and integrated non-Unsafe!
    cand_eids = (
        set(ts.loc[~ts["category"].eq("Unsafe"), "entity_id"].astype(str)) | cleared_all_unsafe_eids
    )
    cand_sub = ts[ts["entity_id"].isin(cand_eids)].copy()

    seq_by_tid: dict[str, str] = {}
    seq_by_gid: dict[str, list[str]] = {}
    for g_id, t_id, seq in zip(seq_tbl["gene_id"], seq_tbl["transcript_id"], seq_tbl["sequence"]):
        seq_by_tid[str(t_id).split(".")[0]] = str(seq)
        seq_by_gid.setdefault(str(g_id), []).append(str(seq))

    cand_9mers_by_eid: dict[str, set[str]] = {}
    cand_tryptic_by_eid: dict[str, set[str]] = {}
    all_base_exact_unique: set[str] = set()
    all_int_exact_unique: set[str] = set()

    for r in cand_sub.itertuples(index=False):
        eid = str(r.entity_id)
        rep_tids = (
            [t.split(".")[0] for t in json.loads(r.representative_transcripts)]
            if pd.notna(r.representative_transcripts) and str(r.representative_transcripts)
            else []
        )
        seqs = [seq_by_tid[t] for t in rep_tids if t in seq_by_tid]
        if not seqs:
            for g in json.loads(r.member_gene_ids):
                seqs.extend(seq_by_gid.get(g, []))
        kms: set[str] = set()
        tps: set[str] = set()
        for seq in seqs:
            for m in re.finditer(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
                seg = m.group()
                for i in range(len(seg) - 8):
                    kms.add(seg[i : i + 9])
            tps |= ab.tryptic_peptides(seq)
        cand_9mers_by_eid[eid] = kms
        cand_tryptic_by_eid[eid] = tps
        if eid not in baseline_unsafe_eids:
            all_base_exact_unique.update(kms - base_unsafe_9mers)
        if eid not in integrated_unsafe_eids:
            all_int_exact_unique.update(kms - int_unsafe_9mers)

    def _compute_hit_1mm(cand_uniq: set[str], u_9mers: set[str]) -> set[str]:
        by_mask: list[dict[str, list[str]]] = [{} for _ in range(9)]
        for q in cand_uniq:
            for pos in range(9):
                by_mask[pos].setdefault(q[:pos] + q[pos + 1 :], []).append(q)
        hits: set[str] = set()
        for u in u_9mers:
            for pos in range(9):
                ql = by_mask[pos].get(u[:pos] + u[pos + 1 :])
                if ql:
                    hits.update(ql)
        return hits

    base_hit_1mm = _compute_hit_1mm(all_base_exact_unique, base_unsafe_9mers)
    int_hit_1mm = _compute_hit_1mm(all_int_exact_unique, int_unsafe_9mers)

    ns_rows: list[dict] = []
    for r in cand_sub.itertuples(index=False):
        eid = str(r.entity_id)
        kms = cand_9mers_by_eid[eid]
        tps = cand_tryptic_by_eid[eid]

        if eid in baseline_unsafe_eids:
            ex_b, s2_b, pct1_b, sh_tr_b = 0, 0, np.nan, 0
        else:
            ex_set_b = kms - base_unsafe_9mers
            s2_set_b = ex_set_b - base_hit_1mm
            ex_b, s2_b = len(ex_set_b), len(s2_set_b)
            pct1_b = round(100.0 * (ex_b - s2_b) / ex_b, 2) if ex_b else np.nan
            sh_tr_b = sum(1 for p in tps if p.replace("I", "L") in base_unsafe_tryptic_il)

        if eid in integrated_unsafe_eids:
            ex_i, s2_i, pct1_i, sh_tr_i = 0, 0, np.nan, 0
        else:
            ex_set_i = kms - int_unsafe_9mers
            s2_set_i = ex_set_i - int_hit_1mm
            ex_i, s2_i = len(ex_set_i), len(s2_set_i)
            pct1_i = round(100.0 * (ex_i - s2_i) / ex_i, 2) if ex_i else np.nan
            sh_tr_i = sum(1 for p in tps if p.replace("I", "L") in int_unsafe_tryptic_il)

        ns_rows.append(
            {
                "entity_id": eid,
                "protein_group_label": str(r.protein_group_label),
                "total_9mers": len(kms),
                "exact_unique_9mers": ex_b,
                "safe_2mm_9mers": s2_b,
                "unsafe_near_self_1mm_pct": pct1_b,
                "total_tryptic_peptides": len(tps),
                "unsafe_shared_tryptic_peptides": sh_tr_b,
                "exact_unique_9mers_integrated": ex_i,
                "safe_2mm_9mers_integrated": s2_i,
                "unsafe_near_self_1mm_pct_integrated": pct1_i,
                "unsafe_shared_tryptic_peptides_integrated": sh_tr_i,
            }
        )
    ns_df = pd.DataFrame(ns_rows).sort_values("protein_group_label").reset_index(drop=True)

    # Detailed tryptic-MS paralog audit table across all 12 target candidates (against integrated Unsafe)
    ms_audit_labels = (
        "CTAGE1",
        "TGM4",
        "SPATA31A3",
        "MAGEA6",
        "MAGEA2/B",
        "PAGE2",
        "PLAC1",
        "DKKL1",
        "MAGEB6",
        "CTAG2",
        "DUX4",
        "FOXR2",
    )
    ms_rows: list[dict] = []
    for lbl in ms_audit_labels:
        sub_m = cand_sub[cand_sub["protein_group_label"].eq(lbl)]
        if sub_m.empty:
            raise ValueError(f"Missing expected MS audit candidate {lbl!r}")
        eid = str(sub_m.iloc[0]["entity_id"])
        tps = sorted(cand_tryptic_by_eid[eid])
        shared_map: dict[str, list[str]] = {}
        partner_counts: dict[str, int] = {}
        for p in tps:
            il = p.replace("I", "L")
            if il in int_unsafe_tryptic_il:
                plist = sorted(int_unsafe_tryptic_il[il])
                shared_map[p] = plist
                for pl in plist:
                    partner_counts[pl] = partner_counts.get(pl, 0) + 1
        top_partners = [
            f"{k} ({v})"
            for k, v in sorted(partner_counts.items(), key=lambda kv: (-kv[1], kv[0]))[:6]
        ]
        ms_rows.append(
            {
                "entity_id": eid,
                "protein_group_label": lbl,
                "total_tryptic_peptides": len(tps),
                "unsafe_shared_tryptic_peptides": len(shared_map),
                "unsafe_shared_tryptic_pct": round(100.0 * len(shared_map) / len(tps), 2)
                if tps
                else 0.0,
                "top_unsafe_ms_paralogs": "; ".join(top_partners),
                "example_shared_peptides": json.dumps(dict(list(shared_map.items())[:8])),
            }
        )
    ms_df = pd.DataFrame(ms_rows)

    _write_csv_with_receipt(
        adj_df,
        "assay-boundary-adjudications",
        {
            "disposition_counts": adj_df["disposition"].value_counts().to_dict(),
        },
    )
    _write_csv_with_receipt(
        ns_df,
        "candidate-near-self-1mm",
        {
            "baseline_unsafe_9mers_indexed": len(base_unsafe_9mers),
            "integrated_unsafe_9mers_indexed": len(int_unsafe_9mers),
            "candidate_unique_9mers_evaluated": len(all_int_exact_unique),
        },
    )
    _write_csv_with_receipt(
        ms_df,
        "tryptic-ms-paralog-audit",
        {
            "audited_candidates": len(ms_df),
        },
    )

    obs_path = DATA_DIR / "assay-boundary-comparison-observations.json"
    obs_path.write_text(
        json.dumps(
            {
                "schema": "oncoref.assay-boundary-comparison-observations.v1",
                "policy_id": ab.ASSAY_BOUNDARY_POLICY_ID,
                "comparison_id": COMPARISON_ID,
                "candidates": candidates_payload,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    ns_idx = ns_df.set_index("protein_group_label")
    ms_idx = ms_df.set_index("protein_group_label")
    audit_md = f"""# Assay-Boundary Adjudications & 1-Mismatch Near-Self Audit (#768, #769, #779)

Policy ID: `{ab.ASSAY_BOUNDARY_POLICY_ID}`

## 1. Assay-Boundary Adjudications (`#768`, `#779`)

1. **`#713` Non-admixed tumor witness gating (`TCGA-TGCT` Note `R` libraries):**
   `11` post-meiotic spermatid markers whose `>= 10x` tumor evidence rested solely on `library_qc` Note `R` parenchymal-admixture `TCGA-TGCT` libraries (`CETN1`, `PGK2`, `PRM2`, `BOD1L2`, `CCDC70`, `TSSK2`, `PRR30`, `ZNRF4`, `IQCF5`, `BPIFA3`, `CYLC2`) reclassify from `CTA -> Reproductive Antigen only`, whereas **`PRM1` stays `CTA`** via its independent female `BRCA` witness (`58.97 TPM >= 57.11 TPM` under `#764` bleed QC) and **`SPATA19` stays `CTA`** via its `#737` `cSCC` somatic tumor witness (`GSM3567709 = 13.28 TPM >= 10.0 TPM`, matching `CANCER_ONLY` MS evidence).
2. **RNA-zero single-antibody HPA IHC reconciliation (`cta-ihc-unreliable.csv`):**
   **`MAGEB6`** (`0.0118 TPM` somatic normal max, `0.0 nTPM`, `100%` required RNA coverage complete, `22.14 TPM` Testis, `10.44 TPM` `MM`) advances `Unsafe -> CTA`, while `NLRP9` (`2.32 TPM` somatic bulk) and `XAGE2` remain `Unsafe`.
3. **`Adult cornea` (`GSE155683`) multi-UMI cell support floor (`>= 10` cells with `>= 2` UMIs):**
   **`DKKL1`** (`8` multi-UMI cells across `21,343` adult cornea cells; `3.39 TPM` normal max, `1,346.2 TPM` Testis, `61.62 TPM` `BRCA`) advances `Unsafe -> CTA`, while `INSL4` (`64`), `GPHB5` (`24`), `QRFP` (`40`), and `C20orf202` (`42`) remain `Unsafe`.
4. **Protamine-diluted germline transcription-factor qualification:**
   **`FOXR2`** (`135` Testis donors `>= 3 TPM`, `5.2 nTPM` spermatogonia), **`SPATA31A3`** (`23` Testis donors `>= 3 TPM`, `1.2 nTPM` spermatids), and **`PPP4R3C`** (`230` Testis donors `>= 3 TPM`, `22.7 nTPM` germ cells) advance `Cancer Antigen only -> CTA`, while single-donor `DEFB125` (`1` donor) remains `Cancer Antigen only`.
5. **Droplet ambient milk soup & `<= 100 kb` cis-chromatin spillover (`#779`):**
   **`CSN2`** (`ENSG00000135222` on `4q13.3`; `breast c-3` macrophages `683.2 -> 0.0 nTPM` after `single_cell_tissues(qc="all", cleanup="all")`; `salivary gland` serous clusters `30.2 kb` from `STATH`/`HTN3`/`HTN1` at `321,577–896,409 nTPM` with spillover ratio `<= 0.000050 <= 1e-4`, HPA surgical salivary bulk `0.90 TPM < 1.0 TPM`, salivary IHC `Not detected`, `0` CELLxGENE rare-cell flags, `2.99 TPM` clean somatic normal max, `28,460.4 TPM` lactating breast, `99,761.8 TPM` `BRCA`) advances `Unsafe -> CTA` (`CTA-MF`), while `CSN1S1` (`16.9 nTPM` brain neurons), `CSN3` (`261.4 nTPM` + `395.9 TPM` bulk), `ANKRD30A` (`4.9 nTPM` brain neurons), and `MAGEA3` remain strictly `Unsafe`.

## 2. 1-Mismatch Hamming Near-Self & Tryptic-MS Paralog Audit (`#769`)

- **`MAGEA6`:** `{int(ns_idx.loc["MAGEA6", "total_9mers"])}` total 9-mers -> `{int(ns_idx.loc["MAGEA6", "exact_unique_9mers"])}` exact-unique -> **`{int(ns_idx.loc["MAGEA6", "safe_2mm_9mers"])}` safe `>= 2-mm`** (`{ns_idx.loc["MAGEA6", "unsafe_near_self_1mm_pct"]:.2f}%` 1-mm near-self against `Unsafe` `MAGEA3`/`MAGEA12`).
- **`MAGEA2/B`:** `{int(ns_idx.loc["MAGEA2/B", "total_9mers"])}` total -> `{int(ns_idx.loc["MAGEA2/B", "exact_unique_9mers"])}` exact-unique -> **`{int(ns_idx.loc["MAGEA2/B", "safe_2mm_9mers"])}` safe `>= 2-mm`** (`{ns_idx.loc["MAGEA2/B", "unsafe_near_self_1mm_pct"]:.2f}%` 1-mm near-self).
- **`PAGE2`:** `{int(ns_idx.loc["PAGE2", "total_9mers"])}` total -> `{int(ns_idx.loc["PAGE2", "exact_unique_9mers"])}` exact-unique -> **`{int(ns_idx.loc["PAGE2", "safe_2mm_9mers"])}` safe `>= 2-mm`** (`{ns_idx.loc["PAGE2", "unsafe_near_self_1mm_pct"]:.2f}%` 1-mm near-self).
- **`CTAGE1` (`ENSG00000212710`) vs. `MIA2` (`ENSG00000150527`):** `{int(ms_idx.loc["CTAGE1", "unsafe_shared_tryptic_peptides"])} / {int(ms_idx.loc["CTAGE1", "total_tryptic_peptides"])}` (`{ms_idx.loc["CTAGE1", "unsafe_shared_tryptic_pct"]:.2f}%`) canonical tryptic peptides (`I == L`) are shared directly with ubiquitous `Unsafe` `{ms_idx.loc["CTAGE1", "top_unsafe_ms_paralogs"]}`, explaining healthy multi-organ bottom-up MS attributions without requiring `CTAGE1`-specific translation.

### Summary of Controversial Target Dispositions & Near-Self Footprints

| Protein | Promoted Category | Sex Subset | Exact Unique 9-mers | Safe `>= 2-mm` 9-mers | `1-mm` Near-Self % | Shared `Unsafe` Tryptic Peptides | Key Rule / Note |
|---|---|---|---:|---:|---:|---:|---|
| `PLAC1` | `CTA` | `CTA-M` | `{int(ns_idx.loc["PLAC1", "exact_unique_9mers_integrated"])}` | `{int(ns_idx.loc["PLAC1", "safe_2mm_9mers_integrated"])}` | `{ns_idx.loc["PLAC1", "unsafe_near_self_1mm_pct_integrated"]:.2f}%` | `{int(ns_idx.loc["PLAC1", "unsafe_shared_tryptic_peptides_integrated"])}` | `#765` nuclear 5'-UTR intron & `PABIR2` readthrough cleared |
| `DKKL1` | `CTA` | `CTA-MF` | `{int(ns_idx.loc["DKKL1", "exact_unique_9mers_integrated"])}` | `{int(ns_idx.loc["DKKL1", "safe_2mm_9mers_integrated"])}` | `{ns_idx.loc["DKKL1", "unsafe_near_self_1mm_pct_integrated"]:.2f}%` | `{int(ns_idx.loc["DKKL1", "unsafe_shared_tryptic_peptides_integrated"])}` | `#768` sparse adult cornea (`8 / 21,343` multi-UMI cells) cleared |
| `MAGEB6` | `CTA` | `CTA-M` | `{int(ns_idx.loc["MAGEB6", "exact_unique_9mers_integrated"])}` | `{int(ns_idx.loc["MAGEB6", "safe_2mm_9mers_integrated"])}` | `{ns_idx.loc["MAGEB6", "unsafe_near_self_1mm_pct_integrated"]:.2f}%` | `{int(ns_idx.loc["MAGEB6", "unsafe_shared_tryptic_peptides_integrated"])}` | `#768` single-antibody IHC with `0.012 TPM` somatic RNA cleared |
| `CSN2` | `CTA` | `CTA-MF` | `{int(ns_idx.loc["CSN2", "exact_unique_9mers_integrated"])}` | `{int(ns_idx.loc["CSN2", "safe_2mm_9mers_integrated"])}` | `{ns_idx.loc["CSN2", "unsafe_near_self_1mm_pct_integrated"]:.2f}%` | `{int(ns_idx.loc["CSN2", "unsafe_shared_tryptic_peptides_integrated"])}` | `#779` breast ambient milk soup + `4q13.3` SCPP cis-spillover cleared |
| `CTAG2` | `CTA` | `CTA-MF` | `{int(ns_idx.loc["CTAG2", "exact_unique_9mers_integrated"])}` | `{int(ns_idx.loc["CTAG2", "safe_2mm_9mers_integrated"])}` | `{ns_idx.loc["CTAG2", "unsafe_near_self_1mm_pct_integrated"]:.2f}%` | `{int(ns_idx.loc["CTAG2", "unsafe_shared_tryptic_peptides_integrated"])}` | `#766` `LAGE-1a` (`ENST00000369585`) `SLLMWITQC` clinical clearance |
| `DUX4` | `CTA` | `CTA-E` | `{int(ns_idx.loc["DUX4", "exact_unique_9mers_integrated"])}` | `{int(ns_idx.loc["DUX4", "safe_2mm_9mers_integrated"])}` | `{ns_idx.loc["DUX4", "unsafe_near_self_1mm_pct_integrated"]:.2f}%` | `{int(ns_idx.loc["DUX4", "unsafe_shared_tryptic_peptides_integrated"])}` | `#765 + #767` symmetric `D4Z4` repeat-array silence |
| `FOXR2` | `CTA` | `CTA-M` | `{int(ns_idx.loc["FOXR2", "exact_unique_9mers_integrated"])}` | `{int(ns_idx.loc["FOXR2", "safe_2mm_9mers_integrated"])}` | `{ns_idx.loc["FOXR2", "unsafe_near_self_1mm_pct_integrated"]:.2f}%` | `{int(ns_idx.loc["FOXR2", "unsafe_shared_tryptic_peptides_integrated"])}` | `#768` protamine-diluted spermatogonial TF (`135` Testis donors) |
| `SPATA31A3` | `CTA` | `CTA-M` | `{int(ns_idx.loc["SPATA31A3", "exact_unique_9mers_integrated"]):,}` | `{int(ns_idx.loc["SPATA31A3", "safe_2mm_9mers_integrated"]):,}` | `{ns_idx.loc["SPATA31A3", "unsafe_near_self_1mm_pct_integrated"]:.2f}%` | `{int(ns_idx.loc["SPATA31A3", "unsafe_shared_tryptic_peptides_integrated"])}` | `#768` germline TF; shares `14` tryptic peptides with `SPATA31C1/C2` |
| `PRM1` | `CTA` | `CTA-M` | `{int(ns_idx.loc["PRM1", "exact_unique_9mers_integrated"])}` | `{int(ns_idx.loc["PRM1", "safe_2mm_9mers_integrated"])}` | `{ns_idx.loc["PRM1", "unsafe_near_self_1mm_pct_integrated"]:.2f}%` | `{int(ns_idx.loc["PRM1", "unsafe_shared_tryptic_peptides_integrated"])}` | `#764`/`#776` Thyroid & opposite-sex bleed cleaned (`5.71 TPM`); female `BRCA` `58.97 TPM` |
| `PRM2` | `Reproductive Antigen only` | `MRA-only` | `{int(ns_idx.loc["PRM2", "exact_unique_9mers_integrated"])}` | `{int(ns_idx.loc["PRM2", "safe_2mm_9mers_integrated"])}` | `{ns_idx.loc["PRM2", "unsafe_near_self_1mm_pct_integrated"]:.2f}%` | `{int(ns_idx.loc["PRM2", "unsafe_shared_tryptic_peptides_integrated"])}` | `#768`/`#776` `TCGA-TGCT` Note `R` parenchymal admixture demoted; female-safe `MRA-only` |
| `SPATA19` | `CTA` | `CTA-M` | `{int(ns_idx.loc["SPATA19", "exact_unique_9mers_integrated"])}` | `{int(ns_idx.loc["SPATA19", "safe_2mm_9mers_integrated"])}` | `{ns_idx.loc["SPATA19", "unsafe_near_self_1mm_pct_integrated"]:.2f}%` | `{int(ns_idx.loc["SPATA19", "unsafe_shared_tryptic_peptides_integrated"])}` | `#768 + #737 + #776` somatic `cSCC` witness (`13.28 TPM`) + female breast `0.41 TPM` |
| `MAGEA6` | `CTA` | `CTA-MF` | `{int(ns_idx.loc["MAGEA6", "exact_unique_9mers_integrated"])}` | `{int(ns_idx.loc["MAGEA6", "safe_2mm_9mers_integrated"])}` | `{ns_idx.loc["MAGEA6", "unsafe_near_self_1mm_pct_integrated"]:.2f}%` | `{int(ns_idx.loc["MAGEA6", "unsafe_shared_tryptic_peptides_integrated"])}` | `#769` high 1-mm cross-reactivity risk vs. `Unsafe` `MAGEA3`/`MAGEA12` |
| `CTAGE1` | `CTA` | `CTA-M` | `{int(ns_idx.loc["CTAGE1", "exact_unique_9mers_integrated"])}` | `{int(ns_idx.loc["CTAGE1", "safe_2mm_9mers_integrated"])}` | `{ns_idx.loc["CTAGE1", "unsafe_near_self_1mm_pct_integrated"]:.2f}%` | `{int(ns_idx.loc["CTAGE1", "unsafe_shared_tryptic_peptides_integrated"])}` | `#769` `55/138` tryptic peptides shared with `Unsafe` `MIA2` |

## Programmatic Inspection

```python
from oncoref import assay_boundary_adjudications, peptide_subtraction, target_selection

adj = assay_boundary_adjudications.assay_boundary_adjudications()
near_self = assay_boundary_adjudications.candidate_near_self_1mm()
ms_audit = assay_boundary_adjudications.tryptic_ms_paralog_audit()
policy_1mm = assay_boundary_adjudications.hamming_1mm_policy()
```
"""
    (AUDIT_DIR / "index.md").write_text(audit_md, encoding="utf-8")
    print(
        "Built assay-boundary-adjudications, candidate-near-self-1mm, tryptic-ms-paralog-audit, and audit doc."
    )


if __name__ == "__main__":
    main()

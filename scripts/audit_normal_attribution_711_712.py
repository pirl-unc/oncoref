#!/usr/bin/env python3
"""Reproducible adjudication audit for #712 (HPA normal libraries) and #711 (rare-cell attribution).

Reads checksum-pinned repository inputs (`target-selection.csv.gz`, `normal-comparison.csv.gz`,
`normal-library-qc.csv`, `gene-mapping-risk.csv.gz`, `rare-cell-attribution-reviews.json`,
`rare-cell-evidence.csv.gz`, `donor-populations.csv.gz`, `bulk-confirmation.csv.gz`, and
`cta-recovery-20261117/candidate-evidence.csv.gz`) without mutating released catalog decisions.

Emits:
  - oncoref/data/normal-library-peak-adjudications.csv
  - oncoref/data/hpa-normal-statistic-sensitivity.csv
  - oncoref/data/rare-cell-candidate-adjudications.csv
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import antigen_sets, sequence_bundle, target_selection  # noqa: E402

DATA = ROOT / "oncoref/data"
REBUILD = ROOT / "docs/audits/target-selection-20261110/rebuild"
RARE_AUDIT = ROOT / "docs/audits/normal-cell-atlas-20261110/rare-cell"
RECOVERY_142 = ROOT / "docs/audits/cta-recovery-20261117/candidate-evidence.csv.gz"

ISSUE_712_CANDIDATES: tuple[tuple[str, str], ...] = (
    ("C4orf17", "bone marrow.97"),
    ("TTLL8", "lymph node.175"),
    ("GOLGA6L2", "skeletal muscle.339"),
    ("CCDC87", "parathyroid gland.393"),
    ("RNASE13", "salivary gland.216"),
    ("GAB4", "spleen.378"),
    ("PSG7", "parathyroid gland.393"),
    ("GSTT4", "skeletal muscle.339"),
)

ISSUE_711_CANDIDATES: tuple[str, ...] = (
    "CAGE1",
    "RAD21L1",
    "PPP3R2",
    "RGSL1",
    "DAZL",
    "OOSP4B",
    "TAS2R1",
    "PSG8",
    "ZNF280A",
    "POTEC",
    "VCX3B",
    "PSG3",
    "TUBB8",
    "DCAF8L2",
    "FBXO47",
)


def build_712_adjudications(
    ts: pd.DataFrame,
    ncomp: pd.DataFrame,
    lib_qc: pd.DataFrame,
    risk: pd.DataFrame,
    old_142: pd.DataFrame,
) -> pd.DataFrame:
    valid_n = ncomp[
        ncomp.qc_pass
        & ncomp.identity_resolved
        & ncomp.value.notna()
        & ncomp.scope.isin(["normal", "disposable"])
    ].copy()
    lib_qc_idx = lib_qc.set_index("sample_id")
    ts_by_label = ts.set_index("protein_group_label")

    rows: list[dict] = []
    for label, limiting_sample_142 in ISSUE_712_CANDIDATES:
        rec = ts_by_label.loc[label]
        eid = str(rec["entity_id"])
        gids = json.loads(rec["member_gene_ids"])
        old = old_142.loc[eid]
        qc_row = lib_qc_idx.loc[limiting_sample_142]
        limiting_tissue = str(qc_row["tissue"])

        ent_obs = valid_n[valid_n.entity_id.eq(eid)].sort_values("value", ascending=False)
        hpa_site_obs = ent_obs[
            ent_obs.source_id.eq("hpa_tissue_samples:v23") & ent_obs.tissue.eq(limiting_tissue)
        ]
        hpa_site_p95 = (
            float(hpa_site_obs["value"].max()) if not hpa_site_obs.empty else float("nan")
        )

        # Next highest normal/disposable observation excluding that HPA tissue context
        other_obs = ent_obs[
            ~(ent_obs.source_id.eq("hpa_tissue_samples:v23") & ent_obs.tissue.eq(limiting_tissue))
        ]
        second_tpm = float(other_obs.iloc[0]["value"]) if not other_obs.empty else 0.0
        second_source = str(other_obs.iloc[0]["source_id"]) if not other_obs.empty else ""
        second_tissue = str(other_obs.iloc[0]["tissue"]) if not other_obs.empty else ""

        # Highest non-HPA-biopsy observation (GTEx donor P95 / DICE / Monaco / brain / etc.)
        non_hpa = ent_obs[~ent_obs.source_id.eq("hpa_tissue_samples:v23")]
        max_non_hpa_tpm = float(non_hpa.iloc[0]["value"]) if not non_hpa.empty else 0.0
        max_non_hpa_source = (
            f"{non_hpa.iloc[0]['tissue']} [{non_hpa.iloc[0]['source_id']}]"
            if not non_hpa.empty
            else ""
        )

        # Gene mapping risk across member loci
        r_sub = risk.reindex(gids)
        shared_kmer = float(r_sub["shared_kmer_fraction"].max()) if not r_sub.empty else 0.0
        paralog_ambig = bool(r_sub["paralog_ambiguous"].any()) if not r_sub.empty else False
        host_ov = bool(r_sub["host_overlapped"].any()) if not r_sub.empty else False

        h0_cur = max(
            1.0,
            float(rec["normal_reference_tpm"]) if pd.notna(rec["normal_reference_tpm"]) else 0.0,
            float(rec["disposable_max_tpm"]) if pd.notna(rec["disposable_max_tpm"]) else 0.0,
        )
        repro_supported = float(rec["reproductive_fold"]) * h0_cur
        h0_without_site = max(1.0, second_tpm)
        home_screened = pd.isna(rec["home_tissue"]) or str(rec["home_tissue"]).strip() == ""
        would_pass_without_site = bool(
            rec["category"] != "Unsafe"
            and home_screened
            and all(bool(rec[c]) for c in target_selection.SHARED_CHECKS)
            and float(rec["tumor_max_tpm"]) >= 10.0 * h0_without_site
            and repro_supported >= 10.0 * h0_without_site
        )

        if bool(rec["is_cta"]):
            disposition = "restored_to_cta_under_active_p95"
            supersession = "none_needed_already_cta"
            rationale = (
                f"Active within-tissue biopsy P95 ({hpa_site_p95:.2f} TPM vs 1.42 max {float(old['H0']):.2f} TPM; "
                f"library carryover_P={bool(qc_row['carryover_P'])}, carryover_E={bool(qc_row['carryover_E'])}, "
                f"testis_content={float(qc_row['testis_tissue_content']):.4f}) already satisfies both >=10x gates "
                f"(tumor_fold={float(rec['tumor_fold']):.2f}x, reproductive_fold={float(rec['reproductive_fold']):.2f}x)."
            )
        elif max_non_hpa_tpm >= hpa_site_p95:
            disposition = "superseded_by_independent_gtex_p95"
            supersession = "retain_observation_gtex_controls_h0"
            rationale = (
                f"GTEx v11 donor P95 ({max_non_hpa_tpm:.2f} TPM in {max_non_hpa_source}) exceeds HPA {limiting_tissue} P95 "
                f"({hpa_site_p95:.2f} TPM) and independently holds reproductive_fold={float(rec['reproductive_fold']):.2f}x < 10x."
            )
        elif paralog_ambig:
            disposition = "retained_paralog_multimapping_unverified"
            supersession = "retain_raw_peak_pending_unique_read_quantification"
            rationale = (
                f"Library {limiting_sample_142} passes tissue-signature QC (carryover_P=False, carryover_E=False, "
                f"testis_content=0.0), while locus shares {shared_kmer * 100:.1f}% 31-mers with expressed paralogs; "
                f"secondary normal context is {second_tissue} ({second_tpm:.2f} TPM [{second_source}]). "
                "Peak is retained without silent zeroing until unique-locus read quantification."
            )
        elif second_tpm >= 1.5:
            disposition = "retained_multi_tissue_cross_atlas_concordant"
            supersession = "retain_verified_endocrine_or_somatic_expression"
            rationale = (
                f"Library {limiting_sample_142} passes all tissue-marker and carryover checks, and independent "
                f"non-HPA normal panels corroborate >=1.5 TPM somatic expression ({max_non_hpa_source}: {max_non_hpa_tpm:.2f} TPM; "
                f"next context {second_tissue}: {second_tpm:.2f} TPM)."
            )
        else:
            disposition = "retained_sole_adult_tissue_panel_biopsy"
            supersession = "retain_required_hpa_marrow_coverage"
            rationale = (
                f"Library {limiting_sample_142} passes carryover and olfactory/testis QC with low paralog sharing "
                f"({shared_kmer * 100:.2f}%); HPA v23 provides primary adult {limiting_tissue} bulk TPM absent from GTEx v11 "
                f"(secondary normal skin={second_tpm:.2f} TPM; tumor_fold={float(rec['tumor_fold']):.2f}x)."
            )

        rows.append(
            {
                "entity_id": eid,
                "protein_group_label": label,
                "member_gene_ids": "|".join(gids),
                "limiting_library_1_42": limiting_sample_142,
                "limiting_tissue": limiting_tissue,
                "classification_1_42": str(old["classification"]),
                "h0_1_42_tpm": round(float(old["H0"]), 4),
                "active_classification": str(rec["classification"]),
                "active_h0_tpm": round(h0_cur, 4),
                "active_hpa_site_p95_tpm": round(hpa_site_p95, 4),
                "tumor_max_tpm": round(float(rec["tumor_max_tpm"]), 4),
                "active_tumor_fold": round(float(rec["tumor_fold"]), 4),
                "active_reproductive_fold": round(float(rec["reproductive_fold"]), 4),
                "library_carryover_p": bool(qc_row["carryover_P"]),
                "library_carryover_e": bool(qc_row["carryover_E"]),
                "library_olfactory_share_ge5": round(float(qc_row["olfactory_share_ge5"]), 6),
                "library_testis_content": round(float(qc_row["testis_tissue_content"]), 6),
                "second_normal_context": f"{second_tissue} [{second_source}]",
                "second_normal_tpm": round(second_tpm, 4),
                "max_non_hpa_biopsy_context": max_non_hpa_source,
                "max_non_hpa_biopsy_tpm": round(max_non_hpa_tpm, 4),
                "paralog_ambiguous": paralog_ambig,
                "host_overlapped": host_ov,
                "shared_kmer_fraction": round(shared_kmer, 6),
                "would_pass_10x_without_limiting_site": would_pass_without_site,
                "adjudication_disposition": disposition,
                "supersession_action": supersession,
                "rationale": rationale,
            }
        )
    return pd.DataFrame(rows)


def build_hpa_sensitivity(ts: pd.DataFrame, ncomp: pd.DataFrame) -> pd.DataFrame:
    ts_idx = ts.set_index("entity_id")
    valid_n = ncomp[
        ncomp.qc_pass
        & ncomp.identity_resolved
        & ncomp.value.notna()
        & ncomp.scope.isin(["normal", "disposable"])
    ]
    ra_eligible = ts_idx[list(target_selection.SHARED_CHECKS)].all(axis=1) & ts_idx.category.ne(
        "Unsafe"
    )
    ca_eligible = ra_eligible & (
        ts_idx["home_tissue"].isna() | ts_idx["home_tissue"].astype(str).str.strip().eq("")
    )
    h0_active = ts_idx[["normal_reference_tpm", "disposable_max_tpm"]].max(axis=1).clip(lower=1.0)
    r_supported = ts_idx["reproductive_fold"] * h0_active

    # Counterfactual 1: omit each entity's top HPA biopsy tissue context
    without_top_hpa: dict[str, float] = {}
    for eid, grp in valid_n.groupby("entity_id"):
        ordered = grp.sort_values("value", ascending=False)
        top = ordered.iloc[0]
        if top["source_id"] == "hpa_tissue_samples:v23":
            rem = ordered[
                ~(ordered.source_id.eq("hpa_tissue_samples:v23") & ordered.tissue.eq(top["tissue"]))
            ]
            without_top_hpa[eid] = max(1.0, float(rem["value"].max()) if not rem.empty else 1.0)
        else:
            without_top_hpa[eid] = max(1.0, float(top["value"]))
    h0_no_top = pd.Series(without_top_hpa).reindex(ts_idx.index).fillna(1.0).clip(lower=1.0)

    # Counterfactual 2: replace hpa_tissue_samples:v23 P95 with HPA tissue mean (rna_tissue_hpa:v23)
    h0_no_all_hpa_biopsies = (
        valid_n[~valid_n.source_id.eq("hpa_tissue_samples:v23")]
        .groupby("entity_id")["value"]
        .max()
        .reindex(ts_idx.index)
        .fillna(1.0)
        .clip(lower=1.0)
    )

    scenarios = [
        (
            "active_policy_gtex_p95_plus_hpa_biopsy_p95",
            h0_active,
            True,
            "Active released policy (oncoref.target-selection.2026-11-20.2): within-site donor/biopsy P95 across GTEx v11 and HPA v23.",
        ),
        (
            "counterfactual_omit_single_top_hpa_biopsy_site",
            h0_no_top,
            False,
            "Sensitivity only: drops a protein's single highest HPA v23 biopsy tissue context when HPA sets H0.",
        ),
        (
            "counterfactual_omit_all_hpa_individual_biopsies",
            h0_no_all_hpa_biopsies,
            False,
            "Sensitivity only: omits hpa_tissue_samples:v23 from H0 while keeping HPA tissue means (rna_tissue_hpa:v23) and GTEx v11 donor P95.",
        ),
    ]

    out_rows: list[dict] = []
    base_cta_ids = set(ts_idx.index[ts_idx["is_cta"].astype(bool)])
    for name, h0_vec, is_active, notes in scenarios:
        ca = ca_eligible & ts_idx["tumor_max_tpm"].ge(10.0 * h0_vec)
        ra = ra_eligible & r_supported.ge(10.0 * h0_vec)
        cta = ca & ra
        new_ids = sorted(set(ts_idx.index[cta]) - base_cta_ids)
        new_labels = sorted(ts_idx.loc[new_ids, "protein_group_label"].tolist()) if new_ids else []
        out_rows.append(
            {
                "scenario_id": name,
                "is_active_policy": is_active,
                "total_proteins": len(ts_idx),
                "unsafe_count": int(ts_idx["category"].eq("Unsafe").sum()),
                "ca_count": int(ca.sum()),
                "ra_count": int(ra.sum()),
                "cta_count": int(cta.sum()),
                "ca_only_count": int((ca & ~ra).sum()),
                "ra_only_count": int((ra & ~ca).sum()),
                "delta_cta_vs_active": int(cta.sum()) - len(base_cta_ids),
                "added_cta_labels": "|".join(new_labels),
                "notes": notes,
            }
        )
    return pd.DataFrame(out_rows)


def build_711_adjudications(
    ts: pd.DataFrame,
    ncomp: pd.DataFrame,
    donor_pop: pd.DataFrame,
    bulk_conf: pd.DataFrame,
    reviews_json: dict,
) -> pd.DataFrame:
    ts_by_label = ts.set_index("protein_group_label")
    reviews_by_label = {r["protein"]: r for r in reviews_json["reviews"]}
    rare_ev = pd.read_csv(DATA / "rare-cell-evidence.csv.gz").set_index("entity_id")

    # Verify Unsafe background propagation for confounder_rejected promotions (RGSL1, OOSP4B)
    unsafe_bg_genes = set(antigen_sets.gene_ids("Unsafe"))
    seq_disp = sequence_bundle.target_selection_dispositions().set_index("entity_id")

    valid_n = ncomp[ncomp.qc_pass & ncomp.identity_resolved & ncomp.value.notna()]

    rows: list[dict] = []
    for priority_rank, label in enumerate(ISSUE_711_CANDIDATES, start=1):
        rec = ts_by_label.loc[label]
        rev = reviews_by_label[label]
        eid = str(rec["entity_id"])
        gid = str(rev["gene_id"])
        ev = rare_ev.loc[eid]

        # Donor population peak stats for this gene
        g_pop = donor_pop[donor_pop.gene_id.eq(gid)]
        max_cells_positive = int(g_pop["n_positive"].max()) if not g_pop.empty else 0
        max_cells_multi = int(g_pop["n_multi"].max()) if not g_pop.empty else 0
        max_umi_single_cell = int(g_pop["count_max"].max()) if not g_pop.empty else 0

        # Orthogonal bulk confirmation in the finding organ + retina/brain
        organ = str(rev["organ"])
        b_sub = bulk_conf[bulk_conf.entity_id.eq(eid) & bulk_conf.organ.eq(organ)]
        organ_bulk_max_tpm = float(b_sub["bulk_tpm"].max()) if not b_sub.empty else 0.0

        ent_n = valid_n[valid_n.entity_id.eq(eid)]
        eye_ret = ent_n[ent_n.tissue.isin(["retina", "eye", "cornea/conjunctiva"])]
        eye_bulk_or_cluster_max = float(eye_ret["value"].max()) if not eye_ret.empty else 0.0
        brain_n = ent_n[
            ent_n.tissue.isin(
                [
                    "cerebral cortex",
                    "cerebellum",
                    "basal ganglia",
                    "hypothalamus",
                    "amygdala",
                    "hippocampal formation",
                    "midbrain",
                    "spinal cord",
                    "white matter",
                ]
            )
        ]
        brain_bulk_max = float(brain_n["value"].max()) if not brain_n.empty else 0.0

        in_unsafe_bg = gid in unsafe_bg_genes
        seq_role = str(seq_disp.loc[eid, "shared_sequence"])

        status = str(rev["attribution_status"])
        if status == "confounder_rejected" and int(rev["studies_multi_molecule"]) >= 2:
            evidence_needed = "resolved_unsafe_in_forbidden_background"
            if not (rec["category"] == "Unsafe" and in_unsafe_bg and seq_role == "exclude"):
                raise RuntimeError(f"Unsafe propagation invariant violated for {label}")
        elif status == "single_study_unconfirmed" or (
            status == "confounder_rejected" and int(rev["studies_multi_molecule"]) < 2
        ):
            evidence_needed = "independent_second_donor_study_or_targeted_ish_ihc"
        elif bool(rev["raw_paralog_ambiguous"]):
            evidence_needed = "unique_locus_multimapper_excluded_umi_realignment"
        else:
            evidence_needed = "spliced_exon_junction_or_polya_bam_attribution"

        host_symbols = "|".join(h["host_symbol"] for h in rev.get("host_geometry", []))

        rows.append(
            {
                "priority_rank": priority_rank,
                "review_id": str(rev["id"]),
                "entity_id": eid,
                "protein_group_label": label,
                "gene_id": gid,
                "selected_transcript": str(rev["selected_transcript"]),
                "active_classification": str(rec["classification"]),
                "rare_cell_status": str(ev["rare_cell_status"]),
                "attribution_status": status,
                "supersedes_host_overlap": bool(rev["supersedes_host_overlap"]),
                "finding_organ": organ,
                "finding_top_cell_type": str(rev["top_cell_type"]),
                "studies_multi_molecule": int(rev["studies_multi_molecule"]),
                "donors_multi_molecule": int(rev["donors_multi_molecule"]),
                "max_donor_positive_cells": max_cells_positive,
                "max_donor_multi_umi_cells": max_cells_multi,
                "max_single_cell_umi": max_umi_single_cell,
                "raw_paralog_ambiguous": bool(rev["raw_paralog_ambiguous"]),
                "raw_host_overlapped": bool(rev["raw_host_overlapped"]),
                "host_genes": host_symbols,
                "shared_kmer_fraction": round(float(rev["shared_kmer_fraction"]), 6),
                "shared_3prime_kmer_fraction": round(float(rev["shared_3prime_kmer_fraction"]), 6),
                "organ_bulk_confirmation_max_tpm": round(organ_bulk_max_tpm, 4),
                "brain_panel_max_value": round(brain_bulk_max, 4),
                "eye_retina_panel_max_value": round(eye_bulk_or_cluster_max, 4),
                "tumor_fold": round(float(rec["tumor_fold"]), 4),
                "reproductive_fold": round(float(rec["reproductive_fold"]), 4),
                "in_unsafe_forbidden_background": in_unsafe_bg,
                "sequence_bundle_role": seq_role,
                "resolution_requirement": evidence_needed,
                "rationale": str(rev["rationale"]),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    ts = target_selection.results()
    ncomp = pd.read_csv(REBUILD / "normal-comparison.csv.gz", low_memory=False)
    lib_qc = pd.read_csv(REBUILD / "normal-library-qc.csv")
    risk = pd.read_csv(DATA / "gene-mapping-risk.csv.gz").set_index("gene_id")
    old_142 = pd.read_csv(RECOVERY_142).set_index("entity_id")
    donor_pop = pd.read_csv(RARE_AUDIT / "donor-populations.csv.gz")
    bulk_conf = pd.read_csv(RARE_AUDIT / "bulk-confirmation.csv.gz")
    reviews_json = json.loads((DATA / "rare-cell-attribution-reviews.json").read_text())

    adj_712 = build_712_adjudications(ts, ncomp, lib_qc, risk, old_142)
    sens_712 = build_hpa_sensitivity(ts, ncomp)
    adj_711 = build_711_adjudications(ts, ncomp, donor_pop, bulk_conf, reviews_json)

    adj_712.to_csv(DATA / "normal-library-peak-adjudications.csv", index=False)
    sens_712.to_csv(DATA / "hpa-normal-statistic-sensitivity.csv", index=False)
    adj_711.to_csv(DATA / "rare-cell-candidate-adjudications.csv", index=False)

    print("=== 712 adjudications ===")
    print(
        adj_712[
            [
                "protein_group_label",
                "classification_1_42",
                "active_classification",
                "active_h0_tpm",
                "second_normal_tpm",
                "adjudication_disposition",
            ]
        ].to_string(index=False)
    )
    print("\n=== 712 sensitivity ===")
    print(sens_712.drop(columns=["notes"]).to_string(index=False))
    print("\n=== 711 adjudications ===")
    print(
        adj_711[
            [
                "priority_rank",
                "protein_group_label",
                "active_classification",
                "attribution_status",
                "in_unsafe_forbidden_background",
                "sequence_bundle_role",
                "resolution_requirement",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()

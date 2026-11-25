#!/usr/bin/env python3
# Licensed under the Apache License, Version 2.0.
"""Build per-protein clean tumor witness deltas under tightened gDNA library QC (olfactory_share_ge5 <= 0.025).

Scans all baseline-admitted (`qc_pass == True` & `identity_resolved == True`) tumor samples across
eligible human per-sample linear-TPM cohorts in ``oncoref.source_matrices``, identifies libraries
newly failing ``library_qc.BACKGROUND_MAX_SHARE`` (> 2.5% of measured olfactory receptors at
>= 5 TPM), and recomputes exact QC-passing tumor peak TPM, top sample, top cancer code, and >= 10x
witness counts (both overall and non-admixed non-TGCT for `#768` ablation orthogonality) for every
non-Unsafe human catalog entity whose peak sample or witness set included a failed library.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from oncoref import (  # noqa: E402
    assay_boundary_adjudications,
    library_qc,
    normal_library_qc,
    rare_cell_compartment_qc,
    repeat_family_equivalence,
    source_matrices,
    target_selection,
    tgct_localization,
    unmeasured_home_panels,
)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _summarize_obs(
    obs_map: dict[str, tuple[str, set[str], float]], denom: float, thr: float
) -> dict[str, object]:
    if obs_map:
        best_skey, (_, best_ccs, clean_max) = max(obs_map.items(), key=lambda kv: (kv[1][2], kv[0]))
        best_cc = sorted(best_ccs)[0] if best_ccs else ""
    else:
        best_skey, best_cc, clean_max = "", "", 0.0

    passing = [v for v in obs_map.values() if v[2] >= thr]
    by_pat: dict[str, float] = {}
    anon_n = 0
    pass_ccs: set[str] = set()
    for part_id, obs_ccs, val_f in passing:
        pass_ccs |= obs_ccs
        if part_id:
            by_pat[part_id] = max(by_pat.get(part_id, 0.0), val_f / denom)
        else:
            anon_n += 1

    return {
        "max_tpm": round(clean_max, 6),
        "top_sample": best_skey,
        "top_cancer_code": best_cc,
        "qualifies_tenfold": bool(clean_max >= thr),
        "samples_at_tenfold": len(passing),
        "independent_patients_at_tenfold": len(by_pat),
        "anonymous_samples_at_tenfold": anon_n,
        "cancer_types_at_tenfold": len(pass_ccs),
        "patient_folds_json": json.dumps(
            [round(x, 6) for x in sorted(by_pat.values(), reverse=True)]
        ),
    }


def main() -> int:
    # 1. Load baseline-admitted tumor sample universe (`tumor-samples.csv.gz` & `library-qc.csv.gz` with empty `reasons`)
    # so upstream-quarantined/duplicate/unrelated-tissue/already-failed columns never leak into clean_cohort_gene_mats.
    lqc_df = pd.read_csv(
        ROOT / "docs/audits/target-selection-20261110/rebuild/library-qc.csv.gz",
        keep_default_na=False,
    )
    lqc_pass_pairs = {
        (str(r.cancer_code), str(r.sample_id))
        for r in lqc_df.itertuples(index=False)
        if str(r.reasons).strip() == ""
    }
    tsamp_df = pd.read_csv(
        ROOT / "docs/audits/target-selection-20261110/rebuild/tumor-samples.csv.gz",
        keep_default_na=False,
    )
    codes_by_key: dict[str, set[str]] = {}
    pid_by_key: dict[str, str] = {}
    sample_by_bare_id: dict[str, tuple[str, str, set[str]]] = {}
    valid_sids_by_cohort: dict[str, set[str]] = {}
    for r in tsamp_df.itertuples(index=False):
        sc = str(r.source_cohort)
        ns = str(r.namespace)
        sid = str(r.sample_id)
        skey = str(r.sample_key) or f"{ns}:{sid}"
        cc = str(r.cancer_code)
        pid = str(r.participant_id).strip()
        is_valid = (cc, sid) in lqc_pass_pairs
        if is_valid and sid:
            if sc:
                valid_sids_by_cohort.setdefault(sc, set()).add(sid)
            if cc:
                valid_sids_by_cohort.setdefault(cc, set()).add(sid)
        for k in (skey, f"{sc}:{sid}", f"{ns}:{sid}"):
            if cc:
                codes_by_key.setdefault(k, set()).add(cc)
            if pid:
                pid_by_key[k] = pid
        if sid and is_valid:
            prev_ns, prev_pid, prev_ccs = sample_by_bare_id.get(sid, (ns or sc, pid, set()))
            if cc:
                prev_ccs.add(cc)
            sample_by_bare_id[sid] = (prev_ns, prev_pid or pid, prev_ccs)

    # 2. Scan all cached linear-TPM cohorts in source_matrices restricted to valid_sids_by_cohort
    cov = pd.read_csv(
        ROOT / "docs/audits/target-selection-20261110/rebuild/tumor-source-coverage.csv",
        keep_default_na=False,
    )
    excluded_scale_codes = set(
        cov.loc[cov["status"].str.contains("Excluded", na=False), "cancer_code"].astype(str)
    )
    reg = source_matrices.registry()

    failed_by_cohort: dict[str, list[dict]] = {}
    failed_bare_sids: set[str] = set()
    clean_cohort_gene_mats: dict[str, tuple[str, pd.DataFrame]] = {}

    for row in reg.itertuples(index=False):
        code = str(row.cancer_code)
        if code in excluded_scale_codes or code == "NUTM":
            continue
        allowed_sids = valid_sids_by_cohort.get(code, set())
        if not allowed_sids:
            continue
        p = source_matrices.local_path(code)
        if p is None or not Path(p).is_file():
            continue
        mat = pd.read_parquet(p)
        num_cols = [
            c
            for c in mat.columns
            if pd.api.types.is_numeric_dtype(mat[c]) and str(c) in allowed_sids
        ]
        if not num_cols:
            continue
        sig = library_qc.signatures(mat, num_cols)
        fail_mask = (sig["olfactory_receptors_measured"] >= library_qc.BACKGROUND_MIN_RECEPTORS) & (
            sig["olfactory_share_ge5"] > library_qc.BACKGROUND_MAX_SHARE
        )
        if not fail_mask.any():
            continue
        fail_rows = sig.loc[fail_mask]
        failed_by_cohort[code] = [
            {
                "sample_id": str(sid),
                "olfactory_receptors_measured": int(r.olfactory_receptors_measured),
                "olfactory_share_ge5": round(float(r.olfactory_share_ge5), 6),
            }
            for sid, r in fail_rows.iterrows()
        ]
        failed_bare_sids.update(fail_rows.index.astype(str))
        keep_cols = [c for c in num_cols if c not in fail_rows.index]
        if keep_cols:
            gmat = mat[keep_cols].apply(pd.to_numeric, errors="coerce")
            gmat.index = mat["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
            clean_cohort_gene_mats[code] = (
                code,
                gmat.groupby(level=0).sum(min_count=1),
            )

    # 3. Load frozen 2026-11-20.2 baseline + integrated normal N adjustments + witness ledgers
    base_df = pd.read_csv(
        ROOT / "tests/fixtures/release-1.50-inputs/oncoref/data/target-selection.csv.gz",
        compression="gzip",
        low_memory=False,
    )
    cleared_obs: dict[str, dict] = {}
    cleared_obs.update(
        rare_cell_compartment_qc.cleared_candidate_comparison_observations()["candidates"]
    )
    cleared_obs.update(
        repeat_family_equivalence.repeat_family_comparison_observations()["candidates"]
    )
    cleared_obs.update(
        assay_boundary_adjudications.assay_boundary_comparison_observations()["candidates"]
    )

    bleed_by_eid = (
        normal_library_qc.multiplex_bleed_gene_deltas(only_applied=True)
        .set_index("entity_id")
        .to_dict("index")
    )
    uh_by_eid = (
        unmeasured_home_panels.unmeasured_home_protein_adjudications(only_ca_cta_gated=True)
        .set_index("entity_id")
        .to_dict("index")
    )
    tgct_restr = tgct_localization.tgct_admixture_restricted_entity_ids(include_multi_r_tgct=True)
    tgct_sup_by_eid = (
        tgct_localization.tgct_aafl_protein_support().set_index("entity_id").to_dict("index")
    )
    tgct_r_samples = set(
        tgct_localization.tgct_reproductive_localization_audit(r_only=True)["sample_id"].astype(str)
    )

    tw_df = pd.read_csv(
        ROOT / "docs/audits/target-selection-20261110/rebuild/tumor-witnesses.csv.gz",
        usecols=[
            "entity_id",
            "source_id",
            "sample_id",
            "participant_id",
            "value",
            "qc_pass",
            "identity_resolved",
        ],
        keep_default_na=False,
    )
    tw_df = tw_df[
        tw_df["qc_pass"].astype(str).eq("True") & tw_df["identity_resolved"].astype(str).eq("True")
    ].copy()
    tw_df["value"] = pd.to_numeric(tw_df["value"], errors="coerce")
    tw_df = tw_df[tw_df["value"].notna()]

    tw_all_by_eid: dict[str, list[tuple[str, str, str, tuple[str, ...], float]]] = {}
    for r in tw_df.itertuples(index=False):
        skey = f"{r.source_id}:{r.sample_id}"
        ccs = tuple(sorted(codes_by_key.get(skey, set())))
        tw_all_by_eid.setdefault(str(r.entity_id), []).append(
            (str(r.source_id), str(r.sample_id), str(r.participant_id).strip(), ccs, float(r.value))
        )

    obs_737 = pd.read_csv(ROOT / "oncoref/data/human-737-tumor-observations.csv.gz")
    obs_737_by_eid: dict[str, list[dict]] = {}
    for r in obs_737.to_dict("records"):
        obs_737_by_eid.setdefault(str(r["entity_id"]), []).append(r)

    delta_rows: list[dict] = []
    for r in base_df.itertuples(index=False):
        eid = str(r.entity_id)
        is_cleared = eid in cleared_obs
        if not is_cleared and (r.category == "Unsafe" or not bool(r.required_coverage_complete)):
            continue

        if is_cleared:
            cp = cleared_obs[eid]
            if not bool(cp.get("required_coverage_complete")):
                continue
            n_val = float(cp.get("normal_reference_tpm") or 0.0)
            d_val = float(cp.get("disposable_max_tpm") or 0.0)
        else:
            bleed_r = bleed_by_eid.get(eid)
            n_val = (
                float(bleed_r["clean_normal_reference_tpm"])
                if bleed_r is not None
                else (float(r.normal_reference_tpm) if pd.notna(r.normal_reference_tpm) else np.nan)
            )
            d_val = float(r.disposable_max_tpm) if pd.notna(r.disposable_max_tpm) else 0.0
            uh_r = uh_by_eid.get(eid)
            if (
                uh_r is not None
                and bool(uh_r.get("home_bulk_measured"))
                and pd.notna(uh_r.get("home_bulk_p95_tpm"))
            ):
                hb = float(uh_r["home_bulk_p95_tpm"])
                if str(uh_r.get("home_tissue_tier")) == "disposable":
                    d_val = max(d_val, hb)
                elif pd.notna(n_val):
                    n_val = max(n_val, hb)

        if pd.isna(n_val):
            continue
        denom = max(float(n_val), float(d_val), float(target_selection.DETECTION_FLOOR))
        thr = float(target_selection.ENRICHMENT_FOLD) * denom

        old_top_s = str(r.tumor_top_sample).strip() if pd.notna(r.tumor_top_sample) else ""
        old_bare_s = old_top_s.split(":", 1)[-1] if old_top_s else ""
        tw_list = tw_all_by_eid.get(eid, [])
        had_bad_tw = any(item[1] in failed_bare_sids for item in tw_list)
        had_bad_tgct_sup = (
            eid in tgct_restr
            and eid in tgct_sup_by_eid
            and str(tgct_sup_by_eid[eid].get("top_non_tgct_sample") or "").split(":", 1)[-1]
            in failed_bare_sids
        )
        had_bad_bleed_sup = (
            eid in bleed_by_eid
            and str(bleed_by_eid[eid].get("top_non_tgct_sample") or "").split(":", 1)[-1]
            in failed_bare_sids
        )
        had_bad_cleared = is_cleared and any(
            str(t.get("sample_id") or "").split(":", 1)[-1] in failed_bare_sids
            for t in cleared_obs[eid].get("top_tumor_observations", [])
        )
        if not (
            old_bare_s in failed_bare_sids
            or had_bad_tw
            or had_bad_tgct_sup
            or had_bad_bleed_sup
            or had_bad_cleared
        ):
            continue

        # Track both (a) pure gDNA-clean observations (preserving TGCT for orthogonal #768 ablations)
        # and (b) non-admixed non-TGCT gDNA-clean observations (for #768 + gDNA QC joint evaluation).
        clean_obs_all: dict[str, tuple[str, set[str], float]] = {}
        clean_obs_non_tgct: dict[str, tuple[str, set[str], float]] = {}

        def _rec_into(
            target_map: dict[str, tuple[str, set[str], float]],
            skey: str,
            pid: str,
            ccs: set[str],
            val_f: float,
        ) -> None:
            if skey in target_map:
                p_prev, c_prev, v_prev = target_map[skey]
                target_map[skey] = (
                    pid or p_prev,
                    c_prev | {c for c in ccs if c},
                    max(v_prev, val_f),
                )
            else:
                target_map[skey] = (pid, {c for c in ccs if c}, val_f)

        def _rec(
            skey: str,
            pid: str,
            ccs: set[str],
            val_f: float,
            *,
            is_tgct_note_r: bool = False,
            _all: dict[str, tuple[str, set[str], float]] = clean_obs_all,
            _non_tgct: dict[str, tuple[str, set[str], float]] = clean_obs_non_tgct,
            _is_tgct_restr: bool = eid in tgct_restr,
        ) -> None:
            _rec_into(_all, skey, pid, ccs, val_f)
            if not (_is_tgct_restr and is_tgct_note_r):
                _rec_into(_non_tgct, skey, pid, ccs, val_f)

        for src_id, samp_id, part_id, ccs_tup, val_f in tw_list:
            if samp_id in failed_bare_sids:
                continue
            is_tgct_r = samp_id in tgct_r_samples or "TGCT" in ccs_tup or "TGCT" in src_id
            _rec(f"{src_id}:{samp_id}", part_id, set(ccs_tup), val_f, is_tgct_note_r=is_tgct_r)

        if (
            not is_cleared
            and old_bare_s
            and old_bare_s not in failed_bare_sids
            and pd.notna(r.tumor_max_tpm)
        ):
            skey_old = old_top_s if ":" in old_top_s else f"baseline_tumor:{old_top_s}"
            old_ccs = set(codes_by_key.get(skey_old, set()))
            if pd.notna(r.tumor_top_cancer_code) and str(r.tumor_top_cancer_code).strip():
                old_ccs.add(str(r.tumor_top_cancer_code).strip())
            is_old_tgct_r = old_bare_s in tgct_r_samples or "TGCT" in old_ccs or "TGCT" in skey_old
            _rec(
                skey_old,
                pid_by_key.get(skey_old, ""),
                old_ccs,
                float(r.tumor_max_tpm),
                is_tgct_note_r=is_old_tgct_r,
            )

        if eid in tgct_restr and eid in tgct_sup_by_eid and not had_bad_tgct_sup:
            tr_sup = tgct_sup_by_eid[eid]
            nt_val = (
                float(tr_sup["max_non_tgct_tpm"]) if pd.notna(tr_sup["max_non_tgct_tpm"]) else 0.0
            )
            nt_raw = (
                str(tr_sup["top_non_tgct_sample"]).strip()
                if pd.notna(tr_sup["top_non_tgct_sample"])
                else ""
            )
            nt_cc = (
                str(tr_sup["top_non_tgct_cancer_code"]).strip()
                if pd.notna(tr_sup["top_non_tgct_cancer_code"])
                else ""
            )
            if nt_val > 0.0 and nt_raw:
                skey_nt = nt_raw if ":" in nt_raw else f"somatic_witness_768:{nt_raw}"
                _rec(
                    skey_nt,
                    pid_by_key.get(skey_nt, f"somatic:{nt_raw}"),
                    {nt_cc} if nt_cc else set(),
                    nt_val,
                    is_tgct_note_r=False,
                )

        for r737 in obs_737_by_eid.get(eid, []):
            sid737 = str(r737["sample_id"]).split(":", 1)[-1]
            if sid737 in failed_bare_sids:
                continue
            skey737 = f"{r737['source_cohort']}:{r737['sample_id']}"
            p737 = str(r737["participant_id"]).strip() if pd.notna(r737["participant_id"]) else ""
            _rec(skey737, p737, {str(r737["cancer_code"])}, float(r737["value"]))

        if is_cleared:
            for ct in cleared_obs[eid].get("top_tumor_observations", []):
                ct_sid = str(ct.get("sample_id") or "").split(":", 1)[-1]
                if ct_sid in failed_bare_sids:
                    continue
                skey_ct = f"{ct['source_id']}:{ct['sample_id']}"
                pct = str(ct.get("participant_id") or "").strip() or pid_by_key.get(skey_ct, "")
                ct_ccs = set(codes_by_key.get(skey_ct, set()))
                if ct.get("cancer_code"):
                    ct_ccs.add(str(ct["cancer_code"]))
                is_ct_tgct = ct_sid in tgct_r_samples or "TGCT" in ct_ccs
                _rec(skey_ct, pct, ct_ccs, float(ct["value"]), is_tgct_note_r=is_ct_tgct)

        # Inspect baseline-admitted clean samples of affected cohorts so sub-threshold clean maxima are exact
        gids = json.loads(r.member_gene_ids)
        for c_code, (_src_name, gmat) in clean_cohort_gene_mats.items():
            if not all(g in gmat.index for g in gids):
                continue
            summed = gmat.loc[gids].sum(axis=0, min_count=len(gids)).dropna()
            if summed.empty:
                continue
            best_sid = str(summed.idxmax())
            candidates_sids = set(summed[summed >= thr].index.astype(str)) | {best_sid}
            for c_sid in candidates_sids:
                v_c = float(summed.at[c_sid])
                if v_c <= 0.0 or c_sid not in sample_by_bare_id:
                    continue
                ns_c, pid_c, ccs_c = sample_by_bare_id[c_sid]
                skey_c = f"{ns_c}:{c_sid}"
                is_c_tgct = c_code == "TGCT" or c_sid in tgct_r_samples or "TGCT" in ccs_c
                _rec(skey_c, pid_c, set(ccs_c) | {c_code}, v_c, is_tgct_note_r=is_c_tgct)

        s_all = _summarize_obs(clean_obs_all, denom, thr)
        s_ntgct = _summarize_obs(clean_obs_non_tgct, denom, thr)
        old_max = float(r.tumor_max_tpm) if pd.notna(r.tumor_max_tpm) else 0.0
        delta_rows.append(
            {
                "entity_id": eid,
                "protein_group_label": str(r.protein_group_label),
                "applied": True,
                "previous_tumor_max_tpm": round(old_max, 6),
                "previous_tumor_top_sample": old_top_s,
                "previous_tumor_top_cancer_code": (
                    str(r.tumor_top_cancer_code) if pd.notna(r.tumor_top_cancer_code) else ""
                ),
                "clean_tumor_max_tpm": s_all["max_tpm"],
                "clean_tumor_top_sample": s_all["top_sample"],
                "clean_tumor_top_cancer_code": s_all["top_cancer_code"],
                "single_tumor_threshold_tpm": round(thr, 6),
                "qualifies_tenfold_clean": s_all["qualifies_tenfold"],
                "clean_samples_at_tenfold": s_all["samples_at_tenfold"],
                "clean_independent_patients_at_tenfold": s_all["independent_patients_at_tenfold"],
                "clean_anonymous_samples_at_tenfold": s_all["anonymous_samples_at_tenfold"],
                "clean_cancer_types_at_tenfold": s_all["cancer_types_at_tenfold"],
                "clean_patient_folds_json": s_all["patient_folds_json"],
                "clean_non_tgct_max_tpm": s_ntgct["max_tpm"],
                "clean_non_tgct_top_sample": s_ntgct["top_sample"],
                "clean_non_tgct_top_cancer_code": s_ntgct["top_cancer_code"],
                "clean_non_tgct_samples_at_tenfold": s_ntgct["samples_at_tenfold"],
                "clean_non_tgct_independent_patients_at_tenfold": s_ntgct[
                    "independent_patients_at_tenfold"
                ],
                "clean_non_tgct_anonymous_samples_at_tenfold": s_ntgct[
                    "anonymous_samples_at_tenfold"
                ],
                "clean_non_tgct_cancer_types_at_tenfold": s_ntgct["cancer_types_at_tenfold"],
                "clean_non_tgct_patient_folds_json": s_ntgct["patient_folds_json"],
            }
        )

    out_df = pd.DataFrame(delta_rows).sort_values(["entity_id"]).reset_index(drop=True)
    csv_path = ROOT / "oncoref/data/tumor-gdna-background-deltas.csv"
    json_path = ROOT / "oncoref/data/tumor-gdna-background-deltas.json"
    out_df.to_csv(csv_path, index=False, lineterminator="\n")
    csv_sha = _sha256_bytes(csv_path.read_bytes())

    receipt = {
        "schema": "oncoref.tumor-gdna-background-deltas.v1",
        "policy_id": "oncoref.tumor-gdna-background-qc.2026-11-25.1",
        "background_tpm": library_qc.BACKGROUND_TPM,
        "background_min_receptors": library_qc.BACKGROUND_MIN_RECEPTORS,
        "background_max_share": library_qc.BACKGROUND_MAX_SHARE,
        "failed_samples_total": len(failed_bare_sids),
        "failed_cohorts": {k: len(v) for k, v in sorted(failed_by_cohort.items())},
        "affected_entities": len(out_df),
        "lost_tenfold_entities": int(
            (
                (out_df["previous_tumor_max_tpm"] >= out_df["single_tumor_threshold_tpm"])
                & (~out_df["qualifies_tenfold_clean"])
            ).sum()
        ),
        "sha256": csv_sha,
    }
    json_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in receipt.items() if k != "sha256"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

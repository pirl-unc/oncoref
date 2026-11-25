#!/usr/bin/env python3
"""Promote the integrated human target-selection policy to the active packaged catalog (#775).

Builds:
- ``oncoref/data/baseline-2026-11-20-deltas.json`` (compact delta allowing exact
  reconstruction of the 129-column ``2026-11-20.2`` baseline in ``human_screening_inputs``
  and ablation tests without duplicating a 3.7 MB CSV in the PyPI wheel),
- ``oncoref/data/target-selection.csv.gz`` and ``oncoref/data/target-selection.json``
  (``oncoref.target-selection.2026-11-23.1``: 332 CTA, 117 CA only, 133 RA only,
  1,062 Other, 18,249 Unsafe),
- ``oncoref/data/host-exposure.csv.gz`` and ``oncoref/data/host-exposure.json``
  (adding verified 14-compartment reproductive/embryonic + HPA scRNA/IHC exposures
  for newly cleared ``PLAC1``, ``ZAN``, ``DUX4``, ``DKKL1``, ``MAGEB6``, and ``CSN2``),
- ``oncoref/data/species/human/target-selection.json``,
- ``oncoref/data/antigen-source-inventory.csv`` SHA-256 entries, and
- ``oncoref/data/human-integrated-policy-preview.{csv,json}``.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from build_host_exposure import catalog_compartments, protein_evidence  # noqa: E402

from oncoref import (  # noqa: E402
    assay_boundary_adjudications,
    curation,
    host_context,
    human_screening_inputs,
    normal_library_qc,
    rare_cell_compartment_qc,
    repeat_family_equivalence,
    species_target_selection,
    target_selection,
    tgct_localization,
    unmeasured_home_panels,
)

PROMOTED_POLICY_ID = "oncoref.target-selection.2026-11-23.1"
PROMOTED_COMPARISON_ID = "oncoref.human-grch38-integrated-tpm.2026-11-23.1"
FROZEN_150_DIR = ROOT / "tests/fixtures/release-1.50-inputs"


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _val_equal(a: object, b: object) -> bool:
    if pd.isna(a) and pd.isna(b):
        return True
    if pd.isna(a) or pd.isna(b):
        return False
    if isinstance(a, (float, np.floating)) and isinstance(b, (float, np.floating)):
        return bool(math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12))
    return a == b


def _json_scalar(v: object) -> object:
    if pd.isna(v):
        return None
    if isinstance(v, (np.bool_, bool)):
        return bool(v)
    if isinstance(v, (np.integer, int)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        return float(v)
    return str(v)


def _load_frozen_baseline() -> tuple[pd.DataFrame, dict, pd.DataFrame, dict]:
    base_csv = FROZEN_150_DIR / "oncoref/data/target-selection.csv.gz"
    base_json = FROZEN_150_DIR / "oncoref/data/target-selection.json"
    he_csv = FROZEN_150_DIR / "oncoref/data/host-exposure.csv.gz"
    he_json = FROZEN_150_DIR / "oncoref/data/host-exposure.json"
    return (
        pd.read_csv(base_csv, compression="gzip", low_memory=False),
        json.loads(base_json.read_text(encoding="utf-8")),
        pd.read_csv(he_csv, compression="gzip", keep_default_na=False),
        json.loads(he_json.read_text(encoding="utf-8")),
    )


def _build_cleared_host_exposures(
    cleared_obs_by_eid: dict[str, dict],
    member_genes_by_eid: dict[str, list[str]],
) -> dict[str, dict]:
    """Build 14-compartment host-exposure dicts (with HPA v23 scRNA & IHC) for newly cleared entities."""
    identity = pd.read_csv(
        ROOT / "docs/audits/cta-selection-20261028/hpa_rna_consensus-accepted-identities.csv.gz"
    )
    lookup = identity.set_index("source_gene_id")["canonical_gene_id"]
    gene_to_eid = {g: eid for eid, gids in member_genes_by_eid.items() for g in gids}
    sizes = {eid: len(gids) for eid, gids in member_genes_by_eid.items()}

    out: dict[str, dict] = {}
    for eid, payload in sorted(cleared_obs_by_eid.items()):
        comp_by_name = {
            str(item["compartment"]): item for item in payload.get("reproductive_compartments", [])
        }
        exp_map: dict[str, dict] = {}
        for c in host_context.COMPARTMENTS:
            citem = comp_by_name.get(c)
            if citem is not None and citem.get("max_tpm") is not None:
                val = float(citem["max_tpm"])
                is_emb = c == "embryo"
                srcs = (
                    ["petropoulos_E-MTAB-3929"]
                    if is_emb
                    else ["gtex_individual_tissue_tpm:v11", "hpa_tissue_samples:v23"]
                )
                exp_map[c] = {
                    "coverage": "measured",
                    "max_tpm": val,
                    "measured_observations": 1,
                    "unresolved_observations": 0,
                    "qc_review_observations": 0,
                    "scope_complete": not is_emb,
                    "scope": (
                        "early embryo only; later fetal/developmental stages unavailable"
                        if is_emb
                        else "sampled adult tissue/state; no universal donor or cell-type coverage claim"
                    ),
                    "sources": srcs,
                    "max_ntpm": None,
                    "protein_status": "unassessed",
                }
            else:
                exp_map[c] = {
                    "coverage": "missing",
                    "max_tpm": None,
                    "measured_observations": 0,
                    "unresolved_observations": 0,
                    "qc_review_observations": 0,
                    "scope_complete": False,
                    "scope": (
                        "early embryo only; later fetal/developmental stages unavailable"
                        if c == "embryo"
                        else "sampled adult tissue/state; no universal donor or cell-type coverage claim"
                    ),
                    "sources": [],
                    "max_ntpm": None,
                    "protein_status": "unassessed",
                }
        out[eid] = exp_map

    references = curation.manifest()["reference_inputs"]
    milk_marker_gids = {
        "ENSG00000167531",  # LALBA
        "ENSG00000126545",  # CSN1S1
        "ENSG00000135222",  # CSN2
        "ENSG00000171209",  # CSN3
    }
    from oncoref import human_protein_atlas

    sc_clean_full = human_protein_atlas.single_cell_tissues(qc="all", cleanup="all")
    breast_milk = sc_clean_full[
        sc_clean_full["Tissue"].eq("breast")
        & sc_clean_full["Gene"].isin(milk_marker_gids)
        & (sc_clean_full["nTPM"] >= 100.0)
    ]
    lactating_sc_clusters = set(
        breast_milk.groupby("Cluster")["Gene"].nunique().loc[lambda s: s >= 3].index.astype(str)
    )

    for role in ("host_cell_rna", "host_protein"):
        source, version = references[role]["source"], references[role]["version"]
        if role == "host_cell_rna":
            data = sc_clean_full.copy()
        else:
            path = curation.reference_input(role)
            data = pd.read_csv(path, sep="\t", low_memory=False)
        data["entity_id"] = data["Gene"].map(lookup).map(gene_to_eid)
        data = data[data["entity_id"].isin(out)]
        if data.empty:
            continue
        data["compartment"] = data["Tissue"].map(catalog_compartments)
        data = data.explode("compartment").dropna(subset=["compartment"])
        for (entity, compartment), rows in data.groupby(["entity_id", "compartment"], sort=True):
            if source != "hpa_normal_tissue" and entity == "PEG:480888ecf4b591dc9c7f":
                # #767 replaces transcript-only Kallisto EM (lacking unprocessed_pseudogene DUX4L*
                # decoys) with symmetric genome-aligned D4Z4 quantification on DUX4.
                continue
            item = out[entity][compartment]
            item.setdefault("orthogonal_sources", []).append(f"{source}:{version}")
            if source == "hpa_normal_tissue":
                item.update(protein_evidence(rows))
            else:
                if compartment == "resting_breast" and rows["Tissue"].eq("breast").any():
                    rows = rows[~rows["Cluster"].astype(str).isin(lactating_sc_clusters)]
                keys = [c for c in ("Tissue", "Cluster", "Cell type") if c in rows]
                grouped = rows.groupby(keys).agg(
                    total=("nTPM", "sum"), measured=("Gene", "nunique")
                )
                complete = grouped["measured"].eq(sizes[entity])
                if complete.any():
                    comp_tot = grouped.loc[complete, "total"]
                    item["max_ntpm"] = float(comp_tot.max())
                    if compartment == "resting_breast":
                        item["female_max_ntpm"] = float(comp_tot.max())
                        non_gland = ["breast" not in str(idx[-1]).lower() for idx in comp_tot.index]
                        male_tot = comp_tot[non_gland]
                        item["male_max_ntpm"] = float(male_tot.max()) if not male_tot.empty else 0.0
                if not complete.all() and not grouped.empty:
                    item["scope_complete"] = False
                item["cell_rna_scope"] = (
                    "unlocalized mammary state"
                    if rows["Tissue"].eq("breast").any()
                    else "named tissue cell clusters"
                )
    return out


def _clean_and_sex_stratify_host_exposures(
    merged_he: dict[str, dict],
    base_df: pd.DataFrame,
    cleared_obs_by_eid: dict[str, dict],
) -> dict[str, dict]:
    """Apply deconvolved GTEx host_exposure, sex-scoped breast scRNA/TPM, & RNA-zero IHC reconciliation (#776, #779)."""
    from copy import deepcopy

    from pyensembl import EnsemblRelease

    from oncoref import bulk_qc, gtex, human_protein_atlas
    from oncoref.gtex import deconvolution as gtex_deconv

    ens115 = EnsemblRelease(115)
    coding_tx = {t.id for t in ens115.transcripts() if t.biotype == "protein_coding"}
    hpa_ts = pd.read_csv(curation.reference_input("hpa_tissue_samples"), sep="\t", low_memory=False)
    hpa_coding = (
        hpa_ts[hpa_ts["enstid"].isin(coding_tx)].drop(columns=["enstid"]).groupby("ensgid").sum()
    )
    hpa_tpm_cols = {c.removeprefix("TPM."): c for c in hpa_coding.columns if c.startswith("TPM.")}

    members_by_eid = {
        str(r.entity_id): json.loads(r.member_gene_ids) for r in base_df.itertuples(index=False)
    }

    def _hpa_coding_val(eid: str, sid: object, fallback: float) -> float:
        col = hpa_tpm_cols.get(str(sid))
        if col is None:
            return float(fallback)
        gids = members_by_eid.get(eid, [])
        present = [g for g in gids if g in hpa_coding.index]
        return (
            float(hpa_coding.loc[present, col].sum())
            if (gids and len(present) == len(gids))
            else float(fallback)
        )

    out = deepcopy(merged_he)
    wanted_eids = set(out)
    wanted_gids = sorted({g for eid in wanted_eids for g in members_by_eid.get(eid, [])})

    # 1. Refresh HPA breast & compartment scRNA nTPM using ambient-soup-cleaned single_cell_tissues
    identity = pd.read_csv(
        ROOT / "docs/audits/cta-selection-20261028/hpa_rna_consensus-accepted-identities.csv.gz"
    )
    lookup = identity.set_index("source_gene_id")["canonical_gene_id"]
    gene_to_eid = {g: eid for eid in wanted_eids for g in members_by_eid.get(eid, [])}
    sizes = {eid: len(members_by_eid.get(eid, [])) for eid in wanted_eids}

    sc_clean = human_protein_atlas.single_cell_tissues(qc="all", cleanup="all")
    milk_marker_gids = {
        "ENSG00000167531",
        "ENSG00000126545",
        "ENSG00000135222",
        "ENSG00000171209",
    }
    breast_milk = sc_clean[
        sc_clean["Tissue"].eq("breast")
        & sc_clean["Gene"].isin(milk_marker_gids)
        & (sc_clean["nTPM"] >= 100.0)
    ]
    lactating_sc_clusters = set(
        breast_milk.groupby("Cluster")["Gene"].nunique().loc[lambda s: s >= 3].index.astype(str)
    )
    sc_sub = sc_clean.copy()
    sc_sub["entity_id"] = sc_sub["Gene"].map(lookup).map(gene_to_eid)
    sc_sub = sc_sub[sc_sub["entity_id"].isin(wanted_eids)]
    sc_sub["compartment"] = sc_sub["Tissue"].map(catalog_compartments)
    sc_sub = sc_sub.explode("compartment").dropna(subset=["compartment"])
    for (eid, comp), rows in sc_sub.groupby(["entity_id", "compartment"], sort=False):
        if eid == "PEG:480888ecf4b591dc9c7f":
            continue
        item = out[str(eid)].get(str(comp))
        if not item:
            continue
        if comp == "resting_breast" and rows["Tissue"].eq("breast").any():
            rows = rows[~rows["Cluster"].astype(str).isin(lactating_sc_clusters)]
        grouped = rows.groupby(["Tissue", "Cluster", "Cell type"]).agg(
            total=("nTPM", "sum"), measured=("Gene", "nunique")
        )
        complete = grouped["measured"].eq(sizes[str(eid)])
        if complete.any():
            comp_tot = grouped.loc[complete, "total"]
            item["max_ntpm"] = float(comp_tot.max())
            if comp == "resting_breast":
                item["female_max_ntpm"] = float(comp_tot.max())
                non_gland = ["breast" not in str(idx[-1]).lower() for idx in comp_tot.index]
                male_tot = comp_tot[non_gland]
                item["male_max_ntpm"] = float(male_tot.max()) if not male_tot.empty else 0.0

    qual_comps_by_eid: dict[str, set[str]] = {}
    for r in base_df.itertuples(index=False):
        qraw = getattr(r, "reproductive_qualifying_compartments", "")
        if pd.notna(qraw) and str(qraw).strip():
            qual_comps_by_eid[str(r.entity_id)] = {
                c.strip() for c in str(qraw).split("|") if c.strip()
            }
    for ceid, cp in cleared_obs_by_eid.items():
        qset = qual_comps_by_eid.setdefault(str(ceid), set())
        for citem in cp.get("reproductive_compartments", []):
            if bool(citem.get("qualifies_tenfold")):
                qset.add(str(citem["compartment"]))

    # 2. Compute deconvolved GTEx v11 donor host_exposure & clean max across reproductive & mammary shards
    def _entity_gtex_host_exposure(
        site_name: str,
        *,
        include_samples: set[str] | None = None,
        exclude_samples: set[str] | None = None,
        sex_filter: str | None = None,
    ) -> tuple[dict[str, float], dict[str, float]]:
        sdf = gtex.sample_tpm(site_name, genes=wanted_gids, qc="all", cleanup="all", sex=sex_filter)
        s_cols = [
            c
            for c in sdf.columns
            if c not in ("Name", "Description")
            and (include_samples is None or c in include_samples)
            and (exclude_samples is None or c not in exclude_samples)
        ]
        if not s_cols or sdf.empty:
            return {}, {}
        gids = sdf["Name"].astype(str).str.split(".").str[0]
        mat = sdf[s_cols].copy()
        mat.index = gids
        mat = mat.groupby(level=0).sum(min_count=1)
        pw_map = sdf.attrs.get("purity_weights", {})
        p_wts = np.array([float(pw_map.get(s, 1.0)) for s in s_cols], dtype=np.float64)
        donors = [gtex_deconv.parse_gtex_sample_barcode(s).donor_id for s in s_cols]
        d_idx_map: dict[str, list[int]] = {}
        for j_i, d_id in enumerate(donors):
            d_idx_map.setdefault(d_id, []).append(j_i)
        d_wts = np.array([p_wts[idxs].mean() for idxs in d_idx_map.values()], dtype=np.float64)

        res: dict[str, float] = {}
        res_max: dict[str, float] = {}
        for eid in wanted_eids:
            mems = members_by_eid.get(eid, [])
            if not mems or not all(m in mat.index for m in mems):
                continue
            s_vec = mat.loc[mems].sum(axis=0, min_count=len(mems)).to_numpy(dtype=np.float64)
            with np.errstate(all="ignore"):
                d_vec = np.array(
                    [
                        float(np.nanmedian(s_vec[idxs]))
                        if np.isfinite(s_vec[idxs]).any()
                        else np.nan
                        for idxs in d_idx_map.values()
                    ],
                    dtype=np.float64,
                )
            fin = np.isfinite(d_vec)
            if not fin.any():
                continue
            rv, rw = d_vec[fin], d_wts[fin]
            if len(rv) >= 50:
                val = bulk_qc.purity_weighted_harrell_davis_quantile(rv, rw, 0.99)
            elif len(rv) >= 20:
                val = bulk_qc.purity_weighted_harrell_davis_quantile(rv, rw, 0.95)
            else:
                val = float(rv.max())
            res[eid] = float(val)
            res_max[eid] = float(rv.max())
        return res, res_max

    rc = pd.read_csv(
        ROOT / "docs/audits/target-selection-20261110/rebuild/reproductive-comparison.csv.gz",
        low_memory=False,
    )
    lactating_sids = set(rc.loc[rc["compartment"].eq("lactating_breast"), "sample_id"].astype(str))
    # Separate coordinated multi-gene milk-operon induction (postpartum lactation in XX women and
    # hyperprolactinemia / galactorrhea in 8/336 XY ICU/medicated cadaver donors) from resting_breast:
    milk_sdf = gtex.sample_tpm(
        "breast_mammary_tissue",
        genes=sorted(milk_marker_gids),
        qc="all",
        cleanup="all",
    )
    milk_scols = [c for c in milk_sdf.columns if c not in ("Name", "Description")]
    milk_mat = milk_sdf[milk_scols]
    secretory_mammary_sids = lactating_sids | {
        c for c in milk_scols if (milk_mat[c] >= 5.0).sum() >= 3 and milk_mat[c].max() >= 15.0
    }
    gtex_he_by_comp: dict[str, dict[str, float]] = {}
    gtex_max_by_comp: dict[str, dict[str, float]] = {}
    for comp_lbl, site_lbl in (
        ("testis", "testis"),
        ("prostate", "prostate"),
        ("ovary", "ovary"),
        ("uterus", "uterus"),
        ("vagina", "vagina"),
        ("fallopian tube", "fallopian_tube"),
    ):
        he_d, mx_d = _entity_gtex_host_exposure(site_lbl)
        gtex_he_by_comp[comp_lbl] = he_d
        gtex_max_by_comp[comp_lbl] = mx_d
    ecto, ecto_mx = _entity_gtex_host_exposure("cervix_ectocervix")
    endo, endo_mx = _entity_gtex_host_exposure("cervix_endocervix")
    gtex_he_by_comp["cervix"] = {
        eid: max(ecto.get(eid, 0.0), endo.get(eid, 0.0)) for eid in (set(ecto) | set(endo))
    }
    gtex_max_by_comp["cervix"] = {
        eid: max(ecto_mx.get(eid, 0.0), endo_mx.get(eid, 0.0))
        for eid in (set(ecto_mx) | set(endo_mx))
    }
    gtex_rb_xx, _ = _entity_gtex_host_exposure(
        "breast_mammary_tissue", exclude_samples=secretory_mammary_sids, sex_filter="XX"
    )
    gtex_rb_xy, _ = _entity_gtex_host_exposure(
        "breast_mammary_tissue", exclude_samples=secretory_mammary_sids, sex_filter="XY"
    )
    gtex_lb_xx, _ = _entity_gtex_host_exposure(
        "breast_mammary_tissue", include_samples=secretory_mammary_sids, sex_filter="XX"
    )
    large_gtex_comps = {"testis", "prostate", "ovary", "uterus", "vagina"}

    # HPA coding biopsy maxima per (eid, compartment)
    hpa_max_by_ec: dict[tuple[str, str], float] = {}
    hpa_rc = rc[
        rc["source_id"].eq("hpa_tissue_samples:v23")
        & rc["identity_resolved"]
        & rc["qc_pass"].astype(bool)
        & rc["library_qc_review"].astype(object).fillna("").eq("")
        & rc["value"].notna()
    ]
    for (eid_raw, comp_raw), grp in hpa_rc.groupby(["entity_id", "compartment"], sort=False):
        eid_s = str(eid_raw)
        if eid_s not in wanted_eids:
            continue
        vals = [_hpa_coding_val(eid_s, s, float(v)) for s, v in zip(grp["sample_id"], grp["value"])]
        if vals:
            c_key = "resting_breast" if str(comp_raw) == "breast" else str(comp_raw)
            hpa_max_by_ec[(eid_s, c_key)] = max(hpa_max_by_ec.get((eid_s, c_key), 0.0), max(vals))

    for eid in wanted_eids:
        e_map = out[eid]
        q_set = qual_comps_by_eid.get(eid, set())
        for comp, g_dict in gtex_he_by_comp.items():
            item = e_map.get(comp)
            if not item or item.get("coverage") != "measured":
                continue
            g_val = g_dict.get(eid)
            g_mx = gtex_max_by_comp.get(comp, {}).get(eid)
            h_val = hpa_max_by_ec.get((eid, comp))
            ntpm = item.get("max_ntpm")
            pstat = item.get("protein_status")
            orth_corrob = (
                comp in q_set
                or (ntpm is not None and float(ntpm) >= 2.0)
                or pstat in {"positive", "unresolved"}
            )
            if g_val is not None:
                base_g = max(g_val, g_mx or 0.0) if comp in q_set else g_val
                if comp in large_gtex_comps:
                    item["max_tpm"] = (
                        max(base_g, h_val) if (h_val is not None and orth_corrob) else base_g
                    )
                else:
                    item["max_tpm"] = max(base_g, h_val) if h_val is not None else base_g
            elif h_val is not None:
                item["max_tpm"] = h_val

        # Also refresh HPA-only compartments (placenta, endometrium, epididymis, seminal vesicle, lactating_breast)
        for h_only_comp in ("placenta", "endometrium", "epididymis", "seminal vesicle"):
            item = e_map.get(h_only_comp)
            h_val = hpa_max_by_ec.get((eid, h_only_comp))
            if item and item.get("coverage") == "measured" and h_val is not None:
                item["max_tpm"] = h_val

        lb_item = e_map.get("lactating_breast")
        if lb_item and lb_item.get("coverage") == "measured":
            g_lb = gtex_lb_xx.get(eid)
            h_lb = hpa_max_by_ec.get((eid, "lactating_breast"))
            if g_lb is not None:
                lb_item["max_tpm"] = max(g_lb, h_lb) if h_lb is not None else g_lb
            elif h_lb is not None:
                lb_item["max_tpm"] = h_lb

        rb_item = e_map.get("resting_breast")
        if rb_item and rb_item.get("coverage") == "measured":
            g_xx = gtex_rb_xx.get(eid)
            g_xy = gtex_rb_xy.get(eid)
            h_rb = hpa_max_by_ec.get((eid, "resting_breast"))
            f_ntpm = rb_item.get("female_max_ntpm", rb_item.get("max_ntpm"))
            pstat = rb_item.get("protein_status")
            orth_rb = (f_ntpm is not None and float(f_ntpm) >= 2.0) or pstat in {
                "positive",
                "unresolved",
            }
            if g_xx is not None:
                rb_item["female_max_tpm"] = (
                    max(g_xx, h_rb) if (h_rb is not None and orth_rb) else g_xx
                )
            elif "female_max_tpm" not in rb_item:
                rb_item["female_max_tpm"] = float(rb_item["max_tpm"])
            if g_xy is not None:
                rb_item["male_max_tpm"] = g_xy
            elif "male_max_tpm" not in rb_item:
                rb_item["male_max_tpm"] = float(rb_item["max_tpm"])
            rb_item["max_tpm"] = max(
                float(rb_item["female_max_tpm"]), float(rb_item["male_max_tpm"])
            )

        # 3. Reconcile single-antibody 'unresolved' HPA IHC in near-zero RNA compartments (#779)
        for _comp, item in e_map.items():
            if not isinstance(item, dict) or item.get("protein_status") != "unresolved":
                continue
            mt = item.get("max_tpm")
            mn = item.get("max_ntpm")
            levels = {str(f.get("Level", "")) for f in (item.get("protein_findings") or [])}
            if (
                mt is not None
                and float(mt) < 1.2
                and (mn is None or float(mn) < 2.0)
                and "High" not in levels
            ):
                item["protein_status"] = "rna_zero_reconciled"

    return out


def main() -> None:
    base_df, base_meta, base_he_df, base_he_meta = _load_frozen_baseline()
    assert len(base_df) == 19893

    # Collect cleared observation payloads for #765, #767, and #768 former-Unsafe/Review rows
    cleared_obs_by_eid: dict[str, dict] = {}
    cleared_obs_by_eid.update(
        rare_cell_compartment_qc.cleared_candidate_comparison_observations()["candidates"]
    )
    cleared_obs_by_eid.update(
        repeat_family_equivalence.repeat_family_comparison_observations()["candidates"]
    )
    cleared_obs_by_eid.update(
        assay_boundary_adjudications.assay_boundary_comparison_observations()["candidates"]
    )

    # Evaluate the full integrated catalog from the frozen 2026-11-20.2 starting state
    # (temporarily monkeypatching _baseline_2026_11_20_frame and packaged_species_host_exposure)
    base_he_dict = base_he_df.set_index("entity_id")["exposures"].map(json.loads).to_dict()
    missing_he_eids = {
        eid: cleared_obs_by_eid[eid] for eid in cleared_obs_by_eid if eid not in base_he_dict
    }
    member_genes_by_eid = {
        str(r.entity_id): json.loads(r.member_gene_ids)
        for r in base_df.loc[base_df["entity_id"].isin(missing_he_eids)].itertuples(index=False)
    }
    extra_he = _build_cleared_host_exposures(missing_he_eids, member_genes_by_eid)
    merged_he = _clean_and_sex_stratify_host_exposures(
        {**base_he_dict, **extra_he}, base_df, cleared_obs_by_eid
    )

    from copy import deepcopy

    target_selection._baseline_2026_11_20_frame = lambda: base_df.copy()  # type: ignore[attr-defined]
    target_selection._baseline_2026_11_20_host_exposures = lambda: deepcopy(base_he_dict)  # type: ignore[attr-defined]
    species_target_selection.packaged_species_host_exposure = lambda species="human", **kw: (
        deepcopy(merged_he)
    )  # type: ignore[assignment]
    species_target_selection._packaged_species_host_exposure_ref = lambda species="human", **kw: (
        merged_he
    )  # type: ignore[assignment]

    base_eval = human_screening_inputs.evaluate_human_catalog(
        include_home_panels=False,
        include_adjudications=False,
        include_737_cohorts=False,
        include_multiplex_bleed_qc=False,
        include_nuclear_spillover_qc=False,
        include_clinical_clearances=False,
        include_repeat_family_equivalence=False,
        include_assay_boundary_fixes=False,
        policy_id="oncoref.target-selection.2026-11-20.2",
    )
    bc = base_eval.catalog.set_index("entity_id")

    integ_res = human_screening_inputs.evaluate_human_catalog(
        include_home_panels=True,
        include_adjudications=True,
        include_737_cohorts=True,
        include_multiplex_bleed_qc=True,
        include_nuclear_spillover_qc=True,
        include_clinical_clearances=True,
        include_repeat_family_equivalence=True,
        include_assay_boundary_fixes=True,
        policy_id=PROMOTED_POLICY_ID,
    )
    ic = integ_res.catalog.set_index("entity_id")

    # Lookup tables for enriching tumor_top_sample, 135-cohort tumor witness counts, and contexts
    obs_737 = human_screening_inputs.upgraded_737_tumor_observations()
    obs_737_by_eid: dict[str, list[dict]] = {}
    for r in obs_737.to_dict("records"):
        obs_737_by_eid.setdefault(str(r["entity_id"]), []).append(r)

    tsamp_df = pd.read_csv(
        ROOT / "docs/audits/target-selection-20261110/rebuild/tumor-samples.csv.gz",
        keep_default_na=False,
    )
    codes_by_key: dict[str, set[str]] = {}
    pid_by_key: dict[str, str] = {}
    sample_by_bare_id: dict[str, tuple[str, str, set[str]]] = {}
    for r in tsamp_df.itertuples(index=False):
        sc = str(r.source_cohort)
        ns = str(r.namespace)
        sid = str(r.sample_id)
        skey = str(r.sample_key) or f"{ns}:{sid}"
        cc = str(r.cancer_code)
        pid = str(r.participant_id).strip()
        for k in (skey, f"{sc}:{sid}", f"{ns}:{sid}"):
            if cc:
                codes_by_key.setdefault(k, set()).add(cc)
            if pid:
                pid_by_key[k] = pid
        if sid:
            prev_ns, prev_pid, prev_ccs = sample_by_bare_id.get(sid, (ns or sc, pid, set()))
            if cc:
                prev_ccs.add(cc)
            sample_by_bare_id[sid] = (prev_ns, prev_pid or pid, prev_ccs)
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
    tw_by_eid: dict[str, list[tuple[str, str, str, tuple[str, ...], float]]] = {}
    for r in tw_df.itertuples(index=False):
        skey = f"{r.source_id}:{r.sample_id}"
        ccs = tuple(sorted(codes_by_key.get(skey, set())))
        tw_by_eid.setdefault(str(r.entity_id), []).append(
            (str(r.source_id), str(r.sample_id), str(r.participant_id).strip(), ccs, float(r.value))
        )
    tgct_r_samples = set(
        tgct_localization.tgct_reproductive_localization_audit(r_only=True)["sample_id"].astype(str)
    )

    bleed_by_eid = (
        normal_library_qc.multiplex_bleed_gene_deltas(only_applied=True)
        .set_index("entity_id")
        .to_dict("index")
    )
    tgct_restr = tgct_localization.tgct_admixture_restricted_entity_ids(include_multi_r_tgct=True)
    tgct_sup_by_eid = (
        tgct_localization.tgct_aafl_protein_support().set_index("entity_id").to_dict("index")
    )
    uh_by_eid = (
        unmeasured_home_panels.unmeasured_home_protein_adjudications(only_ca_cta_gated=True)
        .set_index("entity_id")
        .to_dict("index")
    )

    promoted = base_df.copy()
    promoted["policy"] = PROMOTED_POLICY_ID

    sync_cols = [
        "category",
        "subtype",
        "classification",
        "normal_screen_pass",
        "required_coverage_complete",
        "identity_translation_pass",
        "evidence_resolved",
        "expression_bounds_pass",
        "sequence_pass",
        "is_cancer_antigen",
        "is_reproductive_antigen",
        "is_cta",
        "cta_tier",
        "is_target_candidate",
        "cancer_evidence",
        "reproductive_evidence",
        "normal_reference_tpm",
        "disposable_max_tpm",
        "disposable_max_ntpm",
        "unsafe_sequence_overlap_pct",
        "tumor_max_tpm",
        "tumor_fold",
        "tumor_top_cancer_code",
        "reproductive_max_sample_tpm",
        "reproductive_supported_tpm",
        "reproductive_fold",
        "reproductive_qualifying_compartments",
        "rare_cell_status",
        "rare_cell_action",
        "rare_cell_basis",
        "Unsafe_reasons",
    ]
    ic_aligned = ic.reindex(promoted["entity_id"]).reset_index(drop=True)
    bc_aligned = bc.reindex(promoted["entity_id"]).reset_index(drop=True)
    for col in sync_cols:
        if col not in ic_aligned.columns:
            continue
        for idx in range(len(promoted)):
            iv = ic_aligned.at[idx, col]
            bv = bc_aligned.at[idx, col] if col in bc_aligned.columns else None
            if not _val_equal(iv, bv):
                if isinstance(iv, str) and promoted[col].dtype.kind in ("f", "i", "u", "b"):
                    promoted[col] = promoted[col].astype(object)
                promoted.at[idx, col] = iv

    # Add provenance & sex-restriction columns directly onto promoted table
    extra_cols = [
        "sex_restriction_class",
        "sex_restriction_code",
        "rare_cell_711_requirement",
        "normal_library_712_disposition",
        "home_panel_735_disposition",
        "multiplex_bleed_764_disposition",
        "nuclear_spillover_765_disposition",
        "nuclear_rna_caution",
        "clinical_clearance_766_disposition",
        "repeat_family_767_disposition",
        "assay_boundary_768_disposition",
    ]
    for col in extra_cols:
        promoted[col] = ic_aligned[col].fillna("").astype(str).to_numpy()

    # Reconcile auxiliary display/provenance columns for modified entities
    for idx, row in promoted.iterrows():
        eid = str(row["entity_id"])
        old_row = base_df.iloc[idx]
        cat_now = str(row["category"])
        cls_now = str(row["classification"])

        if cat_now == "Unsafe":
            # Newly Unsafe (#735 SP7/GFY/PANX3/BPIFB3/GLYATL3 or #711 DAZL)
            if str(old_row["category"]) != "Unsafe":
                promoted.at[idx, "normal_screen_pass"] = False
                promoted.at[idx, "decision_reason"] = str(row["Unsafe_reasons"])
                promoted.at[idx, "expression_bounds_status"] = "unassessed"
                for clear_num in (
                    "normal_reference_tpm",
                    "single_tumor_threshold_tpm",
                    "tumor_max_tpm",
                    "tumor_fold",
                    "samples_at_tenfold",
                    "independent_patients_at_tenfold",
                    "anonymous_samples_at_tenfold",
                    "patient_or_sample_witnesses_at_tenfold",
                    "cancer_types_at_tenfold",
                    "reproductive_supported_tpm",
                    "reproductive_fold",
                ):
                    promoted.at[idx, clear_num] = np.nan
                for clear_str in (
                    "independent_patient_tumor_folds",
                    "tumor_top_sample",
                    "tumor_top_cancer_code",
                ):
                    promoted.at[idx, clear_str] = np.nan
            continue

        if bool(row["is_target_candidate"]):
            for chk_col in target_selection.SHARED_CHECKS:
                promoted.at[idx, chk_col] = True
        elif eid in cleared_obs_by_eid:
            promoted.at[idx, "normal_screen_pass"] = True

        n_tpm = row["normal_reference_tpm"]
        d_tpm = row["disposable_max_tpm"]
        t_tpm = row["tumor_max_tpm"]
        r_tpm = row["reproductive_supported_tpm"]
        if pd.notna(n_tpm):
            d_val = float(d_tpm) if pd.notna(d_tpm) else 0.0
            if pd.isna(d_tpm):
                promoted.at[idx, "disposable_max_tpm"] = d_val
            if pd.isna(promoted.at[idx, "disposable_max_ntpm"]):
                promoted.at[idx, "disposable_max_ntpm"] = 0.0
            if pd.isna(promoted.at[idx, "healthy_maximum_ntpm"]):
                promoted.at[idx, "healthy_maximum_ntpm"] = float(n_tpm)
            bar = 10.0 * max(float(n_tpm), d_val, 1.0)
            promoted.at[idx, "single_tumor_threshold_tpm"] = bar
            peaks = [float(x) for x in (t_tpm, r_tpm) if pd.notna(x)]
            peak_val = max(peaks) if peaks else 0.0
            promoted.at[idx, "target_peak_tpm"] = peak_val
            promoted.at[idx, "expression_bounds_status"] = (
                "pass" if bool(row["expression_bounds_pass"]) else "fail"
            )

        # Populate cleared candidate context strings (PLAC1, ZAN, DUX4, DKKL1, MAGEB6, TENT5D)
        if eid in cleared_obs_by_eid and pd.notna(n_tpm):
            cp = cleared_obs_by_eid[eid]
            promoted.at[idx, "normal_reference_contexts"] = json.dumps(
                [
                    {
                        "source_id": "human_normal_reference_20261123",
                        "tissue": "peripheral_normal_reference",
                        "source_context": "peripheral_normal_reference",
                        "value": float(n_tpm),
                    }
                ]
            )
            r_obs = [
                {
                    "source_id": (
                        "human_embryo_atlases"
                        if c["compartment"] == "embryo"
                        else "human_reproductive_panel"
                    ),
                    "tissue": str(c["compartment"]),
                    "source_context": str(c["compartment"]),
                    "value": float(c["max_tpm"]),
                }
                for c in cp.get("reproductive_compartments", [])
                if bool(c.get("qualifies_tenfold"))
            ]
            if r_obs:
                top_r = max(r_obs, key=lambda x: x["value"])
                promoted.at[idx, "reproductive_supported_observations"] = json.dumps([top_r])
                promoted.at[idx, "reproductive_supported_context"] = str(top_r["source_context"])
            exp_ctxs = [
                {
                    "tissue": str(c["compartment"]),
                    "compartment": str(c["compartment"]),
                    "cell_type": "bulk" if c["compartment"] != "embryo" else "embryo_lineage",
                    "stage": "embryo" if c["compartment"] == "embryo" else "adult",
                    "value": float(c["max_tpm"]),
                }
                for c in cp.get("reproductive_compartments", [])
                if float(c.get("max_tpm") or 0.0) >= 1.0 and c["compartment"] != "embryo"
            ]
            promoted.at[idx, "reproductive_expression_contexts"] = json.dumps(exp_ctxs)
            top_tumors = cp.get("top_tumor_observations", [])
            if top_tumors and eid not in tw_by_eid:
                best_t = max(top_tumors, key=lambda x: float(x["value"]))
                promoted.at[idx, "tumor_top_sample"] = (
                    f"{best_t['source_id']}:{best_t['sample_id']}"
                )

        # Update normal_reference_contexts if #764 bleed or #735 home panel changed N
        if eid in bleed_by_eid and pd.notna(n_tpm):
            br = bleed_by_eid[eid]
            promoted.at[idx, "normal_reference_contexts"] = json.dumps(
                [
                    {
                        "source_id": "gtex_donor_p95_bleed_qc_764:v11",
                        "tissue": str(br["top_clean_gtex_site"]),
                        "source_context": str(br["top_clean_gtex_site"]),
                        "value": float(n_tpm),
                    }
                ]
            )
        elif eid in uh_by_eid and pd.notna(n_tpm) and pd.notna(old_row["normal_reference_tpm"]):
            if float(n_tpm) > float(old_row["normal_reference_tpm"]) + 1e-6:
                uh_r = uh_by_eid[eid]
                promoted.at[idx, "normal_reference_contexts"] = json.dumps(
                    [
                        {
                            "source_id": "unmeasured_home_panels_735",
                            "tissue": str(uh_r["home_tissue"]),
                            "source_context": str(uh_r["home_tissue"]),
                            "value": float(n_tpm),
                        }
                    ]
                )

        # Reconcile exact 135-cohort tumor witness/patient counts whenever threshold, #737, #768, or cleared status applies
        d_curr_val = (
            float(promoted.at[idx, "disposable_max_tpm"])
            if pd.notna(promoted.at[idx, "disposable_max_tpm"])
            else 0.0
        )
        d_old_val = (
            float(old_row["disposable_max_tpm"]) if pd.notna(old_row["disposable_max_tpm"]) else 0.0
        )
        needs_witness_recompute = pd.notna(n_tpm) and (
            eid in obs_737_by_eid
            or eid in tgct_restr
            or eid in bleed_by_eid
            or eid in cleared_obs_by_eid
            or pd.isna(old_row["normal_reference_tpm"])
            or abs(float(n_tpm) - float(old_row["normal_reference_tpm"])) > 1e-9
            or abs(d_curr_val - d_old_val) > 1e-9
        )
        if needs_witness_recompute:
            denom = max(float(n_tpm), d_curr_val, 1.0)
            thr = 10.0 * denom
            merged_obs: dict[str, tuple[str, set[str], float]] = {}

            def _record_obs(
                skey: str,
                pid: str,
                ccs: set[str],
                val_f: float,
                obs_map: dict[str, tuple[str, set[str], float]] = merged_obs,
            ) -> None:
                if skey in obs_map:
                    prev_pid, prev_ccs, prev_val = obs_map[skey]
                    obs_map[skey] = (
                        pid or prev_pid,
                        prev_ccs | {c for c in ccs if c},
                        max(prev_val, val_f),
                    )
                else:
                    obs_map[skey] = (pid, {c for c in ccs if c}, val_f)

            for src_id, samp_id, part_id, ccs_tup, val_f in tw_by_eid.get(eid, []):
                if eid in tgct_restr and (
                    samp_id in tgct_r_samples or "TGCT" in ccs_tup or "TGCT" in src_id
                ):
                    continue
                _record_obs(f"{src_id}:{samp_id}", part_id, set(ccs_tup), val_f)
            if (
                eid not in tgct_restr
                and pd.notna(old_row["tumor_max_tpm"])
                and pd.notna(old_row["tumor_top_sample"])
                and str(old_row["tumor_top_sample"]).strip()
            ):
                old_top = str(old_row["tumor_top_sample"]).strip()
                if ":" not in old_top:
                    old_top = f"baseline_tumor:{old_top}"
                old_ccs = set(codes_by_key.get(old_top, set()))
                if (
                    pd.notna(old_row["tumor_top_cancer_code"])
                    and str(old_row["tumor_top_cancer_code"]).strip()
                ):
                    old_ccs.add(str(old_row["tumor_top_cancer_code"]).strip())
                _record_obs(
                    old_top,
                    pid_by_key.get(old_top, ""),
                    old_ccs,
                    float(old_row["tumor_max_tpm"]),
                )
            if eid in tgct_restr and eid in tgct_sup_by_eid:
                tr_sup = tgct_sup_by_eid[eid]
                nt_val = (
                    float(tr_sup["max_non_tgct_tpm"])
                    if pd.notna(tr_sup["max_non_tgct_tpm"])
                    else 0.0
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
                    if ":" in nt_raw:
                        nt_skey = nt_raw
                        nt_pid = pid_by_key.get(nt_skey, f"somatic:{nt_skey}")
                        nt_ccs = set(codes_by_key.get(nt_skey, set()))
                    elif nt_raw in sample_by_bare_id:
                        nt_ns, nt_pid_lookup, cc_lookup = sample_by_bare_id[nt_raw]
                        nt_skey = f"{nt_ns}:{nt_raw}"
                        nt_pid = nt_pid_lookup or f"somatic:{nt_raw}"
                        nt_ccs = set(cc_lookup)
                    else:
                        nt_skey = f"oncoref_non_admixed_tumor_witness_768:{nt_raw}"
                        nt_pid = f"somatic:{nt_raw}"
                        nt_ccs = set()
                    if nt_cc:
                        nt_ccs.add(nt_cc)
                    _record_obs(nt_skey, nt_pid, nt_ccs, nt_val)
            for r737 in obs_737_by_eid.get(eid, []):
                skey = f"{r737['source_cohort']}:{r737['sample_id']}"
                p737 = (
                    str(r737["participant_id"]).strip() if pd.notna(r737["participant_id"]) else ""
                )
                _record_obs(skey, p737, {str(r737["cancer_code"])}, float(r737["value"]))
            if eid in cleared_obs_by_eid:
                for ct in cleared_obs_by_eid[eid].get("top_tumor_observations", []):
                    skey = f"{ct['source_id']}:{ct['sample_id']}"
                    pct = str(ct.get("participant_id") or "").strip() or pid_by_key.get(skey, "")
                    ct_ccs = set(codes_by_key.get(skey, set()))
                    if ct.get("cancer_code"):
                        ct_ccs.add(str(ct["cancer_code"]))
                    _record_obs(skey, pct, ct_ccs, float(ct["value"]))

            if merged_obs:
                best_skey, (_, best_ccs, best_val) = max(
                    merged_obs.items(), key=lambda kv: (kv[1][2], kv[0])
                )
                promoted.at[idx, "tumor_max_tpm"] = best_val
                promoted.at[idx, "tumor_fold"] = best_val / denom
                promoted.at[idx, "tumor_top_sample"] = best_skey
                if best_ccs:
                    orig_cc = (
                        str(old_row["tumor_top_cancer_code"]).strip()
                        if pd.notna(old_row["tumor_top_cancer_code"])
                        else ""
                    )
                    promoted.at[idx, "tumor_top_cancer_code"] = (
                        orig_cc if orig_cc in best_ccs else sorted(best_ccs)[0]
                    )
                r_curr = promoted.at[idx, "reproductive_supported_tpm"]
                promoted.at[idx, "target_peak_tpm"] = max(
                    best_val, float(r_curr) if pd.notna(r_curr) else 0.0
                )

            passing_obs = [v for v in merged_obs.values() if v[2] >= thr]
            by_patient: dict[str, float] = {}
            anon_cnt = 0
            c_codes: set[str] = set()
            for part_id, obs_ccs, val_f in passing_obs:
                c_codes |= obs_ccs
                if part_id:
                    by_patient[part_id] = max(by_patient.get(part_id, 0.0), val_f / denom)
                else:
                    anon_cnt += 1
            folds_sorted = sorted(by_patient.values(), reverse=True)
            promoted.at[idx, "samples_at_tenfold"] = len(passing_obs)
            promoted.at[idx, "independent_patients_at_tenfold"] = len(folds_sorted)
            promoted.at[idx, "anonymous_samples_at_tenfold"] = anon_cnt
            promoted.at[idx, "patient_or_sample_witnesses_at_tenfold"] = (
                len(folds_sorted) + anon_cnt
            )
            promoted.at[idx, "cancer_types_at_tenfold"] = len(c_codes)
            promoted.at[idx, "independent_patient_tumor_folds"] = json.dumps(folds_sorted)

        # Update human-readable decision_reason on changed rows
        if cls_now != str(old_row["classification"]):
            if cls_now == "CTA":
                promoted.at[idx, "decision_reason"] = (
                    "Qualifies both cancer and reproductive/developmental 10x routes under integrated 2026-11-23 policy"
                )
            elif cls_now == "Cancer Antigen only":
                promoted.at[idx, "decision_reason"] = (
                    "Qualifies cancer 10x route after #735/#737 home-tissue and cohort integration"
                )
            elif cls_now == "Reproductive Antigen only":
                promoted.at[idx, "decision_reason"] = (
                    "Qualifies reproductive 10x route; no non-admixed tumor sample reaches 10x"
                )
            elif cls_now == "Other:Expression":
                promoted.at[idx, "decision_reason"] = (
                    "Neither positive route reaches 10x over the screened adult normal baseline"
                )
            elif cls_now == "Other:Review":
                promoted.at[idx, "decision_reason"] = (
                    "Cleared nuclear snRNA-seq spillover (#765); held pending complete peripheral comparison coverage"
                )

    # Compute compact 2026-11-20.2 delta JSON so _baseline_2026_11_20_frame() restores base_df byte-accurately
    base_cols = list(base_df.columns)
    added_cols = [c for c in promoted.columns if c not in base_df.columns]
    row_deltas: dict[str, dict[str, object]] = {}
    for idx in range(len(base_df)):
        eid = str(base_df.at[idx, "entity_id"])
        diff: dict[str, object] = {}
        for col in base_cols:
            if col == "policy":
                continue
            old_v = base_df.at[idx, col]
            new_v = promoted.at[idx, col]
            if not _val_equal(old_v, new_v):
                diff[col] = _json_scalar(old_v)
        if diff:
            row_deltas[eid] = diff

    he_deltas: dict[str, dict[str, dict]] = {}
    for eid, orig_comps in sorted(base_he_dict.items()):
        new_comps = merged_he.get(eid, {})
        c_diff: dict[str, dict] = {}
        for comp, orig_d in orig_comps.items():
            cur_d = new_comps.get(comp)
            if cur_d is None:
                c_diff[comp] = orig_d
                continue
            if orig_d == cur_d:
                continue
            patch = {k: v for k, v in orig_d.items() if k not in cur_d or cur_d[k] != v}
            drop = [k for k in cur_d if k not in orig_d]
            entry: dict[str, object] = {}
            if patch:
                entry["__patch__"] = patch
            if drop:
                entry["__drop__"] = sorted(drop)
            if entry:
                c_diff[comp] = entry
        if c_diff:
            he_deltas[eid] = c_diff

    sb_manifest = json.loads(
        (ROOT / "docs/audits/sequence-bundle-20261121/manifest.json").read_text(encoding="utf-8")
    )
    delta_doc = {
        "schema": "oncoref.baseline-catalog-deltas.v1",
        "baseline_policy_id": "oncoref.target-selection.2026-11-20.2",
        "promoted_policy_id": PROMOTED_POLICY_ID,
        "baseline_package_dependencies_sha256": dict(sb_manifest["package_dependencies"]),
        "added_columns": added_cols,
        "changed_entity_count": len(row_deltas),
        "baseline_host_exposure_entity_ids": sorted(base_he_df["entity_id"].astype(str).tolist()),
        "baseline_host_exposure_deltas": he_deltas,
        "rows": row_deltas,
    }
    delta_path = ROOT / "oncoref/data/baseline-2026-11-20-deltas.json"
    delta_path.write_text(
        json.dumps(delta_doc, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8"
    )

    # Write promoted target-selection.csv.gz
    ts_csv_path = ROOT / "oncoref/data/target-selection.csv.gz"
    promoted.to_csv(
        ts_csv_path,
        index=False,
        compression={"method": "gzip", "mtime": 0},
        lineterminator="\n",
    )
    ts_sha = _sha256_file(ts_csv_path)
    ts_bytes = ts_csv_path.stat().st_size

    # Write promoted host-exposure.csv.gz (930 entities with #776 sex-stratification & bleed QC)
    he_rows = [
        {"entity_id": eid, "exposures": json.dumps(merged_he[eid], allow_nan=False, sort_keys=True)}
        for eid in sorted(merged_he)
    ]
    he_csv_path = ROOT / "oncoref/data/host-exposure.csv.gz"
    pd.DataFrame(he_rows).to_csv(
        he_csv_path, index=False, compression={"method": "gzip", "mtime": 0}
    )
    he_sha = _sha256_file(he_csv_path)
    he_bytes = he_csv_path.stat().st_size

    # Update host-exposure.json
    new_he_meta = dict(base_he_meta)
    new_he_meta["catalog_policy_id"] = PROMOTED_POLICY_ID
    new_he_meta["catalog_data_sha256"] = ts_sha
    new_he_meta["sha256"] = he_sha
    (ROOT / "oncoref/data/host-exposure.json").write_text(
        json.dumps(new_he_meta, indent=2) + "\n", encoding="utf-8"
    )

    # Update target-selection.json
    new_ts_meta = json.loads(json.dumps(base_meta))
    new_ts_meta["policy_id"] = PROMOTED_POLICY_ID
    new_ts_meta["data_sha256"] = ts_sha
    new_ts_meta["data_bytes"] = ts_bytes
    summ = new_ts_meta["summary"]
    summ["policy_id"] = PROMOTED_POLICY_ID
    summ["categories"] = {k: int(v) for k, v in promoted["category"].value_counts().items()}
    summ["classifications"] = {
        k: int(v) for k, v in promoted["classification"].value_counts().items()
    }
    summ["memberships"] = {
        "is_cancer_antigen": int(promoted["is_cancer_antigen"].astype(bool).sum()),
        "is_reproductive_antigen": int(promoted["is_reproductive_antigen"].astype(bool).sum()),
        "is_cta": int(promoted["is_cta"].astype(bool).sum()),
        "is_target_candidate": int(promoted["is_target_candidate"].astype(bool).sum()),
    }
    ctas_only = promoted[promoted["is_cta"].astype(bool)]
    summ["cta_tiers"] = {k: int(v) for k, v in ctas_only["cta_tier"].value_counts().items()}
    shared_ok = ~promoted["category"].eq("Unsafe") & promoted[
        list(target_selection.SHARED_CHECKS)
    ].astype(bool).all(axis=1)
    summ["shared_eligible"] = int(shared_ok.sum())
    summ["assessed_for_cancer_evidence"] = int(promoted["normal_reference_tpm"].notna().sum())
    summ["rare_cell"]["status"] = {
        k: int(v) for k, v in promoted["rare_cell_status"].value_counts().items()
    }
    summ["rare_cell"]["action"] = {
        k: int(v) for k, v in promoted["rare_cell_action"].value_counts().items()
    }
    pub_mask = promoted["published_nomination"].fillna(False).astype(bool)
    summ["published_nomination_by_category"] = {
        k: int(v) for k, v in promoted.loc[pub_mask, "category"].value_counts().items()
    }
    (ROOT / "oncoref/data/target-selection.json").write_text(
        json.dumps(new_ts_meta, indent=2) + "\n", encoding="utf-8"
    )

    # Update species/human/target-selection.json
    sp_human_path = ROOT / "oncoref/data/species/human/target-selection.json"
    sp_meta = json.loads(sp_human_path.read_text(encoding="utf-8"))
    sp_meta["policy_id"] = PROMOTED_POLICY_ID
    sp_meta["comparison_id"] = PROMOTED_COMPARISON_ID
    sp_meta["data_sha256"] = ts_sha
    sp_meta["data_bytes"] = ts_bytes
    sp_meta["host_exposure_sha256"] = he_sha
    sp_meta["host_exposure_bytes"] = he_bytes
    sp_meta["summary"]["categories"] = summ["categories"]
    sp_meta["summary"]["classifications"] = summ["classifications"]
    sp_meta["summary"]["sex_restriction_codes"] = dict(
        sorted(integ_res.summary["sex_restriction_codes"].items())
    )
    sex_df = integ_res.sex_restricted_antigens()
    sp_meta["summary"]["sex_restriction_classes"] = {
        k: int(v) for k, v in sorted(sex_df["sex_restriction_class"].value_counts().items())
    }
    sp_human_path.write_text(json.dumps(sp_meta, indent=2) + "\n", encoding="utf-8")

    # Update SHA-256 entries in antigen-source-inventory.csv
    inv_path = ROOT / "oncoref/data/antigen-source-inventory.csv"
    inv = pd.read_csv(inv_path, keep_default_na=False)
    inv.loc[inv["dataset"].eq("target-selection.csv.gz"), "sha256"] = ts_sha
    inv.loc[inv["dataset"].eq("host-exposure.csv.gz"), "sha256"] = he_sha
    inv.to_csv(inv_path, index=False, lineterminator="\n")

    # Recompute and write human-integrated-policy-preview.{csv,json}
    changed_df, preview_summary = human_screening_inputs.compute_next_policy_transitions(
        integrated_result=integ_res
    )
    prev_csv_path = ROOT / "oncoref/data/human-integrated-policy-preview.csv"
    prev_json_path = ROOT / "oncoref/data/human-integrated-policy-preview.json"
    changed_df.to_csv(prev_csv_path, index=False, lineterminator="\n")
    prev_json_path.write_text(
        json.dumps({**preview_summary, "sha256": _sha256_file(prev_csv_path)}, indent=2) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "policy_id": PROMOTED_POLICY_ID,
                "categories": summ["categories"],
                "sex_restriction_codes": sp_meta["summary"]["sex_restriction_codes"],
                "delta_rows": len(row_deltas),
                "changed_preview_rows": len(changed_df),
                "host_exposure_entities": len(he_rows),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

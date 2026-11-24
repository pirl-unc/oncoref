#!/usr/bin/env python3
"""Promote the integrated human target-selection policy to the active packaged catalog (#775).

Builds:
- ``oncoref/data/baseline-2026-11-20-deltas.json`` (compact delta allowing exact
  reconstruction of the 129-column ``2026-11-20.2`` baseline in ``human_screening_inputs``
  and ablation tests without duplicating a 3.7 MB CSV in the PyPI wheel),
- ``oncoref/data/target-selection.csv.gz`` and ``oncoref/data/target-selection.json``
  (``oncoref.target-selection.2026-11-23.1``: 331 CTA, 117 CA only, 133 RA only,
  1,062 Other, 18,250 Unsafe),
- ``oncoref/data/host-exposure.csv.gz`` and ``oncoref/data/host-exposure.json``
  (adding verified 14-compartment reproductive/embryonic + HPA scRNA/IHC exposures
  for newly cleared ``PLAC1``, ``ZAN``, ``DUX4``, ``DKKL1``, and ``MAGEB6``),
- ``oncoref/data/species/human/target-selection.json``,
- ``oncoref/data/antigen-source-inventory.csv`` SHA-256 entries, and
- ``oncoref/data/human-integrated-policy-preview.{csv,json}``.
"""

from __future__ import annotations

import gzip
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
    reference_data,
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
    for role in ("host_cell_rna", "host_protein"):
        source, version = references[role]["source"], references[role]["version"]
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
                keys = [c for c in ("Tissue", "Cluster", "Cell type") if c in rows]
                grouped = rows.groupby(keys).agg(
                    total=("nTPM", "sum"), measured=("Gene", "nunique")
                )
                complete = grouped["measured"].eq(sizes[entity])
                if complete.any():
                    item["max_ntpm"] = float(grouped.loc[complete, "total"].max())
                if not complete.all():
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
    """Stratify GTEx resting_breast by donor sex and mask #764 bleed on opposite-sex shards (#776)."""
    from copy import deepcopy

    from pyensembl import EnsemblRelease

    med = pd.read_csv(reference_data.pinned_path("gtex_tissue_rna", "v11"), sep="\t", skiprows=2)
    med["gene_id"] = med["Name"].astype(str).str.split(".").str[0]
    sentinel_ref = normal_library_qc.stoichiometric_sentinel_reference(
        med, panels=normal_library_qc.HOST_EXPOSURE_SENTINEL_PANELS
    )

    sex_y = ["ENSG00000129824", "ENSG00000067048", "ENSG00000183878", "ENSG00000012817"]
    wanted_bytes = {g.encode("ascii") for g in (set(sentinel_ref["gene_id"]) | set(sex_y))}

    def _load_gct(pin_key: str) -> pd.DataFrame:
        p = reference_data.pinned_path("gtex_individual_tissue_tpm", pin_key)
        vals: dict[str, np.ndarray] = {}
        with gzip.open(p, "rb") as fh:
            fh.readline()
            fh.readline()
            samples = pd.Index(
                fh.readline().decode("utf-8").rstrip("\n").split("\t")[2:], name="sample_id"
            )
            for line in fh:
                rn, _, rest = line.partition(b"\t")
                bid = rn.split(b".", 1)[0]
                if bid in wanted_bytes:
                    gid = bid.decode("ascii")
                    vals[gid] = vals.get(gid, 0.0) + np.asarray(
                        rest.partition(b"\t")[2].split(), dtype=float
                    )
        return pd.DataFrame.from_dict(vals, orient="index", columns=samples)

    pins = [
        k
        for k in reference_data.REFERENCE_SOURCES["gtex_individual_tissue_tpm"]["pins"]
        if any(
            x in k
            for x in ("breast", "ovary", "uterus", "vagina", "cervix", "fallopian", "prostate")
        )
    ]
    male_xy_samples: set[str] = set()
    tp_bleed_samples: set[str] = set()
    for pk in sorted(pins):
        mat = _load_gct(pk)
        y_sum = mat.reindex(sex_y).fillna(0.0).sum(axis=0)
        if "prostate" not in pk:
            male_xy_samples |= set(y_sum[y_sum >= 5.0].index)
        fl = normal_library_qc.detect_stoichiometric_multiplex_bleed(
            mat,
            site=pk.removeprefix("v11-"),
            sentinel_reference=sentinel_ref,
            max_bleed_fraction=0.05,
        )
        tp_bleed_samples |= set(
            fl[fl["source_key"].isin(["testis", "testis_chromatin", "prostate"])]["sample_id"]
        )

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
    ab_tf_eids = set(
        assay_boundary_adjudications.assay_boundary_adjudications()
        .loc[
            lambda d: d["disposition"].eq("germline_tf_protamine_diluted_qualified"),
            "entity_id",
        ]
        .astype(str)
    )
    male_origin_eids = ab_tf_eids | {
        str(r.entity_id)
        for r in base_df.itertuples(index=False)
        if {
            c.strip()
            for c in str(r.reproductive_qualifying_compartments or "").split("|")
            if c.strip()
        }
        & {"testis", "epididymis", "prostate", "seminal vesicle"}
    }
    for ceid, cp in cleared_obs_by_eid.items():
        if any(
            c.get("qualifies_tenfold")
            and str(c.get("compartment")) in {"testis", "epididymis", "prostate", "seminal vesicle"}
            for c in cp.get("reproductive_compartments", [])
        ):
            male_origin_eids.add(str(ceid))

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
    rc = pd.read_csv(
        ROOT / "docs/audits/target-selection-20261110/rebuild/reproductive-comparison.csv.gz",
        low_memory=False,
    )
    female_gtex_comps = {
        "ovary",
        "uterus",
        "vagina",
        "cervix",
        "fallopian tube",
        "lactating_breast",
    }

    for eid_raw, grp in rc.groupby("entity_id", sort=True):
        eid = str(eid_raw)
        if eid not in out:
            continue
        grp = grp.copy()
        hm = grp["source_id"].eq("hpa_tissue_samples:v23")
        if hm.any():
            grp.loc[hm, "value"] = [
                _hpa_coding_val(eid, s, float(v))
                for s, v in zip(grp.loc[hm, "sample_id"], grp.loc[hm, "value"])
            ]
        for comp in female_gtex_comps:
            item = out[eid].get(comp)
            if not item or item.get("coverage") != "measured":
                continue
            c_rows = grp[
                grp["compartment"].eq(comp) & grp["identity_resolved"] & grp["value"].notna()
            ]
            c_rows = c_rows[~c_rows["sample_id"].isin(male_xy_samples)]
            ntpm, pstat = item.get("max_ntpm"), item.get("protein_status")
            hpa_sub = c_rows[c_rows["source_id"].eq("hpa_tissue_samples:v23")]
            hpa_max = float(hpa_sub["value"].max()) if len(hpa_sub) else 0.0
            if (
                (ntpm is None or float(ntpm) < 0.5)
                and pstat not in {"positive", "unresolved"}
                and hpa_max < 0.5
                and eid in male_origin_eids
            ):
                clean = c_rows[~c_rows["sample_id"].isin(tp_bleed_samples)]
                if len(clean) >= 5:
                    c_rows = clean
            if len(c_rows):
                item["max_tpm"] = float(c_rows["value"].max())

        rb_item = out[eid].get("resting_breast")
        if rb_item and rb_item.get("coverage") == "measured":
            rb_rows = grp[
                grp["compartment"].isin(["resting_breast", "breast"])
                & grp["identity_resolved"]
                & grp["value"].notna()
            ]
            ntpm, pstat = rb_item.get("max_ntpm"), rb_item.get("protein_status")
            hpa_rb = rb_rows[rb_rows["source_id"].eq("hpa_tissue_samples:v23")]
            hpa_rb_max = float(hpa_rb["value"].max()) if len(hpa_rb) else 0.0
            orth_ok = (
                (ntpm is None or float(ntpm) < 0.5)
                and pstat not in {"positive", "unresolved"}
                and hpa_rb_max < 2.0
            )
            fem_rows = rb_rows[~rb_rows["sample_id"].isin(male_xy_samples)]
            male_rows = rb_rows[
                rb_rows["sample_id"].isin(male_xy_samples)
                | ~rb_rows["source_id"].eq("gtex_individual_tissue_tpm:v11")
            ]
            if orth_ok and eid in male_origin_eids:
                fc = fem_rows[~fem_rows["sample_id"].isin(tp_bleed_samples)]
                if len(fc) >= 5:
                    fem_rows = fc
                mc = male_rows[~male_rows["sample_id"].isin(tp_bleed_samples)]
                if len(mc) >= 3:
                    male_rows = mc
            rb_item["female_max_tpm"] = (
                float(fem_rows["value"].max()) if len(fem_rows) else float(rb_item["max_tpm"])
            )
            rb_item["male_max_tpm"] = (
                float(male_rows["value"].max()) if len(male_rows) else float(rb_item["max_tpm"])
            )
            rb_item["max_tpm"] = max(rb_item["female_max_tpm"], rb_item["male_max_tpm"])
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

    # Lookup tables for enriching tumor_top_sample, normal_reference_contexts, reproductive contexts
    obs_737 = human_screening_inputs.upgraded_737_tumor_observations()
    obs_737_by_eid: dict[str, list[dict]] = {}
    for r in obs_737.to_dict("records"):
        obs_737_by_eid.setdefault(str(r["entity_id"]), []).append(r)

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
        "samples_at_tenfold",
        "independent_patients_at_tenfold",
        "independent_patient_tumor_folds",
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
                    "reproductive_supported_tpm",
                    "reproductive_fold",
                ):
                    promoted.at[idx, clear_num] = np.nan
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
            if top_tumors:
                best_t = max(top_tumors, key=lambda x: float(x["value"]))
                promoted.at[idx, "tumor_top_sample"] = (
                    f"{best_t['source_id']}:{best_t['sample_id']}"
                )
                if pd.notna(promoted.at[idx, "single_tumor_threshold_tpm"]):
                    thr = float(promoted.at[idx, "single_tumor_threshold_tpm"])
                    qual_t = [x for x in top_tumors if float(x["value"]) >= thr]
                    promoted.at[idx, "cancer_types_at_tenfold"] = len(
                        {str(x["cancer_code"]) for x in qual_t if x.get("cancer_code")}
                    )
                    anon_n = sum(
                        1 for x in qual_t if not str(x.get("participant_id") or "").strip()
                    )
                    promoted.at[idx, "anonymous_samples_at_tenfold"] = anon_n
                    promoted.at[idx, "patient_or_sample_witnesses_at_tenfold"] = (
                        int(row["independent_patients_at_tenfold"] or 0) + anon_n
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

        # Update tumor_top_sample for #768 TGCT admixture demotions or #737 upgrades
        if eid in tgct_restr and eid in tgct_sup_by_eid:
            tsup = tgct_sup_by_eid[eid]
            if eid in obs_737_by_eid and pd.notna(t_tpm):
                best737 = max(obs_737_by_eid[eid], key=lambda x: float(x["value"]))
                if float(best737["value"]) >= float(tsup.get("max_non_tgct_tpm") or 0.0):
                    promoted.at[idx, "tumor_top_sample"] = (
                        f"{best737['source_cohort']}:{best737['sample_id']}"
                    )
                else:
                    promoted.at[idx, "tumor_top_sample"] = str(
                        tsup.get("top_non_tgct_sample") or ""
                    )
            else:
                promoted.at[idx, "tumor_top_sample"] = str(tsup.get("top_non_tgct_sample") or "")
                if cls_now == "Reproductive Antigen only":
                    promoted.at[idx, "cancer_types_at_tenfold"] = 0
                    promoted.at[idx, "anonymous_samples_at_tenfold"] = 0
                    promoted.at[idx, "patient_or_sample_witnesses_at_tenfold"] = 0
        elif eid in obs_737_by_eid and pd.notna(t_tpm):
            old_t = float(old_row["tumor_max_tpm"]) if pd.notna(old_row["tumor_max_tpm"]) else -1.0
            if float(t_tpm) > old_t + 1e-6:
                best737 = max(obs_737_by_eid[eid], key=lambda x: float(x["value"]))
                promoted.at[idx, "tumor_top_sample"] = (
                    f"{best737['source_cohort']}:{best737['sample_id']}"
                )

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
        c_diff = {
            comp: orig_comps[comp] for comp in orig_comps if orig_comps[comp] != new_comps.get(comp)
        }
        if c_diff:
            he_deltas[eid] = c_diff

    delta_doc = {
        "schema": "oncoref.baseline-catalog-deltas.v1",
        "baseline_policy_id": "oncoref.target-selection.2026-11-20.2",
        "promoted_policy_id": PROMOTED_POLICY_ID,
        "baseline_package_dependencies_sha256": {
            "target-selection.csv.gz": _sha256_file(
                FROZEN_150_DIR / "oncoref/data/target-selection.csv.gz"
            ),
            "target-selection.json": _sha256_file(
                FROZEN_150_DIR / "oncoref/data/target-selection.json"
            ),
        },
        "added_columns": added_cols,
        "changed_entity_count": len(row_deltas),
        "baseline_host_exposure_entity_ids": sorted(base_he_df["entity_id"].astype(str).tolist()),
        "baseline_host_exposure_deltas": he_deltas,
        "rows": row_deltas,
    }
    delta_path = ROOT / "oncoref/data/baseline-2026-11-20-deltas.json"
    delta_path.write_text(json.dumps(delta_doc, indent=2, sort_keys=True) + "\n", encoding="utf-8")

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

    print(
        json.dumps(
            {
                "policy_id": PROMOTED_POLICY_ID,
                "categories": summ["categories"],
                "sex_restriction_codes": sp_meta["summary"]["sex_restriction_codes"],
                "delta_rows": len(row_deltas),
                "host_exposure_entities": len(he_rows),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

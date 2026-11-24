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
        c_diff = {
            comp: orig_comps[comp] for comp in orig_comps if orig_comps[comp] != new_comps.get(comp)
        }
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

# Licensed under the Apache License, Version 2.0.
"""Unified human SpeciesScreeningSpec adapter & v1.49 evidence integration (#755).

Adapts human Ensembl 115 GRCh38 protein-equivalence groups, normal tissue/protein/literature
gates, donor-resolved rare-cell RNA outcomes, #711/#712 normal adjudications, #735 unmeasured
home-tissue panels, and #737 linear-TPM tumor cohorts (THYMCA, BCC, cSCC) directly into
``species_target_selection.evaluate_species_catalog(human_spec, ...)``.

Default packaged catalog calls (``target_selection.results()``, ``human.target_selection()``)
and default ``evaluate_human_catalog()`` calls both evaluate the promoted integrated policy
``oncoref.target-selection.2026-11-25.1``, while explicit ``include_*=False`` flags reproduce
the frozen ``oncoref.target-selection.2026-11-20.2`` baseline and single-rule ablations.
"""

from __future__ import annotations

import hashlib
import json
from importlib.resources import files
from pathlib import Path

import numpy as np
import pandas as pd

from . import (
    assay_boundary_adjudications,
    clinical_organ_clearance,
    library_qc,
    normal_adjudications,
    normal_library_qc,
    rare_cell_compartment_qc,
    repeat_family_equivalence,
    species_target_selection,
    target_selection,
    tgct_localization,
    unmeasured_home_panels,
)
from .load_dataset import get_data
from .reference_identity import ReferenceIdentity
from .species import resolve_species
from .species_policies import species_policy
from .species_target_selection import (
    SpeciesScreeningResult,
    SpeciesScreeningSpec,
    evaluate_species_catalog,
    packaged_species_metadata,
)

HUMAN_BASELINE_COMPARISON_ID = "oncoref.native-gene-tpm-normal-reference.2026-11-20.1"
HUMAN_INTEGRATED_COMPARISON_ID = "oncoref.human-grch38-integrated-tpm.2026-11-25.1"
HUMAN_INTEGRATED_POLICY_PREVIEW_ID = "oncoref.target-selection.2026-11-25.1"

REQUIRED_PERIPHERAL_CONTEXT = "peripheral_normal_reference"


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def human_screening_spec(
    *,
    policy_id: str | None = None,
    include_home_panels: bool = True,
    include_737_cohorts: bool = True,
    include_multiplex_bleed_qc: bool = True,
    include_nuclear_spillover_qc: bool = True,
    include_clinical_clearances: bool = True,
    include_repeat_family_equivalence: bool = True,
    include_assay_boundary_fixes: bool = True,
    include_tumor_gdna_qc: bool | None = None,
) -> SpeciesScreeningSpec:
    """Construct the declarative ``SpeciesScreeningSpec`` for human Ensembl 115 GRCh38."""
    pol = species_policy("human")
    meta = packaged_species_metadata("human")
    ref = (
        ReferenceIdentity.from_dict(meta["reference"])
        if meta is not None
        else pol.primary_reference
    )
    base_any = (
        include_home_panels
        or include_737_cohorts
        or include_multiplex_bleed_qc
        or include_nuclear_spillover_qc
        or include_clinical_clearances
        or include_repeat_family_equivalence
        or include_assay_boundary_fixes
    )
    eff_gdna = base_any if include_tumor_gdna_qc is None else bool(include_tumor_gdna_qc)
    all_layers = (
        include_home_panels
        and include_737_cohorts
        and include_multiplex_bleed_qc
        and include_nuclear_spillover_qc
        and include_clinical_clearances
        and include_repeat_family_equivalence
        and include_assay_boundary_fixes
        and eff_gdna
    )
    any_layer = base_any or eff_gdna
    if policy_id is not None:
        effective_policy = policy_id
    elif all_layers:
        effective_policy = target_selection.POLICY_ID
    elif not any_layer:
        effective_policy = target_selection.BASELINE_2026_11_20_POLICY_ID
    else:
        effective_policy = "oncoref.human-integrated-ablation.2026-11-25.1"
    comp_id = HUMAN_INTEGRATED_COMPARISON_ID if any_layer else HUMAN_BASELINE_COMPARISON_ID
    notes = [
        "Evaluated through species_target_selection.evaluate_species_catalog (#755).",
    ]
    if include_home_panels:
        notes.append(
            "Includes #735 unit-separated adult normal panels for unmeasured home tissues."
        )
    if include_737_cohorts:
        notes.append("Includes #737 QC-passing linear-TPM tumor cohorts (THYMCA, BCC, cSCC).")
    if include_multiplex_bleed_qc:
        notes.append(
            "Includes #764 stoichiometric multiplex index-hopping / cross-sample bleed QC."
        )
    if include_nuclear_spillover_qc:
        notes.append(
            "Includes #765 nuclear snRNA-seq compartment and transcriptional spillover QC."
        )
    if include_clinical_clearances:
        notes.append(
            "Includes #766 epitope/isoform-scoped clinical organ-safety clearances & Xq28 attribution."
        )
    if include_repeat_family_equivalence:
        notes.append(
            "Includes #767 closed homologous repeat-family sequence equivalence & symmetric quantification."
        )
    if include_assay_boundary_fixes:
        notes.append(
            "Includes #768 TGCT admixture witness, RNA-zero HPA IHC, adult cornea multi-UMI, and germline TF rules."
        )
    if eff_gdna:
        notes.append("Includes tightened tumor gDNA library QC (olfactory_share_ge5 <= 2.5%).")
    return SpeciesScreeningSpec(
        species=resolve_species("human"),
        reference=ref,
        policy_id=effective_policy,
        comparison_id=comp_id,
        comparability_reference="docs/audits/target-selection-20261110/rebuild/manifest.json",
        comparison_unit="TPM",
        required_normal_contexts=frozenset({REQUIRED_PERIPHERAL_CONTEXT}),
        anatomy=pol.anatomy,
        rare_cell_provenance_resolver=target_selection.resolve_rare_cell_provenance,
        disposable_cap=pol.thresholds.disposable_cap,
        enrichment_fold=pol.thresholds.enrichment_fold,
        detection_floor=pol.thresholds.detection_floor,
        notes=tuple(notes),
    )


def upgraded_737_tumor_observations(*, cancer_code: str | None = None) -> pd.DataFrame:
    """Return SHA-256-verified QC-passing linear-TPM per-sample observations for #737 cohorts."""
    data_root = files("oncoref").joinpath("data")
    csv_res = data_root.joinpath("human-737-tumor-observations.csv.gz")
    meta_res = data_root.joinpath("human-737-tumor-observations.json")
    raw_bytes = csv_res.read_bytes()
    receipt = json.loads(meta_res.read_text(encoding="utf-8"))
    if _sha256_bytes(raw_bytes) != receipt["sha256"]:
        raise ValueError("human-737-tumor-observations.csv.gz does not match its SHA-256 receipt")
    df = get_data("human-737-tumor-observations")
    if cancer_code is not None:
        df = df[df["cancer_code"].eq(cancer_code)].reset_index(drop=True)
    return df


def build_human_screening_inputs(
    *,
    include_home_panels: bool = True,
    include_adjudications: bool = True,
    include_737_cohorts: bool = True,
    include_multiplex_bleed_qc: bool = True,
    include_nuclear_spillover_qc: bool = True,
    include_clinical_clearances: bool = True,
    include_repeat_family_equivalence: bool = True,
    include_assay_boundary_fixes: bool = True,
    include_tumor_gdna_qc: bool | None = None,
    extra_tumor_observations: pd.DataFrame | None = None,
    policy_id: str | None = None,
) -> dict:
    """Materialize human screening inputs in ``evaluate_species_catalog`` format."""
    base_any = (
        include_home_panels
        or include_737_cohorts
        or include_multiplex_bleed_qc
        or include_nuclear_spillover_qc
        or include_clinical_clearances
        or include_repeat_family_equivalence
        or include_assay_boundary_fixes
        or extra_tumor_observations is not None
    )
    eff_gdna = base_any if include_tumor_gdna_qc is None else bool(include_tumor_gdna_qc)
    spec = human_screening_spec(
        policy_id=policy_id,
        include_home_panels=include_home_panels,
        include_737_cohorts=include_737_cohorts or extra_tumor_observations is not None,
        include_multiplex_bleed_qc=include_multiplex_bleed_qc,
        include_nuclear_spillover_qc=include_nuclear_spillover_qc,
        include_clinical_clearances=include_clinical_clearances,
        include_repeat_family_equivalence=include_repeat_family_equivalence,
        include_assay_boundary_fixes=include_assay_boundary_fixes,
        include_tumor_gdna_qc=eff_gdna,
    )
    ts = target_selection._baseline_2026_11_20_frame()
    exposures = (
        species_target_selection.packaged_species_host_exposure("human")
        if include_multiplex_bleed_qc
        else target_selection._baseline_2026_11_20_host_exposures()
    )

    cleared_rare_eids: set[str] = set()
    cleared_obs_map: dict[str, dict] = {}
    admixture_non_tgct_map: dict[str, tuple[float, str, str]] = {}
    germline_tf_eids: set[str] = set()
    if include_nuclear_spillover_qc:
        nuc_cleared = rare_cell_compartment_qc.nuclear_spillover_adjudications(only_cleared=True)
        cleared_rare_eids = set(
            nuc_cleared.loc[
                nuc_cleared["disposition"].eq("nuclear_spillover_cleared"), "entity_id"
            ].astype(str)
        )
        cleared_obs_map = dict(
            rare_cell_compartment_qc.cleared_candidate_comparison_observations()["candidates"]
        )
    if include_clinical_clearances:
        coc_doc = clinical_organ_clearance.ledger()
        ctag2_attr = coc_doc.get("locus_attributions", {}).get("CTAG2")
        if ctag2_attr and any(
            e["entity_id"] == ctag2_attr["entity_id"]
            and e["status"] == "cleared"
            and "lung" in e.get("cleared_organs", [])
            for e in coc_doc["entries"]
        ):
            cleared_rare_eids.add(str(ctag2_attr["entity_id"]))
    if include_repeat_family_equivalence:
        rf_cleared = repeat_family_equivalence.repeat_family_adjudications(only_cleared=True)
        rf_obs = repeat_family_equivalence.repeat_family_comparison_observations()["candidates"]
        nuc_by_eid = (
            rare_cell_compartment_qc.nuclear_spillover_adjudications()
            .set_index("entity_id")
            .to_dict("index")
            if include_nuclear_spillover_qc
            else {}
        )
        for rf_row in rf_cleared.itertuples(index=False):
            rf_eid = str(rf_row.entity_id)
            rc_status = str(rf_row.baseline_rare_cell_status or "below_screen")
            can_clear = False
            if rc_status in {"below_screen", "dismissed"}:
                can_clear = True
            elif include_nuclear_spillover_qc and rf_eid in nuc_by_eid:
                nr = nuc_by_eid[rf_eid]
                re_eval = rare_cell_compartment_qc.evaluate_nuclear_rare_cell_disposition(
                    baseline_rare_cell_status=str(nr["baseline_rare_cell_status"]),
                    nucleus_only_signal=bool(nr["nucleus_only_signal"]),
                    finding_bulk_tpm=float(nr["finding_bulk_tpm"]),
                    somatic_ihc_positive=bool(nr["somatic_ihc_positive"]),
                    literature_effect=str(nr["literature_effect"]),
                    other_unsafe_gates=bool(nr["other_unsafe_gates"]),
                    transcriptional_spillover_risk=bool(nr["transcriptional_spillover_risk"]),
                    paralog_ambiguous=bool(nr["paralog_ambiguous"]),
                    host_overlapped=bool(nr["host_overlapped"]),
                    closed_repeat_family_cleared=True,
                )
                can_clear = bool(re_eval["cleared"])
            if can_clear:
                cleared_rare_eids.add(rf_eid)
                if rf_eid in rf_obs:
                    cleared_obs_map[rf_eid] = rf_obs[rf_eid]
    if include_assay_boundary_fixes:
        ab_df = assay_boundary_adjudications.assay_boundary_adjudications()
        ab_obs = assay_boundary_adjudications.assay_boundary_comparison_observations()["candidates"]
        for ab_row in ab_df.itertuples(index=False):
            ab_eid = str(ab_row.entity_id)
            disp = str(ab_row.disposition)
            if disp in {
                "ihc_rna_zero_cleared",
                "adult_cornea_sparse_umi_cleared",
                "sc_soup_and_cis_spillover_cleared",
            }:
                cleared_rare_eids.add(ab_eid)
                if ab_eid in ab_obs:
                    cleared_obs_map[ab_eid] = ab_obs[ab_eid]
            elif disp == "germline_tf_protamine_diluted_qualified":
                germline_tf_eids.add(ab_eid)
        tgct_restr = tgct_localization.tgct_admixture_restricted_entity_ids(
            include_multi_r_tgct=True
        )
        tgct_sup = tgct_localization.tgct_aafl_protein_support()
        for tr in tgct_sup.itertuples(index=False):
            teid = str(tr.entity_id)
            if teid in tgct_restr:
                admixture_non_tgct_map[teid] = (
                    float(tr.max_non_tgct_tpm) if pd.notna(tr.max_non_tgct_tpm) else 0.0,
                    str(tr.top_non_tgct_sample) if pd.notna(tr.top_non_tgct_sample) else "",
                    str(tr.top_non_tgct_cancer_code)
                    if pd.notna(tr.top_non_tgct_cancer_code)
                    else "",
                )

    readthrough = ts["readthrough_locus"].fillna(False).astype(bool)
    overlap_ge95 = (
        pd.to_numeric(ts["unsafe_sequence_overlap_pct"], errors="coerce").fillna(0.0) >= 95.0
    )
    seq_pass = ts["sequence_pass"].astype(bool)
    nested_unqualified = ~seq_pass & ~readthrough & ~overlap_ge95

    # Baseline unmeasured-home flag from the active catalog
    baseline_unmeasured_home = ts["classification"].eq("Other:Unscreened tissue") | (
        ts["is_reproductive_antigen"].astype(bool)
        & ~ts["is_cancer_antigen"].astype(bool)
        & ts["cancer_evidence"].eq("qualifies")
    )
    home_unmeasured = baseline_unmeasured_home.copy()

    gate_rows: list[dict] = []
    for r in ts.itertuples(index=False):
        if r.entity_id in cleared_rare_eids:
            continue
        if r.category == "Unsafe":
            gate_rows.append(
                {
                    "entity_id": r.entity_id,
                    "layer": "baseline_normal_safety",
                    "unsafe": True,
                    "hold": False,
                    "reason": str(r.Unsafe_reasons)
                    if pd.notna(r.Unsafe_reasons) and str(r.Unsafe_reasons)
                    else "Normal-tissue exclusion",
                }
            )
        elif not bool(r.normal_screen_pass):
            gate_rows.append(
                {
                    "entity_id": r.entity_id,
                    "layer": "baseline_normal_hold",
                    "unsafe": False,
                    "hold": True,
                    "reason": str(r.decision_reason)
                    if pd.notna(r.decision_reason) and str(r.decision_reason)
                    else "Normal-screen review hold",
                }
            )

    if include_adjudications:
        rc_adj = normal_adjudications.rare_cell_candidate_adjudications()
        for r in rc_adj.itertuples(index=False):
            if r.attribution_status == "confounder_rejected":
                gate_rows.append(
                    {
                        "entity_id": r.entity_id,
                        "layer": "rare_cell_attribution_711",
                        "unsafe": True,
                        "hold": False,
                        "reason": f"#711 confounder rejected ({r.finding_organ} {r.finding_top_cell_type})",
                    }
                )
            else:
                gate_rows.append(
                    {
                        "entity_id": r.entity_id,
                        "layer": "rare_cell_attribution_711",
                        "unsafe": False,
                        "hold": True,
                        "reason": f"#711 pending {r.resolution_requirement}",
                    }
                )

    home_bulk_observations: dict[str, tuple[str, str, float]] = {}
    if include_home_panels:
        uh = unmeasured_home_panels.unmeasured_home_protein_adjudications(only_ca_cta_gated=True)
        eid_to_idx = {eid: idx for idx, eid in enumerate(ts["entity_id"])}
        for r in uh.itertuples(index=False):
            idx = eid_to_idx.get(r.entity_id)
            if idx is None:
                continue
            disp = str(r.home_panel_disposition)
            if disp == "home_expressed_above_tier_ceiling":
                cross = (
                    f"; {r.cross_panel_ceiling_violation}"
                    if pd.notna(r.cross_panel_ceiling_violation)
                    and str(r.cross_panel_ceiling_violation).strip()
                    else ""
                )
                gate_rows.append(
                    {
                        "entity_id": r.entity_id,
                        "layer": "home_tissue_panel_735",
                        "unsafe": True,
                        "hold": False,
                        "reason": (
                            f"#735 adult {r.home_tissue} ({r.home_tissue_tier}) exceeds ceiling"
                            f"{cross}"
                        ),
                    }
                )
                home_unmeasured.iat[idx] = False
            elif disp in {
                "home_cleared_panel_passes_tenfold",
                "home_disposable_follicle_expressed_tenfold",
                "home_expressed_insufficient_tumor_fold",
            }:
                home_unmeasured.iat[idx] = False
                if bool(r.home_bulk_measured) and pd.notna(r.home_bulk_p95_tpm):
                    scope = "disposable" if str(r.home_tissue_tier) == "disposable" else "normal"
                    home_bulk_observations[r.entity_id] = (
                        str(r.home_tissue),
                        scope,
                        float(r.home_bulk_p95_tpm),
                    )
            elif disp == "home_undetected_in_sampled_panel":
                # Conservative fail-closed rule: incomplete anatomical sub-niches stay unscreened
                home_unmeasured.iat[idx] = True

    proteins = pd.DataFrame(
        {
            "entity_id": ts["entity_id"],
            "protein_group_label": ts["protein_group_label"],
            "member_gene_ids": ts["member_gene_ids"],
            "representative_gene_id": ts["representative_gene_id"],
            "member_count": ts["member_gene_ids"].map(lambda s: len(json.loads(s))),
            "representative_sequence_sha256s": ts["representative_sequence_sha256s"].fillna("[]"),
            "identity_translation_pass": ts["identity_translation_pass"].astype(bool),
            "required_coverage_complete": ts["required_coverage_complete"].astype(bool),
            "unresolved_evidence": ~ts["evidence_resolved"].astype(bool),
            "readthrough_locus": readthrough,
            "nested_unqualified": nested_unqualified,
            "unsafe_sequence_overlap_pct": pd.to_numeric(
                ts["unsafe_sequence_overlap_pct"], errors="coerce"
            ),
            "disposable_max_ntpm": pd.to_numeric(ts["disposable_max_ntpm"], errors="coerce"),
            "developmental_mean_tpm": pd.to_numeric(ts["developmental_mean_tpm"], errors="coerce"),
            "home_tissue": ts["home_tissue"].fillna("").astype(str),
            "home_tissue_unmeasured": home_unmeasured,
            "published_nomination": ts["published_nomination"].fillna(False).astype(bool),
        }
    )
    if cleared_obs_map:
        eid_to_pidx = {eid: idx for idx, eid in enumerate(proteins["entity_id"])}
        for ceid, cpayload in cleared_obs_map.items():
            pidx = eid_to_pidx.get(ceid)
            if pidx is None:
                continue
            proteins.at[pidx, "required_coverage_complete"] = bool(
                cpayload["required_coverage_complete"]
            )
            proteins.at[pidx, "unsafe_sequence_overlap_pct"] = float(
                cpayload["unsafe_sequence_overlap_pct"]
            )
            proteins.at[pidx, "disposable_max_ntpm"] = float(cpayload["disposable_max_ntpm"])
            proteins.at[pidx, "nested_unqualified"] = False

    bleed_clean_by_eid: dict[str, dict] = {}
    if include_multiplex_bleed_qc:
        bleed_df = normal_library_qc.multiplex_bleed_gene_deltas(only_applied=True)
        for br in bleed_df.to_dict("records"):
            bleed_clean_by_eid[str(br["entity_id"])] = br

    gdna_clean_by_eid: dict[str, dict] = {}
    if eff_gdna:
        gdna_df = library_qc.tumor_gdna_background_deltas(only_applied=True)
        for gr in gdna_df.to_dict("records"):
            gdna_clean_by_eid[str(gr["entity_id"])] = gr

    comp_id = spec.comparison_id
    unit = spec.comparison_unit
    normal_rows: list[dict] = []
    tumor_rows: list[dict] = []
    repro_rows: list[dict] = []

    for r in ts.itertuples(index=False):
        eid = r.entity_id
        if eid in cleared_obs_map:
            cpayload = cleared_obs_map[eid]
            if not bool(cpayload["required_coverage_complete"]):
                continue
            if cpayload.get("normal_reference_tpm") is not None:
                normal_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": "cleared_candidate_normal_reference_765",
                        "context_id": REQUIRED_PERIPHERAL_CONTEXT,
                        "tissue": REQUIRED_PERIPHERAL_CONTEXT,
                        "scope": "normal",
                        "value": float(cpayload["normal_reference_tpm"]),
                        "comparison_id": comp_id,
                        "unit": unit,
                        "qc_pass": True,
                        "identity_resolved": True,
                    }
                )
            if cpayload.get("disposable_max_tpm") is not None:
                normal_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": "cleared_candidate_normal_reference_765",
                        "context_id": "disposable_reference",
                        "tissue": "disposable_reference",
                        "scope": "disposable",
                        "value": float(cpayload["disposable_max_tpm"]),
                        "comparison_id": comp_id,
                        "unit": unit,
                        "qc_pass": True,
                        "identity_resolved": True,
                    }
                )
            if eid in gdna_clean_by_eid:
                grec = gdna_clean_by_eid[eid]
                c_tmax = float(grec["clean_tumor_max_tpm"])
                if c_tmax > 0.0 and str(grec.get("clean_tumor_top_sample") or "").strip():
                    c_tsmp = str(grec["clean_tumor_top_sample"]).strip()
                    tumor_rows.append(
                        {
                            "entity_id": eid,
                            "source_id": "oncoref_tumor_gdna_qc_20261125",
                            "sample_id": c_tsmp,
                            "participant_id": (
                                f"clean:{c_tsmp}"
                                if int(grec.get("clean_independent_patients_at_tenfold") or 0) > 0
                                else ""
                            ),
                            "cancer_code": str(grec.get("clean_tumor_top_cancer_code") or ""),
                            "value": c_tmax,
                            "comparison_id": comp_id,
                            "unit": unit,
                            "qc_pass": True,
                            "identity_resolved": True,
                        }
                    )
            else:
                for tobs in cpayload.get("top_tumor_observations", []):
                    tumor_rows.append(
                        {
                            "entity_id": eid,
                            "source_id": str(tobs["source_id"]),
                            "sample_id": str(tobs["sample_id"]),
                            "participant_id": str(tobs.get("participant_id") or ""),
                            "cancer_code": str(tobs.get("cancer_code") or ""),
                            "value": float(tobs["value"]),
                            "comparison_id": comp_id,
                            "unit": unit,
                            "qc_pass": True,
                            "identity_resolved": True,
                        }
                    )
            supp_val = cpayload.get("reproductive_supported_tpm")
            he_eid = exposures.get(eid, {})
            for citem in cpayload.get("reproductive_compartments", []):
                comp = str(citem["compartment"])
                he_comp = he_eid.get(comp, {})
                is_emb = citem.get("scope") == "developmental" or comp == "embryo"
                cqual = bool(citem.get("qualifies_tenfold", False))
                raw_cval = float(citem["max_tpm"])
                he_cval = float(
                    he_comp["max_tpm"] if he_comp.get("max_tpm") is not None else raw_cval
                )
                cval = max(raw_cval, he_cval) if cqual else he_cval
                c_ntpm = he_comp.get("max_ntpm") if "max_ntpm" in he_comp else citem.get("max_ntpm")
                c_prot = (
                    he_comp.get("protein_status")
                    if "protein_status" in he_comp
                    else citem.get("protein_status", "not_assessed")
                )
                if is_emb:
                    repro_rows.append(
                        {
                            "entity_id": eid,
                            "source_id": "human_embryo_atlases",
                            "sample_id": f"repro:{comp}",
                            "participant_id": "embryo_d1" if cqual else "",
                            "scope": "developmental",
                            "compartment": comp,
                            "source_context": "embryo_lineage",
                            "developmental_stage": "preimplantation",
                            "cells": 25 if cqual else 20,
                            "positive_cells": 10 if cqual else 1,
                            "value": cval,
                            "host_max_tpm": he_cval,
                            "compartment_qualifying_value": (
                                float(supp_val) if cqual and supp_val is not None else 0.0
                            ),
                            "max_ntpm": c_ntpm,
                            "protein_status": c_prot,
                            "scope_complete": False,
                            "comparison_id": comp_id,
                            "unit": unit,
                            "qc_pass": True,
                            "identity_resolved": True,
                        }
                    )
                else:
                    repro_rows.append(
                        {
                            "entity_id": eid,
                            "source_id": "human_reproductive_panel",
                            "sample_id": f"repro:{comp}",
                            "participant_id": f"repro_donor:{comp}",
                            "scope": "reproductive",
                            "compartment": comp,
                            "source_context": comp,
                            "developmental_stage": "adult",
                            "cells": 0,
                            "positive_cells": 0,
                            "value": cval,
                            "host_max_tpm": he_cval,
                            "compartment_qualifying_value": cval if cqual else min(cval, 0.0),
                            "female_max_tpm": he_comp.get("female_max_tpm"),
                            "male_max_tpm": he_comp.get("male_max_tpm"),
                            "female_max_ntpm": he_comp.get("female_max_ntpm"),
                            "male_max_ntpm": he_comp.get("male_max_ntpm"),
                            "max_ntpm": c_ntpm,
                            "protein_status": c_prot,
                            "scope_complete": True,
                            "comparison_id": comp_id,
                            "unit": unit,
                            "qc_pass": True,
                            "identity_resolved": True,
                        }
                    )
            continue
        if r.category == "Unsafe" or not bool(r.required_coverage_complete):
            continue
        bleed_rec = bleed_clean_by_eid.get(eid)
        norm_val = (
            float(bleed_rec["clean_normal_reference_tpm"])
            if bleed_rec is not None
            else (float(r.normal_reference_tpm) if pd.notna(r.normal_reference_tpm) else np.nan)
        )
        if pd.notna(norm_val):
            normal_rows.append(
                {
                    "entity_id": eid,
                    "source_id": (
                        "human_normal_reference_bleed_qc_764"
                        if bleed_rec is not None
                        else "human_normal_reference_20261120"
                    ),
                    "context_id": REQUIRED_PERIPHERAL_CONTEXT,
                    "tissue": REQUIRED_PERIPHERAL_CONTEXT,
                    "scope": "normal",
                    "value": norm_val,
                    "comparison_id": comp_id,
                    "unit": unit,
                    "qc_pass": True,
                    "identity_resolved": True,
                }
            )
        if pd.notna(r.disposable_max_tpm):
            normal_rows.append(
                {
                    "entity_id": eid,
                    "source_id": "human_normal_reference_20261120",
                    "context_id": "disposable_reference",
                    "tissue": "disposable_reference",
                    "scope": "disposable",
                    "value": float(r.disposable_max_tpm),
                    "comparison_id": comp_id,
                    "unit": unit,
                    "qc_pass": True,
                    "identity_resolved": True,
                }
            )
        if eid in home_bulk_observations:
            htissue, hscope, hval = home_bulk_observations[eid]
            normal_rows.append(
                {
                    "entity_id": eid,
                    "source_id": "unmeasured_home_panels_735",
                    "context_id": f"home_735:{htissue}",
                    "tissue": htissue,
                    "scope": hscope,
                    "value": hval,
                    "comparison_id": comp_id,
                    "unit": unit,
                    "qc_pass": True,
                    "identity_resolved": True,
                }
            )

        if eid in gdna_clean_by_eid and eid in admixture_non_tgct_map:
            grec = gdna_clean_by_eid[eid]
            c_tmax = float(grec["clean_non_tgct_max_tpm"])
            c_tsmp = str(grec.get("clean_non_tgct_top_sample") or "").strip()
            if c_tmax > 0.0 and c_tsmp:
                tumor_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": "oncoref_non_admixed_tumor_witness_768",
                        "sample_id": c_tsmp,
                        "participant_id": (
                            f"somatic:{c_tsmp}"
                            if int(grec.get("clean_non_tgct_independent_patients_at_tenfold") or 0)
                            > 0
                            else ""
                        ),
                        "cancer_code": str(grec.get("clean_non_tgct_top_cancer_code") or ""),
                        "value": c_tmax,
                        "comparison_id": comp_id,
                        "unit": unit,
                        "qc_pass": True,
                        "identity_resolved": True,
                    }
                )
        elif eid in gdna_clean_by_eid:
            grec = gdna_clean_by_eid[eid]
            c_tmax = float(grec["clean_tumor_max_tpm"])
            c_tsmp = str(grec.get("clean_tumor_top_sample") or "").strip()
            if c_tmax > 0.0 and c_tsmp:
                tumor_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": "oncoref_tumor_gdna_qc_20261125",
                        "sample_id": c_tsmp,
                        "participant_id": (
                            f"clean:{c_tsmp}"
                            if int(grec.get("clean_independent_patients_at_tenfold") or 0) > 0
                            else ""
                        ),
                        "cancer_code": str(grec.get("clean_tumor_top_cancer_code") or ""),
                        "value": c_tmax,
                        "comparison_id": comp_id,
                        "unit": unit,
                        "qc_pass": True,
                        "identity_resolved": True,
                    }
                )
        elif eid in admixture_non_tgct_map:
            nt_val, nt_sample, nt_code = admixture_non_tgct_map[eid]
            if nt_val > 0 and nt_sample:
                tumor_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": "oncoref_non_admixed_tumor_witness_768",
                        "sample_id": nt_sample,
                        "participant_id": f"somatic:{nt_sample}",
                        "cancer_code": nt_code,
                        "value": nt_val,
                        "comparison_id": comp_id,
                        "unit": unit,
                        "qc_pass": True,
                        "identity_resolved": True,
                    }
                )
        elif pd.notna(r.tumor_max_tpm):
            n_known = (
                int(r.independent_patients_at_tenfold or 0)
                if pd.notna(r.independent_patients_at_tenfold)
                else 0
            )
            top_sample = (
                str(r.tumor_top_sample)
                if pd.notna(r.tumor_top_sample) and str(r.tumor_top_sample)
                else "baseline_peak_sample"
            )
            top_code = (
                str(r.tumor_top_cancer_code)
                if pd.notna(r.tumor_top_cancer_code) and str(r.tumor_top_cancer_code)
                else ""
            )
            tumor_rows.append(
                {
                    "entity_id": eid,
                    "source_id": "oncoref_tumor_union_20261120",
                    "sample_id": top_sample,
                    "participant_id": f"baseline:{top_sample}" if n_known > 0 else "",
                    "cancer_code": top_code,
                    "value": float(r.tumor_max_tpm),
                    "comparison_id": comp_id,
                    "unit": unit,
                    "qc_pass": True,
                    "identity_resolved": True,
                }
            )
        if (
            bleed_rec is not None
            and eid not in gdna_clean_by_eid
            and eid not in admixture_non_tgct_map
            and pd.notna(bleed_rec.get("max_non_tgct_tumor_tpm"))
            and float(bleed_rec["max_non_tgct_tumor_tpm"]) > 0
            and pd.notna(bleed_rec.get("top_non_tgct_sample"))
            and str(bleed_rec.get("top_non_tgct_sample")).strip()
        ):
            nt_sample = str(bleed_rec["top_non_tgct_sample"]).strip()
            if nt_sample != (str(r.tumor_top_sample) if pd.notna(r.tumor_top_sample) else ""):
                nt_cc_raw = bleed_rec.get("top_non_tgct_cancer_code")
                tumor_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": "oncoref_somatic_tumor_witness_764",
                        "sample_id": nt_sample,
                        "participant_id": f"somatic:{nt_sample}",
                        "cancer_code": str(nt_cc_raw).strip() if pd.notna(nt_cc_raw) else "",
                        "value": float(bleed_rec["max_non_tgct_tumor_tpm"]),
                        "comparison_id": comp_id,
                        "unit": unit,
                        "qc_pass": True,
                        "identity_resolved": True,
                    }
                )

        q_raw_val = r.reproductive_qualifying_compartments
        qual_comps = (
            {c.strip() for c in str(q_raw_val).split("|") if c.strip()}
            if pd.notna(q_raw_val) and str(q_raw_val).strip()
            else set()
        )
        rec_exp = exposures.get(eid, {})
        supp_tpm = (
            float(r.reproductive_supported_tpm)
            if pd.notna(r.reproductive_supported_tpm)
            else np.nan
        )
        added_comps: set[str] = set()
        for comp, item in sorted(rec_exp.items()):
            if not isinstance(item, dict) or item.get("coverage") != "measured":
                continue
            val = item.get("max_tpm")
            if val is None or not np.isfinite(val):
                continue
            added_comps.add(comp)
            is_embryo = comp == "embryo"
            comp_qualifies = comp in qual_comps
            if is_embryo:
                repro_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": "human_embryo_atlases",
                        "sample_id": f"repro:{comp}",
                        "participant_id": "embryo_d1" if comp_qualifies else "",
                        "scope": "developmental",
                        "compartment": comp,
                        "source_context": "embryo_lineage",
                        "developmental_stage": "preimplantation",
                        "cells": 25 if comp_qualifies else 20,
                        "positive_cells": 10 if comp_qualifies else 1,
                        "value": float(val),
                        "compartment_qualifying_value": (
                            supp_tpm if comp_qualifies and pd.notna(supp_tpm) else 0.0
                        ),
                        "max_ntpm": item.get("max_ntpm"),
                        "protein_status": item.get("protein_status", "not_assessed"),
                        "scope_complete": bool(item.get("scope_complete", False)),
                        "comparison_id": comp_id,
                        "unit": unit,
                        "qc_pass": True,
                        "identity_resolved": True,
                    }
                )
            else:
                raw_val = float(val)
                peak_val = (
                    max(raw_val, supp_tpm) if (comp_qualifies and pd.notna(supp_tpm)) else raw_val
                )
                if (
                    eid in germline_tf_eids
                    and comp == "testis"
                    and pd.notna(r.reproductive_max_sample_tpm)
                ):
                    peak_val = max(peak_val, float(r.reproductive_max_sample_tpm))
                repro_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": "human_reproductive_panel",
                        "sample_id": f"repro:{comp}",
                        "participant_id": f"repro_donor:{comp}",
                        "scope": "reproductive",
                        "compartment": comp,
                        "source_context": comp,
                        "developmental_stage": "adult",
                        "cells": 0,
                        "positive_cells": 0,
                        "value": peak_val,
                        "host_max_tpm": raw_val,
                        "compartment_qualifying_value": (
                            peak_val if comp_qualifies else min(raw_val, 0.0)
                        ),
                        "germline_protamine_diluted_qualified": bool(
                            eid in germline_tf_eids and comp == "testis"
                        ),
                        "female_max_tpm": item.get("female_max_tpm"),
                        "male_max_tpm": item.get("male_max_tpm"),
                        "female_max_ntpm": item.get("female_max_ntpm"),
                        "male_max_ntpm": item.get("male_max_ntpm"),
                        "max_ntpm": item.get("max_ntpm"),
                        "protein_status": item.get("protein_status", "not_assessed"),
                        "scope_complete": bool(item.get("scope_complete", True)),
                        "comparison_id": comp_id,
                        "unit": unit,
                        "qc_pass": True,
                        "identity_resolved": True,
                    }
                )
        if eid in germline_tf_eids and "testis" not in added_comps:
            g_max = (
                float(r.reproductive_max_sample_tpm)
                if pd.notna(r.reproductive_max_sample_tpm)
                else 0.0
            )
            if g_max >= 5.0:
                repro_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": "human_reproductive_panel_768",
                        "sample_id": "repro:testis",
                        "participant_id": "repro_donor:testis",
                        "scope": "reproductive",
                        "compartment": "testis",
                        "source_context": "testis",
                        "developmental_stage": "adult",
                        "cells": 0,
                        "positive_cells": 0,
                        "value": g_max,
                        "compartment_qualifying_value": g_max,
                        "germline_protamine_diluted_qualified": True,
                        "max_ntpm": None,
                        "protein_status": "not_assessed",
                        "scope_complete": True,
                        "comparison_id": comp_id,
                        "unit": unit,
                        "qc_pass": True,
                        "identity_resolved": True,
                    }
                )
        # Ensure any qualifying compartment lacking a host_exposure entry is still represented
        for comp in sorted(qual_comps - added_comps):
            if pd.notna(supp_tpm):
                repro_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": "human_reproductive_panel",
                        "sample_id": f"repro_qual:{comp}",
                        "participant_id": f"repro_donor:{comp}",
                        "scope": "reproductive",
                        "compartment": comp,
                        "source_context": comp,
                        "developmental_stage": "adult",
                        "cells": 0,
                        "positive_cells": 0,
                        "value": supp_tpm,
                        "compartment_qualifying_value": supp_tpm,
                        "max_ntpm": None,
                        "protein_status": "not_assessed",
                        "scope_complete": comp != "embryo",
                        "comparison_id": comp_id,
                        "unit": unit,
                        "qc_pass": True,
                        "identity_resolved": True,
                    }
                )

    if include_737_cohorts:
        obs_737 = upgraded_737_tumor_observations()
        valid_entities = set(
            ts.loc[
                ~ts["category"].eq("Unsafe") & ts["required_coverage_complete"].astype(bool),
                "entity_id",
            ]
        ) | {
            ceid
            for ceid, cpayload in cleared_obs_map.items()
            if bool(cpayload.get("required_coverage_complete"))
        }
        sub_737 = obs_737[obs_737["entity_id"].isin(valid_entities)]
        for r in sub_737.itertuples(index=False):
            tumor_rows.append(
                {
                    "entity_id": str(r.entity_id),
                    "source_id": str(r.source_cohort),
                    "sample_id": str(r.sample_id),
                    "participant_id": str(r.participant_id) if pd.notna(r.participant_id) else "",
                    "cancer_code": str(r.cancer_code),
                    "value": float(r.value),
                    "comparison_id": comp_id,
                    "unit": unit,
                    "qc_pass": True,
                    "identity_resolved": True,
                }
            )

    if extra_tumor_observations is not None and not extra_tumor_observations.empty:
        extra = extra_tumor_observations.copy()
        extra["comparison_id"] = comp_id
        extra["unit"] = unit
        tumor_rows.extend(extra.to_dict("records"))

    return {
        "spec": spec,
        "proteins": proteins,
        "normal_gates": pd.DataFrame(gate_rows),
        "comparison_normals": pd.DataFrame(normal_rows),
        "comparison_tumors": pd.DataFrame(tumor_rows),
        "comparison_reproductive": pd.DataFrame(repro_rows),
    }


def evaluate_human_catalog(
    *,
    include_home_panels: bool = True,
    include_adjudications: bool = True,
    include_737_cohorts: bool = True,
    include_multiplex_bleed_qc: bool = True,
    include_nuclear_spillover_qc: bool = True,
    include_clinical_clearances: bool = True,
    include_repeat_family_equivalence: bool = True,
    include_assay_boundary_fixes: bool = True,
    include_tumor_gdna_qc: bool | None = None,
    extra_tumor_observations: pd.DataFrame | None = None,
    policy_id: str | None = None,
) -> SpeciesScreeningResult:
    """Evaluate the 19,893-protein human catalog through ``evaluate_species_catalog``."""
    inputs = build_human_screening_inputs(
        include_home_panels=include_home_panels,
        include_adjudications=include_adjudications,
        include_737_cohorts=include_737_cohorts,
        include_multiplex_bleed_qc=include_multiplex_bleed_qc,
        include_nuclear_spillover_qc=include_nuclear_spillover_qc,
        include_clinical_clearances=include_clinical_clearances,
        include_repeat_family_equivalence=include_repeat_family_equivalence,
        include_assay_boundary_fixes=include_assay_boundary_fixes,
        include_tumor_gdna_qc=include_tumor_gdna_qc,
        extra_tumor_observations=extra_tumor_observations,
        policy_id=policy_id,
    )
    spec = inputs.pop("spec")
    result = evaluate_species_catalog(spec, **inputs)

    base_any_eval = (
        include_home_panels
        or include_737_cohorts
        or include_multiplex_bleed_qc
        or include_nuclear_spillover_qc
        or include_clinical_clearances
        or include_repeat_family_equivalence
        or include_assay_boundary_fixes
        or extra_tumor_observations is not None
    )
    eff_gdna_eval = base_any_eval if include_tumor_gdna_qc is None else bool(include_tumor_gdna_qc)
    is_integrated = base_any_eval or eff_gdna_eval
    # Attach #711, #712, #735, #764-#769 adjudication & near-self provenance columns
    cat = clinical_organ_clearance.annotate_catalog(result.catalog, enabled=True)
    cat = repeat_family_equivalence.annotate_catalog(cat, enabled=True)
    cat = assay_boundary_adjudications.annotate_catalog(cat, enabled=True, integrated=is_integrated)
    object.__setattr__(result, "catalog", cat)
    rc_map = (
        normal_adjudications.rare_cell_candidate_adjudications()
        .set_index("entity_id")["resolution_requirement"]
        .to_dict()
    )
    nl_map = (
        normal_adjudications.normal_library_peak_adjudications()
        .drop_duplicates("entity_id")
        .set_index("entity_id")["adjudication_disposition"]
        .to_dict()
    )
    uh_map = (
        unmeasured_home_panels.unmeasured_home_protein_adjudications()
        .set_index("entity_id")["home_panel_disposition"]
        .to_dict()
    )
    bleed_set = (
        set(normal_library_qc.multiplex_bleed_gene_deltas(only_applied=True)["entity_id"])
        if include_multiplex_bleed_qc
        else set()
    )
    rf_disp_map: dict[str, str] = {}
    rf_cleared_set: set[str] = set()
    if include_repeat_family_equivalence:
        rf_df = repeat_family_equivalence.repeat_family_adjudications()
        for rr in rf_df.itertuples(index=False):
            eid_str = str(rr.entity_id)
            rf_disp_map[eid_str] = str(rr.disposition)
            if bool(rr.closed_repeat_family_cleared):
                rf_cleared_set.add(eid_str)

    nuc_disp_map: dict[str, str] = {}
    nuc_caution_map: dict[str, str] = {}
    nuc_cleared_rare: set[str] = set()
    if include_nuclear_spillover_qc:
        nuc_df = rare_cell_compartment_qc.nuclear_spillover_adjudications()
        for nr in nuc_df.itertuples(index=False):
            eid_str = str(nr.entity_id)
            disp_str = str(nr.disposition)
            caut_str = (
                str(nr.nuclear_rna_caution)
                if pd.notna(nr.nuclear_rna_caution) and str(nr.nuclear_rna_caution)
                else ""
            )
            is_clr = bool(nr.nuclear_spillover_cleared) and disp_str == "nuclear_spillover_cleared"
            if not is_clr and eid_str in rf_cleared_set:
                re_eval = rare_cell_compartment_qc.evaluate_nuclear_rare_cell_disposition(
                    baseline_rare_cell_status=str(nr.baseline_rare_cell_status),
                    nucleus_only_signal=bool(nr.nucleus_only_signal),
                    finding_bulk_tpm=float(nr.finding_bulk_tpm),
                    somatic_ihc_positive=bool(nr.somatic_ihc_positive),
                    literature_effect=str(nr.literature_effect),
                    other_unsafe_gates=bool(nr.other_unsafe_gates),
                    transcriptional_spillover_risk=bool(nr.transcriptional_spillover_risk),
                    paralog_ambiguous=bool(nr.paralog_ambiguous),
                    host_overlapped=bool(nr.host_overlapped),
                    closed_repeat_family_cleared=True,
                )
                if bool(re_eval["cleared"]):
                    disp_str = str(re_eval["disposition"])
                    caut_str = str(re_eval["nuclear_rna_caution"])
                    is_clr = True
            nuc_disp_map[eid_str] = disp_str
            if caut_str:
                nuc_caution_map[eid_str] = caut_str
            if is_clr:
                nuc_cleared_rare.add(eid_str)

    clin_disp_map: dict[str, str] = {}
    if include_clinical_clearances:
        coc_doc = clinical_organ_clearance.ledger()
        ctag2_attr = coc_doc.get("locus_attributions", {}).get("CTAG2")
        if ctag2_attr and any(
            e["entity_id"] == ctag2_attr["entity_id"]
            and e["status"] == "cleared"
            and "lung" in e.get("cleared_organs", [])
            for e in coc_doc["entries"]
        ):
            clin_disp_map[str(ctag2_attr["entity_id"])] = "clinical_organ_and_xq28_cleared"

    ab_disp_map: dict[str, str] = {}
    if include_assay_boundary_fixes:
        ab_df = assay_boundary_adjudications.assay_boundary_adjudications()
        tgct_restr_eval = tgct_localization.tgct_admixture_restricted_entity_ids(
            include_multi_r_tgct=True
        )
        cat_by_eid = cat.set_index("entity_id")
        for ar in ab_df.itertuples(index=False):
            aeid = str(ar.entity_id)
            if str(ar.rule_family) == "tgct_admixture_witness_713" and aeid in cat_by_eid.index:
                crow = cat_by_eid.loc[aeid]
                n_v = (
                    float(crow["normal_reference_tpm"])
                    if pd.notna(crow["normal_reference_tpm"])
                    else 0.0
                )
                d_v = (
                    float(crow["disposable_max_tpm"])
                    if pd.notna(crow["disposable_max_tpm"])
                    else 0.0
                )
                t_v = float(crow["tumor_max_tpm"]) if pd.notna(crow["tumor_max_tpm"]) else 0.0
                dyn_ev = assay_boundary_adjudications.evaluate_non_admixed_tumor_witness(
                    entity_id=aeid,
                    max_non_tgct_tpm=t_v,
                    single_tumor_threshold_tpm=10.0 * max(n_v, d_v, 1.0),
                    admixture_restricted_ids=tgct_restr_eval,
                )
                ab_disp_map[aeid] = str(dyn_ev["disposition"])
            else:
                ab_disp_map[aeid] = str(ar.disposition)

    cat["rare_cell_711_requirement"] = [rc_map.get(eid, "") for eid in cat["entity_id"]]
    cat["normal_library_712_disposition"] = [nl_map.get(eid, "") for eid in cat["entity_id"]]
    cat["home_panel_735_disposition"] = [uh_map.get(eid, "") for eid in cat["entity_id"]]
    cat["multiplex_bleed_764_disposition"] = [
        "multiplex_bleed_cleaned" if eid in bleed_set else "" for eid in cat["entity_id"]
    ]
    cat["nuclear_spillover_765_disposition"] = [
        nuc_disp_map.get(eid, "") for eid in cat["entity_id"]
    ]
    cat["nuclear_rna_caution"] = [nuc_caution_map.get(eid, "") for eid in cat["entity_id"]]
    cat["clinical_clearance_766_disposition"] = [
        clin_disp_map.get(eid, "") for eid in cat["entity_id"]
    ]
    cat["repeat_family_767_disposition"] = [rf_disp_map.get(eid, "") for eid in cat["entity_id"]]
    cat["assay_boundary_768_disposition"] = [ab_disp_map.get(eid, "") for eid in cat["entity_id"]]
    if nuc_cleared_rare:
        mask_nuc = cat["entity_id"].isin(nuc_cleared_rare)
        cat.loc[mask_nuc, "rare_cell_status"] = "nuclear_only"
        cat.loc[mask_nuc, "rare_cell_action"] = "caution"
    return result


_PREVIEW_STRING_COLS = (
    "entity_id",
    "protein_group_label",
    "category_2026_11_20",
    "classification_2026_11_20",
    "home_tissue",
    "tumor_top_cancer_code_2026_11_20",
    "reproductive_qualifying_compartments",
    "category_integrated",
    "classification_integrated",
    "tumor_top_cancer_code_integrated",
    "sex_restriction_code",
    "home_panel_735_disposition",
    "rare_cell_711_requirement",
    "normal_library_712_disposition",
    "multiplex_bleed_764_disposition",
    "nuclear_spillover_765_disposition",
    "nuclear_rna_caution",
    "clinical_clearance_766_disposition",
    "repeat_family_767_disposition",
    "assay_boundary_768_disposition",
    "closed_homolog_family_id",
    "co_silenced_coding_sibling_ids",
    "closed_self_paralog_ids",
    "clinical_clearance_scope",
    "clinical_cleared_organs",
    "clinical_cleared_epitopes",
    "clinical_hla_alleles",
    "clinical_cleared_transcripts",
    "clinical_trial_source_ids",
    "transition_driver",
    "transition",
)


def compute_next_policy_transitions(
    *,
    integrated_result: SpeciesScreeningResult | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Recompute exact protein-level transitions from folding #735, #711/#712, #737, and #764-#768 into human screening."""
    baseline = target_selection._baseline_2026_11_20_frame()
    integrated_res = (
        integrated_result
        if integrated_result is not None
        else evaluate_human_catalog(
            include_home_panels=True,
            include_adjudications=True,
            include_737_cohorts=True,
            include_multiplex_bleed_qc=True,
            include_nuclear_spillover_qc=True,
            include_clinical_clearances=True,
            include_repeat_family_equivalence=True,
            include_assay_boundary_fixes=True,
        )
    )
    integrated = integrated_res.catalog

    merged = baseline[
        [
            "entity_id",
            "protein_group_label",
            "category",
            "classification",
            "home_tissue",
            "normal_reference_tpm",
            "tumor_max_tpm",
            "tumor_top_cancer_code",
            "reproductive_qualifying_compartments",
        ]
    ].merge(
        integrated[
            [
                "entity_id",
                "category",
                "classification",
                "normal_reference_tpm",
                "tumor_max_tpm",
                "tumor_top_cancer_code",
                "sex_restriction_code",
                "home_panel_735_disposition",
                "rare_cell_711_requirement",
                "normal_library_712_disposition",
                "multiplex_bleed_764_disposition",
                "nuclear_spillover_765_disposition",
                "nuclear_rna_caution",
                "clinical_clearance_766_disposition",
                "repeat_family_767_disposition",
                "assay_boundary_768_disposition",
                "closed_homolog_family_id",
                "co_silenced_coding_sibling_ids",
                "closed_self_paralog_ids",
                "exact_unique_9mers",
                "safe_2mm_9mers",
                "unsafe_near_self_1mm_pct",
                "unsafe_shared_tryptic_peptides",
                "clinical_clearance_scope",
                "clinical_cleared_organs",
                "clinical_cleared_epitopes",
                "clinical_hla_alleles",
                "clinical_cleared_transcripts",
                "clinical_trial_source_ids",
            ]
        ],
        on="entity_id",
        suffixes=("_2026_11_20", "_integrated"),
        validate="one_to_one",
    )
    changed = merged[
        merged["classification_2026_11_20"].ne(merged["classification_integrated"])
    ].copy()

    rc_rejected = set(
        normal_adjudications.rare_cell_candidate_adjudications().loc[
            lambda df: df["attribution_status"].eq("confounder_rejected"), "entity_id"
        ]
    )
    gdna_lost = set(
        library_qc.tumor_gdna_background_deltas(only_applied=True).loc[
            lambda df: ~df["qualifies_tenfold_clean"].astype(bool), "entity_id"
        ]
    )
    drivers = []
    for r in changed.itertuples(index=False):
        reasons = []
        if r.entity_id in rc_rejected:
            reasons.append("#711:confounder_rejected")
        if (
            r.home_panel_735_disposition
            and r.home_panel_735_disposition != "already_excluded_upstream"
        ):
            reasons.append(f"#735:{r.home_panel_735_disposition}")
        if r.multiplex_bleed_764_disposition:
            reasons.append(f"#764:{r.multiplex_bleed_764_disposition}")
        if r.nuclear_spillover_765_disposition == "nuclear_spillover_cleared":
            reasons.append(f"#765:{r.nuclear_spillover_765_disposition}")
        if r.clinical_clearance_766_disposition:
            reasons.append(f"#766:{r.clinical_clearance_766_disposition}")
        if r.repeat_family_767_disposition == "closed_repeat_family_cleared":
            reasons.append(f"#767:{r.repeat_family_767_disposition}")
        if r.assay_boundary_768_disposition in {
            "reclassified_ra_only_admixed_tgct_witness",
            "ihc_rna_zero_cleared",
            "adult_cornea_sparse_umi_cleared",
            "germline_tf_protamine_diluted_qualified",
            "sc_soup_and_cis_spillover_cleared",
        }:
            reasons.append(f"#768:{r.assay_boundary_768_disposition}")
        if r.entity_id in gdna_lost:
            reasons.append("library_qc:tumor_gdna_background_excluded")
        if (
            pd.notna(r.tumor_max_tpm_integrated)
            and (
                pd.isna(r.tumor_max_tpm_2026_11_20)
                or float(r.tumor_max_tpm_integrated) > float(r.tumor_max_tpm_2026_11_20) + 1e-6
            )
            and str(r.tumor_top_cancer_code_integrated) in {"THYMCA", "BCC", "cSCC"}
        ):
            reasons.append(f"#737:{r.tumor_top_cancer_code_integrated}")
        drivers.append(";".join(reasons) or "integrated_evaluation")
    changed["transition_driver"] = drivers
    changed["transition"] = (
        changed["classification_2026_11_20"] + " -> " + changed["classification_integrated"]
    )
    for col in _PREVIEW_STRING_COLS:
        if col in changed.columns:
            changed[col] = changed[col].astype("object").fillna("").astype(str)
    changed = changed.sort_values(
        ["transition", "protein_group_label"], ascending=[True, True]
    ).reset_index(drop=True)

    summary = {
        "schema": "oncoref.human-integrated-policy-preview.v1",
        "baseline_policy_id": target_selection.BASELINE_2026_11_20_POLICY_ID,
        "preview_policy_id": target_selection.POLICY_ID,
        "total_proteins": len(integrated),
        "baseline_categories": baseline["category"].value_counts().to_dict(),
        "integrated_categories": integrated["category"].value_counts().to_dict(),
        "baseline_classifications": baseline["classification"].value_counts().to_dict(),
        "integrated_classifications": integrated["classification"].value_counts().to_dict(),
        "integrated_sex_restriction_codes": integrated_res.summary["sex_restriction_codes"],
        "changed_proteins": len(changed),
        "transition_counts": changed["transition"].value_counts().to_dict(),
    }
    return changed, summary


def next_policy_preview(*, recompute: bool = False) -> pd.DataFrame:
    """Return protein transitions under integrated v1.49 (#735 + #711/#712 + #737) evaluation."""
    if recompute:
        changed, _ = compute_next_policy_transitions()
        return changed
    df = get_data("human-integrated-policy-preview").copy()
    for col in _PREVIEW_STRING_COLS:
        if col in df.columns:
            df[col] = df[col].astype("object").fillna("").astype(str)
    return df


def next_policy_summary() -> dict:
    """Return summary counts and SHA-256 provenance for the integrated human policy preview."""
    meta_path = files("oncoref").joinpath("data/human-integrated-policy-preview.json")
    csv_path = Path(str(files("oncoref").joinpath("data/human-integrated-policy-preview.csv")))
    payload = json.loads(meta_path.read_text(encoding="utf-8"))
    if _sha256_bytes(csv_path.read_bytes()) != payload["sha256"]:
        raise ValueError("human-integrated-policy-preview.csv does not match its SHA-256 receipt")
    return payload


__all__ = [
    "HUMAN_BASELINE_COMPARISON_ID",
    "HUMAN_INTEGRATED_COMPARISON_ID",
    "HUMAN_INTEGRATED_POLICY_PREVIEW_ID",
    "build_human_screening_inputs",
    "compute_next_policy_transitions",
    "evaluate_human_catalog",
    "human_screening_spec",
    "next_policy_preview",
    "next_policy_summary",
    "upgraded_737_tumor_observations",
]

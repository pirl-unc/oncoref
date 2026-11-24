# Licensed under the Apache License, Version 2.0.
"""Unified human SpeciesScreeningSpec adapter & v1.49 evidence integration (#755).

Adapts human Ensembl 115 GRCh38 protein-equivalence groups, normal tissue/protein/literature
gates, donor-resolved rare-cell RNA outcomes, #711/#712 normal adjudications, #735 unmeasured
home-tissue panels, and #737 linear-TPM tumor cohorts (THYMCA, BCC, cSCC) directly into
``species_target_selection.evaluate_species_catalog(human_spec, ...)``.

Default packaged catalog calls (``target_selection.results()``, ``human.target_selection()``)
remain deterministically pinned to ``oncoref.target-selection.2026-11-20.2``, while
``evaluate_human_catalog(...)`` and ``next_policy_preview()`` evaluate both exact baseline
parity and integrated v1.49 candidate transitions through the shared multi-species engine.
"""

from __future__ import annotations

import hashlib
import json
from importlib.resources import files
from pathlib import Path

import numpy as np
import pandas as pd

from . import (
    clinical_organ_clearance,
    normal_adjudications,
    normal_library_qc,
    rare_cell_compartment_qc,
    target_selection,
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
    packaged_species_host_exposure,
    packaged_species_metadata,
)

HUMAN_BASELINE_COMPARISON_ID = "oncoref.native-gene-tpm-normal-reference.2026-11-20.1"
HUMAN_INTEGRATED_COMPARISON_ID = "oncoref.human-grch38-integrated-tpm.2026-11-22.1"
HUMAN_INTEGRATED_POLICY_PREVIEW_ID = "oncoref.human-integrated-preview.2026-11-22.1"

REQUIRED_PERIPHERAL_CONTEXT = "peripheral_normal_reference"


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def human_screening_spec(
    *,
    policy_id: str | None = None,
    include_home_panels: bool = False,
    include_737_cohorts: bool = False,
    include_multiplex_bleed_qc: bool = False,
    include_nuclear_spillover_qc: bool = False,
    include_clinical_clearances: bool = False,
) -> SpeciesScreeningSpec:
    """Construct the declarative ``SpeciesScreeningSpec`` for human Ensembl 115 GRCh38."""
    pol = species_policy("human")
    meta = packaged_species_metadata("human")
    ref = (
        ReferenceIdentity.from_dict(meta["reference"])
        if meta is not None
        else pol.primary_reference
    )
    integrated = (
        include_home_panels
        or include_737_cohorts
        or include_multiplex_bleed_qc
        or include_nuclear_spillover_qc
        or include_clinical_clearances
    )
    effective_policy = policy_id or (
        HUMAN_INTEGRATED_POLICY_PREVIEW_ID if integrated else target_selection.POLICY_ID
    )
    comp_id = HUMAN_INTEGRATED_COMPARISON_ID if integrated else HUMAN_BASELINE_COMPARISON_ID
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
    include_home_panels: bool = False,
    include_adjudications: bool = True,
    include_737_cohorts: bool = False,
    include_multiplex_bleed_qc: bool = False,
    include_nuclear_spillover_qc: bool = False,
    include_clinical_clearances: bool = False,
    extra_tumor_observations: pd.DataFrame | None = None,
    policy_id: str | None = None,
) -> dict:
    """Materialize human screening inputs in ``evaluate_species_catalog`` format."""
    spec = human_screening_spec(
        policy_id=policy_id,
        include_home_panels=include_home_panels,
        include_737_cohorts=include_737_cohorts or extra_tumor_observations is not None,
        include_multiplex_bleed_qc=include_multiplex_bleed_qc,
        include_nuclear_spillover_qc=include_nuclear_spillover_qc,
        include_clinical_clearances=include_clinical_clearances,
    )
    ts = target_selection.results()
    exposures = packaged_species_host_exposure("human")

    cleared_rare_eids: set[str] = set()
    cleared_obs_map: dict[str, dict] = {}
    if include_nuclear_spillover_qc:
        nuc_cleared = rare_cell_compartment_qc.nuclear_spillover_adjudications(only_cleared=True)
        cleared_rare_eids = set(
            nuc_cleared.loc[
                nuc_cleared["disposition"].eq("nuclear_spillover_cleared"), "entity_id"
            ].astype(str)
        )
        cleared_obs_map = rare_cell_compartment_qc.cleared_candidate_comparison_observations()[
            "candidates"
        ]
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
            for citem in cpayload.get("reproductive_compartments", []):
                comp = str(citem["compartment"])
                is_emb = citem.get("scope") == "developmental" or comp == "embryo"
                cqual = bool(citem.get("qualifies_tenfold", False))
                cval = float(citem["max_tpm"])
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
                            "compartment_qualifying_value": (
                                float(supp_val) if cqual and supp_val is not None else 0.0
                            ),
                            "max_ntpm": None,
                            "protein_status": "not_assessed",
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
                            "compartment_qualifying_value": cval if cqual else min(cval, 0.0),
                            "max_ntpm": None,
                            "protein_status": "not_assessed",
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

        if pd.notna(r.tumor_max_tpm):
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
            and pd.notna(bleed_rec.get("max_non_tgct_tumor_tpm"))
            and float(bleed_rec["max_non_tgct_tumor_tpm"]) > 0
            and str(bleed_rec.get("top_non_tgct_sample") or "")
        ):
            nt_sample = str(bleed_rec["top_non_tgct_sample"])
            if nt_sample != (str(r.tumor_top_sample) if pd.notna(r.tumor_top_sample) else ""):
                tumor_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": "oncoref_somatic_tumor_witness_764",
                        "sample_id": nt_sample,
                        "participant_id": f"somatic:{nt_sample}",
                        "cancer_code": str(bleed_rec.get("top_non_tgct_cancer_code") or ""),
                        "value": float(bleed_rec["max_non_tgct_tumor_tpm"]),
                        "comparison_id": comp_id,
                        "unit": unit,
                        "qc_pass": True,
                        "identity_resolved": True,
                    }
                )

        qual_comps = {
            c.strip()
            for c in str(r.reproductive_qualifying_compartments or "").split("|")
            if c.strip()
        }
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
                        "value": float(val),
                        "compartment_qualifying_value": (
                            float(val) if comp_qualifies else min(float(val), 0.0)
                        ),
                        "max_ntpm": item.get("max_ntpm"),
                        "protein_status": item.get("protein_status", "not_assessed"),
                        "scope_complete": bool(item.get("scope_complete", True)),
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
    include_home_panels: bool = False,
    include_adjudications: bool = True,
    include_737_cohorts: bool = False,
    include_multiplex_bleed_qc: bool = False,
    include_nuclear_spillover_qc: bool = False,
    include_clinical_clearances: bool = False,
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
        extra_tumor_observations=extra_tumor_observations,
        policy_id=policy_id,
    )
    spec = inputs.pop("spec")
    result = evaluate_species_catalog(spec, **inputs)

    # Attach #711, #712, #735, #764, #765, and #766 adjudication provenance columns for inspection
    cat = clinical_organ_clearance.annotate_catalog(result.catalog, enabled=True)
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
    nuc_disp_map: dict[str, str] = {}
    nuc_caution_map: dict[str, str] = {}
    nuc_cleared_rare: set[str] = set()
    if include_nuclear_spillover_qc:
        nuc_df = rare_cell_compartment_qc.nuclear_spillover_adjudications()
        for nr in nuc_df.itertuples(index=False):
            eid_str = str(nr.entity_id)
            nuc_disp_map[eid_str] = str(nr.disposition)
            if pd.notna(nr.nuclear_rna_caution) and str(nr.nuclear_rna_caution):
                nuc_caution_map[eid_str] = str(nr.nuclear_rna_caution)
            if (
                bool(nr.nuclear_spillover_cleared)
                and str(nr.disposition) == "nuclear_spillover_cleared"
            ):
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
    """Recompute exact protein-level transitions from folding #735, #711/#712, #737, #764, #765, and #766 into human screening."""
    baseline = target_selection.results()
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
        "baseline_policy_id": target_selection.POLICY_ID,
        "preview_policy_id": HUMAN_INTEGRATED_POLICY_PREVIEW_ID,
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

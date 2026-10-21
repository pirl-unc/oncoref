# Licensed under the Apache License, Version 2.0.
"""Absolute normal-RNA ceilings conditioned on tissue-specific negative HPA IHC.

This is a human HPA v23 research screen. The ceilings are project parameters,
not clinical safety limits. Positive IHC and other admission gates are applied
separately; this module never certifies antigen or peptide-HLA absence.
"""

from __future__ import annotations

import json
from functools import lru_cache

import pandas as pd

from .normal_tissues import (
    SCREENED_NORMAL_TISSUE_TIERS,
    normal_tissue_screening_policy,
    normal_tissue_screening_tier,
    resolve_normal_tissue,
)

_RELIABILITY = {"Uncertain": 0, "Approved": 1, "Supported": 2, "Enhanced": 3}
_AUDIT_COLUMNS = [
    "observation_index",
    "Ensembl_Gene_ID",
    "source_name",
    "source_version",
    "screening_policy_id",
    "taxonomy_id",
    "tissue",
    "tissue_concept",
    "screening_tier",
    "rna_ntpm",
    "baseline_ceiling_ntpm",
    "rna_ceiling_ntpm",
    "ihc_evidence_status",
    "ihc_reliability_used",
    "required_ihc_context_count",
    "observed_ihc_context_count",
    "missing_ihc_contexts",
    "supporting_ihc",
    "relaxed",
]


@lru_cache(maxsize=512)
def _tissue(tissue, source, version):
    # Resolve the anatomical tissue separately from its sampled cell types.
    # Kidney RNA needs every routine kidney cell context. Caudate staining
    # cannot clear the broader basal-ganglia RNA context.
    return resolve_normal_tissue(tissue, species="human", source_id=source, source_version=version)


@lru_cache(maxsize=8)
def _routine_by_tissue(contexts, source, version):
    result = {}
    for tissue, cell in contexts:
        annotation = _tissue(tissue, source, version)
        if annotation.mapping_status == "mapped":
            result.setdefault(annotation.concept_id, set()).add((tissue, cell))
    return result


def adaptive_normal_rna_limits(
    observations, *, routine_ihc_contexts, reviewed_unreliable=False, species="human"
):
    """Return one auditable ceiling per RNA row for one gene's human HPA v23 data.

    ``observations`` contains source_name/version, modality, tissue, cell_type,
    rna_ntpm, ihc_level/reliability, valid_measurement and Ensembl_Gene_ID.
    Rows must have unique indices. ``routine_ihc_contexts`` is an independent,
    complete inventory of expected (source tissue, source cell type) pairs,
    never inferred from this gene's available observations. Missing or invalid
    evidence retains the baseline. Callers must separately reject positive IHC,
    incomplete RNA, unresolved positive anatomy and reviewed-unreliable IHC.

    The returned index matches the input RNA rows. ``rna_ceiling_ntpm`` is NaN
    for permitted or unresolved anatomy. An unknown tissue never receives a
    permissive numerical cutoff. Human nTPM rules do not transfer to other species.
    """
    if species != "human":
        raise ValueError("Adaptive numerical RNA ceilings support human HPA v23 only")
    config = normal_tissue_screening_policy("ihc_adaptive_cta")
    if not observations.index.is_unique:
        raise ValueError("Observation indices must be unique")
    if observations.Ensembl_Gene_ID.nunique() > 1:
        raise ValueError("Supply observations for exactly one gene at a time")
    for modality, source, version in (
        ("bulk_rna", config["rna_source"], config["rna_source_version"]),
        ("ihc", config["ihc_source"], config["ihc_source_version"]),
    ):
        rows = observations.loc[observations.modality.eq(modality)]
        if not rows.source_name.eq(source).all() or not rows.source_version.eq(version).all():
            raise ValueError(f"Unsupported {modality} source or version")
    if not observations.modality.isin(["bulk_rna", "ihc"]).all():
        raise ValueError("Only bulk RNA and IHC observations are supported")
    pairs = tuple(routine_ihc_contexts)
    if any(
        not isinstance(pair, (tuple, list))
        or len(pair) != 2
        or not all(isinstance(v, str) for v in pair)
        for pair in pairs
    ):
        raise ValueError("Routine IHC contexts must be (tissue, cell type) string pairs")
    contexts = tuple(sorted({tuple(pair) for pair in pairs}))
    routine = _routine_by_tissue(contexts, config["ihc_source"], config["ihc_source_version"])
    ihc = observations.loc[observations.modality.eq("ihc")].copy()
    ihc["tissue_concept"] = [
        _tissue(row.tissue, row.source_name, row.source_version).concept_id
        for row in ihc.itertuples(index=False)
    ]
    grouped = dict(iter(ihc.groupby("tissue_concept")))
    records = []
    for index, row in observations.loc[observations.modality.eq("bulk_rna")].iterrows():
        annotation = _tissue(row.tissue, row.source_name, row.source_version)
        tier = normal_tissue_screening_tier(annotation, policy="ihc_adaptive_cta")
        required = routine.get(annotation.concept_id, set())
        local = grouped.get(annotation.concept_id, ihc.iloc[:0])
        seen = set(zip(local.tissue, local.cell_type))
        missing = required - seen
        baseline = (
            config[tier]["maximum_rna"] if tier in SCREENED_NORMAL_TISSUE_TIERS else float("nan")
        )
        ceiling, reliability, status = baseline, "Missing", "unavailable"
        if tier in {"permitted", "unresolved"}:
            status = tier
        elif reviewed_unreliable:
            status = "reviewed_unreliable"
        elif not required:
            status = "no_routine_ihc_contexts"
        elif local.empty:
            status = "unavailable"
        elif not local.valid_measurement.eq(True).all():
            status = "invalid_observation"
        elif local.ihc_level.isin(config["positive_ihc_levels"]).any():
            status = "positive_ihc"
        elif missing:
            status = "incomplete_ihc_contexts"
        elif not local.ihc_level.eq("Not detected").all():
            status = "not_all_negative"
        elif not local.ihc_reliability.isin(_RELIABILITY).all():
            status = "unknown_reliability"
        else:
            reliability = min(local.ihc_reliability, key=_RELIABILITY.get)
            ceiling = config["negative_ihc_rna_limits"][reliability][tier]
            status = "negative_" + reliability.lower()
        records.append(
            {
                "observation_index": index,
                "Ensembl_Gene_ID": row.Ensembl_Gene_ID,
                "source_name": row.source_name,
                "source_version": row.source_version,
                "screening_policy_id": config["policy_id"],
                "taxonomy_id": annotation.taxonomy_id,
                "tissue": row.tissue,
                "tissue_concept": annotation.concept_id,
                "screening_tier": tier,
                "rna_ntpm": row.rna_ntpm,
                "baseline_ceiling_ntpm": baseline,
                "rna_ceiling_ntpm": ceiling,
                "ihc_evidence_status": status,
                "ihc_reliability_used": reliability,
                "required_ihc_context_count": len(required),
                "observed_ihc_context_count": len(seen),
                "missing_ihc_contexts": json.dumps(sorted(missing)),
                "supporting_ihc": json.dumps(
                    local[["tissue", "cell_type", "ihc_level", "ihc_reliability"]].to_dict(
                        orient="records"
                    )
                ),
                "relaxed": bool(ceiling > baseline),
            }
        )
    return pd.DataFrame(records, columns=_AUDIT_COLUMNS).set_index("observation_index")

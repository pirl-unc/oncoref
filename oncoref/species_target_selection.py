# Licensed under the Apache License, Version 2.0.
"""Cross-species target-selection pipeline over pluggable native data sources.

Human, dog and mouse share the same three-stage screening architecture while
keeping species, reference identities, anatomical compartments and measurement
scales strictly separate:

1. Reference proteome & protein-equivalence grouping (identical amino-acid
   translations within one verified ReferenceIdentity collapse into PEG:<sha256>;
   9-mer containment against Unsafe sequences uses the species' own proteome).
2. Normal-tissue safety & shared eligibility gates (pluggable bulk/cell/protein
   safety layers plus donor-resolved rare-cell UMI grading), followed by
   symmetric 10 x max(N, D, floor) cancer (CA) and reproductive/developmental (RA)
   enrichment on one declared, reference-compatible expression scale.
3. Downstream host-anatomy protection (HostContext parameterized by comparative
   mammalian anatomy from normal-tissue-taxonomy.yaml), sex-restricted antigen
   subsets (MRA/CTA-M, FRA/CTA-F, ERA/CTA-E, SRA/CTA-MF), and peptide subtraction.

No human expression, safety call or threshold calibration is transferred across
species by orthology. Where a species lacks a reference-reconciled normal or
tumor comparison panel in the wheel, its packaged catalog rows remain explicitly
unassessed (Other:Review) with machine-readable missing-input diagnostics.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Iterable
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from . import cta_selection_policy, cta_sequence_overlap, normal_tissues, target_selection
from .expression_contract import ExpressionContract
from .reference_identity import ReferenceIdentity
from .species import Species, resolve_species

if TYPE_CHECKING:
    from .species_expression import ExpressionDataset

COMPARATIVE_COMPARISON_UNITS = ("TPM", "nTPM", "pTPM", "TMM-CPM")

SEX_RESTRICTION_CLASSES = (
    "Male-restricted",
    "Female-restricted",
    "Embryo / placenta-only",
    "Shared male & female",
)
SEX_SETS: dict[str, tuple[str, bool]] = {
    "MRA": ("Male-restricted", False),
    "FRA": ("Female-restricted", False),
    "ERA": ("Embryo / placenta-only", False),
    "SRA": ("Shared male & female", False),
    "CTA-M": ("Male-restricted", True),
    "CTA-F": ("Female-restricted", True),
    "CTA-E": ("Embryo / placenta-only", True),
    "CTA-MF": ("Shared male & female", True),
}


def _default_nonhuman_rare_cell_provenance(source_id: str, donor_id: str) -> tuple[str, str]:
    return source_id, f"{source_id}:{donor_id}"


@dataclass(frozen=True)
class SpeciesAnatomyProfile:
    """Species-applicable reproductive, mammary, thymic and disposable compartments."""

    taxon_id: int
    species: str
    core_reproductive_compartments: tuple[str, ...]
    male_accessory_compartments: tuple[str, ...]
    female_tract_compartments: tuple[str, ...]
    mammary_compartments: tuple[str, ...]
    disposable_compartments: tuple[str, ...]
    thymic_compartments: tuple[str, ...]
    unmeasured_home_tissues: tuple[str, ...]

    @property
    def extended_reproductive_compartments(self) -> tuple[str, ...]:
        return (
            *self.male_accessory_compartments,
            *self.female_tract_compartments,
            *self.mammary_compartments,
        )

    @property
    def reproductive_compartments(self) -> tuple[str, ...]:
        return (*self.core_reproductive_compartments, *self.extended_reproductive_compartments)

    @property
    def lactating_mammary_compartment(self) -> str:
        return self.mammary_compartments[0]

    @property
    def resting_mammary_compartment(self) -> str:
        return self.mammary_compartments[-1]

    @property
    def male_sex_compartments(self) -> tuple[str, ...]:
        gonad = ("testis",) if "testis" in self.core_reproductive_compartments else ()
        return (*gonad, *self.male_accessory_compartments)

    @property
    def female_sex_compartments(self) -> tuple[str, ...]:
        gonad = ("ovary",) if "ovary" in self.core_reproductive_compartments else ()
        return (*gonad, *self.female_tract_compartments, self.lactating_mammary_compartment)

    @property
    def male_adult_protected_compartments(self) -> tuple[str, ...]:
        return (*self.male_sex_compartments, self.resting_mammary_compartment)

    @property
    def female_adult_protected_compartments(self) -> tuple[str, ...]:
        return (*self.female_sex_compartments, self.resting_mammary_compartment)

    def normalize_compartments(self, values: tuple[str, ...]) -> tuple[str, ...]:
        if not isinstance(values, tuple):
            raise ValueError("Tissue choices must be tuples of names")
        aliases: dict[str, tuple[str, ...]] = {
            c.replace("_", " "): (c,) for c in self.reproductive_compartments
        }
        for c in self.reproductive_compartments:
            aliases[c] = (c,)
        if "uterine_horn" in self.female_tract_compartments:
            aliases["uterus"] = ("uterus", "uterine_horn", "endometrium")
            aliases["uterine horn"] = ("uterine_horn",)
        else:
            aliases["uterus"] = ("uterus", "endometrium")
        aliases["oviduct"] = ("fallopian tube",)
        aliases["mammary"] = self.mammary_compartments
        aliases["mammary gland"] = self.mammary_compartments
        aliases["mammary tissue"] = self.mammary_compartments
        aliases["breast"] = self.mammary_compartments
        lactating = next((c for c in self.mammary_compartments if c.startswith("lactating_")), None)
        resting = next((c for c in self.mammary_compartments if c.startswith("resting_")), None)
        if lactating:
            aliases["lactating breast"] = (lactating,)
            aliases["lactating mammary"] = (lactating,)
        if resting:
            aliases["resting breast"] = (resting,)
            aliases["resting mammary"] = (resting,)
        found: set[str] = set()
        for name in values:
            if not isinstance(name, str):
                raise ValueError("Tissue names must be strings")
            key = " ".join(name.lower().replace("_", " ").split())
            if key not in aliases:
                raise ValueError(
                    f"Unknown compartment {name!r} for species {self.species!r}; "
                    f"valid compartments: {self.reproductive_compartments}"
                )
            found.update(aliases[key])
        return tuple(c for c in self.reproductive_compartments if c in found)

    def as_dict(self) -> dict:
        return {
            **asdict(self),
            "extended_reproductive_compartments": list(self.extended_reproductive_compartments),
            "reproductive_compartments": list(self.reproductive_compartments),
            "male_sex_compartments": list(self.male_sex_compartments),
            "female_sex_compartments": list(self.female_sex_compartments),
            "male_adult_protected_compartments": list(self.male_adult_protected_compartments),
            "female_adult_protected_compartments": list(self.female_adult_protected_compartments),
        }


def anatomy_for_species(species: Species | str | int = "human") -> SpeciesAnatomyProfile:
    """Derive species-valid screening and host-protection anatomy from the central taxonomy."""
    resolved = resolve_species(species)
    scopes = normal_tissues.normal_tissue_scopes(species=resolved)
    repro = scopes["reproductive"]
    disposable_scope = scopes["permitted_normal"]

    core = tuple(c for c in ("embryo", "placenta", "testis", "ovary") if c in repro)
    male_candidates = (
        ("epididymis", "epididymis"),
        ("prostate", "prostate"),
        ("coagulating_gland", "coagulating_gland"),
        ("seminal_vesicle", "seminal vesicle"),
    )
    male_accessory = tuple(label for cid, label in male_candidates if cid in repro)
    female_candidates = (
        ("fallopian_tube", "fallopian tube"),
        ("uterus", "uterus"),
        ("uterine_horn", "uterine_horn"),
        ("endometrium", "endometrium"),
        ("cervix", "cervix"),
        ("vagina", "vagina"),
    )
    female_tract = tuple(label for cid, label in female_candidates if cid in repro)
    mammary = (
        ("lactating_breast", "resting_breast")
        if resolved.taxon_id == 9606
        else ("lactating_mammary", "resting_mammary")
    )
    disposable_candidates = (
        ("appendix", "appendix"),
        ("gallbladder", "gallbladder"),
        ("tonsil", "tonsil"),
        ("hair_follicle", "hair follicle"),
        ("cutaneous_melanocyte", "cutaneous_melanocyte"),
    )
    doc = normal_tissues.normal_tissue_taxonomy()
    applicable = {
        cid
        for cid, concept in doc["concepts"].items()
        if resolved.taxon_id in concept.get("taxon_ids", doc["taxon_ids"])
    }
    disposable = tuple(
        label
        for cid, label in disposable_candidates
        if cid in disposable_scope or (cid in {"gallbladder", "tonsil"} and cid in applicable)
    )
    thymic = ("thymus",) if "thymus" in applicable else ()
    return SpeciesAnatomyProfile(
        taxon_id=resolved.taxon_id,
        species=resolved.name,
        core_reproductive_compartments=core,
        male_accessory_compartments=male_accessory,
        female_tract_compartments=female_tract,
        mammary_compartments=mammary,
        disposable_compartments=disposable,
        thymic_compartments=thymic,
        unmeasured_home_tissues=target_selection.UNMEASURED_HOME_TISSUES,
    )


def group_reference_proteins(
    records: pd.DataFrame | Iterable[dict],
    *,
    species: Species | str | int,
    reference: ReferenceIdentity,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collapse identical qualified translations into protein-equivalence groups."""
    resolved = resolve_species(species)
    if reference.taxon_id != resolved.taxon_id:
        raise ValueError("ReferenceIdentity taxon_id does not match requested species")
    frame = records.copy() if isinstance(records, pd.DataFrame) else pd.DataFrame(list(records))
    if frame.empty or "gene_id" not in frame:
        raise ValueError("Protein grouping requires at least one record with gene_id")
    for col, default in (
        ("symbol", ""),
        ("transcript_id", ""),
        ("protein_id", ""),
        ("transcript_biotype", "protein_coding"),
        ("sequence", ""),
    ):
        if col not in frame:
            frame[col] = default
        frame[col] = frame[col].fillna("").astype(str)
    if "is_qualified" not in frame:
        frame["is_qualified"] = frame.sequence.str.len().gt(0)
    else:
        frame["is_qualified"] = frame["is_qualified"].astype(bool) & frame.sequence.str.len().gt(0)

    rep_by_gene: dict[str, tuple[str, str, bool]] = {}
    symbols_by_gene: dict[str, str] = {}
    for gid, group in frame.groupby("gene_id", sort=True):
        if not isinstance(gid, str) or not gid.strip():
            raise ValueError("gene_id must be a nonempty string")
        sym = next((s for s in group.symbol if s), gid)
        symbols_by_gene[gid] = sym
        qual = group[group.is_qualified].copy()
        if qual.empty:
            rep_by_gene[gid] = ("", "", False)
        else:
            qual["_seq_len"] = qual.sequence.str.len()
            best = qual.sort_values(
                ["_seq_len", "transcript_id", "protein_id"], ascending=[False, True, True]
            ).iloc[0]
            digest = hashlib.sha256(best.sequence.encode("utf-8")).hexdigest()
            rep_by_gene[gid] = (best.sequence, digest, True)

    grouped_genes: dict[str, list[str]] = {}
    rep_digest_by_entity: dict[str, list[str]] = {}
    qual_by_entity: dict[str, bool] = {}
    for gid, (_seq, digest, ok) in sorted(rep_by_gene.items()):
        if ok:
            eid = "PEG:" + digest[:20]
            grouped_genes.setdefault(eid, []).append(gid)
            rep_digest_by_entity[eid] = [digest]
            qual_by_entity[eid] = True
        else:
            eid = f"UNTRANSLATED:{resolved.taxon_id}:{gid}"
            grouped_genes[eid] = [gid]
            rep_digest_by_entity[eid] = []
            qual_by_entity[eid] = False

    gene_to_entity = {gid: eid for eid, gids in grouped_genes.items() for gid in gids}
    group_rows = []
    for eid, gids in sorted(grouped_genes.items()):
        gids = sorted(gids)
        labels = sorted(dict.fromkeys(symbols_by_gene[g] for g in gids))
        group_rows.append(
            {
                "taxon_id": resolved.taxon_id,
                "species": resolved.name,
                "reference_id": reference.reference_id,
                "entity_id": eid,
                "protein_group_label": "/".join(labels),
                "representative_gene_id": gids[0],
                "member_gene_ids": json.dumps(gids),
                "member_count": len(gids),
                "representative_sequence_sha256s": json.dumps(rep_digest_by_entity[eid]),
                "identity_translation_pass": qual_by_entity[eid],
            }
        )
    groups_df = pd.DataFrame(group_rows)

    seq_rows = []
    for row in frame[frame.sequence.str.len().gt(0)].itertuples(index=False):
        eid = gene_to_entity[row.gene_id]
        digest = hashlib.sha256(row.sequence.encode("utf-8")).hexdigest()
        reps = set(rep_digest_by_entity[eid])
        seq_rows.append(
            {
                "entity_id": eid,
                "gene_id": row.gene_id,
                "transcript_id": row.transcript_id or row.gene_id,
                "protein_id": row.protein_id,
                "transcript_biotype": row.transcript_biotype or "protein_coding",
                "sequence": row.sequence,
                "sequence_sha256": digest,
                "matches_representative_sequence": (pd.NA if not reps else digest in reps),
            }
        )
    seq_df = pd.DataFrame(
        seq_rows,
        columns=[
            "entity_id",
            "gene_id",
            "transcript_id",
            "protein_id",
            "transcript_biotype",
            "sequence",
            "sequence_sha256",
            "matches_representative_sequence",
        ],
    )
    seq_df["matches_representative_sequence"] = seq_df["matches_representative_sequence"].astype(
        "boolean"
    )
    seq_df.attrs["taxon_id"] = resolved.taxon_id
    seq_df.attrs["species"] = resolved.name
    seq_df.attrs["reference_id"] = reference.reference_id
    return groups_df, seq_df


def group_proteome(
    proteome,
    *,
    species: Species | str | int,
    reference: ReferenceIdentity,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build protein-equivalence groups directly from an ``orfdata.Proteome``."""
    rows = [
        {
            "gene_id": entry.gene_id,
            "symbol": getattr(entry, "gene_symbol", None)
            or getattr(entry, "gene_name", None)
            or entry.gene_id,
            "transcript_id": getattr(entry, "transcript_id", "") or "",
            "protein_id": entry.accession,
            "transcript_biotype": getattr(entry, "biotype", "protein_coding") or "protein_coding",
            "sequence": entry.sequence,
            "is_qualified": getattr(entry, "kind", "annotated") == "annotated"
            and bool(entry.sequence),
        }
        for entry in proteome
        if getattr(entry, "gene_id", None)
    ]
    return group_reference_proteins(rows, species=species, reference=reference)


def _project_complete_group_matrix(
    dataset: ExpressionDataset,
    proteins: pd.DataFrame,
    *,
    feature_to_gene: dict[str, str] | None = None,
    allow_unresolved_reference: bool = False,
) -> pd.DataFrame:
    contract: ExpressionContract = dataset.contract
    if set(proteins["taxon_id"].unique()) != {contract.taxon_id}:
        raise ValueError("ExpressionDataset taxon_id disagrees with protein inventory")
    ref_ids = set(proteins["reference_id"].unique())
    if contract.reference is not None:
        if ref_ids != {contract.reference.reference_id} and feature_to_gene is None:
            raise ValueError(
                "ExpressionDataset reference differs from protein inventory; supply an explicit feature_to_gene crosswalk"
            )
    elif not allow_unresolved_reference and feature_to_gene is None:
        raise ValueError(
            f"ExpressionDataset {contract.source_id!r} has unresolved quantification reference; "
            "reconcile reference identity or pass an explicit feature_to_gene mapping"
        )
    values = dataset.values
    if feature_to_gene is not None:
        mapped_index = values.index.map(feature_to_gene)
        keep = mapped_index.notna()
        values = values.loc[keep].copy()
        values.index = mapped_index[keep]
        if values.index.duplicated().any():
            raise ValueError("feature_to_gene maps multiple source features to the same gene_id")
    sample_cols = list(values.columns)
    gene_to_entity = {
        g: r.entity_id for r in proteins.itertuples() for g in json.loads(r.member_gene_ids)
    }
    expected_members = proteins.set_index("entity_id").member_count
    sub = values[values.index.isin(gene_to_entity)].copy()
    sub["entity_id"] = [gene_to_entity[g] for g in sub.index]
    summed = sub.groupby("entity_id")[sample_cols].sum(min_count=1)
    present_counts = sub.groupby("entity_id")[sample_cols].count()
    complete = summed.where(present_counts.eq(expected_members.reindex(summed.index), axis=0))
    return complete.reindex(proteins.entity_id)


def observations_from_expression_dataset(
    dataset: ExpressionDataset,
    proteins: pd.DataFrame,
    *,
    comparison_id: str,
    feature_to_gene: dict[str, str] | None = None,
    normal_stat: str = "p95",
    allow_unresolved_reference: bool = False,
) -> dict[str, pd.DataFrame]:
    """Project a native ``ExpressionDataset`` onto protein-group comparison tables."""
    contract: ExpressionContract = dataset.contract
    if contract.measurement.feature_level != "gene":
        raise ValueError(
            f"Gene-level comparison requires gene measurements, got {contract.measurement.feature_level!r}"
        )
    if contract.measurement.unit not in COMPARATIVE_COMPARISON_UNITS:
        raise ValueError(
            f"ExpressionDataset unit {contract.measurement.unit!r} is not in {COMPARATIVE_COMPARISON_UNITS}"
        )
    complete = _project_complete_group_matrix(
        dataset,
        proteins,
        feature_to_gene=feature_to_gene,
        allow_unresolved_reference=allow_unresolved_reference,
    )
    species_obj = resolve_species(contract.taxon_id)
    anatomy = anatomy_for_species(species_obj)
    samples = dataset.samples
    sample_cols = list(complete.columns)

    normal_rows, repro_rows, tumor_rows = [], [], []
    if contract.context == "tumor":
        for lib_id in sample_cols:
            s = samples.loc[lib_id]
            donor = (
                f"{contract.source_id}:{s['donor_id']}"
                if s["pooling"] == "single"
                and isinstance(s.get("donor_id"), str)
                and s["donor_id"].strip()
                else ""
            )
            col = complete[lib_id]
            for eid, val in col.items():
                tumor_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": contract.source_id,
                        "sample_id": str(s["library_id"]),
                        "participant_id": donor,
                        "cancer_code": str(s.get("cancer_code") or s.get("tissue") or "TUMOR"),
                        "value": val,
                        "comparison_id": comparison_id,
                        "unit": contract.measurement.unit,
                        "qc_pass": pd.notna(val),
                        "identity_resolved": pd.notna(val),
                    }
                )
    elif contract.context == "normal":
        tissue_buckets: dict[tuple[str, str], list[str]] = {}
        for lib_id in sample_cols:
            s = samples.loc[lib_id]
            stage_raw = s.get("stage") if "stage" in s else None
            raw_tissue = str(s["tissue"]).strip()
            phys_state = (
                str(s["physiological_state"]).strip().lower()
                if "physiological_state" in s and pd.notna(s["physiological_state"])
                else None
            )
            lookup_tissue = raw_tissue
            raw_lower = " ".join(raw_tissue.lower().replace("_", " ").split())
            if raw_lower in {"lactating mammary", "lactating mammary gland", "lactating breast"}:
                lookup_tissue = (
                    "breast"
                    if anatomy.lactating_mammary_compartment == "lactating_breast"
                    else "mammary gland"
                )
                phys_state = phys_state or "lactating"
            elif raw_lower in {"resting mammary", "resting mammary gland", "resting breast"}:
                lookup_tissue = (
                    "breast"
                    if anatomy.resting_mammary_compartment == "resting_breast"
                    else "mammary gland"
                )
                phys_state = phys_state or "resting"
            ann = normal_tissues.resolve_normal_tissue(
                lookup_tissue,
                species=species_obj,
                stage=stage_raw,
                source_id="anatomy",
                context="normal",
            )
            cid = ann.concept_id or raw_tissue.lower()
            label = cid.replace("_", " ")
            if ann.stage in {"embryonic", "fetal"} or "embryonic_or_fetal" in ann.scopes:
                donor = (
                    f"{contract.source_id}:{s['donor_id']}"
                    if s["pooling"] == "single"
                    and isinstance(s.get("donor_id"), str)
                    and s["donor_id"].strip()
                    else ""
                )
                for eid, val in complete[lib_id].items():
                    repro_rows.append(
                        {
                            "entity_id": eid,
                            "source_id": contract.source_id,
                            "sample_id": str(s["library_id"]),
                            "participant_id": donor,
                            "scope": "developmental",
                            "compartment": "embryo",
                            "developmental_stage": str(stage_raw or ann.stage),
                            "source_context": raw_tissue,
                            "cells": int(s.get("cells", 20)),
                            "positive_cells": int(s.get("positive_cells", 20)),
                            "value": val,
                            "comparison_id": comparison_id,
                            "unit": contract.measurement.unit,
                            "qc_pass": pd.notna(val),
                            "identity_resolved": pd.notna(val),
                        }
                    )
            elif "reproductive" in ann.scopes:
                comp = (
                    anatomy.lactating_mammary_compartment
                    if phys_state == "lactating" or ann.physiological_state == "lactating"
                    else anatomy.resting_mammary_compartment
                    if cid in {"mammary_tissue", "mammary_glandular", "breast"}
                    else cid
                    if cid in anatomy.reproductive_compartments
                    else label
                )
                donor = (
                    f"{contract.source_id}:{s['donor_id']}"
                    if s["pooling"] == "single"
                    and isinstance(s.get("donor_id"), str)
                    and s["donor_id"].strip()
                    else ""
                )
                for eid, val in complete[lib_id].items():
                    repro_rows.append(
                        {
                            "entity_id": eid,
                            "source_id": contract.source_id,
                            "sample_id": str(s["library_id"]),
                            "participant_id": donor,
                            "scope": "reproductive",
                            "compartment": comp,
                            "tissue": label,
                            "value": val,
                            "comparison_id": comparison_id,
                            "unit": contract.measurement.unit,
                            "qc_pass": pd.notna(val),
                            "identity_resolved": pd.notna(val),
                        }
                    )
            elif cid == "thymus":
                tissue_buckets.setdefault(("thymus", "thymic"), []).append(lib_id)
            elif label in anatomy.disposable_compartments or cid in anatomy.disposable_compartments:
                tissue_buckets.setdefault((label, "disposable"), []).append(lib_id)
            elif ann.mapping_status == "mapped":
                tissue_buckets.setdefault((label, "normal"), []).append(lib_id)

        for (context_id, scope), libs in sorted(tissue_buckets.items()):
            mat = complete[libs].T
            donor_keys = [
                str(samples.loc[lib, "donor_id"])
                if samples.loc[lib, "pooling"] == "single"
                and isinstance(samples.loc[lib].get("donor_id"), str)
                and samples.loc[lib, "donor_id"].strip()
                else f"__lib__:{lib}"
                for lib in libs
            ]
            by_donor = mat.groupby(donor_keys, sort=False).median()
            summarized = target_selection.summarize_normal_tissue_tpm(by_donor, method=normal_stat)
            for eid, val in summarized.items():
                normal_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": contract.source_id,
                        "context_id": context_id,
                        "tissue": context_id,
                        "scope": scope,
                        "value": val,
                        "comparison_id": comparison_id,
                        "unit": contract.measurement.unit,
                        "qc_pass": pd.notna(val),
                        "identity_resolved": pd.notna(val),
                    }
                )
    else:
        raise ValueError(f"Unsupported ExpressionDataset context {contract.context!r}")

    return {
        "normals": pd.DataFrame(normal_rows),
        "reproductive": pd.DataFrame(repro_rows),
        "tumors": pd.DataFrame(tumor_rows),
    }


def assess_normal_rna_dataset_gates(
    dataset: ExpressionDataset,
    proteins: pd.DataFrame,
    *,
    tier_ceilings: dict[str, float] | None = None,
    feature_to_gene: dict[str, str] | None = None,
    normal_stat: str = "p95",
    allow_unresolved_reference: bool = False,
) -> pd.DataFrame:
    """Evaluate tiered normal somatic RNA ceilings directly from a species ``ExpressionDataset``."""
    from . import tissue_tiers

    contract = dataset.contract
    if contract.context != "normal":
        raise ValueError("Normal RNA safety gates require a normal ExpressionDataset")
    if contract.measurement.unit not in COMPARATIVE_COMPARISON_UNITS:
        raise ValueError(
            f"Unsupported unit {contract.measurement.unit!r} for normal RNA safety gate"
        )
    ceilings = tier_ceilings or {
        "vital": cta_selection_policy.TIER_BASES["vital"],
        "protected": cta_selection_policy.TIER_BASES["protected"],
        "manageable_hold": cta_selection_policy.TIER_BASES["manageable"],
        "manageable_unsafe": cta_selection_policy.MANAGEABLE_UNSAFE_FLOOR,
    }
    complete = _project_complete_group_matrix(
        dataset,
        proteins,
        feature_to_gene=feature_to_gene,
        allow_unresolved_reference=allow_unresolved_reference,
    )
    species_obj = resolve_species(contract.taxon_id)
    anatomy = anatomy_for_species(species_obj)
    samples = dataset.samples

    organ_libs: dict[tuple[str, str], list[str]] = {}
    for lib_id in complete.columns:
        s = samples.loc[lib_id]
        stage_raw = s.get("stage") if "stage" in s else None
        ann = normal_tissues.resolve_normal_tissue(
            s["tissue"], species=species_obj, stage=stage_raw, source_id="anatomy", context="normal"
        )
        if ann.mapping_status != "mapped" or ann.stage in {"embryonic", "fetal"}:
            continue
        if "reproductive" in ann.scopes or ann.concept_id == "thymus":
            continue
        label = (ann.concept_id or str(s["tissue"]).strip().lower()).replace("_", " ")
        if (
            label in anatomy.disposable_compartments
            or ann.concept_id in anatomy.disposable_compartments
        ):
            continue
        tier = tissue_tiers.consequences(label, None, stage="adult", resolution="bulk")["tier"]
        if tier in {"vital", "protected", "manageable"}:
            organ_libs.setdefault((label, tier), []).append(lib_id)

    rows = []
    unit = contract.measurement.unit
    for (organ, tier), libs in sorted(organ_libs.items()):
        mat = complete[libs].T
        donor_keys = [
            str(samples.loc[lib, "donor_id"])
            if samples.loc[lib, "pooling"] == "single"
            and isinstance(samples.loc[lib].get("donor_id"), str)
            and samples.loc[lib, "donor_id"].strip()
            else f"__lib__:{lib}"
            for lib in libs
        ]
        by_donor = mat.groupby(donor_keys, sort=False).median()
        summarized = target_selection.summarize_normal_tissue_tpm(by_donor, method=normal_stat)
        for eid, val in summarized.items():
            if pd.isna(val):
                continue
            v = float(val)
            if tier in {"vital", "protected"} and v > ceilings[tier]:
                rows.append(
                    {
                        "entity_id": eid,
                        "layer": f"{contract.source_id} normal bulk RNA",
                        "unsafe": True,
                        "hold": False,
                        "reason": f"{organ} ({tier}): {v:g} {unit} > {ceilings[tier]:g}",
                    }
                )
            elif tier == "manageable" and v > ceilings["manageable_unsafe"]:
                rows.append(
                    {
                        "entity_id": eid,
                        "layer": f"{contract.source_id} normal bulk RNA",
                        "unsafe": True,
                        "hold": False,
                        "reason": f"{organ} (manageable): {v:g} {unit} > {ceilings['manageable_unsafe']:g}",
                    }
                )
            elif tier == "manageable" and v > ceilings["manageable_hold"]:
                rows.append(
                    {
                        "entity_id": eid,
                        "layer": f"{contract.source_id} normal bulk RNA",
                        "unsafe": False,
                        "hold": True,
                        "reason": f"{organ} (manageable hold): {v:g} {unit} > {ceilings['manageable_hold']:g}",
                    }
                )
    return pd.DataFrame(rows, columns=["entity_id", "layer", "unsafe", "hold", "reason"])


def assess_pooled_protein_gates(
    dataset: ExpressionDataset,
    proteins: pd.DataFrame,
    *,
    feature_to_gene: dict[str, str],
    min_abundance: float = 1e5,
    hold_only: bool = False,
) -> pd.DataFrame:
    """Evaluate affirmative normal-protein detections from pooled proteomics (e.g. canine iBAQ).

    Pooled tissue channels lack individual-donor and cell-type resolution: positive
    Vital/Protected detection supplies an exclusion (or review hold when ``hold_only=True``),
    while absent/zero values never grant negative-IHC RNA ceiling multipliers.
    """
    from . import tissue_tiers

    contract = dataset.contract
    if contract.context != "normal":
        raise ValueError("Pooled protein safety gates require a normal ExpressionDataset")
    complete = _project_complete_group_matrix(
        dataset, proteins, feature_to_gene=feature_to_gene, allow_unresolved_reference=True
    )
    species_obj = resolve_species(contract.taxon_id)
    anatomy = anatomy_for_species(species_obj)
    samples = dataset.samples
    rows = []
    for lib_id in complete.columns:
        s = samples.loc[lib_id]
        ann = normal_tissues.resolve_normal_tissue(
            s["tissue"], species=species_obj, source_id="anatomy", context="normal"
        )
        if (
            ann.mapping_status != "mapped"
            or "reproductive" in ann.scopes
            or ann.concept_id == "thymus"
        ):
            continue
        label = (ann.concept_id or str(s["tissue"]).strip().lower()).replace("_", " ")
        if (
            label in anatomy.disposable_compartments
            or ann.concept_id in anatomy.disposable_compartments
        ):
            continue
        tier = tissue_tiers.consequences(label, None, stage="adult", resolution="bulk")["tier"]
        if tier not in {"vital", "protected"}:
            continue
        for eid, val in complete[lib_id].items():
            if pd.notna(val) and float(val) >= min_abundance:
                rows.append(
                    {
                        "entity_id": eid,
                        "layer": f"{contract.source_id} pooled normal protein",
                        "unsafe": not hold_only,
                        "hold": hold_only,
                        "reason": (
                            f"{label} ({tier} pooled {contract.measurement.unit}): "
                            f"{float(val):g} >= {min_abundance:g}"
                        ),
                    }
                )
    return pd.DataFrame(rows, columns=["entity_id", "layer", "unsafe", "hold", "reason"])


def classify_sex_restriction(
    table: pd.DataFrame | None = None,
    *,
    species: Species | str | int = "human",
    host_evidence: dict | None = None,
    max_tpm: float = 4.0,
    max_ntpm: float = 4.0,
) -> pd.DataFrame:
    """Partition Reproductive Antigens (RA) and CTAs into MRA/FRA/ERA/SRA and CTA-* subsets.

    Uses the species' own comparative reproductive and mammary anatomy from
    ``anatomy_for_species(species)``. Missing opposite-sex organ coverage fails
    closed (never assuming an unmeasured female or male organ is silent).
    """
    from . import host_context, target_filters

    if table is None:
        resolved = resolve_species(species)
        frame = packaged_species_results(resolved)
    else:
        frame = table.copy()
        inferred = frame.attrs.get(
            "species",
            str(frame["species"].iloc[0]) if "species" in frame and len(frame) else species,
        )
        resolved = resolve_species(inferred)

    anatomy = anatomy_for_species(resolved)
    cols = [
        "entity_id",
        "protein_group_label",
        "classification",
        "is_cta",
        "sex_restriction_class",
        "sex_restriction_code",
        "female_all_adult_eligible",
        "male_all_adult_eligible",
        "female_sex_organ_only_eligible",
        "male_sex_organ_only_eligible",
        "male_sex_organ_max_tpm",
        "female_sex_organ_max_tpm",
        "resting_breast_max_tpm",
        "tumor_max_tpm",
        "tumor_top_cancer_code",
        "reproductive_qualifying_compartments",
    ]
    ra = frame[frame.is_reproductive_antigen.astype(bool)].copy()
    if ra.empty:
        return pd.DataFrame(columns=cols)

    ra.attrs["taxon_id"] = resolved.taxon_id
    ra.attrs["species"] = resolved.name
    exposures = packaged_species_host_exposure(resolved) if host_evidence is None else host_evidence

    f_all_host = host_context.from_legacy(
        "female-patient",
        protect_tissues=anatomy.female_adult_protected_compartments,
        lactating=True,
        species=resolved,
    )
    f_sex_host = host_context.from_legacy(
        "female-patient",
        protect_tissues=anatomy.female_sex_compartments,
        lactating=True,
        species=resolved,
    )
    m_all_host = host_context.from_legacy(
        "male-patient",
        protect_tissues=anatomy.male_adult_protected_compartments,
        species=resolved,
    )
    m_sex_host = host_context.from_legacy(
        "male-patient",
        protect_tissues=anatomy.male_sex_compartments,
        species=resolved,
    )
    if (max_tpm, max_ntpm) != (4.0, 4.0):
        from dataclasses import replace

        f_all_host = replace(f_all_host, max_tpm=max_tpm, max_ntpm=max_ntpm)
        f_sex_host = replace(f_sex_host, max_tpm=max_tpm, max_ntpm=max_ntpm)
        m_all_host = replace(m_all_host, max_tpm=max_tpm, max_ntpm=max_ntpm)
        m_sex_host = replace(m_sex_host, max_tpm=max_tpm, max_ntpm=max_ntpm)

    f_all = target_filters.apply(
        target_filters.TargetFilter(host=f_all_host), ra, host_evidence=exposures
    ).set_index("entity_id")
    f_sex = target_filters.apply(
        target_filters.TargetFilter(host=f_sex_host), ra, host_evidence=exposures
    ).set_index("entity_id")
    m_all = target_filters.apply(
        target_filters.TargetFilter(host=m_all_host), ra, host_evidence=exposures
    ).set_index("entity_id")
    m_sex = target_filters.apply(
        target_filters.TargetFilter(host=m_sex_host), ra, host_evidence=exposures
    ).set_index("entity_id")

    def max_comp_tpm(eid: str, compartments: tuple[str, ...]) -> float:
        rec = exposures.get(eid, {})
        vals = [
            rec[c]["max_tpm"]
            for c in compartments
            if c in rec and rec[c].get("max_tpm") is not None
        ]
        return float(max(vals)) if vals else float("nan")

    rows = []
    for r in ra.itertuples():
        eid = r.entity_id
        fs_all = bool(f_all.at[eid, "host_status"] == "eligible")
        ms_all = bool(m_all.at[eid, "host_status"] == "eligible")
        fs_sex = bool(f_sex.at[eid, "host_status"] == "eligible")
        ms_sex = bool(m_sex.at[eid, "host_status"] == "eligible")
        if fs_all and not ms_all:
            sex_class = "Male-restricted"
            code = "CTA-M" if r.is_cta else "MRA-only"
        elif ms_all and not fs_all:
            sex_class = "Female-restricted"
            code = "CTA-F" if r.is_cta else "FRA-only"
        elif fs_all and ms_all:
            sex_class = "Embryo / placenta-only"
            code = "CTA-E" if r.is_cta else "ERA-only"
        else:
            sex_class = "Shared male & female"
            code = "CTA-MF" if r.is_cta else "SRA-only"
        rows.append(
            {
                "entity_id": eid,
                "protein_group_label": r.protein_group_label,
                "classification": r.classification,
                "is_cta": bool(r.is_cta),
                "sex_restriction_class": sex_class,
                "sex_restriction_code": code,
                "female_all_adult_eligible": fs_all,
                "male_all_adult_eligible": ms_all,
                "female_sex_organ_only_eligible": fs_sex,
                "male_sex_organ_only_eligible": ms_sex,
                "male_sex_organ_max_tpm": max_comp_tpm(eid, anatomy.male_sex_compartments),
                "female_sex_organ_max_tpm": max_comp_tpm(eid, anatomy.female_sex_compartments),
                "resting_breast_max_tpm": max_comp_tpm(eid, (anatomy.resting_mammary_compartment,)),
                "tumor_max_tpm": getattr(r, "tumor_max_tpm", np.nan),
                "tumor_top_cancer_code": getattr(r, "tumor_top_cancer_code", ""),
                "reproductive_qualifying_compartments": getattr(
                    r, "reproductive_qualifying_compartments", ""
                ),
            }
        )
    out = pd.DataFrame(rows, columns=cols).sort_values(
        ["is_cta", "sex_restriction_class", "tumor_max_tpm", "protein_group_label"],
        ascending=[False, True, False, True],
    )
    out.attrs["taxon_id"] = resolved.taxon_id
    out.attrs["species"] = resolved.name
    return out.reset_index(drop=True)


@dataclass(frozen=True)
class SpeciesScreeningSpec:
    """Complete, source-pluggable specification for running target selection in a species."""

    species: Species
    reference: ReferenceIdentity
    policy_id: str = target_selection.POLICY_ID
    comparison_id: str = ""
    comparability_reference: str = ""
    comparison_unit: str = "TPM"
    required_normal_contexts: frozenset[str] = frozenset()
    anatomy: SpeciesAnatomyProfile | None = None
    rare_cell_provenance_resolver: Callable[[str, str], tuple[str, str]] | None = None
    disposable_cap: float = target_selection.DISPOSABLE_CAP
    enrichment_fold: float = target_selection.ENRICHMENT_FOLD
    detection_floor: float = target_selection.DETECTION_FLOOR
    notes: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self):
        resolved = resolve_species(self.species)
        object.__setattr__(self, "species", resolved)
        if self.reference.taxon_id != resolved.taxon_id:
            raise ValueError("SpeciesScreeningSpec reference taxon_id disagrees with species")
        if self.comparison_unit not in COMPARATIVE_COMPARISON_UNITS:
            raise ValueError(
                f"comparison_unit must be one of {COMPARATIVE_COMPARISON_UNITS}, got {self.comparison_unit!r}"
            )
        if self.anatomy is None:
            object.__setattr__(self, "anatomy", anatomy_for_species(resolved))
        elif self.anatomy.taxon_id != resolved.taxon_id:
            raise ValueError("SpeciesAnatomyProfile taxon_id disagrees with species")
        if self.rare_cell_provenance_resolver is None and resolved.taxon_id != 9606:
            object.__setattr__(
                self, "rare_cell_provenance_resolver", _default_nonhuman_rare_cell_provenance
            )


@dataclass(frozen=True)
class SpeciesScreeningResult:
    """Evaluated species catalog, host-exposure ledger, sequence table and provenance."""

    spec: SpeciesScreeningSpec
    catalog: pd.DataFrame
    host_exposure: dict[str, dict[str, dict]]
    sequences: pd.DataFrame
    summary: dict

    def sex_restricted_antigens(self) -> pd.DataFrame:
        """Partition this species' evaluated RAs and CTAs into MRA/FRA/ERA/SRA and CTA-*."""
        return classify_sex_restriction(
            self.catalog,
            species=self.spec.species,
            host_evidence=self.host_exposure,
        )

    def select(self, group: str = "CTA", *, mode: str = "all", host=None, **overrides):
        """Enter the shared AntigenSelection / TargetFilter / PeptideIndex workflow."""
        from . import antigen_sets

        return antigen_sets.select(
            group,
            mode=mode,
            host=host,
            table=self.catalog,
            host_evidence=self.host_exposure,
            species=self.spec.species,
            sequences_table=self.sequences,
            **overrides,
        )


def evaluate_species_catalog(
    spec: SpeciesScreeningSpec,
    proteins: pd.DataFrame,
    *,
    sequences: pd.DataFrame | None = None,
    normal_gates: pd.DataFrame | None = None,
    rare_cell_counts: pd.DataFrame | None = None,
    mapping_risk: pd.DataFrame | None = None,
    confirmed_rare_cell_organs: dict[str, set[str]] | None = None,
    literature_rare_cell_organs: dict[str, set[str]] | None = None,
    comparison_normals: pd.DataFrame | None = None,
    comparison_tumors: pd.DataFrame | None = None,
    comparison_reproductive: pd.DataFrame | None = None,
) -> SpeciesScreeningResult:
    """Run the 3-stage OncoRef target-selection process for any species."""
    required_cols = {
        "entity_id",
        "protein_group_label",
        "member_gene_ids",
        "identity_translation_pass",
    }
    if missing := required_cols - set(proteins):
        raise ValueError(f"proteins missing required columns: {sorted(missing)}")
    if not proteins.entity_id.is_unique:
        raise ValueError("proteins must have unique entity_id rows")

    d = proteins.copy(deep=True).reset_index(drop=True)
    d["taxon_id"] = spec.species.taxon_id
    d["species"] = spec.species.name
    d["reference_id"] = spec.reference.reference_id
    if "representative_gene_id" not in d:
        d["representative_gene_id"] = d.member_gene_ids.map(lambda s: json.loads(s)[0])
    if "member_count" not in d:
        d["member_count"] = d.member_gene_ids.map(lambda s: len(json.loads(s)))
    if "representative_sequence_sha256s" not in d:
        d["representative_sequence_sha256s"] = "[]"
    if "home_tissue" not in d:
        d["home_tissue"] = ""
    d["home_tissue"] = d["home_tissue"].fillna("").astype(str)
    d["home_tissue_unmeasured"] = d["home_tissue"].map(
        lambda text: any(
            part.strip() in spec.anatomy.unmeasured_home_tissues
            for part in text.split(";")
            if part.strip()
        )
    )

    gate_unsafe = pd.Series(False, index=d.index)
    gate_hold = pd.Series(False, index=d.index)
    unsafe_reasons_by_entity: dict[str, list[str]] = {eid: [] for eid in d.entity_id}
    hold_reasons_by_entity: dict[str, list[str]] = {eid: [] for eid in d.entity_id}
    if normal_gates is not None and not normal_gates.empty:
        need_gate = {"entity_id", "layer", "unsafe", "hold"}
        if missing := need_gate - set(normal_gates):
            raise ValueError(f"normal_gates missing columns: {sorted(missing)}")
        for row in normal_gates.itertuples(index=False):
            if row.entity_id not in unsafe_reasons_by_entity:
                continue
            label = str(getattr(row, "reason", "") or row.layer)
            if bool(row.unsafe):
                unsafe_reasons_by_entity[row.entity_id].append(label)
            elif bool(row.hold):
                hold_reasons_by_entity[row.entity_id].append(label)
        gate_unsafe = d.entity_id.map(lambda e: bool(unsafe_reasons_by_entity[e]))
        gate_hold = d.entity_id.map(lambda e: bool(hold_reasons_by_entity[e]))

    rare_status = pd.Series("unassessed", index=d.index)
    rare_action = pd.Series("none", index=d.index)
    rare_basis = pd.Series("", index=d.index)
    if rare_cell_counts is not None and not rare_cell_counts.empty:
        by_gene = dict(iter(rare_cell_counts.groupby("gene_id")))
        risk_by_gene = (
            mapping_risk.set_index("gene_id")["mapping_risk"].astype(bool).to_dict()
            if mapping_risk is not None and not mapping_risk.empty
            else {}
        )
        confirmed_map = confirmed_rare_cell_organs or {}
        lit_map = literature_rare_cell_organs or {}
        resolver = spec.rare_cell_provenance_resolver
        rank = {
            "established": 0,
            "unresolved": 1,
            "dismissed": 2,
            "below_screen": 3,
            "unassessed": 4,
        }
        for idx, row in enumerate(d.itertuples(index=False)):
            receipts = []
            for gid in json.loads(row.member_gene_ids):
                sub = by_gene.get(gid)
                if sub is None:
                    continue
                rec = target_selection.assess_rare_cell_rna(
                    sub,
                    mapping_risk=bool(risk_by_gene.get(gid, False)),
                    confirmed_tissues=confirmed_map.get(row.entity_id, set()),
                    literature_tissues=lit_map.get(row.entity_id, set()),
                    provenance_resolver=resolver,
                )
                if rec["status"] == "unassessed":
                    rec["status"] = "below_screen"
                receipts.append(rec)
            if receipts:
                receipts.sort(
                    key=lambda r: (rank[r["status"]], -(r["max_multi_molecule_pct"] or 0))
                )
                top = receipts[0]
                rare_status.iloc[idx] = top["status"]
                rare_action.iloc[idx] = top["action"]
                if top["findings"]:
                    rare_basis.iloc[idx] = top["findings"][0]["basis"]
                if top["action"] == "unsafe":
                    unsafe_reasons_by_entity[row.entity_id].append(
                        "Established rare normal-cell RNA"
                    )
                elif top["action"] == "review":
                    hold_reasons_by_entity[row.entity_id].append("Unresolved rare normal-cell RNA")

    d["rare_cell_status"] = rare_status
    d["rare_cell_action"] = rare_action
    d["rare_cell_basis"] = rare_basis
    d["unsafe"] = gate_unsafe | rare_action.eq("unsafe")
    d["Unsafe_reasons"] = d.entity_id.map(
        lambda e: " | ".join(dict.fromkeys(unsafe_reasons_by_entity[e]))
    )
    d["normal_screen_pass"] = ~d["unsafe"] & ~gate_hold & ~rare_action.eq("review")

    seq_df = (
        sequences.copy()
        if sequences is not None
        else pd.DataFrame(
            columns=[
                "entity_id",
                "gene_id",
                "transcript_id",
                "protein_id",
                "transcript_biotype",
                "sequence",
                "sequence_sha256",
                "matches_representative_sequence",
            ]
        )
    )
    if not seq_df.empty:
        seqs_by_entity = seq_df.groupby("entity_id")["sequence"].agg(set).to_dict()
        ninemers = [
            frozenset(
                seq[i : i + cta_sequence_overlap.K]
                for seq in seqs_by_entity.get(eid, ())
                for i in range(len(seq) - cta_sequence_overlap.K + 1)
            )
            for eid in d.entity_id
        ]
        overlap_fn = cta_sequence_overlap.unsafe_overlap(ninemers, d.index)
        d["unsafe_sequence_overlap_pct"] = (100.0 * overlap_fn(d["unsafe"])).round(2)
    else:
        d["unsafe_sequence_overlap_pct"] = np.nan

    nested = d["nested_unqualified"].astype(bool) if "nested_unqualified" in d else False
    readthrough = d["readthrough_locus"].astype(bool) if "readthrough_locus" in d else False
    contained = d["unsafe_sequence_overlap_pct"].fillna(0.0).ge(95.0) & ~d["unsafe"]
    d["sequence_pass"] = ~(nested | readthrough | contained)
    d["evidence_resolved"] = (
        ~d["unresolved_evidence"].astype(bool) if "unresolved_evidence" in d else True
    )

    has_panel = (
        bool(spec.required_normal_contexts)
        and bool(spec.comparison_id.strip())
        and comparison_normals is not None
        and not comparison_normals.empty
    )
    norm_by_entity = (
        dict(iter(comparison_normals.groupby("entity_id")))
        if comparison_normals is not None and not comparison_normals.empty
        else {}
    )
    tum_by_entity = (
        dict(iter(comparison_tumors.groupby("entity_id")))
        if comparison_tumors is not None and not comparison_tumors.empty
        else {}
    )
    rep_by_entity = (
        dict(iter(comparison_reproductive.groupby("entity_id")))
        if comparison_reproductive is not None and not comparison_reproductive.empty
        else {}
    )

    empty_tumor = pd.DataFrame(
        columns=[
            "entity_id",
            "source_id",
            "sample_id",
            "participant_id",
            "value",
            "comparison_id",
            "unit",
            "qc_pass",
            "identity_resolved",
        ]
    )
    empty_repro = empty_tumor.assign(scope=pd.Series(dtype=str))

    req_cov = []
    ca_status, ra_status = [], []
    norm_ref, disp_max = [], []
    norm_contexts_json = []
    tum_max, tum_top_code, tum_fold, samples_10x, patients_10x, patient_folds_json = (
        [],
        [],
        [],
        [],
        [],
        [],
    )
    rep_max, rep_supp, rep_fold, rep_qual_comps = [], [], [], []
    bounds_pass = []
    host_exposure: dict[str, dict[str, dict]] = {}

    for row in d.itertuples(index=False):
        eid = row.entity_id
        host_exposure[eid] = {}
        n_df = norm_by_entity.get(eid)
        eligible_for_comparison = (
            not row.unsafe
            and row.normal_screen_pass
            and bool(row.identity_translation_pass)
            and row.evidence_resolved
            and row.sequence_pass
        )
        if not has_panel or n_df is None:
            req_cov.append(False)
            ca_status.append("unassessed")
            ra_status.append("unassessed")
            norm_ref.append(np.nan)
            disp_max.append(np.nan)
            norm_contexts_json.append("[]")
            tum_max.append(np.nan)
            tum_top_code.append("")
            tum_fold.append(np.nan)
            samples_10x.append(0)
            patients_10x.append(0)
            patient_folds_json.append("[]")
            rep_max.append(np.nan)
            rep_supp.append(np.nan)
            rep_fold.append(np.nan)
            rep_qual_comps.append("")
            bounds_pass.append(True)
            continue

        valid_norm = n_df[
            n_df.scope.eq("normal") & n_df.qc_pass & n_df.identity_resolved & n_df.value.notna()
        ]
        covered = set(valid_norm.context_id)
        panel_complete = set(spec.required_normal_contexts) <= covered
        req_cov.append(panel_complete)
        if not panel_complete or not eligible_for_comparison:
            ca_status.append("unassessed" if not panel_complete else "does_not_qualify")
            ra_status.append("unassessed" if not panel_complete else "does_not_qualify")
            norm_ref.append(np.nan)
            disp_max.append(np.nan)
            norm_contexts_json.append("[]")
            tum_max.append(np.nan)
            tum_top_code.append("")
            tum_fold.append(np.nan)
            samples_10x.append(0)
            patients_10x.append(0)
            patient_folds_json.append("[]")
            rep_max.append(np.nan)
            rep_supp.append(np.nan)
            rep_fold.append(np.nan)
            rep_qual_comps.append("")
            bounds_pass.append(True)
            continue

        t_df = tum_by_entity.get(eid, empty_tumor.assign(entity_id=eid))
        if t_df.empty:
            t_df = pd.DataFrame(
                [
                    {
                        "entity_id": eid,
                        "source_id": "none",
                        "sample_id": "none",
                        "participant_id": "",
                        "value": np.nan,
                        "comparison_id": spec.comparison_id,
                        "unit": spec.comparison_unit,
                        "qc_pass": False,
                        "identity_resolved": False,
                    }
                ]
            )
        ca_rec = target_selection.assess_cancer_rna(
            t_df,
            n_df,
            required_normal_contexts=set(spec.required_normal_contexts),
            comparison_id=spec.comparison_id,
            comparability_reference=spec.comparability_reference or spec.comparison_id,
            unit=spec.comparison_unit,
        )
        r_df = rep_by_entity.get(eid, empty_repro)
        if r_df.empty:
            r_eval = pd.DataFrame(
                [
                    {
                        "entity_id": eid,
                        "source_id": "none",
                        "sample_id": "none",
                        "participant_id": "",
                        "scope": "reproductive",
                        "value": np.nan,
                        "comparison_id": spec.comparison_id,
                        "unit": spec.comparison_unit,
                        "qc_pass": False,
                        "identity_resolved": False,
                    }
                ]
            )
        else:
            r_eval = r_df
        ra_rec = target_selection.assess_reproductive_rna(
            r_eval,
            n_df,
            required_normal_contexts=set(spec.required_normal_contexts),
            comparison_id=spec.comparison_id,
            comparability_reference=spec.comparability_reference or spec.comparison_id,
            unit=spec.comparison_unit,
        )

        n_max = float(ca_rec["normal_maximum"])
        d_max = float(ca_rec["disposable_maximum"])
        denom = max(n_max, d_max, spec.detection_floor)
        valid_t = t_df[t_df.qc_pass & t_df.identity_resolved & t_df.value.notna()]
        t_peak = float(valid_t.value.max()) if not valid_t.empty else np.nan
        t_code = (
            str(valid_t.sort_values("value", ascending=False).iloc[0].get("cancer_code", ""))
            if not valid_t.empty
            else ""
        )
        t_ratio = t_peak / denom if pd.notna(t_peak) else np.nan

        valid_r = (
            r_df[r_df.qc_pass & r_df.identity_resolved & r_df.value.notna()]
            if not r_df.empty
            else r_df
        )
        r_peak = float(valid_r.value.max()) if not valid_r.empty else np.nan
        r_supp = ra_rec["supported_maximum"]
        r_ratio = float(r_supp) / denom if r_supp is not None and pd.notna(r_supp) else np.nan

        qual_comps: list[str] = []
        if not valid_r.empty and "compartment" in valid_r:
            threshold = spec.enrichment_fold * denom
            for comp, cgroup in valid_r.groupby("compartment", sort=True):
                cmax = float(cgroup.value.max())
                is_dev = cgroup.scope.eq("developmental").any()
                host_exposure[eid][str(comp)] = {
                    "coverage": "measured",
                    "max_tpm": cmax,
                    "max_ntpm": None,
                    "protein_status": "not_assessed",
                    "scope_complete": not is_dev,
                }
                if cmax >= threshold:
                    qual_comps.append(str(comp))

        peaks = [v for v in (t_peak, r_supp) if v is not None and pd.notna(v)]
        target_peak = max(peaks) if peaks else 0.0
        b_ok = d_max <= spec.disposable_cap and (
            (ca_rec["status"] != "qualifies" and ra_rec["status"] != "qualifies")
            or (target_peak >= spec.detection_floor and target_peak > max(n_max, d_max))
        )

        ca_status.append(
            "qualifies"
            if ca_rec["status"] == "qualifies"
            else (
                "unassessed" if valid_t.empty and comparison_tumors is None else "does_not_qualify"
            )
        )
        ra_status.append(
            "qualifies"
            if ra_rec["status"] == "qualifies"
            else (
                "unassessed"
                if valid_r.empty and comparison_reproductive is None
                else "does_not_qualify"
            )
        )
        norm_ref.append(n_max)
        disp_max.append(d_max)
        norm_contexts_json.append(
            json.dumps(
                [
                    {
                        "tissue": r.get("tissue", r["context_id"]),
                        "context_id": r["context_id"],
                        "value": r["value"],
                    }
                    for r in ca_rec["normal_maximum_observations"]
                ]
            )
        )
        tum_max.append(t_peak)
        tum_top_code.append(t_code)
        tum_fold.append(t_ratio)
        samples_10x.append(ca_rec["samples_at_tenfold"])
        patients_10x.append(ca_rec["independent_participants_at_tenfold"])
        patient_folds_json.append(json.dumps(ca_rec["independent_participant_folds"]))
        rep_max.append(r_peak)
        rep_supp.append(r_supp if r_supp is not None else np.nan)
        rep_fold.append(r_ratio)
        rep_qual_comps.append(
            " | ".join(c for c in spec.anatomy.reproductive_compartments if c in qual_comps)
            or " | ".join(sorted(set(qual_comps)))
        )
        bounds_pass.append(bool(b_ok))

    d["required_coverage_complete"] = req_cov
    d["expression_bounds_pass"] = bounds_pass
    d["cancer_evidence"] = ca_status
    d["reproductive_evidence"] = ra_status
    d["normal_reference_tpm"] = norm_ref
    d["normal_reference_contexts"] = norm_contexts_json
    d["disposable_max_tpm"] = disp_max
    d["tumor_max_tpm"] = tum_max
    d["tumor_top_cancer_code"] = tum_top_code
    d["tumor_fold"] = tum_fold
    d["samples_at_tenfold"] = samples_10x
    d["independent_patients_at_tenfold"] = patients_10x
    d["independent_patient_tumor_folds"] = patient_folds_json
    d["reproductive_max_sample_tpm"] = rep_max
    d["reproductive_supported_tpm"] = rep_supp
    d["reproductive_fold"] = rep_fold
    d["reproductive_qualifying_compartments"] = rep_qual_comps

    classified = target_selection.classify(d)
    classified["policy"] = spec.policy_id
    classified["target_policy_id"] = spec.policy_id

    core_set = set(spec.anatomy.core_reproductive_compartments)
    ext_set = set(spec.anatomy.extended_reproductive_compartments)
    tiers = []
    for row in classified.itertuples(index=False):
        if not row.is_cta:
            tiers.append("")
            continue
        comps = {c for c in str(row.reproductive_qualifying_compartments).split(" | ") if c}
        has_core = bool(comps & core_set)
        ext_max = max(
            (host_exposure[row.entity_id].get(c, {}).get("max_tpm", 0.0) or 0.0 for c in ext_set),
            default=0.0,
        )
        if not has_core:
            tiers.append("extended_only")
        elif ext_max > spec.disposable_cap or (
            pd.notna(row.tumor_max_tpm) and ext_max > row.tumor_max_tpm / spec.enrichment_fold
        ):
            tiers.append("core_with_extended_exposure")
        else:
            tiers.append("core")
    classified["cta_tier"] = tiers

    classified.attrs.update(
        taxon_id=spec.species.taxon_id,
        species=spec.species.name,
        reference_id=spec.reference.reference_id,
        policy_id=spec.policy_id,
        comparison_id=spec.comparison_id,
        comparison_unit=spec.comparison_unit,
    )
    sex_df = classify_sex_restriction(classified, species=spec.species, host_evidence=host_exposure)
    sex_class_by_eid = (
        sex_df.set_index("entity_id")["sex_restriction_class"].to_dict() if not sex_df.empty else {}
    )
    sex_code_by_eid = (
        sex_df.set_index("entity_id")["sex_restriction_code"].to_dict() if not sex_df.empty else {}
    )
    classified["sex_restriction_class"] = classified["entity_id"].map(sex_class_by_eid).fillna("")
    classified["sex_restriction_code"] = classified["entity_id"].map(sex_code_by_eid).fillna("")

    if not seq_df.empty:
        cat_by_eid = classified.set_index("entity_id")["category"].to_dict()
        seq_df["catalog_category"] = seq_df["entity_id"].map(cat_by_eid).fillna("")
        seq_df["selection_category"] = seq_df["catalog_category"]
        seq_df["set"] = seq_df["catalog_category"]
        seq_df["background_reason"] = ""
    seq_df.attrs.update(
        policy_id=spec.policy_id,
        taxon_id=spec.species.taxon_id,
        species=spec.species.name,
        reference=spec.reference.reference_id,
        missing_gene_ids=sorted(
            {
                g
                for r in classified.itertuples(index=False)
                for g in json.loads(r.member_gene_ids)
                if g not in set(seq_df["gene_id"])
            }
        ),
        sequence_scope="all annotated protein translations",
    )

    summary = {
        "taxon_id": spec.species.taxon_id,
        "species": spec.species.name,
        "reference_id": spec.reference.reference_id,
        "policy_id": spec.policy_id,
        "comparison_id": spec.comparison_id,
        "comparison_unit": spec.comparison_unit,
        "proteins": len(classified),
        "categories": classified.category.value_counts().to_dict(),
        "classifications": classified.classification.value_counts().to_dict(),
        "sex_restriction_codes": (
            classified.loc[classified.sex_restriction_code.ne(""), "sex_restriction_code"]
            .value_counts()
            .to_dict()
        ),
    }
    return SpeciesScreeningResult(
        spec=spec,
        catalog=classified,
        host_exposure=deepcopy(host_exposure),
        sequences=seq_df,
        summary=summary,
    )


def packaged_species_metadata(species: Species | str | int = "human") -> dict | None:
    """Return packaged target-selection metadata for an evaluated species, or None if unassessed."""
    from importlib.resources import files

    resolved = resolve_species(species)
    meta_path = files("oncoref").joinpath(f"data/species/{resolved.name}/target-selection.json")
    if not meta_path.is_file():
        return None
    return json.loads(meta_path.read_text(encoding="utf-8"))


def _resolve_packaged_species_asset(
    resolved_name: str, meta: dict, path_key: str, default_stem: str
):
    from importlib.resources import files
    from pathlib import Path

    from . import reference_data

    rel = meta.get(path_key)
    if rel:
        direct = Path(str(files("oncoref").joinpath("data", rel)))
        if direct.is_file():
            return direct
        return reference_data.bundled_or_managed_path(Path(rel).name.removesuffix(".csv.gz"))
    return reference_data.bundled_or_managed_path(default_stem)


def _load_packaged_species_catalog_cached(species_name: str) -> pd.DataFrame:
    from . import reference_data

    resolved = resolve_species(species_name)
    meta = packaged_species_metadata(resolved)
    if meta is None:
        raise FileNotFoundError(f"No packaged target-selection metadata for {resolved.name!r}")
    path = _resolve_packaged_species_asset(
        resolved.name, meta, "data_path", f"{resolved.name}-target-selection"
    )
    if reference_data._sha256(path) != meta["data_sha256"]:
        raise ValueError(
            f"{resolved.name} target-selection CSV does not match its metadata SHA-256"
        )
    df = pd.read_csv(path, low_memory=False)
    if "species" not in df.columns:
        df.insert(0, "species", resolved.name)
    if "taxon_id" not in df.columns:
        df.insert(0, "taxon_id", resolved.taxon_id)
    if "policy" not in df.columns:
        df["policy"] = meta["policy_id"]
    if "target_policy_id" not in df.columns:
        df["target_policy_id"] = meta["policy_id"]
    for str_col in (
        "home_tissue",
        "Unsafe_reasons",
        "decision_reason",
        "reproductive_qualifying_compartments",
        "sex_restriction_class",
        "sex_restriction_code",
        "cta_tier",
        "tumor_top_sample",
        "rare_cell_basis",
        "subtype",
    ):
        if str_col in df.columns:
            df[str_col] = df[str_col].fillna("").astype(str)
    for bool_col in (
        "identity_translation_pass",
        "home_tissue_unmeasured",
        "unsafe",
        "normal_screen_pass",
        "sequence_pass",
        "evidence_resolved",
        "required_coverage_complete",
        "expression_bounds_pass",
        "is_cancer_antigen",
        "is_reproductive_antigen",
        "is_cta",
        "is_target_candidate",
        "published_nomination",
    ):
        if bool_col in df.columns:
            df[bool_col] = df[bool_col].fillna(False).astype(bool)
    df.attrs.update(
        taxon_id=resolved.taxon_id,
        species=resolved.name,
        reference_id=meta["reference"]["reference_id"],
        policy_id=meta["policy_id"],
        comparison_id=meta["comparison_id"],
        comparison_unit=meta["comparison_unit"],
    )
    return df


_SPECIES_CATALOG_CACHE: dict[str, pd.DataFrame] = {}
_SPECIES_HOST_EXPOSURE_CACHE: dict[str, dict] = {}
_SPECIES_SEQUENCES_CACHE: dict[str, pd.DataFrame] = {}
_SPECIES_SYMBOL_INDEX_CACHE: dict[str, dict[str, tuple[str, ...]]] = {}


def _symbol_index_for_species(species_name: str, dataset_name: str) -> dict[str, tuple[str, ...]]:
    cached = _SPECIES_SYMBOL_INDEX_CACHE.get(species_name)
    if cached is None:
        from .gene_ids import unversioned
        from .load_dataset import get_data

        cat = get_data(dataset_name, copy=False).fillna("")
        grouped: dict[str, list[str]] = {}
        for sym, gid in zip(cat["symbol"], cat["ensembl_gene_id"]):
            if isinstance(sym, str) and sym and isinstance(gid, str) and gid:
                grouped.setdefault(sym, []).append(unversioned(gid))
        cached = {k: tuple(dict.fromkeys(v)) for k, v in grouped.items()}
        _SPECIES_SYMBOL_INDEX_CACHE[species_name] = cached
    return cached


def _packaged_species_host_exposure_ref(species: Species | str | int = "human") -> dict:
    from importlib.resources import files

    from . import host_context, reference_data

    resolved = resolve_species(species)
    meta = packaged_species_metadata(resolved)
    if meta is None:
        return {}
    cached = _SPECIES_HOST_EXPOSURE_CACHE.get(resolved.name)
    if cached is None:
        path = _resolve_packaged_species_asset(
            resolved.name, meta, "host_exposure_path", f"{resolved.name}-host-exposure"
        )
        actual_sha = reference_data._sha256(path)
        if actual_sha != meta["host_exposure_sha256"]:
            raise ValueError(
                f"{resolved.name} host-exposure table does not match its metadata SHA-256"
            )
        if receipt_rel := meta.get("host_exposure_receipt_path"):
            receipt = json.loads(
                files("oncoref").joinpath("data", receipt_rel).read_text(encoding="utf-8")
            )
            if (
                receipt["policy_id"] != host_context.POLICY_ID
                or receipt["catalog_policy_id"] != meta["policy_id"]
            ):
                raise ValueError("Host exposure belongs to another policy")
            if receipt.get("catalog_data_sha256") != meta["data_sha256"]:
                raise ValueError("Host exposure belongs to another catalog snapshot")
            if receipt["sha256"] != actual_sha:
                raise ValueError("Host exposure does not match its receipt")
        cached = (
            pd.read_csv(path, keep_default_na=False)
            .set_index("entity_id")
            .exposures.map(json.loads)
            .to_dict()
        )
        _SPECIES_HOST_EXPOSURE_CACHE[resolved.name] = cached
    return cached


def packaged_species_host_exposure(species: Species | str | int = "human") -> dict:
    """Return detached per-compartment host exposure observations for ``species``."""
    return deepcopy(_packaged_species_host_exposure_ref(species))


def packaged_species_sequences(species: Species | str | int = "human") -> pd.DataFrame:
    """Return all verified reference proteome translations for an evaluated species."""
    from . import antigen_sets, reference_data

    resolved = resolve_species(species)
    meta = packaged_species_metadata(resolved)
    if meta is None:
        raise FileNotFoundError(f"No packaged proteome sequence table for {resolved.name!r}")
    if meta.get("sequence_provider") == "sequence_bundle":
        return antigen_sets._load_bundle_sequences_table()
    cached = _SPECIES_SEQUENCES_CACHE.get(resolved.name)
    if cached is None:
        path = reference_data.bundled_or_managed_path(f"{resolved.name}-proteome-sequences")
        if reference_data._sha256(path) != meta["sequences_sha256"]:
            raise ValueError(
                f"{resolved.name} proteome sequence table does not match its metadata SHA-256"
            )
        seq_df = pd.read_csv(path, low_memory=False)
        seq_df["matches_representative_sequence"] = seq_df[
            "matches_representative_sequence"
        ].astype("boolean")
        seq_df.attrs.update(
            taxon_id=resolved.taxon_id,
            species=resolved.name,
            reference_id=meta["reference"]["reference_id"],
            reference=meta["reference"]["reference_id"],
            sequence_source="packaged_proteome",
        )
        _SPECIES_SEQUENCES_CACHE[resolved.name] = seq_df
        cached = seq_df
    out = cached.copy()
    out.attrs.update(cached.attrs)
    return out


def packaged_species_screening_sources(species: Species | str | int = "mouse") -> pd.DataFrame:
    """Return the SHA-256-pinned primary input manifest for an evaluated species."""
    from importlib.resources import files

    resolved = resolve_species(species)
    path = files("oncoref").joinpath(f"data/species/{resolved.name}/screening-sources.csv")
    if not path.is_file():
        return pd.DataFrame(
            columns=["source_id", "pillar", "filename", "unit", "strain", "bytes", "sha256"]
        )
    return pd.read_csv(path)


def packaged_species_tumor_cohorts(species: Species | str | int = "mouse") -> pd.DataFrame:
    """Return the strain-explicit tumor cohort registry for an evaluated species."""
    from importlib.resources import files

    resolved = resolve_species(species)
    path = files("oncoref").joinpath(f"data/species/{resolved.name}/tumor-cohorts.csv")
    if not path.is_file():
        return pd.DataFrame(
            columns=[
                "cancer_code",
                "label",
                "model_kind",
                "strain",
                "mhc_haplotype",
                "source_accession",
                "n_samples",
                "native_feature_space",
                "quantification",
            ]
        )
    return pd.read_csv(path)


def packaged_species_results(
    species: Species | str | int = "human",
    gene: str | None = None,
    *,
    category: str | None = None,
) -> pd.DataFrame:
    """Return the packaged target-selection catalog for human, dog or mouse."""
    from .species_antigens import cta_snapshot

    resolved = resolve_species(species)
    if category is not None and category not in target_selection.CATEGORIES:
        raise ValueError(
            f"Unknown category {category!r}; expected one of {target_selection.CATEGORIES}"
        )
    meta = packaged_species_metadata(resolved)
    if meta is not None:
        cached = _SPECIES_CATALOG_CACHE.get(resolved.name)
        if cached is None:
            cached = _load_packaged_species_catalog_cached(resolved.name)
            _SPECIES_CATALOG_CACHE[resolved.name] = cached
        evaluated = cached.copy()
        evaluated.attrs.update(cached.attrs)
        if gene is not None:
            alias_gene_ids: set[str] = set()
            for cand in cta_snapshot(resolved).candidates():
                if gene in {
                    cand["id"],
                    cand["label"],
                    f"nomination:{resolved.taxon_id}:{cand['id']}",
                }:
                    alias_gene_ids.update(cand["gene_ids"])
            if sym_ds := meta.get("symbol_index"):
                alias_gene_ids.update(
                    _symbol_index_for_species(resolved.name, sym_ds).get(gene, ())
                )
            target_tokens = {gene, *alias_gene_ids}
            evaluated = evaluated[
                evaluated.protein_group_label.isin(target_tokens)
                | evaluated.protein_group_label.str.split("/").map(
                    lambda parts: bool({p for p in parts if len(p) >= 3} & target_tokens)
                )
                | evaluated.entity_id.isin(target_tokens)
                | evaluated.member_gene_ids.map(
                    lambda ids: bool(set(json.loads(ids)) & target_tokens)
                )
            ]
        if category is not None:
            evaluated = evaluated[evaluated.category.eq(category)]
        out = evaluated.reset_index(drop=True)
        out.attrs.update(cached.attrs)
        return out

    snapshot = cta_snapshot(resolved)
    native = snapshot.as_dict()
    ref_dict = native["references"][0]["identity"]
    reference = ReferenceIdentity.from_dict(ref_dict)
    rows = []
    for cand in snapshot.candidates():
        gids = list(cand["gene_ids"])
        mapping = cand["mapping"]
        complete_id = mapping["status"] == "unique" and bool(mapping["complete"]) and bool(gids)
        rows.append(
            {
                "entity_id": f"nomination:{resolved.taxon_id}:{cand['id']}",
                "protein_group_label": cand["label"],
                "representative_gene_id": gids[0] if gids else cand["id"],
                "member_gene_ids": json.dumps(gids),
                "member_count": len(gids),
                "identity_translation_pass": complete_id,
                "home_tissue": "",
                "decision_reason": "; ".join(cand["evidence_coverage"]["missing"]),
                "published_nomination": cand["nomination_kind"] == "published_nomination",
            }
        )
    spec = SpeciesScreeningSpec(
        species=resolved,
        reference=reference,
        notes=tuple(snapshot.coverage["limitations"]),
    )
    evaluated = evaluate_species_catalog(spec, pd.DataFrame(rows)).catalog
    if gene is not None:
        evaluated = evaluated[
            evaluated.protein_group_label.eq(gene)
            | evaluated.entity_id.eq(gene)
            | evaluated.entity_id.str.endswith(f":{gene}")
            | evaluated.member_gene_ids.map(lambda ids: gene in json.loads(ids))
        ]
    if category is not None:
        evaluated = evaluated[evaluated.category.eq(category)]
    return evaluated.reset_index(drop=True)


def packaged_species_coverage(species: Species | str | int = "human") -> dict:
    """Structured readiness, policy parameters, anatomy and category counts."""
    from .expression_registry import expression_sources
    from .species_antigens import cta_snapshot
    from .species_assets import registered_references
    from .species_policies import species_policy

    resolved = resolve_species(species)
    anatomy = anatomy_for_species(resolved)
    pol = species_policy(resolved)
    try:
        refs = [r["id"] for r in registered_references(resolved)]
    except FileNotFoundError:
        refs = []
    if not refs:
        refs = list(pol.registered_reference_aliases)
    pkg_meta = packaged_species_metadata(resolved)
    if pkg_meta is not None:
        return {
            "taxon_id": resolved.taxon_id,
            "species": resolved.name,
            "selection_status": pkg_meta["selection_status"],
            "policy_id": pkg_meta["policy_id"],
            "comparison_id": pkg_meta["comparison_id"],
            "comparison_unit": pkg_meta["comparison_unit"],
            "registered_references": refs,
            "reference": pkg_meta["reference"],
            "anatomy": anatomy.as_dict(),
            "thresholds": asdict(pol.thresholds),
            "counts": pkg_meta["summary"]["categories"],
            "sex_restricted_counts": pkg_meta["summary"].get("sex_restriction_codes", {}),
            "populated_inputs": [b.role for b in pol.datasets if b.status == "active_in_catalog"],
            "missing_inputs": [
                f"{b.role}: {b.source_id} ({b.status})"
                for b in pol.datasets
                if b.status != "active_in_catalog"
            ],
        }

    snap = cta_snapshot(resolved)
    cat = packaged_species_results(resolved)
    expr = expression_sources(species=resolved)
    ready_expr = [s.id for s in expr if s.acquisition_status == "available"]
    pending_expr = [s.id for s in expr if s.acquisition_status != "available"]
    missing = [
        f"{b.role}: {b.source_id} ({b.status})"
        for b in pol.datasets
        if b.status != "active_in_catalog"
    ]
    return {
        "taxon_id": resolved.taxon_id,
        "species": resolved.name,
        "selection_status": "pending_reconciled_panels",
        "policy_id": pol.policy_id,
        "registered_references": refs,
        "reference": snap.as_dict()["references"][0]["identity"],
        "anatomy": anatomy.as_dict(),
        "thresholds": asdict(pol.thresholds),
        "counts": {c: int(cat.category.eq(c).sum()) for c in target_selection.CATEGORIES},
        "sex_restricted_counts": {},
        "populated_inputs": [
            f"registered_references ({', '.join(refs)})",
            f"curated_cta_nominations ({len(cat)})",
            *([f"native_expression_matrices ({', '.join(ready_expr)})"] if ready_expr else []),
            *(
                [f"metadata_or_blocked_expression ({', '.join(pending_expr)})"]
                if pending_expr
                else []
            ),
        ],
        "missing_inputs": missing,
    }

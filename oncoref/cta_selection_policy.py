"""Current CTA selection rules, separate from frozen release policies.

RNA has one tissue/IHC ceiling. Crossing it supplies exclusion evidence in
Vital/Protected compartments; crossing it in Manageable compartments prevents
CTA admission but does not establish either Unsafe or Neither. Neither requires
an affirmative expression profile. Unresolved results have a null category.
These are research selection rules, not estimates of clinical toxicity.
"""

from __future__ import annotations

from collections.abc import Collection, Sequence

import numpy as np
import pandas as pd

from .cell_rna_policy import EXPECTED_UNITS, GRADES, NEUTRAL_REPRODUCTIVE_IMMUNE, POSITIVE
from .compartment_rna_policy import IHC_RANKS
from .cta_funnel import EvidenceLayer, FunnelResult

POLICY_ID = "oncoref.cta-selection.2026-10-28"
TIER_BASES = {"vital": 4.0, "protected": 8.0, "manageable": 16.0}
APPRECIABLE_RNA = 2.0
MIN_MANAGEABLE_TISSUES = 2
UNIT_SCALE = {"TPM": 1.0, "nTPM": 1.0, "CP10K": 0.01}
BOWEL_EPITHELIAL_CELLS = frozenset(
    {
        "distal enterocytes",
        "intestinal goblet cells",
        "paneth cells",
        "undifferentiated cells",
        "glandular cells",
    }
)


def compartment_tiers(contexts: pd.DataFrame) -> pd.Series:
    """Use the same identified cell-compartment tier for RNA and protein.

    Whole-organ bulk retains its anatomical tier. This policy treats identified
    fibroblasts and large-bowel epithelial populations as Manageable; it does
    not describe them as disposable or establish tolerability of cell killing.
    """
    tiers = contexts.category.replace({"renewing_mucosa": "manageable"}).copy()
    mapped = contexts.mapping_status.eq("mapped") & contexts.category.isin(
        [*TIER_BASES, "renewing_mucosa"]
    )
    cellular = contexts.source.isin(["hpa_single_cell_tissue", "hpa_normal_tissue"])
    cells, tissues = contexts.cell_type.str.lower(), contexts.tissue.str.lower()
    manageable = cells.eq("fibroblasts") | (
        tissues.isin(["colon", "rectum"]) & cells.isin(BOWEL_EPITHELIAL_CELLS)
    )
    tiers.loc[mapped & cellular & manageable] = "manageable"
    return tiers


def match_ihc(
    contexts: pd.DataFrame,
    ihc: pd.DataFrame,
    *,
    positive_controls: Collection[str],
    conflicting_genes: Collection[str] = (),
    crosswalk: pd.DataFrame,
) -> pd.DataFrame:
    """Match adult HPA RNA to exact or explicitly reviewed IHC compartments.

    The crosswalk is source/version/tissue/cell specific and carries provenance.
    A parent-compartment match is never labelled exact. Unreviewed same-tissue
    evidence remains visible without changing the ceiling. Sorted immune cells,
    stimulated states, other species and developmental samples cannot borrow
    adult HPA staining. Group negatives and controls must already cover all loci.
    """
    stale = [
        "Level",
        "Reliability",
        "ihc_positive",
        "ihc_exact_match",
        "negative_ihc_grade",
        "ihc_match_scope",
        "ihc_cell_type",
        "ihc_match_source",
    ]
    d = contexts.drop(columns=stale, errors="ignore").copy()
    h = ihc.rename(
        columns={"Gene": "gene_id", "Tissue": "tissue", "Cell type": "ihc_cell_type"}
    ).copy()
    for frame, cells in [(d, "cell_type"), (h, "ihc_cell_type")]:
        for name in ["tissue", cells]:
            frame[name] = frame[name].fillna("").str.strip().str.lower()
    keys = ["source", "version", "tissue", "cell_type"]
    if crosswalk.duplicated(keys).any():
        raise ValueError("Ambiguous IHC compartment crosswalk")
    if h.duplicated(["gene_id", "tissue", "ihc_cell_type"]).any():
        raise ValueError("Repeated IHC gene/tissue/cell observation")
    d = d.merge(crosswalk, on=keys, how="left", validate="many_to_one")
    d["ihc_cell_type"] = d.ihc_cell_type.fillna(d.cell_type)
    d = d.merge(
        h[["gene_id", "tissue", "ihc_cell_type", "Level", "Reliability"]],
        on=["gene_id", "tissue", "ihc_cell_type"],
        how="left",
        validate="many_to_one",
    )
    allowed = (
        d.source.eq("hpa_single_cell_tissue")
        & d.taxonomy_id.eq(9606)
        & d.stage.eq("adult")
        & d.state.isin(["normal", "unspecified", "unstimulated"])
    )
    matched = allowed & d.Level.notna()
    exact = matched & d.ihc_cell_type.eq(d.cell_type)
    same_tissue = pd.MultiIndex.from_frame(d[["gene_id", "tissue"]]).isin(
        pd.MultiIndex.from_frame(h[["gene_id", "tissue"]])
    )
    d["ihc_match_scope"] = np.select(
        [exact, matched, allowed & same_tissue],
        ["exact_cell", "reviewed_compartment", "same_tissue_unresolved_cell"],
        default="unmatched",
    )
    d["ihc_exact_match"] = exact
    d["ihc_positive"] = matched & d.Level.isin(POSITIVE) & d.Reliability.isin(GRADES)
    usable = (
        matched
        & d.Level.eq("Not detected")
        & d.Reliability.isin(GRADES)
        & d.gene_id.isin(positive_controls)
        & ~d.gene_id.isin(conflicting_genes)
    )
    d["negative_ihc_grade"] = d.Reliability.where(usable, "Missing/uncertain")
    return d


def assess_contexts(contexts: pd.DataFrame) -> pd.DataFrame:
    """Apply a single inclusive RNA ceiling to supported native measurements.

    Equality passes. Source-specific support statistics are supplied by adapters.
    CP10K has a nominal numerical /100 policy, not a biochemical TPM conversion.
    Fibroblasts and bowel epithelium use Manageable; organ bulk retains its tier.
    Protein observations are assessed in their own layer, avoiding double counts.
    """
    required = {
        "source",
        "taxonomy_id",
        "stage",
        "tissue",
        "cell_type",
        "category",
        "mapping_status",
        "negative_ihc_grade",
        "value",
        "unit",
    }
    if not required <= set(contexts):
        raise ValueError(f"Missing context fields: {sorted(required - set(contexts))}")
    if not contexts.taxonomy_id.eq(9606).all() or not contexts.stage.eq("adult").all():
        raise ValueError("This measurement policy requires adult human contexts")
    assay = contexts.source.where(~contexts.source.str.startswith("dice_"), "dice")
    if not contexts.unit.eq(assay.map(EXPECTED_UNITS | {"hpa_rna_consensus": "nTPM"})).all():
        raise ValueError("Source and native unit do not match")
    if not contexts.negative_ihc_grade.isin(IHC_RANKS).all():
        raise ValueError("Unknown IHC grade")
    known_tiers = {
        *TIER_BASES,
        "renewing_mucosa",
        "reproductive",
        "developmental",
        "disposable",
        "unresolved",
        "permitted",
        NEUTRAL_REPRODUCTIVE_IMMUNE,
    }
    if not contexts.category.isin(known_tiers).all():
        raise ValueError("Unknown anatomical category")
    if (
        "ihc_positive" in contexts
        and (contexts.ihc_positive.fillna(False) & contexts.negative_ihc_grade.isin(GRADES)).any()
    ):
        raise ValueError("Positive IHC cannot receive negative-IHC credit")
    d = contexts.drop(
        columns=[
            "strong_boundary",
            "reference_cutoff",
            "reference_review",
            "cta_cutoff",
            "cta_excluded",
            "local_exclusion",
        ],
        errors="ignore",
    ).copy()
    d["value"] = pd.to_numeric(d.value, errors="raise")
    if d.value.lt(0).any() or np.isinf(d.value).any():
        raise ValueError("Require finite nonnegative RNA or missing")
    d["selection_tier"] = compartment_tiers(d)
    cells = d.cell_type.fillna("").str.lower()
    stroma = d.source.eq("hpa_single_cell_tissue") & cells.eq("fibroblasts")
    normal = d.category.isin([*TIER_BASES, "renewing_mucosa"]) & d.mapping_status.eq("mapped")
    d.loc[stroma & normal, "selection_tier"] = "manageable"
    d["compartment_basis"] = np.where(
        stroma & normal, "identified_fibroblast", "anatomical_default"
    )
    d.loc[d.selection_tier.ne(d.category) & ~stroma, "compartment_basis"] = "large_bowel_epithelium"
    d["candidate_ceiling"] = (
        d.selection_tier.map(TIER_BASES)
        * 2.0 ** d.negative_ihc_grade.map(IHC_RANKS)
        * d.unit.map(UNIT_SCALE)
    )
    d["candidate_exclusion"] = normal & d.value.gt(d.candidate_ceiling)
    d["unsafe_evidence"] = d.candidate_exclusion & d.selection_tier.isin(["vital", "protected"])
    d["unresolved_context"] = ~d.mapping_status.eq("mapped")
    d["candidate_exclusion"] |= d.unresolved_context
    d["selection_policy"] = POLICY_ID
    return d


def neither_profiles(facts: pd.DataFrame) -> pd.DataFrame:
    """Identify positive Neither profiles; absence of data never counts as low.

    Low-expression: all measured normal, reproductive, embryonic and tumor RNA
    below the shared 2-unit expression floor, with complete required RNA panels,
    measured tumor and developmental coverage, and no detected protein.
    Diffuse-Manageable: no reproductive/developmental association, appreciable
    expression in >=2 distinct Manageable tissue concepts, complete required RNA
    coverage, and no above-ceiling or accepted positive Vital/Protected evidence.
    Optional DICE/cornea findings count when present; missing optional assays are
    disclosed, not treated as measurements of zero. Protein uncertainty blocks
    assignment until adjudicated. The final classifier gives Unsafe precedence.
    """
    d = facts.copy()
    complete = d.required_coverage & ~d.unresolved_evidence
    d["neither_low_expression"] = (
        complete
        & d.tumor_measured
        & d.developmental_measured
        & d.normal_maximum_floor_ratio.lt(1)
        & d.reproductive_ntpm.lt(APPRECIABLE_RNA)
        & d.developmental_all_lineage_max_tpm.lt(APPRECIABLE_RNA)
        & d.tumor_maximum_floor_ratio.lt(1)
        & ~d.any_detected_protein
    )
    d["neither_diffuse_manageable"] = (
        complete
        & ~d.positive_context
        & ~d.vital_protected_exclusion
        & d.manageable_tissue_count.ge(MIN_MANAGEABLE_TISSUES)
    )
    return d


def classify(
    facts: pd.DataFrame, layers: Sequence[EvidenceLayer], *, entity_column: str
) -> FunnelResult:
    """Classify every entity using all layers and affirmative Neither profiles.

    Null category means no supported biological assignment; assessment_status
    explains incomplete evidence versus an observed profile outside the three
    definitions. The funnel partitions all entities at every stage. Neither is
    shown only after its complete-profile evaluation, never as a synonym for hold.
    """
    if facts[entity_column].isna().any() or not facts[entity_column].is_unique:
        raise ValueError("Require unique entity identities")
    if not layers or len({layer[0] for layer in layers}) != len(layers):
        raise ValueError("Require uniquely named evidence layers")
    for _, hold, unsafe in layers:
        for mask in (hold, unsafe):
            if not mask.index.equals(facts.index) or mask.isna().any() or mask.dtype != bool:
                raise ValueError("Complete aligned boolean evidence required")
    holds = pd.DataFrame({name: hold | unsafe for name, hold, unsafe in layers})
    exclusions = pd.DataFrame({name: unsafe for name, _, unsafe in layers})
    held, unsafe = holds.any(axis=1), exclusions.any(axis=1)
    profiles = facts.neither_low_expression | facts.neither_diffuse_manageable
    neither = held & profiles & ~unsafe
    eligible = ~(held | unsafe)
    result = facts.drop(columns=[c for c in facts if c.startswith(("hold:", "unsafe:"))]).copy()
    result["category"] = pd.Series(pd.NA, index=facts.index, dtype="string")
    for label, mask in [("CTA", eligible), ("Neither", neither), ("Unsafe", unsafe)]:
        result.loc[mask, "category"] = label
    result["assessment_status"] = np.select(
        [
            result.category.notna(),
            facts.unresolved_evidence
            | ~facts.required_coverage
            | ~facts.tumor_measured
            | ~facts.developmental_measured,
        ],
        ["screened", "incomplete_evidence"],
        default="outside_category_definitions",
    )
    result["neither_profile"] = np.select(
        [neither & facts.neither_low_expression, neither & facts.neither_diffuse_manageable],
        ["low_expression", "diffuse_manageable"],
        default="",
    )
    result["exclusion_reasons"] = holds.apply(lambda row: " | ".join(row.index[row]), axis=1)
    result["Unsafe_reasons"] = exclusions.apply(lambda row: " | ".join(row.index[row]), axis=1)
    running_hold = pd.Series(False, index=facts.index)
    running_unsafe = running_hold.copy()
    sequential, independent = [], []
    for name, _, _ in layers:
        before = ~(running_hold | running_unsafe)
        previous_unsafe = running_unsafe.copy()
        running_hold |= holds[name]
        running_unsafe |= exclusions[name]
        sequential.append(
            {
                "layer": name,
                "new_candidate_removals": int((before & running_hold).sum()),
                "new_Unsafe": int((running_unsafe & ~previous_unsafe).sum()),
                "CTA_eligible_so_far": int((~(running_hold | running_unsafe)).sum()),
                "Unsafe": int(running_unsafe.sum()),
                "Neither": 0,
                "Unassigned": int((running_hold & ~running_unsafe).sum()),
            }
        )
        independent.append(
            {
                "layer": name,
                "candidate_exclusions_anywhere": int(holds[name].sum()),
                "Unsafe_evidence_anywhere": int(exclusions[name].sum()),
                "CTA_added_if_layer_removed": int(
                    (held & ~holds.drop(columns=name).any(axis=1)).sum()
                ),
                "Unsafe_reclassified_if_layer_removed": int(
                    (unsafe & ~exclusions.drop(columns=name).any(axis=1)).sum()
                ),
            }
        )
        result["hold:" + name], result["unsafe:" + name] = holds[name], exclusions[name]
    sequential.append(
        {
            "layer": "Neither expression profiles",
            "new_candidate_removals": 0,
            "new_Unsafe": 0,
            "CTA_eligible_so_far": int(eligible.sum()),
            "Unsafe": int(unsafe.sum()),
            "Neither": int(neither.sum()),
            "Unassigned": int(result.category.isna().sum()),
        }
    )
    result["policy"] = POLICY_ID
    return FunnelResult(result, pd.DataFrame(sequential), pd.DataFrame(independent))

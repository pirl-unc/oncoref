"""Current CTA selection rules, separate from frozen release policies.

RNA has one tissue/IHC ceiling. Crossing it supplies exclusion evidence in
Vital/Protected compartments; crossing it in Manageable compartments prevents
candidate admission and remains a normal-screen hold unless it exceeds the
declared Unsafe floor. Other includes assessed noncandidates and held/incomplete
cases; evidence status is independent of category and subtype.
These are research selection rules, not estimates of clinical toxicity.

Normal-tissue safety, cancer expression and reproductive /
developmental association are distinct questions. Only safety evidence can
establish Unsafe. Cancer requires one qualifying measured sample, without a
tumor-versus-normal ratio. Association requires positive expression, without a
reproductive-versus-normal maximum veto. Profile diffuseness is assessed
separately from normal-tissue safety.
"""

from __future__ import annotations

from collections.abc import Callable, Collection, Sequence

import numpy as np
import pandas as pd

from .cell_rna_policy import EXPECTED_UNITS, GRADES, NEUTRAL_REPRODUCTIVE_IMMUNE, POSITIVE
from .compartment_rna_policy import IHC_RANKS
from .cta_funnel import EvidenceLayer, FunnelResult
from .expression_diffuseness import MAX_DIFFUSENESS
from .tissue_hazards import context_tiers

POLICY_ID = "oncoref.cta-selection.2026-11-05.1"
TIER_BASES = {"vital": 4.0, "protected": 8.0, "manageable": 16.0}
APPRECIABLE_RNA = 2.0
# Manageable-tier RNA above the highest Manageable ceiling (16 x 2^3 = 128 in
# native TPM/nTPM units) is unsafe evidence; lower excess only blocks CTA.
MANAGEABLE_UNSAFE_FLOOR = 128.0
# A would-be CTA with at least this fraction of its 9-mers also present in Unsafe
# proteins is Other:Contained. Unsafe sequence is blacklisted, so such a protein
# retains almost no usable peptide.
CONTAINED_FLOOR = 0.95
OTHER_SUBTYPE_RULES = {
    "Nonreproductive": "No qualifying reproductive, placental or embryonic/fetal evidence establishes CTA identity. Evidence completeness and normal holds are recorded independently; this is not proof of absence.",
    "Reproductive": "Reproductive presence is supported but cancer presence is unestablished. Inspect assessment status for independent normal-evidence gaps or holds.",
    "Diffuse": "Reproductive presence and cancer presence pass, but the reproductive concentration gate fails. Inspect assessment status independently.",
    "Held": "Target gates pass, but an unresolved or marginal normal-tissue signal blocks CTA without establishing Unsafe.",
    "Nested": "A reference-CDS microprotein (under 50 residues) nested inside a same-strand protein-coding gene and not an established ncORF: its RNA signal belongs to the host, so it cannot be a CTA.",
    "Translation": "The selected primary-locus reference translation needs qualification; this is not normal-tissue exclusion evidence.",
    "Contained": "Every CTA gate passes, but at least 95% of its 9-mers also occur in Unsafe proteins, so almost no peptide survives the Unsafe sequence blacklist.",
}
CANDIDATE_ONLY_LAYERS = frozenset(
    {
        "Reproductive / embryonic association",
        "Tumor sample ≥2 TPM/pTPM",
        "Cancer expression",
        "Required RNA coverage",
        "Translation qualification",
        "Reproductive concentration",
    }
)
UNIT_SCALE = {"TPM": 1.0, "nTPM": 1.0, "CP10K": 0.01}
ADDITIONAL_BULK_UNITS = {"hpa_brain": "nTPM", "gtex_tissue_rna": "TPM"}
BOWEL_EPITHELIAL_CELLS = frozenset(
    {
        "distal enterocytes",
        "intestinal goblet cells",
        "paneth cells",
        "undifferentiated cells",
        "glandular cells",
    }
)


def association_evidence(facts: pd.DataFrame) -> pd.DataFrame:
    """Assess positive reproductive/developmental expression independently.

    Association uses reproductive >=5 nTPM, supported developmental >=5 TPM, or
    supported reproductive-lineage cell RNA >=5 nTPM (>=20 reads), or
    reliable positive reproductive IHC or reviewed human reproductive protein.
    Protein findings must be mapped and support every
    member locus, without an unresolved antibody concern. It supplies presence
    only: normal safety, cancer presence and concentration remain separate gates.
    Developmental cell-support criteria are applied upstream.
    Normal-tissue expression is not an input to this presence assessment;
    normal safety and whole-profile concentration answer separate questions.
    """
    d = facts.copy()
    d["reproductive_rna_association"] = d.reproductive_ntpm.ge(5)
    protein = d.get("reproductive_protein_association", pd.Series(False, index=d.index))
    if protein.isna().any() or protein.dtype != bool:
        raise ValueError("Reproductive protein evidence must be explicitly boolean")
    d["reproductive_protein_association"] = protein
    cell = d.get("reproductive_cell_association", pd.Series(False, index=d.index))
    if cell.isna().any() or cell.dtype != bool:
        raise ValueError("Reproductive cell evidence must be explicitly boolean")
    d["reproductive_cell_association"] = cell
    literature = d.get("knowledge_reproductive_support", pd.Series(False, index=d.index))
    if literature.isna().any() or literature.dtype != bool:
        raise ValueError("Reviewed reproductive knowledge must be explicitly boolean")
    d["reproductive_association"] = d.reproductive_rna_association | protein | cell | literature
    d["developmental_association"] = d.developmental_mean_tpm.ge(5)
    d["positive_context"] = d.reproductive_association | d.developmental_association
    d["reproductive_route"] = d.reproductive_association
    d["developmental_route"] = d.developmental_association
    return d


def reproductive_protein_evidence(
    facts: pd.DataFrame, observations: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Admit mapped HPA positives as association, retaining all supporting rows.

    Low/Medium/High with Approved/Supported/Enhanced reliability may attest to
    reproductive presence. Uncertain, unmapped, foreign-species and foreign-
    version rows cannot. Every locus in a protein group must have affirmative
    evidence; evidence is never copied from a related paralog. No normal-safety
    or RNA-ceiling decision is made here.
    """
    import json

    d, evidence = facts.copy(), observations.copy()
    eligible = (
        evidence.source.eq("hpa_normal_tissue")
        & evidence.version.eq("v23")
        & evidence.taxonomy_id.eq(9606)
        & evidence.mapping_status.eq("mapped")
        & evidence.category.eq("reproductive")
        & evidence.Level.isin(POSITIVE)
        & evidence.Reliability.isin(GRADES)
    )
    positive_ids = set(evidence.loc[eligible, "gene_id"])
    conflicts = d.get("Antibody concern", pd.Series(False, index=d.index))
    if conflicts.isna().any() or conflicts.dtype != bool:
        raise ValueError("Antibody concern must be explicitly boolean")
    members = d.member_gene_ids.map(json.loads)
    if members.map(lambda ids: not ids or len(set(ids)) != len(ids)).any():
        raise ValueError("Require nonempty unique member loci")
    d["reproductive_protein_association"] = (
        members.map(lambda ids: set(ids) <= positive_ids) & ~conflicts
    )
    evidence["qualifying_reproductive_positive"] = eligible
    return d, evidence.loc[evidence.category.eq("reproductive")].copy()


def concentration_evidence(facts: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """Apply the fixed HPA adult/placental identity contrast, never a safety veto.

    A complete same-source profile above 10% blocks the adult/placental route.
    An RNA-, cell- or IHC-supported adult reference may be assessed in its matched RNA
    panel; zero/missing RNA references remain explicitly unassessed rather than
    being replaced with an invented abundance. Supported embryonic evidence remains a separate
    route with developmental concentration explicitly unassessed. Missing HPA
    comparisons are reported, not imputed or substituted with a GTEx ratio.
    """
    if (
        scores.entity_id.duplicated().any()
        or not (
            scores.source.eq("hpa_rna_consensus")
            & scores.version.eq("v23")
            & scores.unit.eq("nTPM")
            & scores.taxonomy_id.eq(9606)
            & scores.summary.eq("consensus")
        ).all()
    ):
        raise ValueError("Require one human HPA v23 consensus score per protein")
    d = facts.copy()
    indexed = scores.set_index("entity_id")
    key = d.protein_equivalence_group_id
    d["hpa_diffuseness"] = key.map(indexed.diffuseness)
    d["hpa_diffuseness_reference"] = key.map(indexed.reference_value)
    complete = (
        key.map(indexed.reference_complete).fillna(False).astype(bool)
        & np.isclose(key.map(indexed.normal_coverage), 1)
        & d.hpa_diffuseness.notna()
    )
    # The association and contrast must refer to exactly the same reproductive RNA.
    comparable = (
        complete
        & d.reproductive_ntpm.gt(0)
        & np.isclose(d.hpa_diffuseness_reference, d.reproductive_ntpm)
    )
    d["adult_profile_status"] = np.select(
        [
            ~d.reproductive_association,
            ~comparable,
            d.hpa_diffuseness.gt(MAX_DIFFUSENESS),
        ],
        ["no_qualifying_adult_reference", "unassessed", "diffuse"],
        default="concentrated",
    )
    d["developmental_profile_status"] = np.where(
        d.developmental_association,
        "unassessed_comparable_panel_missing",
        "no_qualifying_reference",
    )
    d["reproductive_concentration_hold"] = (
        d.adult_profile_status.eq("diffuse") & ~d.developmental_association
    )
    d["reproductive_route"] = d.reproductive_association & ~d.adult_profile_status.eq("diffuse")
    return d


def compartment_tiers(contexts: pd.DataFrame) -> pd.Series:
    """Central, cited tissue/cell policy shared by RNA, IHC and rare-cell screens."""
    return context_tiers(contexts)


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
    Adult human tissue and cell tiers come from the central hazard registry.
    Protein observations are assessed in their own layer, avoiding double counts.
    Reviewed TCR-T evidence may add one RNA ladder step, capped at Enhanced.
    This never relabels IHC, supplies coverage, or overrides protein/clinical injury.
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
    if not contexts.unit.eq(
        assay.map(EXPECTED_UNITS | {"hpa_rna_consensus": "nTPM"} | ADDITIONAL_BULK_UNITS)
    ).all():
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
    normal = d.selection_tier.isin(TIER_BASES) & d.mapping_status.eq("mapped")
    d["compartment_basis"] = np.where(
        d.selection_tier.ne(d.category) & normal,
        "central_tissue_hazard_registry",
        "anatomical_default",
    )
    d["baseline_ceiling"] = (
        d.selection_tier.map(TIER_BASES)
        * 2.0 ** d.negative_ihc_grade.map(IHC_RANKS)
        * d.unit.map(UNIT_SCALE)
    )
    if "clinical_rna_credit" not in d:
        d["clinical_rna_credit"] = False
    if not d.clinical_rna_credit.map(lambda value: isinstance(value, (bool, np.bool_))).all():
        raise ValueError("Clinical RNA credit must be explicitly boolean")
    credit = d.clinical_rna_credit & normal
    if credit.any() and (
        "clinical_product_ids" not in d
        or d.loc[credit, "clinical_product_ids"].fillna("").eq("").any()
    ):
        raise ValueError("Clinical RNA credit requires reviewed product provenance")
    rank = d.negative_ihc_grade.map(IHC_RANKS)
    d["candidate_ceiling"] = (
        d.selection_tier.map(TIER_BASES)
        * 2.0 ** (rank + credit.astype(int)).clip(upper=3)
        * d.unit.map(UNIT_SCALE)
    )
    d["clinical_credit_applied"] = credit & d.candidate_ceiling.gt(d.baseline_ceiling)
    d["baseline_exclusion"] = normal & d.value.gt(d.baseline_ceiling)
    d["candidate_exclusion"] = normal & d.value.gt(d.candidate_ceiling)
    manageable_excess = (
        normal
        & d.selection_tier.eq("manageable")
        & d.value.div(d.unit.map(UNIT_SCALE)).gt(MANAGEABLE_UNSAFE_FLOOR)
    )
    d["unsafe_evidence"] = (
        d.candidate_exclusion & d.selection_tier.isin(["vital", "protected"])
    ) | manageable_excess
    d["unresolved_context"] = ~d.mapping_status.eq("mapped")
    d["candidate_exclusion"] |= d.unresolved_context
    d["selection_policy"] = POLICY_ID
    return d


def other_profiles(facts: pd.DataFrame) -> pd.DataFrame:
    """Describe normal profiles; final Other requires all normal gates to pass.

    Low observed normal expression uses complete required normal panels and no
    detected normal protein. Tumor abundance or missing developmental evidence
    cannot disqualify a normal-pass protein from Other. The adult Disposable
    pattern uses only adult tissue observations. Developmental expression and
    missing developmental measurements never change this adult designation.
    It is a descriptive attribute, not an admission or safety override.
    Optional DICE/cornea findings count when present; missing optional assays are
    disclosed, not treated as measurements of zero. Protein uncertainty blocks
    assignment until adjudicated. The final classifier gives Unsafe precedence.
    """
    d = facts.copy()
    complete = d.required_coverage & ~d.unresolved_evidence
    d["other_low_expression"] = (
        complete & d.normal_maximum_floor_ratio.lt(1) & ~d.any_detected_protein
    )
    d["other_disposable"] = (
        complete
        & (
            d.get("disposable_maximum_floor_ratio", pd.Series(float("nan"), index=d.index)).ge(1)
            | d.get("disposable_protein_positive", pd.Series(False, index=d.index))
        )
        & d.get("outside_disposable_floor_ratio", pd.Series(float("nan"), index=d.index)).lt(1)
        & ~d.get("protein_outside_disposable", pd.Series(True, index=d.index))
    )
    return d


def tumor_restricted_evidence(facts: pd.DataFrame) -> pd.DataFrame:
    """Separate non-reproductive tumor candidates from weak cancer RNA leads.

    Require >=10 native TPM in >=2 distinct tumor samples with resolved identity,
    or reviewed reliable tumor protein evidence. Restriction requires complete
    required normal coverage, all measured normal RNA below the 2-unit floor,
    no detected normal protein, and measured developmental expression below 2
    TPM. The floor includes permitted compartments; CP10K retains its separately
    specified nominal floor, not a TPM conversion. Final safety gates still apply.
    Missing sample counts never qualify. A reviewed normal-coverage assessment
    must address relevant rare cells and states before assigning restriction;
    numerical matches remain leads until then. This is panel-defined restriction,
    not a claim of absence in every normal cell or physiological state.
    """
    d = facts.copy()
    samples = pd.to_numeric(d.tumor_samples_ge10, errors="raise")
    if samples.dropna().lt(0).any() or samples.dropna().mod(1).ne(0).any():
        raise ValueError("Require nonnegative integer tumor sample counts or missing")
    for column in [
        "tumor_identity_resolved",
        "reliable_tumor_protein",
        "normal_restriction_reviewed",
    ]:
        if d[column].isna().any() or d[column].dtype != bool:
            raise ValueError(f"Explicit boolean required: {column}")
    d["tumor_restricted_rna_support"] = samples.ge(2) & d.tumor_identity_resolved
    d["tumor_restricted_normal_low"] = (
        d.required_coverage
        & ~d.unresolved_evidence
        & d.normal_maximum_floor_ratio.lt(1)
        & ~d.any_detected_protein
        & d.developmental_measured
        & d.developmental_all_lineage_max_tpm.lt(APPRECIABLE_RNA)
    )
    d["tumor_restricted_profile"] = (
        (d.tumor_restricted_rna_support | d.reliable_tumor_protein)
        & d.tumor_restricted_normal_low
        & ~d.positive_context
    )
    d["tumor_restricted_identity"] = d.tumor_restricted_profile & d.normal_restriction_reviewed
    d["tumor_restricted_lead"] = d.tumor_restricted_profile & ~d.normal_restriction_reviewed
    return d


def translation_admission(facts: pd.DataFrame) -> pd.Series:
    """Require both entry support and qualification of the exact selected record."""
    columns = ["translation_support_eligible", "translation_qualified"]
    for column in columns:
        if facts[column].isna().any() or facts[column].dtype != bool:
            raise ValueError("Complete boolean translation evidence required")
    return facts[columns].all(axis=1)


def classify(
    facts: pd.DataFrame,
    layers: Sequence[EvidenceLayer],
    *,
    entity_column: str,
    unsafe_overlap: Callable[[pd.Series], pd.Series] | None = None,
) -> FunnelResult:
    """Classify every entity with normal exclusions taking precedence.

    Every entity receives CTA, Unsafe or Other. CTA means every gate passes.
    Unsafe requires positive exclusion evidence, including accepted normal-cell
    findings and Manageable-tier excess above 128. Other is
    everything else; assessment_status records whether its normal screen was
    complete, incomplete or held, and the Held subtype marks proteins whose
    target gates pass but whose normal signal is unresolved or marginal.
    Cancer presence, specificity, developmental/reproductive association and
    coverage can block CTA admission but can never supply Unsafe evidence.

    ``unsafe_overlap`` maps the Unsafe mask to each entity's fraction of 9-mers
    that also occur in an Unsafe entity. A would-be CTA at or above
    ``CONTAINED_FLOOR`` becomes ``Other:Contained``: Unsafe sequence is
    blacklisted, so a protein almost wholly inside it has no usable peptide.
    """
    if facts[entity_column].isna().any() or not facts[entity_column].is_unique:
        raise ValueError("Require unique entity identities")
    if not layers or len({layer[0] for layer in layers}) != len(layers):
        raise ValueError("Require uniquely named evidence layers")
    for name, hold, unsafe in layers:
        for mask in (hold, unsafe):
            if not mask.index.equals(facts.index) or mask.isna().any() or mask.dtype != bool:
                raise ValueError("Complete aligned boolean evidence required")
        if name in CANDIDATE_ONLY_LAYERS and unsafe.any():
            raise ValueError(f"{name} is an admission requirement, not safety evidence")
    holds = pd.DataFrame({name: hold | unsafe for name, hold, unsafe in layers})
    exclusions = pd.DataFrame({name: unsafe for name, _, unsafe in layers})
    held, unsafe = holds.any(axis=1), exclusions.any(axis=1)
    # Other includes held and incomplete cases. Its first identity-failure
    # subtype must never substitute for the independent assessment status.
    concentration_hold = holds.get(
        "Reproductive concentration", pd.Series(False, index=facts.index)
    )
    safety_hold = holds[[n for n in holds if n not in CANDIDATE_ONLY_LAYERS]].any(axis=1)
    # Every protein ends as CTA, Unsafe or Other. Unsafe needs positive exclusion
    # evidence; unresolved or marginal normal evidence blocks CTA and is
    # disclosed in assessment_status, but is never a separate unassessed class.
    translation_hold = holds.get("Translation qualification", pd.Series(False, index=facts.index))
    complete = facts.required_coverage & ~facts.unresolved_evidence
    eligible = ~(held | unsafe) & complete & ~safety_hold
    overlap = pd.Series(np.nan, index=facts.index)
    contained = pd.Series(False, index=facts.index)
    if unsafe_overlap is not None:
        overlap = unsafe_overlap(unsafe)
        if not overlap.index.equals(facts.index):
            raise ValueError("Unsafe sequence overlap must align with the evidence index")
        contained = eligible & overlap.ge(CONTAINED_FLOOR)
        eligible = eligible & ~contained
    other = ~unsafe & ~eligible
    result = facts.drop(columns=[c for c in facts if c.startswith(("hold:", "unsafe:"))]).copy()
    result["category"] = pd.Series(pd.NA, index=facts.index, dtype="string")
    for label, mask in [
        ("CTA", eligible),
        ("Other", other),
        ("Unsafe", unsafe),
    ]:
        result.loc[mask, "category"] = label
    result["assessment_status"] = np.select(
        [
            unsafe | eligible | (other & complete & ~safety_hold & ~translation_hold),
            translation_hold & ~unsafe,
            facts.unresolved_evidence | ~facts.required_coverage,
        ],
        ["screened", "translation_review", "incomplete_evidence"],
        default="normal_evidence_hold",
    )
    association_failed = holds.get("Reproductive / embryonic association", ~facts.positive_context)
    cancer_failed = holds.get("Cancer expression", pd.Series(False, index=facts.index))
    result["unsafe_sequence_overlap_pct"] = (100 * overlap).round(2)
    nested = (
        facts["nested_microprotein"].astype(bool)
        if "nested_microprotein" in facts
        else pd.Series(False, index=facts.index)
    )
    result["subtype"] = np.select(
        [
            contained,
            other & nested,
            other & translation_hold,
            other & association_failed,
            other & cancer_failed,
            other & concentration_hold,
        ],
        ["Contained", "Nested", "Translation", "Nonreproductive", "Reproductive", "Diffuse"],
        default=np.where(other, "Held", ""),
    )
    result["classification"] = result.category.astype(str)
    result.loc[other, "classification"] = "Other:" + result.loc[other, "subtype"]
    if result.loc[other, "subtype"].eq("").any():
        raise ValueError("Every Other protein requires a categorical subtype")
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
                "Other": 0,
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
    if unsafe_overlap is not None:
        remaining = ~(running_hold | running_unsafe)
        sequential.append(
            {
                "layer": "Unsafe sequence containment",
                "new_candidate_removals": int(contained.sum()),
                "new_Unsafe": 0,
                "CTA_eligible_so_far": int((remaining & ~contained).sum()),
                "Unsafe": int(running_unsafe.sum()),
                "Other": 0,
                # Contained proteins are assigned Other at the final stage, like holds.
                "Unassigned": int((running_hold & ~running_unsafe).sum() + contained.sum()),
            }
        )
        independent.append(
            {
                "layer": "Unsafe sequence containment",
                "candidate_exclusions_anywhere": int((overlap.ge(CONTAINED_FLOOR) & ~unsafe).sum()),
                "Unsafe_evidence_anywhere": 0,
                "CTA_added_if_layer_removed": int(contained.sum()),
                "Unsafe_reclassified_if_layer_removed": 0,
            }
        )
    sequential.append(
        {
            "layer": "CTA and Other assignments",
            "new_candidate_removals": 0,
            "new_Unsafe": 0,
            "CTA_eligible_so_far": int(eligible.sum()),
            "Unsafe": int(unsafe.sum()),
            "Other": int(other.sum()),
            "Unassigned": int(result.category.isna().sum()),
        }
    )
    result["policy"] = POLICY_ID
    return FunnelResult(result, pd.DataFrame(sequential), pd.DataFrame(independent))

"""Current CTA selection rules, separate from frozen release policies.

RNA has one tissue/IHC ceiling. Crossing it supplies exclusion evidence in
Vital/Protected compartments; crossing it in Manageable compartments prevents
candidate and Neither admission and remains a normal-screen hold. Neither is
the adequately assessed normal-pass fallback. Unresolved results have a null category.
These are research selection rules, not estimates of clinical toxicity.

Normal-tissue safety, cancer expression and reproductive /
developmental association are distinct questions. Only safety evidence can
establish Unsafe. Cancer requires one qualifying measured sample, without a
tumor-versus-normal ratio. Association requires positive expression, without a
reproductive-versus-normal maximum veto. Profile diffuseness is assessed
separately from normal-tissue safety.
"""

from __future__ import annotations

from collections.abc import Collection, Sequence

import numpy as np
import pandas as pd

from .cell_rna_policy import EXPECTED_UNITS, GRADES, NEUTRAL_REPRODUCTIVE_IMMUNE, POSITIVE
from .compartment_rna_policy import IHC_RANKS
from .cta_funnel import EvidenceLayer, FunnelResult
from .expression_diffuseness import MAX_DIFFUSENESS

POLICY_ID = "oncoref.cta-selection.2026-10-30.1"
TIER_BASES = {"vital": 4.0, "protected": 8.0, "manageable": 16.0}
APPRECIABLE_RNA = 2.0
MIN_MANAGEABLE_TISSUES = 2
CANDIDATE_ONLY_LAYERS = frozenset(
    {
        "Reproductive / embryonic association",
        "Tumor sample ≥2 TPM/pTPM",
        "Required RNA coverage",
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
    reliable positive reproductive IHC. IHC must be mapped and support every
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
    d["reproductive_association"] = d.reproductive_rna_association | protein
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
    An RNA- or IHC-supported adult reference may be assessed in its matched RNA
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
    cells = d.cell_type.fillna("").str.lower()
    stroma = d.source.eq("hpa_single_cell_tissue") & cells.eq("fibroblasts")
    normal = d.category.isin([*TIER_BASES, "renewing_mucosa"]) & d.mapping_status.eq("mapped")
    d.loc[stroma & normal, "selection_tier"] = "manageable"
    d["compartment_basis"] = np.where(
        stroma & normal, "identified_fibroblast", "anatomical_default"
    )
    d.loc[d.selection_tier.ne(d.category) & ~stroma, "compartment_basis"] = "large_bowel_epithelium"
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
    d["unsafe_evidence"] = d.candidate_exclusion & d.selection_tier.isin(["vital", "protected"])
    d["unresolved_context"] = ~d.mapping_status.eq("mapped")
    d["candidate_exclusion"] |= d.unresolved_context
    d["selection_policy"] = POLICY_ID
    return d


def neither_profiles(facts: pd.DataFrame) -> pd.DataFrame:
    """Describe normal profiles; final Neither requires all normal gates to pass.

    Low observed normal expression uses complete required normal panels and no
    detected normal protein. Tumor abundance or missing developmental evidence
    cannot disqualify a normal-pass protein from Neither. The historical breadth
    description is descriptive only; it cannot override a Manageable ceiling.
    Optional DICE/cornea findings count when present; missing optional assays are
    disclosed, not treated as measurements of zero. Protein uncertainty blocks
    assignment until adjudicated. The final classifier gives Unsafe precedence.
    """
    d = facts.copy()
    complete = d.required_coverage & ~d.unresolved_evidence
    d["neither_low_expression"] = (
        complete & d.normal_maximum_floor_ratio.lt(1) & ~d.any_detected_protein
    )
    d["neither_diffuse_manageable"] = (
        complete
        & ~d.positive_context
        & ~d.vital_protected_exclusion
        & d.manageable_tissue_count.ge(MIN_MANAGEABLE_TISSUES)
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


def classify(
    facts: pd.DataFrame, layers: Sequence[EvidenceLayer], *, entity_column: str
) -> FunnelResult:
    """Classify every entity with normal exclusions taking precedence.

    Null category means no supported biological assignment; assessment_status
    explains incomplete normal evidence versus an unresolved normal-screen hold.
    Neither is the normal-pass noncandidate fallback; it is never a synonym for
    missing normal evidence. Tumor-restricted is an additional candidate subtype.
    Cancer presence, specificity, developmental/reproductive association and
    coverage can block CTA admission but can never supply Unsafe evidence.
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
    # A measured diffuse identity can establish Neither only when normal safety
    # has no unresolved hold. Tumor/developmental presence is not required for Neither.
    concentration_hold = holds.get(
        "Reproductive concentration", pd.Series(False, index=facts.index)
    )
    safety_hold = holds[[n for n in holds if n not in CANDIDATE_ONLY_LAYERS]].any(axis=1)
    diffuse = (
        concentration_hold & facts.required_coverage & ~facts.unresolved_evidence & ~safety_hold
    )
    normal_pass = facts.required_coverage & ~facts.unresolved_evidence & ~safety_hold & ~unsafe
    eligible = ~(held | unsafe) & normal_pass
    restricted = (
        facts.get("tumor_restricted_identity", pd.Series(False, index=facts.index))
        & normal_pass
        & ~eligible
    )
    neither = normal_pass & ~eligible & ~restricted
    result = facts.drop(columns=[c for c in facts if c.startswith(("hold:", "unsafe:"))]).copy()
    result["category"] = pd.Series(pd.NA, index=facts.index, dtype="string")
    for label, mask in [
        ("CTA", eligible),
        ("Tumor-restricted", restricted),
        ("Neither", neither),
        ("Unsafe", unsafe),
    ]:
        result.loc[mask, "category"] = label
    result["assessment_status"] = np.select(
        [
            result.category.notna(),
            facts.unresolved_evidence | ~facts.required_coverage,
        ],
        ["screened", "incomplete_evidence"],
        default="normal_evidence_hold",
    )
    result["neither_profile"] = np.select(
        [
            neither & diffuse,
            neither & facts.neither_low_expression,
            neither & facts.neither_diffuse_manageable,
            neither,
        ],
        [
            "diffuse_reproductive",
            "low_normal_expression",
            "manageable_expression",
            "normal_pass_noncandidate",
        ],
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
                "Tumor-restricted": 0,
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
            "layer": "Candidate subtypes and Neither fallback",
            "new_candidate_removals": 0,
            "new_Unsafe": 0,
            "CTA_eligible_so_far": int(eligible.sum()),
            "Unsafe": int(unsafe.sum()),
            "Neither": int(neither.sum()),
            "Tumor-restricted": int(restricted.sum()),
            "Unassigned": int(result.category.isna().sum()),
        }
    )
    result["policy"] = POLICY_ID
    return FunnelResult(result, pd.DataFrame(sequential), pd.DataFrame(independent))

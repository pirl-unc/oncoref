# Licensed under the Apache License, Version 2.0.
"""Disjoint CTA, Unsafe and Neither research sets for genes and exact proteins.

Only Unsafe supplies the default forbidden-protein background. Neither is not
a target set or evidence of absence. Historical expression screens remain dated
inputs; failing their selectivity or admission rules is not itself unsafe.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .cta_screen import BASELINES, RELIABLE
from .load_dataset import get_data
from .normal_tissues import SCREENED_NORMAL_TISSUE_TIERS
from .sequence_groups import protein_sequence_groups

POLICY_ID = "oncoref.target-sets.2026-10-23.1"
GROUPS = ("CTA", "Unsafe", "Neither")
MINIMUM_RNA_WITHOUT_POSITIVE_IHC = 5.0


def results():
    """Current genome-wide gene decisions, with CTA nomination and evidence scope."""
    return get_data("target-gene-sets").copy()


def gene_ids(group):
    """Return one disjoint group; Unsafe alone is the forbidden background."""
    if group not in GROUPS:
        raise ValueError(f"Unknown target group: {group}")
    frame = results()
    return set(frame.loc[frame.Category.eq(group), "Gene ID"])


def gene_names(group):
    """Known symbols only; gene_ids retains unnamed loci."""
    if group not in GROUPS:
        raise ValueError(f"Unknown target group: {group}")
    frame = results()
    return set(frame.loc[frame.Category.eq(group), "Gene"].dropna()) - {""}


def classify_genes(
    categories, screen, observations, *, safety_holds=(), unreliable_genes=(), donor_results=None
):
    """Classify the pinned human HPA v23 genome and ordered-screen candidates.

    ``categories`` is the checksum-pinned 21 October genome-wide expression
    result, with its original absolute-ceiling and protein exclusions. Those
    absolute ceilings are unchanged; the ordered policy only removes relaxation
    when no reliable reproductive positive control exists. Rechecking baseline
    ceilings in that case therefore exactly recovers the current absolute gate.
    ``observations`` contains all mapped RNA and IHC rows, including tissues
    without RNA. Missing IHC and measured negative IHC remain distinct.
    """
    if categories.Ensembl_Gene_ID.duplicated().any() or screen["Gene ID"].duplicated().any():
        raise ValueError("Gene inputs must have unique identities")
    if not categories.category_policy_id.eq(
        "oncoref.normal-expression-categories.2026-10-21.2"
    ).all():
        raise ValueError("Supply the pinned 21 October genome-wide expression policy")
    if not screen.Policy.eq("oncoref.cta-screen.2026-10-23.1").all():
        raise ValueError("Supply the pinned ordered CTA screen")
    if not observations.source_version.eq("v23").all():
        raise ValueError("Only human HPA v23 measurements support these RNA limits")
    for modality, source in [("bulk_rna", "hpa_rna_consensus"), ("ihc", "hpa_normal_tissue")]:
        if not observations.loc[observations.modality.eq(modality), "source_name"].eq(source).all():
            raise ValueError("Unexpected expression source")
    if not observations.modality.isin(["bulk_rna", "ihc"]).all():
        raise ValueError("Only mapped bulk RNA and IHC are supported")
    if "taxonomy_id" not in observations or not observations.taxonomy_id.eq(9606).all():
        raise ValueError("Do not transfer human HPA ceilings to another species")
    frame = categories.set_index("Ensembl_Gene_ID").copy()
    current = screen.set_index("Gene ID").reindex(frame.index)
    valid = observations.valid_measurement
    rna = observations.loc[observations.modality.eq("bulk_rna") & valid].copy()
    ihc = observations.loc[observations.modality.eq("ihc") & valid]
    positives = ihc.loc[ihc.strict_positive]
    mapped = positives.mapping_status.eq("mapped")
    controls = positives.loc[
        mapped
        & positives.tissue_category.eq("reproductive")
        & positives.ihc_reliability.isin(RELIABLE)
        & ~positives.Ensembl_Gene_ID.isin(unreliable_genes)
    ]
    control_ids = set(controls.Ensembl_Gene_ID)
    positive_ids = set(positives.Ensembl_Gene_ID)
    normal_protein = positives.loc[
        mapped
        & positives.screening_tier.isin(SCREENED_NORMAL_TISSUE_TIERS)
        & ~positives.Ensembl_Gene_ID.isin(unreliable_genes)
    ]
    rna["baseline"] = rna.screening_tier.map({tier: pair[0] for tier, pair in BASELINES.items()})
    baseline_high = rna.loc[
        ~rna.Ensembl_Gene_ID.isin(control_ids) & rna.rna_ntpm.gt(np.nextafter(rna.baseline, np.inf))
    ]
    # Baseline absolute counterevidence remains decisive. Selectivity is separate.
    rna_high = frame.rna_exclusions.ne("") | frame.index.isin(baseline_high.Ensembl_Gene_ID)
    protein_high = frame.index.isin(normal_protein.Ensembl_Gene_ID)
    held = frame.index.isin(safety_holds)
    variable = current.Category.eq("Variable") & ~protein_high & ~held
    donor_high_ids = set()
    if donor_results is not None and not donor_results.empty:
        high = donor_results["Highest RNA (nTPM)"].gt(
            np.nextafter(donor_results["Absolute ceiling (nTPM)"], np.inf)
        )
        donor_high_ids = set(donor_results.loc[high, "Gene ID"])
    donor_high = frame.index.isin(donor_high_ids)
    max_rna = rna.groupby("Ensembl_Gene_ID").rna_ntpm.max().reindex(frame.index)
    ihc_count = ihc.groupby("Ensembl_Gene_ID").size().reindex(frame.index, fill_value=0)
    positive_count = positives.groupby("Ensembl_Gene_ID").size().reindex(frame.index, fill_value=0)
    reproductive = pd.to_numeric(frame.reproductive_rna_max_ntpm, errors="coerce")
    usable = frame.quality_flags.eq("") & frame.rna_routine_missing.eq(0) & max_rna.notna()
    no_positive = ~frame.index.isin(positive_ids)
    weak = no_positive & reproductive.lt(MINIMUM_RNA_WITHOUT_POSITIVE_IHC)
    unsafe = protein_high | held | ((rna_high | donor_high) & ~variable)
    cta = (
        current.Category.eq("Retained") & ~weak & ~unsafe & ~variable & usable & reproductive.gt(0)
    )
    result = pd.DataFrame(index=frame.index)
    result["Gene"] = frame.Symbol
    result["Gene ID"] = frame.index
    result["Category"] = np.select([unsafe, cta], ["Unsafe", "CTA"], default="Neither")
    result["Detail"] = "Uncurated"
    result.loc[current.Category.eq("Excluded"), "Detail"] = "Screening"
    result.loc[frame.selectivity_exclusions.ne(""), "Detail"] = "Selectivity"
    result.loc[frame.other_subgroup.eq("Disposable"), "Detail"] = "Disposable"
    result.loc[weak & max_rna.lt(5), "Detail"] = "Low"
    result.loc[~usable | reproductive.isna(), "Detail"] = "Unassessed"
    result.loc[variable, "Detail"] = "Variable"
    result.loc[(rna_high | donor_high) & unsafe, "Detail"] = "RNA"
    result.loc[protein_high, "Detail"] = "Protein"
    result.loc[held, "Detail"] = "Safety"
    result.loc[cta, "Detail"] = "Pass"
    result["Target"] = result.Category.eq("CTA")
    result["Forbidden"] = result.Category.eq("Unsafe")
    result["Reproductive RNA (nTPM)"] = reproductive
    result["Maximum RNA (nTPM)"] = max_rna
    result["IHC measurements"] = ihc_count
    result["Positive IHC measurements"] = positive_count
    result["IHC coverage"] = np.select(
        [ihc_count.eq(0), positive_count.eq(0)],
        ["Missing", "Negative where measured"],
        default="Positive somewhere",
    )
    result["Normal protein positive"] = protein_high
    result["Normal RNA above ceiling"] = rna_high
    result["Normal donor RNA above ceiling"] = donor_high
    result["Reproductive positive control"] = frame.index.isin(control_ids)
    result["Previous CTA screen"] = current.Category.fillna("Not nominated")
    result["Previous expression group"] = frame.get("category", "")
    result["Disposable expression"] = frame.get("disposable_positive", False)
    result["Below-threshold subgroup"] = frame.other_subgroup
    result["Donor review"] = current["Donor review"].fillna("Not assessed")
    result["Donor findings"] = current["Donor findings"].fillna("")
    result["Clinical evidence"] = current["Clinical evidence"].fillna("")
    result["Previous reason"] = current.Reason.fillna("")
    explanations = {
        "Pass": "Passes the curated CTA screen and the 5 nTPM/no-positive-IHC admission rule",
        "Low": "No positive IHC anywhere measured; RNA below 5 nTPM in all measured tissues; normal RNA ceilings also pass",
        "Variable": "Donor-variable vital RNA with negative local IHC; excluded from targets and forbidden background pending review",
        "RNA": "Measured normal RNA exceeds its applicable category/IHC ceiling",
        "Protein": "Positive IHC in a Vital, Protected or Manageable compartment",
        "Safety": "Independent reviewed normal-tissue safety hold",
        "Selectivity": "Fails reproductive selectivity without above-ceiling normal RNA or a normal-protein veto",
        "Unassessed": "Evidence or anatomical coverage is incomplete/unusable; neither target nor default forbidden",
        "Disposable": "Below-threshold normal expression with Disposable expression; no curated CTA selection",
        "Uncurated": "No curated CTA selection and no decisive above-threshold normal-expression exclusion",
        "Screening": "Other CTA admission, identity or curation rule fails; not automatically forbidden",
    }
    result["Reason"] = result.Detail.map(explanations)
    additional = pd.Series(
        {
            gid: "; ".join(
                f"{r.tissue}: {r.rna_ntpm:g} nTPM (baseline {r.baseline:g})"
                for r in rows.itertuples()
            )
            for gid, rows in baseline_high.groupby("Ensembl_Gene_ID", sort=False)
        },
        dtype=str,
    )
    result["RNA evidence"] = frame.rna_exclusions + additional.reindex(frame.index).fillna("").map(
        lambda value: "; " + value if value else ""
    )
    result["IHC evidence"] = frame.ihc_exclusions
    result["Coverage gaps"] = frame.quality_flags
    result["Policy"] = POLICY_ID
    return result.reset_index(drop=True).sort_values(["Category", "Gene ID"])


def partition_proteins(occurrences, genes, *, reference, expression_reference_id):
    """Assign exact complete proteins once: Unsafe > Neither > CTA.

    Neither blocks target assignment when a full sequence is shared, but never
    supplies the forbidden background. Unclassified genes are Neither. Keep
    source occurrences; gene-level expression does not allocate individual isoforms.
    """
    if expression_reference_id != reference.reference_id:
        raise ValueError("Expression and protein reference identities must match")
    if genes["Gene ID"].duplicated().any() or not genes.Category.isin(GROUPS).all():
        raise ValueError("Gene sets must have unique IDs and recognized categories")
    indexed = genes.set_index("Gene ID").Category
    proteins = protein_sequence_groups(occurrences, reference=reference)
    for protein in proteins:
        members = sorted({row["gene_id"] for row in protein["occurrences"]})
        assignments = {gene: indexed.get(gene, "Neither") for gene in members}
        group = next(g for g in ("Unsafe", "Neither", "CTA") if g in assignments.values())
        protein.update(
            group=group,
            member_genes=members,
            member_groups=assignments,
            mixed_gene_groups=len(set(assignments.values())) > 1,
            target=group == "CTA",
            forbidden=group == "Unsafe",
        )
    return proteins


def repartition_protein_table(proteins, genes, *, reference_id):
    """Reassign a checksum-verified exact-sequence inventory without losing members."""
    if (
        proteins["Protein sequence"].duplicated().any()
        or not proteins.Reference.eq(reference_id).all()
    ):
        raise ValueError("Protein inventory must be unique and use the declared reference")
    if genes["Gene ID"].duplicated().any() or not genes.Category.isin(GROUPS).all():
        raise ValueError("Gene sets must have unique IDs and recognized categories")
    indexed = genes.set_index("Gene ID").Category
    result = proteins.copy()
    result["Previous group"] = result.Group
    groups, mixed = [], []
    for value in result.Genes:
        members = value.split(";")
        if not all(members):
            raise ValueError("A complete protein must retain its encoding gene identities")
        assignments = {indexed.get(gene, "Neither") for gene in members}
        groups.append(next(g for g in ("Unsafe", "Neither", "CTA") if g in assignments))
        mixed.append(len(assignments) > 1)
    result["Group"] = groups
    result["Mixed gene groups"] = mixed
    result["Target"] = result.Group.eq("CTA")
    result["Forbidden"] = result.Group.eq("Unsafe")
    result = result.rename(columns={"Exclude by default": "Previous exclusion"})
    return result

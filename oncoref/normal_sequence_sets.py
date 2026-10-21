# Licensed under the Apache License, Version 2.0.
"""Disjoint expression groups for genes and exact reference protein sequences."""

from __future__ import annotations

import pandas as pd

from .sequence_groups import protein_sequence_groups

GROUPS = ("Normal", "Unassessed", "Reproductive", "Other")
SET_POLICY_ID = "oncoref.normal-sequence-sets.v1"


def partition_genes(categories, *, normal_evidence_genes=()):
    """Preserve screening groups; explicit external normal evidence takes precedence.

    ``categories`` is the genome-wide expression report. Normal/Unassessed form
    the default exclusion background. Below-threshold observed normal expression
    remains a separate flag, available to a stricter downstream comparison.
    """
    if categories.Ensembl_Gene_ID.duplicated().any() or categories.Ensembl_Gene_ID.eq("").any():
        raise ValueError("Gene categories require unique nonempty identities")
    if not categories.category.isin(GROUPS).all():
        raise ValueError("Unknown expression group")
    result = (
        categories[
            [
                "Ensembl_Gene_ID",
                "Symbol",
                "category",
                "other_subgroup",
                "normal_expression_observed",
            ]
        ]
        .rename(
            columns={
                "Ensembl_Gene_ID": "gene_id",
                "Symbol": "symbol",
                "category": "expression_group",
                "normal_expression_observed": "observed_normal_expression",
            }
        )
        .copy()
    )
    result["group"] = result.expression_group
    result["reason"] = "Expression screening rules"
    held = result.gene_id.isin(normal_evidence_genes)
    result.loc[held, "group"] = "Normal"
    result.loc[held, "reason"] = "External normal-tissue evidence"
    result["exclude_by_default"] = result.group.isin(["Normal", "Unassessed"])
    return result.sort_values("gene_id").reset_index(drop=True)


def partition_proteins(occurrences, genes, *, reference, expression_reference_id):
    """Assign each complete exact sequence once, retaining every source occurrence.

    All annotated complete translations are considered, not only each gene's
    longest isoform. Normal takes precedence over Unassessed, then Reproductive,
    then Other. An absent gene classification is Unassessed. Same-sequence
    occurrences in different genes therefore cannot enter competing sets.

    This inherits gene-level expression; it does not infer isoform-specific RNA
    allocation, translation, PTMs or k-mer specificity. Different full proteins
    can still share k-mers. Incomplete/invalid translations must be retained by
    the caller for review, with their gene's expression group; they cannot be
    supplied here as complete proteins.
    """
    if expression_reference_id != reference.reference_id:
        raise ValueError("Expression and sequence reference identities must match")
    if genes.gene_id.duplicated().any() or not genes.group.isin(GROUPS).all():
        raise ValueError("Gene groups must be unique and recognized")
    indexed = genes.set_index("gene_id")
    result = protein_sequence_groups(occurrences, reference=reference)
    for protein in result:
        members = sorted({row["gene_id"] for row in protein["occurrences"]})
        assignments = {
            gene: indexed.loc[gene, "group"] if gene in indexed.index else "Unassessed"
            for gene in members
        }
        group = next(group for group in GROUPS if group in assignments.values())
        protein.update(
            group=group,
            member_genes=members,
            member_groups=assignments,
            mixed_gene_groups=len(set(assignments.values())) > 1,
            exclude_by_default=group in {"Normal", "Unassessed"},
            observed_normal_expression=any(
                bool(indexed.loc[gene, "observed_normal_expression"])
                for gene in members
                if gene in indexed.index
            ),
        )
    return result


def gene_group_counts(genes):
    """Readable count table in the same order used for sequence conflict resolution."""
    return pd.DataFrame({"Group": GROUPS, "Genes": [int(genes.group.eq(g).sum()) for g in GROUPS]})

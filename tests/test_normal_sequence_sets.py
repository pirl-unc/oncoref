"""Shared isoforms must not leak from normal or unresolved genes into targets."""

import pandas as pd
import pytest

from oncoref.normal_sequence_sets import partition_genes, partition_proteins
from oncoref.reference_identity import ReferenceIdentity


@pytest.fixture
def reference():
    return ReferenceIdentity(
        taxon_id=9606,
        assembly_accession="fixture.1",
        annotation_source="fixture",
        annotation_release="1",
        source_version="1",
        identifier_namespaces={"gene": "fixture", "transcript": "fixture", "protein": "fixture"},
        content_provenance={"proteins": "a" * 64},
    )


def genes():
    return partition_genes(
        pd.DataFrame(
            [
                {
                    "Ensembl_Gene_ID": gene,
                    "Symbol": gene.upper(),
                    "category": group,
                    "other_subgroup": "",
                    "normal_expression_observed": False,
                }
                for gene, group in [
                    ("r", "Reproductive"),
                    ("n", "Normal"),
                    ("u", "Unassessed"),
                    ("o", "Other"),
                ]
            ]
        )
    )


def occurrence(reference, gene, sequence, suffix=""):
    return {
        "reference_id": reference.reference_id,
        "occurrence_id": gene + suffix,
        "gene_id": gene,
        "transcript_id": gene + "tx" + suffix,
        "protein_id": gene + "p" + suffix,
        "sequence": sequence,
        "complete": True,
        "source": "fixture",
    }


def partition(reference, records):
    return partition_proteins(
        records, genes(), reference=reference, expression_reference_id=reference.reference_id
    )


def test_all_isoforms_and_normal_shared_sequences_are_kept(reference):
    rows = [
        occurrence(reference, "r", "MACDEFGHIK"),
        occurrence(reference, "r", "MQRSTVWY", "2"),
        occurrence(reference, "n", "MQRSTVWY"),
    ]
    groups = partition(reference, rows)
    by_sequence = {p["sequence"]: p for p in groups}
    assert by_sequence["MACDEFGHIK"]["group"] == "Reproductive"
    assert by_sequence["MQRSTVWY"]["group"] == "Normal"
    assert by_sequence["MQRSTVWY"]["mixed_gene_groups"]
    assert len(by_sequence["MQRSTVWY"]["occurrences"]) == 2


@pytest.mark.parametrize(
    "other,expected", [("u", "Unassessed"), ("missing", "Unassessed"), ("o", "Reproductive")]
)
def test_missing_evidence_blocks_shared_target_assignment(reference, other, expected):
    groups = partition(
        reference, [occurrence(reference, "r", "MALI"), occurrence(reference, other, "MALI")]
    )
    assert len(groups) == 1 and groups[0]["group"] == expected


def test_disjoint_proteins_can_still_share_kmers(reference):
    groups = partition(
        reference,
        [occurrence(reference, "r", "MACDEFGHIK"), occurrence(reference, "n", "TACDEFGHIK")],
    )
    assert {p["group"] for p in groups} == {"Reproductive", "Normal"}
    assert len({p["protein_sequence_id"] for p in groups}) == 2
    assert {p["sequence"][1:] for p in groups} == {"ACDEFGHIK"}


def test_reference_mismatch_and_incomplete_translation_fail(reference):
    row = occurrence(reference, "r", "MALI")
    with pytest.raises(ValueError, match="reference identities"):
        partition_proteins([row], genes(), reference=reference, expression_reference_id="other")
    with pytest.raises(ValueError, match="incomplete translation"):
        partition(reference, [{**row, "complete": False}])


def test_external_normal_evidence_is_visible_without_rewriting_expression_group():
    categories = pd.DataFrame(
        [
            {
                "Ensembl_Gene_ID": "MAGEA12",
                "Symbol": "MAGEA12",
                "category": "Reproductive",
                "other_subgroup": "",
                "normal_expression_observed": False,
            }
        ]
    )
    row = partition_genes(categories, normal_evidence_genes=["MAGEA12"]).iloc[0]
    assert row.expression_group == "Reproductive" and row.group == "Normal"
    assert row.exclude_by_default and row.reason == "External normal-tissue evidence"
    assert categories.iloc[0].category == "Reproductive"

"""Species/reference identity, exact sequences and lossless occurrence provenance."""

import copy
import json

import pytest

from oncoref import protein_sequence_groups


@pytest.fixture
def scope():
    return {
        "taxon_id": 9615,
        "assembly_accession": "GCF_fixture.1",
        "annotation_release": "fixture-106",
        "content_provenance": {"annotation.gtf": "a" * 64, "proteins.fa": "b" * 64},
    }


@pytest.fixture
def occurrence(scope):
    return {
        **scope,
        "occurrence_id": "a",
        "gene_id": "g1.1",
        "transcript_id": "t1.2",
        "protein_id": "p1.3",
        "source": "synthetic",
        "complete": True,
        "sequence": " maLi *\n",
        "coordinates": {"contig": "1", "start": 10, "end": 22, "strand": "+"},
        "evidence": [{"source_record": "CDS1", "validated": True}],
    }


def test_exact_sequences_retain_all_loci_and_transcripts(scope, occurrence):
    second = dict(occurrence, occurrence_id="b", gene_id="g2", transcript_id="t2")
    third = dict(occurrence, occurrence_id="c", transcript_id="t1.alternative")
    distinct = dict(occurrence, occurrence_id="d", sequence="MALL")
    rows = [occurrence, second, third, distinct, occurrence]
    original = copy.deepcopy(rows)
    groups = protein_sequence_groups(rows, **scope)
    assert rows == original
    assert {g["sequence"] for g in groups} == {"MALI", "MALL"}
    shared = next(g for g in groups if g["sequence"] == "MALI")
    assert shared["occurrences"] == rows[:3]
    assert shared["occurrences"][0]["sequence"] == " maLi *\n"
    assert shared["occurrences"][0]["protein_id"] == "p1.3"
    shared["occurrences"][0]["evidence"][0]["validated"] = False
    assert rows == original


def test_order_independence_and_json_round_trip(scope, occurrence):
    rows = [
        dict(occurrence, occurrence_id="z"),
        occurrence,
        dict(occurrence, occurrence_id="b", sequence="MALL"),
    ]
    groups = protein_sequence_groups(rows, **scope)
    serialized = json.dumps(groups, sort_keys=True, allow_nan=False)
    assert json.loads(serialized) == groups
    assert (
        json.dumps(protein_sequence_groups(reversed(rows), **scope), sort_keys=True) == serialized
    )


@pytest.mark.parametrize("taxon_id", [9606, 10090])
def test_human_mouse_are_explicit_and_reference_scoped(scope, occurrence, taxon_id):
    dog = protein_sequence_groups([occurrence], **scope)[0]
    other_scope = dict(scope, taxon_id=taxon_id)
    other = protein_sequence_groups([dict(occurrence, taxon_id=taxon_id)], **other_scope)[0]
    assert other["protein_sequence_id"] == dog["protein_sequence_id"]
    assert other["group_id"] != dog["group_id"]
    assert other["reference_id"] != dog["reference_id"]
    assert other["taxon_id"] == taxon_id


@pytest.mark.parametrize(
    "field,value",
    [
        ("gene_id", "another_gene"),
        ("transcript_id", "another_tx"),
        ("protein_id", "other"),
        ("sequence", "MALL"),
        ("source", "different-source"),
        ("coordinates", {"contig": "2"}),
    ],
)
def test_contradictory_occurrence_ids_are_rejected(scope, occurrence, field, value):
    with pytest.raises(ValueError, match="Occurrence 'a': contradictory reuse"):
        protein_sequence_groups([occurrence, dict(occurrence, **{field: value})], **scope)


@pytest.mark.parametrize("left,right", [(True, 1), (1, 1.0), (0.0, -0.0)])
def test_duplicate_provenance_cannot_lose_json_types(scope, occurrence, left, right):
    rows = [dict(occurrence, evidence={"value": value}) for value in (left, right)]
    for ordered in (rows, list(reversed(rows))):
        with pytest.raises(ValueError, match="Occurrence 'a': contradictory reuse"):
            protein_sequence_groups(ordered, **scope)


@pytest.mark.parametrize(
    "change",
    [
        {"taxon_id": 9606},
        {"taxon_id": 9615.0},
        {"assembly_accession": "GCF_fixture.2"},
        {"annotation_release": "107"},
        {"content_provenance": {"proteins.fa": "c" * 64}},
    ],
)
def test_mixed_references_are_rejected(scope, occurrence, change):
    with pytest.raises(ValueError, match="Occurrence 'a': conflicting"):
        protein_sequence_groups([occurrence, dict(occurrence, **change)], **scope)


@pytest.mark.parametrize(
    "sequence", [None, 123, "", "*", "MAX", "MAB", "MAZ", "MAJ", "MAU", "MAO", "MA*I", "MALI**"]
)
def test_invalid_translations_have_recorded_reason(scope, occurrence, sequence):
    with pytest.raises(ValueError, match=r"Occurrence 'a': .*sequence"):
        protein_sequence_groups([dict(occurrence, sequence=sequence)], **scope)


@pytest.mark.parametrize("complete", [None, False, 1, "true"])
def test_completeness_is_explicit(scope, occurrence, complete):
    with pytest.raises(ValueError, match="Occurrence 'a': incomplete translation"):
        protein_sequence_groups([dict(occurrence, complete=complete)], **scope)


@pytest.mark.parametrize(
    "field", ["occurrence_id", "gene_id", "transcript_id", "protein_id", "source", "sequence"]
)
def test_missing_source_fields_are_rejected(scope, occurrence, field):
    del occurrence[field]
    with pytest.raises(ValueError):
        protein_sequence_groups([occurrence], **scope)


@pytest.mark.parametrize(
    "change",
    [
        {"taxon_id": True},
        {"taxon_id": 0},
        {"taxon_id": "9615"},
        {"assembly_accession": " "},
        {"annotation_release": None},
        {"content_provenance": {}},
        {"content_provenance": {"file": "not-a-digest"}},
    ],
)
def test_explicit_reference_contract(scope, change):
    with pytest.raises(ValueError):
        protein_sequence_groups([], **dict(scope, **change))


def test_content_changes_reference_identity(scope, occurrence):
    original = protein_sequence_groups([occurrence], **scope)[0]
    updated = dict(scope, content_provenance={"annotation.gtf": "c" * 64, "proteins.fa": "b" * 64})
    new = protein_sequence_groups(
        [dict(occurrence, content_provenance=updated["content_provenance"])], **updated
    )[0]
    assert original["protein_sequence_id"] == new["protein_sequence_id"]
    assert original["group_id"] != new["group_id"]


def test_digest_case_does_not_change_reference_identity(scope, occurrence):
    original = protein_sequence_groups([occurrence], **scope)[0]
    uppercase = {key: digest.upper() for key, digest in scope["content_provenance"].items()}
    new = protein_sequence_groups(
        [dict(occurrence, content_provenance=uppercase)],
        **dict(scope, content_provenance=uppercase),
    )[0]
    assert new["group_id"] == original["group_id"]
    assert new["occurrences"][0]["content_provenance"] == uppercase


@pytest.mark.parametrize("value", [float("nan"), float("inf"), {"x"}, (1, 2), {1: "bad-key"}])
def test_source_fields_must_round_trip_losslessly(scope, occurrence, value):
    with pytest.raises(ValueError, match="finite JSON-native"):
        protein_sequence_groups([dict(occurrence, extra=value)], **scope)

"""Species/reference identity, exact sequences and lossless occurrence provenance."""

import copy
import json

import pytest

from oncoref import (
    ReferenceIdentity,
    normalize_protein_sequence,
    protein_sequence_groups,
    protein_sequence_id,
)


def group(rows, **scope):
    return protein_sequence_groups(rows, reference=ReferenceIdentity(**scope))


@pytest.fixture
def scope():
    return {
        "taxon_id": 9615,
        "assembly_accession": "GCF_fixture.1",
        "annotation_release": "fixture-106",
        "annotation_source": "example.org",
        "source_version": "snapshot-1",
        "identifier_namespaces": {
            "gene": "example:gene",
            "transcript": "example:transcript",
            "protein": "example:protein",
        },
        "content_provenance": {"annotation.gtf": "a" * 64, "proteins.fa": "b" * 64},
    }


@pytest.fixture
def occurrence(scope):
    return {
        **scope,
        "reference_id": ReferenceIdentity(**scope).reference_id,
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
    groups = group(rows, **scope)
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
    groups = group(rows, **scope)
    serialized = json.dumps(groups, sort_keys=True, allow_nan=False)
    assert json.loads(serialized) == groups
    assert json.dumps(group(reversed(rows), **scope), sort_keys=True) == serialized


@pytest.mark.parametrize("taxon_id", [9606, 10090])
def test_human_mouse_are_explicit_and_reference_scoped(scope, occurrence, taxon_id):
    dog = group([occurrence], **scope)[0]
    other_scope = dict(scope, taxon_id=taxon_id)
    other = group(
        [
            dict(
                occurrence,
                taxon_id=taxon_id,
                reference_id=ReferenceIdentity(**other_scope).reference_id,
            )
        ],
        **other_scope,
    )[0]
    assert other["protein_sequence_id"] == dog["protein_sequence_id"]
    assert other["group_id"] != dog["group_id"]
    assert other["reference_id"] != dog["reference_id"]
    assert other["reference"]["taxon_id"] == taxon_id


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
        group([occurrence, dict(occurrence, **{field: value})], **scope)


@pytest.mark.parametrize("left,right", [(True, 1), (1, 1.0), (0.0, -0.0)])
def test_duplicate_provenance_cannot_lose_json_types(scope, occurrence, left, right):
    rows = [dict(occurrence, evidence={"value": value}) for value in (left, right)]
    for ordered in (rows, list(reversed(rows))):
        with pytest.raises(ValueError, match="Occurrence 'a': contradictory reuse"):
            group(ordered, **scope)


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
        group([occurrence, dict(occurrence, **change)], **scope)


@pytest.mark.parametrize(
    "sequence", [None, 123, "", "*", "MAX", "MAB", "MAZ", "MAJ", "MAU", "MAO", "MA*I", "MALI**"]
)
def test_invalid_translations_have_recorded_reason(scope, occurrence, sequence):
    with pytest.raises(ValueError, match=r"Occurrence 'a': .*sequence"):
        group([dict(occurrence, sequence=sequence)], **scope)


@pytest.mark.parametrize("complete", [None, False, 1, "true"])
def test_completeness_is_explicit(scope, occurrence, complete):
    with pytest.raises(ValueError, match="Occurrence 'a': incomplete translation"):
        group([dict(occurrence, complete=complete)], **scope)


@pytest.mark.parametrize(
    "field", ["occurrence_id", "gene_id", "transcript_id", "protein_id", "source", "sequence"]
)
def test_missing_source_fields_are_rejected(scope, occurrence, field):
    del occurrence[field]
    with pytest.raises(ValueError):
        group([occurrence], **scope)


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
        group([], **dict(scope, **change))


def test_content_changes_reference_identity(scope, occurrence):
    original = group([occurrence], **scope)[0]
    updated = dict(scope, content_provenance={"annotation.gtf": "c" * 64, "proteins.fa": "b" * 64})
    new = group(
        [
            dict(
                occurrence,
                content_provenance=updated["content_provenance"],
                reference_id=ReferenceIdentity(**updated).reference_id,
            )
        ],
        **updated,
    )[0]
    assert original["protein_sequence_id"] == new["protein_sequence_id"]
    assert original["group_id"] != new["group_id"]


def test_digest_case_does_not_change_reference_identity(scope, occurrence):
    original = group([occurrence], **scope)[0]
    uppercase = {key: digest.upper() for key, digest in scope["content_provenance"].items()}
    new = group(
        [dict(occurrence, content_provenance=uppercase)],
        **dict(scope, content_provenance=uppercase),
    )[0]
    assert new["group_id"] == original["group_id"]
    assert new["occurrences"][0]["content_provenance"] == uppercase


@pytest.mark.parametrize("value", [float("nan"), float("inf"), {"x"}, (1, 2), {1: "bad-key"}])
def test_source_fields_must_round_trip_losslessly(scope, occurrence, value):
    with pytest.raises(ValueError, match="finite JSON-native"):
        group([dict(occurrence, extra=value)], **scope)


@pytest.mark.parametrize(
    "sequence", ["Mß", "Mı", "Mﬀ", "Mſ", "ＭＡ", "ΜA", "M\u00a0A", "M\u2003A", "M\u200bA"]
)
def test_unicode_sequences_rejected_before_normalization(scope, occurrence, sequence):
    for make in (normalize_protein_sequence, protein_sequence_id):
        with pytest.raises(ValueError, match="invalid protein sequence characters"):
            make(sequence)
    with pytest.raises(ValueError, match="Occurrence 'a': invalid protein sequence characters"):
        group([dict(occurrence, sequence=sequence)], **scope)


def test_every_non_ascii_character_is_rejected_before_case_conversion():
    # Exhaustive Unicode input boundary: expanding case mappings, homoglyphs,
    # surrogates and Unicode whitespace cannot become accepted AA20 sequence.
    for code in range(128, 0x110000):
        with pytest.raises(ValueError):
            normalize_protein_sequence("M" + chr(code))


def test_ascii_normalization_and_protein_id_contract():
    assert normalize_protein_sequence(" maLi \t\r\n\v\f*") == "MALI"
    assert protein_sequence_id("mali*") == protein_sequence_id("MALI")
    assert protein_sequence_id("MALI") != protein_sequence_id("MALL")
    assert protein_sequence_id("MALI").startswith("openvax:protein-sequence:v1:sha256:")


def test_unversioned_or_missing_reference_id_rejected(scope, occurrence):
    for value in (None, "a" * 64):
        with pytest.raises(ValueError, match="reference_id"):
            group([dict(occurrence, reference_id=value)], **scope)


def test_cyclic_source_provenance_has_clear_error(scope, occurrence):
    occurrence["cycle"] = occurrence
    with pytest.raises(ValueError, match="acyclic"):
        group([occurrence], **scope)


def test_required_identifier_namespaces(scope):
    scope["identifier_namespaces"] = {"gene": "example:gene"}
    with pytest.raises(ValueError, match="identifier namespaces"):
        group([], **scope)

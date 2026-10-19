"""Portable reference schema and explicit consumer migration contracts."""

import copy
import json
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

from oncoref import ReferenceIdentity, adapt_reference_occurrences, protein_sequence_groups


@pytest.fixture
def fields():
    return {
        "taxon_id": 9606,
        "assembly_accession": "GCF_000001405.40",
        "annotation_source": "NCBI RefSeq",
        "annotation_release": "110",
        "source_version": "2022-04",
        "identifier_namespaces": {"gene": "NCBI Gene", "transcript": "RefSeq", "protein": "RefSeq"},
        "content_provenance": {"annotation": "a" * 64, "protein": "b" * 64},
    }


def test_reference_round_trip_immutable_and_detached(fields):
    reference = ReferenceIdentity(**fields)
    before = reference.as_dict()
    assert hash(reference) == hash(ReferenceIdentity.from_dict(before))
    assert ReferenceIdentity.from_dict(json.loads(json.dumps(before))) == reference
    fields["content_provenance"]["protein"] = "c" * 64
    fields["identifier_namespaces"]["protein"] = "another"
    saved = reference.as_dict()
    saved["content_provenance"]["protein"] = "d" * 64
    assert reference.as_dict() == before
    with pytest.raises(FrozenInstanceError):
        reference.taxon_id = 9615
    with pytest.raises(TypeError):
        reference.content_provenance["protein"] = "c" * 64


@pytest.mark.parametrize(
    "field,value",
    [
        ("taxon_id", 9615),
        ("assembly_accession", "GCF_000001405.39"),
        ("annotation_source", "Ensembl"),
        ("annotation_release", "111"),
        ("source_version", "2023-04"),
        ("identifier_namespaces", {"gene": "Ensembl"}),
        ("content_provenance", {"protein": "c" * 64}),
    ],
)
def test_every_declared_scope_component_changes_identity(fields, field, value):
    assert (
        ReferenceIdentity(**fields).reference_id
        != ReferenceIdentity(**{**fields, field: value}).reference_id
    )


def test_asset_order_and_hash_case_do_not_change_identity(fields):
    original = ReferenceIdentity(**fields)
    fields["content_provenance"] = {"protein": "B" * 64, "annotation": "A" * 64}
    fields["identifier_namespaces"] = dict(reversed(list(fields["identifier_namespaces"].items())))
    assert ReferenceIdentity(**fields).reference_id == original.reference_id


@pytest.mark.parametrize(
    "field,value",
    [
        ("taxon_id", True),
        ("taxon_id", 0),
        ("taxon_id", 9606.0),
        ("assembly_accession", ""),
        ("annotation_source", " unknown "),
        ("annotation_release", None),
        ("source_version", "v1\n"),
        ("identifier_namespaces", {}),
        ("identifier_namespaces", {"gene": 123}),
        ("content_provenance", {}),
        ("content_provenance", {"protein": "g" * 64}),
        ("content_provenance", {"protein": "a" * 63}),
    ],
)
def test_invalid_scope_declarations_fail(fields, field, value):
    with pytest.raises(ValueError):
        ReferenceIdentity(**{**fields, field: value})


@pytest.mark.parametrize(
    "change",
    [
        {"schema_version": 2},
        {"schema_version": True},
        {"schema": "consumer.reference"},
        {"reference_id": "a" * 64},
        {"taxon_id": 9615},
        {"unexpected": "metadata"},
    ],
)
def test_reader_rejects_tampering_and_unrecognized_contract(fields, change):
    record = ReferenceIdentity(**fields).as_dict()
    with pytest.raises(ValueError):
        ReferenceIdentity.from_dict({**record, **change})


def test_reader_requires_stored_id_and_schema(fields):
    record = ReferenceIdentity(**fields).as_dict()
    for missing in record:
        with pytest.raises(ValueError):
            ReferenceIdentity.from_dict({k: v for k, v in record.items() if k != missing})
    with pytest.raises(ValueError):
        ReferenceIdentity.from_dict(fields)


@pytest.fixture
def migration(fields):
    source = {
        "species": fields["taxon_id"],
        "assembly": fields["assembly_accession"],
        "release": fields["annotation_release"],
        "files": fields["content_provenance"],
        "upstream_version": fields["source_version"],
        "consumer_note": "retain this",
    }
    row = {
        "occurrence_id": "one",
        "gene_id": "gene.1",
        "transcript_id": "tx.2",
        "protein_id": "protein.3",
        "source": "example",
        "sequence": "mali*",
        "complete": True,
        "old_reference": "opaque:ref-1",
        "coordinates": {"contig": "1", "start": 2, "end": 14},
        "evidence": [{"caller": "example"}],
    }
    options = {
        "reference": ReferenceIdentity(**fields),
        "source_reference": source,
        "source_reference_id": "opaque:ref-1",
        "reference_id_field": "old_reference",
        "field_map": {
            "taxon_id": "species",
            "assembly_accession": "assembly",
            "annotation_release": "release",
            "content_provenance": "files",
            "source_version": "upstream_version",
        },
        "added_fields": {
            "annotation_source": fields["annotation_source"],
            "identifier_namespaces": fields["identifier_namespaces"],
        },
        "evidence": {"reviewed_source_manifest": "d" * 64},
    }
    return row, options


@pytest.mark.parametrize("taxon_id", [9606, 9615, 10090])
def test_explicit_adapter_is_species_and_consumer_independent(migration, fields, taxon_id):
    row, options = migration
    fields["taxon_id"] = taxon_id
    options["reference"] = ReferenceIdentity(**fields)
    options["source_reference"]["species"] = taxon_id
    original = copy.deepcopy(row)
    result = adapt_reference_occurrences([row], **options)
    assert result["source_reference"] == options["source_reference"]
    assert result["occurrences"][0]["source_record"] == original
    assert result["occurrences"][0]["old_reference"] == "opaque:ref-1"
    assert result["occurrences"][0]["reference_migration_id"] == result["migration_id"]
    assert json.loads(json.dumps(result)) == result
    groups = protein_sequence_groups(result["occurrences"], reference=options["reference"])
    assert groups[0]["occurrences"][0]["source_record"] == original
    assert groups[0]["sequence"] == "MALI"
    assert row == original
    result["occurrences"][0]["source_record"]["evidence"][0]["caller"] = "changed"
    assert row == original


@pytest.mark.parametrize(
    "field,value",
    [
        ("species", 9615),
        ("assembly", "other.1"),
        ("release", "other"),
        ("files", {"protein": "f" * 64}),
    ],
)
def test_adapter_rejects_cross_reference_relabeling(migration, field, value):
    row, options = migration
    options["source_reference"][field] = value
    with pytest.raises(ValueError, match="does not match"):
        adapt_reference_occurrences([row], **options)


def test_adapter_requires_scope_and_content_to_come_from_source(migration):
    row, options = migration
    for field in ("taxon_id", "assembly_accession", "annotation_release", "content_provenance"):
        modified = {
            **options,
            "field_map": dict(options["field_map"]),
            "added_fields": dict(options["added_fields"]),
        }
        old = modified["field_map"].pop(field)
        modified["added_fields"][field] = options["source_reference"][old]
        with pytest.raises(ValueError, match="Map source scope/content"):
            adapt_reference_occurrences([row], **modified)


@pytest.mark.parametrize(
    "change",
    [
        {"old_reference": "wrong"},
        {"reference_id": "existing"},
        {"source_record": {}},
        {"reference_migration_id": "existing"},
    ],
)
def test_adapter_rejects_mixed_or_already_migrated_rows(migration, change):
    row, options = migration
    with pytest.raises(ValueError):
        adapt_reference_occurrences([{**row, **change}], **options)


def test_adapter_preserves_source_id_when_field_name_already_matches(migration):
    row, options = migration
    row["reference_id"] = row.pop("old_reference")
    options["reference_id_field"] = "reference_id"
    result = adapt_reference_occurrences([row], **options)
    assert result["occurrences"][0]["source_record"]["reference_id"] == "opaque:ref-1"
    assert result["occurrences"][0]["reference_id"] == options["reference"].reference_id


def test_adapter_requires_explicit_evidence(migration):
    row, options = migration
    with pytest.raises(ValueError, match="migration evidence"):
        adapt_reference_occurrences([row], **{**options, "evidence": {}})


def test_portable_identity_golden_vectors():
    from oncoref import protein_sequence_id

    fixture = Path(__file__).parent / "fixtures" / "reference-identity-v1.json"
    for vector in json.loads(fixture.read_text()):
        reference = ReferenceIdentity.from_dict(vector["reference"])
        assert protein_sequence_id(vector["sequence"]) == vector["protein_sequence_id"]
        row = {
            "occurrence_id": "a",
            "gene_id": "g",
            "transcript_id": "t",
            "protein_id": "p",
            "source": "fixture",
            "sequence": vector["sequence"],
            "complete": True,
            "reference_id": reference.reference_id,
        }
        group = protein_sequence_groups([row], reference=reference)[0]
        assert group["group_id"] == vector["group_id"]


@pytest.mark.parametrize("field", ["taxon_id", "species"])
def test_adapter_rejects_conflicting_per_occurrence_scope(migration, field):
    row, options = migration
    row[field] = 9615
    with pytest.raises(ValueError, match=r"Occurrence.*conflicting reference field"):
        adapt_reference_occurrences([row], **options)


@pytest.mark.parametrize("key", ["/tmp/protein.fa", "https://example.org/file", r"C:\file.fa"])
def test_local_paths_and_urls_are_not_identity_asset_roles(fields, key):
    fields["content_provenance"] = {key: "a" * 64}
    with pytest.raises(ValueError, match="logical asset roles"):
        ReferenceIdentity(**fields)

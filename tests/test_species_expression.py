# Licensed under the Apache License, Version 2.0.
import json
from dataclasses import replace

import pytest

import oncoref
from oncoref import catalog, cli
from oncoref.expression_contract import (
    ExpressionContract,
    ExpressionMeasurement,
    MeasurementTransform,
    ReferenceCrosswalk,
    validate_expression_join,
)
from oncoref.expression_registry import ExpressionSource, expression_sources
from oncoref.reference_identity import ReferenceIdentity


@pytest.fixture
def contract():
    reference = ReferenceIdentity(
        taxon_id=9615,
        assembly_accession="GCF_fixture.1",
        annotation_source="NCBI RefSeq",
        annotation_release="106",
        source_version="fixture-v1",
        identifier_namespaces={"gene": "NCBI Gene"},
        content_provenance={"gtf": "a" * 64},
    )
    return ExpressionContract(
        taxon_id=9615,
        source_id="fixture",
        source_version="v1",
        context="normal",
        reference=reference,
        reference_status="verified",
        feature_namespace="NCBI Gene",
        measurement=ExpressionMeasurement("gene", "TPM", "linear", False),
        content_provenance={"matrix": "b" * 64},
        provenance_notes="Synthetic test fixture",
    )


@pytest.mark.parametrize("value", ["dog", "canine", "Canis lupus familiaris", 9615, "9615"])
def test_species_boundary(value):
    dog = oncoref.for_species(value)
    assert dog.species.taxon_id == 9615
    assert {g.name for g in dog.groups} >= {"references", "expression", "cancer", "antigens"}
    assert dog.cancer.datasets() == ()
    assert dog.antigens.coverage()["coverage"] == "partial"
    assert all(s.taxon_id == 9615 for s in dog.expression.sources())
    assert all(r["taxon_id"] == 9615 for r in catalog.inventory(species=value))


@pytest.mark.parametrize("value", [None, True, 9615.0, "wolf", "", "all"])
def test_unknown_species_cannot_fall_back(value):
    with pytest.raises(ValueError, match="Unknown species"):
        oncoref.for_species(value)


def test_human_compatibility_and_group_partition():
    human = oncoref.for_species("human")
    assert human.expression.sources() == expression_sources()
    rows = human.inventory()
    assert rows and all(r["taxon_id"] == 9606 for r in rows)
    names = [r["name"] for g in human.groups for r in g.datasets()]
    assert sorted(names) == sorted(r["name"] for r in rows)
    assert len(names) == len(set(names))
    assert human.references.datasets()
    assert human.expression.sources(context="normal") == ()
    for source in human.expression.sources():
        c = source.expression_contract
        assert c.reference_status == "legacy_human"
        assert c.reference is None
        with pytest.raises(ValueError, match="legacy_human"):
            c.require_reference()


def test_species_cli_and_no_human_fetch(capsys, monkeypatch):
    assert cli.main(["species", "dog"]) == 0
    assert all(r["taxon_id"] == 9615 for r in json.loads(capsys.readouterr().out))
    monkeypatch.setattr(catalog, "fetch", lambda *a, **kw: pytest.fail("human fallback"))
    assert cli.main(["data", "fetch", "hpa", "--species", "dog"]) == 1
    assert "no human fallback" in capsys.readouterr().err
    assert cli.main(["data", "list", "--species", "dog"]) == 0
    assert "hpa_normal_tissue" not in capsys.readouterr().out


def test_contract_roundtrip_tampering_and_fingerprint(contract):
    record = json.loads(json.dumps(contract.as_dict()))
    assert ExpressionContract.from_dict(record).as_dict() == record
    for key, value in [("taxon_id", 9606), ("source_version", "v2"), ("extra", 1)]:
        with pytest.raises(ValueError):
            ExpressionContract.from_dict({**record, key: value})
    for changed in [
        replace(contract, source_version="v2"),
        replace(contract, context="tumor"),
        replace(contract, reference=replace(contract.reference, annotation_release="107")),
        replace(contract, measurement=replace(contract.measurement, unit="FPKM")),
        replace(contract, content_provenance={"matrix": "c" * 64}),
    ]:
        assert changed.contract_id != contract.contract_id
    with pytest.raises(TypeError):
        contract.content_provenance["matrix"] = "d" * 64


@pytest.mark.parametrize(
    "unit", ["counts", "TMM-CPM", "ERCC-normalized counts", "FPKM", "TPM", "pseudo-TPM", "iBAQ"]
)
def test_native_measurements_stay_distinct(contract, unit):
    other = replace(contract, measurement=replace(contract.measurement, unit=unit))
    assert ExpressionContract.from_dict(other.as_dict()).measurement.unit == unit
    if unit != "TPM":
        with pytest.raises(ValueError, match="measurement semantics"):
            validate_expression_join(contract, other, combine_values=True)


def test_transform_retains_assumptions_and_denominator_loss(contract):
    transform = MeasurementTransform(
        "length proxy", "TMM-CPM", "pseudo-TPM", ("gene lengths incomplete",), 37149, 16172
    )
    measurement = ExpressionMeasurement("gene", "pseudo-TPM", "linear", True, (transform,))
    record = replace(contract, measurement=measurement).as_dict()
    restored = ExpressionContract.from_dict(record)
    assert restored.measurement == measurement
    assert restored.measurement.transformations[0].output_features == 16172
    with pytest.raises(ValueError, match="output must match"):
        replace(measurement, unit="TPM")


@pytest.mark.parametrize("status", ["unresolved", "legacy_human"])
def test_missing_reference_stays_unknown(contract, status):
    unknown = replace(
        contract,
        reference=None,
        reference_status=status,
        taxon_id=9606 if status == "legacy_human" else 9615,
    )
    assert ExpressionContract.from_dict(unknown.as_dict()).reference is None
    with pytest.raises(ValueError, match=status):
        validate_expression_join(unknown, unknown)
    with pytest.raises(ValueError, match="Legacy human"):
        replace(unknown, taxon_id=9615, reference_status="legacy_human")


def test_join_and_derived_summary_guards(contract):
    validate_expression_join(contract, contract, combine_values=True)
    human = replace(contract, taxon_id=9606, reference=replace(contract.reference, taxon_id=9606))
    with pytest.raises(ValueError, match="taxon mismatch"):
        validate_expression_join(contract, human)
    other = replace(
        contract,
        reference=replace(contract.reference, annotation_source="Ensembl", annotation_release="98"),
        feature_namespace="Ensembl",
    )
    with pytest.raises(ValueError, match="reference/namespace mismatch"):
        validate_expression_join(contract, other)
    crosswalk = ReferenceCrosswalk(
        contract.reference,
        other.reference,
        "gene",
        "NCBI Gene",
        "Ensembl",
        "fixture-v1",
        {"reviewed_mapping": "c" * 64},
        (("gene1", "a"), ("gene1", "b"), ("gene2", "b")),
    )
    assert ReferenceCrosswalk.from_dict(json.loads(json.dumps(crosswalk.as_dict()))).targets(
        "gene1"
    ) == ("a", "b")
    assert crosswalk.targets("missing") == ()
    validate_expression_join(contract, other, crosswalk=crosswalk)
    with pytest.raises(ValueError, match="abundance allocation"):
        validate_expression_join(contract, other, crosswalk=crosswalk, combine_values=True)
    with pytest.raises(ValueError, match="directed"):
        validate_expression_join(other, contract, crosswalk=crosswalk)
    for level in ("transcript", "protein", "protein_group"):
        with pytest.raises(ValueError, match="cannot establish"):
            contract.require_feature_level(level)


def test_source_contract_cannot_disagree(contract):
    kwargs = {
        "id": "fixture",
        "category": "expression",
        "cancer_codes": (),
        "source_type": "fixture",
        "taxon_id": 9615,
        "context": "normal",
        "unit": "TPM",
        "source_version": "v1",
        "source_scale_class": "linear",
        "tpm_proxy": False,
    }
    source = ExpressionSource(**kwargs, contract=contract)
    assert source.expression_contract == contract
    with pytest.raises(ValueError, match="explicit contract"):
        ExpressionSource(**kwargs)
    with pytest.raises(ValueError, match="disagree"):
        replace(source, unit="FPKM")

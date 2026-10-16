"""One evidence flow across source formats, taxa and frozen archive versions."""

import copy
import json
from dataclasses import replace
from importlib.resources import files

import pytest

import oncoref
from oncoref.cta_evidence_adapters import normalize_cta_snapshot
from oncoref.expression_contract import ExpressionContract, ExpressionMeasurement
from oncoref.reference_identity import ReferenceIdentity
from oncoref.species_antigens import CTAEvidenceSnapshot


@pytest.fixture(scope="module")
def human():
    return oncoref.human.cta_snapshot()


@pytest.fixture(scope="module")
def dog():
    return oncoref.canine.cta_snapshot()


def reseal(record):
    record = copy.deepcopy(record)
    record.pop("snapshot_id", None)
    return CTAEvidenceSnapshot.create(record)


@pytest.mark.parametrize("record", [[], None, 1, True, "not a snapshot"])
def test_nonrecord_snapshot_inputs_raise_value_error(record):
    with pytest.raises(ValueError):
        CTAEvidenceSnapshot.from_dict(record)


def test_common_fact_shape_preserves_human_and_canine_somatic_signal(human, dog):
    human_row = human.evidence(candidate_id="ENSG00000126890")[0]  # CTAG2
    human_fact = next(f for f in human_row["facts"] if f["source_field"] == "rna_heart_max_ntpm")
    dog_row = next(
        e
        for e in dog.evidence(candidate_id="CTCFL")
        if e["tissue"] == "kidney" and e["locator"].startswith("Table 2")
    )
    dog_fact = dog_row["facts"][0]
    assert set(human_fact) == set(dog_fact)
    assert human_fact["observation"]["value"] == 5.1
    assert dog_fact["observation"]["value"] == 2.2
    assert human_fact["measurement"]["unit"] == "nTPM"
    assert dog_fact["measurement"]["unit"] == "percent of testis-normalized expression"
    assert human_fact["tissue_role"] == dog_fact["tissue_role"] == "somatic"
    assert human_fact["statistic"] == "tissue_maximum"
    assert dog_fact["statistic"] == "native_observation"
    assert human_row["curation"][0]["decision"] == "excluded"
    assert dog_row["curation"] == []
    for row in (human_row, dog_row):
        coverage = row["candidates"][0]["evidence_coverage"]
        assert "quantification_reference" in coverage["missing"]
        assert coverage["independent_donor_count"] is None
        assert row["source"]["sha256"] == row["contract"]["content_provenance"]["source"]


def test_every_historical_human_decision_and_source_field_is_preserved(human):
    evidence = human.evidence()
    original = json.loads(oncoref.cta_evidence().to_json(orient="records", double_precision=15))
    assert [e["source_fields"] for e in evidence] == original
    included = {
        e["candidate_ids"][0] for e in evidence if e["curation"][0]["decision"] == "included"
    }
    assert included == oncoref.cta_gene_ids()
    for row in evidence:
        policy = row["curation"][0]
        assert policy["policy_id"] == "oncoref-human-cta-core"
        assert policy["policy_version"] == row["source"]["version"]
        assert "rule" in policy
    missing = next(e for e in evidence if e["source_fields"]["rna_testis_ntpm"] is None)
    fact = next(f for f in missing["facts"] if f["source_field"] == "rna_testis_ntpm")
    assert fact["observation"]["status"] == "not_reported"
    assert fact["observation"]["value"] is None


def test_human_summary_zero_and_protein_flags_are_not_measured_negatives(human):
    row = human.evidence(candidate_id="ENSG00000185686")[0]  # PRAME
    zero = next(f for f in row["facts"] if f["source_field"] == "rna_placenta_ntpm")
    assert zero["observation"]["value"] == 0.0
    assert zero["statistic"] == "tissue_summary"
    assert not any(f["observation"]["status"] == "no_amplification" for f in row["facts"])
    assert (
        "measured_negative_or_censored_observations"
        in row["candidates"][0]["evidence_coverage"]["missing"]
    )
    protein = [f for f in row["facts"] if f["statistic"] == "tissue_detection_summary"]
    assert [f["tissue"] for f in protein] == ["testis"]
    assert protein[0]["measurement"]["feature_level"] == "protein_group"
    assert protein[0]["observation"]["status"] == "qualitative"
    assert protein[0]["observation"]["value"] is None
    assert not any(f["source_field"] in {"protein_ovary", "protein_placenta"} for f in row["facts"])


def test_native_censoring_counts_lineage_and_group_scope_are_lossless(dog):
    for row in dog.evidence():
        assert len(row["facts"]) == 1
        fact = row["facts"][0]
        assert fact["observation"] == row["observation"]
        assert fact["measurement"] == row["contract"]["measurement"]
    rows = dog.evidence(candidate_id="NT5C1B")
    censored = next(
        r for r in rows if r["tissue"] == "kidney" and r["locator"].startswith("Table 2")
    )
    assert censored["facts"][0]["observation"]["status"] == "below_limit"
    assert censored["facts"][0]["observation"]["value"] == 0.01
    assert censored["technical_replicates"] == 3
    assert censored["donor_count"] is None
    grouped = next(r for r in dog.evidence(candidate_id="AKR1E2") if len(r["candidate_ids"]) == 3)
    assert grouped["candidate_ids"] == ["AKR1E2", "SPECC1", "TPX2"]
    assert len(grouped["facts"]) == 1
    assert len(grouped["candidates"]) == 3
    same = next(r for r in dog.evidence(candidate_id="TPX2") if r["id"] == grouped["id"])
    assert same["facts"] == grouped["facts"]
    mage = next(r for r in dog.candidates() if r["id"] == "MAGE")
    assert "unique_complete_gene_mapping" in mage["evidence_coverage"]["missing"]
    assert len(mage["gene_ids"]) == 2


def test_capabilities_describe_real_gaps_without_taxon_admission_rules(human, dog):
    h = human.evidence(candidate_id="ENSG00000185686")[0]["candidates"][0]
    d = dog.evidence(candidate_id="PRAME")[0]["candidates"][0]
    assert "sequence_reference_assets" in h["evidence_coverage"]["missing"]
    assert "transcript_protein_pairs" in h["evidence_coverage"]["missing"]
    assert "sequence_reference_assets" not in d["evidence_coverage"]["missing"]
    assert "transcript_protein_pairs" not in d["evidence_coverage"]["missing"]
    for candidate in (h, d):
        coverage = candidate["evidence_coverage"]
        assert coverage["status"] == "partial"
        assert coverage["unreported_tissues"] is None
        assert "complete_tissue_panel" in coverage["missing"]
        assert "independent_sample_total" in coverage["missing"]
        assert "independent_donor_total" in coverage["missing"]
    mouse = oncoref.mouse.cta_snapshot()
    assert mouse.as_dict()["schema_version"] == 2
    assert mouse.as_dict()["facts"] == mouse.as_dict()["curation"] == []
    assert mouse.coverage["status"] == "unavailable"


@pytest.mark.parametrize("species", [oncoref.human, oncoref.canine, oncoref.mouse])
def test_v2_replay_uses_frozen_facts_without_adapters_or_current_data(
    species, tmp_path, monkeypatch
):
    import oncoref.cta_evidence_adapters as adapters

    snapshot = species.cta_snapshot()
    expected = snapshot.evidence()
    path = tmp_path / "snapshot.json"
    snapshot.write(path)
    monkeypatch.setattr(
        adapters, "normalize_cta_snapshot", lambda *_: pytest.fail("Reinterpreted archive")
    )
    replay = species.cta_snapshot(path=path, expected_snapshot_id=snapshot.snapshot_id)
    assert replay.as_dict() == snapshot.as_dict()
    assert replay.evidence() == expected
    if expected:
        expected[0]["facts"].clear()
        assert replay.evidence()[0]["facts"]


def test_v1_archives_are_not_silently_upgraded(tmp_path, dog):
    raw = json.loads(files("oncoref").joinpath("data/species/dog/cta-evidence.json").read_text())
    old = CTAEvidenceSnapshot.from_dict(raw)
    path = tmp_path / "v1.json"
    old.write(path)
    replay = oncoref.canine.cta_snapshot(path=path, expected_snapshot_id=old.snapshot_id)
    assert replay.as_dict() == raw
    assert replay.evidence()[0]["facts"] == []
    assert "normalized_observations" in replay.candidates()[0]["evidence_coverage"]["missing"]
    assert normalize_cta_snapshot(replay).as_dict() == dog.as_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: d["facts"][0].update(evidence_id="missing"),
        lambda d: d["facts"].append(copy.deepcopy(d["facts"][0])),
        lambda d: d["facts"][0]["observation"].update(status="not_reported", value=0),
        lambda d: d["facts"][0]["measurement"].update(unit="TPM"),
        lambda d: d["facts"][0].update(tissue_role="reproductive"),
        lambda d: d["facts"][0].update(extra="unknown"),
    ],
)
def test_facts_cannot_contradict_parent_even_with_recomputed_digest(dog, mutation):
    record = dog.as_dict()
    mutation(record)
    with pytest.raises(ValueError):
        reseal(record)


def test_historical_decision_cannot_be_changed_independently_of_source(human):
    record = human.as_dict()
    policy = record["curation"][0]
    policy["decision"] = "excluded" if policy["decision"] == "included" else "included"
    with pytest.raises(ValueError, match="decision disagrees"):
        reseal(record)


def equivalent_native_snapshot(taxon):
    reference = ReferenceIdentity(
        taxon,
        "fixture-assembly",
        "fixture-provider",
        "1",
        "fixture-v1",
        {"gene": "fixture gene"},
        {"annotation": "a" * 64},
    )
    contract = ExpressionContract(
        taxon,
        "fixture-source",
        "1",
        "normal",
        reference,
        "verified",
        "fixture gene",
        ExpressionMeasurement("gene", "nTPM", "source_native", False),
        {"source": "b" * 64},
        "Synthetic test observation with an exact reference.",
    )
    source = {
        "id": "fixture-source",
        "version": "1",
        "url": "https://example.org/fixture",
        "doi": None,
        "pmid": None,
        "sha256": "b" * 64,
        "scope": "Test fixture",
        "study_id": "fixture-study",
        "lineage": ["fixture-study"],
    }
    mapping = {
        "reference_id": reference.reference_id,
        "gene_namespace": "fixture gene",
        "status": "unique",
        "method": "Exact reference fixture",
        "complete": True,
        "targets": [{"id": "1", "symbol": "FIXTURE", "sequences": []}],
        "notes": "Gene only",
    }
    candidate = {
        "id": "candidate",
        "label": "FIXTURE",
        "nomination_kind": "published_nomination",
        "derivation": [
            {
                "kind": "published_nomination",
                "source_id": source["id"],
                "locator": "Fixture",
                "note": "Fixture nomination",
            }
        ],
        "evidence_status": "limited",
        "mappings": [mapping],
        "evidence_ids": ["observation"],
        "notes": "Test fixture",
    }
    evidence = {
        "id": "observation",
        "candidate_ids": ["candidate"],
        "source_id": source["id"],
        "study_id": source["study_id"],
        "contract_id": contract.contract_id,
        "locator": "Fixture",
        "assay": "RNA sequencing",
        "tissue": "kidney",
        "tissue_role": "somatic",
        "sample_ids": ["sample1"],
        "sample_count": 1,
        "donor_count": 1,
        "technical_replicates": None,
        "observation": {
            "status": "reported",
            "value": 5.0,
            "reported": "5.0 nTPM",
            "positive_count": None,
        },
        "source_fields": {},
        "notes": "Fixture only",
    }
    return normalize_cta_snapshot(
        CTAEvidenceSnapshot.create(
            {
                "schema": "oncoref.cta-evidence",
                "schema_version": 1,
                "taxon_id": taxon,
                "source_version": "fixture-v1",
                "coverage": {"status": "partial", "limitations": ["One tissue"]},
                "references": [{"id": "fixture", "identity": reference.as_dict()}],
                "sources": [source],
                "contracts": [contract.as_dict()],
                "candidates": [candidate],
                "evidence": [evidence],
            }
        )
    )


def test_equivalent_evidence_has_identical_shared_policy_result_in_all_species():
    results = []
    coverages = []
    for taxon in (9606, 9615, 10090):
        snapshot = equivalent_native_snapshot(taxon)
        # A hypothetical consumer rule; no taxon checks or source-column parsing.
        facts = [f for e in snapshot.evidence() for f in e["facts"]]
        results.append(
            any(
                f["tissue_role"] == "somatic"
                and f["measurement"]["unit"] == "nTPM"
                and f["observation"]["status"] == "reported"
                and f["observation"]["value"] >= 1
                for f in facts
            )
        )
        coverages.append(snapshot.candidates()[0]["evidence_coverage"])
        ref = ReferenceIdentity.from_dict(snapshot.candidates()[0]["reference"])
        with pytest.raises(ValueError, match="taxon"):
            snapshot.evidence(reference=replace(ref, taxon_id=10090 if taxon != 10090 else 9606))
        with pytest.raises(ValueError, match="annotation"):
            snapshot.evidence(reference=replace(ref, content_provenance={"annotation": "c" * 64}))
    assert results == [True, True, True]
    assert coverages[0] == coverages[1] == coverages[2]


def test_verified_measurement_on_another_reference_is_not_an_identity_bridge():
    snapshot = equivalent_native_snapshot(9615)
    record = snapshot.as_dict()
    original = ExpressionContract.from_dict(record["contracts"][0])
    changed = replace(original, reference=replace(original.reference, annotation_release="2"))
    record["contracts"] = [changed.as_dict()]
    record["evidence"][0]["contract_id"] = changed.contract_id
    replay = reseal(record)
    assert (
        "matching_quantification_reference"
        in replay.candidates()[0]["evidence_coverage"]["missing"]
    )

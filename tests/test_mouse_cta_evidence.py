"""P1A is a strain-qualified nomination, not a measured mouse-wide CTA panel."""

import hashlib
import json
from dataclasses import replace
from importlib.resources import files

import pytest

from oncoref import mouse
from oncoref.reference_identity import ReferenceIdentity
from oncoref.species_antigens import CTAEvidenceSnapshot


def test_source_identity_and_checked_single_gene_reference():
    snapshot = mouse.cta_snapshot()
    (candidate,) = snapshot.candidates()
    assert candidate["id"] == "P1A"
    assert candidate["gene_ids"] == ["22037"]
    assert candidate["mapping"]["targets"] == [
        {
            "id": "22037",
            "symbol": "Trap1a",
            "sequences": [{"transcript_id": "NM_011635.2", "protein_id": "NP_035765.1"}],
        }
    ]
    reference = ReferenceIdentity.from_dict(candidate["reference"])
    assert reference.assembly_accession == "GCF_000001635.27"
    assert reference.annotation_release == "RS_2024_02"
    assert "single-gene subset" in reference.source_version
    assert mouse.cta_gene_ids(reference="p1a-refseq-rs2024-02") == {"22037"}
    receipt = json.loads(
        files("oncoref").joinpath("data/species/mouse/cta-curation.json").read_text()
    )
    audit = receipt["identity_audit"]
    assert audit["source_protein"] == "AAA40475.1"
    assert audit["source_protein_length"] == audit["target_protein_length"] == 224
    assert audit["source_protein_sha256"] == audit["target_protein_sha256"]
    assert audit["cds_identical"] is False and audit["rna_identical"] is False
    assert audit["refseq_translation_matches_protein"] is True
    assert len(receipt["inputs"]) == 7
    for wrong in (
        replace(reference, assembly_accession="GCF_001624185.1"),
        replace(reference, source_version="DBA/2 source sequence"),
        replace(reference, annotation_release="RS_2020_09"),
        replace(reference, content_provenance={"rna": "a" * 64, "protein": "b" * 64}),
        replace(reference, taxon_id=9606),
    ):
        with pytest.raises(ValueError):
            snapshot.evidence(reference=wrong)


def test_positive_normal_rna_keeps_assay_limits_strain_and_unknown_counts():
    evidence = mouse.cta_evidence()
    assert len(evidence) == 6
    normal = {e["tissue"]: e for e in evidence if e["contract"]["context"] == "normal"}
    assert set(normal) == {"testis", "liver", "lung", "spleen", "thymus"}
    for row in evidence:
        assert row["source"]["pmid"] == "10049945"
        assert row["study_id"] == "PMID:10049945"
        assert row["sample_count"] is row["donor_count"] is row["technical_replicates"] is None
        assert row["sample_ids"] == []
        assert row["contract"]["reference_status"] == "unresolved"
        assert row["contract"]["reference"] is None
        assert row["source_fields"]["expression_origin"] == "endogenous"
        assert row["source_fields"]["measured_isoform"] is None
        identity = row["source_fields"]["nomination_identity"]
        assert identity["source_rna"] == "M36387.1"
        assert identity["source_protein"] == "AAA40475.1"
        assert identity["source_protein_sha256"] == identity["target_protein_sha256"]
        assert identity["cds_identical"] is identity["rna_identical"] is False
        receipt = files("oncoref").joinpath("data/species/mouse/cta-curation.json").read_bytes()
        pins = row["contract"]["content_provenance"]
        assert pins["curation_receipt"] == hashlib.sha256(receipt).hexdigest()
        native = next(r for r in json.loads(receipt)["inputs"] if r["name"] == "M36387.1.gb")
        assert pins["nomination_sequence"] == native["sha256"]
        (fact,) = row["facts"]
        assert fact["measurement"]["feature_level"] == "gene"
        assert fact["measurement"]["unit"] is None
        assert fact["observation"]["status"] == "qualitative"
        assert fact["observation"]["value"] is fact["observation"]["positive_count"] is None
    assert all(e["source_fields"]["strain"] == "BALB/c" for e in normal.values())
    assert all(e["source_fields"]["substrain"] is None for e in normal.values())
    assert normal["testis"]["tissue_role"] == "reproductive"
    assert all(normal[t]["tissue_role"] == "somatic" for t in ("lung", "spleen", "liver", "thymus"))
    assert "some experiments" in normal["thymus"]["observation"]["reported"]
    (tumor,) = [e for e in evidence if e["tissue_role"] == "tumor"]
    assert tumor["source_fields"]["cell_line"] == "J558"
    assert tumor["source_fields"]["strain"] is None


def test_sparse_coverage_selection_inventory_and_frozen_replay_agree(tmp_path):
    snapshot = mouse.cta_snapshot()
    assert snapshot.coverage["status"] == "partial"
    coverage = snapshot.candidates()[0]["evidence_coverage"]
    assert coverage["unreported_tissues"] is None
    assert coverage["independent_sample_count"] is coverage["independent_donor_count"] is None
    assert {"complete_tissue_panel", "quantification_reference"} <= set(coverage["missing"])
    assert "sequence_reference_assets" not in coverage["missing"]
    path = tmp_path / "mouse.json"
    snapshot.write(path)
    assert CTAEvidenceSnapshot.read(path).evidence() == snapshot.evidence()
    frame = mouse.cta_selection()
    assert len(frame) == 1 and frame.category.isna().all()
    assert frame.assessment_status.eq("unassessed").all()
    assert mouse.cta_selection(category="CTA").empty
    assert mouse.cta_selection_coverage()["evidence_status"] == "partial"
    sources = mouse.cta_sources()
    assert set(sources.query("status == 'native_observations'").source_id) == {
        "sarma1999",
        "vandeneynde1991",
    }
    audit = mouse.cta_selection_audit()
    assert set(audit.memberships.entity_id) == {"nomination:10090:P1A"}
    assert audit.memberships.gene_ids.tolist() == [("22037",), ("22037",)]
    assert audit.explain("22037").decision.eq("unassessed").all()
    assert audit.funnel.eligible.isna().all()

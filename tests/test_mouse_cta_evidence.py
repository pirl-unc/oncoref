"""P1A is a strain-qualified nomination mapped to a complete Ensembl 115 mouse host reference."""

import gzip
import hashlib
import json
from dataclasses import replace
from importlib.resources import files

import pytest

from oncoref import mouse, species_assets
from oncoref.reference_identity import ReferenceIdentity
from oncoref.species_antigens import CTAEvidenceSnapshot


def _complete_host_roles(identity: dict) -> tuple[str, ...]:
    namespaces = identity["identifier_namespaces"]
    if namespaces == {
        "gene": "Ensembl gene",
        "transcript": "Ensembl transcript",
        "protein": "Ensembl protein",
    }:
        roles = ("protein", "gtf", "rna", "ncrna")
    elif namespaces in (
        {"gene": "NCBI Gene", "transcript": "NCBI RefSeq RNA", "protein": "NCBI RefSeq protein"},
        {"gene": "NCBI Gene", "transcript": "RefSeq", "protein": "RefSeq"},
    ):
        roles = ("protein", "feature_table")
    else:
        raise ValueError("Unsupported CTA reference identifier namespaces")
    if not set(roles) <= set(identity["content_provenance"]):
        raise ValueError(
            "CTA evidence requires a complete annotated protein inventory reference; "
            "candidate subsets are insufficient"
        )
    return roles


def test_source_identity_and_complete_host_reference():
    snapshot = mouse.cta_snapshot()
    (candidate,) = snapshot.candidates()
    assert candidate["id"] == "P1A"
    assert candidate["gene_ids"] == ["ENSMUSG00000051257"]
    assert candidate["mapping"]["complete"] is True
    assert candidate["mapping"]["targets"] == [
        {
            "id": "ENSMUSG00000051257",
            "symbol": "Trap1a",
            "sequences": [
                {
                    "transcript_id": "ENSMUST00000054534.4",
                    "protein_id": "ENSMUSP00000059448.4",
                }
            ],
        }
    ]
    reference = ReferenceIdentity.from_dict(candidate["reference"])
    assert reference.assembly_accession == "GCA_000001635.9"
    assert reference.annotation_source == "Ensembl"
    assert reference.annotation_release == "115"
    assert reference.source_version == "Ensembl-115"
    assert set(reference.content_provenance) == {"gtf", "rna", "ncrna", "protein"}
    assert _complete_host_roles(reference.as_dict()) == ("protein", "gtf", "rna", "ncrna")
    (registered,) = mouse.references.registered_references()
    assert registered["id"] == "ensembl115"
    assert registered["identity"] == reference.as_dict()
    assert mouse.cta_gene_ids(reference="ensembl115") == {"ENSMUSG00000051257"}

    receipt = json.loads(
        files("oncoref").joinpath("data/species/mouse/cta-curation.json").read_text()
    )
    audit = receipt["identity_audit"]
    qual = receipt["host_reference_qualification"]
    assert audit["source_protein"] == "AAA40475.1"
    assert audit["source_protein_length"] == audit["target_protein_length"] == 224
    assert audit["source_protein_sha256"] == audit["target_protein_sha256"]
    assert audit["cds_identical"] is False and audit["rna_identical"] is False
    assert audit["refseq_translation_matches_protein"] is True
    assert qual["reference_id"] == reference.reference_id
    assert qual["gene_id"] == "ENSMUSG00000051257"
    assert qual["transcript_id"] == "ENSMUST00000054534.4"
    assert qual["protein_id"] == "ENSMUSP00000059448.4"
    assert qual["protein_length"] == 224
    assert qual["protein_sha256"] == audit["source_protein_sha256"]
    assert qual["matches_native_aaa40475_protein"] is True
    assert qual["matches_refseq_np035765_protein"] is True
    assert qual["ncbi_gene_cross_reference"] == "22037"
    assert qual["annotated_protein_count"] == 66668
    assert len(receipt["inputs"]) == 7

    for wrong in (
        "p1a-refseq-rs2024-02",
        replace(reference, assembly_accession="GCF_000001635.27"),
        replace(reference, source_version="DBA/2 source sequence"),
        replace(reference, annotation_release="114"),
        replace(reference, content_provenance={"rna": "a" * 64, "protein": "b" * 64}),
        replace(reference, taxon_id=9606),
    ):
        with pytest.raises(ValueError):
            snapshot.evidence(reference=wrong)


def test_candidate_subset_references_rejected_by_complete_host_contract():
    single_gene_refseq = ReferenceIdentity(
        taxon_id=10090,
        assembly_accession="GCF_000001635.27",
        annotation_source="NCBI RefSeq",
        annotation_release="RS_2024_02",
        source_version="NM_011635.2/NP_035765.1 single-gene subset",
        identifier_namespaces={
            "gene": "NCBI Gene",
            "transcript": "NCBI RefSeq RNA",
            "protein": "NCBI RefSeq protein",
        },
        content_provenance={"rna": "a" * 64, "protein": "b" * 64},
    )
    incomplete_ensembl = ReferenceIdentity(
        taxon_id=10090,
        assembly_accession="GCA_000001635.9",
        annotation_source="Ensembl",
        annotation_release="115",
        source_version="Ensembl-115-incomplete",
        identifier_namespaces={
            "gene": "Ensembl gene",
            "transcript": "Ensembl transcript",
            "protein": "Ensembl protein",
        },
        content_provenance={"protein": "a" * 64, "gtf": "b" * 64},
    )
    for subset in (single_gene_refseq, incomplete_ensembl):
        with pytest.raises(ValueError, match="complete annotated protein inventory"):
            _complete_host_roles(subset.as_dict())


def test_offline_mouse_host_proteome_construction_from_registered_assets(tmp_path, monkeypatch):
    orfdata = pytest.importorskip("orfdata")
    pyensembl = pytest.importorskip("pyensembl")
    import requests

    p1a_seq = (
        "MSDNKKPDKAHSGSGGDGDGNRCNLLHRYSLEEILPYLGWLVFAVVTTSFLALQMFIDAL"
        "YEEQYERDVAWIARQSKRMSSVDEDEDDEDDEDDYYDDEDDDDDAFYDDEDDEEEELENL"
        "MDDESEDEAEEEMSVEMGAGAEEMGAGANCACVPGHHLRKNEVKCRMIYFFHDPNFLVSI"
        "PVNPKEQMECRCENADEEVAMEEEEEEEEEEEEEEMGNPDGFSP"
    )
    assert len(p1a_seq) == 224
    cds_seq = "ATG" + "GCT" * 223 + "TAA"
    housekeeping_seq = "M" + "A" * 63
    housekeeping_cds = "ATG" + "GCT" * 63 + "TAA"

    gtf_text = "\n".join(
        [
            "#!genome-build GRCm39",
            "#!genome-build-accession GCA_000001635.9",
            'X\tensembl\tgene\t138234432\t138238918\t.\t+\t.\tgene_id "ENSMUSG00000051257"; gene_version "4"; gene_name "Trap1a"; gene_biotype "protein_coding";',
            'X\tensembl\ttranscript\t138234432\t138238918\t.\t+\t.\tgene_id "ENSMUSG00000051257"; gene_version "4"; transcript_id "ENSMUST00000054534"; transcript_version "4"; gene_name "Trap1a"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";',
            'X\tensembl\texon\t138234432\t138235106\t.\t+\t.\tgene_id "ENSMUSG00000051257"; gene_version "4"; transcript_id "ENSMUST00000054534"; transcript_version "4"; exon_number "1";',
            'X\tensembl\tCDS\t138234432\t138235103\t.\t+\t0\tgene_id "ENSMUSG00000051257"; gene_version "4"; transcript_id "ENSMUST00000054534"; transcript_version "4"; exon_number "1"; protein_id "ENSMUSP00000059448"; protein_version "4";',
            '1\tensembl\tgene\t1000\t1200\t.\t+\t.\tgene_id "ENSMUSG00000000001"; gene_version "5"; gene_name "Gnai3"; gene_biotype "protein_coding";',
            '1\tensembl\ttranscript\t1000\t1200\t.\t+\t.\tgene_id "ENSMUSG00000000001"; gene_version "5"; transcript_id "ENSMUST00000000001"; transcript_version "5"; gene_name "Gnai3"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";',
            '1\tensembl\texon\t1000\t1194\t.\t+\t.\tgene_id "ENSMUSG00000000001"; gene_version "5"; transcript_id "ENSMUST00000000001"; transcript_version "5"; exon_number "1";',
            '1\tensembl\tCDS\t1000\t1191\t.\t+\t0\tgene_id "ENSMUSG00000000001"; gene_version "5"; transcript_id "ENSMUST00000000001"; transcript_version "5"; exon_number "1"; protein_id "ENSMUSP00000000001"; protein_version "5";',
            '2\tensembl\tgene\t2000\t2100\t.\t+\t.\tgene_id "ENSMUSG00000099999"; gene_version "1"; gene_name "Gm99999"; gene_biotype "lncRNA";',
            '2\tensembl\ttranscript\t2000\t2100\t.\t+\t.\tgene_id "ENSMUSG00000099999"; gene_version "1"; transcript_id "ENSMUST00000099999"; transcript_version "1"; gene_name "Gm99999"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";',
            "",
        ]
    )
    rna_text = (
        f">ENSMUST00000054534.4 cdna chromosome:GRCm39:X:138234432:138238918:1 gene:ENSMUSG00000051257.4 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Trap1a\n{cds_seq}\n"
        f">ENSMUST00000000001.5 cdna chromosome:GRCm39:1:1000:1200:1 gene:ENSMUSG00000000001.5 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Gnai3\n{housekeeping_cds}\n"
    )
    ncrna_text = (
        ">ENSMUST00000099999.1 ncrna chromosome:GRCm39:2:2000:2100:1 gene:ENSMUSG00000099999.1 gene_biotype:lncRNA transcript_biotype:lncRNA gene_symbol:Gm99999\n"
        + ("ACGT" * 25)
        + "\n"
    )
    pep_text = (
        f">ENSMUSP00000059448.4 pep chromosome:GRCm39:X:138234432:138238918:1 gene:ENSMUSG00000051257.4 transcript:ENSMUST00000054534.4 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Trap1a\n{p1a_seq}\n"
        f">ENSMUSP00000000001.5 pep chromosome:GRCm39:1:1000:1200:1 gene:ENSMUSG00000000001.5 transcript:ENSMUST00000000001.5 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Gnai3\n{housekeeping_seq}\n"
    )

    raw_by_role = {"gtf": gtf_text, "rna": rna_text, "ncrna": ncrna_text, "protein": pep_text}
    paths = {}
    provenance = {}
    for role, content in raw_by_role.items():
        target = tmp_path / f"{role}.gz"
        target.write_bytes(gzip.compress(content.encode("utf-8"), mtime=0))
        paths[role] = target
        provenance[role] = hashlib.sha256(target.read_bytes()).hexdigest()

    original = species_assets._manifest("mouse")
    manifest = json.loads(json.dumps(original))
    registration = manifest["references"][0]
    fixture_ref = replace(
        ReferenceIdentity.from_dict(registration["identity"]),
        source_version="Ensembl-115-fixture",
        content_provenance=provenance,
    )
    registration["identity"] = fixture_ref.as_dict()
    overrides = {}
    for asset in manifest["assets"]:
        role = next(r for r, aid in registration["assets"].items() if aid == asset["id"])
        asset["sha256"] = provenance[role]
        asset["size_bytes"] = paths[role].stat().st_size
        overrides["10090:" + asset["id"]] = str(paths[role])

    species_dir = tmp_path / "species"
    (species_dir / "mouse").mkdir(parents=True)
    (species_dir / "mouse" / "assets.json").write_text(json.dumps(manifest))
    monkeypatch.setattr(species_assets, "_MANIFEST_ROOT", species_dir)
    monkeypatch.setattr(
        requests, "get", lambda *a, **kw: pytest.fail("Network call on offline mouse host replay")
    )

    resolved_ref, resolved_paths = mouse.references.ensure_reference(
        "ensembl115", offline=True, max_bytes=0, local_overrides=overrides
    )
    assert resolved_ref == fixture_ref and resolved_paths == paths

    context = orfdata.ReferenceContext.from_genome(
        pyensembl.EnsemblRelease(int(resolved_ref.annotation_release), species="mouse")
    )
    assert (
        context.taxon_id,
        context.assembly_accession,
        context.annotation_name.casefold(),
        context.annotation_version,
    ) == (
        resolved_ref.taxon_id,
        resolved_ref.assembly_accession,
        resolved_ref.annotation_source.casefold(),
        resolved_ref.annotation_release,
    )
    proteome = orfdata.Proteome.from_ensembl(
        resolved_paths["protein"],
        resolved_paths["gtf"],
        resolved_paths["rna"],
        resolved_paths["ncrna"],
        reference=context,
        expected_sha256=dict(resolved_ref.content_provenance),
    )
    layers = [
        layer
        for layer in proteome.layers
        if layer.get("kind") == "annotated" and layer.get("role") == "host"
    ]
    assert len(layers) == 1
    assert dict(layers[0]["content_provenance"]) == dict(resolved_ref.content_provenance)
    assert layers[0]["identifier_namespaces"] == dict(resolved_ref.identifier_namespaces)
    annotated = [entry for entry in proteome if entry.kind == "annotated"]
    assert len(annotated) == 2 and layers[0]["added"] == 2
    by_acc = {entry.accession: entry for entry in annotated}
    assert by_acc["ENSMUSP00000059448.4"].gene_id == "ENSMUSG00000051257"
    assert by_acc["ENSMUSP00000059448.4"].transcript_id == "ENSMUST00000054534.4"
    assert by_acc["ENSMUSP00000059448.4"].sequence == p1a_seq


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
    assert audit.memberships.gene_ids.tolist() == [
        ("ENSMUSG00000051257",),
        ("ENSMUSG00000051257",),
    ]
    assert audit.explain("ENSMUSG00000051257").decision.eq("unassessed").all()
    assert audit.funnel.eligible.isna().all()

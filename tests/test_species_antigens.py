"""Species CTA evidence never substitutes human policy or invents missing facts."""

import copy
import gzip
import importlib.util
import json
from dataclasses import FrozenInstanceError, replace
from pathlib import Path

import pytest

import oncoref
from oncoref import catalog, cli
from oncoref.expression_contract import ExpressionContract
from oncoref.reference_identity import ReferenceIdentity, _identity
from oncoref.species_antigens import CTAEvidenceSnapshot


@pytest.fixture
def snapshot():
    return oncoref.canine.antigens.cta_snapshot()


def reseal(record):
    payload = {k: v for k, v in record.items() if k != "snapshot_id"}
    return {**payload, "snapshot_id": _identity("cta-evidence", record["schema_version"], payload)}


def test_species_collection_and_reference_boundary(snapshot):
    dog = oncoref.for_species(9615)
    assert dog.antigens.coverage()["coverage"] == "partial"
    rows = dog.antigens.cta_candidates(reference="refseq106")
    assert len(rows) == 14
    assert {r["taxon_id"] for r in rows} == {9615}
    reference = ReferenceIdentity.from_dict(rows[0]["reference"])
    assert reference.annotation_release == "106"
    assert snapshot.candidates(reference=reference) == rows
    assert snapshot.gene_ids(reference=reference) == dog.antigens.cta_gene_ids(
        reference="refseq106"
    )
    with pytest.raises(TypeError, match="reference"):
        dog.antigens.cta_gene_ids()
    for wrong in (
        "GRCh38",
        "ensembl98",
        replace(reference, annotation_release="107"),
        replace(reference, content_provenance={"gtf": "a" * 64}),
    ):
        with pytest.raises(ValueError):
            snapshot.candidates(reference=wrong)
    with pytest.raises(ValueError, match="taxon"):
        snapshot.evidence(reference=replace(reference, taxon_id=9606))
    with pytest.raises(ValueError, match="antigens group"):
        dog.expression.cta_snapshot()


def test_partial_mouse_and_unchanged_human_defaults(snapshot):
    mouse = oncoref.for_species("mouse").antigens
    assert mouse.coverage()["coverage"] == "partial"
    missing = mouse.cta_snapshot()
    assert missing.coverage["status"] == "partial"
    assert len(missing.as_dict()["candidates"]) == 1
    with pytest.raises(ValueError, match="annotation"):
        mouse.cta_candidates(reference="refseq106")
    # Historical flat human APIs still return their own human gene universe.
    human_ids = oncoref.cta_gene_ids()
    assert human_ids and all(g.startswith("ENSG") for g in human_ids)
    assert not set(human_ids) & snapshot.gene_ids(reference="refseq106")
    assert oncoref.human.antigens.datasets()
    assert oncoref.human.cta_snapshot().coverage["status"] == "partial"
    assert len(oncoref.mouse.cta_candidates()) == 1
    assert len(oncoref.mouse.cta_evidence()) == 6


def test_source_backed_counterevidence_and_native_units(snapshot):
    evidence = snapshot.evidence(reference="refseq106")
    kidney = {
        e["candidate_ids"][0]: e
        for e in evidence
        if e["tissue"] == "kidney" and e["locator"].startswith("Table 2")
    }
    for gene, expected in [("CTCFL", 2.2), ("DDX53", 4.9)]:
        row = kidney[gene]
        assert row["observation"]["value"] == expected
        assert row["tissue_role"] == "somatic"
        assert row["contract"]["measurement"]["unit"] == "percent of testis-normalized expression"
        assert row["donor_count"] is None and row["sample_count"] is None
        assert row["technical_replicates"] == 3
        assert row["source"]["pmid"] == "40359694"
    assert kidney["NT5C1B"]["observation"]["status"] == "below_limit"
    assert kidney["NT5C1B"]["observation"]["value"] == 0.01
    blank = kidney["PRAME"]["observation"]
    assert blank["status"] == "no_amplification" and blank["value"] is None
    placenta = next(
        e for e in evidence if e["tissue"] == "placenta" and e["candidate_ids"] == ["PAGE3-like"]
    )
    assert placenta["tissue_role"] == "reproductive"
    assert placenta["observation"]["value"] == 15.2
    pbmc = next(e for e in evidence if e["candidate_ids"] == ["IGF2BP3"])
    assert pbmc["source"]["pmid"] == "30854786"
    assert pbmc["contract"]["measurement"]["unit"] is None
    assert pbmc["observation"]["value"] is None
    assert pbmc["sample_count"] is None
    assert pbmc["tissue"] == "peripheral blood mononuclear cells"
    # The abstract doesn't say which of the three genes was high in each tissue.
    grouped = next(
        e
        for e in evidence
        if e["candidate_ids"] == ["AKR1E2", "SPECC1", "TPX2"] and e["tissue_role"] == "somatic"
    )
    assert grouped["tissue"] is None and grouped["observation"]["value"] is None
    for row in evidence:
        contract = ExpressionContract.from_dict(row["contract"])
        assert contract.taxon_id == 9615
        assert contract.content_provenance["source"] == row["source"]["sha256"]
        assert row["study_id"] in row["source"]["lineage"]
        with pytest.raises(ValueError, match="unresolved"):
            contract.require_reference()


def test_family_scope_ambiguous_mapping_and_all_sequence_associations(snapshot):
    candidates = {r["id"]: r for r in snapshot.candidates(reference="refseq106")}
    mage = candidates["MAGE"]
    assert mage["mapping"]["status"] == "ambiguous"
    assert not mage["mapping"]["complete"]
    assert set(mage["gene_ids"]) == {"403710", "492173"}
    assert candidates["MAGEA13P"]["nomination_kind"] == "assay_followup"
    assert candidates["MAGEB10"]["nomination_kind"] == "assay_followup"
    prame = candidates["PRAME"]["mapping"]["targets"][0]
    assert prame["id"] == "477552"
    assert len(prame["sequences"]) == 6
    assert {s["transcript_id"] for s in prame["sequences"]} == {
        f"XM_03857575{i}.1" for i in range(5, 10)
    } | {"XM_038575760.1"}
    mage_evidence = snapshot.evidence(reference="refseq106", candidate_id="MAGE")
    ihc = [
        e for e in mage_evidence if e["locator"] == "Table 3" and e["tissue"].endswith(" melanoma")
    ]
    assert sum(e["observation"]["positive_count"] for e in ihc) == 4
    assert sum(e["sample_count"] for e in ihc) == 33
    assert all(e["contract"]["measurement"]["feature_level"] == "protein_group" for e in ihc)
    assert all(
        e["contract"]["measurement"]["feature_level"] == "unknown"
        for e in mage_evidence
        if e["locator"].startswith("Table 2")
    )
    # A multi-target observation never becomes a separate observation per gene.
    assert all(e["candidate_ids"] == ["MAGE"] for e in mage_evidence)


def test_conflicting_reports_and_nomination_provenance_are_retained(snapshot):
    rows = snapshot.candidates()
    assert not any(r["nomination_kind"] == "expression_derived" for r in rows)
    assert sum(r["nomination_kind"] == "published_nomination" for r in rows) == 12
    with pytest.raises(ValueError, match="candidate"):
        snapshot.evidence(reference="refseq106", candidate_id="human-only")
    lmck = next(
        e
        for e in snapshot.evidence(reference="refseq106", candidate_id="MAGE")
        if e["sample_ids"] == ["LMCK"]
    )
    assert lmck["observation"]["value"] == 0.9
    assert lmck["observation"]["reported"] == "- 0.9"
    assert any("threshold" in s for s in snapshot.coverage["limitations"])


def test_offline_export_replay_and_immutable_results(snapshot, tmp_path, monkeypatch):
    import oncoref.species_antigens as module

    path = tmp_path / "cta.json"
    snapshot.write(path)
    monkeypatch.setattr(
        module, "files", lambda *a: pytest.fail("Offline replay read installed assets")
    )
    replay = oncoref.canine.antigens.cta_snapshot(
        path=path, expected_snapshot_id=snapshot.snapshot_id
    )
    assert replay.as_dict() == snapshot.as_dict()
    assert replay.evidence(reference="refseq106") == snapshot.evidence(reference="refseq106")
    rows = replay.candidates(reference="refseq106")
    rows[0]["mapping"]["targets"].clear()
    assert replay.candidates(reference="refseq106")[0]["gene_ids"]
    with pytest.raises(FrozenInstanceError):
        replay._serialized = "{}"
    with pytest.raises(ValueError, match="taxon"):
        oncoref.for_species("mouse").antigens.cta_snapshot(path=path)
    with pytest.raises(ValueError, match="Unexpected"):
        CTAEvidenceSnapshot.read(path, expected_snapshot_id="another-version")
    altered = json.loads(path.read_text())
    altered["source_version"] += "-tampered"
    path.write_text(json.dumps(altered))
    with pytest.raises(ValueError, match="digest"):
        CTAEvidenceSnapshot.read(path)


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: d.update(schema_version=3),
        lambda d: d.update(taxon_id=9606),
        lambda d: d.update(extra="unknown"),
        lambda d: d["candidates"][0]["mappings"][0].update(gene_namespace="Ensembl gene"),
        lambda d: d["candidates"][0]["mappings"][0].update(status="unresolved"),
        lambda d: d["candidates"][0]["evidence_ids"].append("missing"),
        lambda d: d["evidence"][0].update(study_id="independent-study"),
        lambda d: d["evidence"][0]["observation"].update(value=0),
        lambda d: d["evidence"][0].update(donor_count=-1),
        lambda d: d["evidence"][0].update(contract_id="missing"),
        lambda d: d["sources"][0].update(sha256="b" * 64),
        lambda d: d["candidates"][0].update(nomination_kind="expression_derived"),
        lambda d: d["candidates"][0].update(nomination_kind="legacy_human"),
        lambda d: d["coverage"].update(status="unavailable"),
    ],
)
def test_invalid_semantics_rejected_even_with_recomputed_digest(snapshot, mutation):
    record = snapshot.as_dict()
    mutation(record)
    with pytest.raises(ValueError):
        CTAEvidenceSnapshot.from_dict(reseal(record))


def test_ambiguous_and_unresolved_mapping_does_not_pick_one_gene(snapshot):
    record = snapshot.as_dict()
    candidate = record["candidates"][0]
    original = candidate["mappings"][0]["targets"][0]["id"]
    second = copy.deepcopy(candidate["mappings"][0]["targets"][0])
    second["id"] = "ambiguous-fixture"
    candidate["mappings"][0]["targets"].append(second)
    candidate["mappings"][0]["status"] = "ambiguous"
    replay = CTAEvidenceSnapshot.from_dict(reseal(record))
    assert original not in replay.gene_ids(reference="refseq106")
    assert set(replay.candidates()[0]["gene_ids"]) == {original, "ambiguous-fixture"}
    candidate["mappings"][0].update(targets=[], status="unresolved", complete=False)
    replay = CTAEvidenceSnapshot.from_dict(reseal(record))
    assert replay.candidates(reference="refseq106")[0]["gene_ids"] == []
    assert replay.evidence(reference="refseq106", candidate_id=candidate["id"])


def test_cli_inventory_includes_bundled_evidence_with_cache_override(tmp_path, capsys):
    assert (
        cli.main(
            [
                "data",
                "list",
                "--species",
                "dog",
                "--group",
                "antigens",
                "--cache-root",
                str(tmp_path),
                "--json",
            ]
        )
        == 0
    )
    rows = json.loads(capsys.readouterr().out)
    assert rows == catalog.inventory(species="dog", group="antigens")
    assert rows[0]["held"] == "wheel" and rows[0]["available"]


def test_builder_keeps_many_to_many_sequence_pairs_and_rejects_changed_inputs(tmp_path):
    spec = importlib.util.spec_from_file_location(
        "cta_builder", Path(__file__).parents[1] / "scripts/build_canine_cta_evidence.py"
    )
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    source = tmp_path / "features.gz"
    with gzip.open(source, "wt") as stream:
        stream.write("# feature\tGeneID\tsymbol\tproduct_accession\trelated_accession\n")
        stream.write("mRNA\t1\tTEST\tRNA1.1\tPROTEIN1.1\n")
        stream.write("CDS\t1\tTEST\tPROTEIN1.1\tRNA1.1\n")
        stream.write("mRNA\t1\tTEST\tRNA2.1\tPROTEIN1.1\n")
        stream.write("mRNA\t2\tTEST\tRNA3.1\tPROTEIN2.1\n")
    genes = builder.reference_targets(source)
    assert len(genes) == 2
    assert genes["1"]["sequences"] == [
        {"transcript_id": "RNA1.1", "protein_id": "PROTEIN1.1"},
        {"transcript_id": "RNA2.1", "protein_id": "PROTEIN1.1"},
    ]
    with pytest.raises(ValueError, match="checksum"):
        builder.build(source, source, source)


def test_simple_species_methods_have_identical_signatures_and_record_shapes():
    import inspect

    for method in ("cta_candidates", "cta_evidence", "cta_gene_ids", "cta_snapshot"):
        signatures = [
            inspect.signature(getattr(species, method))
            for species in (oncoref.human, oncoref.canine, oncoref.mouse)
        ]
        assert signatures[0] == signatures[1] == signatures[2]
        for species in (oncoref.human, oncoref.canine, oncoref.mouse):
            assert getattr(species, method) == getattr(species.antigens, method)
    assert set(oncoref.human.cta_candidates()[0]) == set(oncoref.canine.cta_candidates()[0])
    assert set(oncoref.human.cta_evidence()[0]) == set(oncoref.canine.cta_evidence()[0])
    with pytest.raises(TypeError, match="reference"):
        oncoref.canine.cta_gene_ids()
    with pytest.raises(ValueError, match="explicit reference"):
        oncoref.human.cta_gene_ids(reference=None)


def test_human_adapter_preserves_existing_facts_and_offline_replay(tmp_path):
    original = json.loads(oncoref.cta_evidence().to_json(orient="records", double_precision=15))
    evidence = oncoref.human.cta_evidence()
    assert [e["source_fields"] for e in evidence] == original
    assert all(e["contract"]["reference_status"] == "legacy_human" for e in evidence)
    snapshot = oncoref.human.cta_snapshot()
    assert snapshot.taxon_id == 9606
    assert oncoref.human.cta_gene_ids(reference="ensembl115") == {
        r["Ensembl_Gene_ID"] for r in original
    }
    with pytest.raises(ValueError, match="annotation"):
        oncoref.human.cta_gene_ids(reference="refseq106")
    path = tmp_path / "human.json"
    snapshot.write(path)
    assert oncoref.human.cta_snapshot(path=path).as_dict() == snapshot.as_dict()
    # Historical default is a curated subset, not redefined by candidate enumeration.
    assert oncoref.cta_gene_ids() < oncoref.human.cta_gene_ids(reference="ensembl115")


def test_rbmy_corrected_transcripts_propagate_without_merging_genes():
    expected = {
        "ENSG00000169800": "ENST00000303766",
        "ENSG00000226941": "ENST00000250831",
    }
    candidates = {row["id"]: row for row in oncoref.human.cta_candidates(reference="ensembl115")}
    evidence = {
        row["source_fields"]["Ensembl_Gene_ID"]: row["source_fields"]
        for row in oncoref.human.cta_evidence()
    }
    assert set(expected) <= oncoref.cta_gene_ids()
    for gid, tid in expected.items():
        assert candidates[gid]["gene_ids"] == [gid]
        (target,) = candidates[gid]["mapping"]["targets"]
        assert target["id"] == gid
        assert target["symbol"] == evidence[gid]["Symbol"]
        assert any(
            s["transcript_id"].startswith(tid + ".") and s["protein_id"]
            for s in target["sequences"]
        )
        assert len(target["sequences"]) > 1
        assert evidence[gid]["Canonical_Transcript_ID"] == tid

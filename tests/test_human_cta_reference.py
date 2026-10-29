"""Exact sequence references stay separate from historical human assay evidence."""

import gzip
import hashlib
import importlib.util
import json
from dataclasses import replace
from pathlib import Path

import pytest

import oncoref
from oncoref.human_cta_reference import human_cta_reference
from oncoref.reference_identity import ReferenceIdentity
from oncoref.species_antigens import CTAEvidenceSnapshot
from oncoref.species_assets import registered_references

FIXTURE = Path(__file__).parent / "fixtures" / "human-cta-reference"
SPEC = importlib.util.spec_from_file_location(
    "build_human_cta_reference", Path(__file__).parents[1] / "scripts/build_human_cta_reference.py"
)
builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(builder)


@pytest.fixture
def inputs():
    manifest = json.loads((FIXTURE / "manifest.json").read_text())
    reference = replace(
        ReferenceIdentity.from_dict(manifest["full_reference"]),
        source_version="Ensembl-115-real-fixture-slices",
        content_provenance=manifest["fixture_sha256"],
    )
    return (
        reference,
        {r: FIXTURE / (r + ".gz") for r in reference.content_provenance},
        manifest["genes"],
    )


def mutate(inputs, tmp_path, role, transform):
    reference, paths, genes = inputs
    path = tmp_path / (role + ".gz")
    path.write_bytes(gzip.compress(transform(gzip.decompress(paths[role].read_bytes())), mtime=0))
    reference = replace(
        reference,
        content_provenance={
            **reference.content_provenance,
            role: hashlib.sha256(path.read_bytes()).hexdigest(),
        },
    )
    return reference, {**paths, role: path}, genes


def test_all_versioned_isoforms_keep_distinct_rbmy_genes(inputs):
    result = builder.build(*inputs)
    ref, live = human_cta_reference()
    assert ref == ReferenceIdentity.from_dict(registered_references("human")[0]["identity"])
    assert set(ref.content_provenance) == {"gtf", "rna", "ncrna", "protein"}
    expected = {
        "ENSG00000169800": ("ENST00000303766.12", "ENSP00000307155.7", 3),
        "ENSG00000226941": ("ENST00000250831.6", "ENSP00000250831.5", 2),
    }
    for gid, (tid, pid, count) in expected.items():
        mapping = result["mappings"][gid]
        assert mapping["complete"] and mapping["status"] == "unique"
        (target,) = mapping["targets"]
        assert target["id"] == gid and len(target["sequences"]) == count
        assert {"transcript_id": tid, "protein_id": pid} in target["sequences"]
        assert any(s["protein_id"] is None for s in target["sequences"])
        assert target == live[gid]["targets"][0]
        assert result["sequence_audit"][pid]["transcript"] == tid
        assert result["sequence_audit"][tid]["gene"].startswith(gid + ".")


def test_no_symbol_or_missing_gene_fallback(inputs):
    ref, paths, genes = inputs
    result = builder.build(ref, paths, [*genes, "PRAME", "ENSG00000000000"])
    for gid in ("PRAME", "ENSG00000000000"):
        assert result["mappings"][gid]["status"] == "unresolved"
        assert result["mappings"][gid]["targets"] == []
        assert result["mappings"][gid]["complete"] is False


def test_source_hash_and_wrong_species_fail_before_parsing(inputs, tmp_path):
    ref, paths, genes = inputs
    bad = tmp_path / "changed.gz"
    bad.write_bytes(paths["rna"].read_bytes() + b"changed")
    with pytest.raises(ValueError, match="checksum"):
        builder.build(ref, {**paths, "rna": bad}, genes)
    with pytest.raises(ValueError, match="complete Ensembl"):
        builder.build(replace(ref, taxon_id=9615), paths, genes)


@pytest.mark.parametrize("field", [b"gene:ENSG00000169800.14", b"transcript:ENST00000303766.12"])
def test_header_identity_conflicts_are_not_resolved_by_symbol(inputs, tmp_path, field):
    changed = mutate(inputs, tmp_path, "protein", lambda b: b.replace(field, field + b"9"))
    with pytest.raises(ValueError, match="identity conflict"):
        builder.build(*changed)


def test_conflicting_duplicate_sequence_identity_fails(inputs, tmp_path):
    def duplicate(data):
        record = b">" + data.split(b">")[1]
        header, sequence = record.split(b"\n", 1)
        return data + header + b"\nA" + sequence[1:]

    with pytest.raises(ValueError, match="conflicting FASTA"):
        builder.build(*mutate(inputs, tmp_path, "protein", duplicate))


def test_missing_sequence_retains_annotation_and_incomplete_status(inputs, tmp_path):
    def remove(data):
        return b"".join(
            b">" + r for r in data.split(b">")[1:] if not r.startswith(b"ENSP00000307155.7 ")
        )

    result = builder.build(*mutate(inputs, tmp_path, "protein", remove))
    mapping = result["mappings"]["ENSG00000169800"]
    assert mapping["complete"] is False
    assert mapping["status"] == "unique"
    assert {"transcript_id": "ENST00000303766.12", "protein_id": "ENSP00000307155.7"} in mapping[
        "targets"
    ][0]["sequences"]
    assert result["sequence_audit"]["ENSP00000307155.7"] is None
    assert "ENSP00000307155.7" in mapping["notes"]


def test_multiple_annotation_loci_remain_ambiguous(inputs, tmp_path):
    def duplicate_locus(data):
        gene = next(
            line
            for line in data.splitlines()
            if b"\tgene\t" in line and b'gene_id "ENSG00000169800"' in line
        )
        fields = gene.split(b"\t")
        fields[0] = b"synthetic-alternative-locus"
        return data + b"\t".join(fields) + b"\n"

    result = builder.build(*mutate(inputs, tmp_path, "gtf", duplicate_locus))
    mapping = result["mappings"]["ENSG00000169800"]
    assert mapping["status"] == "ambiguous"
    assert "synthetic-alternative-locus" in mapping["notes"]
    assert len(mapping["targets"][0]["sequences"]) == 3


def test_live_snapshot_keeps_historical_facts_and_exact_reference_on_replay(tmp_path, monkeypatch):
    snapshot = oncoref.human.cta_snapshot()
    reference = ReferenceIdentity.from_dict(snapshot.candidates()[0]["reference"])
    historical = json.loads(oncoref.cta_evidence().to_json(orient="records", double_precision=15))
    assert [e["source_fields"] for e in snapshot.evidence()] == historical
    assert len(snapshot.candidates()) == 2533
    assert all(e["contract"]["reference_status"] == "legacy_human" for e in snapshot.evidence())
    for wrong in (
        replace(reference, taxon_id=9615),
        replace(reference, content_provenance={"rna": "a" * 64}),
        replace(reference, annotation_release="116"),
    ):
        with pytest.raises(ValueError):
            snapshot.candidates(reference=wrong)
    path = tmp_path / "snapshot.json"
    snapshot.write(path)
    import oncoref.human_cta_snapshot as live

    monkeypatch.setattr(
        live, "human_cta_reference", lambda: pytest.fail("Live reference read on replay")
    )
    replay = CTAEvidenceSnapshot.read(path)
    assert replay.as_dict() == snapshot.as_dict()
    assert replay.candidates(reference=reference) == snapshot.candidates(reference=reference)
    assert "quantification_reference" in replay.candidates()[0]["evidence_coverage"]["missing"]
    assert "sequence_reference_assets" not in replay.candidates()[0]["evidence_coverage"]["missing"]


def test_unknown_live_candidate_has_explicit_unresolved_reference(monkeypatch):
    import oncoref.human_cta_snapshot as live

    reference, mappings = human_cta_reference()
    gid = "ENSG00000169800"
    monkeypatch.setattr(
        live,
        "human_cta_reference",
        lambda: (reference, {k: v for k, v in mappings.items() if k != gid}),
    )
    snapshot = live.human_cta_snapshot.__wrapped__()
    mapping = next(c for c in snapshot.candidates() if c["id"] == gid)["mapping"]
    assert mapping["status"] == "unresolved" and mapping["targets"] == []
    assert mapping["complete"] is False


def test_human_assets_use_common_cli_and_catalog_contracts(inputs, tmp_path, monkeypatch, capsys):
    from oncoref import catalog, cli, species_assets

    reference, paths, _ = inputs
    registration = registered_references("human")[0]
    original = species_assets._manifest("human")
    manifest = json.loads(json.dumps(original))
    manifest["references"][0]["identity"] = reference.as_dict()
    overrides = {}
    for asset in manifest["assets"]:
        role = next(r for r, aid in registration["assets"].items() if aid == asset["id"])
        asset["sha256"] = reference.content_provenance[role]
        asset["size_bytes"] = paths[role].stat().st_size
        overrides["9606:" + asset["id"]] = str(paths[role])
    root = tmp_path / "species"
    (root / "human").mkdir(parents=True)
    (root / "human" / "assets.json").write_text(json.dumps(manifest))
    monkeypatch.setattr(species_assets, "_MANIFEST_ROOT", root)
    override_path = tmp_path / "overrides.json"
    override_path.write_text(json.dumps(overrides))
    monkeypatch.setenv("ONCOREF_ASSET_OVERRIDES", str(override_path))
    import requests

    monkeypatch.setattr(
        requests, "get", lambda *a, **kw: pytest.fail("Network on offline reference replay")
    )
    result_ref, result_paths = oncoref.human.references.ensure_reference(
        "ensembl115", offline=True, max_bytes=0
    )
    assert result_ref == reference and result_paths == paths
    for action in ("metadata", "status", "path", "fetch", "versions", "provenance"):
        command = ["data", action, "ensembl115-protein", "--species", "human"]
        if action == "fetch":
            command += ["--offline", "--max-bytes", "0"]
        assert cli.main(command) == 0
        assert capsys.readouterr().out
    assert (
        cli.main(
            [
                "data",
                "list",
                "--species",
                "human",
                "--group",
                "references",
                "--json",
                "--asset-overrides",
                str(override_path),
            ]
        )
        == 0
    )
    rows = json.loads(capsys.readouterr().out)
    assert {a["id"] for a in manifest["assets"]} <= {r["name"] for r in rows}
    assert catalog.versions("ensembl115-protein") == ("Ensembl-115",)
    assert (
        catalog.provenance("ensembl115-protein", verify_content=True)["verification_state"]
        == "available"
    )
    assert catalog.path("ensembl115-protein") == paths["protein"]
    assert catalog.ensure("ensembl115-protein") == paths["protein"]
    assert catalog.fetch("ensembl115-protein") == []
    assert catalog.status("ensembl115-protein")[0]["present"]

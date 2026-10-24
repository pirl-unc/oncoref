"""Pinned release provenance through real PyEnsembl lookup and old JSON records."""

import gzip
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from pyensembl import EnsemblRelease

from oncoref import gene_identity

from .test_gene_identity import PRAME, Annotation, gene


def _snapshot(root):
    result = {}
    for entry in root.rglob("*"):
        state = entry.stat()
        result[str(entry.relative_to(root))] = (
            state.st_ino,
            state.st_mode,
            state.st_size,
            state.st_mtime_ns,
            state.st_ctime_ns,
            entry.read_bytes() if entry.is_file() else None,
        )
    return result


@pytest.mark.parametrize(
    "release,options,assembly",
    [
        (75, {}, "GRCh37"),
        (93, {}, "GRCh38"),
        ("2013_09", {"assembly": "GRCh37"}, "GRCh37"),
        ("2026_04", {}, "GRCh38"),
        ("2026_04", {"assembly_accession": "GCA_000001405.29"}, "GRCh38"),
        ("2026_04", {"reference_name": "study annotation"}, "GRCh38"),
    ],
)
@pytest.mark.skipif(
    not hasattr(EnsemblRelease, "assembly"),
    reason="Resolved assembly selections require PyEnsembl 2.29+ (Python 3.11+)",
)
def test_real_annotation_lookup_and_json_preserve_release_and_cache(
    tmp_path, monkeypatch, http_downloads, release, options, assembly
):
    cache = tmp_path / "cache"
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(cache))
    genome = EnsemblRelease(
        release, species="human", server="https://fixture.invalid/ensembl", **options
    )
    # Tiny synthetic data exercises the real downloader, index and gene lookup;
    # it is not a claim about this gene's coordinates in an upstream annotation.
    gtf = f'22\ttest\tgene\t1\t9\t.\t+\t.\tgene_id "{PRAME}"; gene_name "PRAME";\n'
    bodies = {genome.gtf_url: gzip.compress(gtf.encode())}
    bodies.update(
        (url, gzip.compress(b">synthetic\nATGAAATAA\n")) for url in genome.transcript_fasta_urls
    )
    bodies.update((url, gzip.compress(b">synthetic\nMK\n")) for url in genome.protein_fasta_urls)
    http_downloads(bodies.__getitem__)
    genome.download()
    genome.index()
    encoded_genome = genome.to_json()
    genome.close()
    before = _snapshot(cache)
    reopened = EnsemblRelease.from_json(encoded_genome)
    with reopened:
        record = gene_identity.resolve_gene_identity(PRAME, genome=reopened)
        assert record.verified and record.canonical_gene_id == PRAME
        assert record.annotation_assembly == assembly
        assert record.annotation_release == release
        assert type(record.annotation_release) is type(release)
        encoded = json.loads(json.dumps(record.as_dict()))
        restored = gene_identity.GeneIdentity(**encoded)
        assert json.loads(json.dumps(restored.as_dict())) == encoded
        assert type(restored.annotation_release) is type(release)
        assert reopened.reference_name == genome.reference_name
        assert reopened.assembly_accession == genome.assembly_accession
    assert _snapshot(cache) == before


def test_numbered_record_keeps_legacy_fields_and_reference_hashes():
    path = Path(__file__).parent / "fixtures/gene-identity/numbered-identity-oncoref-1.19.0.json"
    archived = json.loads(path.read_text())
    restored = gene_identity.GeneIdentity(**archived)
    assert json.loads(json.dumps(restored.as_dict())) == archived
    assert type(restored.annotation_release) is int
    current = gene_identity.resolve_gene_identity(PRAME, genome=Annotation([gene(PRAME)]))
    observed = json.loads(json.dumps(current.as_dict()))
    observed["oncoref_version"] = archived["oncoref_version"]
    assert observed == archived


@pytest.mark.parametrize(
    "release",
    [
        None,
        True,
        False,
        0,
        -1,
        1.5,
        "93",
        "latest",
        "2026_4",
        "2026-04",
        "2026_00",
        "2026_13",
        "0000_04",
        "2026_04\n",
        " 2026_04",
        "２０２６_04",
        [],
        {},
    ],
)
def test_invalid_release_fails_before_lookup(release):
    genome = Annotation([])
    genome.release = release
    with pytest.raises(ValueError, match="release"):
        gene_identity.resolve_gene_identity(PRAME, genome=genome)


@pytest.mark.parametrize("assembly", [None, "", "GRCm39", "GCA_000001405.29", [], {}])
def test_unknown_resolved_assembly_does_not_fall_back_to_cache_label(assembly):
    genome = Annotation([])
    genome.assembly = assembly
    with pytest.raises(ValueError, match="GRCh37/GRCh38"):
        gene_identity.resolve_gene_identity(PRAME, genome=genome)


@pytest.mark.parametrize(
    "species", [None, SimpleNamespace(), SimpleNamespace(latin_name="canis_familiaris")]
)
def test_missing_or_nonhuman_species_rejected(species):
    genome = Annotation([])
    genome.species = species
    genome.release = "2026_04"
    with pytest.raises(ValueError, match="human"):
        gene_identity.resolve_gene_identity(PRAME, genome=genome)

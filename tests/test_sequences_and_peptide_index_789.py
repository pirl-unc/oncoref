# Licensed under the Apache License, Version 2.0.
"""Regression tests for #789: dog Ensembl 98 pyensembl path, offline auto, and ~10x build perf."""

from __future__ import annotations

import gzip
import hashlib
import json
import sqlite3
import time
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from oncoref import (
    antigen_sets,
    canine,
    reference_data,
    sequence_bundle,
    species_target_selection,
)
from oncoref.peptide_subtraction import MatchPolicy, PeptideIndex


def _make_sqlite_pyensembl_genome(
    tmp_path: Path,
    *,
    release: int,
    latin_name: str,
    cds_rows: list[tuple[str, str, str, str, str]],
    fasta_entries: dict[str, str],
):
    db_path = tmp_path / f"{latin_name}_{release}.db"
    conn = sqlite3.connect(db_path)
    conn.execute(
        "CREATE TABLE CDS ("
        "gene_id TEXT NOT NULL, "
        "transcript_id TEXT NOT NULL, "
        "protein_id TEXT NOT NULL, "
        "transcript_biotype TEXT NOT NULL, "
        "protein_version TEXT NOT NULL)"
    )
    conn.execute("CREATE INDEX CDS_index_gene_id ON CDS (gene_id)")
    conn.executemany("INSERT INTO CDS VALUES (?, ?, ?, ?, ?)", cds_rows)
    conn.commit()

    class _FakeDB:
        # Real pyensembl.Database exposes column_exists + run_sql_query (not table_exists)
        def column_exists(self, table: str, col: str) -> bool:
            return table == "CDS" and col in {
                "gene_id",
                "transcript_id",
                "protein_id",
                "transcript_biotype",
                "protein_version",
            }

        def run_sql_query(self, sql: str, query_params=()):
            return conn.execute(sql, query_params).fetchall()

    stripped = {k.rsplit(".", 1)[0]: k for k in fasta_entries if "." in k}
    pseq = SimpleNamespace(fasta_dictionary=fasta_entries, _stripped_index=stripped)
    return SimpleNamespace(
        release=release,
        species=SimpleNamespace(latin_name=latin_name),
        db=_FakeDB(),
        protein_sequences=pseq,
        _conn=conn,
    )


def test_dog_ensembl_98_pyensembl_path_and_release_guardrails(tmp_path, monkeypatch):
    dog_cta = canine.target_selection(category="CTA").head(2).copy()
    gids = [gid for raw in dog_cta["member_gene_ids"] for gid in json.loads(raw)]
    cds_rows = []
    fasta_entries = {}
    for idx, gid in enumerate(gids):
        tid = f"ENSCAFT99000{idx:04d}"
        pid = f"ENSCAFP99000{idx:04d}"
        seq = f"MACDEFGHIKLMNPQRSTVWY{'A' * idx}"
        # Two exon rows per transcript verify duplicate-CDS collapse preserves GTF order once.
        cds_rows.append((gid, tid, pid, "protein_coding", "1"))
        cds_rows.append((gid, tid, pid, "protein_coding", "1"))
        fasta_entries[f"{pid}.1"] = seq

    dog_genome = _make_sqlite_pyensembl_genome(
        tmp_path,
        release=98,
        latin_name="canis_familiaris",
        cds_rows=cds_rows,
        fasta_entries=fasta_entries,
    )
    try:
        seq_df = antigen_sets._sequences(
            dog_cta,
            "CTA",
            species="dog",
            sequence_source="pyensembl",
            genome=dog_genome,
        )
        assert len(seq_df) == len(gids)
        assert seq_df.attrs["reference"] == "canis_familiaris:CanFam3.1:Ensembl:98"
        assert seq_df.attrs["sequence_source"] == "pyensembl"
        assert seq_df.attrs["missing_gene_ids"] == []

        # Wrong release for dog (115 instead of 98) fails closed with accurate message
        wrong_dog = SimpleNamespace(
            release=115,
            species=SimpleNamespace(latin_name="canis_familiaris"),
        )
        with pytest.raises(ValueError, match="require dog Ensembl 98"):
            antigen_sets._sequences(
                dog_cta, "CTA", species="dog", sequence_source="pyensembl", genome=wrong_dog
            )

        # Wrong release for human (98 instead of 115) still fails closed
        wrong_human = SimpleNamespace(
            release=98,
            species=SimpleNamespace(latin_name="homo_sapiens"),
        )
        with pytest.raises(ValueError, match="require human Ensembl 115"):
            antigen_sets._sequences(
                antigen_sets.proteins("CTA").head(1),
                "CTA",
                species="human",
                sequence_source="pyensembl",
                genome=wrong_human,
            )

        # Default EnsemblRelease instantiation when genome=None uses release 98 & species='dog'
        captured = {}

        def _fake_ensembl_release(release, species):
            captured["release"] = release
            captured["species"] = species
            return dog_genome

        import pyensembl

        monkeypatch.setattr(pyensembl, "EnsemblRelease", _fake_ensembl_release)
        auto_py_df = antigen_sets._sequences(
            dog_cta, "CTA", species="dog", sequence_source="pyensembl"
        )
        assert captured == {"release": 98, "species": "dog"}
        assert len(auto_py_df) == len(gids)
    finally:
        dog_genome._conn.close()


def test_sequence_source_auto_never_downloads_and_sequence_bundle_is_cached(tmp_path, monkeypatch):
    # 1. sequence_bundle.is_cached() returns False on empty dir without creating folders
    empty_cache = tmp_path / "no_such_bundle_cache"
    assert sequence_bundle.is_cached(cache_dir=empty_cache) is False
    assert not empty_cache.exists()

    # 2. packaged_species_sequences("human", allow_download=False) raises ReferenceDataError (no AttributeError)
    monkeypatch.setenv("ONCOREF_SEQUENCE_CACHE_DIR", str(empty_cache))
    with pytest.raises(reference_data.ReferenceDataError, match="not cached locally"):
        species_target_selection.packaged_species_sequences("human", allow_download=False)

    # 3. sequence_source="auto" passes allow_download=False and falls back to local pyensembl on ReferenceDataError
    dog_cta = canine.target_selection(category="CTA").head(1).copy()
    gid = json.loads(dog_cta.iloc[0]["member_gene_ids"])[0]
    seen_allow_download = []

    def _no_download_packaged(species_arg, *, allow_download=True):
        seen_allow_download.append(allow_download)
        if not allow_download:
            raise reference_data.ReferenceDataError("uncached offline test")
        raise AssertionError("sequence_source='auto' must not call allow_download=True")

    monkeypatch.setattr(
        species_target_selection, "packaged_species_sequences", _no_download_packaged
    )
    fallback_genome = _make_sqlite_pyensembl_genome(
        tmp_path,
        release=98,
        latin_name="canis_familiaris",
        cds_rows=[(gid, "ENSCAFT_FB1", "ENSCAFP_FB1", "protein_coding", "")],
        fasta_entries={"ENSCAFP_FB1.2": "MKTIIALSYIFCLVFA"},
    )
    import pyensembl

    monkeypatch.setattr(pyensembl, "EnsemblRelease", lambda rel, species: fallback_genome)
    try:
        res = antigen_sets._sequences(dog_cta, "CTA", species="dog", sequence_source="auto")
        assert seen_allow_download == [False]
        assert res.attrs["sequence_source"] == "pyensembl"
        assert list(res["sequence"]) == ["MKTIIALSYIFCLVFA"]
    finally:
        fallback_genome._conn.close()


def test_bulk_sqlite_full_scan_for_large_gene_sets_and_version_fallbacks(tmp_path):
    # >900 genes exercises the single-table scan branch in _pyensembl_translations_by_gene
    n_genes = 950
    protein_rows = pd.DataFrame(
        {
            "entity_id": [f"E{i}" for i in range(n_genes)],
            "representative_sequence_sha256s": [
                json.dumps([hashlib.sha256(b"ACDEFGHIK").hexdigest()]) for _ in range(n_genes)
            ],
            "member_gene_ids": [json.dumps([f"ENSG{i:08d}"]) for i in range(n_genes)],
            "category": ["Unsafe"] * n_genes,
        }
    )
    cds_rows = []
    fasta_entries = {}
    for i in range(n_genes - 1):
        gid = f"ENSG{i:08d}"
        tid = f"ENST{i:08d}"
        if i % 3 == 0:
            pid, pver = f"ENSP{i:08d}", "2"
            fasta_entries[f"{pid}.2"] = "ACDEFGHIK"
        elif i % 3 == 1:
            pid, pver = f"ENSP{i:08d}.4", ""
            fasta_entries[f"ENSP{i:08d}"] = "ACDEFGHIK"
        else:
            pid, pver = f"ENSP{i:08d}", ""
            fasta_entries[f"ENSP{i:08d}.7"] = "LMNPQRSTV"
        cds_rows.append((gid, tid, pid, "protein_coding", pver))

    genome = _make_sqlite_pyensembl_genome(
        tmp_path,
        release=115,
        latin_name="homo_sapiens",
        cds_rows=cds_rows,
        fasta_entries=fasta_entries,
    )
    try:
        out = antigen_sets._sequences(
            protein_rows, "Unsafe", species="human", sequence_source="pyensembl", genome=genome
        )
        assert len(out) == n_genes - 1
        assert out.attrs["missing_gene_ids"] == [f"ENSG{n_genes - 1:08d}"]
        assert out["matches_representative_sequence"].iloc[0] == True  # noqa: E712
        assert out["matches_representative_sequence"].iloc[2] == False  # noqa: E712
    finally:
        genome._conn.close()


def test_peptide_index_fast_grouping_and_multispecies_from_sequence_bundle(tmp_path):
    # 1. Verify PeptideIndex fast grouping preserves all duplicate-sequence occurrences & perf
    seqs = [f"ACDEFGHIKLMNPQ{i % 5000:05d}" for i in range(25000)]
    # Replace digits with canonical AAs so every sequence is canonical AA20
    aa_map = str.maketrans("0123456789", "ACDEFGHIKL")
    seqs = [s.translate(aa_map) for s in seqs]
    records = pd.DataFrame(
        {
            "entity_id": [f"PEG:{i}" for i in range(len(seqs))],
            "gene_id": [f"G{i}" for i in range(len(seqs))],
            "transcript_id": [f"T{i}" for i in range(len(seqs))],
            "sequence": seqs,
        }
    )
    t0 = time.perf_counter()
    idx = PeptideIndex(records, policy=MatchPolicy(max_distance=0), lengths=[9])
    elapsed = time.perf_counter() - t0
    assert elapsed < 3.0, f"PeptideIndex init on 25k records took {elapsed:.2f}s"
    probe = seqs[0][-9:]
    m = idx.matches(probe)
    # 25000 / 5000 = 5 identical sequence occurrences retained
    assert len(m) == 5
    assert set(m["gene_id"]) == {f"G{k}" for k in range(0, 25000, 5000)}

    # 2. PeptideIndex.from_sequence_bundle supports dog packaged proteome and custom non-human bundle
    dog_idx = PeptideIndex.from_sequence_bundle(species="dog", background="Other", lengths=[9])
    assert dog_idx.metadata["species"] == "dog"
    assert dog_idx.metadata["taxon_id"] == 9615

    with pytest.raises(ValueError, match="requires a SequenceBundle"):
        PeptideIndex.from_sequence_bundle(species="dog", roles=("exclude",), lengths=[9])

    # Custom synthetic canine SequenceBundle (taxon_id=9615) infers dog species instead of hardcoding human
    dog_other = canine.target_selection(category="Other")
    raw_records = []
    for row in dog_other.itertuples(index=False):
        seq = "ACDEFGHIKLMNPQRST"
        for gid in json.loads(row.member_gene_ids):
            raw_records.append(
                {
                    "record_id": f"REC_{gid}",
                    "reference_id": "dog:synthetic:98",
                    "sequence": seq,
                    "sequence_sha256": sequence_bundle.digest(seq.encode()),
                    "valid_segments": [
                        {
                            "start": 0,
                            "end": len(seq),
                            "sequence": seq,
                            "sha256": sequence_bundle.digest(seq.encode()),
                        }
                    ],
                    "shared_sequence": "review",
                    "occurrences": [
                        {
                            "occurrence_id": f"OCC_{gid}",
                            "gene_id": gid,
                            "transcript_id": f"T_{gid}",
                            "transcript_biotype": "protein_coding",
                            "shared_sequence": "review",
                        }
                    ],
                }
            )
    folder = tmp_path / "dog_syn_bundle"
    folder.mkdir()
    files = {
        "records.jsonl.gz": gzip.compress(
            ("\n".join(json.dumps(r) for r in raw_records) + "\n").encode(), mtime=0
        ),
        "coverage.json": json.dumps(
            {"reference_id": "dog:synthetic:98", "missing_translations": []}
        ).encode(),
        "decisions.csv.gz": gzip.compress(b"gene,role\n", mtime=0),
        "orfdata-proteome.json.gz": gzip.compress(b"{}", mtime=0),
    }
    manifest = {
        "schema": sequence_bundle.SCHEMA,
        "schema_version": 1,
        "reference": {"reference_id": "dog:synthetic:98", "taxon_id": 9615},
        "selection_policy_id": dog_other.attrs["policy_id"],
        "role_policy_id": sequence_bundle.ROLE_POLICY_ID,
        "record_count": len(raw_records),
        "occurrence_count": len(raw_records),
        "package_dependencies": {},
        "payloads": {
            n: {"bytes": len(v), "sha256": sequence_bundle.digest(v)} for n, v in files.items()
        },
    }
    manifest["bundle_id"] = sequence_bundle.content_id(manifest)
    for name, value in files.items():
        (folder / name).write_bytes(value)
    (folder / "manifest.json").write_text(json.dumps(manifest))
    loaded_dog_bundle = sequence_bundle.load(
        folder, expected_id=manifest["bundle_id"], require_current=False
    )
    syn_dog_idx = PeptideIndex.from_sequence_bundle(
        loaded_dog_bundle, background="Other", lengths=[9]
    )
    assert syn_dog_idx.metadata["species"] == "dog"
    assert not syn_dog_idx.matches("ACDEFGHIK").empty

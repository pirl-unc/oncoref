# Licensed under the Apache License, Version 2.0.
"""Regression tests for #791: HPA cancer prevalence caching + symbol, oncoref.dog, and representative_transcripts."""

from __future__ import annotations

import importlib
import json

import pandas as pd
import pytest

import oncoref
from oncoref import hpa_cancer, load_dataset, reference_data, species_target_selection
from oncoref.reference_identity import ReferenceIdentity


def test_hpa_cancer_prevalence_caching_copy_detachment_and_canonical_symbol(monkeypatch, tmp_path):
    ihc_path = tmp_path / "ihc.tsv"
    rna_path = tmp_path / "rna.tsv"
    pd.DataFrame(
        [
            {
                "Gene": "ENSG00000141510",  # TP53
                "Gene name": "TP53_HPA",
                "Cancer": "lung cancer",
                "High": 3,
                "Medium": 2,
                "Low": 1,
                "Not detected": 4,
            },
            {
                "Gene": "ENSG_UNMAPPED_TEST",
                "Gene name": "UNMAPPED_HPA",
                "Cancer": "lung cancer",
                "High": 0,
                "Medium": 0,
                "Low": 0,
                "Not detected": 5,
            },
        ]
    ).to_csv(ihc_path, sep="\t", index=False)

    pd.DataFrame(
        [
            {
                "gene_id": "ENSG00000141510",
                "cancer_code": "LUAD",
                "cancer": "LUAD (TCGA)",
                "cohort": "TCGA",
                "samples": 100,
                "sum_ptpm": 500.0,
                "expressed_samples_ptpm_ge_1": 80,
            },
            {
                "gene_id": "ENSG00000141510",
                "cancer_code": "LUSC",
                "cancer": "LUSC (TCGA)",
                "cohort": "TCGA",
                "samples": 50,
                "sum_ptpm": 250.0,
                "expressed_samples_ptpm_ge_1": 40,
            },
            {
                "gene_id": "ENSG00000147381",  # MAGEA4 (RNA-only in this synthetic slice)
                "cancer_code": "LUAD",
                "cancer": "LUAD (TCGA)",
                "cohort": "TCGA",
                "samples": 100,
                "sum_ptpm": 300.0,
                "expressed_samples_ptpm_ge_1": 35,
            },
            {
                "gene_id": "ENSG00000147381",
                "cancer_code": "LUSC",
                "cancer": "LUSC (TCGA)",
                "cohort": "TCGA",
                "samples": 50,
                "sum_ptpm": 200.0,
                "expressed_samples_ptpm_ge_1": 25,
            },
        ]
    ).to_csv(rna_path, sep="\t", index=False)

    monkeypatch.setattr(
        reference_data,
        "ensure",
        lambda name, *a, **kw: ihc_path if name == "hpa_cancer_ihc" else rna_path,
    )
    hpa_cancer._clear_cache()

    read_csv_calls: list[str] = []
    real_read_csv = pd.read_csv

    def _spy_read_csv(filepath_or_buffer, *args, **kwargs):
        p_str = str(filepath_or_buffer)
        if p_str in {str(ihc_path), str(rna_path)}:
            read_csv_calls.append(p_str)
        return real_read_csv(filepath_or_buffer, *args, **kwargs)

    monkeypatch.setattr(pd, "read_csv", _spy_read_csv)

    # First call parses once and attaches canonical symbol right after gene_id
    ihc1 = hpa_cancer.hpa_cancer_ihc_prevalence()
    assert list(ihc1.columns[:3]) == ["gene_id", "symbol", "hpa_gene_name"]
    assert ihc1.set_index("gene_id").loc["ENSG00000141510", "symbol"] == "TP53"
    assert ihc1.set_index("gene_id").loc["ENSG_UNMAPPED_TEST", "symbol"] == ""

    rna1 = hpa_cancer.hpa_cancer_rna_prevalence()
    assert list(rna1.columns[:2]) == ["gene_id", "symbol"]
    assert (
        rna1.set_index(["gene_id", "cancer_code"]).loc[("ENSG00000141510", "LUAD"), "symbol"]
        == "TP53"
    )
    assert (
        rna1.set_index(["gene_id", "cancer_code"]).loc[("ENSG00000147381", "LUAD"), "symbol"]
        == "MAGEA4"
    )
    assert read_csv_calls == [str(ihc_path), str(rna_path)]

    # Mutating returned DataFrames never corrupts the in-process cache
    ihc1.loc[:, "symbol"] = "MUTATED"
    rna1.loc[:, "symbol"] = "MUTATED"

    # Subsequent calls hit the in-process cache (0 additional pd.read_csv calls!)
    ihc2 = hpa_cancer.hpa_cancer_ihc_prevalence("ENSG00000141510")
    rna2 = hpa_cancer.hpa_cancer_rna_prevalence("ENSG00000147381", cohort="TCGA")
    assert read_csv_calls == [str(ihc_path), str(rna_path)]
    assert list(ihc2["symbol"]) == ["TP53"]
    assert list(rna2["symbol"]) == ["MAGEA4", "MAGEA4"]

    # include_symbol=False omits symbol column cleanly
    ihc_no_sym = hpa_cancer.hpa_cancer_ihc_prevalence(include_symbol=False)
    rna_no_sym = hpa_cancer.hpa_cancer_rna_prevalence(include_symbol=False)
    assert "symbol" not in ihc_no_sym.columns
    assert "symbol" not in rna_no_sym.columns

    # Comparison API also attaches canonical symbol (including RNA-only rows) without column collisions
    comp = hpa_cancer.hpa_cancer_rna_ihc_comparison(include_missing_ihc=True)
    assert list(comp.columns[:2]) == ["gene_id", "symbol"]
    comp_idx = comp.set_index(["gene_id", "cancer"])
    assert comp_idx.loc[("ENSG00000147381", "lung cancer"), "symbol"] == "MAGEA4"
    assert "symbol" not in hpa_cancer.hpa_cancer_rna_ihc_comparison(include_symbol=False).columns

    # Global derived-cache clear also clears HPA cancer prevalence caches
    load_dataset._clear_cache()
    _ = hpa_cancer.hpa_cancer_rna_prevalence("ENSG00000141510")
    assert len(read_csv_calls) == 3

    with pytest.raises(ValueError, match="include_symbol"):
        hpa_cancer.hpa_cancer_rna_prevalence(include_symbol="yes")  # type: ignore[arg-type]


def test_oncoref_dog_top_level_and_submodule_symmetry():
    import oncoref.dog as dog_imported
    from oncoref import canine, dog

    assert hasattr(oncoref, "dog")
    assert "dog" in oncoref.__all__
    assert oncoref.dog is canine
    assert dog is canine
    assert dog_imported is canine
    assert importlib.import_module("oncoref.dog") is canine
    assert oncoref.dog.species.name == "dog"
    assert oncoref.dog.species.taxon_id == 9615
    assert (
        oncoref.dog.target_selection("PRAME").iloc[0]["entity_id"]
        == canine.target_selection("PRAME").iloc[0]["entity_id"]
    )


def test_representative_transcripts_populated_across_dog_and_mouse_catalogs(monkeypatch):
    # Block proteome sequence downloads to prove representative_transcripts loads offline from wheel sidecar
    species_target_selection._SPECIES_CATALOG_CACHE.clear()
    with monkeypatch.context() as m_ctx:
        m_ctx.setattr(
            species_target_selection,
            "packaged_species_sequences",
            lambda *a, **kw: (_ for _ in ()).throw(AssertionError("Should not need proteome CSV")),
        )
        dog_cat = oncoref.dog.target_selection()
        mouse_cat = oncoref.mouse.target_selection()
    species_target_selection._SPECIES_CATALOG_CACHE.clear()

    assert len(dog_cat) == 20082
    assert "representative_transcripts" in dog_cat.columns
    dog_parsed = dog_cat["representative_transcripts"].map(json.loads)
    dog_translated = dog_cat["identity_translation_pass"]
    assert dog_parsed[dog_translated].map(len).ge(1).all()
    assert dog_parsed[~dog_translated].map(len).eq(0).all()
    assert all(t.startswith("ENSCAFT") for tids in dog_parsed[dog_translated] for t in tids)

    assert len(mouse_cat) == 21433
    assert "representative_transcripts" in mouse_cat.columns
    mouse_parsed = mouse_cat["representative_transcripts"].map(json.loads)
    mouse_translated = mouse_cat["identity_translation_pass"]
    assert mouse_parsed[mouse_translated].map(len).ge(1).all()
    assert mouse_parsed[~mouse_translated].map(len).eq(0).all()
    assert all(t.startswith("ENSMUST") for tids in mouse_parsed[mouse_translated] for t in tids)

    # Verify every selected representative transcript matches a representative sequence in proteome CSVs
    for sp, cat in (("dog", dog_cat), ("mouse", mouse_cat)):
        seq = species_target_selection.packaged_species_sequences(sp)
        valid_pairs = set(
            zip(
                seq.loc[seq["matches_representative_sequence"].fillna(False), "entity_id"],
                seq.loc[seq["matches_representative_sequence"].fillna(False), "transcript_id"],
            )
        )
        for eid, raw_tids in zip(cat["entity_id"], cat["representative_transcripts"]):
            for tid in json.loads(raw_tids):
                assert (eid, tid) in valid_pairs

    # group_reference_proteins also emits representative_transcripts directly on live evaluations
    ref = ReferenceIdentity(
        taxon_id=9615,
        assembly_accession="GCA_000002285.2",
        annotation_source="Ensembl",
        annotation_release="98",
        source_version="Ensembl-98",
        identifier_namespaces={
            "gene": "Ensembl gene",
            "protein": "Ensembl protein",
            "transcript": "Ensembl transcript",
        },
        content_provenance={"gtf": "a" * 64, "protein": "b" * 64},
    )
    groups_df, _seq_df = species_target_selection.group_reference_proteins(
        [
            {
                "gene_id": "ENSCAFG1",
                "symbol": "G1",
                "transcript_id": "ENSCAFT2",
                "protein_id": "ENSCAFP2",
                "sequence": "MACDEFGHIK",
            },
            {
                "gene_id": "ENSCAFG1",
                "symbol": "G1",
                "transcript_id": "ENSCAFT1",
                "protein_id": "ENSCAFP1",
                "sequence": "MACDEFGHIK",
            },
            {
                "gene_id": "ENSCAFG2",
                "symbol": "G2",
                "transcript_id": "",
                "protein_id": "",
                "sequence": "",
            },
        ],
        species="dog",
        reference=ref,
    )
    by_label = groups_df.set_index("protein_group_label")
    assert json.loads(by_label.loc["G1", "representative_transcripts"]) == ["ENSCAFT1"]
    assert json.loads(by_label.loc["G2", "representative_transcripts"]) == []

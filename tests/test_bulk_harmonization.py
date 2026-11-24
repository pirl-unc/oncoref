# Licensed under the Apache License, Version 2.0.
"""Tests for harmonized bulk & scRNA ``(qc=, cleanup=, normalize=)`` across modalities (#780)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oncoref import (
    barkbase,
    bgee_srp009687,
    bgee_srp181218,
    bulk_qc,
    doga,
    gse304066,
    hpa,
    human_protein_atlas,
    library_qc,
    pfizer_canine_rna,
    source_matrices,
    species_expression,
)
from oncoref.expression_contract import ExpressionContract, ExpressionMeasurement
from oncoref.gene_ids import canonical_gene_space
from oncoref.reference_identity import ReferenceIdentity
from oncoref.species_expression import ExpressionDataset


def test_bulk_processing_spec_all_alias_and_cache_key_identity():
    qc_all = bulk_qc.expand_ops(
        "all",
        allowed=library_qc.ALL_TUMOR_QC,
        aliases=library_qc.TUMOR_QC_ALIASES,
        param_name="qc",
    )
    qc_explicit = bulk_qc.expand_ops(
        list(reversed(library_qc.ALL_TUMOR_QC)),
        allowed=library_qc.ALL_TUMOR_QC,
        aliases=library_qc.TUMOR_QC_ALIASES,
        param_name="qc",
    )
    assert qc_all == qc_explicit == frozenset(library_qc.ALL_TUMOR_QC)

    cl_all = bulk_qc.expand_ops(
        "all",
        allowed=library_qc.ALL_TUMOR_CLEANUP,
        aliases=library_qc.TUMOR_CLEANUP_ALIASES,
        param_name="cleanup",
    )
    cl_explicit = bulk_qc.expand_ops(
        set(library_qc.ALL_TUMOR_CLEANUP),
        allowed=library_qc.ALL_TUMOR_CLEANUP,
        aliases=library_qc.TUMOR_CLEANUP_ALIASES,
        param_name="cleanup",
    )
    assert cl_all == cl_explicit == frozenset(library_qc.ALL_TUMOR_CLEANUP)

    spec_a = bulk_qc.BulkProcessingSpec(
        domain="tumor_bulk",
        source_id="LUAD",
        source_version="5.22.8",
        qc=qc_all,
        cleanup=cl_all,
        normalize="clean_tpm",
    )
    spec_b = bulk_qc.BulkProcessingSpec(
        domain="tumor_bulk",
        source_id="LUAD",
        source_version="5.22.8",
        qc=list(qc_explicit),  # type: ignore[arg-type]
        cleanup=tuple(cl_explicit),  # type: ignore[arg-type]
        normalize="tpm_clean",
    )
    assert spec_a == spec_b
    assert hash(spec_a) == hash(spec_b)
    assert spec_a.cache_key() == spec_b.cache_key()

    with pytest.raises(ValueError, match="Unknown"):
        bulk_qc.resolve_normalize_option("zscore")


def test_hardcoded_carryover_and_immune_ensg_ids_match_canonical_gene_space():
    cspace = canonical_gene_space()
    gid_to_sym = dict(zip(cspace["ensembl_gene_id"], cspace["symbol"]))
    expected = {
        "ENSG00000163736": "PPBP",
        "ENSG00000163737": "PF4",
        "ENSG00000169704": "GP9",
        "ENSG00000101162": "TUBB1",
        "ENSG00000005961": "ITGA2B",
        "ENSG00000185245": "GP1BA",
        "ENSG00000244734": "HBB",
        "ENSG00000206172": "HBA1",
        "ENSG00000188536": "HBA2",
        "ENSG00000223609": "HBD",
        "ENSG00000158578": "ALAS2",
        "ENSG00000004939": "SLC4A1",
    }
    for gid, sym in expected.items():
        assert gid_to_sym.get(gid) == sym


def test_hpa_normal_bulk_qc_cleanup_nan_alignment_and_tonsil_protection(monkeypatch, tmp_path):
    synthetic_hpa = pd.DataFrame(
        [
            # Liver native markers + 2% splash into placenta
            ("ENSG00000132693", "CRP", "liver", 30000.0, 32000.0, 28000.0),
            ("ENSG00000171557", "FGG", "liver", 22000.0, 24000.0, 20000.0),
            ("ENSG00000171560", "FGA", "liver", 18000.0, 19000.0, 17000.0),
            ("ENSG00000132693", "CRP", "placenta", 600.0, 620.0, 580.0),
            ("ENSG00000171557", "FGG", "placenta", 450.0, 460.0, 440.0),
            ("ENSG00000171560", "FGA", "placenta", 360.0, 370.0, 350.0),
            # Endogenous squamous epithelium in tonsil (CRNN / KRT13 / KRT10) -> MUST NOT be masked!
            ("ENSG00000143536", "CRNN", "tonsil", 222.3, 250.0, 234.4),
            ("ENSG00000171401", "KRT13", "tonsil", 800.0, 900.0, 850.0),
            ("ENSG00000186395", "KRT10", "tonsil", 12.0, 14.0, 12.9),
            # Non-coding gene with NaN pTPM/nTPM must not crash melt or vanish!
            ("ENSG00000222222", "LINC_TEST", "placenta", 5.0, np.nan, np.nan),
            ("ENSG00000222222", "LINC_TEST", "tonsil", np.nan, np.nan, np.nan),
            # Placenta & retina endogenous genes + ribosomal/biological for clean_tpm
            ("ENSG00000170965", "PLAC1", "placenta", 180.0, 190.0, 200.0),
            ("ENSG00000170965", "PLAC1", "choroid plexus", 0.0, 0.0, 0.0),
            ("ENSG00000170965", "PLAC1", "retina", 0.0, 0.0, 0.0),
            ("ENSG00000142937", "RPS8", "placenta", 4000.0, 4400.0, 4000.0),
            ("ENSG00000142937", "RPS8", "choroid plexus", 4000.0, 4400.0, 4000.0),
            ("ENSG00000142937", "RPS8", "retina", 4000.0, 4400.0, 4000.0),
            ("ENSG00000142937", "RPS8", "tonsil", 4000.0, 4400.0, 4000.0),
            ("ENSG00000142937", "RPS8", "liver", 4000.0, 4400.0, 4000.0),
            ("ENSG00000111640", "GAPDH", "placenta", 1200.0, 1320.0, 1200.0),
            ("ENSG00000111640", "GAPDH", "choroid plexus", 1200.0, 1320.0, 1200.0),
            ("ENSG00000111640", "GAPDH", "retina", 1200.0, 1320.0, 1200.0),
            ("ENSG00000111640", "GAPDH", "tonsil", 1200.0, 1320.0, 1200.0),
            ("ENSG00000111640", "GAPDH", "liver", 1200.0, 1320.0, 1200.0),
        ],
        columns=["Gene", "Gene name", "Tissue", "TPM", "pTPM", "nTPM"],
    )
    monkeypatch.setattr(hpa, "_read_hpa", lambda name, ver=None: synthetic_hpa.copy())
    monkeypatch.setattr(
        hpa.reference_data,
        "local_path",
        lambda name, ver=None: tmp_path / f"{name}-{ver}.tsv",
    )

    assert human_protein_atlas.rna_consensus is hpa.hpa_rna_consensus

    # 1. duplicate_gtex_tissues drops 'liver' while retaining 'placenta', 'choroid plexus', 'retina', 'tonsil'
    dedup = human_protein_atlas.rna_tissue_hpa(qc="duplicate_gtex_tissues", use_cache=False)
    assert "liver" not in set(dedup["Tissue"])
    assert {"placenta", "choroid plexus", "retina", "tonsil"} <= set(dedup["Tissue"])

    # 2. qc="splash" masks CRP in placenta across ALL 3 columns (TPM, pTPM, nTPM) simultaneously,
    #    while preserving endogenous CRNN in tonsil and LINC_TEST partial-NaN row!
    splashed = human_protein_atlas.rna_tissue_hpa(qc="splash", use_cache=False)
    plac_crp = splashed.loc[
        splashed["Gene name"].eq("CRP") & splashed["Tissue"].eq("placenta")
    ].iloc[0]
    assert np.isnan(plac_crp["TPM"]) and np.isnan(plac_crp["pTPM"]) and np.isnan(plac_crp["nTPM"])
    tonsil_crnn = splashed.loc[
        splashed["Gene name"].eq("CRNN") & splashed["Tissue"].eq("tonsil"), "TPM"
    ].iloc[0]
    assert tonsil_crnn == pytest.approx(222.3)
    linc_plac = splashed.loc[
        splashed["Gene name"].eq("LINC_TEST") & splashed["Tissue"].eq("placenta")
    ].iloc[0]
    assert linc_plac["TPM"] == pytest.approx(5.0) and np.isnan(linc_plac["pTPM"])

    # 3. Brain subregion resolution recognizes corpus callosum / cerebellar cortex as native 'brain'
    assert "brain" in hpa._resolve_hpa_native_panels(
        "corpus callosum, genu (white matter)",
        brain_subregions=frozenset({"corpus callosum, genu (white matter)"}),
    )
    assert "striated_muscle" not in hpa._resolve_hpa_native_panels("smooth muscle")

    # 4. nTPM-only consensus rejects normalize='clean_tpm'
    synthetic_cons = synthetic_hpa[["Gene", "Gene name", "Tissue", "nTPM"]].copy()
    monkeypatch.setattr(hpa, "_hpa_rna_consensus_for_version", lambda ver: synthetic_cons.copy())
    with pytest.raises(ValueError, match="linear 'TPM'"):
        human_protein_atlas.rna_consensus(normalize="clean_tpm", use_cache=False)

    # 5. Cache & alias equivalence
    clean_all = human_protein_atlas.rna_tissue_hpa(
        qc="duplicate_gtex_tissues", cleanup="all", normalize="clean_tpm", use_cache=True
    )
    clean_explicit = human_protein_atlas.rna_tissue_hpa(
        qc=["duplicate_gtex_tissues"],
        cleanup=["deconvolution", "index_hopping"],
        normalize="clean_tpm",
        use_cache=True,
    )
    pd.testing.assert_frame_equal(clean_all, clean_explicit)


def test_sorted_immune_bulk_and_single_cell_ambient_soup(monkeypatch, tmp_path):
    synthetic_imm = pd.DataFrame(
        [
            # Ex-vivo activated T-cell + unfractionated total PBMC + 5% platelet-contaminated NK sort -> dropped by qc="all"
            ("ENSG00000163736", "PPBP", "Naive CD4 T-cell activated", 0.0),
            ("ENSG00000163737", "PF4", "Naive CD4 T-cell activated", 0.0),
            ("ENSG00000244734", "HBB", "Naive CD4 T-cell activated", 0.0),
            ("ENSG00000206172", "HBA1", "Naive CD4 T-cell activated", 0.0),
            ("ENSG00000167286", "CD3D", "Naive CD4 T-cell activated", 950.0),
            ("ENSG00000163736", "PPBP", "total PBMC", 120.0),
            ("ENSG00000163737", "PF4", "total PBMC", 90.0),
            ("ENSG00000244734", "HBB", "total PBMC", 9000.0),
            ("ENSG00000206172", "HBA1", "total PBMC", 5600.0),
            ("ENSG00000167286", "CD3D", "total PBMC", 400.0),
            ("ENSG00000163736", "PPBP", "bad platelet-soaked NK", 500.0),
            ("ENSG00000169704", "GP9", "bad platelet-soaked NK", 25.0),
            ("ENSG00000005961", "ITGA2B", "bad platelet-soaked NK", 11.0),
            ("ENSG00000244734", "HBB", "bad platelet-soaked NK", 0.0),
            ("ENSG00000206172", "HBA1", "bad platelet-soaked NK", 0.0),
            ("ENSG00000167286", "CD3D", "bad platelet-soaked NK", 0.0),
            # FACS-sorted naive B-cell with 10% reticulocyte globin carryover -> KEPT by qc="all", CLEANED by sort_carryover!
            ("ENSG00000163736", "PPBP", "naive B-cell", 42.5),
            ("ENSG00000169704", "GP9", "naive B-cell", 2.1),
            ("ENSG00000005961", "ITGA2B", "naive B-cell", 0.9),
            ("ENSG00000244734", "HBB", "naive B-cell", 4500.0),
            ("ENSG00000188536", "HBA2", "naive B-cell", 2600.0),
            ("ENSG00000206172", "HBA1", "naive B-cell", 550.0),
            ("ENSG00000177455", "CD19", "naive B-cell", 720.0),
            ("ENSG00000167286", "CD3D", "naive B-cell", 0.0),
            # Resting naive CD4 T-cell with 0.5% platelet carryover + classical monocyte with endogenous NRGN
            ("ENSG00000163736", "PPBP", "naive CD4 T-cell", 42.5),
            ("ENSG00000169704", "GP9", "naive CD4 T-cell", 2.1),
            ("ENSG00000005961", "ITGA2B", "naive CD4 T-cell", 0.9),
            ("ENSG00000244734", "HBB", "naive CD4 T-cell", 0.0),
            ("ENSG00000206172", "HBA1", "naive CD4 T-cell", 0.0),
            ("ENSG00000177455", "CD19", "naive CD4 T-cell", 0.0),
            ("ENSG00000167286", "CD3D", "naive CD4 T-cell", 880.0),
            ("ENSG00000154146", "NRGN", "classical monocyte", 154.2),
            ("ENSG00000163736", "PPBP", "classical monocyte", 1.7),
            ("ENSG00000169704", "GP9", "classical monocyte", 0.0),
        ],
        columns=["Gene", "Gene name", "Immune cell", "TPM"],
    )
    synthetic_sc = pd.DataFrame(
        [
            # Breast c-16 alveolar milk dominator (CSN2 = 346,158.4 nTPM) + c-3 macrophage soup (683.2 nTPM)
            (
                "ENSG00000135222",
                "CSN2",
                "breast",
                "c-16",
                "breast glandular cells",
                3557723,
                346158.4,
            ),
            ("ENSG00000135222", "CSN2", "breast", "c-3", "macrophages", 6414, 683.2),
            # Misannotated breast c-24 't-cells' with HBB = 597,119.7 nTPM -> dropped by qc='misannotated_erythroid'
            ("ENSG00000244734", "HBB", "breast", "c-24", "t-cells", 900000, 597119.7),
            ("ENSG00000135222", "CSN2", "breast", "c-24", "t-cells", 500, 120.0),
            # Negative control CSN1S1 in brain excitatory neurons (16.9 nTPM, no brain dominator) -> UNTOUCHED!
            ("ENSG00000126545", "CSN1S1", "brain", "c-31", "excitatory neurons", 670, 16.9),
        ],
        columns=["Gene", "Gene name", "Tissue", "Cluster", "Cell type", "Read count", "nTPM"],
    )

    def _mock_table(src, genes=None, version=None):
        if src == "hpa_single_cell_tissue":
            return synthetic_sc.copy()
        return synthetic_imm.copy()

    monkeypatch.setattr(human_protein_atlas, "table", _mock_table)
    monkeypatch.setattr(
        human_protein_atlas.reference_data,
        "local_path",
        lambda name, ver=None: tmp_path / f"{name}-{ver}.tsv",
    )

    qc_df = human_protein_atlas.immune_cells("monaco", qc="all", use_cache=False)
    surviving_cells = set(qc_df["Immune cell"])
    # Activated T-cells & total PBMC dropped, while FACS-sorted naive B-cell is RETAINED!
    assert surviving_cells == {"naive B-cell", "naive CD4 T-cell", "classical monocyte"}
    assert qc_df.loc[
        qc_df["Gene name"].eq("HBB") & qc_df["Immune cell"].eq("naive B-cell"), "TPM"
    ].iloc[0] == pytest.approx(4500.0)

    cleaned = human_protein_atlas.immune_cells(
        "monaco",
        qc="all",
        cleanup="all",
        genes=["PPBP", "HBB", "CD19", "CD3D", "NRGN"],
        use_cache=False,
    )
    bcell = cleaned[cleaned["Immune cell"].eq("naive B-cell")].set_index("Gene name")["TPM"]
    tcell = cleaned[cleaned["Immune cell"].eq("naive CD4 T-cell")].set_index("Gene name")["TPM"]
    mono = cleaned[cleaned["Immune cell"].eq("classical monocyte")].set_index("Gene name")["TPM"]
    assert bcell["HBB"] < 1.0
    assert bcell["PPBP"] < 1.0
    assert bcell["CD19"] == pytest.approx(720.0)
    assert tcell["PPBP"] < 1.0
    assert tcell["CD3D"] == pytest.approx(880.0)
    assert mono["NRGN"] == pytest.approx(154.2)

    # Single-cell tissue QC & ambient soup subtraction
    sc_clean = human_protein_atlas.single_cell_tissues(qc="all", cleanup="all", use_cache=False)
    assert "c-24" not in set(sc_clean["Cluster"])  # misannotated erythroid cluster dropped
    c16_csn2 = sc_clean.loc[sc_clean["Cluster"].eq("c-16"), "nTPM"].iloc[0]
    c3_csn2 = sc_clean.loc[sc_clean["Cluster"].eq("c-3"), "nTPM"].iloc[0]
    brain_csn1s1 = sc_clean.loc[sc_clean["Gene name"].eq("CSN1S1"), "nTPM"].iloc[0]
    assert c16_csn2 == pytest.approx(346158.4)
    assert c3_csn2 == pytest.approx(0.0)
    assert brain_csn1s1 == pytest.approx(16.9)
    with pytest.raises(ValueError, match="linear 'TPM'"):
        human_protein_atlas.single_cell_tissues(normalize="clean_tpm", use_cache=False)


@pytest.mark.parametrize(
    "taxon_id,or_prefix",
    [
        (9615, "OR1A"),  # Canine (Canis lupus familiaris)
        (10090, "Olfr"),  # Mouse (Mus musculus)
    ],
)
def test_nonhuman_species_expression_dataset_process_and_contract(taxon_id, or_prefix, tmp_path):
    # Verify canine tissue substring disambiguation ('adrenal' != kidney, 'kidneycortex' != brain, 'smooth muscle' != striated)
    assert species_expression._tissue_to_species_panels("adrenal") == frozenset({"adrenal_gland"})
    assert species_expression._tissue_to_species_panels("kidneycortex") == frozenset({"kidney"})
    assert species_expression._tissue_to_species_panels("smooth muscle") == frozenset()

    ref = ReferenceIdentity(
        taxon_id=taxon_id,
        assembly_accession=f"GCF_{taxon_id}.1",
        annotation_source="Ensembl",
        annotation_release="112",
        source_version="v1",
        identifier_namespaces={"gene": "Ensembl"},
        content_provenance={"gtf": "a" * 64},
    )
    contract = ExpressionContract(
        taxon_id=taxon_id,
        source_id=f"atlas_{taxon_id}",
        source_version="v1",
        context="normal",
        reference=ref,
        reference_status="verified",
        feature_namespace="Ensembl",
        measurement=ExpressionMeasurement("gene", "TPM", "linear", False),
        content_provenance={"matrix": "b" * 64},
        provenance_notes="Synthetic non-human multi-organ bulk fixture",
    )

    or_genes = [f"FEAT_OR_{i}" for i in range(25)]
    or_syms = [f"{or_prefix}{i + 1}" for i in range(25)]
    base_feats = [
        ("F_ALB", "ALB"),
        ("F_FGG", "FGG"),
        ("F_FGA", "FGA"),
        ("F_PRM1", "PRM1"),
        ("F_PRM2", "PRM2"),
        ("F_TNP1", "TNP1"),
        ("F_GFAP", "GFAP"),
        ("F_MBP", "MBP"),
        ("F_RPL8", "RPL8"),
        ("F_MT1", "MT-CO1"),
        ("F_ACTB", "ACTB"),
        *zip(or_genes, or_syms),
    ]
    feat_ids = [f for f, _ in base_feats]
    features = pd.DataFrame(
        {"source_feature_id": feat_ids, "GeneSymbol": [s for _, s in base_feats]},
        index=pd.Index(feat_ids, name="feature_id"),
    )

    samples = pd.DataFrame(
        [
            {
                "sample_id": "liver_1",
                "taxon_id": taxon_id,
                "library_id": "L1",
                "run_ids": ["R1"],
                "donor_id": "D1",
                "donor_scope": "S",
                "pooling": "single",
                "tissue": "liver",
            },
            {
                "sample_id": "liver_2",
                "taxon_id": taxon_id,
                "library_id": "L2",
                "run_ids": ["R2"],
                "donor_id": "D2",
                "donor_scope": "S",
                "pooling": "single",
                "tissue": "liver",
            },
            {
                "sample_id": "testis_1",
                "taxon_id": taxon_id,
                "library_id": "T1",
                "run_ids": ["R3"],
                "donor_id": "D1",
                "donor_scope": "S",
                "pooling": "single",
                "tissue": "testis",
            },
            {
                "sample_id": "brain_splash",
                "taxon_id": taxon_id,
                "library_id": "B1",
                "run_ids": ["R4"],
                "donor_id": "D1",
                "donor_scope": "S",
                "pooling": "single",
                "tissue": "brain",
            },
            {
                "sample_id": "brain_hop",
                "taxon_id": taxon_id,
                "library_id": "B2",
                "run_ids": ["R5"],
                "donor_id": "D2",
                "donor_scope": "S",
                "pooling": "single",
                "tissue": "brain",
            },
            {
                "sample_id": "brain_swap",
                "taxon_id": taxon_id,
                "library_id": "B3",
                "run_ids": ["R6"],
                "donor_id": "D3",
                "donor_scope": "S",
                "pooling": "single",
                "tissue": "brain",
            },
            {
                "sample_id": "pooled_lib",
                "taxon_id": taxon_id,
                "library_id": "P1",
                "run_ids": ["R7"],
                "donor_id": None,
                "donor_scope": "S",
                "pooling": "pooled",
                "tissue": "brain",
            },
            {
                "sample_id": "gdna_lib",
                "taxon_id": taxon_id,
                "library_id": "G1",
                "run_ids": ["R8"],
                "donor_id": "D4",
                "donor_scope": "S",
                "pooling": "single",
                "tissue": "brain",
            },
        ]
    ).set_index("sample_id")

    mat = pd.DataFrame(0.0, index=features.index, columns=samples.index)
    for col in ("liver_1", "liver_2"):
        mat.loc[["F_ALB", "F_FGG", "F_FGA"], col] = [40000.0, 20000.0, 15000.0]
    mat.loc[["F_PRM1", "F_PRM2", "F_TNP1"], "testis_1"] = [25000.0, 20000.0, 10000.0]
    for col in ("brain_splash", "brain_hop", "pooled_lib", "gdna_lib"):
        mat.loc[["F_GFAP", "F_MBP"], col] = [8000.0, 6000.0]
    mat.loc[["F_PRM1", "F_PRM2", "F_TNP1"], "brain_splash"] = [500.0, 400.0, 200.0]
    mat.loc[["F_ALB", "F_FGG", "F_FGA"], "brain_hop"] = [40.0, 20.0, 15.0]
    mat.loc[["F_ALB", "F_FGG", "F_FGA"], "brain_swap"] = [32000.0, 16000.0, 12000.0]
    mat.loc[or_genes[:20], "gdna_lib"] = 8.0
    mat.loc["F_RPL8", :] = 5000.0
    mat.loc["F_MT1", :] = 3000.0
    mat.loc["F_ACTB", :] = 2000.0

    ds = ExpressionDataset(mat, features, samples, contract)
    qc_ds = ds.process(qc="all")
    assert set(qc_ds.samples.index) == {
        "liver_1",
        "liver_2",
        "testis_1",
        "brain_splash",
        "brain_hop",
    }
    assert np.isnan(qc_ds.values.loc["F_PRM1", "brain_splash"])
    assert qc_ds.values.loc["F_GFAP", "brain_splash"] == pytest.approx(8000.0)
    assert qc_ds.values.loc["F_PRM1", "testis_1"] == pytest.approx(25000.0)

    full_ds = ds.process(qc="all", cleanup="all", normalize="clean_tpm")
    explicit_ds = ds.process(
        qc=("splash", "swaps", "broad_background", "pooled_libraries"),
        cleanup=("deconvolution", "index_hopping"),
        normalize="clean_tpm",
    )
    pd.testing.assert_frame_equal(full_ds.values, explicit_ds.values)
    assert full_ds.contract == explicit_ds.contract
    assert full_ds.contract.measurement.unit == "clean-TPM"
    assert full_ds.values.loc["F_ALB", "brain_hop"] < 1e-6
    assert full_ds.values["brain_hop"].sum() == pytest.approx(1_000_000.0)

    # Empty sample survival save/load & tissue_summary schema roundtrip must not crash
    only_pooled = ExpressionDataset(
        mat[["pooled_lib"]], features, samples.loc[["pooled_lib"]], contract
    )
    empty_after_qc = only_pooled.process(qc="pooled_libraries")
    assert empty_after_qc.values.shape[1] == 0
    empty_sum = empty_after_qc.tissue_summary()
    assert list(empty_sum.columns) == [
        "feature_id",
        "tissue",
        "native_mean",
        "native_max",
        "measured_libraries",
        "total_libraries",
        "known_individual_donors_measured",
        "libraries_with_unresolved_or_pooled_donors",
        "unit",
    ]
    out_dir = tmp_path / f"empty_ds_{taxon_id}"
    empty_after_qc.save(out_dir)
    reloaded = ExpressionDataset.load(out_dir)
    assert reloaded.values.shape == (len(features), 0)

    # Protein-group iBAQ dataset must reject RNA bulk cleanup
    prot_contract = ExpressionContract(
        taxon_id=taxon_id,
        source_id=f"prot_{taxon_id}",
        source_version="v1",
        context="normal",
        reference=ref,
        reference_status="verified",
        feature_namespace="Ensembl",
        measurement=ExpressionMeasurement("protein_group", "iBAQ", "linear", False),
        content_provenance={"matrix": "b" * 64},
        provenance_notes="Synthetic proteomics fixture",
    )
    prot_ds = ExpressionDataset(mat, features, samples, prot_contract)
    with pytest.raises(ValueError, match="gene/transcript"):
        prot_ds.process(cleanup="all")


def test_named_canine_modules_accept_harmonized_kwargs_and_validate():
    for mod in (barkbase, bgee_srp009687, bgee_srp181218, pfizer_canine_rna, doga, gse304066):
        with pytest.raises(ValueError, match="Unknown species bulk qc"):
            mod.read(qc="not_a_real_qc_step")
        with pytest.raises(ValueError, match="Unknown species bulk cleanup"):
            mod.read(cleanup="not_a_real_cleanup_step")


def test_tumor_bulk_qc_cleanup_and_strict_reproductive_cta_protection(monkeypatch, tmp_path):
    """Tumor cleanup must subtract somatic bleed/stroma while NEVER subtracting CTAs, CGA, MAGEA3, or oncofetal genes."""
    or_rows = [(f"ENSG0000099{i:04d}", f"OR1A{i}") for i in range(1, 115)]
    marker_rows = [
        # Testis / germline markers (Note R) + CTAs + unsafe CTAs (MAGEA3) + trophoblast/oncofetal (CGA, AFP, CLDN6)
        ("ENSG00000175646", "PRM1"),
        ("ENSG00000122304", "PRM2"),
        ("ENSG00000118245", "TNP1"),
        ("ENSG00000010318", "PHF7"),
        ("ENSG00000203784", "LELP1"),
        ("ENSG00000163206", "SMCP"),
        ("ENSG00000147381", "MAGEA4"),
        ("ENSG00000221867", "MAGEA3"),
        ("ENSG00000184033", "CTAG1B"),
        ("ENSG00000185686", "PRAME"),
        ("ENSG00000170965", "PLAC1"),
        ("ENSG00000135346", "CGA"),
        ("ENSG00000081051", "AFP"),
        ("ENSG00000184697", "CLDN6"),
        # Lung native markers (cohort is LUAD)
        ("ENSG00000185303", "SFTPA2"),
        ("ENSG00000122852", "SFTPA1"),
        ("ENSG00000168484", "SFTPC"),
        # Foreign somatic panels: liver (including ALB from gtex-deconvolution-anchors), pancreas, thyroid, kidney
        ("ENSG00000163631", "ALB"),
        ("ENSG00000132693", "CRP"),
        ("ENSG00000257017", "HP"),
        ("ENSG00000229314", "ORM1"),
        ("ENSG00000171557", "FGG"),
        ("ENSG00000171560", "FGA"),
        ("ENSG00000171564", "FGB"),
        ("ENSG00000254647", "INS"),
        ("ENSG00000275896", "PRSS2"),
        ("ENSG00000204983", "PRSS1"),
        ("ENSG00000091704", "CPA1"),
        ("ENSG00000137392", "CLPS"),
        ("ENSG00000142789", "CELA3A"),
        ("ENSG00000042832", "TG"),
        ("ENSG00000009765", "IYD"),
        ("ENSG00000165409", "TSHR"),
        ("ENSG00000169344", "UMOD"),
        ("ENSG00000167580", "AQP2"),
        ("ENSG00000074803", "SLC12A1"),
        ("ENSG00000149452", "SLC22A8"),
        # Clean-TPM compartments
        ("ENSG00000142937", "RPS8"),
        ("ENSG00000198804", "MT-CO1"),
        ("ENSG00000111640", "GAPDH"),
    ]
    all_genes = [*marker_rows, *or_rows]
    # Include 6 clean LUAD tumors so cohort median for foreign liver/pancreas/thyroid/kidney is strictly 0.0!
    clean_luad_cols = [f"TUMOR_CLEAN_{i}" for i in range(1, 7)]
    sample_names = [
        "TUMOR_CTA_HIGH",
        "TUMOR_NOTE_R_MET",
        *clean_luad_cols,
        "TUMOR_GDNA_FAIL",
        "TUMOR_MULTI_TISSUE",
        "TUMOR_LIVER_SWAP",
    ]
    df = pd.DataFrame(
        {
            "Ensembl_Gene_ID": [g for g, _ in all_genes],
            "Symbol": [s for _, s in all_genes],
            **{s: np.zeros(len(all_genes), dtype=float) for s in sample_names},
        }
    ).set_index("Symbol", drop=False)

    for s in ("TUMOR_CTA_HIGH", *clean_luad_cols, "TUMOR_GDNA_FAIL", "TUMOR_MULTI_TISSUE"):
        df.loc[["SFTPA2", "SFTPA1", "SFTPC"], s] = [600.0, 500.0, 450.0]
        df.loc[["RPS8", "MT-CO1", "GAPDH"], s] = [4000.0, 2500.0, 1500.0]

    df.loc[["PRM1", "PRM2", "TNP1", "PHF7", "LELP1", "SMCP"], "TUMOR_CTA_HIGH"] = [
        120.0,
        110.0,
        55.0,
        25.0,
        18.0,
        15.0,
    ]
    df.loc[
        ["MAGEA4", "MAGEA3", "CTAG1B", "PRAME", "PLAC1", "CGA", "AFP", "CLDN6"], "TUMOR_CTA_HIGH"
    ] = [
        250.0,
        190.0,
        180.0,
        310.0,
        95.0,
        140.0,
        85.0,
        65.0,
    ]
    df.loc[["ALB", "CRP", "HP", "ORM1", "FGG", "FGA", "FGB"], "TUMOR_CTA_HIGH"] = [
        450.0,
        332.0,
        322.0,
        240.0,
        232.0,
        174.0,
        167.0,
    ]
    # Note R sample resected from liver metastasis (high liver stroma + no native lung surfactant)
    # MUST NOT be dropped by normal_swaps because Note R never rejects a tumor library!
    df.loc[["PRM1", "PRM2", "TNP1", "PHF7", "LELP1", "SMCP"], "TUMOR_NOTE_R_MET"] = [
        120.0,
        110.0,
        55.0,
        25.0,
        18.0,
        15.0,
    ]
    df.loc[["ALB", "CRP", "HP", "ORM1", "FGG", "FGA", "FGB"], "TUMOR_NOTE_R_MET"] = [
        15000.0,
        18000.0,
        16000.0,
        12000.0,
        11000.0,
        9000.0,
        8500.0,
    ]
    df.loc[["RPS8", "MT-CO1", "GAPDH"], "TUMOR_NOTE_R_MET"] = [4000.0, 2500.0, 1500.0]

    df.loc[[s for _, s in or_rows[:40]], "TUMOR_GDNA_FAIL"] = 9.0
    df.loc[["CRP", "HP", "ORM1", "FGG", "FGA", "FGB"], "TUMOR_MULTI_TISSUE"] = 400.0
    df.loc[["INS", "PRSS2", "PRSS1", "CPA1", "CLPS", "CELA3A"], "TUMOR_MULTI_TISSUE"] = 1200.0
    df.loc[["TG", "IYD", "TSHR"], "TUMOR_MULTI_TISSUE"] = 40.0
    df.loc[["UMOD", "AQP2", "SLC12A1", "SLC22A8"], "TUMOR_MULTI_TISSUE"] = 10.0
    df.loc[["CRP", "HP", "ORM1", "FGG", "FGA", "FGB"], "TUMOR_LIVER_SWAP"] = [
        18000.0,
        16000.0,
        12000.0,
        11000.0,
        9000.0,
        8500.0,
    ]
    df.loc[["RPS8", "MT-CO1", "GAPDH"], "TUMOR_LIVER_SWAP"] = [4000.0, 2500.0, 1500.0]
    df = df.reset_index(drop=True)

    pq_file = tmp_path / "v5.22.8" / "LUAD.parquet"
    pq_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(pq_file, index=False)
    monkeypatch.setenv("CANCERDATA_SOURCE_MATRICES", str(tmp_path))
    monkeypatch.setattr(source_matrices, "resolve_version", lambda code, version=None: "5.22.8")

    sig = library_qc.signatures(df, sample_names)
    assert library_qc.notes(sig.loc["TUMOR_CTA_HIGH"]) == ["R"]
    assert library_qc.notes(sig.loc["TUMOR_NOTE_R_MET"]) == ["R"]

    cleaned = source_matrices.read(
        "LUAD", auto_fetch=False, qc="all", cleanup="all", normalize=None
    )
    kept_cols = [c for c in cleaned.columns if c not in {"Ensembl_Gene_ID", "Symbol"}]
    assert kept_cols == ["TUMOR_CTA_HIGH", "TUMOR_NOTE_R_MET", *clean_luad_cols]
    assert set(cleaned.attrs["bulk_spec"]["dropped_samples"]) == {
        "TUMOR_GDNA_FAIL",
        "TUMOR_MULTI_TISSUE",
        "TUMOR_LIVER_SWAP",
    }

    by_sym = cleaned.set_index("Symbol")
    for protected_sym, expected_tpm in [
        ("PRM1", 120.0),
        ("PRM2", 110.0),
        ("TNP1", 55.0),
        ("MAGEA4", 250.0),
        ("MAGEA3", 190.0),
        ("CTAG1B", 180.0),
        ("PRAME", 310.0),
        ("PLAC1", 95.0),
        ("CGA", 140.0),
        ("AFP", 85.0),
        ("CLDN6", 65.0),
        ("SFTPA2", 600.0),
    ]:
        assert by_sym.loc[protected_sym, "TUMOR_CTA_HIGH"] == pytest.approx(expected_tpm), (
            protected_sym
        )

    assert by_sym.loc["ALB", "TUMOR_CTA_HIGH"] < 350.0
    assert by_sym.loc["CRP", "TUMOR_CTA_HIGH"] < 200.0
    assert by_sym.loc["FGG", "TUMOR_CTA_HIGH"] < 150.0

    # Verify PAAD native pancreas protection even if >50% of PDAC samples lack islet INS
    paad_df = (
        df[["Ensembl_Gene_ID", "Symbol", *clean_luad_cols]].copy().set_index("Symbol", drop=False)
    )
    paad_df.loc[["SFTPA2", "SFTPA1", "SFTPC"], clean_luad_cols] = 0.0
    paad_df.loc[["INS", "PRSS1", "CPA1"], "TUMOR_CLEAN_1"] = [80000.0, 95000.0, 40000.0]
    paad_proc = library_qc.apply_tumor_bulk_spec(
        paad_df.reset_index(drop=True), cancer_code="PAAD", qc="all", cleanup="all"
    )
    paad_sym = paad_proc.set_index("Symbol")
    assert paad_sym.loc["PRSS1", "TUMOR_CLEAN_1"] == pytest.approx(95000.0)
    assert paad_sym.loc["INS", "TUMOR_CLEAN_1"] == pytest.approx(80000.0)

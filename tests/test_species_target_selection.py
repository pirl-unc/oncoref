"""Cross-species target-selection, comparative HostContext and peptide subtraction."""

from __future__ import annotations

import hashlib
import json
from types import SimpleNamespace

import pandas as pd
import pytest

from oncoref import canine, cli, host_context, human, mouse, target_selection
from oncoref.expression_contract import ExpressionContract, ExpressionMeasurement
from oncoref.reference_identity import ReferenceIdentity
from oncoref.species_expression import ExpressionDataset
from oncoref.species_target_selection import (
    SpeciesScreeningSpec,
    anatomy_for_species,
    evaluate_species_catalog,
    group_proteome,
    group_reference_proteins,
    observations_from_expression_dataset,
)


def _make_seq(tag: str, length: int = 40) -> str:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    digest = hashlib.sha256(tag.encode("utf-8")).digest()
    return "".join(
        alphabet[(digest[i % len(digest)] + 7 * i) % len(alphabet)] for i in range(length)
    )


def test_comparative_anatomy_profiles_and_host_contexts():
    h = anatomy_for_species("human")
    d = anatomy_for_species("dog")
    m = anatomy_for_species("mouse")

    assert "seminal vesicle" in h.male_accessory_compartments
    assert "seminal vesicle" not in d.male_accessory_compartments
    assert "coagulating_gland" in m.male_accessory_compartments
    assert "coagulating_gland" not in h.male_accessory_compartments
    assert "coagulating_gland" not in d.male_accessory_compartments

    assert "uterine_horn" in d.female_tract_compartments
    assert "uterine_horn" in m.female_tract_compartments
    assert "uterine_horn" not in h.female_tract_compartments

    assert "appendix" in h.disposable_compartments
    assert "appendix" not in d.disposable_compartments
    assert "appendix" not in m.disposable_compartments
    assert "tonsil" in d.disposable_compartments
    assert "tonsil" not in m.disposable_compartments

    with pytest.raises(ValueError, match="Unknown compartment 'seminal vesicle' for species 'dog'"):
        host_context.HostContext(protect=("seminal vesicle",), species="dog")

    dog_host = host_context.HostContext(protect=("mammary", "uterus"), species="dog")
    assert dog_host.protected == (
        "uterus",
        "uterine_horn",
        "endometrium",
        "lactating_mammary",
        "resting_mammary",
    )

    male_dog = host_context.from_legacy("male-patient", species="dog")
    assert "uterine_horn" in male_dog.absent
    assert "lactating_mammary" in male_dog.absent
    assert "seminal vesicle" not in male_dog.states()

    female_mouse = host_context.from_legacy("female-patient", lactating=True, species="mouse")
    assert "coagulating_gland" in female_mouse.absent
    assert "seminal vesicle" in female_mouse.absent
    assert "lactating_mammary" in female_mouse.protected


def test_end_to_end_canine_pipeline_from_expression_datasets_to_peptide_subtraction():
    dog_ref = ReferenceIdentity.from_dict(canine.references.registered_references()[0]["identity"])
    unsafe_seq = "M" + "ACDEFGHIKLMNPQRSTVWY" * 2
    contained_seq = unsafe_seq[:39]
    cta_seq = _make_seq("DOG_CTA", 45)
    ca_seq = _make_seq("DOG_CA", 45)
    ra_seq = _make_seq("DOG_RA", 45)
    bone_seq = _make_seq("DOG_SP7", 45)
    low_seq = _make_seq("DOG_LOW", 45)

    proteins, seqs = group_reference_proteins(
        [
            {
                "gene_id": "400001",
                "symbol": "CTAG_A",
                "transcript_id": "XM_1.1",
                "protein_id": "XP_1.1",
                "sequence": cta_seq,
            },
            {
                "gene_id": "400002",
                "symbol": "CTAG_B",
                "transcript_id": "XM_2.1",
                "protein_id": "XP_2.1",
                "sequence": cta_seq,
            },
            {
                "gene_id": "400003",
                "symbol": "DOG_CA",
                "transcript_id": "XM_3.1",
                "protein_id": "XP_3.1",
                "sequence": ca_seq,
            },
            {
                "gene_id": "400004",
                "symbol": "DOG_RA",
                "transcript_id": "XM_4.1",
                "protein_id": "XP_4.1",
                "sequence": ra_seq,
            },
            {
                "gene_id": "400005",
                "symbol": "DOG_UNSAFE",
                "transcript_id": "XM_5.1",
                "protein_id": "XP_5.1",
                "sequence": unsafe_seq,
            },
            {
                "gene_id": "400006",
                "symbol": "DOG_CONTAINED",
                "transcript_id": "XM_6.1",
                "protein_id": "XP_6.1",
                "sequence": contained_seq,
            },
            {
                "gene_id": "400007",
                "symbol": "DOG_BONE",
                "transcript_id": "XM_7.1",
                "protein_id": "XP_7.1",
                "sequence": bone_seq,
            },
            {
                "gene_id": "400008",
                "symbol": "DOG_LOW",
                "transcript_id": "XM_8.1",
                "protein_id": "XP_8.1",
                "sequence": low_seq,
            },
            {
                "gene_id": "400009",
                "symbol": "DOG_UNTRANS",
                "transcript_id": "XR_9.1",
                "protein_id": "",
                "sequence": "",
                "is_qualified": False,
            },
        ],
        species="dog",
        reference=dog_ref,
    )
    assert len(proteins) == 8
    proteins.loc[proteins.protein_group_label.eq("DOG_BONE"), "home_tissue"] = "bone"

    gene_ids = [f"40000{i}" for i in range(1, 10)]
    normal_samples = pd.DataFrame(
        [
            {
                "taxon_id": 9615,
                "library_id": "L_liv1",
                "run_ids": ["R1"],
                "donor_id": "D1",
                "donor_scope": "study",
                "pooling": "single",
                "tissue": "liver",
            },
            {
                "taxon_id": 9615,
                "library_id": "L_hrt1",
                "run_ids": ["R2"],
                "donor_id": "D1",
                "donor_scope": "study",
                "pooling": "single",
                "tissue": "heart",
            },
            {
                "taxon_id": 9615,
                "library_id": "L_thy1",
                "run_ids": ["R3"],
                "donor_id": "D1",
                "donor_scope": "study",
                "pooling": "single",
                "tissue": "thymus",
            },
            {
                "taxon_id": 9615,
                "library_id": "L_gal1",
                "run_ids": ["R4"],
                "donor_id": "D1",
                "donor_scope": "study",
                "pooling": "single",
                "tissue": "gallbladder",
            },
            {
                "taxon_id": 9615,
                "library_id": "L_tes1",
                "run_ids": ["R5"],
                "donor_id": "D2",
                "donor_scope": "study",
                "pooling": "single",
                "tissue": "testis",
            },
            {
                "taxon_id": 9615,
                "library_id": "L_uth1",
                "run_ids": ["R6"],
                "donor_id": "D3",
                "donor_scope": "study",
                "pooling": "single",
                "tissue": "uterine horn",
            },
        ],
        index=["L_liv1", "L_hrt1", "L_thy1", "L_gal1", "L_tes1", "L_uth1"],
    )
    normal_vals = pd.DataFrame(
        {
            "L_liv1": [0.2, 0.2, 0.5, 0.3, 80.0, 0.1, 0.2, 0.5, 0.1],
            "L_hrt1": [0.1, 0.1, 0.4, 0.2, 90.0, 0.1, 0.1, 0.4, 0.1],
            "L_thy1": [25.0, 25.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0],
            "L_gal1": [1.0, 0.5, 0.2, 0.2, 20.0, 0.1, 0.2, 0.2, 0.1],
            "L_tes1": [25.0, 15.0, 1.0, 2.0, 30.0, 25.0, 1.0, 2.0, 20.0],
            "L_uth1": [0.0, 0.0, 0.5, 35.0, 15.0, 0.0, 0.0, 1.0, 0.0],
        },
        index=gene_ids,
    )
    normal_contract = ExpressionContract(
        9615,
        "synthetic-barkbase",
        "1.0",
        "normal",
        dog_ref,
        "verified",
        "NCBI Gene",
        ExpressionMeasurement("gene", "TMM-CPM", "linear", False),
        {"matrix": "a" * 64},
        "Synthetic canine normal TMM-CPM matrix.",
    )
    normal_ds = ExpressionDataset(
        normal_vals,
        pd.DataFrame({"symbol": gene_ids}, index=gene_ids),
        normal_samples,
        normal_contract,
    )

    tumor_samples = pd.DataFrame(
        [
            {
                "taxon_id": 9615,
                "library_id": "T_osa1",
                "run_ids": ["RT1"],
                "donor_id": "P1",
                "donor_scope": "study",
                "pooling": "single",
                "tissue": "osteosarcoma",
                "cancer_code": "OSA",
            },
            {
                "taxon_id": 9615,
                "library_id": "T_osa_pool",
                "run_ids": ["RT2"],
                "donor_id": "POOL",
                "donor_scope": "study",
                "pooling": "pooled",
                "tissue": "osteosarcoma",
                "cancer_code": "OSA",
            },
        ],
        index=["T_osa1", "T_osa_pool"],
    )
    tumor_vals = pd.DataFrame(
        {
            "T_osa1": [12.0, 8.0, 25.0, 1.0, 50.0, 20.0, 45.0, 3.0, 30.0],
            "T_osa_pool": [15.0, 10.0, 30.0, 0.5, 40.0, 18.0, 50.0, 2.0, 25.0],
        },
        index=gene_ids,
    )
    tumor_contract = ExpressionContract(
        9615,
        "synthetic-canine-osa",
        "1.0",
        "tumor",
        dog_ref,
        "verified",
        "NCBI Gene",
        ExpressionMeasurement("gene", "TMM-CPM", "linear", False),
        {"matrix": "b" * 64},
        "Synthetic canine osteosarcoma TMM-CPM matrix.",
    )
    tumor_ds = ExpressionDataset(
        tumor_vals,
        pd.DataFrame({"symbol": gene_ids}, index=gene_ids),
        tumor_samples,
        tumor_contract,
    )

    norm_obs = observations_from_expression_dataset(
        normal_ds, proteins, comparison_id="canine-tmmcpm-v1"
    )
    tum_obs = observations_from_expression_dataset(
        tumor_ds, proteins, comparison_id="canine-tmmcpm-v1"
    )

    unsafe_eid = proteins.loc[proteins.protein_group_label.eq("DOG_UNSAFE"), "entity_id"].iloc[0]
    gates = pd.DataFrame(
        [
            {
                "entity_id": unsafe_eid,
                "layer": "Canine normal bulk RNA",
                "unsafe": True,
                "hold": False,
                "reason": "heart 90 TMM-CPM",
            }
        ]
    )

    spec = SpeciesScreeningSpec(
        species="dog",
        reference=dog_ref,
        comparison_id="canine-tmmcpm-v1",
        comparability_reference="synthetic-canine-receipt",
        comparison_unit="TMM-CPM",
        required_normal_contexts=frozenset({"liver", "heart"}),
    )
    res = canine.evaluate_target_catalog(
        spec,
        proteins,
        sequences=seqs,
        normal_gates=gates,
        comparison_normals=norm_obs["normals"],
        comparison_reproductive=norm_obs["reproductive"],
        comparison_tumors=tum_obs["tumors"],
    )
    by_label = res.catalog.set_index("protein_group_label")
    assert by_label.loc["CTAG_A/CTAG_B", "classification"] == "CTA"
    assert by_label.loc["CTAG_A/CTAG_B", "cta_tier"] == "core"
    assert by_label.loc["CTAG_A/CTAG_B", "samples_at_tenfold"] == 2
    assert by_label.loc["CTAG_A/CTAG_B", "independent_patients_at_tenfold"] == 1
    assert by_label.loc["DOG_CA", "classification"] == "Cancer Antigen only"
    assert by_label.loc["DOG_RA", "classification"] == "Reproductive Antigen only"
    assert by_label.loc["DOG_RA", "reproductive_qualifying_compartments"] == "uterine_horn"
    assert by_label.loc["DOG_UNSAFE", "classification"] == "Unsafe"
    assert by_label.loc["DOG_CONTAINED", "classification"] == "Other:Sequence"
    assert by_label.loc["DOG_CONTAINED", "unsafe_sequence_overlap_pct"] >= 95.0
    assert by_label.loc["DOG_BONE", "classification"] == "Other:Unscreened tissue"
    assert by_label.loc["DOG_LOW", "classification"] == "Other:Expression"
    assert by_label.loc["DOG_UNTRANS", "classification"] == "Other:Review"

    sel_all = res.select("CTA")
    assert sel_all.proteins().protein_group_label.tolist() == ["CTAG_A/CTAG_B"]
    assert (sel_all.decisions().papers_counted == 0).all()
    sel_protect_testis = res.select(
        "CTA", host=host_context.HostContext(protect=("testis",), species="dog")
    )
    assert sel_protect_testis.proteins().empty
    bg_index = sel_protect_testis.background_index("Unsafe")
    assert bg_index.subtract(cta_seq, lengths=[9]).empty
    assert not sel_all.background_index("Unsafe").subtract(cta_seq, lengths=[9]).empty


def test_mouse_rare_cell_umi_proteome_adapter_and_species_boundary_enforcement():
    mouse_ref = ReferenceIdentity.from_dict(mouse.references.registered_references()[0]["identity"])
    fake_proteome = [
        SimpleNamespace(
            gene_id="ENSMUSG00000051257",
            gene_symbol="Trap1a",
            transcript_id="ENSMUST1.1",
            accession="ENSMUSP1.1",
            biotype="protein_coding",
            sequence=_make_seq("TRAP1A", 50),
            kind="annotated",
        ),
        SimpleNamespace(
            gene_id="ENSMUSG00000099991",
            gene_symbol="Magea_clean",
            transcript_id="ENSMUST2.1",
            accession="ENSMUSP2.1",
            biotype="protein_coding",
            sequence=_make_seq("MAGEA_CLEAN", 50),
            kind="annotated",
        ),
    ]
    proteins, seqs = group_proteome(fake_proteome, species="mouse", reference=mouse_ref)
    umi = pd.DataFrame(
        [
            {
                "gene_id": "ENSMUSG00000051257",
                "source_id": "tabula_muris_10x",
                "donor_id": "M1",
                "tissue": "lung",
                "cell_type": "type II pneumocyte",
                "tier": "vital",
                "n_cells": 2000,
                "n_positive": 80,
                "n_multi": 45,
            },
            {
                "gene_id": "ENSMUSG00000051257",
                "source_id": "mouse_cell_atlas",
                "donor_id": "M2",
                "tissue": "lung",
                "cell_type": "type II pneumocyte",
                "tier": "vital",
                "n_cells": 2000,
                "n_positive": 70,
                "n_multi": 40,
            },
            {
                "gene_id": "ENSMUSG00000099991",
                "source_id": "tabula_muris_10x",
                "donor_id": "M1",
                "tissue": "lung",
                "cell_type": "type II pneumocyte",
                "tier": "vital",
                "n_cells": 2000,
                "n_positive": 0,
                "n_multi": 0,
            },
        ]
    )
    eid_clean = proteins.loc[proteins.protein_group_label.eq("Magea_clean"), "entity_id"].iloc[0]
    eid_trap = proteins.loc[proteins.protein_group_label.eq("Trap1a"), "entity_id"].iloc[0]
    normals = pd.DataFrame(
        [
            {
                "entity_id": eid,
                "source_id": "mouse-norm",
                "context_id": ctx,
                "tissue": ctx,
                "scope": "normal",
                "value": 0.2,
                "comparison_id": "mouse-tpm-v1",
                "unit": "TPM",
                "qc_pass": True,
                "identity_resolved": True,
            }
            for eid in (eid_clean, eid_trap)
            for ctx in ("lung", "liver")
        ]
    )
    repro = pd.DataFrame(
        [
            {
                "entity_id": eid_clean,
                "source_id": "mouse-norm",
                "sample_id": "S_cg1",
                "participant_id": "M3",
                "scope": "reproductive",
                "compartment": "coagulating_gland",
                "value": 30.0,
                "comparison_id": "mouse-tpm-v1",
                "unit": "TPM",
                "qc_pass": True,
                "identity_resolved": True,
            },
        ]
    )
    tumors = pd.DataFrame(
        [
            {
                "entity_id": eid_clean,
                "source_id": "mouse-tum",
                "sample_id": "P815_1",
                "participant_id": "DBA2_1",
                "value": 45.0,
                "comparison_id": "mouse-tpm-v1",
                "unit": "TPM",
                "qc_pass": True,
                "identity_resolved": True,
            },
        ]
    )
    spec = SpeciesScreeningSpec(
        species="mouse",
        reference=mouse_ref,
        comparison_id="mouse-tpm-v1",
        comparability_reference="synthetic-mouse-receipt",
        comparison_unit="TPM",
        required_normal_contexts=frozenset({"lung", "liver"}),
    )
    res = evaluate_species_catalog(
        spec,
        proteins,
        sequences=seqs,
        rare_cell_counts=umi,
        comparison_normals=normals,
        comparison_reproductive=repro,
        comparison_tumors=tumors,
    )
    by_label = res.catalog.set_index("protein_group_label")
    assert by_label.loc["Trap1a", "classification"] == "Unsafe"
    assert by_label.loc["Trap1a", "rare_cell_status"] == "established"
    assert by_label.loc["Magea_clean", "classification"] == "CTA"
    assert by_label.loc["Magea_clean", "cta_tier"] == "extended_only"
    assert (
        by_label.loc["Magea_clean", "reproductive_qualifying_compartments"] == "coagulating_gland"
    )

    with pytest.raises(ValueError, match="taxon_id"):
        SpeciesScreeningSpec(species="dog", reference=mouse_ref)


def test_uniform_species_module_surface_cli_and_explicit_unassessed_nominations(capsys):
    for mod, expected_taxon in ((human, 9606), (canine, 9615), (mouse, 10090)):
        cov = mod.target_selection_coverage()
        assert cov["taxon_id"] == expected_taxon
        assert set(cov["counts"]) == set(target_selection.CATEGORIES)
        assert mod.host_compartments() == tuple(cov["anatomy"]["reproductive_compartments"])
        table = mod.target_selection()
        assert (table.taxon_id == expected_taxon).all()
        assert set(target_selection.SHARED_CHECKS) <= set(table.columns)

    dog_rows = canine.target_selection()
    assert len(dog_rows) == 14
    assert dog_rows.classification.eq("Other:Review").all()
    assert not dog_rows.is_target_candidate.any()

    mouse_rows = mouse.target_selection("ENSMUSG00000051257")
    assert len(mouse_rows) == 1
    assert mouse_rows.iloc[0].classification == "Other:Review"
    assert mouse.target_selection(category="CTA").empty

    assert cli.main(["antigens", "funnel", "--species", "dog", "--format", "json"]) == 0
    funnel_rows = json.loads(capsys.readouterr().out)
    assert funnel_rows[0]["remaining"] == 14
    assert funnel_rows[-1]["remaining"] == 0

    assert cli.main(["antigens", "explain", "P1A", "--species", "mouse", "--format", "json"]) == 0
    exp = json.loads(capsys.readouterr().out)
    assert len(exp) == 1 and exp[0]["classification"] == "Other:Review"

    assert (
        cli.main(["antigens", "list", "--set", "Other", "--species", "dog", "--format", "json"])
        == 0
    )
    listed = json.loads(capsys.readouterr().out)
    assert len(listed) == 14

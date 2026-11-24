"""Cross-species policies, target-selection, sex-restricted subsets & peptide subtraction."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from oncoref import antigen_sets, canine, cli, host_context, human, mouse, target_selection
from oncoref.expression_contract import ExpressionContract, ExpressionMeasurement
from oncoref.species_expression import ExpressionDataset
from oncoref.species_policies import (
    ROLE_IDS,
    all_species_policies,
    role_comparison_table,
    species_policy,
)
from oncoref.species_target_selection import (
    SpeciesScreeningSpec,
    anatomy_for_species,
    assess_normal_rna_dataset_gates,
    assess_pooled_protein_gates,
    classify_sex_restriction,
    evaluate_species_catalog,
    group_proteome,
    group_reference_proteins,
    observations_from_expression_dataset,
)

ROOT = Path(__file__).resolve().parents[1]


def _make_seq(tag: str, length: int = 40) -> str:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    digest = hashlib.sha256(tag.encode("utf-8")).digest()
    return "".join(
        alphabet[(digest[i % len(digest)] + 7 * i) % len(alphabet)] for i in range(length)
    )


def test_declarative_species_policies_cover_all_eight_roles_and_comparative_anatomy():
    policies = all_species_policies()
    assert [p.species.taxon_id for p in policies] == [9606, 9615, 10090]
    for pol in policies:
        table = pol.dataset_table()
        assert set(table.role) == set(ROLE_IDS)
        serialized = pol.as_dict()
        assert serialized["taxon_id"] == pol.species.taxon_id
        assert set(serialized["sex_subset_compartments"]) == {
            "MRA / CTA-M",
            "FRA / CTA-F",
            "ERA / CTA-E",
            "SRA / CTA-MF",
        }
    matrix = role_comparison_table()
    assert len(matrix) == len(ROLE_IDS)
    assert species_policy("human").thresholds.negative_ihc_multipliers == {
        "Approved": 2.0,
        "Supported": 4.0,
        "Enhanced": 8.0,
    }
    assert species_policy("dog").thresholds.negative_ihc_multipliers == {}
    assert species_policy("mouse").thresholds.negative_ihc_multipliers == {}


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


def test_human_sex_restricted_subsets_match_published_csv_and_selection_api(capsys):
    # Frozen 2026-11-10 rebuild CSV matches exact reconstructed 2026-11-20.2 baseline
    frozen_published = pd.read_csv(
        ROOT / "docs/audits/target-selection-20261110/rebuild/summary/sex-restricted-antigens.csv"
    )
    frozen_computed = classify_sex_restriction(
        target_selection._baseline_2026_11_20_frame(),
        species="human",
        host_evidence=target_selection._baseline_2026_11_20_host_exposures(),
    )
    pd.testing.assert_frame_equal(frozen_computed, frozen_published)

    # Live promoted catalog (#775 + #776) matches docs/antigen-whitepaper/sex-restricted-antigens.csv
    published = pd.read_csv(ROOT / "docs/antigen-whitepaper/sex-restricted-antigens.csv")
    computed = human.sex_restricted_antigens()
    pd.testing.assert_frame_equal(computed, published)

    by_lbl = computed.set_index("protein_group_label")
    for rescued_cta_m in (
        "PRM1",
        "SPATA19",
        "CT83",
        "TGM4",
        "SEMG2",
        "ODF1",
        "ADAM2",
        "PLAC1",
        "MAGEB6",
        "FOXR2",
        "SPATA31A3",
    ):
        assert by_lbl.loc[rescued_cta_m, "sex_restriction_code"] == "CTA-M"
    assert by_lbl.loc["DUX4", "sex_restriction_code"] == "CTA-E"
    for shared_cta_mf in ("PRAME", "LALBA", "SCGB1D4", "CTAG2", "DKKL1"):
        assert by_lbl.loc[shared_cta_mf, "sex_restriction_code"] == "CTA-MF"

    cta_m = antigen_sets.proteins("CTA-M")
    cta_e = antigen_sets.proteins("CTA-E")
    cta_f = antigen_sets.proteins("CTA-F")
    cta_mf = antigen_sets.proteins("CTA-MF")
    assert len(cta_m) + len(cta_e) + len(cta_f) + len(cta_mf) == len(antigen_sets.proteins("CTA"))
    assert set(cta_m.entity_id) == set(
        published.loc[published.sex_restriction_code.eq("CTA-M"), "entity_id"]
    )
    assert set(cta_e.entity_id) == set(
        published.loc[published.sex_restriction_code.eq("CTA-E"), "entity_id"]
    )
    assert set(cta_f.entity_id) == set(
        published.loc[published.sex_restriction_code.eq("CTA-F"), "entity_id"]
    )
    assert len(antigen_sets.proteins("MRA")) == int(
        published.sex_restriction_class.eq("Male-restricted").sum()
    )
    assert cli.main(["antigens", "list", "--set", "CTA-E", "--ids", "--format", "json"]) == 0
    cli_ids = {r["gene_id"] for r in json.loads(capsys.readouterr().out)}
    assert cli_ids == antigen_sets.gene_ids("CTA-E")


def test_end_to_end_canine_pipeline_gates_sex_subsets_and_peptide_subtraction():
    dog_pol = canine.screening_policy()
    dog_ref = dog_pol.primary_reference
    unsafe_seq = "M" + "ACDEFGHIKLMNPQRSTVWY" * 2
    contained_seq = unsafe_seq[:39]
    cta_m_seq = _make_seq("DOG_CTA_M", 45)
    cta_e_seq = _make_seq("DOG_CTA_E", 45)
    ca_seq = _make_seq("DOG_CA", 45)
    fra_seq = _make_seq("DOG_FRA", 45)
    prot_unsafe_seq = _make_seq("DOG_PFIZER_HIT", 45)
    bone_seq = _make_seq("DOG_SP7", 45)
    low_seq = _make_seq("DOG_LOW", 45)

    proteins, seqs = group_reference_proteins(
        [
            {
                "gene_id": "400001",
                "symbol": "CTAG_A",
                "transcript_id": "XM_1.1",
                "protein_id": "XP_1.1",
                "sequence": cta_m_seq,
            },
            {
                "gene_id": "400002",
                "symbol": "CTAG_B",
                "transcript_id": "XM_2.1",
                "protein_id": "XP_2.1",
                "sequence": cta_m_seq,
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
                "symbol": "DOG_FRA",
                "transcript_id": "XM_4.1",
                "protein_id": "XP_4.1",
                "sequence": fra_seq,
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
            {
                "gene_id": "400010",
                "symbol": "DOG_CTA_E",
                "transcript_id": "XM_10.1",
                "protein_id": "XP_10.1",
                "sequence": cta_e_seq,
            },
            {
                "gene_id": "400011",
                "symbol": "DOG_PFIZER_HIT",
                "transcript_id": "XM_11.1",
                "protein_id": "XP_11.1",
                "sequence": prot_unsafe_seq,
            },
        ],
        species="dog",
        reference=dog_ref,
    )
    assert len(proteins) == 10
    proteins.loc[proteins.protein_group_label.eq("DOG_BONE"), "home_tissue"] = "bone"

    gene_ids = [f"40000{i}" for i in range(1, 10)] + ["400010", "400011"]
    # Sample all canine male and female adult reproductive compartments plus placenta so
    # opposite-sex clearance is measured rather than guessed.
    repro_tissues = [
        ("L_tes1", "D2", "testis"),
        ("L_epi1", "D2", "epididymis"),
        ("L_pro1", "D2", "prostate"),
        ("L_ova1", "D3", "ovary"),
        ("L_ft1", "D3", "fallopian tube"),
        ("L_ut1", "D3", "uterus"),
        ("L_uth1", "D3", "uterine horn"),
        ("L_end1", "D3", "endometrium"),
        ("L_cer1", "D3", "cervix"),
        ("L_vag1", "D3", "vagina"),
        ("L_lmam1", "D3", "lactating mammary"),
        ("L_rmam1", "D3", "mammary gland"),
        ("L_pla1", "D4", "placenta"),
    ]
    base_samples = [
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
        *(
            {
                "taxon_id": 9615,
                "library_id": lib,
                "run_ids": [f"R_{lib}"],
                "donor_id": donor,
                "donor_scope": "study",
                "pooling": "single",
                "tissue": tissue,
            }
            for lib, donor, tissue in repro_tissues
        ),
    ]
    normal_samples = pd.DataFrame(base_samples).set_index("library_id", drop=False)

    # 11 genes: 400001(CTAG_A), 400002(CTAG_B), 400003(CA), 400004(FRA), 400005(UNSAFE),
    # 400006(CONTAINED), 400007(BONE), 400008(LOW), 400009(UNTRANS), 400010(CTA_E), 400011(PFIZER_HIT)
    cols_data: dict[str, list[float]] = {
        "L_liv1": [0.2, 0.2, 0.5, 0.3, 80.0, 0.1, 0.2, 0.5, 0.1, 0.1, 0.2],
        "L_hrt1": [0.1, 0.1, 0.4, 0.2, 90.0, 0.1, 0.1, 0.4, 0.1, 0.1, 0.2],
        "L_thy1": [25.0, 25.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "L_gal1": [1.0, 0.5, 0.2, 0.2, 20.0, 0.1, 0.2, 0.2, 0.1, 0.2, 0.2],
    }
    for lib, _donor, tissue in repro_tissues:
        if tissue == "testis":
            cols_data[lib] = [25.0, 15.0, 1.0, 0.1, 30.0, 25.0, 1.0, 2.0, 20.0, 0.2, 25.0]
        elif tissue == "uterine horn":
            cols_data[lib] = [0.1, 0.1, 0.5, 35.0, 15.0, 0.0, 0.0, 1.0, 0.0, 0.1, 0.1]
        elif tissue == "placenta":
            cols_data[lib] = [0.0, 0.0, 0.2, 0.2, 5.0, 0.0, 0.0, 0.2, 0.0, 42.0, 0.0]
        else:
            cols_data[lib] = [0.1, 0.1, 0.2, 0.1, 5.0, 0.0, 0.0, 0.2, 0.0, 0.1, 0.1]

    normal_vals = pd.DataFrame(cols_data, index=gene_ids)
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

    # Pooled canine proteomics dataset catches DOG_PFIZER_HIT in heart despite low mRNA!
    prot_samples = pd.DataFrame(
        [
            {
                "taxon_id": 9615,
                "library_id": "P_hrt_pool",
                "run_ids": ["MS1"],
                "donor_id": None,
                "donor_scope": "study",
                "pooling": "pooled",
                "tissue": "heart",
            }
        ]
    ).set_index("library_id", drop=False)
    prot_vals = pd.DataFrame({"P_hrt_pool": [0.0, 5e6]}, index=["PG_CLEAN", "PG_HIT"])
    prot_contract = ExpressionContract(
        9615,
        "synthetic-pfizer-dog-protein",
        "1.0",
        "normal",
        None,
        "unresolved",
        "UniProt/FASTA group",
        ExpressionMeasurement("protein_group", "iBAQ", "linear", False),
        {"matrix": "c" * 64},
        "Synthetic pooled canine TMT-iBAQ matrix.",
    )
    prot_ds = ExpressionDataset(
        prot_vals,
        pd.DataFrame({"symbol": ["PG_CLEAN", "PG_HIT"]}, index=["PG_CLEAN", "PG_HIT"]),
        prot_samples,
        prot_contract,
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
            "T_osa1": [12.0, 8.0, 25.0, 1.0, 50.0, 20.0, 45.0, 3.0, 30.0, 38.0, 30.0],
            "T_osa_pool": [15.0, 10.0, 30.0, 0.5, 40.0, 18.0, 50.0, 2.0, 25.0, 40.0, 28.0],
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

    rna_gates = assess_normal_rna_dataset_gates(normal_ds, proteins)
    prot_gates = assess_pooled_protein_gates(
        prot_ds,
        proteins,
        feature_to_gene={"PG_CLEAN": "400003", "PG_HIT": "400011"},
        min_abundance=1e5,
    )
    all_gates = pd.concat([rna_gates, prot_gates], ignore_index=True)

    norm_obs = observations_from_expression_dataset(
        normal_ds, proteins, comparison_id="canine-tmmcpm-v1"
    )
    tum_obs = observations_from_expression_dataset(
        tumor_ds, proteins, comparison_id="canine-tmmcpm-v1"
    )

    spec = dog_pol.to_screening_spec(
        comparison_id="canine-tmmcpm-v1",
        comparability_reference="synthetic-canine-receipt",
        comparison_unit="TMM-CPM",
        required_normal_contexts=frozenset({"liver", "heart"}),
    )
    res = canine.evaluate_target_catalog(
        spec,
        proteins,
        sequences=seqs,
        normal_gates=all_gates,
        comparison_normals=norm_obs["normals"],
        comparison_reproductive=norm_obs["reproductive"],
        comparison_tumors=tum_obs["tumors"],
    )
    by_label = res.catalog.set_index("protein_group_label")
    assert by_label.loc["CTAG_A/CTAG_B", "classification"] == "CTA"
    assert by_label.loc["CTAG_A/CTAG_B", "sex_restriction_code"] == "CTA-M"
    assert by_label.loc["DOG_CTA_E", "classification"] == "CTA"
    assert by_label.loc["DOG_CTA_E", "sex_restriction_code"] == "CTA-E"
    assert by_label.loc["DOG_FRA", "classification"] == "Reproductive Antigen only"
    assert by_label.loc["DOG_FRA", "sex_restriction_code"] == "FRA-only"
    assert by_label.loc["DOG_CA", "classification"] == "Cancer Antigen only"
    assert by_label.loc["DOG_UNSAFE", "classification"] == "Unsafe"
    assert by_label.loc["DOG_PFIZER_HIT", "classification"] == "Unsafe"
    assert "pooled iBAQ" in by_label.loc["DOG_PFIZER_HIT", "Unsafe_reasons"]
    assert by_label.loc["DOG_CONTAINED", "classification"] == "Other:Sequence"
    assert by_label.loc["DOG_BONE", "classification"] == "Other:Unscreened tissue"
    assert by_label.loc["DOG_LOW", "classification"] == "Other:Expression"
    assert by_label.loc["DOG_UNTRANS", "classification"] == "Other:Review"

    # Sex-specific selections work identically on the canine result!
    assert res.select("CTA-M").proteins().protein_group_label.tolist() == ["CTAG_A/CTAG_B"]
    assert res.select("CTA-E").proteins().protein_group_label.tolist() == ["DOG_CTA_E"]
    assert res.select("FRA").proteins().protein_group_label.tolist() == ["DOG_FRA"]

    # Protecting canine testis excludes CTAG_A/CTAG_B while retaining adult-absent DOG_CTA_E,
    # and moves CTAG_A/CTAG_B into the canine peptide-subtraction background!
    sel_protect_testis = res.select(
        "CTA", host=host_context.HostContext(protect=("testis",), species="dog")
    )
    assert sel_protect_testis.proteins().protein_group_label.tolist() == ["DOG_CTA_E"]
    bg_index = sel_protect_testis.background_index("Unsafe")
    assert bg_index.subtract(cta_m_seq, lengths=[9]).empty
    assert not bg_index.subtract(cta_e_seq, lengths=[9]).empty


def test_mouse_rare_cell_umi_fail_closed_sex_subsets_and_species_boundaries():
    mouse_pol = mouse.screening_policy()
    mouse_ref = mouse_pol.primary_reference
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
    # Only coagulating_gland is measured here (female organs unmeasured) -> must fail closed to CTA-MF!
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
    spec = mouse_pol.to_screening_spec(
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
    # Unmeasured female mouse organs fail closed (never claimed female-safe CTA-M without measurements):
    assert by_label.loc["Magea_clean", "sex_restriction_code"] == "CTA-MF"
    assert res.select("CTA-M").proteins().empty

    with pytest.raises(ValueError, match="taxon_id"):
        SpeciesScreeningSpec(species="dog", reference=mouse_ref)


def test_uniform_species_module_surface_cli_and_explicit_unassessed_nominations(capsys):
    for mod, expected_taxon in ((human, 9606), (canine, 9615), (mouse, 10090)):
        cov = mod.target_selection_coverage()
        assert cov["taxon_id"] == expected_taxon
        assert set(cov["counts"]) == set(target_selection.CATEGORIES)
        assert mod.host_compartments() == tuple(cov["anatomy"]["reproductive_compartments"])
        assert set(mod.screening_datasets().role) == set(ROLE_IDS)
        table = mod.target_selection()
        assert (table.taxon_id == expected_taxon).all()
        assert {"policy", "target_policy_id", *target_selection.SHARED_CHECKS} <= set(table.columns)

    # Bidirectional HostContext species guard rejects cross-species host contexts in both directions:
    with pytest.raises(
        ValueError, match="HostContext species 'dog' disagrees with catalog species 'human'"
    ):
        antigen_sets.select(
            "CTA", host=host_context.HostContext(protect=("testis",), species="dog")
        )
    with pytest.raises(
        ValueError, match="HostContext species 'human' disagrees with catalog species 'dog'"
    ):
        antigen_sets.select(
            "CTA",
            species="dog",
            host=host_context.HostContext(protect=("testis",), species="human"),
        )

    dog_rows = canine.target_selection()
    assert len(dog_rows) == 14
    assert dog_rows.classification.eq("Other:Review").all()
    assert not dog_rows.is_target_candidate.any()
    assert canine.sex_restricted_antigens().empty
    assert not classify_sex_restriction(species="mouse").empty

    mouse_rows = mouse.target_selection("ENSMUSG00000051257")
    assert len(mouse_rows) == 1
    assert mouse_rows.iloc[0].classification in {"Other:Review", "Unsafe"}
    assert not mouse_rows.iloc[0].is_cta
    assert not mouse.target_selection(category="CTA").empty

    assert cli.main(["antigens", "funnel", "--species", "dog", "--format", "json"]) == 0
    funnel_rows = json.loads(capsys.readouterr().out)
    assert funnel_rows[0]["remaining"] == 14
    assert funnel_rows[-1]["remaining"] == 0

    assert cli.main(["antigens", "explain", "P1A", "--species", "mouse", "--format", "json"]) == 0
    exp = json.loads(capsys.readouterr().out)
    assert len(exp) == 1 and exp[0]["classification"] in {"Other:Review", "Unsafe"}

    assert (
        cli.main(["antigens", "list", "--set", "Other", "--species", "dog", "--format", "json"])
        == 0
    )
    listed = json.loads(capsys.readouterr().out)
    assert len(listed) == 14

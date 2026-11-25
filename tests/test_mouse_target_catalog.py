"""End-to-end tests for the quantitative Ensembl 115 GRCm39 murine antigen catalog (#734)."""

from __future__ import annotations

import json

import pandas as pd

from oncoref import (
    antigen_sets,
    cli,
    data_manifest,
    host_context,
    mouse,
    reference_data,
    target_selection,
)
from oncoref.species_policies import ROLE_IDS


def test_mouse_catalog_genome_wide_partition_and_eight_active_pillars():
    cov = mouse.target_selection_coverage()
    assert cov["taxon_id"] == 10090
    assert cov["species"] == "mouse"
    assert cov["selection_status"] == "evaluated"
    assert cov["comparison_unit"] == "TPM"
    assert cov["thresholds"]["negative_ihc_multipliers"] == {}
    assert set(cov["populated_inputs"]) == set(ROLE_IDS)
    assert cov["missing_inputs"] == []

    ds = mouse.screening_datasets()
    assert set(ds.role) == set(ROLE_IDS)
    assert ds.status.eq("active_in_catalog").all()

    cat = mouse.target_selection()
    assert len(cat) > 21000
    assert cat.entity_id.is_unique
    loci = cat.member_gene_ids.map(json.loads).explode()
    assert loci.is_unique
    assert len(loci) > 21700
    assert all(str(g).startswith("ENSMUSG") for g in loci)

    ca = cat.is_cancer_antigen.astype(bool)
    ra = cat.is_reproductive_antigen.astype(bool)
    cta = cat.is_cta.astype(bool)
    cand = cat.is_target_candidate.astype(bool)
    assert (ca & ra).equals(cta)
    assert (ca | ra).equals(cand)
    assert int(cta.sum()) > 0
    assert int((ca & ~ra).sum()) > 0
    assert int((ra & ~ca).sum()) > 0
    assert int(cat.category.eq("Unsafe").sum()) > 10000
    assert {c: int(cat.category.eq(c).sum()) for c in target_selection.CATEGORIES} == cov["counts"]


def test_mouse_sources_and_strain_explicit_tumor_cohorts():
    sources = mouse.screening_sources()
    assert len(sources) >= 24
    assert sources.source_id.is_unique
    assert sources.sha256.str.fullmatch(r"[0-9a-f]{64}").all()
    assert (sources["bytes"] > 0).all()

    cohorts = mouse.tumor_cohorts()
    assert set(cohorts.cancer_code) == {
        "MC38",
        "B16F10",
        "CT26",
        "LLC",
        "4T1",
        "KP_LUAD",
        "BMPC_PRAD",
        "P815",
    }
    assert set(cohorts.strain) == {"C57BL/6", "BALB/c", "FVB", "DBA/2"}
    assert int(cohorts.n_samples.sum()) == 36


def test_p1a_trap1a_alias_resolution_and_legacy_cta_snapshot_continuity():
    # Frozen historical cta_selection() / cta_snapshot() surface remains unchanged
    legacy_sel = mouse.cta_selection()
    assert len(legacy_sel) == 1
    assert legacy_sel.assessment_status.eq("unassessed").all()
    assert mouse.cta_snapshot().coverage["status"] == "partial"

    # Quantitative target_selection resolves both P1A alias and ENSMUSG00000051257
    by_alias = mouse.target_selection("P1A")
    by_ens = mouse.target_selection("ENSMUSG00000051257")
    by_sym = mouse.target_selection("Trap1a")
    assert len(by_alias) == len(by_ens) == len(by_sym) == 1
    pd.testing.assert_frame_equal(by_alias, by_ens)
    pd.testing.assert_frame_equal(by_alias, by_sym)

    row = by_alias.iloc[0]
    assert row.published_nomination
    assert not row.is_cta
    assert row.classification in {"Unsafe", "Other:Review"}


def test_mouse_sixteen_compartment_host_context_and_sex_restricted_subsets():
    anat = mouse.screening_policy().anatomy
    assert len(anat.reproductive_compartments) == 16
    assert "coagulating_gland" in anat.male_accessory_compartments
    assert "seminal vesicle" in anat.male_accessory_compartments
    assert "uterine_horn" in anat.female_tract_compartments
    assert "lactating_mammary" in anat.mammary_compartments

    exp = mouse.host_exposure()
    measured_union = {
        comp
        for comps in exp.values()
        for comp, rec in comps.items()
        if rec.get("coverage") == "measured"
    }
    assert measured_union == set(anat.reproductive_compartments)

    sex_df = mouse.sex_restricted_antigens()
    assert len(sex_df) == int(mouse.target_selection().is_reproductive_antigen.sum())
    cta_m = antigen_sets.proteins("CTA-M", species="mouse")
    cta_e = antigen_sets.proteins("CTA-E", species="mouse")
    cta_f = antigen_sets.proteins("CTA-F", species="mouse")
    cta_mf = antigen_sets.proteins("CTA-MF", species="mouse")
    total_cta = antigen_sets.proteins("CTA", species="mouse")
    assert len(cta_m) + len(cta_e) + len(cta_f) + len(cta_mf) == len(total_cta)
    assert {"Pbsn", "Sbp", "Tex19.1", "Plac1"} <= set(cta_m.protein_group_label)
    assert {"Serpinb9e", "Pramel13"} <= set(cta_e.protein_group_label)
    assert {"Gzmf", "Gzmg", "Gzmd"} <= set(cta_f.protein_group_label)

    # Sorted & activated immune gates (ImmGen ULI/MNP + LPS BMDM) exclude cytokines & myeloid markers
    by_lbl = mouse.target_selection().set_index("protein_group_label")
    for imm_unsafe in (
        "Csf3",
        "Cxcl2",
        "Il23a",
        "Acod1",
        "Il10",
        "Tnfsf4",
        "Tnfsf8",
        "Zfp658",
        "Glrp1",
        "Utf1",
        "Peg12",
        "Cpvl",
        "Dppa2",
        "Qrfp",
        "Otos",
    ):
        assert by_lbl.loc[imm_unsafe, "classification"] == "Unsafe"

    # Pure HostContext protection of non-pregnant, non-lactating female adult organs keeps all CTA-M eligible,
    # while TargetFilter additionally requires 10x tumor fold over max(H0, protected female TPM)
    f_host = host_context.from_legacy(
        "female-patient",
        protect_tissues=anat.female_adult_protected_compartments,
        lactating=False,
        species="mouse",
    )
    assert host_context.assess(cta_m, f_host, evidence=exp)["pass_host"].all()
    female_sel = antigen_sets.select("CTA", species="mouse", host=f_host).proteins()
    assert not female_sel.empty
    assert set(female_sel.entity_id) <= (set(cta_m.entity_id) | set(cta_e.entity_id))

    cta_p = antigen_sets.proteins("CTA-P", species="mouse")
    ra_p = antigen_sets.proteins("RA-P", species="mouse")
    assert len(cta_p) > 0
    assert set(cta_p.entity_id) <= set(total_cta.entity_id)
    assert set(cta_p.entity_id) <= set(ra_p.entity_id)
    assert set(ra_p.entity_id) == set(antigen_sets.proteins("PRA", species="mouse").entity_id)

    protect_testis = antigen_sets.select(
        "CTA",
        species="mouse",
        host=host_context.HostContext(protect=("testis",), species="mouse"),
    ).proteins()
    assert len(protect_testis) < len(total_cta)


def test_mouse_packaged_proteome_sequences_and_peptide_subtraction(capsys):
    cta_seqs = antigen_sets.sequences("CTA", species="mouse")
    assert not cta_seqs.empty
    assert cta_seqs.attrs["missing_gene_ids"] == []
    assert cta_seqs.attrs["sequence_source"] == "packaged_proteome"
    assert cta_seqs["matches_representative_sequence"].any()

    # Build a PeptideIndex over a compact slice of mouse Unsafe proteins and verify subtraction
    cat = mouse.target_selection()
    unsafe_slice = cat[cat.category.eq("Unsafe")].head(25)
    seq_table = antigen_sets._sequences(unsafe_slice, "Unsafe", species="mouse")
    from oncoref.peptide_subtraction import PeptideIndex

    idx = PeptideIndex(seq_table, lengths=[9])
    sample_seq = str(seq_table.iloc[0].sequence)
    assert idx.subtract(sample_seq, lengths=[9]).empty

    assert (
        cli.main(
            [
                "antigens",
                "list",
                "--set",
                "CTA-M",
                "--species",
                "mouse",
                "--ids",
                "--format",
                "json",
            ]
        )
        == 0
    )
    cli_ids = {r["gene_id"] for r in json.loads(capsys.readouterr().out)}
    assert cli_ids == antigen_sets.gene_ids("CTA-M", species="mouse")


def test_mouse_heavy_tables_are_externalized_and_pinned_in_reference_releases():
    import hashlib
    from pathlib import Path

    import oncoref
    from oncoref import species_target_selection

    audit = data_manifest.package_data_size_audit().set_index("dataset_stem")
    meta = species_target_selection.packaged_species_metadata("mouse")
    assert meta is not None
    expected_shas = {
        "mouse-target-selection": meta["data_sha256"],
        "mouse-host-exposure": meta["host_exposure_sha256"],
        "mouse-proteome-sequences": meta["sequences_sha256"],
    }
    data_dir = Path(oncoref.__file__).resolve().parent / "data"
    for stem, expected_sha in expected_shas.items():
        assert stem in data_manifest.EXTERNALIZED_ORIGINATED
        assert stem in reference_data.REFERENCE_SOURCES
        assert not bool(audit.loc[stem, "shipped_in_wheel"])
        prov = reference_data.provenance(stem)
        assert prov["sha256"] == expected_sha
        assert "reference-mouse-20261125" in prov["url"]
        p = data_dir / f"{stem}.csv.gz"
        if p.exists():
            assert hashlib.sha256(p.read_bytes()).hexdigest() == expected_sha

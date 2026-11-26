"""End-to-end tests for the quantitative Ensembl 98 CanFam3.1 canine antigen catalog."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import oncoref
from oncoref import (
    antigen_sets,
    canine,
    cli,
    data_manifest,
    host_context,
    reference_data,
    species_target_selection,
    target_selection,
)
from oncoref.species_policies import ROLE_IDS


def test_dog_catalog_genome_wide_partition_and_seven_active_pillars():
    cov = canine.target_selection_coverage()
    assert cov["taxon_id"] == 9615
    assert cov["species"] == "dog"
    assert cov["selection_status"] == "evaluated"
    assert cov["comparison_unit"] == "TPM"
    assert cov["thresholds"]["negative_ihc_multipliers"] == {}
    assert len(set(cov["populated_inputs"])) == 7
    assert any("canine_normal_scrna_gap" in m for m in cov["missing_inputs"])

    ds = canine.screening_datasets()
    assert set(ds.role) == set(ROLE_IDS)
    assert int(ds.status.eq("active_in_catalog").sum()) >= 7

    cat = canine.target_selection()
    assert len(cat) == 20082
    assert cat.entity_id.is_unique
    loci = cat.member_gene_ids.map(json.loads).explode()
    assert loci.is_unique
    assert len(loci) == 20199
    assert all(str(g).startswith("ENSCAFG") for g in loci)

    ca = cat.is_cancer_antigen.astype(bool)
    ra = cat.is_reproductive_antigen.astype(bool)
    cta = cat.is_cta.astype(bool)
    cand = cat.is_target_candidate.astype(bool)
    assert (ca & ra).equals(cta)
    assert (ca | ra).equals(cand)
    assert int(cta.sum()) >= 30
    assert int((ca & ~ra).sum()) >= 50
    assert int((ra & ~ca).sum()) >= 200
    assert int(cat.category.eq("Unsafe").sum()) > 15000
    assert {c: int(cat.category.eq(c).sum()) for c in target_selection.CATEGORIES} == cov["counts"]


def test_dog_sources_cohorts_and_all_fourteen_curated_nominations_resolve():
    sources = canine.screening_sources()
    assert len(sources) >= 30
    assert sources.source_id.is_unique
    assert sources.sha256.str.fullmatch(r"[0-9a-f]{64}").all()
    assert (sources["bytes"] > 0).all()

    cohorts = canine.tumor_cohorts()
    assert set(cohorts.cancer_code) == {
        "OSA_NCI",
        "OSA_COTC030",
        "OSA_CCOGC",
        "CMT",
        "HSA",
        "LYM",
        "DLBCL",
        "GLIOMA",
    }
    assert int(cohorts.n_samples.sum()) > 500

    # Frozen legacy cta_selection() surface remains unchanged while target_selection(cand) resolves all 14
    legacy = canine.cta_selection()
    assert len(legacy) == 14
    assert legacy.assessment_status.eq("unassessed").all()

    for cand in canine.cta_candidates():
        resolved = canine.target_selection(cand["id"])
        assert not resolved.empty, f"Failed to resolve canine nomination {cand['id']!r}"
        assert resolved["published_nomination"].any()

    assert canine.target_selection("PRAME").iloc[0].classification == "CTA"
    assert canine.target_selection("MAGEB10").iloc[0].classification == "Reproductive Antigen only"
    assert canine.target_selection("SSX6").iloc[0].classification == "Reproductive Antigen only"
    for unsafe_nom in ("AKR1E2", "IGF2BP3", "NT5C1B", "PAGE3-like", "SPECC1", "TPX2"):
        assert canine.target_selection(unsafe_nom).iloc[0].classification == "Unsafe"
    for retro_or_unsafe in (
        "ROMO1",
        "RPL10L",
        "ENSCAFG00000017286",
        "ENSCAFG00000007385",
        "ENSCAFG00000014111",
    ):
        assert canine.target_selection(retro_or_unsafe).iloc[0].classification != "CTA"
    for review_nom in ("CAGE1", "CTCFL", "DDX53"):
        assert canine.target_selection(review_nom).iloc[0].classification == "Other:Review"


def test_dog_comparative_anatomy_host_context_and_four_way_sex_subsets(capsys):
    anat = canine.screening_policy().anatomy
    assert "seminal vesicle" not in anat.male_accessory_compartments
    assert "uterine_horn" in anat.female_tract_compartments

    exp = canine.host_exposure()
    sample_eid = canine.target_selection("PRAME").iloc[0].entity_id
    assert exp[sample_eid]["fallopian tube"]["coverage"] == "proteome_only"
    assert exp[sample_eid]["uterine_horn"]["coverage"] == "measured"

    # Direct HostContext protecting fallopian tube fails closed ("unresolved") on proteome_only
    ft_res = host_context.assess(
        canine.target_selection("PRAME"),
        host_context.HostContext(protect=("fallopian tube",), species="dog"),
        evidence=exp,
    )
    assert ft_res.iloc[0].host_status == "unresolved"

    # Four-way sex-restriction partition sums identically to total canine CTAs;
    # unmeasured male canine resting_mammary fails closed (CSN2/CSN3 -> CTA-MF; ovary/tract -> RA-F)
    cta_m = antigen_sets.proteins("CTA-M", species="dog")
    cta_p = antigen_sets.proteins("CTA-P", species="dog")
    cta_f = antigen_sets.proteins("CTA-F", species="dog")
    cta_mf = antigen_sets.proteins("CTA-MF", species="dog")
    ra_f = antigen_sets.proteins("RA-F", species="dog")
    total_cta = antigen_sets.proteins("CTA", species="dog")
    assert len(cta_m) + len(cta_p) + len(cta_f) + len(cta_mf) == len(total_cta) == 54
    assert {"PRAME", "NAA11", "PLD6", "PSMA8"} <= set(cta_m.protein_group_label)
    assert {"CSN2", "CSN3"} <= set(cta_mf.protein_group_label)
    assert {"BMP15", "ZP3", "KPNA7"} <= set(ra_f.protein_group_label)
    assert {"IGF2BP1", "GPR50"} <= set(cta_p.protein_group_label)
    assert set(cta_p.entity_id) == set(antigen_sets.proteins("CTA-E", species="dog").entity_id)

    # Spayed female dog HostContext retains male-restricted and pregnancy-only CTAs
    spayed_bitch = host_context.HostContext(
        absent=(
            "testis",
            "epididymis",
            "prostate",
            "ovary",
            "uterus",
            "uterine_horn",
            "endometrium",
            "lactating_mammary",
        ),
        protect=("cervix", "vagina", "resting_mammary"),
        species="dog",
    )
    assert host_context.assess(cta_m, spayed_bitch, evidence=exp)["pass_host"].all()

    # Packaged canine sequences and PeptideIndex subtraction
    cta_seqs = antigen_sets.sequences("CTA", species="dog")
    assert not cta_seqs.empty
    assert cta_seqs.attrs["missing_gene_ids"] == []
    assert cta_seqs.attrs["sequence_source"] == "packaged_proteome"

    from oncoref.peptide_subtraction import PeptideIndex

    unsafe_slice = canine.target_selection(category="Unsafe").head(25)
    seq_table = antigen_sets._sequences(unsafe_slice, "Unsafe", species="dog")
    idx = PeptideIndex(seq_table, lengths=[9])
    assert idx.subtract(str(seq_table.iloc[0].sequence), lengths=[9]).empty

    assert (
        cli.main(
            ["antigens", "list", "--set", "CTA-M", "--species", "dog", "--ids", "--format", "json"]
        )
        == 0
    )
    cli_ids = {r["gene_id"] for r in json.loads(capsys.readouterr().out)}
    assert cli_ids == antigen_sets.gene_ids("CTA-M", species="dog")


def test_dog_heavy_tables_are_externalized_and_pinned_in_reference_releases():
    audit = data_manifest.package_data_size_audit().set_index("dataset_stem")
    meta = species_target_selection.packaged_species_metadata("dog")
    assert meta is not None
    expected_shas = {
        "dog-target-selection": meta["data_sha256"],
        "dog-host-exposure": meta["host_exposure_sha256"],
        "dog-proteome-sequences": meta["sequences_sha256"],
    }
    data_dir = Path(oncoref.__file__).resolve().parent / "data"
    for stem, expected_sha in expected_shas.items():
        assert stem in data_manifest.EXTERNALIZED_ORIGINATED
        assert stem in reference_data.REFERENCE_SOURCES
        assert not bool(audit.loc[stem, "shipped_in_wheel"])
        prov = reference_data.provenance(stem)
        assert prov["sha256"] == expected_sha
        assert "reference-dog-20261126" in prov["url"]
        p = data_dir / f"{stem}.csv.gz"
        if p.exists():
            assert hashlib.sha256(p.read_bytes()).hexdigest() == expected_sha

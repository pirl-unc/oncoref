import hashlib
import json
from types import SimpleNamespace

import pandas as pd
import pytest

from oncoref import antigen_sets


def test_current_memberships_and_forbidden_background():
    sets = {g: set(antigen_sets.proteins(g).entity_id) for g in antigen_sets.SETS}
    assert sets["CTA"] == sets["CA"] & sets["RA"]
    assert sets["CA|RA"] == sets["CA"] | sets["RA"]
    assert not sets["CA|RA"] & sets["Unsafe"]
    assert sets["Unsafe|Other"] == sets["Unsafe"] | sets["Other"]
    assert not sets["Unsafe"] & sets["Other"]
    assert not sets["CA|RA"] & sets["Unsafe|Other"]
    assert sets["CA|RA"] | sets["Unsafe|Other"] == set(
        antigen_sets.target_selection.results().entity_id
    )
    assert antigen_sets.gene_ids("CTA") == {
        g for ids in antigen_sets.proteins("CTA").member_gene_ids for g in json.loads(ids)
    }
    with pytest.raises(ValueError):
        antigen_sets.proteins("Neither")


def test_sequences_include_every_annotated_translation_and_report_gaps(monkeypatch):
    monkeypatch.setattr(
        antigen_sets,
        "proteins",
        lambda group: pd.DataFrame(
            {
                "entity_id": ["p1", "p2"],
                "representative_sequence_sha256s": [
                    json.dumps([hashlib.sha256(b"ACDEFGHIK").hexdigest()]),
                    "[]",
                ],
                "member_gene_ids": ['["g1"]', '["g2"]'],
            }
        ),
    )
    transcripts = [
        SimpleNamespace(id="t1", biotype="protein_coding", protein_sequence="ACDEFGHIK"),
        SimpleNamespace(id="t2", biotype="protein_coding_LoF", protein_sequence="ACDEFGHIK"),
        SimpleNamespace(id="t3", biotype="nonsense_mediated_decay", protein_sequence="LMNPQRSTV"),
        SimpleNamespace(id="t4", biotype="retained_intron", protein_sequence=None),
    ]
    genome = SimpleNamespace(
        release=115,
        species=SimpleNamespace(latin_name="homo_sapiens"),
        gene_by_id=lambda gene: SimpleNamespace(transcripts=transcripts if gene == "g1" else []),
    )
    result = antigen_sets.sequences("Unsafe", genome=genome)
    assert result.matches_representative_sequence.tolist() == [True, True, False]
    assert list(result.transcript_id) == ["t1", "t2", "t3"]
    assert set(result.transcript_biotype) == {
        "protein_coding",
        "protein_coding_LoF",
        "nonsense_mediated_decay",
    }
    assert result.attrs["missing_gene_ids"] == ["g2"]
    assert result.attrs["sequence_scope"] == "all annotated protein translations"
    assert result.sequence_sha256.iloc[0] == result.sequence_sha256.iloc[1]


def test_loss_context_is_independent_of_membership_and_retains_organ_scope():
    # A Review hold still has measured normal exposure; bulk testis does not
    # establish selective germ-cell injury or waive endocrine consequences.
    assert antigen_sets.target_selection.results("CAGE1").category.iloc[0] == "Other"
    records = antigen_sets.loss_consequences("CAGE1")
    testis = next(r for r in records if r["expression_context"]["tissue"] == "testis")
    assert testis["matched_level"] == "bulk" and testis["tier"] == "protected"
    assert testis["evidence"] and "androgen" in testis["persistent_loss"]
    # RPS4Y2 (Cancer Antigen only) fails RA's 10x fold requirement, while CTAG2 clears to CTA;
    # both retain organ-level reproductive loss consequences independently of category.
    assert not antigen_sets.target_selection.results("RPS4Y2").is_reproductive_antigen.iloc[0]
    assert antigen_sets.loss_consequences("RPS4Y2")
    assert antigen_sets.target_selection.results("CTAG2").is_cta.iloc[0]
    assert antigen_sets.loss_consequences("CTAG2")
    lalba = antigen_sets.loss_consequences("LALBA")
    assert any(r["expression_context"].get("compartment") == "lactating_breast" for r in lalba)
    with pytest.raises(ValueError, match="Unknown"):
        antigen_sets.loss_consequences("not_a_gene")


def test_proteins_and_gene_ids_accept_host_profiles_and_reject_backgrounds():
    selected = antigen_sets.select("CTA", host_profile="female-patient")
    rows = antigen_sets.proteins("CTA", profile="female-patient")
    assert set(rows.entity_id) == set(selected.proteins().entity_id)
    assert antigen_sets.gene_ids("CTA", profile="female-patient") == selected.gene_ids()
    assert rows.is_cta.all() and rows.category.eq("CTA").all()
    with pytest.raises(ValueError, match="Host profiles"):
        antigen_sets.proteins("Unsafe", profile="female-patient")


def test_proteins_uniform_schema_and_cross_species_wiring():
    import oncoref

    cta = antigen_sets.proteins("CTA")
    cta_male = antigen_sets.proteins("CTA-M")
    cta_peri = antigen_sets.proteins("CTA-P")
    ra_peri = antigen_sets.proteins("RA-P")
    assert list(cta.columns) == list(cta_male.columns) == list(cta_peri.columns)
    assert set(cta_peri.entity_id) <= set(cta.entity_id)
    assert set(cta_peri.entity_id) <= set(ra_peri.entity_id)
    assert set(ra_peri.entity_id) == set(antigen_sets.proteins("PRA").entity_id)
    assert "catalog_category" in cta.columns
    assert cta.attrs["set"] == "CTA"

    dog_ras = antigen_sets.proteins("RA", species="canine")
    assert len(dog_ras) == 431 and "catalog_category" in dog_ras.columns
    assert len(antigen_sets.proteins("CTA", species="canine")) == 118
    assert antigen_sets.proteins("Unsafe", species="canine").attrs["species"] == "dog"
    assert not antigen_sets.proteins("Other", species="canine").empty
    assert antigen_sets.gene_ids("Other", species="canine")
    assert not antigen_sets.proteins("Other", species="mouse").empty

    for mod_name in (
        "antigen_sets",
        "curation",
        "figure_style",
        "host_context",
        "peptide_subtraction",
        "sequence_bundle",
        "species_policies",
        "species_target_selection",
        "target_filters",
        "tcr_trials",
    ):
        assert mod_name in oncoref.__all__
        assert hasattr(oncoref, mod_name)

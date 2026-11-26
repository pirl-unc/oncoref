"""Tests for uniform species-parameterized runtime across human, dog, and mouse (#754)."""

from __future__ import annotations

import json

import pandas as pd

from oncoref import (
    antigen_sets,
    canine,
    cli,
    host_context,
    human,
    mouse,
    species_target_selection,
)


def test_packaged_species_metadata_and_coverage_are_uniform_across_evaluated_species():
    for sp_name, expected_taxon, expected_n in (
        ("human", 9606, 19893),
        ("dog", 9615, 20082),
        ("mouse", 10090, 21433),
    ):
        meta = species_target_selection.packaged_species_metadata(sp_name)
        assert meta is not None
        assert meta["schema"] == "oncoref.species-target-selection-metadata.v1"
        assert meta["taxon_id"] == expected_taxon
        assert meta["species"] == sp_name
        assert meta["selection_status"] == "evaluated"
        assert meta["summary"]["proteins"] == expected_n

        cov = species_target_selection.packaged_species_coverage(sp_name)
        assert cov["taxon_id"] == expected_taxon
        assert cov["selection_status"] == "evaluated"
        assert cov["counts"] == meta["summary"]["categories"]
        assert cov["sex_restricted_counts"] == meta["summary"]["sex_restriction_codes"]


def test_packaged_species_results_and_host_exposure_share_identical_contract():
    h_cat = human.target_selection()
    d_cat = canine.target_selection()
    m_cat = mouse.target_selection()
    for df, sp, tid in ((h_cat, "human", 9606), (d_cat, "dog", 9615), (m_cat, "mouse", 10090)):
        assert {"taxon_id", "species", "entity_id", "protein_group_label", "category"} <= set(
            df.columns
        )
        assert (df["taxon_id"] == tid).all()
        assert (df["species"] == sp).all()
        assert df.attrs["taxon_id"] == tid
        assert df.attrs["species"] == sp

    # Canonical member symbol resolution works through packaged_species_results on human
    ssx2 = human.target_selection("SSX2")
    assert not ssx2.empty
    assert "SSX2/B" in set(ssx2["protein_group_label"])

    # host_context.assess automatically resolves species host exposure when evidence is omitted
    m_cta = mouse.target_selection(category="CTA")
    assessed = host_context.assess(
        m_cta, host_context.HostContext(protect=("testis",), species="mouse")
    )
    assert (assessed["host_status"] != "not_applied").all()
    assert assessed["host_status"].isin(["eligible", "excluded", "unresolved"]).all()


def test_antigen_sets_and_cli_have_uniform_columns_and_behavior(capsys):
    for sp in ("human", "dog", "mouse"):
        cta = antigen_sets.proteins("CTA", species=sp)
        unsafe = antigen_sets.proteins("Unsafe", species=sp)
        other = antigen_sets.proteins("Other", species=sp)
        assert cta.attrs["species"] == unsafe.attrs["species"] == other.attrs["species"] == sp
        assert "taxon_id" in unsafe.columns and "species" in unsafe.columns

    # CLI explain and funnel execute through the identical species-parameterized path
    for sp, gene in (("human", "MAGEA4"), ("mouse", "Trap1a"), ("dog", "AKR1E2")):
        assert cli.main(["antigens", "explain", gene, "--species", sp, "--format", "json"]) == 0
        explained = json.loads(capsys.readouterr().out)
        assert len(explained) == 1
        assert cli.main(["antigens", "funnel", "--species", sp, "--format", "json"]) == 0
        funnel_rows = pd.DataFrame(json.loads(capsys.readouterr().out))
        assert "All assessed proteins" in set(funnel_rows["stage"])

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Driver identities remain source-anchored without optional genome annotations."""

import warnings

import pandas as pd
import pytest

from oncoref import drivers, genome
from oncoref.load_dataset import get_data

TP53 = "ENSG00000141510"
KRAS = "ENSG00000133703"
ZNG1C = "ENSG00000196873"


@pytest.fixture(autouse=True)
def clear_driver_caches():
    for cache in (
        drivers._driver_symbol_index,
        drivers._driver_gene_evidence_frame,
        drivers._driver_variant_evidence_frame,
    ):
        cache.cache_clear()
    yield
    for cache in (
        drivers._driver_symbol_index,
        drivers._driver_gene_evidence_frame,
        drivers._driver_variant_evidence_frame,
    ):
        cache.cache_clear()


def test_complete_driver_tables_and_retired_symbol_queries_without_genome(monkeypatch):
    def unexpected_genome_lookup(*args, **kwargs):
        raise AssertionError("driver identity validation must use the bundled reference")

    monkeypatch.setattr(genome, "canonical_gene_id_and_name", unexpected_genome_lookup)
    with pytest.warns(drivers.DriverGeneIdentityWarning, match="CBWD3.*" + ZNG1C) as caught:
        genes = drivers.driver_gene_evidence_df()
    assert len(caught) == 1
    assert "ZNG1C" in str(caught[0].message)
    pd.testing.assert_frame_equal(genes, get_data("driver-gene-evidence"))
    pd.testing.assert_frame_equal(
        drivers.driver_variant_evidence_df(), get_data("driver-variant-evidence")
    )
    assert (len(genes), len(drivers.driver_variant_evidence_df())) == (739, 579)
    assert len(drivers.driver_legacy_migration_audit_df()) == 1318

    expected = genes.loc[genes.ensembl_gene_id.eq(ZNG1C)].reset_index(drop=True)
    with warnings.catch_warnings(record=True) as repeated:
        warnings.simplefilter("always")
        for query in ("CBWD3", "cbwd3", " ZNG1C ", ZNG1C, ZNG1C + ".1"):
            pd.testing.assert_frame_equal(drivers.driver_gene_evidence(query), expected)
    assert not repeated
    changed = drivers.driver_gene_evidence("CBWD3")
    changed.loc[0, "gene_symbol"] = "changed"
    assert drivers.driver_gene_evidence("CBWD3").iloc[0].gene_symbol == "CBWD3"


@pytest.fixture(params=["gene", "variant"])
def source_table(request, monkeypatch):
    dataset = f"driver-{request.param}-evidence"
    frame = get_data(dataset).iloc[[0]].copy()
    frame["gene_symbol"] = "RETIRED_TEST_SYMBOL"
    frame["ensembl_gene_id"] = TP53
    frame["aliases"] = "TP53"
    monkeypatch.setattr(drivers, "get_data", lambda name, **kwargs: frame)
    return frame, getattr(drivers, f"driver_{request.param}_evidence_df")


@pytest.mark.parametrize("aliases", ["KRAS; TP53", "TP53;KRAS", "UNKNOWN;;KRAS; TP53; TP53"])
def test_all_aliases_are_checked_against_the_recorded_id(source_table, aliases):
    frame, accessor = source_table
    frame["aliases"] = aliases
    before = frame.copy(deep=True)
    with pytest.warns(drivers.DriverGeneIdentityWarning) as caught:
        result = accessor()
    message = str(caught[0].message)
    assert "TP53 -> " + TP53 in message
    assert "KRAS -> " + KRAS in message
    pd.testing.assert_frame_equal(result, before)
    pd.testing.assert_frame_equal(frame, before)


@pytest.mark.parametrize("aliases", [None, float("nan"), "", " ; ", "UNKNOWN", "KRAS"])
def test_missing_or_nonmatching_aliases_do_not_validate_a_retired_symbol(source_table, aliases):
    frame, accessor = source_table
    frame["aliases"] = aliases
    with pytest.raises(ValueError, match="no recorded alias supports this Ensembl ID"):
        accessor()


def test_known_conflicting_primary_symbol_cannot_be_overridden_by_an_alias(source_table):
    frame, accessor = source_table
    frame["gene_symbol"] = "KRAS"
    frame["aliases"] = "TP53"
    with pytest.raises(ValueError, match="symbol/Ensembl mismatch.*" + KRAS):
        accessor()


@pytest.mark.parametrize("symbol", [None, "", " "])
def test_alias_does_not_substitute_for_a_missing_source_symbol(source_table, symbol):
    frame, accessor = source_table
    frame["gene_symbol"] = symbol
    with pytest.raises(ValueError, match="missing its source gene symbol"):
        accessor()


def test_direct_match_does_not_need_alias_fallback(source_table):
    frame, accessor = source_table
    frame["gene_symbol"] = "TP53"
    frame["aliases"] = "KRAS"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        pd.testing.assert_frame_equal(accessor(), frame)
    assert not caught


def test_shared_alias_keeps_all_candidate_ids_and_uses_the_recorded_id(source_table, monkeypatch):
    frame, accessor = source_table
    frame["aliases"] = "SHARED_ALIAS"
    # The correct ID is deliberately second; a first-hit symbol map would lose it.
    monkeypatch.setattr(
        drivers,
        "canonical_gene_space",
        lambda: pd.DataFrame(
            {"symbol": ["SHARED_ALIAS", "SHARED_ALIAS"], "ensembl_gene_id": [KRAS, TP53]}
        ),
    )
    with pytest.warns(drivers.DriverGeneIdentityWarning) as caught:
        pd.testing.assert_frame_equal(accessor(), frame)
    assert TP53 in str(caught[0].message)
    assert KRAS in str(caught[0].message)


def test_ambiguous_primary_symbol_is_warned_and_cannot_be_queried_without_id(
    source_table, monkeypatch
):
    frame, accessor = source_table
    frame["gene_symbol"] = "SHARED_SYMBOL"
    monkeypatch.setattr(
        drivers,
        "canonical_gene_space",
        lambda: pd.DataFrame(
            {"symbol": ["SHARED_SYMBOL", "SHARED_SYMBOL"], "ensembl_gene_id": [KRAS, TP53]}
        ),
    )
    with pytest.warns(drivers.DriverGeneIdentityWarning, match="ambiguous source symbol"):
        pd.testing.assert_frame_equal(accessor(), frame)
    query = (
        drivers.driver_gene_evidence
        if accessor is drivers.driver_gene_evidence_df
        else drivers.driver_variant_evidence
    )
    with pytest.raises(ValueError, match=r"Ambiguous driver gene symbol.*use an Ensembl gene ID"):
        query("SHARED_SYMBOL")
    pd.testing.assert_frame_equal(query(TP53), frame.reset_index(drop=True))


def test_alias_column_is_required_for_auditable_fallback(source_table):
    frame, accessor = source_table
    frame.drop(columns="aliases", inplace=True)
    with pytest.raises(ValueError, match=r"missing required columns.*aliases"):
        accessor()

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Driver identities remain source-anchored without optional genome annotations."""

import warnings

import pandas as pd
import pytest

from oncoref import drivers, genome
from oncoref.gene_ids import canonical_gene_space, ensembl_id_aliases, entrez_gene_mappings
from oncoref.load_dataset import get_data

TP53 = "ENSG00000141510"
KRAS = "ENSG00000133703"
ZNG1C = "ENSG00000196873"


@pytest.fixture(autouse=True)
def clear_driver_caches():
    for cache in (
        drivers._driver_symbol_index,
        drivers._driver_source_symbol_index,
        drivers._driver_gene_evidence_frame,
        drivers._driver_variant_evidence_frame,
    ):
        cache.cache_clear()
    yield
    for cache in (
        drivers._driver_symbol_index,
        drivers._driver_source_symbol_index,
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
    monkeypatch.setattr(
        drivers,
        "get_data",
        lambda name, **kwargs: frame if name == dataset else get_data(name, **kwargs),
    )
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


@pytest.mark.parametrize("kind", ["gene", "variant"])
def test_every_shipped_driver_identifier_roundtrips_without_genome(kind, monkeypatch):
    def unexpected_genome_lookup(*args, **kwargs):
        raise AssertionError("driver queries must use bundled mappings")

    monkeypatch.setattr(genome, "canonical_gene_id_and_name", unexpected_genome_lookup)
    frame = get_data(f"driver-{kind}-evidence")
    query = getattr(drivers, f"driver_{kind}_evidence")
    symbols = canonical_gene_space().set_index("ensembl_gene_id")["symbol"].to_dict()
    expected = frame.groupby("ensembl_gene_id")["record_id"].agg(set).to_dict()
    probes = []
    for row in (
        frame[["gene_symbol", "aliases", "ensembl_gene_id"]]
        .drop_duplicates()
        .itertuples(index=False)
    ):
        tokens = [row.gene_symbol, symbols[row.ensembl_gene_id], row.ensembl_gene_id]
        if pd.notna(row.aliases):
            tokens.extend(alias.strip() for alias in row.aliases.split(";") if alias.strip())
        for token in tokens:
            probes.extend((form, row.ensembl_gene_id) for form in (token, f" {token} "))
        probes.append((row.ensembl_gene_id.lower(), row.ensembl_gene_id))
        probes.append((row.ensembl_gene_id + ".7", row.ensembl_gene_id))
    probes.extend(
        (alias, primary) for alias, primary in ensembl_id_aliases().items() if primary in expected
    )
    probes.extend(
        (str(row.entrez_id), row.canonical_ensembl_gene_id)
        for row in entrez_gene_mappings().itertuples(index=False)
        if row.canonical_ensembl_gene_id in expected
    )
    for token, gid in set(probes):
        assert set(query(token).record_id) == expected[gid], (kind, token, gid)


@pytest.mark.parametrize("kind", ["gene", "variant"])
def test_every_bundled_ncbi_alias_for_driver_genes_roundtrips(kind, monkeypatch):
    def unexpected_genome_lookup(*args, **kwargs):
        raise AssertionError("driver aliases must not fall back to installed annotations")

    monkeypatch.setattr(genome, "canonical_gene_id_and_name", unexpected_genome_lookup)
    frame = get_data(f"driver-{kind}-evidence")
    expected = frame.groupby("ensembl_gene_id")["record_id"].agg(set).to_dict()
    ids_by_symbol = canonical_gene_space().groupby("symbol")["ensembl_gene_id"].agg(set).to_dict()
    query = getattr(drivers, f"driver_{kind}_evidence")
    for row in get_data("ncbi-symbol-synonyms").itertuples(index=False):
        gene_ids = ids_by_symbol.get(row.official_symbol, set())
        matches = gene_ids.intersection(expected)
        if matches:
            assert len(gene_ids) == 1, (row.alias, row.official_symbol, gene_ids)
            gid = next(iter(matches))
            assert set(query(row.alias).record_id) == expected[gid], (kind, row.alias, gid)


@pytest.mark.parametrize("query", [drivers.driver_gene_evidence, drivers.driver_variant_evidence])
def test_exact_alias_case_wins_and_unanchored_case_fold_is_ambiguous(query):
    pd.testing.assert_frame_equal(query("emb"), query("XPO1"))
    assert query("EMB").empty  # EMB is a separate gene, not a driver in these sources.
    with pytest.raises(ValueError, match=r"Ambiguous driver gene symbol.*use an Ensembl gene ID"):
        query("Emb")
    for token, gid in (("alb", "ENSG00000163631"), ("met", "ENSG00000105976")):
        # These folds also collide with aliases of other genes in the snapshot.
        with pytest.raises(ValueError, match="Ambiguous driver gene symbol"):
            query(token)
        pd.testing.assert_frame_equal(query(token.upper()), query(gid))


@pytest.mark.parametrize(
    "old,current,gid",
    [
        ("CBWD3", "ZNG1C", "ENSG00000196873"),
        ("FAM46D", "TENT5D", "ENSG00000174016"),
        ("H3F3A", "H3-3A", "ENSG00000163041"),
        ("H3F3C", "H3-5", "ENSG00000188375"),
        ("HIST1H1C", "H1-2", "ENSG00000187837"),
        ("HIST1H1E", "H1-4", "ENSG00000168298"),
        ("RQCD1", "CNOT9", "ENSG00000144580"),
        ("TCEB1", "ELOC", "ENSG00000154582"),
        ("WHSC1", "NSD2", "ENSG00000109685"),
    ],
)
def test_all_nine_source_symbol_renames(old, current, gid):
    expected = drivers.driver_gene_evidence(gid)
    assert not expected.empty
    assert set(expected.gene_symbol) == {old}
    for token in (old, old.lower(), current, current.lower()):
        pd.testing.assert_frame_equal(drivers.driver_gene_evidence(token), expected)


@pytest.mark.parametrize("alias", ["OLD_SOURCE_NAME", "KRAS"])
def test_recorded_aliases_are_queryable_without_becoming_global_synonyms(source_table, alias):
    frame, accessor = source_table
    frame["gene_symbol"] = "TP53"
    frame["aliases"] = alias
    query = (
        drivers.driver_gene_evidence
        if accessor is drivers.driver_gene_evidence_df
        else drivers.driver_variant_evidence
    )
    if alias == "KRAS":
        # A source alias that collides with a bundled gene cannot silently select it.
        with pytest.raises(ValueError, match="Ambiguous driver gene symbol"):
            query(alias)
    else:
        pd.testing.assert_frame_equal(query(alias), frame.reset_index(drop=True))
        pd.testing.assert_frame_equal(query(f" {alias.lower()} "), frame.reset_index(drop=True))

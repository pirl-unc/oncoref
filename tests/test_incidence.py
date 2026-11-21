# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

import hashlib
import re

import pandas as pd
import pytest

from oncoref import incidence
from oncoref import load_dataset as ld


def test_derived_burden_cache_invalidates_on_clear(monkeypatch):
    # The lru_cached burden maps must be cleared by load_dataset._clear_cache(),
    # so swapping a bundled fixture is reflected (test-isolation guard).
    incidence._tissue_burden_map()  # populate the cache
    real = incidence.get_data

    def fake(name, *a, **k):
        if name == "tissue-burden-map":
            return pd.DataFrame(
                {"primary_tissue": ["unobtanium"], "burden_category": ["xyzzy"], "scope": ["solid"]}
            )
        return real(name, *a, **k)

    monkeypatch.setattr(incidence, "get_data", fake)
    ld._clear_cache()
    try:
        assert incidence._tissue_burden_map() == {"unobtanium": "xyzzy"}
    finally:
        monkeypatch.undo()
        ld._clear_cache()  # restore real data for other tests


def test_burden_map_nonempty_floats():
    mapping = incidence.cancer_burden()
    assert mapping
    assert all(isinstance(v, float) for v in mapping.values())


def test_burden_metric_validation():
    with pytest.raises(ValueError):
        incidence.cancer_burden(metric="not_a_metric")


def test_burden_all_metrics_resolve():
    for metric in incidence._BURDEN_METRICS:
        mapping = incidence.cancer_burden(metric=metric)
        assert mapping


def test_burden_table_exposes_source_provenance_schema():
    df = incidence.cancer_burden_df()
    expected = {
        "us_incidence_count",
        "us_incidence_total",
        "us_mortality_count",
        "us_mortality_total",
        "world_incidence_count",
        "world_incidence_total",
        "world_mortality_count",
        "world_mortality_total",
        "us_source_locator",
        "us_source_locator_status",
        "world_source_locator",
        "world_source_locator_status",
        "source_site_labels",
        "source_site_codes",
        "included_source_sites",
        "excluded_source_sites",
        "derivation_basis",
        "rounding_rule",
        "provenance_notes",
        "aggregation",
        "source_anchor",
    }
    assert expected <= set(df.columns)
    assert set(df["us_source_locator_status"]) <= {"not_extracted", "extracted"}
    assert set(df["world_source_locator_status"]) <= {
        "extracted",
        "extracted_partial_scope",
        "parent_extracted_subtype_split_approximated",
        "not_in_globocan2024_34site_table",
    }
    assert set(df["rounding_rule"]) <= {
        "round_half_even_1dp_over_all_cancers_incl_nmsc",
        "subtype_shares_sum_to_parent_1dp",
        "not_extracted",
    }
    assert set(df["derivation_basis"]) <= {
        "direct_source",
        "not_extracted",
        "sum_of_sites",
        "residual",
        "literature_approximation",
    }
    by_category = df.set_index("burden_category")
    assert by_category.loc["colorectal", "derivation_basis"] == "sum_of_sites"
    assert by_category.loc["other_and_unknown_primary", "derivation_basis"] == "residual"

    # GLOBOCAN 2024 raw-count provenance (#542): 29 rows carry extracted integer counts
    # that reproduce their published 1-decimal percentage shares over All cancers.
    extracted = df[df["world_source_locator_status"].isin({"extracted", "extracted_partial_scope"})]
    assert len(extracted) == 29
    assert extracted["world_incidence_count"].notna().all()
    assert extracted["world_mortality_count"].notna().all()
    assert (extracted["world_incidence_total"] == 20_636_778).all()
    assert (extracted["world_mortality_total"] == 9_762_507).all()
    for row in extracted.itertuples():
        assert round(
            100.0 * int(row.world_incidence_count) / int(row.world_incidence_total), 1
        ) == (row.world_incidence_pct)
        assert round(
            100.0 * int(row.world_mortality_count) / int(row.world_mortality_total), 1
        ) == (row.world_mortality_pct)

    assert int(by_category.loc["lung", "world_mortality_count"]) == 1_861_839
    assert by_category.loc["lung", "world_mortality_pct"] == 19.1
    assert int(by_category.loc["liver", "world_mortality_count"]) == 732_489
    assert by_category.loc["liver", "world_mortality_pct"] == 7.5
    assert int(by_category.loc["cervix", "world_mortality_count"]) == 279_583
    assert by_category.loc["cervix", "world_mortality_pct"] == 2.9

    # Site codes & inclusion/exclusion provenance (#542, #543) populated for all 37 rows.
    assert df["source_site_codes"].astype(str).str.strip().ne("").all()
    assert df["included_source_sites"].astype(str).str.strip().ne("").all()


def test_burden_provenance_is_a_lossless_pirlygenes_replacement():
    df = incidence.cancer_burden_df()
    assert not df["aggregation"].isna().any()
    assert not df["source_anchor"].isna().any()
    assert df["aggregation"].astype(str).str.strip().ne("").all()
    assert df["source_anchor"].astype(str).str.strip().ne("").all()

    # Frozen from pirlygenes' 37 populated burden-category aggregation values.
    ordered = df.sort_values("burden_category")
    aggregation_payload = "\n".join(
        f"{row.burden_category}\t{row.aggregation}" for row in ordered.itertuples()
    )
    assert hashlib.sha256(aggregation_payload.encode()).hexdigest() == (
        "0a1e3af0cf10da3fb43e180ca0d0c6d0f7cd081f53af3a4c58f3bcab48cbd589"
    )

    anchor_pattern = re.compile(r"(?:PMID:\d+|DOI:10\.\d{4,9}/\S+)")
    for source_anchor in df["source_anchor"].astype(str):
        assert all(anchor_pattern.fullmatch(token) for token in source_anchor.split(";"))

    by_category = df.set_index("burden_category")
    assert by_category.loc["adrenal", "source_anchor"] == "PMID:29977421;PMID:33112261"
    assert by_category.loc["soft_tissue_sarcoma", "source_anchor"] == (
        "PMID:38230766;PMID:39498209"
    )
    assert by_category.loc["kaposi_sarcoma", "source_anchor"] == (
        "PMID:42417444;DOI:10.3322/caac.70090;DOI:10.1016/S2214-109X(23)00349-2"
    )

    acs_rows = df["source"].astype(str).str.contains("ACS-CFF2024")
    globocan24_rows = df["source"].astype(str).str.contains("GLOBOCAN2024")
    globocan22_rows = df["source"].astype(str).str.contains("GLOBOCAN2022")
    assert globocan24_rows.sum() == 31
    assert df.loc[acs_rows, "source_anchor"].str.contains("PMID:38230766", regex=False).all()
    assert (
        df.loc[globocan24_rows, "source_anchor"]
        .str.contains("DOI:10.3322/caac.70090", regex=False)
        .all()
    )
    assert (
        df.loc[globocan22_rows, "source_anchor"]
        .str.contains("DOI:10.3322/caac.21834", regex=False)
        .all()
    )

    head_and_neck = by_category.loc["head_and_neck"]
    assert head_and_neck["us_mortality_pct"] == 2.6
    assert int(head_and_neck["us_mortality_count"]) == 16_110
    assert int(head_and_neck["world_mortality_count"]) == 475_065
    assert "C07-C08" in head_and_neck["included_source_sites"]
    assert "downstream 1.8% omitted larynx" in head_and_neck["provenance_notes"]


def test_burden_category_from_primary_tissue():
    assert incidence.burden_category("PRAD") == "prostate"
    assert incidence.burden_category("LUAD") == "lung"


def test_burden_category_sarcoma_bone_vs_soft():
    # Sarcoma family splits on primary_tissue: bone vs soft tissue.
    assert incidence.burden_category("SARC_OS") == "bone_and_joint"
    assert incidence.burden_category("SARC") == "soft_tissue_sarcoma"


def test_burden_category_override_table():
    # cancer-code-burden-map.csv carries the true ontology exceptions.
    overrides = incidence.cancer_code_burden_map()
    assert overrides
    for code, expected in overrides.items():
        assert incidence.burden_category(code) == expected


def test_burden_category_unknown_returns_none():
    assert incidence.burden_category("not_a_real_cancer") is None


def test_every_registry_code_resolves_to_a_burden_category():
    # Drift guard: the primary_tissue/family -> burden maps in incidence.py must
    # cover the registry's vocabulary. A new registry primary_tissue/family that
    # nothing maps would silently yield None here.
    from oncoref.cancer_types import cancer_type_registry

    registry = cancer_type_registry()
    unmapped = [
        code for code in registry["code"].astype(str) if incidence.burden_category(code) is None
    ]
    assert not unmapped, f"registry codes with no burden category: {unmapped}"


def test_family_burden_map_has_no_stale_families():
    # Drift guard for the data-driven family fallback: every family key must be a
    # current registry family. (This caught the stale 'cns' slug left after the
    # registry split cns -> cns-glial/-embryonal/... .)
    from oncoref.cancer_types import cancer_type_registry

    registry_families = set(cancer_type_registry()["family"].dropna().astype(str))
    stale = [f for f in incidence._family_burden_map() if f not in registry_families]
    assert not stale, f"family-burden-map keys not in the registry: {stale}"


def test_mixed_chol_and_salivary_burden_scope():
    # #543: scalar burden_category stays disjoint (no double-counted totals), while
    # burden_categories / burden_scope / cohorts_for_burden_category surface CHOL's
    # mixed intrahepatic (liver C22.1) + extrahepatic (gallbladder_biliary C24)
    # composition and explicit salivary C07-C08 inclusion in head_and_neck.
    assert incidence.burden_category("CHOL") == "gallbladder_biliary"
    assert incidence.burden_categories("CHOL") == ("gallbladder_biliary", "liver")
    chol_scope = incidence.burden_scope("CHOL")
    assert chol_scope["scope_status"] == "mixed_sites"
    assert chol_scope["expression_source_site_counts"] == {"liver": 32, "gallbladder_biliary": 4}

    assert incidence.burden_categories("BTC") == ("gallbladder_biliary", "liver")
    # Narrowed child GBC (gallbladder only) must NOT inherit BTC's mixed liver scope.
    assert incidence.burden_categories("GBC") == ("gallbladder_biliary",)
    assert incidence.burden_scope("GBC")["scope_status"] == "exact"

    # Salivary cohorts inherit SGC's explicit C07-C08 subsite inclusion in head_and_neck.
    for sal_code in ("SGC", "ADCC"):
        scope = incidence.burden_scope(sal_code)
        assert scope["primary_burden_category"] == "head_and_neck"
        assert scope["burden_categories"] == ("head_and_neck",)
        assert scope["scope_status"] == "included_subsite"
        assert "C07-C08" in scope["notes"]

    # Coverage queries include CHOL under liver by default, while disjoint partitioning excludes it.
    assert "CHOL" in incidence.cohorts_for_burden_category("liver", include_mixed=True)
    assert "CHOL" not in incidence.cohorts_for_burden_category("liver", include_mixed=False)
    assert "LIHC" in incidence.cohorts_for_burden_category("liver", include_mixed=False)
    assert "ADCC" in incidence.cohorts_for_burden_category("head_and_neck", include_mixed=False)

    # Disjointness & overlap checks for custom world vs US category sets:
    assert incidence.burden_categories_disjoint(["colorectal", "anus"], metric="us_mortality_pct")
    assert not incidence.burden_categories_disjoint(
        ["colorectal", "anus"], metric="world_mortality_pct"
    )
    assert incidence.burden_category_overlaps(
        ["colorectal", "anus", "lung"], metric="world_mortality_pct"
    ) == (("anus", "colorectal"),)
    assert not incidence.burden_categories_disjoint(
        ["adrenal", "other_and_unknown_primary"], metric="world_incidence_pct"
    )
    assert incidence.burden_categories_disjoint(
        ["lung", "breast", "colorectal", "adrenal", "liver", "head_and_neck"],
        metric="world_mortality_pct",
    )


def test_proteoform_and_expression_commonprefix_operates_on_symbol_strings():
    # #578: slash-joined gene symbol contraction is string prefixing, not path math.
    from oncoref.proteoforms import _contract_members

    assert _contract_members("XAGE1A/XAGE1B") == "XAGE1A/B"
    assert _contract_members("CGB3/CGB5/CGB8") == "CGB3/5/8"
    assert _contract_members("AKAP17A/AKAP17A") == "AKAP17A"

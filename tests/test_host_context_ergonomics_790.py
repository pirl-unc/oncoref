# Licensed under the Apache License, Version 2.0.
"""Regression tests for #790: HostContext.for_host, compartments, SEX_RESTRICTION_CODES, and pregnancy scope."""

from __future__ import annotations

import json
import warnings

import pytest

from oncoref import (
    antigen_sets,
    antigens,
    canine,
    cli,
    host_context,
    mouse,
    target_filters,
    target_selection,
)


def test_for_host_promoted_builder_and_non_deprecated_from_legacy_alias():
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        assert host_context.from_legacy is host_context.for_host
        assert host_context.for_host("all") is None
        assert host_context.HostContext.for_host("all") is None

        h1 = host_context.for_host("female-patient", lactating=True, species="dog")
        h2 = host_context.HostContext.for_host("female-patient", lactating=True, species="dog")
        h3 = host_context.from_legacy("female-patient", lactating=True, species="dog")
        assert h1 == h2 == h3
        assert h1 is not None
        assert h1.species == "dog"
        assert h1.protected == ("lactating_mammary",)
        assert "testis" in h1.absent and "prostate" in h1.absent


def test_public_cross_species_compartments_and_sex_restriction_codes():
    expected_codes = (
        "CTA-M",
        "CTA-F",
        "CTA-MF",
        "CTA-P",
        "RA-M-only",
        "RA-F-only",
        "RA-MF-only",
        "RA-P-only",
    )
    assert expected_codes == host_context.SEX_RESTRICTION_CODES
    assert expected_codes == antigens.SEX_RESTRICTION_CODES

    # Every non-empty sex_restriction_code across all 3 species catalogs belongs to SEX_RESTRICTION_CODES
    for df in (
        target_selection.results(),
        canine.target_selection(),
        mouse.target_selection(),
    ):
        observed = set(df.loc[df["sex_restriction_code"].ne(""), "sex_restriction_code"])
        assert observed <= set(host_context.SEX_RESTRICTION_CODES)

    assert host_context.compartments("human") == host_context.COMPARTMENTS
    mouse_comps = host_context.compartments("mouse")
    dog_comps = host_context.compartments("dog")
    assert "coagulating_gland" in mouse_comps
    assert "uterine_horn" in mouse_comps and "lactating_mammary" in mouse_comps
    assert "uterine_horn" in dog_comps and "lactating_mammary" in dog_comps
    assert "seminal vesicle" not in dog_comps


def test_pregnant_fail_closed_default_and_explicit_allow_incomplete_scope(capsys):
    ctas = target_selection.results(category="CTA")

    # 1. Default pregnant=True fails closed across all 3 species due to embryo scope_complete=False
    preg_default = host_context.for_host("female-patient", pregnant=True)
    assert preg_default is not None and preg_default.allow_incomplete_scope is False
    assessed_default = host_context.assess(ctas, preg_default)
    assert int(assessed_default["pass_host"].sum()) == 0
    assert (assessed_default["host_status"].isin(["excluded", "unresolved"])).all()
    assert antigens.table("CTA", sex="female", pregnant=True).empty
    assert antigens.table("CTA", species="dog", sex="female", pregnant=True).empty
    assert antigens.table("CTA", species="mouse", sex="female", pregnant=True).empty

    # 2. Maternal-fetal interface only ("placenta", "lactating_breast") has scope_complete=True
    maternal_interface = host_context.HostContext(protect=("placenta", "lactating_breast"))
    assessed_interface = host_context.assess(ctas, maternal_interface)
    assert int(assessed_interface["pass_host"].sum()) == 306

    # 3. Explicit allow_incomplete_scope=True evaluates measured early-embryo + placenta + lactating RNA (<8 TPM)
    # while still failing closed on >=8 TPM expression OR unmeasured coverage != "measured"
    preg_opt_in = host_context.HostContext.for_host(
        "female-patient", pregnant=True, allow_incomplete_scope=True
    )
    assert preg_opt_in is not None and preg_opt_in.allow_incomplete_scope is True
    assert any("allow_incomplete_scope=True" in note for note in preg_opt_in.assumptions)
    assessed_opt_in = host_context.assess(ctas, preg_opt_in)
    n_eligible = int(assessed_opt_in["pass_host"].sum())
    assert 100 < n_eligible < 306

    # Missing coverage still fails closed even when allow_incomplete_scope=True
    missing_res = host_context.assess(
        ctas.head(1),
        host_context.HostContext(protect=("embryo",), allow_incomplete_scope=True),
        evidence={
            ctas.iloc[0]["entity_id"]: {
                "embryo": {"coverage": "missing", "max_tpm": None, "scope_complete": False}
            }
        },
    ).iloc[0]
    assert missing_res["host_status"] == "unresolved"
    assert not bool(missing_res["pass_host"])

    # antigens.table, antigen_sets.select, explain, and CLI wire allow_incomplete_scope end-to-end
    fem_preg_tbl = antigens.table("CTA", sex="female", pregnant=True, allow_incomplete_scope=True)
    assert 0 < len(fem_preg_tbl) <= n_eligible
    assert (fem_preg_tbl["effective_tumor_fold"] >= 10.0).all()

    sel = antigen_sets.select(
        "CTA", host_profile="female-patient", pregnant=True, allow_incomplete_scope=True
    )
    assert len(sel.proteins()) > 0

    first_lbl = str(fem_preg_tbl.iloc[0]["protein_group_label"])
    exp_default = antigens.explain(first_lbl, sex="female", pregnant=True)
    exp_opt_in = antigens.explain(
        first_lbl, sex="female", pregnant=True, allow_incomplete_scope=True
    )
    assert exp_default["host"]["pass_host"] is False
    assert "embryo: sampled scope incomplete" in exp_default["host"]["reasons"]
    assert exp_opt_in["host"]["pass_host"] is True

    assert (
        cli.main(
            [
                "antigens",
                "list",
                "--sex",
                "female",
                "--pregnant",
                "--allow-incomplete-scope",
                "--format",
                "json",
            ]
        )
        == 0
    )
    cap = capsys.readouterr()
    cli_rows = json.loads(cap.out)
    assert len(cli_rows) == len(fem_preg_tbl)

    # Non-boolean allow_incomplete_scope fails closed
    with pytest.raises(ValueError, match="allow_incomplete_scope"):
        host_context.HostContext(allow_incomplete_scope="yes")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="allow_incomplete_scope"):
        target_filters.TargetFilter(allow_incomplete_scope=1)  # type: ignore[arg-type]

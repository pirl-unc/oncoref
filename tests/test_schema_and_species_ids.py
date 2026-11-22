# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

from __future__ import annotations

import warnings

import pandas as pd
import pytest

import oncoref
from oncoref import coverage, drivers, gene_ids, normalization, source_matrices


def test_is_ensembl_gene_id_supports_human_canine_and_murine_prefixes():
    assert gene_ids.is_ensembl_gene_id("ENSG00000141510.18")
    assert gene_ids.is_ensembl_gene_id("ENSCAFG00000012345.2")
    assert gene_ids.is_ensembl_gene_id("ENSMUSG00000051257.4")
    assert not gene_ids.is_ensembl_gene_id("TP53")
    assert not gene_ids.is_ensembl_gene_id("NY-ESO-1")
    assert not gene_ids.is_ensembl_gene_id("7157")

    assert gene_ids.is_ensembl_gene_id("ENSCAFG00000012345", species="dog")
    assert not gene_ids.is_ensembl_gene_id("ENSCAFG00000012345", species="human")
    assert gene_ids.is_ensembl_gene_id("ENSMUSG00000051257", species="mouse")
    assert gene_ids.display_gene_name("ENSCAFG00000099999") is None


def test_canonicalize_gene_columns_bridges_schema_dialects():
    lower = pd.DataFrame(
        {
            "ensembl_gene_id": ["ENSG00000141510"],
            "gene_symbol": ["TP53"],
            "Cancer_Type": ["LUAD"],
        }
    )
    expr_style = oncoref.canonicalize_gene_columns(
        lower, style="expression", include_cancer_code=True
    )
    assert list(expr_style.columns) == ["Ensembl_Gene_ID", "Symbol", "cancer_code"]

    back = oncoref.canonicalize_gene_columns(expr_style, style="lower")
    assert list(back.columns) == ["ensembl_gene_id", "symbol", "cancer_code"]


def test_coverage_panel_supports_non_human_ensembl_gene_ids(tmp_path):
    panel = tmp_path / "mouse_dog_panel.csv"
    panel.write_text("gene_id,gene_symbol\nENSMUSG00000051257.4,Trap1a\nENSCAFG00000012345,PRAME\n")
    label, ids = coverage.resolve_gene_set(str(panel))
    assert label == "mouse_dog_panel"
    assert ids == {"ENSMUSG00000051257", "ENSCAFG00000012345"}

    bare = tmp_path / "bare_mouse.txt"
    bare.write_text("ENSMUSG00000051257.4\nENSCAFG00000012345.1\n")
    _, bare_ids = coverage.resolve_gene_set(str(bare))
    assert bare_ids == {"ENSMUSG00000051257", "ENSCAFG00000012345"}


def test_driver_gene_evidence_supports_canonical_cancer_code_filter():
    codes = drivers.driver_cancer_codes()
    assert "LUAD" in codes and "BRCA" in codes

    luad = drivers.driver_gene_evidence(
        "KRAS", cancer_code="lung adenocarcinoma", include_cancer_code=True
    )
    assert len(luad) == 1
    assert luad.iloc[0]["cancer_code"] == "LUAD"
    assert luad.iloc[0]["source_scope"] == "LUAD"

    pancan = drivers.driver_gene_evidence("TP53", source_scope="PANCAN", include_cancer_code=True)
    assert len(pancan) == 1
    assert pd.isna(pancan.iloc[0]["cancer_code"])


def test_clean_tpm_records_unfilled_compartments_and_audits_sclc():
    assert "SCLC" in source_matrices.partial_technical_compartment_cohorts()
    gene_table = pd.DataFrame({"Ensembl_Gene_ID": ["ENSG00000141510", "ENSG00000108821"]})
    values = pd.DataFrame({"s1": [10.0, 20.0], "s2": [5.0, 15.0]})

    with pytest.warns(RuntimeWarning, match="cannot fill"):
        warned = normalization.clean_tpm(values, gene_table=gene_table)
    assert "unfilled_compartments" in warned.attrs

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        recorded = normalization.clean_tpm(
            values,
            gene_table=gene_table,
            on_empty_compartment="record",
        )
    assert set(recorded.attrs["unfilled_compartments"]) == {"ribosomal", "technical"}
    pd.testing.assert_frame_equal(warned, recorded)

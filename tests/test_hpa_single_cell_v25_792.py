# Licensed under the Apache License, Version 2.0.
"""Regression tests for #792: HPA v25.1 single-cell nCPM column support and provenance unit."""

from __future__ import annotations

import pytest

from oncoref import hpa, human_protein_atlas, normal_cell_data, reference_data


def _seed(cache_root, name, version, filename, text):
    path = cache_root / "sources" / name / version / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def test_hpa_single_cell_v25_1_ncpm_and_v23_ntpm_both_work(monkeypatch, tmp_path):
    monkeypatch.setenv("CANCERDATA_DATA_DIR", str(tmp_path))
    monkeypatch.setitem(reference_data.REFERENCE_SOURCES["hpa_single_cell"], "pins", {})
    hpa._hpa_single_cell_for_version.cache_clear()
    hpa.hpa_cell_type_expression.cache_clear()

    _seed(
        tmp_path,
        "hpa_single_cell",
        "v25.1",
        "rna_single_cell_type.tsv",
        "Gene\tGene name\tCell type\tnCPM\n"
        "ENSG00000000003\tTSPAN6\tadipocytes\t2.8\n"
        "ENSG00000000003\tTSPAN6\tspermatocytes\t19.4\n"
        "ENSG00000147381\tMAGEA4\tspermatocytes\t412.0\n",
    )
    _seed(
        tmp_path,
        "hpa_single_cell",
        "v23",
        "rna_single_cell_type.tsv",
        "Gene\tGene name\tCell type\tnTPM\n"
        "ENSG00000000003\tTSPAN6\tAdipocytes\t147.5\n"
        "ENSG00000147381\tMAGEA4\tSpermatocytes\t380.0\n",
    )

    try:
        # 1. Default version is v25.1 (nCPM), and reference-releases.json records unit == 'nCPM'
        assert reference_data.resolve_version("hpa_single_cell") == "v25.1"
        prov_v25 = reference_data.REFERENCE_SOURCES["hpa_single_cell"]
        assert (
            reference_data._release_manifest["sources"]["hpa_single_cell"]["versions"]["v25.1"][
                "unit"
            ]
            == "nCPM"
        )
        assert prov_v25["default_version"] == "v25.1"

        # 2. hpa_cell_type_expression() succeeds on default v25.1 (nCPM) & returns defensive copies
        wide_v25 = hpa.hpa_cell_type_expression()
        assert list(wide_v25.columns) == [
            "Ensembl_Gene_ID",
            "Symbol",
            "adipocytes",
            "spermatocytes",
        ]
        assert wide_v25.attrs["value_column"] == "nCPM"
        by_sym = wide_v25.set_index("Symbol")
        assert by_sym.loc["TSPAN6", "adipocytes"] == pytest.approx(2.8)
        assert by_sym.loc["MAGEA4", "adipocytes"] == pytest.approx(0.0)
        assert by_sym.loc["MAGEA4", "spermatocytes"] == pytest.approx(412.0)

        # Mutating returned frame does not corrupt cache
        wide_v25.loc[:, "Symbol"] = "MUTATED"
        assert "MAGEA4" in set(human_protein_atlas.cell_type_expression()["Symbol"])

        # 3. Historical v23 (nTPM) also works via explicit version="v23"
        wide_v23 = hpa.hpa_cell_type_expression("v23")
        assert list(wide_v23.columns) == [
            "Ensembl_Gene_ID",
            "Symbol",
            "Adipocytes",
            "Spermatocytes",
        ]
        assert wide_v23.attrs["value_column"] == "nTPM"
        assert wide_v23.set_index("Symbol").loc["TSPAN6", "Adipocytes"] == pytest.approx(147.5)

        # 4. gene_cell_type_ntpm() works on both default v25.1 (nCPM) and explicit v23 (nTPM)
        sc_v25 = hpa.gene_cell_type_ntpm("ENSG00000000003.15")
        assert sc_v25 == {"adipocytes": 2.8, "spermatocytes": 19.4}
        sc_v23 = hpa.gene_cell_type_ntpm("ENSG00000000003", version="v23")
        assert sc_v23 == {"adipocytes": 147.5}

        # 5. normal_cell_data.observations("hpa_single_cell") also preserves native unit='nCPM' on v25.1
        obs_v25 = normal_cell_data.observations("hpa_single_cell", genes=["TSPAN6"])
        assert set(obs_v25["unit"]) == {"nCPM"}
        assert obs_v25.set_index("cell_type").loc["adipocytes", "value"] == pytest.approx(2.8)
    finally:
        hpa._hpa_single_cell_for_version.cache_clear()
        hpa.hpa_cell_type_expression.cache_clear()

import numpy as np
import pandas as pd

from oncoref import library_qc as qc


def test_library_screen_flows_through_public_qc_and_respects_scale():
    from oncoref import expression

    panel = qc.marker_panels().query("panel == 'testis'")
    raw = panel[["gene_id", "symbol"]].rename(
        columns={"gene_id": "Ensembl_Gene_ID", "symbol": "Symbol"}
    )
    raw["sample"] = panel.reference_tpm * 0.002
    kwargs = {
        "min_detected_genes": 0,
        "min_housekeeping_detected": 0,
        "max_top_gene_fraction": 1,
        "max_top10_gene_fraction": 1,
    }
    result = expression.sample_expression_qc_from_matrix(
        raw,
        source_metadata={"source_scale_class": "linear_rnaseq_tpm", "linear_tpm_comparable": True},
        **kwargs,
    )
    assert result.iloc[0].sample_qc_status == "fail"
    assert "normal_testis_tissue" in result.iloc[0].sample_qc_reasons
    proxy = expression.sample_expression_qc_from_matrix(
        raw,
        source_metadata={"source_scale_class": "microarray_proxy", "linear_tpm_comparable": False},
        **kwargs,
    )
    assert "normal_testis_tissue" not in proxy.iloc[0].sample_qc_reasons


def test_coordinated_testis_program_and_missing_marker():
    panel = qc.marker_panels().query("panel == 'testis'")
    raw = panel.rename(columns={"gene_id": "Ensembl_Gene_ID", "symbol": "Symbol"})
    raw["coordinated"] = raw.reference_tpm * 0.001
    raw["partial"] = raw.coordinated
    raw.loc[raw.index[0], "partial"] = 0
    signatures = qc.signatures(raw, ["coordinated", "partial"])
    assert "normal_testis_tissue" in qc.reasons(signatures.loc["coordinated"])
    assert not qc.reasons(signatures.loc["partial"])
    missing = qc.signatures(raw.iloc[1:], ["coordinated"])
    assert np.isnan(missing.iloc[0].testis_tissue_content)
    assert not qc.reasons(missing.iloc[0])


def test_background_boundary_and_minimum_coverage():
    row = pd.Series(
        {
            "olfactory_receptors_measured": 100,
            "olfactory_share_ge5": 0.2,
            "testis_tissue_content": np.nan,
            "unrelated_tissues": "",
        }
    )
    assert not qc.reasons(row)
    row.olfactory_share_ge5 = 0.21
    assert qc.reasons(row) == ["broad_background"]
    row.olfactory_receptors_measured = 99
    row.unrelated_tissues = "liver;pancreas;brain;kidney"
    assert qc.reasons(row) == ["several_unrelated_tissues"]

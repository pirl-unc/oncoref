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
    assert result.iloc[0].sample_qc_status == "pass"
    assert result.iloc[0].sample_qc_reasons == ""
    assert result.iloc[0].sample_qc_notes == "R"
    assert result.iloc[0].passes_expression_qc
    assert result.iloc[0].recommended_for_absolute_tpm_floor
    proxy = expression.sample_expression_qc_from_matrix(
        raw,
        source_metadata={"source_scale_class": "microarray_proxy", "linear_tpm_comparable": False},
        **kwargs,
    )
    assert proxy.iloc[0].sample_qc_notes == ""
    # An advisory reproductive note does not override an independent QC failure.
    failed = expression.sample_expression_qc_from_matrix(
        raw,
        source_metadata={"source_scale_class": "linear_rnaseq_tpm", "linear_tpm_comparable": True},
        **(kwargs | {"min_detected_genes": 100}),
    ).iloc[0]
    assert failed.sample_qc_status == "fail"
    assert failed.sample_qc_reasons == "low_detected_genes"
    assert failed.sample_qc_notes == "R"


def test_coordinated_testis_program_and_missing_marker():
    panel = qc.marker_panels().query("panel == 'testis'")
    raw = panel.rename(columns={"gene_id": "Ensembl_Gene_ID", "symbol": "Symbol"})
    raw["coordinated"] = raw.reference_tpm * 0.001
    raw["partial"] = raw.coordinated
    raw.loc[raw.index[0], "partial"] = 0
    signatures = qc.signatures(raw, ["coordinated", "partial"])
    assert not qc.reasons(signatures.loc["coordinated"])
    assert qc.notes(signatures.loc["coordinated"]) == ["R"]
    assert not qc.reasons(signatures.loc["partial"])
    assert not qc.notes(signatures.loc["partial"])
    missing = qc.signatures(raw.iloc[1:], ["coordinated"])
    assert np.isnan(missing.iloc[0].testis_tissue_content)
    assert not qc.reasons(missing.iloc[0])
    assert not qc.notes(missing.iloc[0])


def test_reproductive_notes_survive_strict_pass_selection(monkeypatch):
    from oncoref import expression

    panel = qc.marker_panels().query("panel == 'testis'")
    raw = panel[["gene_id", "symbol"]].rename(
        columns={"gene_id": "Ensembl_Gene_ID", "symbol": "Symbol"}
    )
    raw["noted"] = panel.reference_tpm * 0.002
    raw["failed"] = 0.0
    metrics = expression.sample_expression_qc_from_matrix(
        raw,
        source_metadata={"source_scale_class": "linear_rnaseq_tpm", "linear_tpm_comparable": True},
        min_detected_genes=1,
        min_housekeeping_detected=0,
        max_top_gene_fraction=1,
        max_top10_gene_fraction=1,
    )
    monkeypatch.setattr(expression, "sample_expression_qc", lambda *args, **kwargs: metrics)
    for mode in ("pass", "pass_or_warn"):
        selected = expression._apply_sample_qc_filter(raw, "TGCT", sample_qc=mode, auto_fetch=False)
        assert "noted" in selected
        assert "failed" not in selected


def test_background_boundary_and_minimum_coverage():
    row = pd.Series(
        {
            "olfactory_receptors_measured": 100,
            "olfactory_share_ge5": 0.025,
            "testis_tissue_content": np.nan,
            "unrelated_tissues": "",
        }
    )
    assert not qc.reasons(row)
    row.olfactory_share_ge5 = 0.026
    assert qc.reasons(row) == ["broad_background"]
    row.olfactory_receptors_measured = 99
    row.unrelated_tissues = "liver;pancreas;brain;kidney"
    assert qc.reasons(row) == ["several_unrelated_tissues"]

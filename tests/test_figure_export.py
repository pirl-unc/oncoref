# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

from __future__ import annotations

import matplotlib

matplotlib.use("Agg", force=False)
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.patches import FancyArrowPatch

import oncoref
from oncoref import cli, coverage, cta_curation_plots, figure_export, plots


def test_top_level_plotting_and_export_submodules_are_accessible():
    assert oncoref.figure_export is figure_export
    assert oncoref.plots is plots
    assert oncoref.cta_curation_plots is cta_curation_plots


def test_save_figure_writes_300_dpi_png_and_vector_pdf_with_collision_free_labels(tmp_path):
    fig, ax = plt.subplots(figsize=(6, 4))
    xs = [1.0, 1.01, 1.02, 2.5]
    ys = [3.0, 3.01, 2.99, 4.0]
    ax.scatter(xs, ys, s=[80, 140, 80, 90], linewidths=1.2, edgecolors="black")
    texts = [
        ax.text(x, y, name, fontsize=9 if i == 1 else 7)
        for i, (x, y, name) in enumerate(zip(xs, ys, ["LUAD", "BRCA_TNBC", "LUSC", "SKCM"]))
    ]
    ax.set_xlim(0.5, 3.2)
    ax.set_ylim(2.4, 4.6)
    figure_export.queue_scatter_labels(texts, ax)

    out_png = tmp_path / "scatter.png"
    figure_export.save_figure(fig, out_png, dpi=150)  # clamped up to >=300 dpi
    out_pdf = tmp_path / "scatter.pdf"

    assert out_png.is_file() and out_png.stat().st_size > 0
    assert out_pdf.is_file() and out_pdf.stat().st_size > 0

    renderer = fig.canvas.get_renderer()
    bounds = ax.get_window_extent()
    boxes = [t.get_window_extent(renderer) for t in texts]
    for i, a in enumerate(boxes):
        assert a.x0 >= bounds.x0 and a.x1 <= bounds.x1
        assert a.y0 >= bounds.y0 and a.y1 <= bounds.y1
        for b in boxes[i + 1 :]:
            overlaps = a.x0 < b.x1 and a.x1 > b.x0 and a.y0 < b.y1 and a.y1 > b.y0
            assert not overlaps
    # Crowded cluster produces at least one data-coordinate leader line.
    arrows = [p for p in ax.patches if isinstance(p, FancyArrowPatch)]
    assert len(arrows) >= 1
    plt.close(fig)


def test_plots_save_emits_matching_pdf_and_repel_falls_back_to_queue(tmp_path, monkeypatch):
    import sys

    monkeypatch.setitem(sys.modules, "adjustText", None)
    out = tmp_path / "apd1.png"
    fig = plots.apd1_vs_tmb(annotate=True, save=str(out))
    assert out.is_file()
    assert out.with_suffix(".pdf").is_file()
    plt.close(fig)


def test_threshold_contract_and_percentile_patient_coverage(monkeypatch, tmp_path):
    t_tpm = coverage.Threshold("tpm", 25)
    t_pct = coverage.Threshold("pctile", 95)
    assert t_tpm.slug == "t25" and t_tpm.count_col == "n_gt25"
    assert t_pct.kind == "percentile" and t_pct.slug == "p95" and t_pct.count_col == "n_p95"
    assert list(t_tpm.compare([25.0, 25.1])) == [False, True]
    assert list(t_pct.compare([94.9, 95.0])) == [False, True]

    with pytest.raises(ValueError, match="threshold kind"):
        coverage.Threshold("bogus", 10)
    with pytest.raises(ValueError, match="percentile thresholds"):
        coverage.Threshold("percentile", 0)

    # 20-gene transcriptome so gene GA at 100 TPM in p0/p1 ranks above p95 (>=95).
    bg_ids = [f"ENSG_BG{i}" for i in range(18)]
    fixture = pd.DataFrame(
        {
            "Ensembl_Gene_ID": ["ENSG_A", "ENSG_B", *bg_ids],
            "Symbol": ["GA", "GB", *[f"BG{i}" for i in range(18)]],
            "p0": [100.0, 0.0, *range(1, 19)],
            "p1": [0.0, 100.0, *range(1, 19)],
            "p2": [0.0, 0.0, *range(1, 19)],
        }
    )
    monkeypatch.setattr(coverage, "per_sample_expression", lambda *a, **k: fixture.copy())
    monkeypatch.setattr(
        coverage,
        "_coverage_source_metadata",
        lambda codes, **kw: {
            str(c): {
                "source_cohort": "GEO_MICROARRAY",
                "source_type": "geo-microarray",
                "source_scale_class": "microarray_proxy",
                "linear_tpm_comparable": False,
                "normalization": "tpm_clean",
            }
            for c in codes
        },
    )
    panel = tmp_path / "panel.csv"
    panel.write_text("Ensembl_Gene_ID\nENSG_A\nENSG_B\n")

    # Non-comparable microarray source rejects forced TPM mode.
    with pytest.raises(ValueError, match="non-comparable source scale"):
        coverage.patient_coverage(
            str(panel), codes=["FAKE"], threshold_mode="tpm", proteoform=False
        )

    # Auto mode selects percentile and writes p95 artifacts + PDFs.
    res = coverage.render(
        str(panel),
        codes=["FAKE"],
        out_dir=tmp_path / "pct_out",
        threshold_mode="auto",
        proteoform=False,
    )
    assert res["threshold_mode"] == "percentile"
    assert res["threshold"] == 95
    assert "n_p95" in res["counts"].columns and "pct_p95" in res["counts"].columns
    stacked_png = tmp_path / "pct_out" / "panel_stacked_coverage_p95.png"
    curves_png = tmp_path / "pct_out" / "panel_coverage_curves_p95.png"
    assert stacked_png.is_file() and stacked_png.with_suffix(".pdf").is_file()
    assert curves_png.is_file() and curves_png.with_suffix(".pdf").is_file()


def test_covering_set_parity_and_cli(monkeypatch, tmp_path):
    weights, unmapped = coverage.burden_weights(["LUAD", "LUSC", "UNMAPPED_SYNTHETIC"])
    assert unmapped == ["UNMAPPED_SYNTHETIC"]
    assert weights["UNMAPPED_SYNTHETIC"] == 0.0
    assert weights["LUAD"] > 0 and weights["LUAD"] == pytest.approx(weights["LUSC"])

    fake_ref = pd.DataFrame(
        [
            {
                "cancer_code": "LUAD",
                "source_cohort": "SRC1",
                "n_samples": 50,
                "Ensembl_Gene_ID": "ENSG00000147381",
                "Symbol": "MAGEA4",
                "expression": 40.0,
                "q3": 80.0,
            },
            {
                "cancer_code": "SKCM",
                "source_cohort": "SRC1",
                "n_samples": 50,
                "Ensembl_Gene_ID": "ENSG00000184033",
                "Symbol": "CTAG1A",
                "expression": 50.0,
                "q3": 90.0,
            },
            {
                "cancer_code": "PRAD",
                "source_cohort": "SRC1",
                "n_samples": 50,
                "Ensembl_Gene_ID": "ENSG00000147381",
                "Symbol": "MAGEA4",
                "expression": 2.0,
                "q3": 5.0,
            },
        ]
    )
    panel = tmp_path / "mini_cta.csv"
    panel.write_text("Ensembl_Gene_ID\nENSG00000147381\nENSG00000184033\n")

    res = coverage.render_covering_set(
        str(panel),
        stat="q3",
        threshold=30.0,
        out_dir=tmp_path / "covset",
        df=fake_ref,
    )
    assert res["n_cancer_types"] == 3
    assert res["n_coverable"] == 2
    # Denominator is all 3 eligible cancer types (66.67%), never 2/2 = 100%.
    assert res["steps"]["cum_pct_cancer_types"].iloc[-1] == pytest.approx(200.0 / 3.0)
    for p in res["paths"].values():
        assert pd.notna(p) and (tmp_path / "covset" / p.split("/")[-1]).is_file()

    monkeypatch.setattr(coverage, "cancer_reference_expression", lambda **kw: fake_ref.copy())
    cli_out = tmp_path / "cli_covset"
    rc = cli.main(
        [
            "plot",
            "covering-set",
            "--gene-set",
            str(panel),
            "--stat",
            "q3",
            "--threshold",
            "30",
            "--cohort",
            "LUAD",
            "--cohort",
            "SKCM",
            "--out",
            str(cli_out),
        ]
    )
    assert rc == 0
    assert (cli_out / "mini_cta_covering_set_q3.csv").is_file()
    assert (cli_out / "mini_cta_covering_set_q3.png").is_file()
    assert (cli_out / "mini_cta_covering_set_q3.pdf").is_file()


def test_cta_curation_stage_definitions_attrition_and_active_policy_exports():
    assert cta_curation_plots.STAGES["nominated"] == "source union"
    assert cta_curation_plots.STAGES["family_eligible"] == "Histone / alpha-tubulin exclusion"
    assert (
        cta_curation_plots.STAGES["hpa_restriction"]
        == "Adaptive reproductive RNA + IHC restriction"
    )
    assert cta_curation_plots.STAGES["default_panel"] == "Reviewed default panel"

    defs = cta_curation_plots.stage_definitions("core")
    assert [d["stage_id"] for d in defs] == list(cta_curation_plots.STAGES)
    assert len(defs[0]["ihc_reproductive_tissues"]) > len(defs[0]["rna_reproductive_tissues"])

    for scope in ("core", "extended"):
        attr = cta_curation_plots.final_stage_attrition(tissue_scope=scope)
        stage_map = {label: n for label, n, _ in cta_curation_plots.stage_counts(scope)}
        hpa_n = stage_map[cta_curation_plots.STAGES["hpa_restriction"]]
        default_n = stage_map[cta_curation_plots.STAGES["default_panel"]]
        assert int(attr["genes"].sum()) == hpa_n
        assert int(attr.loc[attr["outcome"] == "retained_default", "genes"].iloc[0]) == default_n

    gates = cta_curation_plots.active_policy_gate_matrix()
    assert len(gates) == 19893
    assert int(gates["is_cta"].sum()) == 332
    assert set(gates["gate_cancer_enrichment_10x"]).issubset({"pass", "fail", "unassessed"})
    assert set(gates["gate_reproductive_enrichment_10x"]).issubset({"pass", "fail", "unassessed"})

    transitions = cta_curation_plots.active_policy_transitions()
    by_policy = transitions.groupby("from_policy")["protein_groups"].sum().to_dict()
    assert by_policy == {
        "oncoref.cta-hpa-v23-core": 19893,
        "oncoref.cta-selection.2026-11-08.1": 19893,
    }

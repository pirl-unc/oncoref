# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

import sys
from pathlib import Path
from types import ModuleType

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from oncoref import cli, data_bundle, plots

pytestmark = pytest.mark.xdist_group("plots")


@pytest.fixture(autouse=True)
def _close_figures_after_test():
    yield
    plt.close("all")


def test_apd1_vs_tmb_renders(tmp_path):
    out = tmp_path / "apd1_vs_tmb.png"
    fig = plots.apd1_vs_tmb(annotate=False, save=str(out))
    assert out.exists() and out.stat().st_size > 0
    assert fig is not None
    # The lookup includes means, unspecified summaries and approximations.
    assert fig.axes[0].get_xlabel() == "Tumor mutational burden estimate (mut/Mb, log scale)"
    assert plots._reference_metric_axis("tmb")[1] == "Tumor mutational burden estimate (mut/Mb)"


def test_scatter_label_adjustment_is_bounded_by_iterations_not_wall_clock(monkeypatch):
    # A time_lim makes label layout depend on machine load: the same command renders
    # a clean figure on an idle box and a pile of overlapping labels on a busy one.
    # The bound must be deterministic, so iterations, never seconds.
    observed = {}
    fake_adjust_text = ModuleType("adjustText")

    def adjust_text(texts, **kwargs):
        observed["texts"] = texts
        observed.update(kwargs)

    fake_adjust_text.adjust_text = adjust_text
    monkeypatch.setitem(sys.modules, "adjustText", fake_adjust_text)

    texts = [object()]
    axis = object()
    plots._repel_labels(axis, texts)

    assert observed["texts"] == texts
    assert observed["ax"] is axis
    assert "time_lim" not in observed
    assert observed["iter_lim"] > 0
    # Leader lines back to the anchor point, so a moved label stays traceable.
    assert observed["arrowprops"]["arrowstyle"] == "-"


def test_scatter_label_adjustment_passes_anchor_points_when_given(monkeypatch):
    observed = {}
    fake_adjust_text = ModuleType("adjustText")

    def adjust_text(texts, **kwargs):
        observed.update(kwargs)

    fake_adjust_text.adjust_text = adjust_text
    monkeypatch.setitem(sys.modules, "adjustText", fake_adjust_text)

    plots._repel_labels(object(), [object()], [1.0], [2.0])
    # Anchors let adjustText spring a label back toward its own point rather than
    # leaving it wherever repulsion pushed it.
    assert observed["x"] == [1.0] and observed["y"] == [2.0]


def test_apd1_vs_tmb_strict_pd1_filters_proxy_targets(tmp_path):
    fig = plots.apd1_vs_tmb(strict_pd1=True, save=str(tmp_path / "strict.png"))
    labels = {t.get_text() for t in fig.axes[0].texts}
    proxy_codes = set(
        plots.cancer_apd1_response_df().query("drug_target != 'PD-1'")["cancer_code"].astype(str)
    )
    assert labels.isdisjoint(proxy_codes)


def test_apd1_orr_bars_renders(tmp_path):
    out = tmp_path / "bars.png"
    plots.apd1_orr_bars(save=str(out))
    assert out.exists() and out.stat().st_size > 0


def test_incidence_vs_mortality_renders(tmp_path):
    out = tmp_path / "inc.png"
    plots.incidence_vs_mortality(region="world", save=str(out))
    assert out.exists() and out.stat().st_size > 0


def test_ici_orr_pooled_forest_renders(tmp_path):
    out = tmp_path / "forest.png"
    fig = plots.ici_orr_pooled_forest(save=str(out))
    assert out.exists() and out.stat().st_size > 0
    # one row per cancer type at its fallback regimen (every covered cancer)
    assert len(fig.axes[0].get_yticklabels()) > 50
    # pinning a regimen also works
    plots.ici_orr_pooled_forest(regimen="PD-1", save=str(tmp_path / "f2.png"))


def test_ici_orr_pooled_forest_keeps_adcc_trial_points_and_ci(monkeypatch):
    from oncoref import ici

    monkeypatch.setattr(ici, "cancer_ici_response", lambda: {"ADCC": 6.0})
    axis = plots.ici_orr_pooled_forest().axes[0]
    # Both directly reported combo trials plus the primary-only pooled diamond.
    assert sorted(float(c.get_offsets()[0, 0]) for c in axis.collections) == [5.0, 6.0, 6.2]
    pooled_ci = [line for line in axis.lines if line.get_linewidth() == 2.6]
    assert len(pooled_ci) == 1
    assert tuple(pooled_ci[0].get_xdata()) == ici._wilson_ci(2, 32)


def test_incidence_bad_region():
    with pytest.raises(ValueError, match="region must be"):
        plots.incidence_vs_mortality(region="moon")


def test_cli_plot(tmp_path):
    out = tmp_path / "cli.png"
    assert cli.main(["plot", "apd1-vs-tmb", "--out", str(out)]) == 0
    assert out.exists()


def test_cli_plot_burden_category_bars(tmp_path):
    out = tmp_path / "cats.png"
    assert cli.main(["plot", "burden-category-bars", "--region", "world", "--out", str(out)]) == 0
    assert out.exists()


def test_patient_coverage_renderer_writes_pirlygenes_style_artifacts(tmp_path, monkeypatch):
    from oncoref import coverage

    fixture = pd.DataFrame(
        {
            "Ensembl_Gene_ID": ["ENSG_A", "ENSG_B", "ENSG_C"],
            "Symbol": ["GA", "GB", "GC"],
            "p0": [10.0, 0.0, 0.0],
            "p1": [10.0, 0.0, 10.0],
            "p2": [0.0, 10.0, 0.0],
            "p3": [0.0, 0.0, 0.0],
        }
    )
    monkeypatch.setattr(coverage, "per_sample_expression", lambda *a, **k: fixture.copy())
    panel = tmp_path / "panel.csv"
    panel.write_text("Ensembl_Gene_ID\nENSG_A\nENSG_B\nENSG_C\n")

    result = coverage.render_patient_coverage(
        str(panel),
        cohorts=["LUAD"],
        threshold=5,
        thresholds=(5,),
        out_dir=tmp_path / "out",
        proteoform=False,
    )

    assert result["n_cohorts"] == 1
    for path in result["paths"].values():
        assert Path(path).exists()
    assert (tmp_path / "out" / "panel_patient_counts.csv").exists()
    assert (tmp_path / "out" / "panel_stacked_coverage_t5.png").exists()
    assert (tmp_path / "out" / "panel_coverage_curves_t5.png").exists()


def test_cli_plot_patient_coverage_delegates(monkeypatch, tmp_path, capsys):
    from oncoref import coverage

    captured = {}

    def fake(gene_set, **kwargs):
        captured["gene_set"] = gene_set
        captured.update(kwargs)
        out = tmp_path / "out"
        out.mkdir(exist_ok=True)
        counts = out / "counts.csv"
        counts.write_text("ok\n")
        return {"paths": {"counts_csv": str(counts)}, "label": "CTA", "n_cohorts": 2}

    monkeypatch.setattr(coverage, "render_patient_coverage", fake)
    rc = cli.main(
        [
            "plot",
            "patient-coverage",
            "--gene-set",
            "cta",
            "--codes",
            "LUAD,SKCM",
            "--threshold",
            "25",
            "--out",
            str(tmp_path / "out"),
        ]
    )
    assert rc == 0
    assert captured["gene_set"] == "cta"
    assert captured["cohorts"] == ["LUAD", "SKCM"]
    assert captured["threshold"] == 25
    assert "counts_csv" in capsys.readouterr().out


def test_cta_curation_source_counts_partition():
    from oncoref import cta_curation_plots as ccp

    rows = ccp._per_source_counts(ccp._evidence())
    assert rows
    for row in rows:
        assert row["kept_confident"] + row["kept_weak"] + row["excluded"] == row["total"]
        assert row["total"] > 0
    totals = [row["total"] for row in rows]
    assert totals == sorted(totals, reverse=True)


def test_cta_curation_tag_sets_cover_primary_sources():
    from oncoref import cta_curation_plots as ccp

    df = ccp._evidence()
    sets = ccp._tag_sets(df)
    assert set(sets) == set(ccp.PRIMARY_SOURCES)
    assert sets["CTpedia"]
    assert sets["CTexploreR"]
    assert sets["daSilva2017_protein"]

    # The da Silva primary source is the mass-spec-validated subset only. The broad
    # ``daSilva2017`` cross-reference tag (the full 1,103-gene 0.9-threshold set)
    # shares a prefix with it and must not leak into the source figures.
    tagged = {
        str(ensg)
        for ensg, raw in zip(df["Ensembl_Gene_ID"], df["source_databases"].fillna(""))
        if "daSilva2017_protein" in {t.strip() for t in str(raw).split(";")}
    }
    assert sets["daSilva2017_protein"] == tagged


def test_cta_curation_stage_counts_are_monotonic():
    from oncoref import cta_curation_plots as ccp

    stages = ccp.stage_counts()
    assert next(label for label, _, _ in stages) == "source union"
    remaining = [n for _, n, _ in stages]
    assert remaining == sorted(remaining, reverse=True)
    # Each stage's drop is exactly the step down from the previous stage.
    for (_, prev, _), (_, cur, dropped) in zip(stages, stages[1:]):
        assert prev - cur == dropped


def test_cta_curation_stage_funnel_lands_on_shipped_set():
    from oncoref import cta
    from oncoref import cta_curation_plots as ccp

    # The funnel's last stage must equal what cta_gene_ids() actually returns,
    # not the raw passes_filters column (which is one adjudication short).
    assert ccp.stage_counts()[-1][1] == len(cta.cta_gene_ids())


def test_cta_curation_renderer_writes_every_figure(tmp_path):
    from oncoref import cta_curation_plots as ccp

    result = ccp.render(out_dir=tmp_path)
    assert set(result["paths"]) == set(ccp.FILENAMES)
    assert result["n_genes"] > 0
    for path in result["paths"].values():
        assert path.exists() and path.stat().st_size > 0


def test_cli_plot_cta_curation_delegates(monkeypatch, tmp_path, capsys):
    from oncoref import cta_curation_plots

    captured = {}

    def fake(*, out_dir):
        captured["out_dir"] = out_dir
        out = tmp_path / "cur"
        out.mkdir(exist_ok=True)
        fig = out / "cta-source-venn.png"
        fig.write_text("ok\n")
        return {"n_genes": 3, "paths": {"source_venn": fig}}

    monkeypatch.setattr(cta_curation_plots, "render", fake)
    rc = cli.main(["plot", "cta-curation", "--out", str(tmp_path / "cur")])
    assert rc == 0
    assert captured["out_dir"] == str(tmp_path / "cur")
    assert "source_venn" in capsys.readouterr().out


def test_cli_plot_coverage_stacked_needs_codes(capsys):
    # The coverage plots require --codes; missing it is a clean error, not a crash.
    assert cli.main(["plot", "cta-coverage-stacked", "--out", "x.png"]) == 1


def test_cli_plot_threshold_tpm_is_opt_in(monkeypatch, tmp_path):
    # Without --threshold-tpm the CLI must NOT inject a value, so the patient heatmap keeps
    # its within-sample p95 default (it never silently reverts to a flat TPM cut); when the
    # flag is given it is forwarded.
    captured = {}

    def fake(*, save, **kwargs):
        captured.clear()
        captured.update(kwargs)
        open(save, "wb").close()

    monkeypatch.setattr(plots, "cta_patient_count_heatmap", fake)
    assert cli.main(["plot", "cta-patient-heatmap", "--out", str(tmp_path / "a.png")]) == 0
    assert "threshold_tpm" not in captured  # -> p95 default preserved
    assert (
        cli.main(
            [
                "plot",
                "cta-patient-heatmap",
                "--out",
                str(tmp_path / "b.png"),
                "--threshold-tpm",
                "25",
            ]
        )
        == 0
    )
    assert captured["threshold_tpm"] == 25.0


# ---- CTA expression heatmap (needs the expression bundle / percentile data) ----


def test_cta_expression_heatmap_renders(tmp_path):
    assert data_bundle.item_is_local("cancer-reference-expression-percentiles"), (
        "tests require a staged percentile artifact"
    )

    out = tmp_path / "cta.png"
    oncoref = __import__("oncoref")
    local = oncoref.locally_available_percentile_cohorts(include_recomputable=False)
    availability = oncoref.cancer_reference_expression_availability(local)
    cohorts = availability.loc[availability["available"], "cancer_code"].astype(str).tolist()[:6]
    assert cohorts
    fig = plots.cta_expression_heatmap(
        cohorts=cohorts, n_cohorts=4, n_ctas=8, proteoform=False, save=str(out)
    )
    assert out.exists() and out.stat().st_size > 0
    assert fig is not None


def test_cta_expression_heatmap_bad_stat():
    with pytest.raises(ValueError, match="stat must be one of"):
        plots.cta_expression_heatmap(stat="mean", cohorts=["X"])


def test_cta_expression_heatmap_proteoform_reads_collapsed_vectors(tmp_path, monkeypatch):
    # proteoform=True reads the proteoform-summed percentile vectors and labels CTA
    # columns by the proteoform symbol (NY-ESO-1). Stub the per-cohort vector.
    import pandas as pd

    from oncoref import cta
    from oncoref.proteoforms import gene_to_proteoform_id

    # Use two real CTAs and their real proteoform keys, so the heatmap's CTA-key filter
    # (derived from the same CTA set) keeps both rows. At least one is a collapsed group.
    grouped = next(g for g in sorted(cta.cta_gene_ids()) if "/" in gene_to_proteoform_id([g])[g])
    singleton = next(
        g for g in sorted(cta.cta_gene_ids()) if "/" not in gene_to_proteoform_id([g])[g]
    )
    keys = gene_to_proteoform_id([grouped, singleton])

    def fake_pct(code, *, as_tpm=True, proteoform=False):
        assert proteoform  # the heatmap must request the collapsed variant
        return pd.DataFrame(
            {
                "proteoform_key": [keys[grouped], keys[singleton]],
                "Ensembl_Gene_ID": [grouped, singleton],
                "Symbol": [keys[grouped], "OTHER_CTA"],
                "proteoform_members": ["CTAG1A/CTAG1B", "OTHER_CTA"],
                "p25": [5.0, 1.0],
                "p50": [40.0, 2.0],
                "p75": [80.0, 3.0],
            }
        )

    monkeypatch.setattr(plots, "cohort_gene_percentiles", fake_pct)
    monkeypatch.setattr(plots, "available_percentile_cohorts", lambda *, proteoform=False: ["LUAD"])
    out = tmp_path / "cta_pf.png"
    fig = plots.cta_expression_heatmap(proteoform=True, n_cohorts=2, n_ctas=4, save=str(out))
    assert out.exists() and fig is not None


def test_cta_expression_heatmap_proteoform_missing_bundle():
    with pytest.raises(ValueError, match="proteoform-summed percentile vector"):
        plots.cta_expression_heatmap(proteoform=True, cohorts=[])


def test_cta_expression_heatmap_no_cohorts():
    with pytest.raises(ValueError, match="no cohorts"):
        plots.cta_expression_heatmap(cohorts=[])


def test_cta_expression_heatmap_skips_unusable_cohorts():
    # A code with no percentile vector is skipped with a warning, not a crash; with
    # only unusable codes the matrix is empty -> clean ValueError.
    with (
        pytest.warns(UserWarning, match="without a percentile vector"),
        pytest.raises(ValueError, match="no CTA expression data"),
    ):
        plots.cta_expression_heatmap(cohorts=["NOT_A_REAL_COHORT"])


# ---- CTA addressable burden (P2) + 9mer scaffold (P3) ----


def test_cta_addressable_burden_renders(tmp_path, monkeypatch):
    # Hermetic: stub the within-sample prevalence so the test doesn't need the
    # bundle. Real codes map to a burden category + incidence, so bars render.
    monkeypatch.setattr(
        plots,
        "_cta_prevalence_by_cohort",
        lambda threshold: {"LUAD": 0.6, "SKCM": 0.4, "BRCA": 0.2},
    )
    out = tmp_path / "burden.png"
    fig = plots.cta_addressable_burden(n=10, save=str(out))
    assert out.exists() and out.stat().st_size > 0
    assert fig is not None


def test_cta_addressable_burden_labels_mortality_metric(tmp_path, monkeypatch):
    monkeypatch.setattr(
        plots,
        "_cta_prevalence_by_cohort",
        lambda threshold: {"LUAD": 0.6, "SKCM": 0.4, "BRCA": 0.2},
    )
    out = tmp_path / "burden_mortality.png"
    fig = plots.cta_addressable_burden(metric="world_mortality_pct", n=10, save=str(out))
    assert out.exists() and fig is not None
    ax = fig.axes[0]
    # The metric rides the axis label. The in-figure title is stripped by the
    # publication style (figure_style minimal mode) and restored when it is off.
    assert "WORLD mortality share" in ax.get_xlabel()
    assert ax.get_title() == ""


def test_cta_addressable_burden_no_within_sample(monkeypatch):
    monkeypatch.setattr(plots, "_cta_prevalence_by_cohort", lambda threshold: {})
    with pytest.raises(ValueError, match="no CTA prevalence available"):
        plots.cta_addressable_burden()


def test_cta_addressable_burden_no_mapped_cohorts(monkeypatch):
    # A code that resolves to no burden category yields no bars -> clean ValueError.
    monkeypatch.setattr(plots, "_cta_prevalence_by_cohort", lambda threshold: {"NOT_A_CODE": 0.5})
    with pytest.raises(ValueError, match="no cohort mapped"):
        plots.cta_addressable_burden()


def test_cta_specific_9mer_load_renders(tmp_path, monkeypatch):
    from oncoref import peptides

    monkeypatch.setattr(plots, "_cached_per_sample_cohorts", lambda: ["LUAD", "SKCM", "MM"])
    monkeypatch.setattr(plots, "cancer_tmb", lambda: {"LUAD": 6.9, "SKCM": 13.0, "MM": 1.5})
    monkeypatch.setattr(
        peptides,
        "cta_specific_9mer_load",
        lambda code, **k: {"LUAD": 120.0, "SKCM": 540.0, "MM": 60.0}[code],
    )
    out = tmp_path / "9mer.png"
    fig = plots.cta_specific_9mer_load(against="tmb", save=str(out))
    assert out.exists() and fig is not None


def test_cta_specific_9mer_load_collapses_crc_msi_source_scope(tmp_path, monkeypatch):
    from oncoref import peptides, source_matrices

    monkeypatch.setattr(plots, "_cached_per_sample_cohorts", lambda: ["COAD_MSI", "READ_MSI"])
    monkeypatch.setattr(plots, "cancer_tmb", lambda: {"CRC_MSI": 46.0})
    monkeypatch.setattr(
        peptides,
        "cta_specific_9mer_load",
        lambda code, **k: {"COAD_MSI": 10.0, "READ_MSI": 30.0}[code],
    )
    monkeypatch.setattr(
        source_matrices,
        "cohort_info",
        lambda code: {"n_samples": {"COAD_MSI": 3, "READ_MSI": 1}[code]},
    )

    out = tmp_path / "crc_msi_9mer.png"
    fig = plots.cta_specific_9mer_load(against="tmb", save=str(out))
    offsets = [tuple(coll.get_offsets()[0]) for coll in fig.axes[0].collections]
    assert len(offsets) == 1
    assert offsets[0][0] == pytest.approx(15.0)
    assert offsets[0][1] == pytest.approx(46.0)
    assert [t.get_text() for t in fig.axes[0].texts] == ["CRC_MSI"]


def test_cta_specific_9mer_load_bad_against():
    with pytest.raises(ValueError, match="against must be"):
        plots.cta_specific_9mer_load(against="nonsense")


def test_cta_specific_9mer_load_no_data(monkeypatch):
    monkeypatch.setattr(plots, "_cached_per_sample_cohorts", lambda: [])
    with pytest.raises(ValueError, match="no cohort with both"):
        plots.cta_specific_9mer_load()


# ---- per-patient plots (consume the per-sample matrices via coverage.py) ----


def test_cta_addressable_burden_per_sample_source(tmp_path, monkeypatch):
    # source="per_sample" pulls the faithful union from coverage; stub it hermetically.
    from oncoref import coverage

    monkeypatch.setattr(
        coverage,
        "addressable_fraction_by_cohort",
        lambda **k: __import__("pandas").Series({"LUAD": 0.95, "SKCM": 0.99, "BRCA": 0.5}),
    )
    out = tmp_path / "burden_faithful.png"
    fig = plots.cta_addressable_burden(source="per_sample", n=10, save=str(out))
    assert out.exists() and out.stat().st_size > 0
    assert fig is not None


def test_cta_addressable_burden_bad_source():
    with pytest.raises(ValueError, match="source must be"):
        plots.cta_addressable_burden(source="nonsense")


def test_cta_patient_count_heatmap_renders(tmp_path, monkeypatch):
    import pandas as pd

    def fake_ws(code, *, threshold, proteoform, scope):
        # within-sample p95 prevalence: two cohorts, three CTA proteoforms
        base = {"LUAD": [0.7, 0.3, 0.1], "SKCM": [0.9, 0.2, 0.5]}[code]
        return pd.DataFrame(
            {
                "proteoform_key": ["E1", "E2", "E3"],
                "Symbol": ["GA", "GB", "GC"],
                "frac_samples_top5pct": base,
                "n_samples": [100, 100, 100],
            }
        )

    from oncoref import proteoforms

    monkeypatch.setattr(
        plots,
        "locally_available_within_sample_cohorts",
        lambda **kwargs: ["LUAD", "SKCM"],
    )
    monkeypatch.setattr(plots, "within_sample_top_fraction", fake_ws)
    monkeypatch.setattr(plots, "cta_gene_ids", lambda: ["E1", "E2", "E3"])
    monkeypatch.setattr(proteoforms, "gene_to_proteoform_id", lambda ids: {i: i for i in ids})
    out = tmp_path / "patient_heatmap.png"
    fig = plots.cta_patient_count_heatmap(n_ctas=3, save=str(out))
    assert out.exists() and out.stat().st_size > 0
    assert fig is not None


def test_cta_patient_count_heatmap_no_cohorts(monkeypatch):
    monkeypatch.setattr(plots, "locally_available_within_sample_cohorts", lambda **kwargs: [])
    with pytest.raises(ValueError, match="no local proteoform within-sample"):
        plots.cta_patient_count_heatmap()


def test_cta_patient_count_heatmap_absolute_mode_uses_only_cached_matrices(monkeypatch):
    monkeypatch.setattr(plots, "_cached_per_sample_cohorts", lambda: [])
    monkeypatch.setattr(
        plots,
        "locally_available_within_sample_cohorts",
        lambda **kwargs: pytest.fail("absolute mode must not plan from summary shards"),
    )
    with pytest.raises(ValueError, match="no cached per-sample matrices"):
        plots.cta_patient_count_heatmap(threshold_tpm=50)


def test_cta_patient_count_heatmap_duplicate_symbols(tmp_path, monkeypatch):
    # Distinct proteoforms that happen to share a display Symbol must not crash the
    # cohort×CTA frame alignment — keying is on the unique proteoform_key.
    import pandas as pd

    def fake_ws(code, *, threshold, proteoform, scope):
        return pd.DataFrame(
            {
                "proteoform_key": ["K1", "K2", "K3"],  # unique keys
                "Symbol": ["GA", "GA", "GB"],  # GA duplicated as a display label
                "frac_samples_top5pct": [0.7, 0.3, 0.1],
                "n_samples": [100, 100, 100],
            }
        )

    from oncoref import proteoforms

    monkeypatch.setattr(
        plots,
        "locally_available_within_sample_cohorts",
        lambda **kwargs: ["LUAD", "SKCM"],
    )
    monkeypatch.setattr(plots, "within_sample_top_fraction", fake_ws)
    monkeypatch.setattr(plots, "cta_gene_ids", lambda: ["K1", "K2", "K3"])
    monkeypatch.setattr(proteoforms, "gene_to_proteoform_id", lambda ids: {i: i for i in ids})
    out = tmp_path / "dup.png"
    fig = plots.cta_patient_count_heatmap(save=str(out))  # must not raise
    assert out.exists() and fig is not None


def test_cta_coverage_curves_renders(tmp_path, monkeypatch):
    import pandas as pd

    from oncoref import coverage

    def fake_greedy(code, *, threshold_tpm, max_genes):
        return pd.DataFrame(
            {
                "rank": [1, 2, 3],
                "Ensembl_Gene_ID": ["E1", "E2", "E3"],
                "Symbol": ["GA", "GB", "GC"],
                "marginal_patients": [60, 20, 10],
                "marginal_fraction": [0.6, 0.2, 0.1],
                "cumulative_patients": [60, 80, 90],
                "cumulative_fraction": [0.6, 0.8, 0.9],
            }
        )

    monkeypatch.setattr(coverage, "greedy_coverage", fake_greedy)
    out = tmp_path / "curves.png"
    fig = plots.cta_coverage_curves(["LUAD", "SKCM"], save=str(out))
    assert out.exists() and out.stat().st_size > 0
    assert fig is not None


def test_cta_coverage_curves_empty_raises(monkeypatch):
    import pandas as pd

    from oncoref import coverage

    monkeypatch.setattr(coverage, "greedy_coverage", lambda *a, **k: pd.DataFrame())
    with pytest.raises(ValueError, match="no coverage curve"):
        plots.cta_coverage_curves(["LUAD"])


def _percentile_coverage_table():
    rows = []
    for code, final in (("LUAD", 0.8), ("SKCM", 1.0)):
        rows.extend(
            [
                {
                    "cancer_code": code,
                    "within_sample_percentile": 0.95,
                    "n_patients": 10,
                    "rank": 1,
                    "Ensembl_Gene_ID": "E1",
                    "Symbol": "MAGEA4",
                    "marginal_patients": 6,
                    "marginal_fraction": 0.6,
                    "cumulative_patients": 6,
                    "cumulative_fraction": 0.6,
                    "proteoform_key": "E1",
                    "proteoform_members": "MAGEA4",
                },
                {
                    "cancer_code": code,
                    "within_sample_percentile": 0.95,
                    "n_patients": 10,
                    "rank": 2,
                    "Ensembl_Gene_ID": "E2",
                    "Symbol": "CTAG1B",
                    "marginal_patients": round((final - 0.6) * 10),
                    "marginal_fraction": final - 0.6,
                    "cumulative_patients": round(final * 10),
                    "cumulative_fraction": final,
                    "proteoform_key": "E2",
                    "proteoform_members": "CTAG1B",
                },
            ]
        )
    return pd.DataFrame(rows)


def test_cta_within_sample_percentile_coverage_plots_render(tmp_path):
    coverage_table = _percentile_coverage_table()
    curves = tmp_path / "percentile-curves.png"
    stacked = tmp_path / "percentile-stacked.png"

    curve_figure = plots.cta_within_sample_percentile_coverage_curves(
        ["lung adenocarcinoma", "SKCM"], coverage_table=coverage_table, save=str(curves)
    )
    stacked_figure = plots.cta_within_sample_percentile_coverage_stacked_bars(
        ["LUAD", "SKCM"], coverage_table=coverage_table, save=str(stacked)
    )

    assert curves.stat().st_size > 0 and stacked.stat().st_size > 0
    assert "within-sample p95" in curve_figure.axes[0].get_ylabel()
    assert "within-sample p95" in stacked_figure.axes[0].get_xlabel()


def test_cta_within_sample_percentile_addressable_burden_uses_joint_union(tmp_path, monkeypatch):
    from oncoref import coverage

    supplied = _percentile_coverage_table()
    captured = {}

    def fake_fraction(cohorts, **kwargs):
        captured.update(kwargs)
        return pd.Series({"LUAD": 0.8, "SKCM": 1.0})

    monkeypatch.setattr(
        coverage, "within_sample_percentile_addressable_fraction_by_cohort", fake_fraction
    )
    out = tmp_path / "percentile-burden.png"
    figure = plots.cta_within_sample_percentile_addressable_burden(
        cohorts=["LUAD", "SKCM"], coverage_table=supplied, save=str(out)
    )

    assert captured["coverage"] is supplied
    assert captured["percentile"] == 0.95
    assert out.stat().st_size > 0 and figure is not None


def test_antigen_family_grouping():
    assert plots._antigen_family("MAGEA4") == "MAGE-A"
    assert plots._antigen_family("MAGEB2") == "MAGE-B"
    assert plots._antigen_family("CTAG1B") == "CTAG/NY-ESO"
    assert plots._antigen_family("LAGE1") == "CTAG/NY-ESO"  # same family as CTAG
    assert plots._antigen_family("SSX2") == "SSX"
    assert plots._antigen_family("CT83") == "other"


def test_cta_coverage_stacked_bars_renders(tmp_path, monkeypatch):
    import pandas as pd

    from oncoref import coverage

    def fake_greedy(code, *, threshold_tpm, max_genes):
        return pd.DataFrame(
            {
                "rank": [1, 2, 3],
                "Ensembl_Gene_ID": ["E1", "E2", "E3"],
                "Symbol": ["MAGEA4", "CTAG1B", "SSX2"],
                "marginal_patients": [60, 20, 10],
                "marginal_fraction": [0.6, 0.2, 0.1],
                "cumulative_patients": [60, 80, 90],
                "cumulative_fraction": [0.6, 0.8, 0.9],
            }
        )

    monkeypatch.setattr(coverage, "greedy_coverage", fake_greedy)
    out = tmp_path / "stacked.png"
    fig = plots.cta_coverage_stacked_bars(["LUAD", "SKCM"], save=str(out))
    assert out.exists() and out.stat().st_size > 0
    assert fig is not None


def test_cta_coverage_stacked_bars_empty_raises(monkeypatch):
    import pandas as pd

    from oncoref import coverage

    monkeypatch.setattr(coverage, "greedy_coverage", lambda *a, **k: pd.DataFrame())
    with pytest.raises(ValueError, match="no coverage to plot"):
        plots.cta_coverage_stacked_bars(["LUAD"])


def test_burden_category_bars_renders(tmp_path):
    out = tmp_path / "cats.png"
    fig = plots.burden_category_bars(region="us", n=8, save=str(out))
    assert out.exists() and out.stat().st_size > 0
    assert fig is not None


def test_burden_category_bars_bad_region():
    with pytest.raises(ValueError, match="region must be"):
        plots.burden_category_bars(region="moon")


def test_cta_burden_vs_response_renders(tmp_path, monkeypatch):
    from oncoref import coverage

    monkeypatch.setattr(plots, "_cached_per_sample_cohorts", lambda: ["LUAD", "SKCM", "MM"])
    monkeypatch.setattr(
        plots, "cancer_apd1_response", lambda: {"LUAD": 19.0, "SKCM": 42.0, "MM": 3.0}
    )
    monkeypatch.setattr(
        coverage,
        "mean_antigens_per_patient",
        lambda code, **k: {"LUAD": 1.5, "SKCM": 4.0, "MM": 8.0}[code],
    )
    out = tmp_path / "load.png"
    fig = plots.cta_burden_vs_response(against="apd1", save=str(out))
    assert out.exists() and fig is not None


def test_cta_burden_vs_response_supports_burden_axis(tmp_path, monkeypatch):
    from oncoref import coverage

    monkeypatch.setattr(plots, "_cached_per_sample_cohorts", lambda: ["LUAD", "SKCM"])
    monkeypatch.setattr(
        coverage,
        "mean_antigens_per_patient",
        lambda code, **k: {"LUAD": 1.5, "SKCM": 4.0}[code],
    )
    out = tmp_path / "burden_axis.png"
    fig = plots.cta_burden_vs_response(against="us_incidence", save=str(out))
    assert out.exists() and fig is not None
    assert fig.axes[0].get_ylabel() == "US incidence share (%)"


def test_cta_burden_vs_ici_includes_ici_only_anchor(tmp_path, monkeypatch):
    from oncoref import coverage

    assert "THYM" not in plots.cancer_apd1_response()
    monkeypatch.setattr(plots, "_cached_per_sample_cohorts", lambda: ["THYM"])
    monkeypatch.setattr(
        coverage,
        "mean_antigens_per_patient",
        lambda code, **k: {"THYM": 2.5}[code],
    )
    out = tmp_path / "ici_only.png"
    fig = plots.cta_burden_vs_response(against="ici", save=str(out))
    assert out.exists() and fig is not None
    assert [t.get_text() for t in fig.axes[0].texts] == ["THYM"]


def test_cta_burden_vs_burden_axis_keeps_direct_msi_cohorts(tmp_path, monkeypatch):
    from oncoref import coverage

    monkeypatch.setattr(plots, "_cached_per_sample_cohorts", lambda: ["COAD_MSI", "READ_MSI"])
    monkeypatch.setattr(
        coverage,
        "mean_antigens_per_patient",
        lambda code, **k: {"COAD_MSI": 10.0, "READ_MSI": 30.0}[code],
    )
    out = tmp_path / "burden_msi.png"
    fig = plots.cta_burden_vs_response(against="us_incidence", save=str(out))
    assert out.exists() and fig is not None
    assert {t.get_text() for t in fig.axes[0].texts} == {"COAD_MSI", "READ_MSI"}
    xs = sorted(float(coll.get_offsets()[0][0]) for coll in fig.axes[0].collections)
    assert xs == [10.0, 30.0]


def test_cta_specific_9mer_burden_axis_keeps_direct_msi_cohorts(tmp_path, monkeypatch):
    from oncoref import peptides

    monkeypatch.setattr(plots, "_cached_per_sample_cohorts", lambda: ["COAD_MSI", "READ_MSI"])
    monkeypatch.setattr(
        peptides,
        "cta_specific_9mer_load",
        lambda code, **k: {"COAD_MSI": 100.0, "READ_MSI": 300.0}[code],
    )
    out = tmp_path / "9mer_burden_msi.png"
    fig = plots.cta_specific_9mer_load(against="us_incidence", save=str(out))
    assert out.exists() and fig is not None
    assert {t.get_text() for t in fig.axes[0].texts} == {"COAD_MSI", "READ_MSI"}
    xs = sorted(float(coll.get_offsets()[0][0]) for coll in fig.axes[0].collections)
    assert xs == [100.0, 300.0]


def test_cta_burden_vs_response_bad_against():
    with pytest.raises(ValueError, match="against must be"):
        plots.cta_burden_vs_response(against="nonsense")


def test_cta_burden_vs_response_no_data(monkeypatch):
    monkeypatch.setattr(plots, "_cached_per_sample_cohorts", lambda: [])
    with pytest.raises(ValueError, match="no cohort with both"):
        plots.cta_burden_vs_response()


def test_apd1_response_signature_scatter_renders(tmp_path, monkeypatch):
    from oncoref import response_signatures as rs

    monkeypatch.setattr(plots, "_cached_per_sample_cohorts", lambda: ["LUAD", "SKCM", "MM"])
    monkeypatch.setattr(
        plots, "cancer_apd1_response", lambda: {"LUAD": 19.0, "SKCM": 42.0, "MM": 3.0}
    )
    monkeypatch.setattr(
        rs, "signature_score", lambda code, sig, **k: {"LUAD": 3.1, "SKCM": 2.9, "MM": 1.0}[code]
    )
    out = tmp_path / "sig.png"
    fig = plots.apd1_response_signature_scatter("t_cell_inflamed", save=str(out))
    assert out.exists() and fig is not None


def test_apd1_response_signature_scatter_no_data(monkeypatch):
    monkeypatch.setattr(plots, "_cached_per_sample_cohorts", lambda: [])
    with pytest.raises(ValueError, match="no cohort with both"):
        plots.apd1_response_signature_scatter("t_cell_inflamed")


def test_apd1_response_signature_bad_name():
    with pytest.raises(ValueError, match="unknown signature"):
        plots.apd1_response_signature_scatter("not_a_signature")


def test_regenerate_plots_runner_references_real_functions():
    # Guard: every (family, name, fn_attr, kwargs) job in the batch runner must
    # name a real callable on oncoref.plots — catches typos like a removed or
    # renamed figure before a full run does.
    import importlib.util
    from pathlib import Path

    from oncoref import plots

    runner = Path(__file__).resolve().parent.parent / "scripts" / "regenerate_plots.py"
    spec = importlib.util.spec_from_file_location("_regen_plots", runner)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    availability = {
        "per_sample": ("LUAD", "SKCM"),
        "percentile_gene": ("LUAD", "SKCM"),
        "percentile_proteoform": ("LUAD", "SKCM"),
        "within_sample": ("LUAD", "SKCM"),
    }
    jobs = mod._jobs(availability)
    assert jobs, "runner produced no jobs"
    names = {name for _, name, _, _ in jobs}
    assert {"apd1_vs_tmb_ici", "apd1_vs_tmb_strict_pd1"} <= names
    assert {
        "cta_expression_heatmap_q1",
        "cta_expression_heatmap_median",
        "cta_expression_heatmap_q3",
    } <= names
    assert {
        "cta_burden_vs_us_incidence_t25",
        "cta_burden_vs_world_mortality_t50",
        "cta_specific_9mer_load_vs_world_mortality_t50",
    } <= names
    assert {
        "cta_addressable_burden_per_sample_p90_world_mortality",
        "cta_addressable_burden_per_sample_t50_us_mortality",
    } <= names
    assert {
        "cta_coverage_curves_p90",
        "cta_coverage_curves_p95",
        "cta_coverage_stacked_bars_p90",
        "cta_coverage_stacked_bars_p95",
    } <= names
    assert {"cta_patient_count_heatmap_p90", "cta_patient_count_heatmap_t50"} <= names
    for family, name, fn_attr, kwargs in jobs:
        assert family and name and isinstance(kwargs, dict)
        assert callable(getattr(plots, fn_attr, None)), f"{fn_attr} is not a plots function"


def test_regenerate_plots_runner_computes_percentile_coverage_once(monkeypatch):
    import importlib.util

    runner = Path(__file__).resolve().parent.parent / "scripts" / "regenerate_plots.py"
    spec = importlib.util.spec_from_file_location("_regen_plots_shared", runner)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    availability = {
        "per_sample": ("LUAD", "SKCM"),
        "percentile_gene": (),
        "percentile_proteoform": (),
        "within_sample": (),
    }
    jobs = mod._jobs(availability)
    sentinel = _percentile_coverage_table()
    calls = []

    def fake_sweep(cohorts, **kwargs):
        calls.append((cohorts, kwargs))
        return sentinel

    monkeypatch.setattr(mod, "within_sample_percentile_coverage_sweep", fake_sweep)
    prepared, result = mod._attach_shared_percentile_coverage(jobs, availability)

    assert result is sentinel
    assert calls == [(("LUAD", "SKCM"), {"percentiles": (0.9, 0.95), "on_missing": "record"})]
    dependent = [
        kwargs for _, _, function, kwargs in prepared if function in mod._PERCENTILE_COVERAGE_PLOTS
    ]
    assert len(dependent) == 12
    assert all(kwargs["coverage_table"] is sentinel for kwargs in dependent)


def test_regenerate_plots_runner_records_exact_data_inventory():
    import importlib.util

    runner = Path(__file__).resolve().parent.parent / "scripts" / "regenerate_plots.py"
    spec = importlib.util.spec_from_file_location("_regen_plots_inventory", runner)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    availability = {
        "per_sample": ("LUAD", "SKCM"),
        "percentile_gene": ("LUAD",),
        "percentile_proteoform": ("SKCM",),
        "within_sample": ("LUAD", "SKCM"),
    }
    coverage = pd.DataFrame()
    coverage.attrs.update({"audited_cohorts": ["LUAD"], "missing_cohorts": {"SKCM": "not staged"}})

    text = "\n".join(mod._availability_index_lines(availability, coverage))

    assert "Cached per-sample matrices (2): `LUAD`, `SKCM`" in text
    assert "p90/p95 joint coverage cohorts audited (1): `LUAD`" in text
    assert "- `SKCM`: not staged" in text


def test_regenerate_plots_runner_writes_all_figures_pdf(tmp_path):
    import importlib.util
    from pathlib import Path

    runner = Path(__file__).resolve().parent.parent / "scripts" / "regenerate_plots.py"
    spec = importlib.util.spec_from_file_location("_regen_plots", runner)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    rels = ["first/one.png", "second/two.png"]
    for rel in rels:
        path = tmp_path / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        plt.imsave(path, [[[0.2, 0.4, 0.6], [0.9, 0.8, 0.1]]])

    pdf = mod._write_all_figures_pdf(tmp_path, rels)

    assert pdf == tmp_path / "oncoref-all-figures.pdf"
    assert pdf.exists()
    assert pdf.stat().st_size > 0


def test_regenerate_plots_runner_closes_each_returned_figure(tmp_path, monkeypatch):
    import importlib.util
    import sys

    runner = Path(__file__).resolve().parent.parent / "scripts" / "regenerate_plots.py"
    spec = importlib.util.spec_from_file_location("_regen_plots_memory", runner)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    def stub_plot(*, save, label):
        fig, ax = plt.subplots()
        ax.set_title(label)
        fig.savefig(save)
        return fig

    monkeypatch.setattr(mod.plots, "_memory_test_plot", stub_plot, raising=False)
    monkeypatch.setattr(
        mod,
        "_jobs",
        lambda availability: [
            ("memory", "first", "_memory_test_plot", {"label": "first"}),
            ("memory", "second", "_memory_test_plot", {"label": "second"}),
        ],
    )
    monkeypatch.setattr(
        mod,
        "_plot_data_availability",
        lambda: {
            "per_sample": (),
            "percentile_gene": (),
            "percentile_proteoform": (),
            "within_sample": (),
        },
    )
    monkeypatch.setattr(
        mod.cta_curation_plots,
        "render",
        lambda **kwargs: {"n_genes": 0, "paths": {}},
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["regenerate_plots.py", "--out-dir", str(tmp_path), "--no-timestamp"],
    )

    existing_figures = set(plt.get_fignums())
    assert mod.main() == 0
    assert set(plt.get_fignums()) == existing_figures
    assert (tmp_path / "memory" / "first.png").exists()
    assert (tmp_path / "memory" / "second.png").exists()
    assert (tmp_path / "oncoref-all-figures.pdf").exists()


def test_burden_weights_split_a_category_across_its_codes():
    # A burden category's incidence share is divided among the codes that map to
    # it, so a heavily sub-typed cancer is not counted once per subtype.
    from oncoref.incidence import burden_category, cancer_burden

    weights = plots._burden_weights(["LUAD", "LUSC", "SKCM"])
    assert set(weights.index) == {"LUAD", "LUSC", "SKCM"}
    assert (weights > 0).all()
    lung = cancer_burden().get(burden_category("LUAD"))
    assert weights["LUAD"] + weights["LUSC"] == pytest.approx(lung)


def test_greedy_cover_is_deterministic_and_monotonic():
    matrix = pd.DataFrame(
        {"A": [100.0, 100.0, 0.0], "B": [0.0, 100.0, 100.0], "C": [0.0, 0.0, 0.0]},
        index=["X", "Y", "Z"],
    )
    weights = pd.Series({"X": 1.0, "Y": 1.0, "Z": 5.0})
    steps, coverable = plots._greedy_cover(matrix, 30.0, weights)
    assert sorted(coverable) == ["X", "Y", "Z"]
    # B first: it reaches Z (weight 5) plus Y, beating A's X+Y.
    assert [s[0] for s in steps] == ["B", "A"]
    assert [s[3] for s in steps] == [2, 3]
    assert [s[2] for s in steps] == [6.0, 7.0]
    # A gene covering nothing is never selected.
    assert "C" not in [s[0] for s in steps]
    assert plots._greedy_cover(matrix, 30.0, weights)[0] == steps


def test_greedy_cover_ignores_genes_below_threshold():
    matrix = pd.DataFrame({"A": [10.0, 20.0]}, index=["X", "Y"])
    steps, coverable = plots._greedy_cover(matrix, 30.0, pd.Series({"X": 1.0, "Y": 1.0}))
    assert steps == [] and list(coverable) == []


def test_cta_covering_set_renders(tmp_path):
    out = tmp_path / "covering.png"
    fig = plots.cta_covering_set(save=str(out), n_genes=8)
    assert out.exists() and out.stat().st_size > 0
    ax = fig.axes[0]
    assert ax.get_ylabel() == "cumulative coverage (%)"
    assert len(ax.get_xticklabels()) <= 8
    # Both curves are cumulative, so they never decrease.
    for line in ax.get_lines():
        ydata = list(line.get_ydata())
        assert ydata == sorted(ydata)


def test_cta_covering_set_rejects_an_unknown_stat():
    with pytest.raises(ValueError, match="stat must be one of"):
        plots.cta_covering_set(stat="p99")


def test_cli_stat_is_opt_in_so_each_plot_keeps_its_default(monkeypatch, tmp_path):
    # --stat unset must not force "median" onto cta-covering-set, whose own
    # default is q3 (the subset-antigen-appropriate statistic).
    seen = {}
    monkeypatch.setattr(plots, "cta_covering_set", lambda **kw: seen.update(kw))

    assert cli.main(["plot", "cta-covering-set", "--out", str(tmp_path / "a.png")]) == 0
    assert "stat" not in seen

    seen.clear()
    assert (
        cli.main(["plot", "cta-covering-set", "--stat", "q1", "--out", str(tmp_path / "b.png")])
        == 0
    )
    assert seen["stat"] == "q1"


def test_cli_expression_heatmap_still_defaults_to_median(monkeypatch, tmp_path):
    seen = {}
    monkeypatch.setattr(plots, "cta_expression_heatmap", lambda **kw: seen.update(kw))
    assert cli.main(["plot", "cta-expression-heatmap", "--out", str(tmp_path / "h.png")]) == 0
    assert seen["stat"] == "median"


def test_heatmap_missing_cells_are_not_the_figure_ground(tmp_path):
    # On a white ground an unpainted NaN cell reads as magma's pale-yellow top end,
    # inverting "absent" into "highest". It must get a colour outside the ramp.
    grid = pd.DataFrame(
        {"GENE_A": [10.0, float("nan")], "GENE_B": [1.0, 2.0]}, index=["LUAD", "SKCM"]
    )
    fig = plots._cohort_gene_heatmap(
        grid,
        title="t",
        cbar_label="c",
        cmap="magma",
        lognorm=True,
        save=str(tmp_path / "h.png"),
    )
    bad = np.asarray(fig.axes[0].images[0].cmap.get_bad(), dtype=float)
    assert bad[3] > 0, "missing cells must be painted, not left transparent"
    assert not np.allclose(bad[:3], 1.0), "missing cells must not be the white ground"


def test_stack_size_is_bounded_however_many_rows():
    from oncoref import figure_style

    # Unclamped, the per-row spacing applies directly.
    inches, density = figure_style.stack_size(10, per_item=0.42, floor=6)
    assert inches == pytest.approx(4.2 if 4.2 > 6 else 6)
    assert density == 1.0

    # Past the cap, the figure stops growing and reports how far it compressed.
    inches, density = figure_style.stack_size(1000, per_item=0.42, floor=6)
    assert inches == figure_style.MAX_FIGURE_INCHES
    assert 0 < density < 1
    # No row count, however absurd, can exceed the cap.
    for n in (85, 500, 10_000, 1_000_000):
        got, _ = figure_style.stack_size(n, per_item=0.42, floor=6)
        assert got <= figure_style.MAX_FIGURE_INCHES


def test_tick_fontsize_shrinks_with_density_but_has_a_floor():
    from oncoref import figure_style

    assert figure_style.tick_fontsize(1.0, 8.0) == 8.0
    assert figure_style.tick_fontsize(0.5, 8.0) == 4.0
    # Never vanishes, however dense.
    assert figure_style.tick_fontsize(0.0001, 8.0) == figure_style.MIN_TICK_FONTSIZE


def test_ici_regimen_comparison_stays_within_the_page_cap(tmp_path):
    from oncoref import figure_style

    # This is the figure that reached 36 inches: every cancer type on its own row.
    out = tmp_path / "ici.png"
    fig = plots.ici_regimen_comparison(save=str(out), min_regimens=1)
    width, height = fig.get_size_inches()
    assert height <= figure_style.MAX_FIGURE_INCHES
    assert width <= figure_style.MAX_FIGURE_INCHES
    assert out.exists() and out.stat().st_size > 0


def test_every_cta_expression_heatmap_axis_is_capped(tmp_path):
    from oncoref import figure_style

    fig = plots.cta_expression_heatmap(stat="q3", save=str(tmp_path / "h.png"))
    width, height = fig.get_size_inches()
    assert width <= figure_style.MAX_FIGURE_INCHES
    assert height <= figure_style.MAX_FIGURE_INCHES


def test_grouped_barh_bars_are_row_sized_not_figure_sized():
    # The figure-height and bar-height variables must stay distinct: shadowing one
    # with the other drew every bar 18.5 DATA units tall, so the bars detached from
    # their row labels entirely and merged into solid blocks.
    fig = plots.burden_category_bars(region="us")
    ax = fig.axes[0]
    heights = {round(p.get_height(), 3) for p in ax.patches}
    assert heights == {0.4}, heights
    # Two series share a row, so a bar must be at most half a row.
    assert max(heights) <= 0.5
    # And the axis spans the categories, not some multiple of them.
    lo, hi = sorted(ax.get_ylim())
    assert hi - lo < len(ax.get_yticks()) + 4


def test_bar_rows_stay_legible_relative_to_the_figure():
    # "Tiny labels, fat bars" is an aspect-ratio failure: a 37-row chart at 0.5 in/row
    # is an 18in ribbon. Keep rows dense enough that label and bar weight are comparable.
    fig = plots.burden_category_bars(region="us")
    width, height = fig.get_size_inches()
    n_rows = len(fig.axes[0].get_yticks())
    assert height / n_rows < 0.4, "rows too tall; labels will be dwarfed by bars"
    assert height / width < 2.0, "figure too tall and narrow to read as one chart"


@pytest.fixture
def _print_preset():
    from oncoref import figure_style

    yield
    figure_style.use("print")


def test_slide_preset_trims_rows_and_says_so(_print_preset):
    from oncoref import figure_style

    figure_style.use("slide")
    fig = plots.burden_category_bars(region="us")
    ax = fig.axes[0]
    assert len(ax.get_yticks()) == figure_style.max_items()
    # Trimming is disclosed on the axis, never silent.
    assert "top 12 of 37" in ax.get_xlabel()
    # And it keeps the head of the ranking, not an arbitrary slice.
    assert "prostate" in [t.get_text() for t in ax.get_yticklabels()]


def test_slide_preset_fits_a_slide(_print_preset):
    from oncoref import figure_style

    figure_style.use("slide")
    width, height = plots.burden_category_bars(region="us").get_size_inches()
    assert (width, height) <= figure_style.SLIDE_ASPECT
    assert width / height > 1.0, "a slide is landscape"


def test_slide_preset_keeps_every_scatter_point_and_rations_only_labels(_print_preset):
    from oncoref import figure_style

    printed = plots.apd1_vs_tmb()
    n_points = sum(len(c.get_offsets()) for c in printed.axes[0].collections)

    figure_style.use("slide")
    slide = plots.apd1_vs_tmb()
    ax = slide.axes[0]
    # Every point survives — trimming a scatter would misstate the distribution.
    assert sum(len(c.get_offsets()) for c in ax.collections) == n_points
    # But only the top slice is labelled, so the text is readable.
    labels = [t for t in ax.texts if t.get_text()]
    assert 0 < len(labels) <= figure_style.max_items()


def test_slide_labels_drop_underscores_and_print_keeps_them(_print_preset):
    from oncoref import figure_style

    figure_style.use("print")
    assert figure_style.label("non_hodgkin_lymphoma") == "non_hodgkin_lymphoma"
    figure_style.use("slide")
    assert figure_style.label("non_hodgkin_lymphoma") == "non hodgkin lymphoma"


def test_unknown_preset_is_rejected(_print_preset):
    from oncoref import figure_style

    with pytest.raises(ValueError, match="unknown preset"):
        figure_style.use("poster")


def test_font_scale_preserves_the_size_hierarchy(_print_preset):
    from oncoref import figure_style

    # Scaling must multiply, not lift to a floor: a 4pt label in a dense panel and a
    # 9pt row label encode how crowded each panel is, and flattening them collides.
    figure_style.use("print")
    small, large = figure_style.fs(4), figure_style.fs(9)
    figure_style.use("slide")
    assert figure_style.fs(4) / small == pytest.approx(figure_style.fs(9) / large, rel=1e-3)
    assert figure_style.fs(4) < figure_style.fs(9)


def test_row_height_grows_sublinearly_with_the_text_scale(_print_preset):
    from oncoref import figure_style

    figure_style.use("print")
    base = figure_style.row_height(0.30)
    figure_style.use("slide")
    slide = figure_style.row_height(0.30)
    scale = figure_style.font_scale()
    # Taller labels need more room...
    assert slide > base
    # ...but strictly proportional growth would make a slide chart taller than the slide.
    assert slide < base * scale


def test_log_axis_keeps_non_positive_points_via_symlog():
    # A plain log axis deletes every non-positive point without warning. A cohort at
    # 0 mut/Mb belongs in the low-burden quadrant; dropping it inverts the message.
    points = [("A", 0.0, 10.0), ("B", 1.0, 20.0), ("C", 10.0, 30.0)]
    fig = plots._family_scatter(
        points, xlabel="x", ylabel="y", title="t", logx=True, annotate=False
    )
    ax = fig.axes[0]
    assert ax.get_xscale() == "symlog"
    assert sum(len(c.get_offsets()) for c in ax.collections) == 3


def test_log_axis_stays_plain_log_when_every_value_is_positive():
    points = [("A", 1.0, 10.0), ("B", 10.0, 20.0)]
    fig = plots._family_scatter(
        points, xlabel="x", ylabel="y", title="t", logx=True, annotate=False
    )
    assert fig.axes[0].get_xscale() == "log"


@pytest.fixture
def _no_highlight():
    from oncoref import figure_style

    yield
    figure_style.set_highlight(None)


def test_highlight_mutes_context_without_changing_positions(_no_highlight):
    from oncoref import figure_style

    plain = plots.apd1_vs_tmb(annotate=False)
    plain_xy = np.vstack([c.get_offsets() for c in plain.axes[0].collections])

    figure_style.set_highlight("BRCA_Basal")
    hot = plots.apd1_vs_tmb(annotate=False)
    hot_xy = np.vstack([c.get_offsets() for c in hot.axes[0].collections])

    # A highlight is a display decision: same points, same places, nothing filtered.
    assert plain_xy.shape == hot_xy.shape
    assert np.allclose(np.sort(plain_xy, axis=0), np.sort(hot_xy, axis=0))


def test_highlighted_mark_is_emphasised(_no_highlight):
    from oncoref import figure_style

    figure_style.set_highlight("BRCA_Basal")
    ax = plots.apd1_vs_tmb(annotate=False).axes[0]
    sizes = sorted({float(s) for c in ax.collections for s in np.atleast_1d(c.get_sizes())})
    # Exactly one mark is enlarged, and it is bigger than the rest.
    assert len(sizes) == 2 and sizes[1] > sizes[0]


def test_highlighted_point_is_labelled_even_past_the_slide_budget(_no_highlight):
    from oncoref import figure_style

    figure_style.use("slide")
    figure_style.set_highlight("BRCA_Basal")
    try:
        ax = plots.apd1_vs_tmb().axes[0]
        labels = [t.get_text() for t in ax.texts if t.get_text()]
        # BRCA_Basal has a low ORR so it is nowhere near the top-12 by y.
        assert "BRCA_Basal" in labels
    finally:
        figure_style.use("print")


def test_mute_blends_toward_the_page_not_toward_black(_no_highlight):
    from oncoref import figure_style

    muted = figure_style.mute("#2a7f4f")
    # Fading must lighten, so overlapping context marks do not compound into blobs.
    assert all(ch > 0.5 for ch in muted)


def test_no_highlight_leaves_colours_untouched(_no_highlight):
    from oncoref import figure_style

    figure_style.set_highlight(None)
    plain, _ = plots._family_colors(["LUAD", "BRCA_Basal"])
    figure_style.set_highlight("BRCA_Basal")
    hot, _ = plots._family_colors(["LUAD", "BRCA_Basal"])
    assert plain["LUAD"] != hot["LUAD"]
    assert hot["BRCA_Basal"] == figure_style.HIGHLIGHT_COLOR

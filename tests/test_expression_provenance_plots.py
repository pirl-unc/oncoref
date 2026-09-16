# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

import matplotlib.pyplot as plt
import pytest

from oncoref import cli, figure_style
from oncoref import expression_provenance_plots as epp
from oncoref.load_dataset import get_data

pytestmark = pytest.mark.xdist_group("plots")


@pytest.fixture(autouse=True)
def _close_figures_after_test():
    yield
    plt.close("all")


def test_availability_defaults_to_selected_sources_only():
    every = epp.availability(selected_only=False)
    selected = epp.availability()
    assert len(selected) < len(every)
    # One selected source per cancer code, so the figures never double-count samples.
    assert selected["cancer_code"].is_unique


def test_coverage_summary_is_consistent_with_the_registry():
    summary = epp.coverage_summary()
    registry = get_data("cancer-type-registry")
    assert summary["n_codes_total"] == len(registry)
    assert 0 < summary["n_codes_covered"] <= summary["n_codes_total"]
    assert summary["n_samples"] == int(epp.availability()["n_reference_samples"].sum())
    assert summary["n_tpm_exact_codes"] <= summary["n_codes_covered"]


def test_burden_coverage_never_exceeds_the_total_burden():
    for region in ("us", "world"):
        rows, totals = epp.burden_coverage(region)
        assert not rows.empty
        assert 0 < totals["incidence_covered"] <= totals["incidence_total"]
        assert 0 < totals["mortality_covered"] <= totals["mortality_total"]
        # A category counts as covered iff at least one of its codes has a source.
        assert (rows["covered"] == (rows["sum"] > 0)).all()
        assert (rows["sum"] <= rows["size"]).all()


def test_renderer_writes_every_figure(tmp_path):
    result = epp.render(out_dir=tmp_path)
    assert set(result["paths"]) == set(epp.FILENAMES)
    for path in result["paths"].values():
        assert path.exists() and path.stat().st_size > 0
    assert result["summary"]["n_samples"] > 0


def test_renderer_honors_region(tmp_path):
    us = epp.render(out_dir=tmp_path / "us", region="us")["burden"]
    world = epp.render(out_dir=tmp_path / "world", region="world")["burden"]
    assert us != world


def test_cli_plot_expression_provenance_delegates(monkeypatch, tmp_path, capsys):
    captured = {}

    def fake(*, out_dir, region):
        captured["out_dir"] = out_dir
        captured["region"] = region
        out = tmp_path / "prov"
        out.mkdir(exist_ok=True)
        fig = out / "expr-source-samples.png"
        fig.write_text("ok\n")
        return {
            "summary": {
                "n_samples": 10,
                "n_codes_covered": 2,
                "n_codes_total": 4,
                "n_sources": 1,
            },
            "burden": {
                "incidence_covered": 5.0,
                "incidence_total": 10.0,
                "mortality_covered": 4.0,
                "mortality_total": 10.0,
            },
            "paths": {"source_samples": fig},
        }

    monkeypatch.setattr(epp, "render", fake)
    rc = cli.main(["plot", "expression-provenance", "--out", str(tmp_path / "prov")])
    assert rc == 0
    assert captured["out_dir"] == str(tmp_path / "prov")
    assert captured["region"] == "us"
    assert "source_samples" in capsys.readouterr().out


def test_figures_are_saved_on_an_opaque_white_ground(tmp_path):
    import matplotlib.image as mpimg

    path = epp.render(out_dir=tmp_path)["paths"]["source_samples"]
    img = mpimg.imread(path)
    # Top-left corner is figure margin: opaque and white, never transparent or grey.
    corner = img[0, 0]
    assert corner[:3].min() > 0.98
    if len(corner) == 4:
        assert corner[3] == 1.0


def test_minimal_style_strips_titles_and_restores_them():
    figure_style.apply()
    fig, ax = plt.subplots()
    ax.set_title("burned-in title")
    fig.suptitle("burned-in suptitle")

    figure_style.finalize(fig)
    assert ax.get_title() == ""
    assert fig._suptitle.get_text() == ""

    # Opting out keeps the chatty labelling for interactive/notebook use.
    ax.set_title("kept")
    figure_style.finalize(fig, keep_titles=True)
    assert ax.get_title() == "kept"


def test_finalize_removes_top_and_right_spines():
    fig, ax = plt.subplots()
    figure_style.finalize(fig)
    assert not ax.spines["top"].get_visible()
    assert not ax.spines["right"].get_visible()
    assert ax.spines["left"].get_visible()
    assert ax.spines["bottom"].get_visible()

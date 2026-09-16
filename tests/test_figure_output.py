# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

import re

import pytest

from oncoref import figure_output

_RUN = re.compile(r"^run_\d{8}-\d{6}(-\d+)?$")


def test_new_run_dir_is_timestamped_and_linked(tmp_path):
    run = figure_output.new_run_dir(tmp_path)
    assert run.is_dir()
    assert _RUN.match(run.name)
    link = tmp_path / figure_output.LATEST
    assert link.is_symlink() and link.resolve() == run.resolve()


def test_two_runs_in_the_same_second_do_not_collide(tmp_path):
    first = figure_output.new_run_dir(tmp_path, stamp="20260915-120000")
    second = figure_output.new_run_dir(tmp_path, stamp="20260915-120000")
    assert first != second
    assert first.is_dir() and second.is_dir()
    assert (tmp_path / figure_output.LATEST).resolve() == second.resolve()


def test_latest_is_repointed_not_duplicated(tmp_path):
    figure_output.new_run_dir(tmp_path, stamp="20260915-120000")
    second = figure_output.new_run_dir(tmp_path, stamp="20260915-130000")
    link = tmp_path / figure_output.LATEST
    assert link.resolve() == second.resolve()
    runs = sorted(p.name for p in tmp_path.iterdir() if p.is_dir() and not p.is_symlink())
    assert runs == ["run_20260915-120000", "run_20260915-130000"]


def test_explicit_out_is_used_verbatim(tmp_path):
    target = tmp_path / "nested" / "mine.png"
    resolved = figure_output.resolve(str(target))
    assert resolved == target
    assert target.parent.is_dir()
    # An explicit path never gets a run directory grafted onto it.
    assert not (tmp_path / "figures").exists()


def test_resolve_builds_a_named_png_in_a_fresh_run(tmp_path):
    path = figure_output.resolve(None, name="cta-covering-set", base=tmp_path)
    assert path.name == "cta-covering-set.png"
    assert _RUN.match(path.parent.name)


def test_resolve_builds_a_directory_for_multi_figure_families(tmp_path):
    path = figure_output.resolve(None, name="cta-curation", suffix="", base=tmp_path)
    assert path.is_dir() and path.name == "cta-curation"
    assert _RUN.match(path.parent.name)


def test_env_var_overrides_the_figures_root(tmp_path, monkeypatch):
    monkeypatch.setenv(figure_output.FIGURES_DIR_ENV_VAR, str(tmp_path / "elsewhere"))
    run = figure_output.new_run_dir()
    assert run.parent == tmp_path / "elsewhere"


def test_figures_root_defaults_to_figures(monkeypatch):
    monkeypatch.delenv(figure_output.FIGURES_DIR_ENV_VAR, raising=False)
    assert figure_output.figures_root().name == "figures"


@pytest.mark.parametrize(
    "which,expect_dir",
    [("cta-curation", True), ("expression-provenance", True), ("apd1-vs-tmb", False)],
)
def test_cli_plot_without_out_writes_into_a_run_dir(
    which, expect_dir, tmp_path, monkeypatch, capsys
):
    from oncoref import cli

    monkeypatch.setenv(figure_output.FIGURES_DIR_ENV_VAR, str(tmp_path / "figures"))
    assert cli.main(["plot", which]) == 0
    capsys.readouterr()
    runs = [p for p in (tmp_path / "figures").iterdir() if p.is_dir() and not p.is_symlink()]
    assert len(runs) == 1 and _RUN.match(runs[0].name)
    produced = list(runs[0].iterdir())
    assert produced
    assert produced[0].is_dir() is expect_dir

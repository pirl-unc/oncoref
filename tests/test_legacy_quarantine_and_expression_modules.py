"""Tests for legacy CTA module quarantine (`oncoref._legacy`) and `expression.py` modularization (#756)."""

from __future__ import annotations

import importlib
from pathlib import Path

from oncoref import (
    _legacy,
    expression,
    expression_availability,
    expression_pan_cancer,
    expression_qc,
    expression_universe_deltas,
)


def test_quarantined_legacy_cta_modules_resolve_identically_to_top_level_shims():
    assert len(_legacy.QUARANTINED_CTA_MODULES) == 22
    pkg_dir = Path(expression.__file__).resolve().parent

    for mod_name in _legacy.QUARANTINED_CTA_MODULES:
        shim_path = pkg_dir / f"{mod_name}.py"
        impl_path = pkg_dir / "_legacy" / f"{mod_name}.py"
        assert shim_path.is_file()
        assert impl_path.is_file()
        # Top-level shim is compact (<25 lines)
        assert len(shim_path.read_text(encoding="utf-8").splitlines()) < 25

        top_mod = importlib.import_module(f"oncoref.{mod_name}")
        leg_mod = importlib.import_module(f"oncoref._legacy.{mod_name}")
        # Exact module-object identity preserves monkeypatch.setattr(oncoref.<mod>, ...) semantics
        assert top_mod is leg_mod


def test_active_shared_cta_modules_remain_first_class_top_level_implementations():
    pkg_dir = Path(expression.__file__).resolve().parent
    for active in (
        "cta_literature",
        "cta_sequence_overlap",
        "cta_selection_policy",
        "cta_evidence_model",
        "cta_evidence_schema",
        "cta_evidence_adapters",
    ):
        assert not (pkg_dir / "_legacy" / f"{active}.py").exists()
        mod = importlib.import_module(f"oncoref.{active}")
        assert Path(mod.__file__).resolve() == pkg_dir / f"{active}.py"


def test_expression_submodules_extract_cohesive_surfaces_and_preserve_reexports():
    expr_lines = len(Path(expression.__file__).read_text(encoding="utf-8").splitlines())
    assert expr_lines < 4500

    assert expression.sample_expression_qc is expression_qc.sample_expression_qc
    assert (
        expression.housekeeping_cancer_expression_coverage
        is expression_qc.housekeeping_cancer_expression_coverage
    )
    assert (
        expression.expression_artifact_build_metadata
        is expression_qc.expression_artifact_build_metadata
    )
    assert (
        expression.cancer_reference_expression_availability
        is expression_availability.cancer_reference_expression_availability
    )
    assert (
        expression.umbrella_expression_coverage
        is expression_availability.umbrella_expression_coverage
    )
    assert (
        expression.expression_artifact_gene_universe_deltas
        is expression_universe_deltas.expression_artifact_gene_universe_deltas
    )
    assert expression.pan_cancer_expression is expression_pan_cancer.pan_cancer_expression

"""Relaxed cell sampling must not silently relax independent embryo support."""

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/audit_cta_recovery.py"
SPEC = importlib.util.spec_from_file_location("recovery_audit", SCRIPT)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


@pytest.mark.parametrize(
    "donors,stages,expected",
    [
        (["one", "one"], ["3 day", "3 day"], 0),
        (["one", ""], ["3 day", "3 day"], 0),
        (["one", " "], ["3 day", "3 day"], 0),
        (["one", "two"], ["3 day", "4 day"], 0),
        (["one", "two"], ["3 day", "3 day"], 1),
    ],
)
def test_three_cell_counterfactual_keeps_stage_specific_embryo_replication(
    tmp_path, monkeypatch, donors, stages, expected
):
    monkeypatch.setattr(AUDIT, "OUT", tmp_path)
    catalog = pd.DataFrame(
        [
            {
                "entity_id": "A",
                "R_supported": 0.0,
                "H0": 1.0,
                "shared_pass": True,
                "tumor_fold": 20.0,
                "classification": "Cancer Antigen only",
            }
        ]
    ).set_index("entity_id")
    rows = pd.DataFrame(
        [
            {
                "entity_id": "A",
                "scope": "developmental",
                "source_id": "petropoulos_E-MTAB-3929",
                "cells": 3,
                "positive_cells": 2,
                "value": 30.0,
                "identity_resolved": True,
                "participant_id": donor,
                "developmental_stage": stage,
                "source_context": "lineage",
            }
            for donor, stage in zip(donors, stages)
        ]
    )
    before = catalog.copy(deep=True)
    assert len(AUDIT.embryo_sensitivity(catalog, rows)) == expected
    pd.testing.assert_frame_equal(catalog, before)

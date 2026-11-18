"""An IHC counterfactual must not turn missing coverage into negative evidence."""

import importlib.util
from pathlib import Path

import pandas as pd

PATH = Path(__file__).resolve().parents[1] / "scripts/audit_reproductive_ihc_route.py"
SPEC = importlib.util.spec_from_file_location("ihc_counterfactual", PATH)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def observations():
    return pd.DataFrame(
        [
            ["a", "testis", "spermatogonia", "High", "Enhanced", True],
            ["a", "brain", "neurons", "Not detected", "Enhanced", False],
        ],
        columns=[
            "canonical_gene_id",
            "Tissue",
            "Cell type",
            "Level",
            "Reliability",
            "reproductive",
        ],
    )


def test_observed_negatives_do_not_claim_complete_coverage():
    result = AUDIT.assess_group(observations(), ["a"], {"brain", "heart"})
    assert result["ihc_route"]
    assert not result["full_source_tissue_panel_complete"]
    assert not AUDIT.assess_group(observations(), ["a", "b"], {"brain"})["ihc_route"]


def test_positive_unknown_or_unreliable_normal_row_blocks_negative_claim():
    for level in ["Low", "Medium", "High", "N/A", ""]:
        rows = observations()
        rows.loc[1, "Level"] = level
        assert not AUDIT.assess_group(rows, ["a"], {"brain"})["ihc_route"]
    for i in [0, 1]:
        rows = observations()
        rows.loc[i, "Reliability"] = "Uncertain"
        assert not AUDIT.assess_group(rows, ["a"], {"brain"})["ihc_route"]
    assert not AUDIT.assess_group(observations().iloc[:1], ["a"], {"brain"})["ihc_route"]

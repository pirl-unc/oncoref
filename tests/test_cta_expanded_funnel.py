"""Boundary and evidence-integrity checks for the expanded research funnel."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def script(name):
    path = Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_tumor_maximum_is_an_actual_sample_and_keeps_missing_distinct_from_zero():
    module = script("cta_tumor_union")
    frame = pd.DataFrame(
        {
            "Ensembl_Gene_ID": ["boundary", "low", "missing", "zero"],
            "a": [2.0, 1.9999, np.nan, 0.0],
            "b": [1.0, np.nan, np.nan, 0.0],
        }
    )
    result = module.summarize_matrix(frame, ["a", "b"]).set_index("gene_id")
    assert result.loc["boundary", "samples_ge2"] == 1
    assert result.loc["boundary", "maximum_sample"] == "a"
    assert result.loc["low", "samples_ge2"] == 0
    assert pd.isna(result.loc["missing", "maximum_tpm"])
    assert result.loc["missing", "maximum_sample"] == ""
    assert result.loc["zero", "maximum_tpm"] == 0
    frame.loc[0, "a"] = -1
    with pytest.raises(ValueError, match="nonnegative"):
        module.summarize_matrix(frame, ["a", "b"])


def test_optional_panels_do_not_waive_required_bulk_or_immune_coverage():
    module = script("cta_expanded_funnel")
    row = {
        "normal_bulk_tissues": 50,
        "cell_records": 426,
        "hpa_immune_cell_records": 19,
        "hpa_immune_monaco_records": 30,
        "reproductive_ntpm": 5,
        "adult_bulk_max_ntpm": 1,
        "hpa_immune_cell_max_nTPM": 0,
        "hpa_immune_monaco_max_TPM": 0,
        "cell_supported_records": 0,
        "cell_max_ntpm": np.nan,
        "dice_sources": np.nan,
        "cornea_records": np.nan,
        "dice_third_max_tpm": np.nan,
    }
    facts = pd.DataFrame(
        [row, row | {"normal_bulk_tissues": 41}, row | {"hpa_immune_cell_records": 18}]
    )
    assert module.coverage(facts).tolist() == [True, False, False]
    assert not module.coverage(facts, optional=False).any()


def test_streamed_winners_match_full_merge_and_preserve_source_identity(tmp_path):
    module = script("cta_tumor_union")
    frames = []
    paths = []
    for code, values in [("Z", [2, np.nan, 4]), ("A", [2, 0, 3]), ("M", [1, 1, 5])]:
        frame = pd.DataFrame(
            {
                "gene_id": ["tie", "missing_in_first", "different_winner"],
                "maximum_tpm": values,
                "cancer_code": code,
                "maximum_sample": f"{code}-sample",
                "source_identity": f"{code}-identity",
            }
        )
        path = tmp_path / f"{code}.csv"
        frame.to_csv(path, index=False)
        paths.append(path)
        frames.append(pd.read_csv(path))
    expected = (
        pd.concat(frames, ignore_index=True)
        .sort_values(
            ["maximum_tpm", "cancer_code", "maximum_sample"], ascending=[False, True, True]
        )
        .drop_duplicates("gene_id")
        .reset_index(drop=True)
    )
    result = module.winning_observations(paths, entity_column="gene_id")
    pd.testing.assert_frame_equal(result, expected)
    assert result.set_index("gene_id").loc["tie", "source_identity"] == "A-identity"

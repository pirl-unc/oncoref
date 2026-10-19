"""Reject incompatible source scales before they become plausible TPM output."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from oncoref import expression, expression_builders, normalization
from oncoref.expression_engine import (
    canonicalize_source_gene_matrix,
    validate_nonnegative_expression,
)


@pytest.mark.parametrize("value", [-0.01, np.inf, -np.inf, "broken"])
@pytest.mark.parametrize("unit", ["TPM", "FPKM", "RPKM", "raw_counts", "log2(TPM+1)"])
def test_declared_scale_rejects_invalid_values_without_clipping(value, unit):
    frame = pd.DataFrame({"gene_id": ["g1", "g2"], "sample": [value, 100.0]})
    before = frame.copy(deep=True)
    with pytest.raises(ValueError, match=r"negative or infinite|must be numeric"):
        expression_builders.normalize_source_matrix_to_tpm(
            frame, unit=unit, row_id_col="gene_id", gene_lengths_kb={"g1": 1.0, "g2": 1.0}
        )
    pd.testing.assert_frame_equal(frame, before)


@pytest.mark.parametrize("unit", ["TPM", "FPKM", "RPKM", "raw_counts", "log2(TPM+1)"])
def test_unit_conversion_preserves_unmeasured_values(unit):
    frame = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "partial": [2.0, np.nan],
            "missing": [np.nan, np.nan],
            "zero": [0.0, 0.0],
        }
    )
    out = expression_builders.normalize_source_matrix_to_tpm(
        frame, unit=unit, row_id_col="gene_id", gene_lengths_kb={"g1": 1.0, "g2": 1.0}
    )
    assert pd.isna(out.loc[1, "partial"])
    assert out["missing"].isna().all()
    assert out["zero"].eq(0).all()
    assert out.loc[0, "partial"] == (3.0 if unit == "log2(TPM+1)" else 1_000_000.0)


def test_log_input_is_inverted_before_linear_operations():
    frame = pd.DataFrame({"gene_id": ["g1", "g2"], "sample": [2.0, 4.0]})
    out = expression_builders.normalize_source_matrix_to_tpm(
        frame, unit="log2(TPM+1)", row_id_col="gene_id"
    )
    assert out["sample"].tolist() == [3.0, 15.0]
    assert out["sample"].sum() == 18.0


def test_inverse_transform_rejects_overflow():
    frame = pd.DataFrame({"gene_id": ["g1"], "sample": [2000.0]})
    with pytest.raises(ValueError, match=r"inverse log2.*negative or infinite"):
        expression_builders.normalize_source_matrix_to_tpm(
            frame, unit="log2(TPM+1)", row_id_col="gene_id"
        )


@pytest.mark.parametrize(
    "lengths", [{"g1": 1.0}, {"g1": 1.0, "g2": -2.0}, {"g1": 1.0, "g2": np.inf}]
)
def test_count_conversion_does_not_turn_unknown_gene_lengths_into_zeros(lengths):
    frame = pd.DataFrame({"gene_id": ["g1", "g2"], "sample": [5.0, 10.0]})
    with pytest.raises(ValueError, match="finite positive length for every source row"):
        expression_builders.normalize_source_matrix_to_tpm(
            frame, unit="raw_counts", row_id_col="gene_id", gene_lengths_kb=lengths
        )


@pytest.mark.parametrize("mode", ["tpm_raw", "tpm_clean", "tpm_clean_log1p", "tpm_clean_hk"])
def test_source_reader_validates_before_canonicalization(tmp_path, monkeypatch, mode):
    # Summing duplicate IDs would hide the invalid negative measurement.
    frame = pd.DataFrame(
        {
            "Ensembl_Gene_ID": ["G1", "G1"],
            "Symbol": ["A", "A"],
            "s1": [-2.0, 5.0],
        }
    )
    path = tmp_path / "invalid.parquet"
    frame.to_parquet(path, index=False)

    def must_not_canonicalize(*args, **kwargs):
        pytest.fail("invalid source reached linear row aggregation")

    monkeypatch.setattr(expression, "_canonicalize_gene_rows", must_not_canonicalize)
    with pytest.raises(ValueError, match="negative or infinite"):
        expression._load_per_sample_matrix(str(path), path.stat().st_mtime, mode)


def test_source_validator_accepts_missing_and_zero_without_mutation():
    values = pd.DataFrame({"s1": [np.nan, 0.0, 3.0]})
    before = values.copy()
    validate_nonnegative_expression(values, context="test")
    pd.testing.assert_frame_equal(values, before)


def test_clean_tpm_does_not_normalize_a_transformed_negative_field():
    with pytest.raises(ValueError, match="negative or infinite"):
        normalization.clean_tpm(pd.DataFrame({"s1": [-1.0, 10.0]}))


def test_builder_canonicalization_cannot_hide_negative_duplicate_rows():
    frame = pd.DataFrame(
        {
            "Ensembl_Gene_ID": ["ENSG00000184507", "ENSG00000184507"],
            "Symbol": ["NUTM1", "NUTM1"],
            "sample": [-2.0, 10.0],
        }
    )
    with pytest.raises(ValueError, match=r"source canonicalization.*negative or infinite"):
        canonicalize_source_gene_matrix(frame)


@pytest.mark.parametrize("value", [-1.0, np.inf, "broken"])
def test_gdc_reader_does_not_hide_invalid_tpm(tmp_path, value):
    path = tmp_path / "counts.tsv"
    pd.DataFrame(
        {
            "gene_id": ["ENSG00000184507"],
            "gene_name": ["NUTM1"],
            "gene_type": ["protein_coding"],
            "tpm_unstranded": [value],
        }
    ).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match=r"negative or infinite|must be numeric"):
        expression_builders.read_gdc_star_counts_tpm(path)


def test_gdc_assembly_keeps_unmeasured_genes_distinct_from_zero(tmp_path):
    paths = {}
    for sample, genes, values in [
        (
            "s1",
            [
                "ENSG5.1",
                "ENSG1",
                "ENSG2.1",
                "ENSG2.2",
                "ENSG3.1",
                "ENSG3.2",
                "ENSG4.1",
                "ENSG4.2",
                "ENSG5.2",
            ],
            [2.0, 10.0, np.nan, np.nan, 0.0, np.nan, np.nan, 5.0, 3.0],
        ),
        ("s2", ["ENSG1"], [20.0]),
    ]:
        path = tmp_path / f"{sample}.tsv"
        pd.DataFrame(
            {
                "gene_id": genes,
                "gene_name": [gene.split(".")[0] for gene in genes],
                "gene_type": "protein_coding",
                "tpm_unstranded": values,
            }
        ).to_csv(path, sep="\t", index=False)
        paths[sample] = path
    out = expression_builders._gdc_raw_tpm_matrix(
        pd.DataFrame({"sample_id": ["s1", "s2"], "source_file_id": ["s1", "s2"]}),
        cache_dir=tmp_path,
        file_paths=paths,
    ).set_index("source_gene_id")
    assert pd.isna(out.loc["ENSG2", "s1"])
    assert pd.isna(out.loc["ENSG2", "s2"])
    assert out.loc["ENSG3", "s1"] == 0.0
    assert pd.isna(out.loc["ENSG3", "s2"])
    assert out.loc["ENSG4", "s1"] == 5.0
    assert out.loc["ENSG5", "s1"] == 5.0
    assert out.index.is_unique
    assert out["Symbol"].to_dict() == {f"ENSG{i}": f"ENSG{i}" for i in range(1, 6)}


def test_staging_validates_all_inputs_before_copying_any(tmp_path, monkeypatch):
    path = Path(__file__).parents[1] / "scripts" / "stage_source_matrices.py"
    spec = importlib.util.spec_from_file_location("stage_scale_test", path)
    stage = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = stage
    spec.loader.exec_module(stage)
    cache = tmp_path / "input"
    derived = cache / "TEST" / "derived"
    derived.mkdir(parents=True)
    for code, value in [("A", 10.0), ("B", -1.0)]:
        pd.DataFrame({"Ensembl_Gene_ID": ["G1"], "Symbol": ["A"], "s1": [value]}).to_parquet(
            derived / f"{code}_per_sample_tpm.parquet", index=False
        )
    monkeypatch.setattr(
        stage.sm,
        "registry",
        lambda: pd.DataFrame(
            {
                "cancer_code": ["A", "B"],
                "source_cohort": ["TEST", "TEST"],
                "n_samples": [1, 1],
            }
        ),
    )
    published = tmp_path / "published"
    monkeypatch.setattr(stage.sm, "local_path", lambda code: published / f"{code}.parquet")
    release = tmp_path / "release"
    with pytest.raises(ValueError, match=r"B: linear TPM.*negative or infinite"):
        stage.stage(cache, release_dir=release, codes=None, limit=None)
    assert not published.exists()
    assert not list(release.glob("*.parquet"))


def test_nutm_source_is_not_advertised_as_comparable_linear_tpm():
    metadata = expression.cancer_reference_expression_source_metadata("NUTM")
    assert metadata["source_scale_class"] == "invalid_legacy_tpm"
    assert metadata["linear_tpm_comparable"] is False
    assert metadata["unit"] == "unverified transformed expression"
    assert "2163 negative" in metadata["notes"]

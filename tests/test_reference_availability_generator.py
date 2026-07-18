import pandas as pd
import pytest
from scripts.generate_reference_availability import build_reference_availability


def _write_shard(path, code, source, genes, n_samples):
    pd.DataFrame(
        {
            "Ensembl_Gene_ID": genes,
            "cancer_code": code,
            "source_cohort": source,
            "source_project": "PROJECT",
            "source_version": "v1",
            "tumor_origin": "primary",
            "metastasis_site": pd.NA,
            "processing_pipeline": "rnaseq",
            "notes": "line one\nline two",
            "n_samples": n_samples,
        }
    ).to_csv(path, index=False)


def test_build_reference_availability_is_source_specific_and_selects_richest(tmp_path):
    _write_shard(tmp_path / "small.csv", "X", "SMALL", ["E1", "E2"], 20)
    _write_shard(tmp_path / "rich.csv", "X", "RICH", ["E1", "E2", "E3"], 5)
    _write_shard(tmp_path / "other.csv", "Y", "OTHER", ["E1"], 4)

    table = build_reference_availability(tmp_path, chunksize=1)

    x = table[table["cancer_code"] == "X"].set_index("source_cohort")
    assert x.loc["RICH", "n_reference_genes"] == 3
    assert x.loc["RICH", "selected"]
    assert not x.loc["SMALL", "selected"]
    assert x.loc["SMALL", "n_reference_samples"] == 20
    assert x.loc["SMALL", "notes"] == "line one line two"


def test_build_reference_availability_rejects_split_source_identity(tmp_path):
    _write_shard(tmp_path / "one.csv", "X", "SAME", ["E1"], 1)
    _write_shard(tmp_path / "two.csv", "X", "SAME", ["E2"], 1)

    with pytest.raises(ValueError, match="split across multiple shards"):
        build_reference_availability(tmp_path, chunksize=1)

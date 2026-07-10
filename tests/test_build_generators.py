# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

import glob
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from oncoref import expression_builders
from oncoref.cancer_types import cohort_registry

_SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"

# Optional real-data parity: the per-sample matrices + pirlygenes' shipped
# percentile artifact live only on a maintainer's machine (~22 GB cache), so this
# is gated like the HPA-dependent tests and skips cleanly everywhere else.
_ACC_MATRIX = glob.glob(
    os.path.expanduser("~/.cache/pirlygenes/expression/*/derived/tcga_acc_per_sample_tpm.parquet")
)
_ACC_REF = Path(
    os.path.expanduser(
        "~/code/pirlygenes/pirlygenes/data/cancer-reference-expression-percentiles/ACC.parquet"
    )
)
_PARITY_READY = bool(_ACC_MATRIX) and _ACC_REF.exists()


def _load_script(name):
    spec = importlib.util.spec_from_file_location(name, _SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _matrix(genes, samples, values):
    """genes × samples DataFrame with id columns."""
    df = pd.DataFrame(values, columns=samples)
    df.insert(0, "Symbol", [f"G{i}" for i in range(len(genes))])
    df.insert(0, "Ensembl_Gene_ID", genes)
    return df


# ---------- source-matrix ingestion builders ----------


def test_atomic_write_preserves_existing_artifact_on_failure(tmp_path):
    path = tmp_path / "artifact.csv"
    path.write_text("old\n")

    def _write_then_fail(tmp_path):
        tmp_path.write_text("new\n")
        raise RuntimeError("boom")

    with pytest.raises(RuntimeError, match="boom"):
        expression_builders._atomic_write(path, _write_then_fail)

    assert path.read_text() == "old\n"
    assert not list(tmp_path.glob("*.tmp"))


def test_geo_matrix_builder_writes_canonical_per_sample_matrix_and_sidecars(tmp_path):
    path = tmp_path / "geo.csv"
    pd.DataFrame(
        {
            "GeneID": ["ENSG00000141510.17", "ENSG00000141510", "ENSG00000146648"],
            "Symbol": ["TP53", "TP53", "EGFR"],
            "annotation": ["a", "b", "c"],
            "sample_1": ["2", "3", "5"],
            "sample_2": ["4", "6", "0"],
        }
    ).to_csv(path, index=False)
    source = expression_builders.GeoMatrixSource(
        cancer_code="X",
        source_cohort="TEST_GEO",
        source_project="GEO",
        file_name=path.name,
        unit="FPKM",
        gene_id_col="GeneID",
        symbol_col="Symbol",
        drop_cols=("annotation",),
        sep=",",
    )

    result = expression_builders.build_source_matrices(
        source,
        cache_dir=tmp_path,
        source_path=path,
    )

    out = pd.read_parquet(result.matrix_paths["X"])
    assert list(out.columns) == ["Ensembl_Gene_ID", "Symbol", "sample_1", "sample_2"]
    assert set(out["Ensembl_Gene_ID"]) == {"ENSG00000141510", "ENSG00000146648"}
    assert np.allclose(out[["sample_1", "sample_2"]].sum(axis=0), [1_000_000.0, 1_000_000.0])
    by_id = out.set_index("Ensembl_Gene_ID")
    assert np.isclose(by_id.loc["ENSG00000141510", "sample_1"], 500_000.0)
    assert np.isclose(by_id.loc["ENSG00000141510", "sample_2"], 1_000_000.0)

    stats = result.mapping_audit["mapping_status"].value_counts().to_dict()
    assert stats == {"resolved": 3}
    literal_zero = result.parse_diagnostics.set_index("value_col").loc["sample_2", "n_literal_zero"]
    assert literal_zero == 1
    assert result.sidecar_paths["mapping_audit"].exists()
    assert result.sidecar_paths["parse_diagnostics"].exists()
    assert result.sidecar_paths["X_sample_qc"].exists()
    assert result.sidecar_paths["summary_rows"].exists()
    assert set(result.sample_qc["sample_id"]) == {"sample_1", "sample_2"}
    assert set(result.sample_qc["source_cohort"]) == {"TEST_GEO"}
    summary = result.summary_rows.set_index("Ensembl_Gene_ID")
    assert list(result.summary_rows.columns) == list(
        expression_builders.REFERENCE_EXPRESSION_COLUMNS
    )
    assert set(summary["cancer_code"]) == {"X"}
    assert set(summary["source_cohort"]) == {"TEST_GEO"}
    assert set(summary["source_project"]) == {"GEO"}
    assert set(summary["tumor_origin"]) == {"primary"}
    assert summary.loc["ENSG00000141510", "n_samples"] == 2
    assert summary.loc["ENSG00000141510", "n_detected"] == 2
    assert summary.loc["ENSG00000141510", "TPM_median"] == pytest.approx(750_000.0)
    assert summary.loc["ENSG00000141510", "TPM_clean_median"] == pytest.approx(562_500.0)
    assert (
        summary.loc["ENSG00000141510", "processing_pipeline"]
        == "test_geo_fpkm_to_tpm_oncoref_canonical_clean_tpm_16_9_75"
    )


def test_geo_matrix_builder_reconciles_stale_per_code_artifacts(tmp_path):
    out_dir = tmp_path / "derived"
    out_dir.mkdir()
    stale_matrix = out_dir / "STALE_per_sample_tpm.parquet"
    stale_qc = out_dir / "STALE_sample_qc.csv"
    stale_matrix.write_text("stale matrix")
    stale_qc.write_text("stale qc")
    path = tmp_path / "geo.csv"
    pd.DataFrame(
        {
            "GeneID": ["ENSG00000141510", "ENSG00000146648"],
            "Symbol": ["TP53", "EGFR"],
            "sample_1": ["2", "3"],
        }
    ).to_csv(path, index=False)
    source = expression_builders.GeoMatrixSource(
        cancer_code="LIVE",
        source_cohort="TEST_GEO_STALE",
        source_project="GEO",
        file_name=path.name,
        unit="TPM",
        gene_id_col="GeneID",
        symbol_col="Symbol",
        sep=",",
    )

    result = expression_builders.build_source_matrices(
        source,
        cache_dir=tmp_path,
        source_path=path,
    )

    assert set(result.matrix_paths) == {"LIVE"}
    assert not stale_matrix.exists()
    assert not stale_qc.exists()
    assert (out_dir / "LIVE_per_sample_tpm.parquet").exists()
    assert (out_dir / "LIVE_sample_qc.csv").exists()


def test_geo_matrix_builder_preserves_stale_artifacts_when_no_samples_route(tmp_path):
    out_dir = tmp_path / "derived"
    out_dir.mkdir()
    stale_matrix = out_dir / "STALE_per_sample_tpm.parquet"
    stale_qc = out_dir / "STALE_sample_qc.csv"
    stale_matrix.write_text("stale matrix")
    stale_qc.write_text("stale qc")
    path = tmp_path / "geo.csv"
    pd.DataFrame(
        {
            "GeneID": ["ENSG00000141510"],
            "Symbol": ["TP53"],
            "sample_1": ["2"],
        }
    ).to_csv(path, index=False)
    source = expression_builders.GeoMatrixSource(
        cancer_code=["LIVE"],
        source_cohort="TEST_GEO_EMPTY",
        source_project="GEO",
        file_name=path.name,
        unit="TPM",
        gene_id_col="GeneID",
        symbol_col="Symbol",
        sep=",",
        sample_to_cancer_code=lambda _sample: None,
    )

    with pytest.raises(ValueError, match="no samples were routed"):
        expression_builders.build_source_matrices(
            source,
            cache_dir=tmp_path,
            source_path=path,
        )

    assert stale_matrix.exists()
    assert stale_qc.exists()


def test_geo_matrix_builder_routes_samples_and_reads_transposed_matrix(tmp_path):
    path = tmp_path / "transposed.tsv"
    pd.DataFrame(
        {
            "sample_id": ["tumor_a", "tumor_b"],
            "TP53": ["1", "3"],
            "EGFR": ["1", "1"],
        }
    ).to_csv(path, sep="\t", index=False)
    source = expression_builders.GeoMatrixSource(
        cancer_code=["CODE_A", "CODE_B"],
        source_cohort="TEST_TRANSPOSED",
        file_name=path.name,
        unit="TPM",
        gene_id_col="sample_id",
        transposed=True,
        sample_to_cancer_code=lambda sample: "CODE_A" if sample == "tumor_a" else "CODE_B",
    )

    result = expression_builders.build_source_matrices(
        source,
        cache_dir=tmp_path,
        source_path=path,
    )

    assert set(result.matrix_paths) == {"CODE_A", "CODE_B"}
    code_a = pd.read_parquet(result.matrix_paths["CODE_A"])
    code_b = pd.read_parquet(result.matrix_paths["CODE_B"])
    assert list(code_a.columns) == ["Ensembl_Gene_ID", "Symbol", "tumor_a"]
    assert list(code_b.columns) == ["Ensembl_Gene_ID", "Symbol", "tumor_b"]
    assert set(code_a["Symbol"]) == {"TP53", "EGFR"}
    assert np.isclose(code_a["tumor_a"].sum(), 1_000_000.0)
    assert np.isclose(code_b["tumor_b"].sum(), 1_000_000.0)


def test_geo_matrix_builder_transposed_all_blank_sample_is_missing_not_zero(tmp_path):
    path = tmp_path / "transposed_blank_sample.tsv"
    pd.DataFrame(
        {
            "sample_id": ["tumor_a", "blank_sample", "tumor_b"],
            "TP53": ["1", "", "3"],
            "EGFR": ["1", "", "1"],
        }
    ).to_csv(path, sep="\t", index=False)
    source = expression_builders.GeoMatrixSource(
        cancer_code=["CODE_A", "CODE_B", "CODE_BLANK"],
        source_cohort="TEST_TRANSPOSED_BLANK",
        file_name=path.name,
        unit="TPM",
        gene_id_col="sample_id",
        transposed=True,
        sample_to_cancer_code=lambda sample: {
            "tumor_a": "CODE_A",
            "tumor_b": "CODE_B",
            "blank_sample": "CODE_BLANK",
        }.get(sample),
    )

    result = expression_builders.build_source_matrices(
        source,
        cache_dir=tmp_path,
        source_path=path,
    )

    assert set(result.matrix_paths) == {"CODE_A", "CODE_B"}
    assert "CODE_BLANK" not in result.matrices
    assert "blank_sample" not in set(result.parse_diagnostics["value_col"])
    assert "blank_sample" not in set(result.sample_qc["sample_id"])
    for matrix in result.matrices.values():
        assert "blank_sample" not in matrix.columns
        assert np.allclose(
            matrix[expression_builders.sample_columns(matrix)].sum(axis=0), 1_000_000.0
        )


def test_source_matrix_unit_helpers_validate_raw_counts_lengths():
    df = pd.DataFrame({"gene_id": ["g1", "g2"], "s1": [10.0, 10.0]})
    out = expression_builders.normalize_source_matrix_to_tpm(
        df,
        unit="raw_counts",
        row_id_col="gene_id",
        gene_lengths_kb={"g1": 1.0, "g2": 2.0},
    )

    assert np.isclose(out["s1"].sum(), 1_000_000.0)
    assert np.isclose(out.loc[out["gene_id"] == "g1", "s1"].iloc[0], 666_666.6666666666)
    with pytest.raises(ValueError, match="gene_lengths_kb"):
        expression_builders.normalize_source_matrix_to_tpm(df, unit="raw_counts")


def test_source_matrix_builder_emits_summary_rows_for_raw_counts(tmp_path):
    path = tmp_path / "raw_counts.csv"
    pd.DataFrame(
        {
            "GeneID": ["ENSG00000141510", "ENSG00000146648"],
            "Symbol": ["TP53", "EGFR"],
            "sample_1": [10, 10],
            "sample_2": [30, 10],
        }
    ).to_csv(path, index=False)
    source = expression_builders.GeoMatrixSource(
        cancer_code="RAW",
        source_cohort="TEST_RAW",
        source_project="GEO",
        citation="PMID:1",
        file_name=path.name,
        unit="raw_counts",
        gene_id_col="GeneID",
        symbol_col="Symbol",
        sep=",",
        pipeline_stem="test_raw",
        notes="raw count source notes",
        tumor_origin="metastasis",
        metastasis_site="liver",
    )

    result = expression_builders.build_source_matrices(
        source,
        cache_dir=tmp_path,
        source_path=path,
        gene_lengths_kb={"ENSG00000141510": 1.0, "ENSG00000146648": 2.0},
    )

    summary = result.summary_rows.set_index("Symbol")
    assert set(summary["cancer_code"]) == {"RAW"}
    assert set(summary["notes"]) == {"raw count source notes"}
    assert set(summary["tumor_origin"]) == {"metastasis"}
    assert set(summary["metastasis_site"]) == {"liver"}
    assert set(summary["processing_pipeline"]) == {
        "test_raw_raw_counts_to_tpm_oncoref_canonical_clean_tpm_16_9_75"
    }
    assert summary.loc["TP53", "source_version"].startswith("PMID:1; unit=raw_counts")
    assert summary.loc["TP53", "n_samples"] == 2
    assert summary.loc["TP53", "n_detected"] == 2
    assert summary.loc["TP53", "TPM_median"] > summary.loc["EGFR", "TPM_median"]
    assert pd.read_csv(result.sidecar_paths["summary_rows"]).shape == result.summary_rows.shape


def test_summarize_source_matrix_matches_reference_stat_contract():
    matrix = pd.DataFrame(
        {
            "Ensembl_Gene_ID": ["ENSG00000141510", "ENSG00000146648"],
            "Symbol": ["TP53", "EGFR"],
            "s0": [0.0, 10.0],
            "s1": [1.0, 10.0],
            "s2": [2.0, 10.0],
            "s3": [3.0, 10.0],
            "s4": [4.0, 10.0],
        }
    )
    source = expression_builders.GeoMatrixSource(
        cancer_code="X",
        source_cohort="TEST_STATS",
        file_name="unused.tsv",
        unit="TPM",
    )

    summary = expression_builders.summarize_source_matrix(
        matrix,
        cancer_code="X",
        source=source,
    ).set_index("Symbol")

    assert summary.loc["TP53", "TPM_median"] == 2.0
    assert summary.loc["TP53", "TPM_mean"] == 2.0
    assert summary.loc["TP53", "TPM_q1"] == 1.0
    assert summary.loc["TP53", "TPM_q3"] == 3.0
    assert summary.loc["TP53", "TPM_min"] == 0.0
    assert summary.loc["TP53", "TPM_max"] == 4.0
    assert summary.loc["TP53", "TPM_std"] == pytest.approx(round(np.sqrt(2.5), 6))
    assert summary.loc["TP53", "TPM_p10"] == 0.4
    assert summary.loc["TP53", "TPM_p90"] == 3.6
    assert summary.loc["TP53", "n_samples"] == 5
    assert summary.loc["TP53", "n_detected"] == 4


def test_geo_matrix_source_from_entry_compiles_yaml_filters_and_routing():
    source = expression_builders.geo_matrix_source_from_entry(
        {
            "id": "synthetic",
            "source_type": "geo-matrix",
            "cancer_codes": ["CODE_A", "CODE_B"],
            "source_cohort": "SYNTHETIC",
            "file_url": "https://example.org/source.tsv.gz",
            "file_name": "source.tsv.gz",
            "unit": "log2-TPM",
            "gene_id_col": "",
            "sample_filter": {"include_match": "tumor", "exclude_match": "bad"},
            "sample_to_cancer_code": {
                "rules": [
                    {"match": "^tumor_a", "cancer_code": "CODE_A"},
                    {"match": "^tumor_b", "cancer_code": "CODE_B"},
                ]
            },
            "notes": "source-level summary row notes",
            "pipeline_stem": "synthetic_pipeline",
            "tumor_origin": "metastasis",
            "metastasis_site": "liver",
        }
    )

    assert source.cancer_code == ["CODE_A", "CODE_B"]
    assert source.unit == "log2(TPM+1)"
    assert source.sample_filter(["tumor_a1", "normal_a1", "tumor_bad", "tumor_b1"]) == [
        "tumor_a1",
        "tumor_b1",
    ]
    assert source.sample_to_cancer_code("tumor_a1") == "CODE_A"
    assert source.sample_to_cancer_code("tumor_b1") == "CODE_B"
    assert source.sample_to_cancer_code("normal_a1") is None
    assert source.notes == "source-level summary row notes"
    assert source.pipeline_stem == "synthetic_pipeline"
    assert source.tumor_origin == "metastasis"
    assert source.metastasis_site == "liver"


def test_geo_matrix_source_from_entry_validates_tumor_origin():
    with pytest.raises(ValueError, match="tumor_origin"):
        expression_builders.geo_matrix_source_from_entry(
            {
                "id": "synthetic",
                "source_type": "geo-matrix",
                "cancer_codes": ["CODE_A"],
                "source_cohort": "SYNTHETIC",
                "file_url": "https://example.org/source.tsv.gz",
                "file_name": "source.tsv.gz",
                "unit": "TPM",
                "tumor_origin": "metastatic",
            }
        )


def test_geo_matrix_source_from_registry_loads_packaged_geo_entry():
    source = expression_builders.geo_matrix_source_from_registry("gse328026-sarc-pec")

    assert source.cancer_code == "SARC_PEC"
    assert source.source_cohort == "GSE328026_PECOMA_2026"
    assert source.unit == "TPM"
    assert source.file_name == "GSE328026_TPMs_all_Samples.txt.gz"
    assert source.notes.startswith("n=69 PEComa tumors")
    assert source.pipeline_stem == ""
    assert source.tumor_origin == "primary"
    assert source.metastasis_site is None


def test_recount3_source_from_registry_loads_packaged_routes():
    source = expression_builders.recount3_source_from_registry("gse98894-midnet")

    assert source.source_id == "gse98894-midnet"
    assert source.srp == "SRP107025"
    assert source.cancer_code == ["NET_MIDGUT", "NET_PANCREAS", "NET_RECTAL"]
    assert source.expected_n == {"NET_MIDGUT": 81, "NET_PANCREAS": 113, "NET_RECTAL": 18}
    assert source.tumor_origin == "primary"
    assert source.sample_to_cancer_code({"origin": "ileum"}, "") == "NET_MIDGUT"
    assert source.sample_to_cancer_code({"origin": "pancreas"}, "") == "NET_PANCREAS"
    assert source.sample_to_cancer_code({"origin": "rectal"}, "") == "NET_RECTAL"
    assert source.sample_to_cancer_code({"origin": "lung"}, "") is None


def test_treehouse_source_from_registry_loads_direct_cohort_routes():
    source = expression_builders.treehouse_source_from_registry("treehouse-polya-25-01")

    assert source.source_id == "treehouse-polya-25-01"
    assert source.source_cohort == "TREEHOUSE_POLYA_25_01"
    assert source.tpm_file.startswith("Tumor-25.01-Polya")
    assert source.clinical_file.startswith("clinical_Treehouse")
    assert len(source.cohorts) == 26
    by_code = {cohort.cancer_code: cohort for cohort in source.cohorts}
    assert by_code["SARC_EWS"].disease_label == "Ewing sarcoma"
    assert by_code["SARC_MPNST"].group == "sarc_rare_direct"
    assert by_code["SARC_GIST"].group == "sarc_subtypes"

    rare = expression_builders.treehouse_cohorts_for_group("sarc_rare_direct")
    assert "SARC_MPNST" in {cohort.cancer_code for cohort in rare}

    ribod = expression_builders.treehouse_source_from_registry("treehouse-ribod-25-01")
    assert [cohort.cancer_code for cohort in ribod.cohorts] == ["SARC_CHOR", "RB"]


def test_treehouse_source_from_registry_loads_tcga_subset_routes():
    source = expression_builders.treehouse_source_from_registry("treehouse-polya-25-01-tcga-subset")

    assert source.source_cohort == "TREEHOUSE_POLYA_25_01_TCGA_SUBSET"
    assert source.source_project == "Treehouse (TCGA samples)"
    assert source.pipeline_stem == "treehouse_polya_25_01_tcga_subset"
    assert isinstance(source.cancer_code, list)
    assert len(source.cancer_code) == 30
    assert len(source.cohorts) == 30
    by_code = {cohort.cancer_code: cohort for cohort in source.cohorts}
    assert by_code["BRCA"].disease_label == "breast invasive carcinoma"
    assert by_code["BRCA"].selection == "tcga"
    assert by_code["BRCA"].effective_cache_stem == "tcga_brca"
    assert by_code["UCEC"].disease_label == "uterine corpus endometrioid carcinoma"
    assert "GBM" not in by_code
    assert "LGG" not in by_code
    assert "SARC" not in by_code

    tcga_direct = expression_builders.treehouse_cohorts_for_group(
        "tcga_direct",
        source_id="treehouse-polya-25-01-tcga-subset",
    )
    assert [cohort.cancer_code for cohort in tcga_direct] == source.cancer_code
    assert {cohort.selection for cohort in tcga_direct} == {"tcga"}


def test_treehouse_source_from_registry_loads_brca_pam50_routes():
    source = expression_builders.treehouse_source_from_registry(
        "treehouse-polya-25-01-tcga-brca-pam50"
    )

    assert source.source_cohort == "TREEHOUSE_POLYA_25_01_TCGA_BRCA_PAM50"
    assert source.pipeline_stem == "treehouse_polya_25_01_tcga_brca_pam50"
    assert source.cancer_code == [
        "BRCA_Basal",
        "BRCA_HER2",
        "BRCA_LumA",
        "BRCA_LumB",
        "BRCA_Normal",
    ]
    by_code = {cohort.cancer_code: cohort for cohort in source.cohorts}
    assert by_code["BRCA_Basal"].disease_label == "breast invasive carcinoma"
    assert (
        by_code["BRCA_Basal"].selection
        == "cbio_clinical:brca_tcga_pan_can_atlas_2018:SUBTYPE:BRCA_Basal"
    )
    assert (
        by_code["BRCA_HER2"].selection
        == "cbio_clinical:brca_tcga_pan_can_atlas_2018:SUBTYPE:BRCA_Her2"
    )
    assert by_code["BRCA_HER2"].effective_cache_stem == "tcga_brca_her2"

    pam50 = expression_builders.treehouse_cohorts_for_group(
        "tcga_brca_pam50",
        source_id="treehouse-polya-25-01-tcga-brca-pam50",
    )
    assert [cohort.cancer_code for cohort in pam50] == source.cancer_code


def test_treehouse_source_from_registry_loads_hnsc_hpv_routes():
    source = expression_builders.treehouse_source_from_registry(
        "treehouse-polya-25-01-tcga-hnsc-hpv"
    )

    assert source.source_cohort == "TREEHOUSE_POLYA_25_01_TCGA_HNSC_HPV"
    assert source.pipeline_stem == "treehouse_polya_25_01_tcga_hnsc_hpv"
    assert source.cancer_code == ["HNSC_HPVneg", "HNSC_HPVpos"]
    by_code = {cohort.cancer_code: cohort for cohort in source.cohorts}
    assert by_code["HNSC_HPVneg"].disease_label == "head & neck squamous cell carcinoma"
    assert (
        by_code["HNSC_HPVneg"].selection
        == "cbio_clinical:hnsc_tcga_pan_can_atlas_2018:SUBTYPE:HNSC_HPV-"
    )
    assert (
        by_code["HNSC_HPVpos"].selection
        == "cbio_clinical:hnsc_tcga_pan_can_atlas_2018:SUBTYPE:HNSC_HPV+"
    )
    assert by_code["HNSC_HPVpos"].effective_cache_stem == "tcga_hnsc_hpvpos"

    hpv = expression_builders.treehouse_cohorts_for_group(
        "tcga_hnsc_hpv",
        source_id="treehouse-polya-25-01-tcga-hnsc-hpv",
    )
    assert [cohort.cancer_code for cohort in hpv] == source.cancer_code


def test_treehouse_source_from_registry_loads_luad_mutation_routes():
    source = expression_builders.treehouse_source_from_registry(
        "treehouse-polya-25-01-tcga-luad-mut"
    )

    assert source.source_cohort == "TREEHOUSE_POLYA_25_01_TCGA_LUAD_MUT"
    assert source.pipeline_stem == "treehouse_polya_25_01_tcga_luad_mut"
    assert source.cancer_code == ["LUAD_EGFR", "LUAD_KRAS", "LUAD_STK11"]
    by_code = {cohort.cancer_code: cohort for cohort in source.cohorts}
    assert by_code["LUAD_EGFR"].disease_label == "lung adenocarcinoma"
    assert by_code["LUAD_EGFR"].selection == "cbio_mutation:luad_tcga_pan_can_atlas_2018:EGFR"
    assert (
        by_code["LUAD_STK11"].selection == "cbio_mutation:luad_tcga_pan_can_atlas_2018:STK11,KEAP1"
    )
    assert by_code["LUAD_STK11"].effective_cache_stem == "tcga_luad_stk11"

    luad_mut = expression_builders.treehouse_cohorts_for_group(
        "tcga_luad_mut",
        source_id="treehouse-polya-25-01-tcga-luad-mut",
    )
    assert [cohort.cancer_code for cohort in luad_mut] == source.cancer_code

    registry = cohort_registry()
    assert registry["TREEHOUSE_POLYA_25_01_TCGA_LUAD_MUT"]["n_samples"] == 362


def test_treehouse_source_from_registry_loads_ucec_subtype_routes():
    source = expression_builders.treehouse_source_from_registry(
        "treehouse-polya-25-01-tcga-ucec-subtype"
    )

    assert source.source_cohort == "TREEHOUSE_POLYA_25_01_TCGA_UCEC_SUBTYPE"
    assert source.pipeline_stem == "treehouse_polya_25_01_tcga_ucec_subtype"
    assert source.cancer_code == ["UCEC_POLE", "UCEC_MSI", "UCEC_CNL", "UCEC_CNH"]
    by_code = {cohort.cancer_code: cohort for cohort in source.cohorts}
    assert by_code["UCEC_POLE"].disease_label == "endometrial carcinoma"
    assert (
        by_code["UCEC_POLE"].selection
        == "cbio_clinical:ucec_tcga_pan_can_atlas_2018:SUBTYPE:UCEC_POLE"
    )
    assert (
        by_code["UCEC_CNL"].selection
        == "cbio_clinical:ucec_tcga_pan_can_atlas_2018:SUBTYPE:UCEC_CN_LOW"
    )
    assert (
        by_code["UCEC_CNH"].selection
        == "cbio_clinical:ucec_tcga_pan_can_atlas_2018:SUBTYPE:UCEC_CN_HIGH"
    )
    assert by_code["UCEC_CNH"].effective_cache_stem == "tcga_ucec_cnh"

    subtypes = expression_builders.treehouse_cohorts_for_group(
        "tcga_ucec_subtype",
        source_id="treehouse-polya-25-01-tcga-ucec-subtype",
    )
    assert [cohort.cancer_code for cohort in subtypes] == source.cancer_code

    registry = cohort_registry()
    assert registry["TREEHOUSE_POLYA_25_01_TCGA_UCEC_SUBTYPE"]["n_samples"] == 172


def test_treehouse_source_from_registry_loads_stad_subtype_routes():
    source = expression_builders.treehouse_source_from_registry(
        "treehouse-polya-25-01-tcga-stad-subtype"
    )

    assert source.source_cohort == "TREEHOUSE_POLYA_25_01_TCGA_STAD_SUBTYPE"
    assert source.pipeline_stem == "treehouse_polya_25_01_tcga_stad_subtype"
    assert source.cancer_code == ["STAD_EBV", "STAD_MSI", "STAD_GS", "STAD_CIN"]
    by_code = {cohort.cancer_code: cohort for cohort in source.cohorts}
    assert by_code["STAD_EBV"].disease_label == "stomach adenocarcinoma"
    assert (
        by_code["STAD_EBV"].selection
        == "cbio_clinical:stad_tcga_pan_can_atlas_2018:SUBTYPE:STAD_EBV"
    )
    assert (
        by_code["STAD_CIN"].selection
        == "cbio_clinical:stad_tcga_pan_can_atlas_2018:SUBTYPE:STAD_CIN"
    )
    assert by_code["STAD_CIN"].effective_cache_stem == "tcga_stad_cin"

    subtypes = expression_builders.treehouse_cohorts_for_group(
        "tcga_stad_subtype",
        source_id="treehouse-polya-25-01-tcga-stad-subtype",
    )
    assert [cohort.cancer_code for cohort in subtypes] == source.cancer_code

    registry = cohort_registry()
    assert registry["TREEHOUSE_POLYA_25_01_TCGA_STAD_SUBTYPE"]["n_samples"] == 374


def test_treehouse_source_from_registry_loads_coadread_msi_routes():
    source = expression_builders.treehouse_source_from_registry(
        "treehouse-polya-25-01-tcga-coadread-msi"
    )

    assert source.source_cohort == "TREEHOUSE_POLYA_25_01_TCGA_COADREAD_MSI"
    assert source.pipeline_stem == "treehouse_polya_25_01_tcga_coadread_msi"
    assert source.cancer_code == ["COAD_MSI", "COAD_MSS", "READ_MSI", "READ_MSS"]
    by_code = {cohort.cancer_code: cohort for cohort in source.cohorts}
    assert by_code["COAD_MSI"].disease_label == "colon adenocarcinoma"
    assert by_code["READ_MSI"].disease_label == "rectum adenocarcinoma"
    assert (
        by_code["COAD_MSI"].selection
        == "cbio_sample_clinical:coadread_tcga_pan_can_atlas_2018:MSI_SENSOR_SCORE:>=10"
    )
    assert (
        by_code["COAD_MSS"].selection
        == "cbio_sample_clinical:coadread_tcga_pan_can_atlas_2018:MSI_SENSOR_SCORE:<10"
    )
    assert by_code["READ_MSS"].effective_cache_stem == "tcga_read_mss"

    subtypes = expression_builders.treehouse_cohorts_for_group(
        "tcga_coadread_msi",
        source_id="treehouse-polya-25-01-tcga-coadread-msi",
    )
    assert [cohort.cancer_code for cohort in subtypes] == source.cancer_code

    registry = cohort_registry()
    assert registry["TREEHOUSE_POLYA_25_01_TCGA_COADREAD_MSI"]["n_samples"] == 361


def test_treehouse_source_from_registry_loads_glioma_gdc_project_routes():
    source = expression_builders.treehouse_source_from_registry("treehouse-polya-25-01-tcga-glioma")

    assert source.source_cohort == "TREEHOUSE_POLYA_25_01_TCGA_SUBSET"
    assert source.pipeline_stem == "treehouse_polya_25_01_tcga_glioma_split"
    assert source.cancer_code == ["GBM", "LGG"]
    by_code = {cohort.cancer_code: cohort for cohort in source.cohorts}
    assert by_code["GBM"].disease_label == "glioma"
    assert by_code["GBM"].selection == "gdc_project:TCGA-GBM"
    assert by_code["GBM"].effective_cache_stem == "tcga_gbm"
    assert by_code["LGG"].selection == "gdc_project:TCGA-LGG"

    glioma = expression_builders.treehouse_cohorts_for_group(
        "tcga_glioma",
        source_id="treehouse-polya-25-01-tcga-glioma",
    )
    assert [cohort.cancer_code for cohort in glioma] == ["GBM", "LGG"]


def test_treehouse_sample_ids_filter_disease_and_tcga_selection():
    clinical = pd.DataFrame(
        [
            {"th_dataset_id": "TREEHOUSE-1", "disease": "Ewing sarcoma"},
            {"th_dataset_id": "TREEHOUSE-2", "disease": "ewing sarcoma"},
            {"th_dataset_id": "TCGA-AB-1234-01A", "disease": "Ewing sarcoma"},
            {"th_dataset_id": "TCGA-XY-9999-01A", "disease": "osteosarcoma"},
        ]
    )
    cohort = expression_builders.TreehouseCohort("SARC_EWS", "Ewing sarcoma")
    assert expression_builders.treehouse_sample_ids(clinical, cohort) == [
        "TREEHOUSE-1",
        "TREEHOUSE-2",
        "TCGA-AB-1234-01A",
    ]

    tcga = expression_builders.TreehouseCohort(
        "SARC_EWS",
        "Ewing sarcoma",
        selection="tcga",
    )
    assert expression_builders.treehouse_sample_ids(clinical, tcga) == ["TCGA-AB-1234-01A"]

    unsupported = expression_builders.TreehouseCohort(
        "SARC_EWS",
        "Ewing sarcoma",
        selection="legacy_pam50:BRCA_Basal",
    )
    with pytest.raises(ValueError, match="unsupported Treehouse cohort selection"):
        expression_builders.treehouse_sample_ids(clinical, unsupported)


def test_treehouse_sample_ids_filter_gdc_project_selection():
    clinical = pd.DataFrame(
        [
            {"th_dataset_id": "TCGA-GB-0001-01A", "disease": "glioma"},
            {"th_dataset_id": "TCGA-LG-0002-01A", "disease": "glioma"},
            {"th_dataset_id": "TREEHOUSE-3", "disease": "glioma"},
        ]
    )
    cohort = expression_builders.TreehouseCohort(
        "GBM",
        "glioma",
        selection="gdc_project:TCGA-GBM",
    )
    case_sets = {"gdc_project:TCGA-GBM": {"TCGA-GB-0001"}}

    assert expression_builders.treehouse_sample_ids(
        clinical,
        cohort,
        selection_case_sets=case_sets,
    ) == ["TCGA-GB-0001-01A"]

    with pytest.raises(ValueError, match="requires a precomputed GDC case set"):
        expression_builders.treehouse_sample_ids(clinical, cohort)


def test_treehouse_sample_ids_filter_cbio_clinical_selection():
    clinical = pd.DataFrame(
        [
            {"th_dataset_id": "TCGA-BR-0001-01A", "disease": "breast invasive carcinoma"},
            {"th_dataset_id": "TCGA-BR-0002-01A", "disease": "breast invasive carcinoma"},
            {"th_dataset_id": "TREEHOUSE-3", "disease": "breast invasive carcinoma"},
        ]
    )
    selector = "cbio_clinical:brca_tcga_pan_can_atlas_2018:SUBTYPE:BRCA_Basal"
    cohort = expression_builders.TreehouseCohort(
        "BRCA_Basal",
        "breast invasive carcinoma",
        selection=selector,
    )
    case_sets = {selector: {"TCGA-BR-0001"}}

    assert expression_builders.treehouse_sample_ids(
        clinical,
        cohort,
        selection_case_sets=case_sets,
    ) == ["TCGA-BR-0001-01A"]

    with pytest.raises(ValueError, match="requires a precomputed cBioPortal case set"):
        expression_builders.treehouse_sample_ids(clinical, cohort)


def test_treehouse_sample_ids_filter_cbio_sample_clinical_numeric_selection():
    clinical = pd.DataFrame(
        [
            {"th_dataset_id": "TCGA-CO-0001-01A", "disease": "colon adenocarcinoma"},
            {"th_dataset_id": "TCGA-CO-0002-01A", "disease": "colon adenocarcinoma"},
            {"th_dataset_id": "TCGA-RE-0003-01A", "disease": "rectum adenocarcinoma"},
        ]
    )
    selector = "cbio_sample_clinical:coadread_tcga_pan_can_atlas_2018:MSI_SENSOR_SCORE:>=10"
    cohort = expression_builders.TreehouseCohort(
        "COAD_MSI",
        "colon adenocarcinoma",
        selection=selector,
    )
    case_sets = {selector: {"TCGA-CO-0001", "TCGA-RE-0003"}}

    assert expression_builders.treehouse_sample_ids(
        clinical,
        cohort,
        selection_case_sets=case_sets,
    ) == ["TCGA-CO-0001-01A"]

    with pytest.raises(ValueError, match="requires a precomputed cBioPortal case set"):
        expression_builders.treehouse_sample_ids(clinical, cohort)


def test_treehouse_sample_ids_filter_cbio_mutation_selection():
    clinical = pd.DataFrame(
        [
            {"th_dataset_id": "TCGA-LU-0001-01A", "disease": "lung adenocarcinoma"},
            {"th_dataset_id": "TCGA-LU-0002-01A", "disease": "lung adenocarcinoma"},
            {"th_dataset_id": "TREEHOUSE-3", "disease": "lung adenocarcinoma"},
        ]
    )
    selector = "cbio_mutation:luad_tcga_pan_can_atlas_2018:STK11,KEAP1"
    cohort = expression_builders.TreehouseCohort(
        "LUAD_STK11",
        "lung adenocarcinoma",
        selection=selector,
    )
    case_sets = {selector: {"TCGA-LU-0002"}}

    assert expression_builders.treehouse_sample_ids(
        clinical,
        cohort,
        selection_case_sets=case_sets,
    ) == ["TCGA-LU-0002-01A"]

    with pytest.raises(ValueError, match="requires a precomputed cBioPortal case set"):
        expression_builders.treehouse_sample_ids(clinical, cohort)


def test_build_treehouse_source_matrices_writes_canonical_artifacts(tmp_path):
    clinical_path = tmp_path / "clinical.tsv"
    pd.DataFrame(
        [
            {"th_dataset_id": "SAMPLE_A", "disease": "synthetic tumor"},
            {"th_dataset_id": "SAMPLE_B", "disease": "Synthetic Tumor"},
            {"th_dataset_id": "CONTROL", "disease": "other"},
        ]
    ).to_csv(clinical_path, sep="\t", index=False)
    tpm_path = tmp_path / "treehouse.tsv"
    pd.DataFrame(
        {
            "Gene": ["TP53", "EGFR"],
            "SAMPLE_A": np.log2(np.array([2.0, 1.0]) + 1.0),
            "SAMPLE_B": np.log2(np.array([4.0, 0.0]) + 1.0),
            "CONTROL": np.log2(np.array([100.0, 100.0]) + 1.0),
        }
    ).to_csv(tpm_path, sep="\t", index=False)
    source = expression_builders.TreehouseSource(
        source_id="synthetic-treehouse",
        source_cohort="SYNTHETIC_TREEHOUSE",
        cancer_code="CODE_A",
        tpm_file=tpm_path.name,
        clinical_file=clinical_path.name,
        source_project="Treehouse",
        cohorts=(
            expression_builders.TreehouseCohort(
                "CODE_A",
                "synthetic tumor",
                extra_notes="Synthetic cohort note.",
            ),
        ),
        notes="Synthetic source note.",
        pipeline_stem="synthetic_treehouse",
    )

    result = expression_builders.build_treehouse_source_matrices(source, cache_dir=tmp_path)

    assert set(result.matrix_paths) == {"CODE_A"}
    out = pd.read_parquet(result.matrix_paths["CODE_A"]).set_index("Symbol")
    assert list(out.columns) == ["Ensembl_Gene_ID", "SAMPLE_A", "SAMPLE_B"]
    np.testing.assert_allclose(
        out.loc["TP53", ["SAMPLE_A", "SAMPLE_B"]].astype(float).to_numpy(),
        [2.0, 4.0],
    )
    np.testing.assert_allclose(
        out.loc["EGFR", ["SAMPLE_A", "SAMPLE_B"]].astype(float).to_numpy(),
        [1.0, 0.0],
    )
    assert result.sidecar_paths["mapping_audit"].exists()
    assert result.sidecar_paths["parse_diagnostics"].exists()
    assert result.sidecar_paths["summary_rows"].exists()
    assert result.summary_rows["notes"].str.contains("Synthetic cohort note").all()
    summary = result.summary_rows.set_index("Symbol")
    assert summary.loc["TP53", "TPM_median"] == 3.0
    assert set(result.sample_qc["sample_id"]) == {"SAMPLE_A", "SAMPLE_B"}


def test_build_treehouse_source_matrices_splits_gdc_project_cohorts(
    tmp_path,
    monkeypatch,
):
    clinical_path = tmp_path / "clinical.tsv"
    pd.DataFrame(
        [
            {"th_dataset_id": "TCGA-GB-0001-01A", "disease": "glioma"},
            {"th_dataset_id": "TCGA-LG-0002-01A", "disease": "glioma"},
            {"th_dataset_id": "TREEHOUSE-3", "disease": "glioma"},
        ]
    ).to_csv(clinical_path, sep="\t", index=False)
    tpm_path = tmp_path / "treehouse.tsv"
    pd.DataFrame(
        {
            "Gene": ["TP53", "EGFR"],
            "TCGA-GB-0001-01A": np.log2(np.array([2.0, 1.0]) + 1.0),
            "TCGA-LG-0002-01A": np.log2(np.array([4.0, 0.0]) + 1.0),
            "TREEHOUSE-3": np.log2(np.array([100.0, 100.0]) + 1.0),
        }
    ).to_csv(tpm_path, sep="\t", index=False)
    source = expression_builders.TreehouseSource(
        source_id="synthetic-treehouse-glioma",
        source_cohort="SYNTHETIC_TREEHOUSE_TCGA",
        cancer_code=["GBM", "LGG"],
        tpm_file=tpm_path.name,
        clinical_file=clinical_path.name,
        source_project="Treehouse (TCGA samples)",
        cohorts=(
            expression_builders.TreehouseCohort(
                "GBM",
                "glioma",
                selection="gdc_project:TCGA-GBM",
                cache_stem="tcga_gbm",
            ),
            expression_builders.TreehouseCohort(
                "LGG",
                "glioma",
                selection="gdc_project:TCGA-LGG",
                cache_stem="tcga_lgg",
            ),
        ),
    )

    def fake_case_map(project_ids, *, cache_path=None, force_download=False):
        assert set(project_ids) == {"TCGA-GBM", "TCGA-LGG"}
        assert cache_path is not None
        return pd.DataFrame(
            {
                "submitter_id": ["TCGA-GB-0001", "TCGA-LG-0002"],
                "project_id": ["TCGA-GBM", "TCGA-LGG"],
            }
        )

    monkeypatch.setattr(
        expression_builders,
        "treehouse_gdc_case_project_map",
        fake_case_map,
    )
    result = expression_builders.build_treehouse_source_matrices(source, cache_dir=tmp_path)

    assert set(result.matrix_paths) == {"GBM", "LGG"}
    gbm = pd.read_parquet(result.matrix_paths["GBM"]).set_index("Symbol")
    lgg = pd.read_parquet(result.matrix_paths["LGG"]).set_index("Symbol")
    assert list(gbm.columns) == ["Ensembl_Gene_ID", "TCGA-GB-0001-01A"]
    assert list(lgg.columns) == ["Ensembl_Gene_ID", "TCGA-LG-0002-01A"]
    np.testing.assert_allclose(gbm.loc["TP53", "TCGA-GB-0001-01A"], 2.0)
    np.testing.assert_allclose(lgg.loc["TP53", "TCGA-LG-0002-01A"], 4.0)


def test_build_treehouse_source_matrices_splits_cbio_clinical_cohorts(
    tmp_path,
    monkeypatch,
):
    clinical_path = tmp_path / "clinical.tsv"
    pd.DataFrame(
        [
            {"th_dataset_id": "TCGA-BR-0001-01A", "disease": "breast invasive carcinoma"},
            {"th_dataset_id": "TCGA-BR-0002-01A", "disease": "breast invasive carcinoma"},
            {"th_dataset_id": "TREEHOUSE-3", "disease": "breast invasive carcinoma"},
        ]
    ).to_csv(clinical_path, sep="\t", index=False)
    tpm_path = tmp_path / "treehouse.tsv"
    pd.DataFrame(
        {
            "Gene": ["TP53", "ERBB2"],
            "TCGA-BR-0001-01A": np.log2(np.array([2.0, 1.0]) + 1.0),
            "TCGA-BR-0002-01A": np.log2(np.array([4.0, 8.0]) + 1.0),
            "TREEHOUSE-3": np.log2(np.array([100.0, 100.0]) + 1.0),
        }
    ).to_csv(tpm_path, sep="\t", index=False)
    source = expression_builders.TreehouseSource(
        source_id="synthetic-treehouse-brca-pam50",
        source_cohort="SYNTHETIC_TREEHOUSE_BRCA_PAM50",
        cancer_code=["BRCA_Basal", "BRCA_HER2"],
        tpm_file=tpm_path.name,
        clinical_file=clinical_path.name,
        source_project="Treehouse (TCGA-BRCA) x cBioPortal PAM50",
        cohorts=(
            expression_builders.TreehouseCohort(
                "BRCA_Basal",
                "breast invasive carcinoma",
                selection="cbio_clinical:brca_tcga_pan_can_atlas_2018:SUBTYPE:BRCA_Basal",
                cache_stem="tcga_brca_basal",
            ),
            expression_builders.TreehouseCohort(
                "BRCA_HER2",
                "breast invasive carcinoma",
                selection="cbio_clinical:brca_tcga_pan_can_atlas_2018:SUBTYPE:BRCA_Her2",
                cache_stem="tcga_brca_her2",
            ),
        ),
    )

    def fake_clinical_map(
        study_id,
        attribute_id,
        *,
        clinical_data_type="PATIENT",
        cache_path=None,
        force_download=False,
    ):
        assert study_id == "brca_tcga_pan_can_atlas_2018"
        assert attribute_id == "SUBTYPE"
        assert clinical_data_type == "PATIENT"
        assert cache_path is not None
        return pd.DataFrame(
            {
                "case_id": ["TCGA-BR-0001", "TCGA-BR-0002"],
                "value": ["BRCA_Basal", "BRCA_Her2"],
            }
        )

    monkeypatch.setattr(
        expression_builders,
        "treehouse_cbioportal_clinical_attribute_map",
        fake_clinical_map,
    )
    result = expression_builders.build_treehouse_source_matrices(source, cache_dir=tmp_path)

    assert set(result.matrix_paths) == {"BRCA_Basal", "BRCA_HER2"}
    basal = pd.read_parquet(result.matrix_paths["BRCA_Basal"]).set_index("Symbol")
    her2 = pd.read_parquet(result.matrix_paths["BRCA_HER2"]).set_index("Symbol")
    assert list(basal.columns) == ["Ensembl_Gene_ID", "TCGA-BR-0001-01A"]
    assert list(her2.columns) == ["Ensembl_Gene_ID", "TCGA-BR-0002-01A"]
    np.testing.assert_allclose(basal.loc["TP53", "TCGA-BR-0001-01A"], 2.0)
    np.testing.assert_allclose(her2.loc["ERBB2", "TCGA-BR-0002-01A"], 8.0)


def test_build_treehouse_source_matrices_splits_cbio_sample_numeric_cohorts(
    tmp_path,
    monkeypatch,
):
    clinical_path = tmp_path / "clinical.tsv"
    pd.DataFrame(
        [
            {"th_dataset_id": "TCGA-CO-0001-01A", "disease": "colon adenocarcinoma"},
            {"th_dataset_id": "TCGA-CO-0002-01A", "disease": "colon adenocarcinoma"},
            {"th_dataset_id": "TCGA-RE-0003-01A", "disease": "rectum adenocarcinoma"},
        ]
    ).to_csv(clinical_path, sep="\t", index=False)
    tpm_path = tmp_path / "treehouse.tsv"
    pd.DataFrame(
        {
            "Gene": ["TP53", "MLH1"],
            "TCGA-CO-0001-01A": np.log2(np.array([2.0, 1.0]) + 1.0),
            "TCGA-CO-0002-01A": np.log2(np.array([4.0, 8.0]) + 1.0),
            "TCGA-RE-0003-01A": np.log2(np.array([6.0, 3.0]) + 1.0),
        }
    ).to_csv(tpm_path, sep="\t", index=False)
    source = expression_builders.TreehouseSource(
        source_id="synthetic-treehouse-coadread-msi",
        source_cohort="SYNTHETIC_TREEHOUSE_COADREAD_MSI",
        cancer_code=["COAD_MSI", "COAD_MSS"],
        tpm_file=tpm_path.name,
        clinical_file=clinical_path.name,
        source_project="Treehouse (TCGA-COAD/READ) x cBioPortal MSIsensor",
        cohorts=(
            expression_builders.TreehouseCohort(
                "COAD_MSI",
                "colon adenocarcinoma",
                selection="cbio_sample_clinical:coadread_tcga_pan_can_atlas_2018:MSI_SENSOR_SCORE:>=10",
                cache_stem="tcga_coad_msi",
            ),
            expression_builders.TreehouseCohort(
                "COAD_MSS",
                "colon adenocarcinoma",
                selection="cbio_sample_clinical:coadread_tcga_pan_can_atlas_2018:MSI_SENSOR_SCORE:<10",
                cache_stem="tcga_coad_mss",
            ),
        ),
    )

    def fake_clinical_map(
        study_id,
        attribute_id,
        *,
        clinical_data_type="PATIENT",
        cache_path=None,
        force_download=False,
    ):
        assert study_id == "coadread_tcga_pan_can_atlas_2018"
        assert attribute_id == "MSI_SENSOR_SCORE"
        assert clinical_data_type == "SAMPLE"
        assert cache_path is not None
        return pd.DataFrame(
            {
                "case_id": ["TCGA-CO-0001", "TCGA-CO-0002", "TCGA-RE-0003"],
                "value": ["12.5", "2.0", "14.0"],
            }
        )

    monkeypatch.setattr(
        expression_builders,
        "treehouse_cbioportal_clinical_attribute_map",
        fake_clinical_map,
    )
    result = expression_builders.build_treehouse_source_matrices(source, cache_dir=tmp_path)

    assert set(result.matrix_paths) == {"COAD_MSI", "COAD_MSS"}
    msi = pd.read_parquet(result.matrix_paths["COAD_MSI"]).set_index("Symbol")
    mss = pd.read_parquet(result.matrix_paths["COAD_MSS"]).set_index("Symbol")
    assert list(msi.columns) == ["Ensembl_Gene_ID", "TCGA-CO-0001-01A"]
    assert list(mss.columns) == ["Ensembl_Gene_ID", "TCGA-CO-0002-01A"]
    np.testing.assert_allclose(msi.loc["TP53", "TCGA-CO-0001-01A"], 2.0)
    np.testing.assert_allclose(mss.loc["MLH1", "TCGA-CO-0002-01A"], 8.0)


def test_build_treehouse_source_matrices_splits_cbio_mutation_cohorts(
    tmp_path,
    monkeypatch,
):
    clinical_path = tmp_path / "clinical.tsv"
    pd.DataFrame(
        [
            {"th_dataset_id": "TCGA-LU-0001-01A", "disease": "lung adenocarcinoma"},
            {"th_dataset_id": "TCGA-LU-0002-01A", "disease": "lung adenocarcinoma"},
            {"th_dataset_id": "TCGA-LU-0003-01A", "disease": "lung adenocarcinoma"},
        ]
    ).to_csv(clinical_path, sep="\t", index=False)
    tpm_path = tmp_path / "treehouse.tsv"
    pd.DataFrame(
        {
            "Gene": ["TP53", "EGFR"],
            "TCGA-LU-0001-01A": np.log2(np.array([2.0, 1.0]) + 1.0),
            "TCGA-LU-0002-01A": np.log2(np.array([4.0, 8.0]) + 1.0),
            "TCGA-LU-0003-01A": np.log2(np.array([6.0, 3.0]) + 1.0),
        }
    ).to_csv(tpm_path, sep="\t", index=False)
    source = expression_builders.TreehouseSource(
        source_id="synthetic-treehouse-luad-mut",
        source_cohort="SYNTHETIC_TREEHOUSE_LUAD_MUT",
        cancer_code=["LUAD_EGFR", "LUAD_STK11"],
        tpm_file=tpm_path.name,
        clinical_file=clinical_path.name,
        source_project="Treehouse (TCGA-LUAD) x cBioPortal mutation calls",
        cohorts=(
            expression_builders.TreehouseCohort(
                "LUAD_EGFR",
                "lung adenocarcinoma",
                selection="cbio_mutation:luad_tcga_pan_can_atlas_2018:EGFR",
                cache_stem="tcga_luad_egfr",
            ),
            expression_builders.TreehouseCohort(
                "LUAD_STK11",
                "lung adenocarcinoma",
                selection="cbio_mutation:luad_tcga_pan_can_atlas_2018:STK11,KEAP1",
                cache_stem="tcga_luad_stk11",
            ),
        ),
    )

    def fake_mutation_case_set(
        study_id,
        gene_symbols,
        *,
        molecular_profile_id=None,
        sample_list_id=None,
        cache_path=None,
        force_download=False,
    ):
        assert study_id == "luad_tcga_pan_can_atlas_2018"
        assert cache_path is not None
        genes = tuple(gene_symbols)
        if genes == ("EGFR",):
            return pd.DataFrame(
                {
                    "case_id": ["TCGA-LU-0001"],
                    "sample_id": ["TCGA-LU-0001-01"],
                    "gene_symbol": ["EGFR"],
                    "entrez_gene_id": [1956],
                }
            )
        assert genes == ("STK11", "KEAP1")
        return pd.DataFrame(
            {
                "case_id": ["TCGA-LU-0002", "TCGA-LU-0003"],
                "sample_id": ["TCGA-LU-0002-01", "TCGA-LU-0003-01"],
                "gene_symbol": ["STK11", "KEAP1"],
                "entrez_gene_id": [6794, 9817],
            }
        )

    monkeypatch.setattr(
        expression_builders,
        "treehouse_cbioportal_mutation_case_set",
        fake_mutation_case_set,
    )
    result = expression_builders.build_treehouse_source_matrices(source, cache_dir=tmp_path)

    assert set(result.matrix_paths) == {"LUAD_EGFR", "LUAD_STK11"}
    egfr = pd.read_parquet(result.matrix_paths["LUAD_EGFR"]).set_index("Symbol")
    stk11 = pd.read_parquet(result.matrix_paths["LUAD_STK11"]).set_index("Symbol")
    assert list(egfr.columns) == ["Ensembl_Gene_ID", "TCGA-LU-0001-01A"]
    assert list(stk11.columns) == [
        "Ensembl_Gene_ID",
        "TCGA-LU-0002-01A",
        "TCGA-LU-0003-01A",
    ]
    np.testing.assert_allclose(egfr.loc["TP53", "TCGA-LU-0001-01A"], 2.0)
    np.testing.assert_allclose(stk11.loc["EGFR", "TCGA-LU-0002-01A"], 8.0)


def test_recount3_gene_sums_to_tpm_length_normalizes_and_collapses_versions():
    gene_sums = pd.DataFrame(
        {"S1": [1000.0, 1000.0], "S2": [0.0, 500.0]},
        index=["ENSG00000000001.5", "ENSG00000000002.3"],
    )
    bp_length = pd.Series({"ENSG00000000001": 1000.0, "ENSG00000000002": 2000.0})

    tpm = expression_builders.recount3_gene_sums_to_tpm(gene_sums, bp_length)

    np.testing.assert_allclose(tpm.sum(axis=0).to_numpy(), [1e6, 1e6])
    np.testing.assert_allclose(tpm["S1"].to_numpy(), [2e6 / 3, 1e6 / 3], rtol=1e-9)
    np.testing.assert_allclose(tpm["S2"].to_numpy(), [0.0, 1e6])

    dup = pd.DataFrame(
        {"S1": [600.0, 400.0]},
        index=["ENSG00000000003.1", "ENSG00000000003.1_PAR_Y"],
    )
    collapsed = expression_builders.recount3_gene_sums_to_tpm(
        dup,
        pd.Series({"ENSG00000000003": 1000.0}),
    )
    assert collapsed.index.tolist() == ["ENSG00000000003"]
    np.testing.assert_allclose(collapsed["S1"].to_numpy(), [1e6])


def test_recount3_parse_attributes_and_aggregate_runs_to_samples():
    attrs = expression_builders.parse_recount3_sample_attributes(
        "origin;;pancreas|type;;liver metastasis|n;;1"
    )
    assert attrs == {"origin": "pancreas", "type": "liver metastasis", "n": "1"}
    assert expression_builders.parse_recount3_sample_attributes("") == {}

    gene_sums = pd.DataFrame(
        {"R1": [10.0, 1.0], "R2": [20.0, 3.0], "R3": [5.0, 5.0], "R4": [99.0, 99.0]},
        index=["g1", "g2"],
    )
    meta = pd.DataFrame(
        {
            "external_id": ["R1", "R2", "R3", "R4"],
            "sample_acc": ["A", "A", "B", "C"],
        }
    )

    sample_gs, sample_meta = expression_builders.aggregate_recount3_runs_to_samples(
        gene_sums,
        meta,
        keep_runs={"R1", "R2", "R3"},
    )

    assert set(sample_gs.columns) == {"A", "B"}
    np.testing.assert_allclose(sample_gs.loc["g1", "A"], 30.0)
    np.testing.assert_allclose(sample_gs.loc["g2", "A"], 4.0)
    np.testing.assert_allclose(sample_gs.loc["g1", "B"], 5.0)
    assert list(sample_meta.index) == list(sample_gs.columns)


def test_build_recount3_source_matrices_writes_canonical_artifacts(tmp_path, monkeypatch):
    annotation = pd.DataFrame(
        {
            "bp_length": [1000.0, 2000.0],
            "Symbol": ["TP53", "EGFR"],
        },
        index=["ENSG00000141510", "ENSG00000146648"],
    )
    gene_sums = pd.DataFrame(
        {
            "R1": [100.0, 50.0],
            "R2": [50.0, 25.0],
            "R3": [0.0, 80.0],
            "R4": [999.0, 999.0],
        },
        index=["ENSG00000141510.1", "ENSG00000146648.2"],
    )
    metadata = pd.DataFrame(
        {
            "external_id": ["R1", "R2", "R3", "R4"],
            "sample_acc": ["SAMPLE_A", "SAMPLE_A", "SAMPLE_B", "CONTROL"],
            "sample_attributes": [
                "code;;CODE_A",
                "code;;CODE_A",
                "code;;CODE_B",
                "code;;CONTROL",
            ],
            "sample_title": ["A1", "A2", "B1", "C1"],
        }
    )
    monkeypatch.setattr(
        expression_builders,
        "fetch_recount3_gene_annotation",
        lambda _cache: annotation,
    )
    monkeypatch.setattr(
        expression_builders,
        "fetch_recount3_gene_sums",
        lambda _srp, _cache: gene_sums,
    )
    monkeypatch.setattr(
        expression_builders,
        "fetch_recount3_sample_metadata",
        lambda _srp, _cache: metadata,
    )
    source = expression_builders.Recount3Source(
        source_id="synthetic-recount3",
        srp="SRP000000",
        source_cohort="SYNTHETIC_RECOUNT3",
        cancer_code=["CODE_A", "CODE_B"],
        sample_to_cancer_code=lambda attrs, _title: (
            attrs.get("code") if attrs.get("code") != "CONTROL" else None
        ),
        expected_n={"CODE_A": 1, "CODE_B": 1},
    )
    out_dir = tmp_path / "derived"
    out_dir.mkdir()
    stale_matrix = out_dir / "STALE_per_sample_tpm.parquet"
    stale_qc = out_dir / "STALE_sample_qc.csv"
    stale_matrix.write_text("stale matrix")
    stale_qc.write_text("stale qc")

    result = expression_builders.build_recount3_source_matrices(source, cache_dir=tmp_path)

    assert set(result.matrix_paths) == {"CODE_A", "CODE_B"}
    assert not stale_matrix.exists()
    assert not stale_qc.exists()
    assert result.sidecar_paths["mapping_audit"].exists()
    assert result.sidecar_paths["parse_diagnostics"].exists()
    code_a = pd.read_parquet(result.matrix_paths["CODE_A"])
    code_b = pd.read_parquet(result.matrix_paths["CODE_B"])
    assert list(code_a.columns) == ["Ensembl_Gene_ID", "Symbol", "SAMPLE_A"]
    assert list(code_b.columns) == ["Ensembl_Gene_ID", "Symbol", "SAMPLE_B"]
    assert set(code_a["Ensembl_Gene_ID"]) == {"ENSG00000141510", "ENSG00000146648"}
    assert np.isclose(code_a["SAMPLE_A"].sum(), 1_000_000.0)
    assert np.isclose(code_b["SAMPLE_B"].sum(), 1_000_000.0)
    assert set(result.sample_qc["cancer_code"]) == {"CODE_A", "CODE_B"}


def test_build_geo_matrix_script_uses_registry_config(tmp_path, capsys):
    source_path = tmp_path / "source.tsv"
    pd.DataFrame(
        {
            "gene": ["TP53", "EGFR"],
            "tumor_a": ["1", "3"],
            "normal_a": ["10", "10"],
        }
    ).to_csv(source_path, sep="\t", index=False)
    registry_path = tmp_path / "expression_sources.yaml"
    registry_path.write_text(
        """
sources:
  - id: synthetic-geo
    source_type: geo-matrix
    cancer_codes: [CODE_A]
    file_url: https://example.org/source.tsv.gz
    file_name: source.tsv.gz
    unit: TPM
    gene_id_col: gene
    source_cohort: SYNTHETIC_GEO
    source_project: GEO
    sample_filter:
      include_match: "^tumor_"
""".lstrip()
    )
    mod = _load_script("build_geo_matrix")

    status = mod.main(
        [
            "--source-id",
            "synthetic-geo",
            "--registry",
            str(registry_path),
            "--cache-dir",
            str(tmp_path / "cache"),
            "--output-dir",
            str(tmp_path / "out"),
            "--source-path",
            str(source_path),
        ]
    )

    assert status == 0
    out = pd.read_parquet(tmp_path / "out" / "CODE_A_per_sample_tpm.parquet")
    assert list(out.columns) == ["Ensembl_Gene_ID", "Symbol", "tumor_a"]
    assert set(out["Symbol"]) == {"TP53", "EGFR"}
    assert np.isclose(out["tumor_a"].sum(), 1_000_000.0)
    stdout = capsys.readouterr().out
    assert '"sample_counts": {' in stdout
    assert '"CODE_A": 1' in stdout


def test_build_recount3_script_uses_registry_config(tmp_path, monkeypatch, capsys):
    mod = _load_script("build_recount3_source")
    source = expression_builders.Recount3Source(
        source_id="synthetic-recount3",
        srp="SRP000000",
        source_cohort="SYNTHETIC_RECOUNT3",
        cancer_code="CODE_A",
    )

    def _fake_build(source_obj, *, cache_dir, output_dir=None, **_kwargs):
        assert source_obj is source
        assert Path(cache_dir) == tmp_path / "cache"
        assert output_dir is None
        matrix = pd.DataFrame(
            {
                "Ensembl_Gene_ID": ["ENSG00000141510"],
                "Symbol": ["TP53"],
                "S1": [1_000_000.0],
            }
        )
        return expression_builders.SourceMatrixBuildResult(
            source=source,
            matrices={"CODE_A": matrix},
            matrix_paths={"CODE_A": tmp_path / "CODE_A_per_sample_tpm.parquet"},
            summary_rows=pd.DataFrame(
                columns=list(expression_builders.REFERENCE_EXPRESSION_COLUMNS)
            ),
            mapping_audit=pd.DataFrame(),
            parse_diagnostics=pd.DataFrame(),
            sample_qc=pd.DataFrame(),
            sidecar_paths={"mapping_audit": tmp_path / "mapping_audit.csv"},
        )

    monkeypatch.setattr(
        mod,
        "recount3_source_from_registry",
        lambda source_id, registry_path=None: source,
    )
    monkeypatch.setattr(mod, "build_recount3_source_matrices", _fake_build)

    assert mod.main(["synthetic-recount3", "--cache-dir", str(tmp_path / "cache")]) == 0
    stdout = capsys.readouterr().out
    assert '"source_id": "synthetic-recount3"' in stdout
    assert '"CODE_A": 1' in stdout


def test_build_treehouse_script_uses_registry_config(tmp_path, monkeypatch, capsys):
    mod = _load_script("build_treehouse_source")
    source = expression_builders.TreehouseSource(
        source_id="synthetic-treehouse",
        source_cohort="SYNTHETIC_TREEHOUSE",
        cancer_code="CODE_A",
        tpm_file="treehouse.tsv",
        clinical_file="clinical.tsv",
        cohorts=(expression_builders.TreehouseCohort("CODE_A", "synthetic tumor"),),
    )

    def _fake_build(source_obj, *, cache_dir, output_dir=None, **_kwargs):
        assert source_obj is source
        assert Path(cache_dir) == tmp_path / "cache"
        assert output_dir is None
        matrix = pd.DataFrame(
            {
                "Ensembl_Gene_ID": ["ENSG00000141510"],
                "Symbol": ["TP53"],
                "S1": [10.0],
            }
        )
        return expression_builders.SourceMatrixBuildResult(
            source=source,
            matrices={"CODE_A": matrix},
            matrix_paths={"CODE_A": tmp_path / "CODE_A_per_sample_tpm.parquet"},
            summary_rows=pd.DataFrame(
                columns=list(expression_builders.REFERENCE_EXPRESSION_COLUMNS)
            ),
            mapping_audit=pd.DataFrame(),
            parse_diagnostics=pd.DataFrame(),
            sample_qc=pd.DataFrame(),
            sidecar_paths={"mapping_audit": tmp_path / "mapping_audit.csv"},
        )

    monkeypatch.setattr(
        mod,
        "treehouse_source_from_registry",
        lambda source_id, registry_path=None: source,
    )
    monkeypatch.setattr(mod, "build_treehouse_source_matrices", _fake_build)

    assert mod.main(["synthetic-treehouse", "--cache-dir", str(tmp_path / "cache")]) == 0
    stdout = capsys.readouterr().out
    assert '"source_id": "synthetic-treehouse"' in stdout
    assert '"CODE_A": 1' in stdout


def test_build_geo_matrix_script_requires_gene_lengths_for_raw_counts(tmp_path):
    registry_path = tmp_path / "expression_sources.yaml"
    registry_path.write_text(
        """
sources:
  - id: synthetic-counts
    source_type: geo-matrix
    cancer_codes: [CODE_A]
    file_url: https://example.org/source.tsv.gz
    file_name: source.tsv.gz
    unit: raw_counts
    gene_id_col: gene
    source_cohort: SYNTHETIC_COUNTS
""".lstrip()
    )
    mod = _load_script("build_geo_matrix")

    with pytest.raises(SystemExit, match="raw_counts"):
        mod.main(
            [
                "--source-id",
                "synthetic-counts",
                "--registry",
                str(registry_path),
                "--cache-dir",
                str(tmp_path / "cache"),
            ]
        )


# ---------- cohort_percentile_vectors ----------


def test_percentile_vectors_schema_matches_reader():
    # 26 dense breakpoints p0..p100, plus the two id columns — the exact schema
    # expression.cohort_gene_percentiles reads back.
    genes = ["A", "B"]
    vals = np.array([[1.0, 2.0, 3.0, 4.0], [10.0, 20.0, 30.0, 40.0]])
    out = expression_builders.cohort_percentile_vectors(_matrix(genes, list("wxyz"), vals))
    bp_cols = [c for c in out.columns if c not in ("Ensembl_Gene_ID", "Symbol")]
    assert len(bp_cols) == 26
    assert bp_cols[0] == "p0" and bp_cols[-1] == "p100"
    assert "p50" in bp_cols
    assert len(out) == 2


def test_percentile_vectors_log1p_roundtrip():
    # Stored log1p; expm1 of p0/p50/p100 recovers min/median/max of the gene's
    # across-sample distribution (the reader's as_tpm=True path).
    vals = np.array([[0.0, 10.0, 100.0, 1000.0]])
    out = expression_builders.cohort_percentile_vectors(_matrix(["A"], list("wxyz"), vals))
    row = out.iloc[0]
    assert np.isclose(np.expm1(np.float32(row["p0"])), 0.0, atol=1e-1)
    assert np.isclose(np.expm1(np.float32(row["p100"])), 1000.0, rtol=2e-2)
    # median of [0,10,100,1000] in log1p space, restored
    expected_med = np.median(np.log1p(vals[0]))
    assert np.isclose(np.float32(row["p50"]), expected_med, rtol=2e-2)


def test_percentile_vectors_ignores_nan():
    # A gene unmeasured in some samples: NaN cells are dropped, not treated as 0.
    vals = np.array([[np.nan, 10.0, 10.0, 10.0]])
    out = expression_builders.cohort_percentile_vectors(_matrix(["A"], list("wxyz"), vals))
    assert np.isclose(np.expm1(np.float32(out.iloc[0]["p50"])), 10.0, rtol=2e-2)


def test_percentile_vectors_requires_samples():
    with pytest.raises(ValueError):
        expression_builders.cohort_percentile_vectors(_matrix(["A"], [], np.empty((1, 0))))


# ---------- cohort_medoids ----------


def test_medoids_returns_base_plus_k():
    rng = np.arange(50.0).reshape(5, 10)  # 5 genes × 10 samples
    out = expression_builders.cohort_medoids(
        _matrix([f"g{i}" for i in range(5)], [f"s{j}" for j in range(10)], rng), k=3
    )
    rep_cols = [c for c in out.columns if c not in ("Ensembl_Gene_ID", "Symbol")]
    assert len(rep_cols) == 3
    assert list(out["Ensembl_Gene_ID"]) == [f"g{i}" for i in range(5)]


def test_medoids_small_cohort_keeps_all():
    vals = np.array([[1.0, 2.0]])
    out = expression_builders.cohort_medoids(_matrix(["A"], ["s1", "s2"], vals), k=5)
    assert [c for c in out.columns if c not in ("Ensembl_Gene_ID", "Symbol")] == ["s1", "s2"]


def test_medoids_small_cohort_uses_stable_sample_id_order():
    vals = np.array([[2.0, 1.0]])
    out = expression_builders.cohort_medoids(_matrix(["A"], ["s2", "s1"], vals), k=5)
    assert [c for c in out.columns if c not in ("Ensembl_Gene_ID", "Symbol")] == ["s1", "s2"]


def test_medoids_central_first_then_outlier():
    # 4 near-identical "typical" samples + 1 far outlier. The medoid (pick 1)
    # must come from the dense cluster; the farthest pick (2) must be the outlier.
    genes = [f"g{i}" for i in range(6)]
    typical = np.tile(np.array([[1.0], [2.0], [3.0], [4.0], [5.0], [6.0]]), (1, 4))
    typical = typical + np.array([[0, 0.01, -0.01, 0.0]] * 6)  # tiny jitter
    outlier = np.array([[100.0], [100.0], [100.0], [100.0], [100.0], [100.0]])
    vals = np.hstack([typical, outlier])
    samples = ["t1", "t2", "t3", "t4", "outlier"]
    out = expression_builders.cohort_medoids(_matrix(genes, samples, vals), k=2)
    picks = [c for c in out.columns if c not in ("Ensembl_Gene_ID", "Symbol")]
    assert picks[0] != "outlier"  # central medoid from the dense cluster
    assert picks[1] == "outlier"  # farthest-first grabs the outlier


def test_medoids_preserve_original_tpm():
    # Distance uses log1p internally, but stored values are the original TPM.
    vals = np.array([[7.0, 8.0, 9.0]])
    out = expression_builders.cohort_medoids(_matrix(["A"], ["s1", "s2", "s3"], vals), k=3)
    kept = out[[c for c in out.columns if c not in ("Ensembl_Gene_ID", "Symbol")]].to_numpy()
    assert set(kept.ravel()) == {7.0, 8.0, 9.0}


def test_medoids_can_select_on_biological_view_but_return_full_values():
    full = _matrix(
        ["BIO", "TECH"],
        ["s3", "s1", "s2"],
        np.array(
            [
                [10.0, 0.0, 5.0],
                [1000.0, 1000.0, 0.0],
            ]
        ),
    )
    biological = full[full["Ensembl_Gene_ID"] == "BIO"].reset_index(drop=True)

    out = expression_builders.cohort_medoids(
        full,
        sample_cols=["s3", "s1", "s2"],
        k=1,
        selection_df=biological,
    )

    rep_cols = [c for c in out.columns if c not in ("Ensembl_Gene_ID", "Symbol")]
    assert rep_cols == ["s2"]
    assert list(out["Ensembl_Gene_ID"]) == ["BIO", "TECH"]
    assert out.set_index("Ensembl_Gene_ID").loc["TECH", "s2"] == 0.0


def test_medoids_deterministic():
    rng = (np.arange(60.0) * 1.7 % 11).reshape(6, 10)
    df = _matrix([f"g{i}" for i in range(6)], [f"s{j}" for j in range(10)], rng)
    a = expression_builders.cohort_medoids(df, k=4)
    b = expression_builders.cohort_medoids(df, k=4)
    assert list(a.columns) == list(b.columns)


# ---------- generator scripts (end-to-end on a synthetic per-sample dir) ----------


def _write_cohort(tmp_path, code, genes, samples, values):
    df = _matrix(genes, samples, values)
    path = tmp_path / f"{code}.parquet"
    df.to_parquet(path, index=False)
    return path


def test_percentiles_generator_writes_readable_shards(tmp_path):
    gen = _load_script("generate_cohort_percentiles")
    inp = tmp_path / "in"
    inp.mkdir()
    _write_cohort(
        inp,
        "COHORT_A",
        ["A", "B"],
        list("wxyz"),
        np.array([[1.0, 2.0, 3.0, 4.0], [10.0, 20.0, 30.0, 40.0]]),
    )
    out = tmp_path / "out"
    gen.build(inp, drop_genes=set(), out_dir=out)
    shard = pd.read_parquet(out / "COHORT_A.parquet")
    assert len(shard) == 2
    bp_cols = [c for c in shard.columns if c not in ("Ensembl_Gene_ID", "Symbol")]
    assert len(bp_cols) == 26
    # p100 of gene A (max 4.0), stored log1p
    a = shard[shard["Ensembl_Gene_ID"] == "A"].iloc[0]
    assert np.isclose(np.expm1(np.float32(a["p100"])), 4.0, rtol=2e-2)


def test_percentiles_generator_drops_genes(tmp_path):
    gen = _load_script("generate_cohort_percentiles")
    inp = tmp_path / "in"
    inp.mkdir()
    _write_cohort(
        inp,
        "COHORT_A",
        ["KEEP", "DROPME"],
        ["s1", "s2"],
        np.array([[1.0, 2.0], [9.0, 9.0]]),
    )
    out = tmp_path / "out"
    gen.build(inp, drop_genes={"DROPME"}, out_dir=out)
    shard = pd.read_parquet(out / "COHORT_A.parquet")
    assert list(shard["Ensembl_Gene_ID"]) == ["KEEP"]


def test_representatives_generator_writes_shards_and_provenance(tmp_path):
    gen = _load_script("generate_representatives")
    inp = tmp_path / "in"
    inp.mkdir()
    # 6 samples so k=3 actually selects a subset
    vals = np.array([[float(g * 10 + s) for s in range(6)] for g in range(4)])
    _write_cohort(inp, "COHORT_A", [f"g{i}" for i in range(4)], [f"s{j}" for j in range(6)], vals)
    out = tmp_path / "out"
    gen.build(inp, k=3, out_dir=out)

    shard = pd.read_parquet(out / "COHORT_A.parquet")
    rep_cols = [c for c in shard.columns if c not in ("Ensembl_Gene_ID", "Symbol")]
    assert rep_cols == ["COHORT_A__rep1", "COHORT_A__rep2", "COHORT_A__rep3"]

    prov = pd.read_csv(out / "_provenance.csv")
    assert set(prov["representative_id"]) == set(rep_cols)
    assert (prov["n_cohort_samples"] == 6).all()
    # The reader merges on these exact columns — all must be present (source_project
    # is best-effort and empty for an unregistered synthetic code, but the column
    # must exist so consumers don't KeyError).
    for col in ("representative_id", "source_cohort", "source_project", "n_cohort_samples"):
        assert col in prov.columns
    # Unregistered code -> source_cohort falls back to the code itself.
    assert (prov["source_cohort"] == "COHORT_A").all()


def test_representatives_generator_selects_on_biological_view(tmp_path, monkeypatch):
    gen = _load_script("generate_representatives")
    inp = tmp_path / "in"
    inp.mkdir()
    _write_cohort(
        inp,
        "COHORT_A",
        ["BIO", "TECH"],
        ["sample_a", "sample_b"],
        np.array([[5.0, 10.0], [1000.0, 0.0]]),
    )
    monkeypatch.setattr(gen, "clean_tpm_censored_gene_ids", lambda: {"TECH"})
    seen = {}

    def fake_medoids(df, sample_cols, *, k, selection_df):
        seen["value_ids"] = df["Ensembl_Gene_ID"].tolist()
        seen["selection_ids"] = selection_df["Ensembl_Gene_ID"].tolist()
        seen["sample_cols"] = list(sample_cols)
        return df[["Ensembl_Gene_ID", "Symbol", sample_cols[0]]].copy()

    monkeypatch.setattr(gen, "cohort_medoids", fake_medoids)

    gen.build(inp, k=1, out_dir=tmp_path / "out")

    assert seen == {
        "value_ids": ["BIO", "TECH"],
        "selection_ids": ["BIO"],
        "sample_cols": ["sample_a", "sample_b"],
    }


def _write_rebuild_inputs(tmp_path):
    cache = tmp_path / "cache"
    ref = tmp_path / "ref"
    source_dir = cache / "TEST_SOURCE" / "derived"
    source_dir.mkdir(parents=True)
    ref.mkdir()
    _matrix(
        ["ENSG000001", "ENSG000002", "ENSG000003"],
        ["pass_sample", "warn_sample", "fail_sample"],
        np.array(
            [
                [10.0, 20.0, 30.0],
                [40.0, 50.0, 60.0],
                [70.0, 80.0, 90.0],
            ]
        ),
    ).to_parquet(source_dir / "X_per_sample_tpm.parquet", index=False)
    _write_cohort(ref, "X", ["ENSG000001"], ["reference"], np.array([[1.0]]))
    return cache, ref


def _patch_rebuild_registry(monkeypatch, gen):
    monkeypatch.setattr(
        gen,
        "source_registry",
        lambda: pd.DataFrame({"cancer_code": ["X"], "source_cohort": ["TEST_SOURCE"]}),
    )
    monkeypatch.setattr(gen, "cohort_source_version", lambda code: "test-source-version")
    monkeypatch.setattr(
        gen,
        "sample_expression_qc_from_matrix",
        lambda raw, cancer_type: pd.DataFrame(
            {
                "cancer_code": [cancer_type, cancer_type, cancer_type],
                "sample_id": ["pass_sample", "warn_sample", "fail_sample"],
                "sample_qc_status": ["pass", "warn", "fail"],
                "sample_qc_reasons": [
                    "",
                    "nonlinear_or_proxy_expression_scale",
                    "low_detected_genes",
                ],
            }
        ),
    )


def test_rebuild_expression_artifacts_defaults_to_qc_passing_samples(tmp_path, monkeypatch):
    gen = _load_script("rebuild_expression_artifacts")
    cache, ref = _write_rebuild_inputs(tmp_path)
    _patch_rebuild_registry(monkeypatch, gen)
    out = tmp_path / "out"

    gen.rebuild(cache, ref, out, limit=None, validate=False, sample_qc="pass")

    clean = pd.read_parquet(out / "clean" / "X.parquet")
    assert [c for c in clean.columns if c not in ("Ensembl_Gene_ID", "Symbol")] == ["pass_sample"]

    reps = pd.read_parquet(out / "cancer-reference-expression-representatives" / "X.parquet")
    assert [c for c in reps.columns if c not in ("Ensembl_Gene_ID", "Symbol")] == ["X__rep1"]

    prov = pd.read_csv(out / "cancer-reference-expression-representatives" / "_provenance.csv")
    assert prov.loc[0, "n_source_samples"] == 3
    assert prov.loc[0, "n_cohort_samples"] == 1
    assert prov.loc[0, "sample_qc"] == "pass"
    assert prov.loc[0, "sample_qc_policy_version"] == "sample_expression_qc_v2"
    assert prov.loc[0, "n_qc_pass"] == 1
    assert prov.loc[0, "n_qc_warn"] == 1
    assert prov.loc[0, "n_qc_fail"] == 1

    qc = pd.read_csv(out / "source-matrix-sample-qc.csv")
    assert list(qc["sample_id"]) == ["pass_sample", "warn_sample", "fail_sample"]

    build_meta = pd.read_csv(out / "expression-artifact-build-metadata.csv")
    assert build_meta.loc[0, "n_source_samples"] == 3
    assert build_meta.loc[0, "n_cohort_samples"] == 1
    assert build_meta.loc[0, "sample_qc_policy_version"] == "sample_expression_qc_v2"

    metadata = json.loads((out / "expression-artifact-build-metadata.json").read_text())
    assert metadata["sample_qc"] == "pass"
    assert metadata["sample_qc_manifest"] == "source-matrix-sample-qc.csv"
    assert metadata["n_source_samples"] == 3
    assert metadata["n_cohort_samples"] == 1


def test_rebuild_expression_artifacts_keeps_warn_proxy_source_when_pass_empty(
    tmp_path, monkeypatch
):
    gen = _load_script("rebuild_expression_artifacts")
    cache, ref = _write_rebuild_inputs(tmp_path)
    monkeypatch.setattr(
        gen,
        "source_registry",
        lambda: pd.DataFrame({"cancer_code": ["X"], "source_cohort": ["TEST_SOURCE"]}),
    )
    monkeypatch.setattr(gen, "cohort_source_version", lambda code: "test-source-version")
    monkeypatch.setattr(
        gen,
        "sample_expression_qc_from_matrix",
        lambda raw, cancer_type: pd.DataFrame(
            {
                "cancer_code": [cancer_type, cancer_type, cancer_type],
                "sample_id": ["pass_sample", "warn_sample", "fail_sample"],
                "sample_qc_status": ["warn", "warn", "warn"],
                "sample_qc_reasons": [
                    "nonlinear_or_proxy_expression_scale",
                    "nonlinear_or_proxy_expression_scale",
                    "nonlinear_or_proxy_expression_scale",
                ],
                "source_scale_class": [
                    "microarray_tpm_proxy",
                    "microarray_tpm_proxy",
                    "microarray_tpm_proxy",
                ],
            }
        ),
    )
    out = tmp_path / "out"

    gen.rebuild(cache, ref, out, limit=None, validate=False, sample_qc="pass")

    clean = pd.read_parquet(out / "clean" / "X.parquet")
    assert [c for c in clean.columns if c not in ("Ensembl_Gene_ID", "Symbol")] == [
        "pass_sample",
        "warn_sample",
        "fail_sample",
    ]
    build_meta = pd.read_csv(out / "expression-artifact-build-metadata.csv")
    assert build_meta.loc[0, "sample_qc"] == "pass"
    assert build_meta.loc[0, "sample_qc_effective"] == "pass_or_warn"
    assert build_meta.loc[0, "sample_qc_fallback_reason"] == "no_pass_samples_tpm_proxy_source"
    metadata = json.loads((out / "expression-artifact-build-metadata.json").read_text())
    assert metadata["sample_qc_fallbacks"] == 1


def test_rebuild_expression_artifacts_keeps_concentration_only_source_when_pass_empty(
    tmp_path, monkeypatch
):
    gen = _load_script("rebuild_expression_artifacts")
    cache, ref = _write_rebuild_inputs(tmp_path)
    monkeypatch.setattr(
        gen,
        "source_registry",
        lambda: pd.DataFrame({"cancer_code": ["X"], "source_cohort": ["TEST_SOURCE"]}),
    )
    monkeypatch.setattr(gen, "cohort_source_version", lambda code: "test-source-version")
    monkeypatch.setattr(
        gen,
        "sample_expression_qc_from_matrix",
        lambda raw, cancer_type: pd.DataFrame(
            {
                "cancer_code": [cancer_type, cancer_type, cancer_type],
                "sample_id": ["pass_sample", "warn_sample", "fail_sample"],
                "sample_qc_status": ["fail", "fail", "fail"],
                "sample_qc_reasons": [
                    "high_top10_gene_fraction",
                    "high_top_gene_fraction;high_top10_gene_fraction",
                    "high_top10_gene_fraction",
                ],
            }
        ),
    )
    out = tmp_path / "out"

    gen.rebuild(cache, ref, out, limit=None, validate=False, sample_qc="pass")

    clean = pd.read_parquet(out / "clean" / "X.parquet")
    assert [c for c in clean.columns if c not in ("Ensembl_Gene_ID", "Symbol")] == [
        "pass_sample",
        "warn_sample",
        "fail_sample",
    ]
    build_meta = pd.read_csv(out / "expression-artifact-build-metadata.csv")
    assert build_meta.loc[0, "sample_qc"] == "pass"
    assert build_meta.loc[0, "sample_qc_effective"] == "all"
    assert (
        build_meta.loc[0, "sample_qc_fallback_reason"]
        == "no_pass_samples_high_concentration_source"
    )


def test_rebuild_expression_artifacts_clips_negative_source_values():
    gen = _load_script("rebuild_expression_artifacts")
    df = _matrix(["G1", "G2"], ["s1", "s2"], np.array([[-2.0, 3.0], [4.0, -0.5]]))

    out, n_negative = gen._clip_negative_expression(df, ["s1", "s2"])

    assert n_negative == 2
    assert out[["s1", "s2"]].to_numpy().min() == 0.0
    assert df[["s1", "s2"]].to_numpy().min() < 0.0


def test_rebuild_expression_artifacts_disambiguates_source_by_registry_sample_count(
    tmp_path, monkeypatch
):
    gen = _load_script("rebuild_expression_artifacts")
    cache = tmp_path / "cache"
    ref = tmp_path / "ref"
    ref.mkdir()
    _write_cohort(ref, "X", ["G1", "G2"], ["a", "b"], np.array([[1.0, 2.0], [3.0, 4.0]]))
    source_a = cache / "source-a" / "derived"
    source_b = cache / "source-b" / "derived"
    source_a.mkdir(parents=True)
    source_b.mkdir(parents=True)
    _matrix(["G1", "G2"], ["s1"], np.array([[1.0], [2.0]])).to_parquet(
        source_a / "X_per_sample_tpm.parquet", index=False
    )
    _matrix(["G1", "G2"], ["s1", "s2"], np.array([[1.0, 2.0], [2.0, 3.0]])).to_parquet(
        source_b / "X_per_sample_tpm.parquet", index=False
    )
    monkeypatch.setattr(
        gen,
        "source_registry",
        lambda: pd.DataFrame(
            {"cancer_code": ["X"], "source_cohort": ["AGGREGATE_LABEL"], "n_samples": [2]}
        ),
    )
    monkeypatch.setattr(gen, "cohort_source_version", lambda code: "test-source-version")
    monkeypatch.setattr(
        gen,
        "sample_expression_qc_from_matrix",
        lambda raw, cancer_type: pd.DataFrame(
            {
                "cancer_code": [cancer_type] * (len(raw.columns) - 2),
                "sample_id": [c for c in raw.columns if c not in ("Ensembl_Gene_ID", "Symbol")],
                "sample_qc_status": ["pass"] * (len(raw.columns) - 2),
                "sample_qc_reasons": [""] * (len(raw.columns) - 2),
            }
        ),
    )

    gen.rebuild(cache, ref, tmp_path / "out", limit=None, validate=False, sample_qc="pass")

    build_meta = pd.read_csv(tmp_path / "out" / "expression-artifact-build-metadata.csv")
    assert build_meta.loc[0, "n_source_samples"] == 2
    assert build_meta.loc[0, "source_matrix_path"].endswith(
        "source-b/derived/X_per_sample_tpm.parquet"
    )


def test_rebuild_expression_artifacts_prefers_aggregate_source_directory(tmp_path):
    gen = _load_script("rebuild_expression_artifacts")
    preferred = tmp_path / "geo-heme" / "derived" / "X_per_sample_tpm.parquet"
    duplicate = tmp_path / "gse100026-cml" / "derived" / "X_per_sample_tpm.parquet"
    preferred.parent.mkdir(parents=True)
    duplicate.parent.mkdir(parents=True)
    _matrix(["G1"], ["s1"], np.array([[1.0]])).to_parquet(preferred, index=False)
    _matrix(["G1"], ["s1"], np.array([[1.0]])).to_parquet(duplicate, index=False)

    selected = gen._select_source(
        "X",
        [("gse100026-cml", duplicate), ("geo-heme", preferred)],
        {"X": "GEO_HEME_2022"},
        {"X": 1},
    )

    assert selected == preferred


def test_rebuild_expression_artifacts_selects_representatives_on_biological_view(
    tmp_path, monkeypatch
):
    gen = _load_script("rebuild_expression_artifacts")
    cache = tmp_path / "cache"
    ref = tmp_path / "ref"
    source_dir = cache / "TEST_SOURCE" / "derived"
    source_dir.mkdir(parents=True)
    ref.mkdir()
    _matrix(
        ["BIO", "TECH"],
        ["sample_a", "sample_b"],
        np.array([[5.0, 10.0], [1000.0, 0.0]]),
    ).to_parquet(source_dir / "X_per_sample_tpm.parquet", index=False)
    _write_cohort(ref, "X", ["BIO"], ["reference"], np.array([[1.0]]))

    monkeypatch.setattr(
        gen,
        "source_registry",
        lambda: pd.DataFrame({"cancer_code": ["X"], "source_cohort": ["TEST_SOURCE"]}),
    )
    monkeypatch.setattr(gen, "cohort_source_version", lambda code: "test-source-version")
    monkeypatch.setattr(
        gen,
        "sample_expression_qc_from_matrix",
        lambda raw, cancer_type: pd.DataFrame(
            {
                "cancer_code": [cancer_type, cancer_type],
                "sample_id": ["sample_a", "sample_b"],
                "sample_qc_status": ["pass", "pass"],
                "sample_qc_reasons": ["", ""],
            }
        ),
    )
    monkeypatch.setattr(gen, "clean_tpm", lambda values, gene_table: values.copy())
    monkeypatch.setattr(gen, "clean_tpm_censored_gene_ids", lambda: {"TECH"})

    seen = {}

    def fake_medoids(df, sample_cols, *, k, selection_df):
        seen["value_ids"] = df["Ensembl_Gene_ID"].tolist()
        seen["selection_ids"] = selection_df["Ensembl_Gene_ID"].tolist()
        seen["sample_cols"] = list(sample_cols)
        return df[["Ensembl_Gene_ID", "Symbol", sample_cols[0]]].copy()

    monkeypatch.setattr(gen, "cohort_medoids", fake_medoids)

    gen.rebuild(cache, ref, tmp_path / "out", limit=None, validate=False, sample_qc="pass")

    assert seen == {
        "value_ids": ["BIO", "TECH"],
        "selection_ids": ["BIO"],
        "sample_cols": ["sample_a", "sample_b"],
    }


def test_rebuild_expression_artifacts_keeps_all_samples_when_requested(tmp_path, monkeypatch):
    gen = _load_script("rebuild_expression_artifacts")
    cache, ref = _write_rebuild_inputs(tmp_path)
    _patch_rebuild_registry(monkeypatch, gen)
    out = tmp_path / "out"

    gen.rebuild(cache, ref, out, limit=None, validate=False, sample_qc="all")

    clean = pd.read_parquet(out / "clean" / "X.parquet")
    assert [c for c in clean.columns if c not in ("Ensembl_Gene_ID", "Symbol")] == [
        "pass_sample",
        "warn_sample",
        "fail_sample",
    ]
    build_meta = pd.read_csv(out / "expression-artifact-build-metadata.csv")
    assert build_meta.loc[0, "sample_qc"] == "all"
    assert build_meta.loc[0, "n_source_samples"] == 3
    assert build_meta.loc[0, "n_cohort_samples"] == 3


# ---------- real-data parity (skipped without the maintainer's matrix cache) ----------


@pytest.mark.skipif(not _PARITY_READY, reason="per-sample matrix cache / pirlygenes ref absent")
def test_percentiles_reproduce_pirlygenes_reference():
    # End-to-end on REAL data: raw per-sample matrix -> clean_tpm -> percentile
    # vectors must reproduce pirlygenes' shipped percentile artifact for the same
    # cohort. Proves the generator + oncoref's clean_tpm port are faithful.
    from oncoref import normalization as nz

    raw = pd.read_parquet(_ACC_MATRIX[0])
    samples = [c for c in raw.columns if c not in ("Ensembl_Gene_ID", "Symbol")]
    gene_table = raw[["Symbol", "Ensembl_Gene_ID"]]
    clean = nz.clean_tpm(raw[samples], gene_table=gene_table)
    clean_df = pd.concat([gene_table, clean], axis=1)

    mine = expression_builders.cohort_percentile_vectors(clean_df, samples).set_index(
        "Ensembl_Gene_ID"
    )
    ref = pd.read_parquet(_ACC_REF).set_index("Ensembl_Gene_ID")
    # Column schema is identical.
    assert [c for c in mine.columns if c != "Symbol"] == [c for c in ref.columns if c != "Symbol"]

    common = mine.index.intersection(ref.index)
    assert len(common) > 10_000
    # The deterministic mid/upper percentiles match (expm1 back to TPM); tiny tail
    # deviation at p99 is float16 rounding, so correlation must be essentially 1.
    for col in ("p50", "p95"):
        a = np.expm1(mine.loc[common, col].astype("float32"))
        b = np.expm1(ref.loc[common, col].astype("float32"))
        mask = (a > 0) | (b > 0)
        assert np.corrcoef(a[mask], b[mask])[0, 1] > 0.999

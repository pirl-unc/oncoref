"""HPA cancer denominators, cohort identity, crosswalk and source integrity."""

import hashlib
import io
import zipfile

import numpy as np
import pandas as pd
import pytest
from build_hpa_cancer_reference import complete_gene_chunks

from oncoref import hpa_cancer, reference_data


def ihc_row(cancer="lung cancer", **counts):
    return {
        "Gene": "ENSG1",
        "Gene name": "GENE",
        "Cancer": cancer,
        "High": 1,
        "Medium": 2,
        "Low": 3,
        "Not detected": 4,
        **counts,
    }


def rna_rows(values, cancer="Lung Adenocarcinoma (TCGA)"):
    return pd.DataFrame(
        {
            "Gene": "ENSG1",
            "Sample": [f"S{i}" for i in range(len(values))],
            "Cancer": cancer,
            "pTPM": values,
        }
    )


def test_ihc_counts_are_patients_and_unknown_is_not_negative():
    raw = pd.DataFrame(
        [
            ihc_row("lung cancer"),
            ihc_row("breast cancer", High=0, Medium=0, Low=0),
            ihc_row("renal cancer", High=0, Medium=0, Low=0, **{"Not detected": 0}),
            ihc_row("thyroid cancer", **{"Not detected": np.nan}),
        ]
    )
    result = hpa_cancer._summarize_ihc(raw).set_index("cancer")
    assert result.loc["lung cancer", "total"] == 10
    assert result.loc["lung cancer", "detected"] == 6
    assert result.loc["lung cancer", "prevalence_detected"] == 0.6
    assert result.loc["lung cancer", "prevalence_medium_high"] == 0.3
    assert result.loc["breast cancer", "prevalence_detected"] == 0
    assert result.loc["renal cancer", "total"] == 0
    assert pd.isna(result.loc["renal cancer", "prevalence_detected"])
    assert pd.isna(result.loc["thyroid cancer", "total"])
    assert pd.isna(result.loc["thyroid cancer", "prevalence_detected"])


@pytest.mark.parametrize("value", [-1, 0.5, np.inf, "nonsense"])
def test_invalid_ihc_counts_rejected(value):
    with pytest.raises((ValueError, TypeError)):
        hpa_cancer._summarize_ihc(pd.DataFrame([ihc_row(High=value)]))


def test_duplicate_ihc_rejected():
    with pytest.raises(ValueError, match="unique"):
        hpa_cancer._summarize_ihc(pd.DataFrame([ihc_row(), ihc_row()]))


def test_rna_measured_denominator_and_thresholds():
    rows = rna_rows([0, 0.1, 1, 2, 5, np.nan])
    row = hpa_cancer._summarize_rna_rows(rows, hpa_cancer.hpa_cancer_rna_cohorts()).iloc[0]
    assert row["source_samples"] == 6
    assert row["samples"] == 5
    assert row["missing_samples"] == 1
    assert row["mean_ptpm"] == pytest.approx(8.1 / 5)
    assert row["expressed_samples_ptpm_ge_0_1"] == 4
    assert row["prevalence_ptpm_ge_1"] == 3 / 5
    assert row["prevalence_ptpm_ge_2"] == 2 / 5
    assert row["prevalence_ptpm_ge_5"] == 1 / 5


def test_rna_all_missing_and_measured_zero_are_distinct():
    data = pd.concat([rna_rows([np.nan]), rna_rows([0], "Lung Squamous Cell Carcinoma (TCGA)")])
    rows = hpa_cancer._summarize_rna_rows(data, hpa_cancer.hpa_cancer_rna_cohorts()).set_index(
        "cancer_code"
    )
    assert rows.loc["LUAD", "samples"] == 0
    assert pd.isna(rows.loc["LUAD", "mean_ptpm"])
    assert pd.isna(rows.loc["LUAD", "prevalence_ptpm_ge_1"])
    assert rows.loc["LUSC", "prevalence_ptpm_ge_1"] == 0


@pytest.mark.parametrize("value", [-1, np.inf, -np.inf, "invalid"])
def test_invalid_rna_rejected(value):
    with pytest.raises(ValueError):
        hpa_cancer._summarize_rna_rows(rna_rows([value]), hpa_cancer.hpa_cancer_rna_cohorts())


def test_rna_duplicate_identity_rejected_but_cohorts_stay_separate():
    rows = rna_rows([1, 2])
    with pytest.raises(ValueError, match="unique"):
        hpa_cancer._summarize_rna_rows(pd.concat([rows, rows]), hpa_cancer.hpa_cancer_rna_cohorts())
    both = pd.concat([rows, rows.assign(Cancer="Lung Adenocarcinoma (validation)")])
    result = hpa_cancer._summarize_rna_rows(both, hpa_cancer.hpa_cancer_rna_cohorts())
    assert set(result["cohort"]) == {"TCGA", "validation"}
    assert list(result["samples"]) == [2, 2]
    with pytest.raises(ValueError, match="unrecognized"):
        hpa_cancer._summarize_rna_rows(
            rows.assign(Cancer="Unknown (TCGA)"), hpa_cancer.hpa_cancer_rna_cohorts()
        )


def comparison_rna():
    return pd.DataFrame(
        [
            {
                "gene_id": "ENSG1",
                "cancer_code": code,
                "cancer": label,
                "cohort": cohort,
                "samples": n,
                "sum_ptpm": mean * n,
                "expressed_samples_ptpm_ge_1": pos,
            }
            for code, label, cohort, n, mean, pos in [
                ("LUAD", "LUAD (TCGA)", "TCGA", 100, 2, 40),
                ("LUSC", "LUSC (TCGA)", "TCGA", 10, 9, 9),
                ("LUAD", "LUAD (validation)", "validation", 1, 200, 1),
            ]
        ]
    )


def compare(rna, cohort="TCGA"):
    ihc = hpa_cancer._summarize_ihc(
        pd.DataFrame(
            [ihc_row(c) for c in ["lung cancer", "glioma", "skin cancer", "breast cancer"]]
        )
    )
    return hpa_cancer._compare(
        ihc, rna, hpa_cancer.hpa_cancer_crosswalk(), cohort=cohort, threshold=1
    ).set_index("cancer")


def test_comparison_pools_by_sample_count_and_exposes_mapping_gaps():
    result = compare(comparison_rna())
    row = result.loc["lung cancer"]
    assert row["rna_samples"] == 110
    assert row["rna_mean_ptpm"] == pytest.approx(290 / 110)
    assert row["rna_prevalence"] == pytest.approx(49 / 110)
    assert row["comparison_status"] == "comparable"
    assert not row["cohorts_paired"]
    assert result.loc["glioma", "comparison_status"] == "scope_mismatch"
    assert result.loc["skin cancer", "comparison_status"] == "unmatched"
    assert result.loc["breast cancer", "comparison_status"] == "missing_rna"
    assert pd.isna(result.loc["glioma", "rna_samples"])


def test_partial_type_pool_is_missing_and_validation_is_not_supplemented_with_tcga():
    result = compare(comparison_rna(), "validation")
    row = result.loc["lung cancer"]
    assert row["rna_available_types"] == 1
    assert row["comparison_status"] == "incomplete_rna_types"
    assert pd.isna(row["rna_prevalence"])
    data = comparison_rna()
    data.loc[data.cancer_code == "LUSC", "samples"] = 0
    assert compare(data).loc["lung cancer", "comparison_status"] == "no_measured_rna"


def test_crosswalk_covers_all_types_without_claiming_gbm_is_all_glioma():
    crosswalk = hpa_cancer.hpa_cancer_crosswalk().set_index("cancer")
    assert len(crosswalk) == 20
    assert crosswalk.index.is_unique
    assert crosswalk.loc["colorectal cancer", "rna_codes"] == "COAD|READ"
    assert crosswalk.loc["renal cancer", "rna_codes"] == "KICH|KIRC|KIRP"
    assert crosswalk.loc["lung cancer", "rna_codes"] == "LUAD|LUSC"
    assert crosswalk.loc["glioma", "mapping_status"] == "scope_mismatch"
    assert set(crosswalk.loc[crosswalk.mapping_status == "unmatched"].index) == {
        "carcinoid",
        "lymphoma",
        "skin cancer",
    }
    codes = set("|".join(crosswalk.rna_codes).split("|")) - {""}
    assert codes == set(hpa_cancer.hpa_cancer_rna_cohorts().cancer_code)


def test_builder_keeps_gene_together_across_chunks_and_rejects_repeated_blocks(tmp_path):
    path = tmp_path / "rna.tsv"
    rows = pd.concat([rna_rows([0, 1, 2]), rna_rows([3, 4]).assign(Gene="ENSG2")])
    rows.to_csv(path, sep="\t", index=False)
    chunks = list(complete_gene_chunks(path, chunksize=2))
    assert [len(c) for c in chunks] == [3, 2]
    duplicate = pd.concat([rows, rows.iloc[:1]])
    duplicate.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="gene-contiguous"):
        list(complete_gene_chunks(path, chunksize=2))


def test_duplicate_sample_across_parser_boundary_rejected(tmp_path):
    path = tmp_path / "rna.tsv"
    rows = rna_rows([1, 2])
    pd.concat([rows, rows.iloc[:1]]).to_csv(path, sep="\t", index=False)
    batch = next(complete_gene_chunks(path, chunksize=2))
    with pytest.raises(ValueError, match="unique"):
        hpa_cancer._summarize_rna_rows(batch, hpa_cancer.hpa_cancer_rna_cohorts())


def test_source_identity_and_legacy_release_unknown():
    sources = hpa_cancer.hpa_cancer_sources()
    assert sources["hpa_release"] == "25.1"
    assert sources["ensembl_release"] == 109
    assert sources["scope"].startswith("genome-wide")
    for row in sources["raw_sources"].values():
        assert len(row["sha256"]) == 64
        assert row["acquired_at"]
        assert row["url"].startswith("https://www.proteinatlas.org/")
    for row in sources["legacy_tsarina"].values():
        assert row["hpa_release"] is None
        assert row["release_status"] == "unrecorded"
    sources["hpa_release"] = "changed"
    assert hpa_cancer.hpa_cancer_sources()["hpa_release"] == "25.1"
    cautions = hpa_cancer.hpa_cancer_assay_limitations()
    assert (
        cautions.loc[cautions.gene_id == "ENSG00000176746", "scope"].iloc[0]
        == "normal_tissue_curation_not_cancer_revalidation"
    )


@pytest.mark.parametrize("archive_format", ["zip", "file"])
def test_pinned_download_rejects_wrong_bytes_and_preserves_old_file(
    monkeypatch, tmp_path, archive_format
):
    monkeypatch.setenv("CANCERDATA_DATA_DIR", str(tmp_path))
    content = b"Gene\tCancer\nA\tB\n"
    if archive_format == "zip":
        stream = io.BytesIO()
        with zipfile.ZipFile(stream, "w") as z:
            z.writestr("test.tsv", content)
        archive = stream.getvalue()
    else:
        archive = content
    pin = {
        "bytes": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
        "archive_sha256": hashlib.sha256(archive).hexdigest(),
    }
    monkeypatch.setitem(
        reference_data.REFERENCE_SOURCES,
        "pinned",
        {
            "filename": "test.tsv",
            "description": "test",
            "default_version": "v1",
            "urls": {"v1": "https://example.org/source"},
            "pins": {"v1": pin},
            "archive_format": archive_format,
        },
    )
    monkeypatch.setattr(reference_data.urllib.request, "urlopen", lambda url: io.BytesIO(archive))
    path = reference_data.download("pinned")
    assert path.read_bytes() == content
    assert reference_data.verify("pinned")
    # A matching mutable manifest cannot bless altered bytes under the same pin.
    path.write_bytes(b"X" * len(content))
    manifest = reference_data._read_manifest()
    manifest["pinned@v1"]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    reference_data._write_manifest(manifest)
    assert not reference_data._cached_file_ok("pinned", "v1", path)
    assert not reference_data.verify("pinned")
    assert reference_data.ensure("pinned").read_bytes() == content
    monkeypatch.setattr(reference_data.urllib.request, "urlopen", lambda url: io.BytesIO(b"wrong"))
    with pytest.raises(reference_data.ReferenceDataError, match="archive checksum"):
        reference_data.download("pinned", force=True)
    assert path.read_bytes() == content
    assert not list(path.parent.glob("*.part"))
    # Content pin also matters, independently of the transport/archive hash.
    monkeypatch.setattr(reference_data.urllib.request, "urlopen", lambda url: io.BytesIO(archive))
    pin["sha256"] = "0" * 64
    with pytest.raises(reference_data.ReferenceDataError, match="content checksum"):
        reference_data.download("pinned", force=True)
    assert path.read_bytes() == content


def test_public_apis_select_ids_and_do_not_merge_cohorts(monkeypatch, tmp_path):
    ihc = tmp_path / "ihc.tsv"
    pd.DataFrame([ihc_row(), {**ihc_row(), "Gene": "ENSG2"}]).to_csv(ihc, sep="\t", index=False)
    rna = tmp_path / "rna.tsv"
    comparison_rna().to_csv(rna, sep="\t", index=False)
    monkeypatch.setattr(
        reference_data, "ensure", lambda name: ihc if name == "hpa_cancer_ihc" else rna
    )
    assert list(hpa_cancer.hpa_cancer_ihc_prevalence("ENSG1").gene_id) == ["ENSG1"]
    assert hpa_cancer.hpa_cancer_ihc_prevalence("ABSENT").empty
    assert len(hpa_cancer.hpa_cancer_rna_prevalence(cohort="TCGA")) == 2
    result = hpa_cancer.hpa_cancer_rna_ihc_comparison("ENSG1")
    assert result.iloc[0].rna_samples == 110
    assert hpa_cancer.hpa_cancer_rna_ihc_comparison("ABSENT").empty
    assert hpa_cancer.hpa_cancer_rna_ihc_comparison([]).empty
    pd.testing.assert_frame_equal(result, hpa_cancer.hpa_cancer_rna_ihc_comparison(iter(["ENSG1"])))
    with pytest.raises(ValueError, match="threshold"):
        hpa_cancer.hpa_cancer_rna_ihc_comparison(threshold=3)
    with pytest.raises(ValueError, match="cohort"):
        hpa_cancer.hpa_cancer_rna_ihc_comparison(cohort="all")


def test_opt_in_keeps_rna_only_genes_and_missing_is_not_zero(monkeypatch, tmp_path):
    ihc = tmp_path / "ihc.tsv"
    pd.DataFrame([ihc_row(High=0, Medium=0, Low=0)]).to_csv(ihc, sep="\t", index=False)
    rna = tmp_path / "rna.tsv"
    data = pd.concat([comparison_rna(), comparison_rna().assign(gene_id="RNA_ONLY")])
    data.to_csv(rna, sep="\t", index=False)
    monkeypatch.setattr(
        reference_data, "ensure", lambda name: ihc if name == "hpa_cancer_ihc" else rna
    )
    assert set(hpa_cancer.hpa_cancer_rna_ihc_comparison().gene_id) == {"ENSG1"}
    frame = hpa_cancer.hpa_cancer_rna_ihc_comparison(include_missing_ihc=True)
    assert len(frame) == 40
    assert not frame.duplicated(["gene_id", "cancer"]).any()
    r = frame.set_index(["gene_id", "cancer"])
    measured = r.loc[("ENSG1", "lung cancer")]
    absent = r.loc[("RNA_ONLY", "lung cancer")]
    assert measured.comparison_status == "comparable"
    assert measured.prevalence_detected == 0
    assert absent.measurement_status == absent.comparison_status == "missing_ihc"
    assert absent.rna_samples == 110 and absent.rna_prevalence == pytest.approx(49 / 110)
    assert (
        absent[
            ["total", "detected", "prevalence_detected", "high", "medium", "low", "not_detected"]
        ]
        .isna()
        .all()
    )
    assert r.loc[("RNA_ONLY", "glioma"), "comparison_status"] == "scope_mismatch"
    assert r.loc[("RNA_ONLY", "glioma"), "measurement_status"] == "missing_ihc"
    assert pd.isna(r.loc[("RNA_ONLY", "glioma"), "rna_prevalence"])
    assert r.loc[("RNA_ONLY", "skin cancer"), "comparison_status"] == "unmatched"
    validation = hpa_cancer.hpa_cancer_rna_ihc_comparison(
        "RNA_ONLY", cohort="validation", include_missing_ihc=True
    ).set_index("cancer")
    assert validation.loc["lung cancer", "comparison_status"] == "incomplete_rna_types"
    assert pd.isna(validation.loc["lung cancer", "rna_prevalence"])
    for selection in ([], ["ABSENT"]):
        assert hpa_cancer.hpa_cancer_rna_ihc_comparison(selection, include_missing_ihc=True).empty
    selected = hpa_cancer.hpa_cancer_rna_ihc_comparison(
        iter(["RNA_ONLY"]), include_missing_ihc=True
    )
    assert len(selected) == 20 and set(selected.gene_id) == {"RNA_ONLY"}

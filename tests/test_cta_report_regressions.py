"""Regression cases for the CTA analysis and report review findings."""

import gzip
import json

import cta_decision_guide as guide
import cta_mortality_coverage_report as mortality
import cta_proteoform_report as proteins
import cta_report_common as common
import cta_report_protein_lengths as lengths
import cta_report_render as render
import cta_threshold_report as genes
import numpy as np
import pandas as pd
import pytest


def test_fixed_background_rejects_annotation_padding_effect():
    background = np.array([[0.0], [1.0], [2.0], [4.0], [10.0]])
    original, cutoff = genes.threshold_metrics(background, [3], (90,), background=background)
    padded = np.vstack([background, np.zeros((1000, 1))])
    changed, changed_cutoff = genes.threshold_metrics(padded, [3], (90,), background=background)
    np.testing.assert_equal(cutoff, changed_cutoff)
    assert original[0].n_expressing.equals(changed[0].n_expressing)
    naive, _ = genes.threshold_metrics(padded, [3], (90,))
    assert not naive[0].n_expressing.equals(original[0].n_expressing)


def test_common_background_fails_on_missing_measurement():
    frame = pd.DataFrame({common.ID: ["A", "B"], "P": [1.0, np.nan]})
    with pytest.raises(ValueError, match="missing"):
        common.background_values(frame, ["A", "B"], ["P"])


def test_tiny_cohort_cannot_enter_primary_gene_selection():
    frame = pd.DataFrame(
        {"complete_measurement": [True] * 3, "n_patients": [1, 9, 10], "n_expressing": [1, 9, 10]}
    )
    assert genes.qualifies(frame, 75).tolist() == [False, False, True]


def test_primary_and_metastatic_tumors_are_not_averaged():
    samples = ["TCGA-AA-0001-01A", "TCGA-AA-0001-06A", "TCGA-AA-0002-06A"]
    values, patients, audit = genes.group_patients(
        np.array([[1.0, 100.0, 20.0]]), samples, "TCGA", {}
    )
    assert patients == ["TCGA-AA-0001", "TCGA-AA-0002"]
    assert values.tolist() == [[1.0, 20.0]]
    assert audit.included_in_patient_group.tolist() == [True, False, True]


def test_resume_rejects_changed_inputs_or_corrupt_outputs(tmp_path):
    prefix = tmp_path / "COHORT"
    data = prefix.with_suffix(".parquet")
    data.write_bytes(b"original")
    common.write_checkpoint(prefix, {}, "source-and-policy-v1", [".parquet"])
    assert common.checkpoint_valid(prefix, "source-and-policy-v1", [".parquet"])
    assert not common.checkpoint_valid(prefix, "source-and-policy-v2", [".parquet"])
    data.write_bytes(b"changed")
    assert not common.checkpoint_valid(prefix, "source-and-policy-v1", [".parquet"])


def test_failed_or_stale_analysis_is_not_rendered(tmp_path):
    validation = tmp_path / "validation.json"
    validation.write_text(json.dumps({"status": "failed"}))
    with pytest.raises(ValueError, match="did not pass"):
        common.verify_analysis(tmp_path)
    validation.write_text(json.dumps({"status": "passed"}))
    data = tmp_path / "data.csv"
    data.write_text("current")
    common.seal_stage(tmp_path, "analysis", [], [validation, data])
    common.verify_analysis(tmp_path)
    data.write_text("stale")
    with pytest.raises(ValueError, match="Stale"):
        common.verify_analysis(tmp_path)


def test_readme_start_section_preserves_later_burden_documentation():
    before = "# Panel\n\n## Start here\nold\n\n## Incidence and mortality represented\ncaveats and rebuild command\n"
    after = guide.replace_start_section(before, "new guide")
    assert "old" not in after
    assert after.endswith("## Incidence and mortality represented\ncaveats and rebuild command\n")
    assert guide.replace_start_section(after, "new guide") == after


def test_evidence_sections_follow_metadata_after_page_reordering():
    index = pd.DataFrame(
        {
            "pdf_page": [1, 2, 3],
            "category": ["safety", "methods", "safety"],
            "title": ["lung", "rules", "heart"],
        }
    )
    groups = dict(guide.supporting_sections(index, 3))
    assert groups["safety"].title.tolist() == ["lung", "heart"]
    with pytest.raises(ValueError, match="page index"):
        guide.supporting_sections(index, 4)


def test_new_figure_category_is_included(tmp_path):
    (tmp_path / "plots").mkdir()
    (tmp_path / "plots/safety.png").write_bytes(b"placeholder")
    pd.DataFrame({"name": ["safety"], "title": ["Safety"], "category": ["new-category"]}).to_csv(
        tmp_path / "plot_index.csv", index=False
    )
    assert render.figure_rows(tmp_path).name.tolist() == ["safety"]


def test_missing_or_partial_protein_is_not_false():
    result = proteins.patient_positive_status(
        pd.Series([np.nan, 10.0, 0.0, 10.0]),
        pd.Series([1.0] * 4),
        pd.Series([True, False, True, True]),
    )
    assert result.iloc[:2].isna().all()
    assert result.iloc[2:].tolist() == [False, True]


def test_mortality_pancreas_uses_adenocarcinoma_not_neuroendocrine():
    frame = pd.DataFrame(
        {
            "burden_category": ["pancreas"] * 2,
            "cancer_code": ["NET_PANCREAS", "PAAD"],
            "n_patients": [30, 177],
            "fraction": [0.167, 0.006],
        }
    )
    assert mortality.best_observed(frame, "fraction").cancer_code == "PAAD"
    assert mortality.best_observed(frame.iloc[:1], "fraction") is None


def test_thymoma_is_not_head_neck_evidence():
    frame = pd.DataFrame(
        {
            "burden_category": ["head_and_neck"] * 2,
            "cancer_code": ["THYM", "HNSC"],
            "n_patients": [119, 500],
            "fraction": [0.9, 0.3],
        }
    )
    assert mortality.best_observed(frame, "fraction").cancer_code == "HNSC"


def test_wrong_fasta_cannot_be_stamped_release_112(tmp_path):
    fasta = tmp_path / "release115.fa.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">protein gene:G transcript:T\nMAAA\n")
    with pytest.raises(ValueError, match="pinned Ensembl 112"):
        lengths.build_annotations(pd.DataFrame({lengths.ID: ["G"]}), fasta)


def test_archives_and_gzip_are_byte_reproducible(tmp_path):
    data = pd.DataFrame({"gene": ["G"], "value": [1.25]})
    first, second = tmp_path / "a.csv.gz", tmp_path / "b.csv.gz"
    common.write_csv(data, first)
    common.write_csv(data, second)
    assert first.read_bytes() == second.read_bytes()
    dest = tmp_path / "bundle.zip"
    common.deterministic_zip(tmp_path, dest)
    before = dest.read_bytes()
    first.touch()
    common.deterministic_zip(tmp_path, dest)
    assert before == dest.read_bytes()


def test_pdf_rendering_is_byte_reproducible(tmp_path):
    for filename in ("a.pdf", "b.pdf"):
        report = render.Report(tmp_path / filename)
        report.start("Test")
        report.paragraph("Current measured evidence")
        report.finish("Test", "methods")
        report.save()
    assert (tmp_path / "a.pdf").read_bytes() == (tmp_path / "b.pdf").read_bytes()


def test_rounded_zero_maximum_names_all_ties_and_its_bound():
    from cta_hpa_tissue_report import tissue_maximum

    result = tissue_maximum(pd.Series([0.0, 0.0], index=["lung", "liver"]), 2)
    assert result == {
        "value": 0.0,
        "tissues": "liver;lung",
        "status": "reported_zero",
        "zero_upper_bound": 0.1,
    }
    missing = tissue_maximum(pd.Series([np.nan], index=["lung"]), 2)
    assert missing["status"] == "unavailable"
    assert pd.isna(missing["value"])


def test_resume_keeps_cache_metadata_out_of_scientific_exports(tmp_path):
    prefix = tmp_path / "COHORT"
    prefix.with_suffix(".parquet").write_bytes(b"scientific table")
    original = {"n_patients": 20, "n_common_background_genes": 7876}
    common.write_checkpoint(prefix, original, "input fingerprint", [".parquet"])
    assert common.checkpoint_payload(prefix) == original
    assert common.checkpoint_valid(prefix, "input fingerprint", [".parquet"])

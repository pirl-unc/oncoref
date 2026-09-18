# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Evidence-synthesis contracts: an RNA zero is not negative IHC, a missing
measurement is not a zero, and no negative assay result promotes a gene."""

import pandas as pd

from oncoref import cta, cta_review


def _tables(bulk_rna=None, ihc=None, single_cell_rna=None):
    """Minimal source-shaped tables, so synthesis is tested without the atlas."""
    return {
        "bulk_rna": pd.DataFrame(bulk_rna or [], columns=["Gene", "Gene name", "Tissue", "nTPM"]),
        "ihc": pd.DataFrame(
            ihc or [],
            columns=["Gene", "Gene name", "Tissue", "Cell type", "Level", "Reliability"],
        ),
        "single_cell_rna": pd.DataFrame(
            single_cell_rna or [], columns=["Gene", "Gene name", "Cell type", "nTPM"]
        ),
    }


def _synthesize(**tables):
    universe = pd.DataFrame({"Ensembl_Gene_ID": ["G1"]})
    return cta_review._synthesize_atlas(universe, _tables(**tables)).iloc[0]


def test_grouping_survives_the_mapping_protocol():
    # A groupby's ``keys`` attribute is the grouping label, so dict(groupby)
    # takes the mapping path and calls that string. Regression: synthesis must
    # still route each gene's rows to that gene.
    row = _synthesize(
        bulk_rna=[["G1", "A", "heart muscle", 7.0], ["G2", "B", "heart muscle", 99.0]]
    )
    assert row["heart_rna_max_ntpm"] == 7.0


def test_reported_zero_is_not_negative_protein_evidence():
    row = _synthesize(
        bulk_rna=[["G1", "A", "heart muscle", 0.0]],
        ihc=[["G1", "A", "heart muscle", "myocytes", "Not detected", "Enhanced"]],
    )
    # An RNA zero is an estimate; only the IHC annotation reports non-detection.
    assert row["somatic_rna_status"] == "measured"
    assert row["somatic_rna_max_ntpm"] == 0.0
    assert row["cardiomyocyte_ihc_status"] == "not_detected"
    # Nothing was detected anywhere, so no warning is asserted.
    assert row["atlas_warning_codes"] == ""


def test_absent_measurement_is_unavailable_not_zero():
    row = _synthesize(bulk_rna=[["G1", "A", "testis", 40.0]])
    assert row["somatic_rna_status"] == "unavailable"
    assert pd.isna(row["somatic_rna_max_ntpm"])
    assert row["somatic_ihc_status"] == "unavailable"
    # An unavailable measurement must be reported as a gap, never as evidence.
    assert "somatic_ihc" in row["atlas_evidence_gaps"]
    assert "bulk_rna" not in row["atlas_evidence_gaps"]


def test_reproductive_tissue_is_excluded_from_somatic_summaries():
    row = _synthesize(bulk_rna=[["G1", "A", "testis", 40.0], ["G1", "A", "lung", 1.0]])
    assert row["bulk_rna_max_ntpm"] == 40.0
    assert row["somatic_rna_max_ntpm"] == 1.0


def test_rna_ihc_discordance_is_named_not_resolved():
    row = _synthesize(
        bulk_rna=[["G1", "A", "heart muscle", 5.1]],
        ihc=[["G1", "A", "heart muscle", "myocytes", "Not detected", "Enhanced"]],
    )
    codes = row["atlas_warning_codes"].split(";")
    # The conflict is reported; negative IHC does not cancel the RNA estimate.
    assert "somatic_rna_estimate" in codes
    assert "rna_ihc_discordance" in codes
    assert f"heart_rna_ge_{cta_review.SAFETY_NTPM_THRESHOLD:g}_ntpm" in codes


def test_somatic_protein_detection_is_reported():
    row = _synthesize(
        ihc=[["G1", "A", "liver", "hepatocytes", "Medium", "Supported"]],
    )
    assert row["somatic_ihc_status"] == "detected"
    assert "somatic_protein_detected" in row["atlas_warning_codes"]
    assert row["somatic_ihc_detected_tissues"] == "liver"


def test_cardiomyocyte_ihc_accepts_either_source_label():
    for label in ("myocytes", "cardiomyocytes"):
        row = _synthesize(
            ihc=[["G1", "A", "heart muscle", label, "Not detected", "Enhanced"]],
        )
        assert row["cardiomyocyte_ihc_status"] == "not_detected", label


def test_warning_tier_is_opt_in_and_disjoint_from_strict():
    strict = cta.cta_gene_names()
    warnings = cta.cta_warning_gene_names()
    assert warnings
    assert not (warnings & strict)
    assert cta.cta_gene_names(include_warnings=True) == strict | warnings
    assert len(cta.cta_warning_gene_ids()) == len(warnings)


def test_reviewed_modalities_are_all_consumed():
    # A typo'd modality would load fine and never surface in the summary.
    reviewed = cta_review.cta_reviewed_evidence()
    assert set(reviewed["modality"]) <= set(cta_review.REVIEWED_MODALITIES)


def test_reviewed_evidence_rows_state_their_provenance_and_limits():
    reviewed = cta_review.cta_reviewed_evidence()
    for column in ("assay", "scope", "source_version", "source_anchor", "limitations"):
        assert reviewed[column].str.strip().ne("").all(), column
    # Reviewed rows are supplemental observations about real candidates.
    universe = cta_review._universe()
    assert set(reviewed["Ensembl_Gene_ID"]) <= set(universe["Ensembl_Gene_ID"])


def test_unreviewed_modalities_are_marked_not_reviewed():
    summary = cta_review.cta_evidence_summary()
    reviewed = cta_review.cta_reviewed_evidence()
    for modality in cta_review.REVIEWED_MODALITIES:
        expected = set(reviewed.loc[reviewed["modality"].eq(modality), "Ensembl_Gene_ID"])
        status = summary.set_index("Ensembl_Gene_ID")[f"{modality}_review_status"]
        assert set(status.loc[list(expected)].unique()) <= {"reviewed"}
        assert (status.drop(index=list(expected)) == "not_reviewed").all()


def test_partial_safety_mapping_caveats_a_non_detection():
    summary = cta_review.cta_evidence_summary()
    partial = summary["brain_ihc_mapping_coverage"].ne("complete")
    negative = summary["brain_ihc_status"].ne("detected")
    caveated = summary["atlas_evidence_gaps"].str.contains("brain_ihc_partial_mapping")
    # HPA v23 stains 9 of 14 requested brain regions, so a non-detection there
    # must not read as group-wide absence.
    assert (partial & negative).any()
    assert caveated.eq(partial & negative).all()
    # The regions behind the caveat are named, not merely counted.
    unmapped = summary.loc[partial, "brain_ihc_unmapped_tissues"].iloc[0].split(";")
    assert "spinal cord" in unmapped


def test_detection_does_not_need_the_coverage_caveat():
    # Positive protein evidence stands on its own; an unstained region elsewhere
    # cannot turn a detection into a gap.
    row = _synthesize(ihc=[["G1", "A", "cerebral cortex", "neurons", "High", "Enhanced"]])
    assert row["brain_ihc_status"] == "detected"
    assert row["brain_ihc_mapping_coverage"] == "partial"
    assert "brain_ihc_partial_mapping" not in row["atlas_evidence_gaps"]


def test_versioned_gene_ids_still_join(monkeypatch):
    # A curated row written as ENSG...14 must not merge onto nothing and then
    # read as "not reviewed"; that would hide a review rather than error.
    reviewed = cta_review.cta_reviewed_evidence()
    reviewed["Ensembl_Gene_ID"] = reviewed["Ensembl_Gene_ID"] + ".14"
    monkeypatch.setattr(cta_review, "cta_reviewed_evidence", lambda: reviewed)
    summary = cta_review.cta_evidence_summary()
    status = summary.set_index("Ensembl_Gene_ID")["isoform_review_status"]
    assert status["ENSG00000126890"] == "reviewed"


def test_summary_carries_no_unversioned_measurement_columns():
    summary = cta_review.cta_evidence_summary()
    # Curation tables bring their own HPA numbers of unstated version. Seating
    # those beside the pinned atlas columns produced near-anagram pairs such as
    # rna_heart_max_ntpm vs heart_rna_max_ntpm that silently disagreed on the
    # watchlist candidates.
    for leaked in (
        "rna_heart_max_ntpm",
        "rna_brain_max_ntpm",
        "rna_max_somatic_ntpm",
        "hpa_testis_ntpm",
        "hpa_max_somatic_ntpm",
    ):
        assert leaked not in summary.columns, leaked
    assert set(cta_review._universe().columns) == {
        "Symbol",
        "Ensembl_Gene_ID",
        "candidate_origin",
    }


def test_review_module_is_a_public_facade():
    import oncoref

    assert oncoref.cta_review is cta_review
    assert "cta_review" in oncoref.__all__
    for name in ("cta_evidence_summary", "cta_normal_tissue_evidence", "cta_reviewed_evidence"):
        assert callable(getattr(oncoref.cta_review, name))
    # Newer names stay in the semantic module rather than the flat namespace.
    assert not hasattr(oncoref, "cta_evidence_summary")


def test_every_candidate_gets_a_comparable_summary():
    summary = cta_review.cta_evidence_summary()
    assert summary["Ensembl_Gene_ID"].is_unique
    assert summary["evidence_summary"].str.strip().ne("").all()
    # The caveat travels with every row, including the clean ones.
    assert summary["evidence_summary"].str.contains("do not establish").all()
    assert set(summary["discovery_tier"]) <= {
        "strict",
        "warning",
        "low_expression",
        "candidate",
        "excluded",
    }

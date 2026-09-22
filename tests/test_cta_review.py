# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Evidence-synthesis contracts: an RNA zero is not negative IHC, a missing
measurement is not a zero, and no negative assay result promotes a gene."""

import pandas as pd
import pytest

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


def test_mixed_negative_and_unscored_ihc_is_incomplete():
    # HPA leaves some cell types unscored. A partial negative is not a clean
    # non-detection, so it is reported as a gap rather than as absence.
    row = _synthesize(
        ihc=[
            ["G1", "A", "liver", "hepatocytes", "Not detected", "Enhanced"],
            ["G1", "A", "liver", "bile duct cells", None, "Enhanced"],
        ],
    )
    assert row["somatic_ihc_status"] == "incomplete"
    assert "somatic_ihc" in row["atlas_evidence_gaps"]
    assert row["somatic_ihc_not_detected_rows"] == 1
    assert row["somatic_ihc_unavailable_rows"] == 1


def test_long_form_evidence_keeps_assay_vocabularies_separate():
    df = cta_review.cta_normal_tissue_evidence()
    status = df.groupby("modality")["measurement_status"].agg(set)
    # An RNA zero is an estimate; only IHC can report a non-detection.
    assert "reported_zero" in status["bulk_rna"]
    assert "not_detected" not in status["bulk_rna"]
    assert "not_detected" in status["ihc"]
    assert "reported_zero" not in status["ihc"]
    # Every row names the pinned release it came from.
    assert df["source_version"].eq("v23").all()
    assert df["source_url"].str.startswith("http").all()
    assert set(df["unit"]) == {"nTPM", "IHC category"}
    # Single-cell types are aggregated across organs, so they carry no tissue.
    assert df.loc[df["modality"].eq("single_cell_rna"), "tissue"].isna().all()


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


def test_partial_mapping_is_reported_for_every_candidate():
    summary = cta_review.cta_evidence_summary()
    # HPA v23 maps 8 of the 14 requested brain regions for IHC and 10 for RNA,
    # and routinely surveys 4 and 10 of those. Neither covers thalamus; IHC
    # additionally lacks spinal cord, which RNA does survey.
    for modality in ("rna", "ihc"):
        partial = summary[f"brain_{modality}_mapping_coverage"].ne("complete")
        caveated = summary["atlas_coverage_limits"].str.contains(
            f"brain_{modality}_partial_mapping"
        )
        assert partial.all(), modality
        assert caveated.all(), modality
    # The regions behind the caveat are named, not merely counted.
    assert "spinal cord" in summary["brain_ihc_unmapped_tissues"].iloc[0].split(";")
    assert "thalamus" in summary["brain_rna_unmapped_tissues"].iloc[0].split(";")
    # A completely mapped group stays silent, so the caveat means something.
    assert summary["heart_ihc_mapping_coverage"].eq("complete").all()
    assert not summary["atlas_coverage_limits"].str.contains("heart_ihc_partial_mapping").any()


def test_coverage_caveat_survives_a_positive_detection():
    # Coverage describes what was assayed, not what was concluded. Gating this
    # on the verdict would let one positive region imply the group was covered.
    row = _synthesize(ihc=[["G1", "A", "cerebral cortex", "neurons", "High", "Enhanced"]])
    assert row["brain_ihc_status"] == "detected"
    assert row["brain_ihc_mapping_coverage"] == "partial"
    assert "brain_ihc_partial_mapping" in row["atlas_coverage_limits"]


def test_never_assayed_safety_group_is_named_as_a_gap():
    # 150 candidates have no HPA antibody data at all. "Never assayed" must be
    # distinguishable from "assayed, nothing found" for every group, not just
    # the one that happens to be partially mapped.
    row = _synthesize(bulk_rna=[["G1", "A", "lung", 1.0]])
    gaps = row["atlas_evidence_gaps"].split(";")
    for group in ("lung", "liver", "pancreas", "heart"):
        assert row[f"{group}_ihc_status"] == "unavailable", group
        assert f"{group}_ihc_unavailable" in gaps, group


def test_gradient_ihc_staining_counts_as_protein_detected():
    # HPA's Ascending/Descending levels are gradient staining, i.e. protein
    # present. Treating them as unmeasured would drop a somatic protein
    # detection in the permissive direction.
    for level in ("Ascending", "Descending"):
        row = _synthesize(ihc=[["G1", "A", "liver", "hepatocytes", level, "Supported"]])
        assert row["somatic_ihc_status"] == "detected", level
        assert "somatic_protein_detected" in row["atlas_warning_codes"], level
    # An unusable annotation is neither a detection nor a non-detection.
    row = _synthesize(ihc=[["G1", "A", "liver", "hepatocytes", "Not representative", "Uncertain"]])
    assert row["somatic_ihc_status"] == "unavailable"


def test_placeholder_rows_are_not_measurements():
    # HPA marks "no antibody data" with an all-null tissue/cell-type/level row.
    long_form = cta_review.cta_normal_tissue_evidence()
    ihc = long_form.loc[long_form["modality"].eq("ihc")]
    assert not ihc[["tissue", "cell_type", "ihc_level"]].isna().all(axis=1).any()
    assert ihc["tissue"].notna().all()


def test_evidence_summary_states_warnings_and_gaps():
    summary = cta_review.cta_evidence_summary()
    row = summary.loc[summary["Symbol"].eq("CTAG2")].iloc[0]
    text = row["evidence_summary"]
    # The one-liner must not read as reassuring while sibling columns disagree.
    assert "rna_ihc_discordance" in text
    assert "brain_ihc_partial_mapping" in text
    assert "do not establish" in text
    flagged = summary["atlas_warning_codes"].ne("")
    assert summary.loc[flagged, "evidence_summary"].str.contains("Warnings:").all()
    # A zero RNA estimate is named as an estimate, never as a bare measurement.
    zero = summary["heart_rna_max_ntpm"].eq(0)
    assert zero.any()
    assert summary.loc[zero, "evidence_summary"].str.contains("estimated 0 nTPM").all()


def test_ctag2_clinical_review_reaches_summary_without_promoting_target():
    row = cta_review.cta_evidence_summary().set_index("Symbol").loc["CTAG2"]
    assert row["clinical_review_status"] == "reviewed"
    assert "cardiac arrest" in row["clinical_review"]
    assert "MAGE-A4" in row["clinical_review"]
    assert "titin" in row["clinical_review"]
    assert row["cardiomyocyte_ihc_status"] == "not_detected"
    assert row["heart_rna_max_ntpm"] == 5.1
    assert row["discovery_tier"] == "warning"
    assert "CTAG2" not in cta.cta_gene_names()
    assert "CTAG2" in cta.cta_clinical_target_gene_names()
    warning = cta.cta_warning_references().set_index("Symbol").loc["CTAG2"]
    assert "cardiac arrests" in warning["unresolved_evidence"]


def test_versioned_gene_ids_still_join(monkeypatch):
    # A curated row written as ENSG...14 must not merge onto nothing and then
    # read as "not reviewed"; that would hide a review rather than error.
    reviewed = cta_review.cta_reviewed_evidence()
    reviewed["Ensembl_Gene_ID"] = reviewed["Ensembl_Gene_ID"] + ".14"
    monkeypatch.setattr(cta_review, "cta_reviewed_evidence", lambda: reviewed)
    # The summary is memoized, so a swapped input needs the cache dropped or
    # this test would assert against the real table and pass for free.
    cta_review._evidence_summary_frame.cache_clear()
    try:
        summary = cta_review.cta_evidence_summary()
        status = summary.set_index("Ensembl_Gene_ID")["isoform_review_status"]
        assert status["ENSG00000126890"] == "reviewed"
    finally:
        cta_review._evidence_summary_frame.cache_clear()


def test_derived_frames_are_memoized_and_defensively_copied():
    first = cta_review.cta_evidence_summary()
    second = cta_review.cta_evidence_summary()
    # Callers get independent frames, so in-place mutation cannot leak.
    assert first is not second
    first.loc[0, "evidence_summary"] = "mutated"
    assert cta_review.cta_evidence_summary().loc[0, "evidence_summary"] != "mutated"
    assert cta_review._evidence_summary_frame.cache_info().hits > 0


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
    for name in (
        "cta_atlas_coverage",
        "cta_evidence_summary",
        "cta_normal_tissue_evidence",
        "cta_reviewed_evidence",
    ):
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


def test_partial_observation_of_a_mapped_scope_is_incomplete():
    # A gene measured in 3 of the 5 brain regions the source surveys must not
    # report a brain maximum as though the group had been surveyed. The second
    # gene sets the denominator, so it does not come from the gene under test.
    shared = ("cerebellum", "retina", "hypothalamus")
    rows = [["G1", "A", tissue, 1.5] for tissue in shared]
    rows += [["G2", "B", tissue, 2.0] for tissue in (*shared, "cerebral cortex", "amygdala")]
    row = _synthesize(bulk_rna=rows)
    assert row["brain_rna_measured_rows"] == 3
    assert row["brain_rna_expected_rows"] == 5
    assert row["brain_rna_status"] == "incomplete"
    assert "brain_rna_incomplete" in row["atlas_evidence_gaps"].split(";")


def test_fully_observed_scope_is_measured_not_incomplete():
    summary = cta_review.cta_evidence_summary()
    full = summary["brain_rna_measured_rows"].eq(summary["brain_rna_expected_rows"])
    assert full.any()
    assert summary.loc[full, "brain_rna_status"].eq("measured").all()
    # The incomplete state is reachable on real data, so the check is not inert.
    partial = summary["brain_rna_status"].eq("incomplete")
    assert partial.any()
    assert (
        summary.loc[partial, "brain_rna_measured_rows"]
        < summary.loc[partial, "brain_rna_expected_rows"]
    ).all()


def test_evidence_gaps_are_per_gene_and_coverage_limits_are_per_release():
    summary = cta_review.cta_evidence_summary()
    # A fixed limitation of the release belongs in its own field; repeated into
    # every row it would leave atlas_evidence_gaps never empty and unable to
    # distinguish a gene with missing data from one measured everywhere.
    assert summary["atlas_evidence_gaps"].eq("").any()
    assert summary["atlas_evidence_gaps"].nunique() > 1
    assert summary["atlas_coverage_limits"].nunique() == 1
    assert summary["atlas_coverage_limits"].str.contains("brain_ihc_partial_mapping").all()
    assert not summary["atlas_evidence_gaps"].str.contains("partial_mapping").any()


def test_warning_column_prefix_is_not_doubled():
    summary = cta_review.cta_evidence_summary()
    assert "warning_code" in summary.columns
    assert "warning_warning_code" not in summary.columns


def test_warning_reviews_name_real_candidates_with_agreeing_identifiers():
    """CI guard on curated data: this table alone widens the CTA gene set.

    A typo'd symbol would silently become a reviewed warning-tier CTA, and a
    Symbol/ID disagreement would make cta_warning_gene_names and
    cta_warning_gene_ids describe different genes while the disjointness test's
    equal-cardinality check still passed.
    """
    refs = cta.cta_warning_references()
    universe = cta_review._universe()
    ids = cta_review._strip_version(refs["Ensembl_Gene_ID"])
    assert set(ids) <= set(universe["Ensembl_Gene_ID"])
    assert set(refs["Symbol"]) <= set(universe["Symbol"])
    expected = universe.set_index("Ensembl_Gene_ID")["Symbol"]
    assert refs["Symbol"].tolist() == expected.loc[ids].tolist()


def test_review_level_vocabulary_deliberately_differs_from_the_strict_filter():
    """Pin the divergence so it stays a decision rather than drifting.

    The shipped CTA table scores protein detection with Low/Medium/High only,
    which gates gene inclusion. This module also counts HPA's gradient levels,
    because for reporting a gradient is protein present. The safety-relevant
    path is therefore the narrower one; widening the shared constant would
    change the strict default and require regenerating the table.
    """
    from oncoref.cta_tissues import PROTEIN_DETECTED_LEVELS

    assert set(PROTEIN_DETECTED_LEVELS) == {"Low", "Medium", "High"}
    assert set(cta_review._IHC_DETECTED_LEVELS) == {
        "Low",
        "Medium",
        "High",
        "Ascending",
        "Descending",
    }
    assert "Not representative" not in cta_review._IHC_DETECTED_LEVELS


def test_ntpm_phrase_distinguishes_zero_from_unmeasured():
    # The three cases must not collapse into one number in prose.
    assert cta_review._ntpm_phrase(float("nan")) == "unavailable"
    assert cta_review._ntpm_phrase(0) == "estimated 0 nTPM"
    assert cta_review._ntpm_phrase(5.1) == "5.1 nTPM"


def test_data_derived_denominators_are_not_self_fulfilling():
    """Pin the four scopes whose denominator comes from the supplied tables.

    Those counts are derived from the same fixture that supplies the numerator,
    so a single-gene fixture reports "measured" by construction. Giving a second
    gene a tissue the first lacks makes the denominator independent of the gene
    under test, which is what the real atlas does.
    """
    row = _synthesize(
        bulk_rna=[
            ["G1", "A", "lung", 1.0],
            ["G2", "B", "lung", 1.0],
            ["G2", "B", "liver", 2.0],
        ],
        single_cell_rna=[
            ["G1", "A", "hepatocytes", 1.0],
            ["G2", "B", "hepatocytes", 1.0],
            ["G2", "B", "cardiomyocytes", 3.0],
        ],
    )
    # G1 is measured in 1 of the 2 tissues and cell types the source knows.
    assert row["bulk_rna_measured_rows"] == 1
    assert row["bulk_rna_expected_rows"] == 2
    assert row["bulk_rna_status"] == "incomplete"
    assert row["somatic_rna_status"] == "incomplete"
    assert row["single_cell_rna_status"] == "incomplete"
    gaps = row["atlas_evidence_gaps"].split(";")
    assert "bulk_rna_incomplete" in gaps
    assert "somatic_rna_incomplete" in gaps
    assert "single_cell_rna_incomplete" in gaps


def test_gap_list_covers_every_status_the_record_reports():
    """Every status the synthesis records must reach the gap field.

    Statuses are discovered rather than enumerated, so adding a scope to the
    synthesis brings it under this invariant automatically. The only names
    subtracted are the review columns the merge adds afterwards, and those come
    from the module's own constant, not a list kept by hand here.
    """
    summary = cta_review.cta_evidence_summary()
    # Discovered from the frame, minus the review columns the merge adds after
    # synthesis -- derived from the module's own constant rather than listed by
    # hand, so a new scope is covered automatically. Enumerating the scopes
    # here would be the same hand-kept list that fell behind twice.
    added_by_merge = {f"{modality}_review_status" for modality in cta_review.REVIEWED_MODALITIES}
    statuses = [c for c in summary.columns if c.endswith("_status") and c not in added_by_merge]
    assert len(statuses) > 10
    assert not any(c.endswith("_review_status") for c in statuses)
    for column in statuses:
        prefix = column[: -len("_status")]
        flagged = summary[column].isin(["unavailable", "incomplete"])
        if not flagged.any():
            continue
        for status in ("unavailable", "incomplete"):
            expected = summary[column].eq(status)
            if column.endswith("_ihc_status") and status == "incomplete":
                expected |= summary[column].eq("detected") & (
                    summary[f"{prefix}_measured_tissues"].lt(summary[f"{prefix}_expected_tissues"])
                    | summary[f"{prefix}_unavailable_rows"].gt(0)
                )
            if not expected.any():
                continue
            token = f"{prefix}_{status}"
            named = (
                summary["atlas_evidence_gaps"]
                .str.split(";")
                .apply(lambda gaps, name=token: name in gaps)
            )
            assert named.eq(expected).all(), token


def test_non_detection_requires_a_full_survey_not_an_unreachable_one():
    """The rule is conditional: a partial survey cannot claim a non-detection.

    Pinning "no gene is ever fully surveyed" instead would enshrine an
    arithmetic artifact -- that is what counting HPA's rare special-study
    labels into the denominator produced, and it made a detection the only
    escape from ``incomplete``.
    """
    summary = cta_review.cta_evidence_summary()
    for scope in ("somatic_ihc", "brain_ihc", "heart_ihc"):
        measured = summary[f"{scope}_measured_tissues"]
        expected = summary[f"{scope}_expected_tissues"]
        status = summary[f"{scope}_status"]
        # Reachable in both directions, so neither branch is vacuous.
        assert status.eq("not_detected").any(), scope
        # A clean negative is only ever claimed from a complete survey.
        clean = status.eq("not_detected")
        assert measured.loc[clean].ge(expected.loc[clean]).all(), scope
        # And a short survey is never reported as one.
        short = measured.lt(expected) & summary[f"{scope}_detected_rows"].eq(0) & measured.gt(0)
        assert status.loc[short].eq("incomplete").all(), scope


def test_denominators_count_only_routinely_surveyed_labels():
    """A scope expects what the source runs, not every label it ever used.

    HPA v23 mixes its standard panel with special-study labels measured for a
    handful of genes. Counting those would put every threshold out of reach.
    """
    ihc = cta_review._atlas_tables()["ihc"]
    routine = cta_review._routine_labels("ihc")
    genes = ihc["Gene"].nunique()
    per_label = ihc.groupby("Tissue")["Gene"].nunique()
    rare = set(per_label.index) - routine
    assert rare, "expected HPA to carry special-study labels"
    # The two populations are far apart, so the cut is not delicate.
    assert per_label[sorted(rare)].max() < 0.05 * genes
    assert per_label[sorted(routine)].min() > 0.9 * genes
    summary = cta_review.cta_evidence_summary()
    # Every denominator is reachable; an unreachable one makes a detection the
    # only way out of "incomplete" and inverts the gap field.
    for scope in ("somatic_ihc", "brain_ihc", "heart_ihc", "lung_ihc"):
        measured = summary[f"{scope}_measured_tissues"]
        expected = summary[f"{scope}_expected_tissues"]
        assert measured.max() >= expected.iloc[0], scope


def test_gap_freedom_is_not_predicted_by_protein_detection():
    """Regression: an empty gap field must not mark the worst candidates.

    When a denominator was unreachable, the only escape from an _incomplete IHC
    status was a detection there, so gap-freedom became a precise marker for
    "protein detected in brain". Set inequality would be too weak a check --
    it passes on a one-gene difference -- so compare the detection rate among
    gap-free rows against the base rate instead.
    """
    summary = cta_review.cta_evidence_summary()
    gap_free = summary["atlas_evidence_gaps"].eq("")
    assert gap_free.sum() > 50
    detected = summary["atlas_warning_codes"].str.contains("somatic_protein_detected")
    base_rate = detected.mean()
    among_gap_free = detected.loc[gap_free].mean()
    # Gap-freedom must not concentrate protein detections; under the inversion
    # this ratio was 1/base_rate, every gap-free row being a detection.
    assert among_gap_free < 3 * base_rate
    assert among_gap_free < 0.5
    # And most gap-free rows are ordinary clean candidates.
    assert (~detected.loc[gap_free]).mean() > 0.5


def test_atlas_coverage_rows_reconcile_with_the_requested_groups():
    coverage = cta_review.cta_atlas_coverage()
    # One row per requested tissue, so the levels add up and a region covered
    # by a single substructure is not tallied as fully covered.
    assert set(coverage["modality"]) == {"rna", "ihc"}
    for (group, modality), rows in coverage.groupby(["safety_group", "modality"]):
        assert len(rows) == len(cta_review.SAFETY_TISSUE_GROUPS[group]), (group, modality)
        assert set(rows["coverage_level"]) <= {"complete", "partial", "unavailable"}
    brain_ihc = coverage.loc[
        coverage["safety_group"].eq("brain") & coverage["modality"].eq("ihc")
    ].set_index("requested_tissue")
    assert brain_ihc.loc["basal ganglia", "coverage_level"] == "partial"
    assert brain_ihc.loc["basal ganglia", "source_tissues"] == "caudate"
    assert brain_ihc.loc["midbrain", "source_tissues"] == "dorsal raphe;substantia nigra"
    assert brain_ihc.loc["spinal cord", "coverage_level"] == "unavailable"
    assert brain_ihc.loc["spinal cord", "source_tissues"] == ""
    # The vocabulary matches the summary's column prefixes, so a join works.
    summary = cta_review.cta_evidence_summary()
    assert "brain_rna_mapping_coverage" in summary.columns
    assert coverage["modality"].isin(["rna", "ihc"]).all()


def test_coverage_accessor_and_denominator_agree_on_what_was_surveyed():
    """The two public surfaces must not give opposite answers about a tissue.

    A label can be mapped and still sit outside the panel the source runs for
    most genes. Reporting it as covered here while excluding it from the
    summary's denominator had the accessor advertising 8 of 14 brain regions
    while the status it qualifies rested on 4 labels.
    """
    coverage = cta_review.cta_atlas_coverage()
    summary = cta_review.cta_evidence_summary()
    # RNA scopes count measurements (_expected_rows), IHC counts tissues.
    for modality, suffix in (("rna", "expected_rows"), ("ihc", "expected_tissues")):
        rows = coverage.loc[coverage["modality"].eq(modality)]
        for group in cta_review.SAFETY_TISSUE_GROUPS:
            group_rows = rows.loc[rows["safety_group"].eq(group)]
            surveyed = {
                label
                for joined in group_rows["surveyed_source_tissues"]
                for label in joined.split(";")
                if label
            }
            expected = summary[f"{group}_{modality}_{suffix}"].iloc[0]
            assert len(surveyed) == expected, (group, modality)
    # A mapped-but-not-routine label is reported as such, not as covered.
    brain_ihc = coverage.loc[
        coverage["safety_group"].eq("brain") & coverage["modality"].eq("ihc")
    ].set_index("requested_tissue")
    assert brain_ihc.loc["retina", "coverage_level"] == "complete"
    assert brain_ihc.loc["retina", "survey_state"] == "mapped_not_surveyed"
    assert brain_ihc.loc["retina", "surveyed_coverage_level"] == "unavailable"
    assert brain_ihc.loc["cerebellum", "survey_state"] == "surveyed"
    assert brain_ihc.loc["cerebellum", "surveyed_coverage_level"] == "complete"
    # A region with no source label at all is a different case, not the same
    # one: a single boolean would have put thalamus and retina together.
    assert brain_ihc.loc["thalamus", "survey_state"] == "not_mapped"


def test_unsurveyed_labels_are_reported_once_per_release():
    """A fact identical in every row belongs to the release, not the gene.

    Naming these per group and modality meant ten columns of which nine were
    permanently empty and one held a constant string in all 403 rows -- the
    criterion this module already applied when it split the mapping limits out
    of the per-gene gap field.
    """
    summary = cta_review.cta_evidence_summary()
    assert not [c for c in summary.columns if c.endswith("_unsurveyed_tissues")]
    limits = summary["atlas_coverage_limits"]
    assert limits.str.contains("brain_ihc_unsurveyed_labels").all()
    # Only the scope that has unsurveyed labels says so.
    assert not limits.str.contains("brain_rna_unsurveyed_labels").any()
    assert not limits.str.contains("heart_ihc_unsurveyed_labels").any()
    # The labels themselves live on the accessor, one row per tissue.
    coverage = cta_review.cta_atlas_coverage()
    brain_ihc = coverage.loc[coverage["safety_group"].eq("brain") & coverage["modality"].eq("ihc")]
    # One filter, and exact: the conflated form returned 10 regions here, six
    # of which HPA never maps, so a subset assertion would have passed on it.
    mapped_not_surveyed = brain_ihc.loc[brain_ihc["survey_state"].eq("mapped_not_surveyed")]
    assert set(mapped_not_surveyed["requested_tissue"]) == {
        "choroid plexus",
        "hypothalamus",
        "midbrain",
        "retina",
    }
    # The labels behind them are recoverable, including midbrain's two nuclei,
    # which carry no surveyed label of their own.
    labels = {t for joined in mapped_not_surveyed["source_tissues"] for t in joined.split(";")}
    assert labels == {
        "choroid plexus",
        "dorsal raphe",
        "hypothalamus",
        "retina",
        "substantia nigra",
    }


def test_coverage_accessor_states_both_region_counts():
    # The aggregate verdict is "partial" on either basis for brain, so the
    # counts are what distinguish 8 mapped regions from 4 surveyed ones.
    coverage = cta_review.cta_atlas_coverage()
    brain_ihc = coverage.loc[
        coverage["safety_group"].eq("brain") & coverage["modality"].eq("ihc")
    ].iloc[0]
    assert brain_ihc["group_requested_regions"] == 14
    assert brain_ihc["group_mapped_regions"] == 8
    assert brain_ihc["group_surveyed_regions"] == 4
    assert brain_ihc["group_surveyed_coverage_state"] == "partial"
    heart = coverage.loc[coverage["safety_group"].eq("heart")]
    assert heart["group_mapped_regions"].eq(1).all()
    assert heart["group_surveyed_regions"].eq(1).all()
    assert heart["group_surveyed_coverage_state"].eq("complete").all()
    # Surveyed can never exceed mapped, which can never exceed requested.
    assert (coverage["group_surveyed_regions"] <= coverage["group_mapped_regions"]).all()
    assert (coverage["group_mapped_regions"] <= coverage["group_requested_regions"]).all()


def test_cardiomyocyte_denominator_uses_the_same_routine_intersection():
    # The one scope that was left out of the correction: if heart muscle ever
    # left the routine panel this would silently become unreachable.
    summary = cta_review.cta_evidence_summary()
    routine = cta_review._routine_labels("ihc")
    heart = set(cta_review._resolve("ihc")["heart"].source_tissues)
    assert summary["cardiomyocyte_ihc_expected_tissues"].eq(len(heart & routine)).all()
    assert summary["cardiomyocyte_ihc_status"].eq("not_detected").any()


def test_surveyed_level_downgrades_a_partly_routine_mapping():
    """A region keeps its mapped level only if all its labels are routine.

    Not reachable from v23 -- no mapping there mixes routine and special-study
    labels -- but the rule is what keeps the accessor from overstating a region
    that HPA covers by one routine label and one it rarely runs.
    """

    class _Mapping:
        def __init__(self, level, tissues):
            self.coverage_level = level
            self.source_tissues = tissues

    exact = _Mapping("complete", ("cerebellum", "dorsal raphe"))
    assert cta_review._surveyed_level(exact, ("cerebellum",)) == "partial"
    assert cta_review._surveyed_level(exact, ()) == "unavailable"
    assert cta_review._surveyed_level(exact, ("cerebellum", "dorsal raphe")) == "complete"
    substructure = _Mapping("partial", ("caudate",))
    assert cta_review._surveyed_level(substructure, ("caudate",)) == "partial"


def test_aggregate_level_collapses_per_tissue_coverage():
    """One verdict per group, mirroring how hpa aggregates mapped levels.

    The all-unavailable case is not reachable from v23 -- every safety group
    has at least one routinely surveyed label -- but a release that dropped a
    group's whole panel must report it as unavailable rather than partial.
    """
    assert cta_review._aggregate_level(["complete", "complete"]) == "complete"
    assert cta_review._aggregate_level(["unavailable", "unavailable"]) == "unavailable"
    assert cta_review._aggregate_level(["complete", "unavailable"]) == "partial"
    assert cta_review._aggregate_level(["partial"]) == "partial"
    assert cta_review._aggregate_level(["complete", "partial"]) == "partial"


def test_out_of_scope_labels_cannot_complete_a_survey():
    """A negative must be earned in the labels the scope is counted over.

    Counting a measurement in a rarely-run special-study label toward a
    routine-only denominator let a gene scored in none of the four routine
    brain regions report a clean brain non-detection.
    """
    resolution = cta_review._resolve("ihc")["brain"]
    routine = cta_review._routine_labels("ihc")
    outside = sorted(set(resolution.source_tissues) - routine)
    assert len(outside) >= 2, "expected v23 to map brain labels it rarely stains"
    tables = dict(cta_review._atlas_tables())
    scored_outside = pd.DataFrame(
        [["GX", "X", tissue, "cells", "Not detected", "Enhanced"] for tissue in outside],
        columns=["Gene", "Gene name", "Tissue", "Cell type", "Level", "Reliability"],
    )
    tables["ihc"] = pd.concat([tables["ihc"], scored_outside], ignore_index=True)
    row = cta_review._synthesize_atlas(pd.DataFrame({"Ensembl_Gene_ID": ["GX"]}), tables).iloc[0]
    assert row["brain_ihc_measured_tissues"] == 0
    assert row["brain_ihc_expected_tissues"] == len(set(resolution.source_tissues) & routine)
    assert row["brain_ihc_status"] == "incomplete"
    assert "brain_ihc_incomplete" in row["atlas_evidence_gaps"].split(";")


def test_brain_coverage_differs_by_modality():
    """The two modalities miss different regions; neither stands for the other.

    Pins the claim the docs and the accessor docstring make, which previously
    said a brain result of any kind spoke for neither spinal cord nor thalamus
    — false for RNA, which maps and routinely surveys spinal cord.
    """
    coverage = cta_review.cta_atlas_coverage()
    brain = coverage.loc[coverage["safety_group"].eq("brain")]
    missed = {
        modality: set(rows.loc[rows["survey_state"].ne("surveyed"), "requested_tissue"])
        for modality, rows in brain.groupby("modality")
    }
    shared = missed["rna"] & missed["ihc"]
    assert shared == {"thalamus", "medulla oblongata", "pons", "white matter"}
    assert "spinal cord" in missed["ihc"]
    assert "spinal cord" not in missed["rna"]
    assert missed["rna"] < missed["ihc"]


def test_a_blank_finding_names_its_row_instead_of_crashing_every_gene():
    # Joining a NaN raises TypeError deep in a groupby agg, taking down the
    # whole summary rather than pointing at the row that needs fixing.
    reviewed = cta_review.cta_reviewed_evidence()
    reviewed.loc[reviewed.index[0], "finding"] = None
    with pytest.raises(ValueError, match="no finding"):
        cta_review._review_notes(reviewed, "isoform")


def test_reviewed_evidence_findings_are_present():
    # The one column whose emptiness is fatal was the one left unchecked.
    reviewed = cta_review.cta_reviewed_evidence()
    assert reviewed["finding"].notna().all()
    assert reviewed["finding"].str.strip().ne("").all()


def test_repeated_rna_label_cannot_replace_a_missing_tissue():
    rows = _tables(
        bulk_rna=[
            ["G1", "A", "lung", 0.0],
            ["G1", "A", "lung", 0.0],
            ["G1", "A", "liver", float("inf")],
        ]
    )["bulk_rna"]
    stats = cta_review._rna_stats(rows, frozenset({"lung", "liver"}), "Tissue")
    assert stats["status"] == "incomplete"
    assert stats["measured_rows"] == 1
    assert stats["max_ntpm"] == 0


def test_positive_ihc_does_not_hide_missing_tissues():
    row = _synthesize(
        ihc=[
            ["G1", "A", "lung", "cells", "High", "Enhanced"],
            ["G2", "B", "lung", "cells", "Not detected", "Enhanced"],
            ["G2", "B", "liver", "cells", "Not detected", "Enhanced"],
        ]
    )
    assert row["somatic_ihc_status"] == "detected"
    assert "somatic_ihc_incomplete" in row["atlas_evidence_gaps"].split(";")


def test_atlas_loader_preserves_unscored_tissue_rows(monkeypatch):
    tables = _tables(
        ihc=[
            ["G1", "A", "liver", "hepatocytes", "Not detected", "Enhanced"],
            ["G1", "A", "liver", "bile duct cells", None, "Enhanced"],
            ["G2", "B", None, None, None, "Uncertain"],
        ]
    )
    monkeypatch.setattr(cta_review.hpa, "hpa_normal_tissue", lambda version: tables["ihc"])
    monkeypatch.setattr(cta_review.hpa, "_read_hpa", lambda *args: tables["bulk_rna"])
    ihc = cta_review._atlas_tables.__wrapped__()["ihc"]
    assert len(ihc) == 2
    assert cta_review._ihc_stats(ihc, frozenset({"liver"}))["status"] == "incomplete"


@pytest.mark.parametrize("finding", [None, "", "   "])
def test_review_notes_reject_all_blank_findings(finding):
    reviewed = pd.DataFrame(
        {"Ensembl_Gene_ID": ["G1"], "modality": ["isoform"], "finding": [finding]}
    )
    with pytest.raises(ValueError, match="G1"):
        cta_review._review_notes(reviewed, "isoform")


def test_summary_rejects_missing_atlas_candidates(monkeypatch):
    monkeypatch.setattr(
        cta_review, "_synthesize_atlas", lambda *args: pd.DataFrame({"Ensembl_Gene_ID": []})
    )
    with pytest.raises(ValueError, match="Atlas candidate mismatch: missing="):
        cta_review._evidence_summary_frame.__wrapped__()

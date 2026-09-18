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
    # HPA v23 covers 8 of the 14 requested brain regions by IHC and 10 by RNA,
    # so neither modality speaks for spinal cord or thalamus.
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
    assert ihc["ihc_level"].notna().all()
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

    Scanned against the synthesis output rather than the merged summary, so
    the test and the gaps loop share one scope. Scanning the merged frame
    would sweep in the ``*_review_status`` columns added afterwards and need a
    hand-maintained exclusion -- the fragility the loop was rewritten to drop.
    """
    universe = cta_review._universe()
    summary = cta_review._synthesize_atlas(universe, cta_review._atlas_tables())
    statuses = [c for c in summary.columns if c.endswith("_status")]
    assert len(statuses) > 10
    for column in statuses:
        prefix = column[: -len("_status")]
        flagged = summary[column].isin(["unavailable", "incomplete"])
        if not flagged.any():
            continue
        for status in ("unavailable", "incomplete"):
            expected = summary[column].eq(status)
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
    routine = cta_review._routine_labels(ihc, "Tissue")
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


def test_empty_gaps_does_not_select_the_most_concerning_genes():
    # Regression: when a denominator was unreachable, the only escape from an
    # _incomplete IHC status was a detection there, so an empty gap field
    # became a perfect marker for "protein detected in brain".
    summary = cta_review.cta_evidence_summary()
    gap_free = summary.loc[summary["atlas_evidence_gaps"].eq("")]
    assert len(gap_free) > 50
    detected = set(summary.loc[summary["brain_ihc_status"].eq("detected"), "Symbol"])
    assert set(gap_free["Symbol"]) != detected
    # Most gap-free rows carry no protein-detection warning at all.
    flagged = gap_free["atlas_warning_codes"].str.contains("somatic_protein_detected")
    assert flagged.mean() < 0.5


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

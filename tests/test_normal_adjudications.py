# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Tests for normal-library (#712) and rare-cell attribution (#711) adjudications."""

from __future__ import annotations

import pytest

import oncoref
from oncoref import normal_adjudications, target_selection


def test_normal_library_peak_adjudications_covers_all_eight_712_candidates():
    df = oncoref.normal_library_peak_adjudications()
    assert len(df) == 8
    assert set(df["protein_group_label"]) == {
        "C4orf17",
        "TTLL8",
        "GOLGA6L2",
        "CCDC87",
        "RNASE13",
        "GAB4",
        "PSG7",
        "GSTT4",
    }
    assert set(df["adjudication_disposition"]) <= set(
        normal_adjudications.ADJUDICATION_DISPOSITIONS_712
    )

    # All 6 HPA libraries pass published carryover & testis-content QC
    assert not df["library_carryover_p"].any()
    assert not df["library_carryover_e"].any()
    assert (df["library_testis_content"] == 0.0).all()

    by_label = df.set_index("protein_group_label")
    assert by_label.loc["TTLL8", "active_classification"] == "CTA"
    assert by_label.loc["GOLGA6L2", "active_classification"] == "CTA"
    assert by_label.loc["GAB4", "adjudication_disposition"] == "superseded_by_independent_gtex_p95"
    assert (
        by_label.loc["CCDC87", "adjudication_disposition"]
        == "retained_multi_tissue_cross_atlas_concordant"
    )
    assert bool(by_label.loc["PSG7", "paralog_ambiguous"]) is True
    assert bool(by_label.loc["GSTT4", "paralog_ambiguous"]) is True

    restored = oncoref.normal_library_peak_adjudications(
        disposition="restored_to_cta_under_active_p95"
    )
    assert set(restored["protein_group_label"]) == {"TTLL8", "GOLGA6L2"}

    with pytest.raises(ValueError, match="Unknown disposition"):
        oncoref.normal_library_peak_adjudications(disposition="invalid")


def test_hpa_normal_statistic_sensitivity_preserves_unsafe_precedence():
    sens = oncoref.hpa_normal_statistic_sensitivity().set_index("scenario_id")
    ts = target_selection.results()

    active = sens.loc["active_policy_gtex_p95_plus_hpa_biopsy_p95"]
    assert bool(active["is_active_policy"]) is True
    assert int(active["total_proteins"]) == len(ts)
    assert int(active["cta_count"]) == int(ts["is_cta"].sum())
    assert int(active["unsafe_count"]) == int(ts["category"].eq("Unsafe").sum())
    assert int(active["delta_cta_vs_active"]) == 0

    # Every counterfactual preserves the exact Unsafe count
    assert (sens["unsafe_count"] == int(active["unsafe_count"])).all()
    assert (sens.loc[~sens["is_active_policy"], "delta_cta_vs_active"] >= 0).all()


def test_rare_cell_candidate_adjudications_and_unsafe_propagation():
    reviews = oncoref.rare_cell_attribution_reviews()
    assert reviews["schema"] == "oncoref.rare-cell-attribution-reviews"
    assert len(reviews["reviews"]) == 15

    df = oncoref.rare_cell_candidate_adjudications()
    assert len(df) == 15
    assert set(df["protein_group_label"]) == {
        "CAGE1",
        "RAD21L1",
        "PPP3R2",
        "RGSL1",
        "DAZL",
        "OOSP4B",
        "TAS2R1",
        "PSG8",
        "ZNF280A",
        "POTEC",
        "VCX3B",
        "PSG3",
        "TUBB8",
        "DCAF8L2",
        "FBXO47",
    }

    by_label = df.set_index("protein_group_label")
    # RGSL1 and OOSP4B had host-overlap confounders rejected -> promoted to Unsafe & forbidden background
    for unsafe_lbl in ("RGSL1", "OOSP4B"):
        row = by_label.loc[unsafe_lbl]
        assert row["attribution_status"] == "confounder_rejected"
        assert row["active_classification"] == "Unsafe"
        assert bool(row["in_unsafe_forbidden_background"]) is True
        assert row["sequence_bundle_role"] == "exclude"
        assert row["resolution_requirement"] == "resolved_unsafe_in_forbidden_background"

    # Remaining 13 stay held in Other:Review with explicit evidence requirements
    held_ids = oncoref.unresolved_tenfold_rare_cell_entity_ids()
    assert len(held_ids) == 13
    assert set(df.loc[df["active_classification"].eq("Other:Review"), "entity_id"]) == held_ids

    spliced_needed = oncoref.rare_cell_candidate_adjudications(
        resolution_requirement="spliced_exon_junction_or_polya_bam_attribution"
    )
    assert {"CAGE1", "RAD21L1", "PPP3R2", "TAS2R1"} == set(spliced_needed["protein_group_label"])

    single_study = oncoref.rare_cell_candidate_adjudications(
        attribution_status="single_study_unconfirmed"
    )
    assert {"DCAF8L2", "FBXO47"} == set(single_study["protein_group_label"])

    second_study_needed = oncoref.rare_cell_candidate_adjudications(
        resolution_requirement="independent_second_donor_study_or_targeted_ish_ihc"
    )
    assert {"DAZL", "DCAF8L2", "FBXO47"} == set(second_study_needed["protein_group_label"])

    with pytest.raises(ValueError, match="Unknown attribution_status"):
        oncoref.rare_cell_candidate_adjudications(attribution_status="bad")
    with pytest.raises(ValueError, match="Unknown resolution_requirement"):
        oncoref.rare_cell_candidate_adjudications(resolution_requirement="bad")

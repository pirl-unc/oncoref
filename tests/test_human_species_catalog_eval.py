"""Tests for unified human SpeciesScreeningSpec evaluation and v1.49 layer integration (#755)."""

from __future__ import annotations

import pandas as pd

import oncoref
from oncoref import human, human_screening_inputs, target_selection


def test_human_evaluate_species_catalog_achieves_exact_19893_baseline_parity():
    baseline = target_selection.results()
    res = human.evaluate_integrated_catalog(
        include_home_panels=False,
        include_adjudications=False,
        include_737_cohorts=False,
    )
    cat = res.catalog

    assert len(cat) == len(baseline) == 19893
    assert (cat["taxon_id"] == 9606).all()
    assert (cat["species"] == "human").all()
    assert cat["classification"].equals(baseline["classification"])
    assert cat["category"].equals(baseline["category"])
    assert cat["is_cta"].equals(baseline["is_cta"].astype(bool))
    assert cat["is_cancer_antigen"].equals(baseline["is_cancer_antigen"].astype(bool))
    assert cat["is_reproductive_antigen"].equals(baseline["is_reproductive_antigen"].astype(bool))

    # Exact sex-restriction subset parity across all 455 human RA/CTA proteins
    expected_sex = human.sex_restricted_antigens().set_index("entity_id")["sex_restriction_code"]
    got_sex = res.sex_restricted_antigens().set_index("entity_id")["sex_restriction_code"]
    pd.testing.assert_series_equal(got_sex, expected_sex)

    # Pinned default human catalog remains untouched
    assert human.target_selection().attrs["policy_id"] == target_selection.POLICY_ID
    assert int(human.target_selection()["is_cta"].sum()) == 327


def test_upgraded_737_tumor_observations_and_integrated_policy_preview():
    obs_737 = human.upgraded_737_tumor_observations()
    assert set(obs_737["cancer_code"]) == {"THYMCA", "BCC", "cSCC"}
    assert (obs_737["value"] >= 1.0).all()

    summary = human.next_policy_summary()
    preview = human.next_policy_preview()
    integrated_res = human.evaluate_integrated_catalog(
        include_home_panels=True,
        include_adjudications=True,
        include_737_cohorts=True,
        include_multiplex_bleed_qc=True,
        include_nuclear_spillover_qc=True,
    )
    recomputed, recomputed_summary = human_screening_inputs.compute_next_policy_transitions(
        integrated_result=integrated_res
    )

    assert len(preview) == len(recomputed) == summary["changed_proteins"]
    assert summary["total_proteins"] == 19893
    assert summary["integrated_categories"] == recomputed_summary["integrated_categories"]
    pd.testing.assert_frame_equal(preview, recomputed)

    by_label = preview.set_index("protein_group_label")
    # #735 ceiling violations become Unsafe (SP7 bone >150 TPM, PANX3 bone ceiling, GFY olfactory >300 TPM)
    for unsafe_lbl in ("SP7", "PANX3", "BPIFB3", "GFY", "GLYATL3"):
        assert by_label.loc[unsafe_lbl, "classification_integrated"] == "Unsafe"
        assert (
            "#735:home_expressed_above_tier_ceiling"
            in by_label.loc[unsafe_lbl, "transition_driver"]
        )

    # #735 cleared home panels admit TAS1R2 / KRTAP21-2 / KRTAP9-7 to Cancer Antigen only
    # and promote reproductive-positive OR10J1 from Reproductive Antigen only to CTA
    for ca_lbl in ("TAS1R2", "KRTAP21-2", "KRTAP9-7"):
        assert by_label.loc[ca_lbl, "classification_integrated"] == "Cancer Antigen only"
    assert by_label.loc["OR10J1", "classification_integrated"] == "CTA"

    # #737 linear-TPM cohorts admit BPY2/B/C (cSCC) and promote DEFB107A/B (BCC) to CTA
    assert by_label.loc["BPY2/B/C", "classification_integrated"] == "Cancer Antigen only"
    assert "#737:cSCC" in by_label.loc["BPY2/B/C", "transition_driver"]
    assert by_label.loc["DEFB107A/B", "classification_integrated"] == "CTA"
    assert "#737:BCC" in by_label.loc["DEFB107A/B", "transition_driver"]

    # #711 confounder-rejected DAZL transitions from Other:Review to Unsafe
    assert by_label.loc["DAZL", "classification_integrated"] == "Unsafe"
    assert "#711:confounder_rejected" in by_label.loc["DAZL", "transition_driver"]

    # #735 insufficient 10x fold over home moves OR11A1 and PABPC3 to Other:Expression
    assert by_label.loc["OR11A1", "classification_integrated"] == "Other:Expression"
    assert by_label.loc["PABPC3", "classification_integrated"] == "Other:Expression"

    # #735 incomplete anatomical sub-niche (DSPP soft dental pulp) stays conservatively unchanged
    assert "DSPP" not in by_label.index

    # #711 and #712 provenance annotations are populated on the integrated catalog
    integrated_cat = integrated_res.catalog.set_index("protein_group_label")
    assert integrated_cat.loc["CAGE1", "classification"] == "Other:Review"
    assert (
        integrated_cat.loc["CAGE1", "rare_cell_711_requirement"]
        == "spliced_exon_junction_or_polya_bam_attribution"
    )
    assert (
        integrated_cat.loc["TTLL8", "normal_library_712_disposition"]
        == "restored_to_cta_under_active_p95"
    )
    assert hasattr(oncoref, "human_screening_inputs")

"""Tests for unified human SpeciesScreeningSpec evaluation and v1.49 layer integration (#755)."""

from __future__ import annotations

import pandas as pd

import oncoref
from oncoref import human, human_screening_inputs, target_selection


def test_human_evaluate_species_catalog_achieves_exact_19893_baseline_parity():
    promoted = target_selection.results()
    res_default = human.evaluate_integrated_catalog()
    cat_default = res_default.catalog

    assert len(cat_default) == len(promoted) == 19893
    assert (cat_default["taxon_id"] == 9606).all()
    assert (cat_default["species"] == "human").all()
    assert cat_default["classification"].equals(promoted["classification"])
    assert cat_default["category"].equals(promoted["category"])
    assert cat_default["is_cta"].equals(promoted["is_cta"].astype(bool))
    assert cat_default["is_cancer_antigen"].equals(promoted["is_cancer_antigen"].astype(bool))
    assert cat_default["is_reproductive_antigen"].equals(
        promoted["is_reproductive_antigen"].astype(bool)
    )
    assert int(promoted["is_cta"].sum()) == 322

    # Exact sex-restriction subset parity across all 468 human RA/CTA proteins
    expected_sex = human.sex_restricted_antigens().set_index("entity_id")["sex_restriction_code"]
    got_sex = res_default.sex_restricted_antigens().set_index("entity_id")["sex_restriction_code"]
    pd.testing.assert_series_equal(got_sex.sort_index(), expected_sex.sort_index())

    # Orthogonal single-rule ablation turning off chemosensory GPCR guard restores 325 CTA (2026-11-26.1)
    res_no_gpcr = human.evaluate_integrated_catalog(include_chemosensory_gpcr_guard=False)
    assert int(res_no_gpcr.catalog["is_cta"].sum()) == 325

    # Disabling all integrated layers reproduces the frozen 2026-11-20.2 baseline (327 CTA)
    baseline_20261120 = target_selection._baseline_2026_11_20_frame()
    res_base = human.evaluate_integrated_catalog(
        include_home_panels=False,
        include_adjudications=False,
        include_737_cohorts=False,
        include_multiplex_bleed_qc=False,
        include_nuclear_spillover_qc=False,
        include_clinical_clearances=False,
        include_repeat_family_equivalence=False,
        include_assay_boundary_fixes=False,
    )
    assert res_base.catalog["classification"].equals(baseline_20261120["classification"])
    assert int(res_base.catalog["is_cta"].sum()) == 327

    # Pinned default human catalog matches promoted policy_id (322 CTA)
    assert human.target_selection().attrs["policy_id"] == target_selection.POLICY_ID
    assert int(human.target_selection()["is_cta"].sum()) == 322
    assert set(promoted.loc[promoted["is_cta"], "tumor_activation_mode"]) == {
        "ectopic",
        "ectopic_and_lineage_retained",
        "lineage_retained",
    }
    by_prom = promoted.set_index("protein_group_label")
    assert by_prom.loc["SEMG2", "tumor_activation_mode"] == "lineage_retained"
    assert by_prom.loc["CSN2", "tumor_activation_mode"] == "ectopic_and_lineage_retained"
    assert by_prom.loc["PRAME", "tumor_activation_mode"] == "ectopic_and_lineage_retained"
    assert by_prom.loc["NY-ESO-1", "tumor_activation_mode"] == "ectopic_and_lineage_retained"
    assert cat_default["tumor_activation_mode"].equals(promoted["tumor_activation_mode"])
    assert hasattr(human, "tumor_cohorts") and hasattr(human, "cell_lines")

    # Orthogonal #768 TGCT Note-R single-rule ablation keeps TGCT-witnessed RAs as CTA
    res_no_768 = human.evaluate_integrated_catalog(include_assay_boundary_fixes=False).catalog
    by_no_768 = res_no_768.set_index("protein_group_label")
    for tgct_cta in ("BOD1L2", "CCDC70", "PGK2", "PRR30", "ZNRF4", "TSSK2"):
        assert by_no_768.loc[tgct_cta, "classification"] == "CTA"


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
        include_clinical_clearances=True,
        include_repeat_family_equivalence=True,
        include_assay_boundary_fixes=True,
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

    # #735 cleared home panels admit KRTAP21-2 / KRTAP20-3 to Cancer Antigen only and promote
    # reproductive-positive KRTAP20-4 to CTA, while monogenic chemosensory GPCRs (OR*, TAS1R*,
    # TAS2R*, VN1R*) stay out of CA/CTA (matching 2026-11-20 Other:Unscreened tissue / RA only,
    # with OR52D1 transitioning Other:Unscreened tissue -> Reproductive Antigen only via GSE107746)
    for ca_lbl in ("KRTAP21-2", "KRTAP20-3"):
        assert by_label.loc[ca_lbl, "classification_integrated"] == "Cancer Antigen only"
    for unchanged_vs_20261120 in ("TAS1R2", "OR2B6", "OR51E2", "OR10J1", "OR1N1", "OR4M2"):
        assert unchanged_vs_20261120 not in by_label.index
    assert by_label.loc["OR52D1", "classification_integrated"] == "Reproductive Antigen only"
    assert (
        "chemosensory_gpcr:unscreened_monogenic_sensory_neuron"
        in by_label.loc["OR52D1", "transition_driver"]
    )
    assert by_label.loc["KRTAP20-4", "classification_integrated"] == "CTA"
    assert by_label.loc["KRTAP9-7", "classification_integrated"] == "Other:Expression"
    assert by_label.loc["CAPZA3", "classification_integrated"] == "Reproductive Antigen only"
    assert (
        "library_qc:tumor_gdna_background_excluded" in by_label.loc["CAPZA3", "transition_driver"]
    )

    # #737 linear-TPM cohorts admit BPY2/B/C (cSCC) and promote DEFB107A/B (BCC) to CTA
    assert by_label.loc["BPY2/B/C", "classification_integrated"] == "Cancer Antigen only"
    assert "#737:cSCC" in by_label.loc["BPY2/B/C", "transition_driver"]
    assert by_label.loc["DEFB107A/B", "classification_integrated"] == "CTA"
    assert "#737:BCC" in by_label.loc["DEFB107A/B", "transition_driver"]

    # #779 droplet ambient-soup & <=100 kb SCPP cis-spillover clears CSN2 from Unsafe to CTA
    assert by_label.loc["CSN2", "classification_integrated"] == "CTA"
    assert "#768:sc_soup_and_cis_spillover_cleared" in by_label.loc["CSN2", "transition_driver"]

    # #711 confounder-rejected DAZL transitions from Other:Review to Unsafe
    assert by_label.loc["DAZL", "classification_integrated"] == "Unsafe"
    assert "#711:confounder_rejected" in by_label.loc["DAZL", "transition_driver"]

    # #735 insufficient 10x fold over home moves OR11A1 and PABPC3 to Other:Expression,
    # while high-fold chemosensory GPCRs sit in Other:Unscreened tissue / RA only on integrated_cat
    assert by_label.loc["OR11A1", "classification_integrated"] == "Other:Expression"
    assert by_label.loc["PABPC3", "classification_integrated"] == "Other:Expression"

    # #735 incomplete anatomical sub-niche (DSPP soft dental pulp) stays conservatively unchanged
    assert "DSPP" not in by_label.index

    # #711, #712, and chemosensory GPCR guard annotations are populated on the integrated catalog
    integrated_cat = integrated_res.catalog.set_index("protein_group_label")
    for gpcr_unscreened in ("TAS1R2", "OR7C2", "OR2D2"):
        assert integrated_cat.loc[gpcr_unscreened, "classification"] == "Other:Unscreened tissue"
    for gpcr_ra_only in ("OR1N1", "OR4M2", "OR10J1"):
        assert integrated_cat.loc[gpcr_ra_only, "classification"] == "Reproductive Antigen only"
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

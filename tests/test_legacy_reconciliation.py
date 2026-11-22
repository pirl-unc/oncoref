# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

from __future__ import annotations

import warnings

import pandas as pd

import oncoref
from oncoref import cta, legacy, normal_expression_categories, target_selection


def test_top_level_legacy_cta_accessors_warn_and_point_to_target_selection():
    legacy._WARNED_LEGACY_APIS.clear()
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", DeprecationWarning)
            fn = oncoref.cta_gene_ids
            fn_again = oncoref.cta_gene_ids
            upper_fn = oncoref.CTA_gene_ids

        assert fn is cta.cta_gene_ids
        assert fn_again is cta.cta_gene_ids
        assert upper_fn is cta.CTA_gene_ids
        emitted = [w for w in caught if w.category is DeprecationWarning]
        assert len(emitted) == 2
        msg = str(emitted[0].message)
        assert "oncoref.cta_gene_ids" in msg
        assert "antigen_sets.gene_ids('CTA')" in msg
        assert target_selection.POLICY_ID in msg
    finally:
        legacy._WARNED_LEGACY_APIS.clear()


def test_explicit_cta_submodule_remains_warning_free_for_historical_provenance():
    legacy._WARNED_LEGACY_APIS.clear()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", DeprecationWarning)
        ids = cta.cta_gene_ids()
        assert len(ids) > 500
    emitted = [w for w in caught if w.category is DeprecationWarning]
    assert not emitted


def test_legacy_api_migrations_covers_all_cta_generations():
    table = oncoref.legacy_api_migrations()
    assert not table.empty
    assert set(table.columns) == {
        "legacy_surface",
        "generation",
        "current_status",
        "replacement_surface",
        "active_policy_id",
        "notes",
    }
    assert oncoref.legacy_api_migration("cta_gene_ids") is not None
    assert oncoref.legacy_api_migration("oncoref.target_sets")["active_policy_id"] == (
        target_selection.POLICY_ID
    )
    assert oncoref.legacy_api_migration("oncoref.normal_expression_categories") is not None
    assert oncoref.legacy_api_migration("nonexistent_module") is None


def test_reconcile_with_target_selection_preserves_unsafe_precedence_and_disposable_cap():
    catalog = target_selection.results()
    reconciled = normal_expression_categories.reconcile_with_target_selection(catalog)
    assert len(reconciled) == len(catalog)
    assert reconciled.attrs["policy_id"] == target_selection.POLICY_ID
    assert set(reconciled["reconciled_normal_profile"]) <= set(
        normal_expression_categories.RECONCILED_NORMAL_PROFILES
    )

    # Every catalog Unsafe row stays unsafe_normal regardless of Disposable values.
    unsafe_mask = catalog["category"].eq("Unsafe")
    assert (reconciled.loc[unsafe_mask, "reconciled_normal_profile"] == "unsafe_normal").all()

    # Review / Sequence / Unscreened tissue never collapse into below_threshold_silent.
    unscreened_mask = catalog["classification"].isin(
        ["Other:Review", "Other:Sequence", "Other:Unscreened tissue"]
    )
    assert (
        reconciled.loc[unscreened_mask, "reconciled_normal_profile"] == "review_or_unscreened"
    ).all()

    # Every Other row remains in the default sequence-subtraction background.
    other_mask = catalog["category"].eq("Other")
    assert reconciled.loc[other_mask, "in_default_subtraction_background"].all()

    disposable = normal_expression_categories.adult_disposable_only_proteins(catalog)
    assert (disposable["target_selection_classification"] == "Other:Expression").all()
    assert (
        (disposable["disposable_max_tpm"] >= target_selection.DETECTION_FLOOR)
        | (disposable["disposable_max_ntpm"] >= target_selection.DETECTION_FLOOR)
    ).all()
    assert (disposable["disposable_max_tpm"] <= target_selection.DISPOSABLE_CAP).all()
    assert (disposable["disposable_max_ntpm"] <= target_selection.DISPOSABLE_CAP).all()
    assert (disposable["normal_reference_tpm"] < target_selection.DETECTION_FLOOR).all()
    assert (disposable["healthy_maximum_ntpm"] < target_selection.DETECTION_FLOOR).all()

    # Synthetic check for disposable_above_cap vs unsafe_normal precedence.
    synthetic = catalog.iloc[:2].copy()
    synthetic.loc[synthetic.index[0], "category"] = "Unsafe"
    synthetic.loc[synthetic.index[0], "classification"] = "Unsafe"
    synthetic.loc[synthetic.index[0], "disposable_max_tpm"] = 10.0
    synthetic.loc[synthetic.index[1], "category"] = "Other"
    synthetic.loc[synthetic.index[1], "classification"] = "Other:Expression"
    synthetic.loc[synthetic.index[1], "is_target_candidate"] = False
    for col in (
        "normal_screen_pass",
        "required_coverage_complete",
        "identity_translation_pass",
        "evidence_resolved",
        "sequence_pass",
    ):
        synthetic.loc[synthetic.index[1], col] = True
    synthetic.loc[synthetic.index[1], "expression_bounds_pass"] = False
    synthetic.loc[synthetic.index[1], "disposable_max_tpm"] = 50.0
    synthetic.loc[synthetic.index[1], "disposable_max_ntpm"] = 0.0
    synthetic.loc[synthetic.index[1], "normal_reference_tpm"] = 0.0
    synthetic.loc[synthetic.index[1], "healthy_maximum_ntpm"] = 0.0

    syn_rec = normal_expression_categories.reconcile_with_target_selection(synthetic)
    assert list(syn_rec["reconciled_normal_profile"]) == [
        "unsafe_normal",
        "disposable_above_cap",
    ]
    assert pd.Series(syn_rec["in_default_subtraction_background"]).all()

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Tests for quantitative unmeasured home-tissue normal panels and adjudications (#735)."""

from __future__ import annotations

import pytest

import oncoref
from oncoref import data_manifest, reference_data, target_selection, unmeasured_home_panels


def test_unmeasured_home_tissue_panels_covers_all_tissues_and_unit_contracts():
    panels = oncoref.unmeasured_home_tissue_panels()
    assert len(panels) == 11
    assert set(panels["home_tissue"]) == set(unmeasured_home_panels.HOME_TISSUES_735)
    assert set(panels["modality"]) == set(unmeasured_home_panels.PANEL_MODALITIES_735)

    # Every SHA-256 is a lowercase 64-hex digest
    assert panels["sha256"].str.fullmatch(r"[0-9a-f]{64}").all()

    # Unit separation: bulk panels are strictly TPM, scRNA-seq strictly UMI, qPCR strictly TaqMan_copies
    assert set(panels.loc[panels["modality"].eq("bulk_rnaseq_tpm"), "native_unit"]) == {"TPM"}
    assert set(panels.loc[panels["modality"].eq("scrnaseq_umi"), "native_unit"]) == {"UMI"}
    assert set(panels.loc[panels["modality"].eq("targeted_qpcr_copies"), "native_unit"]) == {
        "TaqMan_copies"
    }

    hf = oncoref.unmeasured_home_tissue_panels(home_tissue="hair follicle", modality="scrnaseq_umi")
    assert set(hf["accession"]) == {"GSE129611", "GSE193269"}
    assert int(hf["n_donors"].sum()) == 8

    with pytest.raises(ValueError, match="Unknown home_tissue"):
        oncoref.unmeasured_home_tissue_panels(home_tissue="spleen")
    with pytest.raises(ValueError, match="Unknown modality"):
        oncoref.unmeasured_home_tissue_panels(modality="microarray")


def test_unmeasured_home_protein_adjudications_matches_catalog_and_key_outcomes():
    ts = target_selection.results()
    expected_home = ts[ts["home_tissue"].fillna("").ne("")]
    adj = oncoref.unmeasured_home_protein_adjudications()
    assert len(adj) == len(expected_home) == 571
    assert set(adj["entity_id"]) == set(expected_home["entity_id"])
    assert set(adj["home_panel_disposition"]) <= set(
        unmeasured_home_panels.HOME_PANEL_DISPOSITIONS_735
    )

    gated = oncoref.unmeasured_home_protein_adjudications(only_ca_cta_gated=True)
    assert len(gated) == 95
    assert not gated["home_panel_disposition"].eq("already_excluded_upstream").any()

    by_label = gated.set_index("protein_group_label")

    # Vital/Manageable ceiling exclusions, including PANX3 cross-panel Vital bone detection
    above = oncoref.unmeasured_home_protein_adjudications(
        disposition="home_expressed_above_tier_ceiling",
        only_ca_cta_gated=True,
    )
    assert set(above["protein_group_label"]) == {"SP7", "PANX3", "BPIFB3", "GFY", "GLYATL3"}
    assert float(by_label.loc["SP7", "home_bulk_p95_tpm"]) > 150.0
    assert "bone:" in str(by_label.loc["PANX3", "cross_panel_ceiling_violation"])
    assert float(by_label.loc["GFY", "home_bulk_p95_tpm"]) > 300.0

    # Insufficient 10x tumor-over-home fold
    insuff = oncoref.unmeasured_home_protein_adjudications(
        disposition="home_expressed_insufficient_tumor_fold",
        only_ca_cta_gated=True,
    )
    assert "PABPC3" in set(insuff["protein_group_label"])
    assert float(by_label.loc["PABPC3", "tumor_over_home_bulk_fold"]) < 10.0

    # Disposable hair follicle single-cell UMI confirmation
    disp_hf = oncoref.unmeasured_home_protein_adjudications(
        disposition="home_disposable_follicle_expressed_tenfold",
        only_ca_cta_gated=True,
    )
    assert set(disp_hf["protein_group_label"]) == {"KRTAP21-2", "KRTAP9-7"}

    # Cleared below ceiling and >=10x tumor fold (including hair-bearing skin P99<=0.50 & max<4.0 KRTAP20-4)
    assert by_label.loc["TAS1R2", "home_panel_disposition"] == "home_cleared_panel_passes_tenfold"
    assert by_label.loc["OR10J1", "home_panel_disposition"] == "home_cleared_panel_passes_tenfold"
    assert (
        by_label.loc["KRTAP20-4", "home_panel_disposition"] == "home_cleared_panel_passes_tenfold"
    )

    # Mineralized odontoblast-front DSPP remains explicitly undetected in soft pulp biopsies
    assert by_label.loc["DSPP", "home_panel_disposition"] == "home_undetected_in_sampled_panel"

    with pytest.raises(ValueError, match="Unknown disposition"):
        oncoref.unmeasured_home_protein_adjudications(disposition="unknown")


def test_unmeasured_home_gene_expression_markers_and_externalized_registration():
    assert "unmeasured-home-gene-expression" in data_manifest.EXTERNALIZED_ORIGINATED
    prov = reference_data.provenance("unmeasured-home-gene-expression")
    assert prov["sha256"] == "28c1a03739035ec7185847be276fc2c788c6ba8b5329476a761b43717c619fbe"
    assert "reference-unmeasured-homes-20261122" in prov["url"]

    krt85 = oncoref.unmeasured_home_gene_expression(gene="KRT85")
    assert len(krt85) == 1
    assert int(krt85["hair_sc_studies_ge2"].iloc[0]) == 2
    assert int(krt85["hair_sc_donors_ge2"].iloc[0]) >= 6
    assert float(krt85["hair_sc_max_pct_ge2"].iloc[0]) > 15.0

    alpl = oncoref.unmeasured_home_gene_expression(gene="ALPL")
    assert float(alpl["bone_bulk_p95_tpm"].iloc[0]) > 400.0
    assert float(alpl["bone_sc_max_pct_ge2"].iloc[0]) > 25.0

    summary = oncoref.unmeasured_home_tissue_summary()
    assert len(summary) == 6
    assert int(summary["total_home_proteins"].sum()) == 571
    assert int(summary["ca_cta_gated_proteins"].sum()) == 95

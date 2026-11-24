# Licensed under the Apache License, Version 2.0.
"""Regression tests for #765 nuclear snRNA-seq compartment & transcriptional spillover QC."""

from __future__ import annotations

import pytest

from oncoref import human_screening_inputs, rare_cell_compartment_qc


def test_spillover_geometry_and_intronless_cds_contract():
    # PLAC1: 197.7 kb non-coding 5'-UTR span (0.32% CDS) + 5.2 kb downstream of PABIR2 (54 TPM) + 1 CDS exon
    plac1 = rare_cell_compartment_qc.assess_transcriptional_spillover_geometry(
        same_strand_upstream_bp=5211,
        same_strand_partner_tpm=54.05,
        noncoding_5prime_span_bp=197673,
        cds_span_ratio=0.003219,
        n_cds_exons=1,
    )
    assert plac1["transcriptional_spillover_risk"] is True
    assert "same_strand_upstream_readthrough" in plac1["spillover_mechanisms"]
    assert "long_noncoding_5prime_intron_span" in plac1["spillover_mechanisms"]
    assert plac1["is_intronless_cds"] is True
    # Guardrail 2: single-CDS-exon ORFs accept CDS-exonic reads rather than demanding splice junctions
    assert plac1["accepted_cds_confirmation_assays"] == ("cds_exonic_reads",)

    # NUTM1: 135 bp head-to-head bidirectional CpG promoter with NOP10 (4,501 TPM in neutrophils)
    nutm1 = rare_cell_compartment_qc.assess_transcriptional_spillover_geometry(
        bidirectional_promoter_bp=135,
        bidirectional_partner_tpm=4501.1,
        n_cds_exons=8,
    )
    assert nutm1["transcriptional_spillover_risk"] is True
    assert nutm1["spillover_mechanisms"] == ("bidirectional_promoter_prompt",)
    assert nutm1["is_intronless_cds"] is False
    assert "spliced_exon_junction" in nutm1["accepted_cds_confirmation_assays"]
    assert "cds_exonic_reads" in nutm1["accepted_cds_confirmation_assays"]


def test_magea12_bulk_dilution_and_literature_guardrails_stay_unsafe():
    # Guardrail 1: MAGEA12 (0.30 TPM brain bulk, literature_effect='confirm_detection') must NEVER clear,
    # even if a caller passes hypothetical spillover geometry.
    magea12_lit = rare_cell_compartment_qc.evaluate_nuclear_rare_cell_disposition(
        baseline_rare_cell_status="established",
        nucleus_only_signal=True,
        finding_bulk_tpm=0.01,
        somatic_ihc_positive=False,
        literature_effect="confirm_detection",
        other_unsafe_gates=False,
        transcriptional_spillover_risk=True,
        paralog_ambiguous=False,
        host_overlapped=False,
    )
    assert magea12_lit["cleared"] is False
    assert magea12_lit["effective_rare_cell_action"] == "unsafe"

    # Even without literature_effect, 0.30 TPM bulk (>= 0.05 TPM ceiling) blocks clearance
    diluted_rare = rare_cell_compartment_qc.evaluate_nuclear_rare_cell_disposition(
        baseline_rare_cell_status="established",
        nucleus_only_signal=True,
        finding_bulk_tpm=0.30,
        somatic_ihc_positive=False,
        literature_effect="not_reviewed",
        other_unsafe_gates=False,
        transcriptional_spillover_risk=True,
        paralog_ambiguous=False,
        host_overlapped=False,
    )
    assert diluted_rare["cleared"] is False
    assert diluted_rare["disposition"] == "retained_established_cytoplasmic_or_bulk"

    # NaN / None inputs must fail closed rather than silently clearing via IEEE-754 NaN >= 0.05 == False
    geom_nan = rare_cell_compartment_qc.assess_transcriptional_spillover_geometry(
        bidirectional_promoter_bp=float("nan"),
        same_strand_upstream_bp=float("nan"),
    )
    assert geom_nan["transcriptional_spillover_risk"] is False

    for bad_bulk in (float("nan"), None, -0.01):
        nan_disp = rare_cell_compartment_qc.evaluate_nuclear_rare_cell_disposition(
            baseline_rare_cell_status="established",
            nucleus_only_signal=True,
            finding_bulk_tpm=bad_bulk,
            somatic_ihc_positive=False,
            literature_effect="not_reviewed",
            other_unsafe_gates=False,
            transcriptional_spillover_risk=True,
            paralog_ambiguous=False,
            host_overlapped=False,
        )
        assert nan_disp["cleared"] is False
        assert nan_disp["effective_rare_cell_action"] == "unsafe"


def test_packaged_nuclear_spillover_tables_and_catalog_transitions():
    comp = rare_cell_compartment_qc.census_collection_compartments()
    assert len(comp) == 93
    assert set(comp["assay_compartment"]) <= set(rare_cell_compartment_qc.ASSAY_COMPARTMENTS)

    adj = rare_cell_compartment_qc.nuclear_spillover_adjudications()
    by_label = adj.set_index("protein_group_label")

    # MAGEA12 negative-control anchor remains strictly retained/Unsafe
    assert bool(by_label.loc["MAGEA12", "nuclear_spillover_cleared"]) is False
    assert by_label.loc["MAGEA12", "disposition"] == "retained_literature_or_ihc_supported"

    # NUTM1 135 bp NOP10 bidirectional promoter spillover recorded
    assert bool(by_label.loc["NUTM1", "nuclear_spillover_cleared"]) is True
    assert by_label.loc["NUTM1", "bidirectional_promoter_bp"] == pytest.approx(135.0)
    assert "NOP10:" in str(by_label.loc["NUTM1", "bidirectional_partner"])

    # PLAC1 clears nuclear snRNA-seq spillover and advances Unsafe -> CTA with nuclear_only caution
    assert bool(by_label.loc["PLAC1", "nuclear_spillover_cleared"]) is True
    assert float(by_label.loc["PLAC1", "finding_bulk_tpm"]) < 0.01

    res = human_screening_inputs.evaluate_human_catalog(include_nuclear_spillover_qc=True)
    cat = res.catalog.set_index("protein_group_label")
    assert cat.loc["MAGEA12", "classification"] == "Unsafe"
    assert cat.loc["PLAC1", "classification"] == "CTA"
    assert cat.loc["PLAC1", "rare_cell_status"] == "nuclear_only"
    assert cat.loc["PLAC1", "nuclear_rna_caution"] == "snrnaseq_intronic_spillover_cleared"
    assert float(cat.loc["PLAC1", "tumor_max_tpm"]) == pytest.approx(152.613, rel=1e-4)
    assert float(cat.loc["PLAC1", "reproductive_supported_tpm"]) == pytest.approx(
        187.021784, rel=1e-4
    )
    assert cat.loc["TENT5D", "classification"] == "Reproductive Antigen only"
    # DUX4 remains Other:Review until #767 clears its closed D4Z4 paralog ambiguity
    assert cat.loc["DUX4", "classification"] == "Other:Review"

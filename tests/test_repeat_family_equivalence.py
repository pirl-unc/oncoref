"""Tests for closed homologous repeat-family sequence equivalence & symmetric quantification (#767)."""

from __future__ import annotations

import json

from oncoref import (
    antigen_sets,
    human,
    repeat_family_equivalence,
    sequence_bundle,
    target_selection,
)


def test_partition_sequence_partners_separates_closed_arrays_from_unsafe_paralogs():
    # Closed D4Z4 repeat array: all partners >=85% CDS 31-mers, unprocessed_pseudogene, zero Unsafe
    dux_res = repeat_family_equivalence.partition_sequence_partners(
        partners=[
            {
                "gene_id": "ENSG00000260832",
                "symbol": "DUX4L2",
                "biotype": "unprocessed_pseudogene",
                "cds_31mer_recall": 0.96,
                "is_unsafe": False,
            },
            {
                "gene_id": "ENSG00000259958",
                "symbol": "DUX4L4",
                "biotype": "unprocessed_pseudogene",
                "cds_31mer_recall": 0.95,
                "is_unsafe": False,
            },
        ],
        unsafe_sequence_overlap_pct=0.72,
        family_somatic_max_tpm=1.72,
    )
    assert dux_res["closed_repeat_family_cleared"] is True
    assert dux_res["disposition"] == "closed_repeat_family_cleared"
    assert dux_res["p_ext"] == ()
    assert len(dux_res["p_self"]) == 2

    # Negative control 1: PSG8 shares sequence with Unsafe PSG1 (unsafe_sequence_overlap_pct = 74.07%)
    psg8_res = repeat_family_equivalence.partition_sequence_partners(
        partners=[
            {
                "gene_id": "ENSG00000231924",
                "symbol": "PSG1",
                "biotype": "protein_coding",
                "cds_31mer_recall": 0.88,
                "is_unsafe": True,
            }
        ],
        unsafe_sequence_overlap_pct=74.07,
        family_somatic_max_tpm=0.0,
    )
    assert psg8_res["closed_repeat_family_cleared"] is False
    assert psg8_res["disposition"] == "retained_external_unsafe_or_somatic_paralog"

    # Negative control 2: processed_pseudogene or lncRNA host partner blocks closed-family clearance
    proc_res = repeat_family_equivalence.partition_sequence_partners(
        partners=[
            {
                "gene_id": "ENSG00000248542",
                "symbol": "RPL10P9",
                "biotype": "processed_pseudogene",
                "cds_31mer_recall": 0.92,
                "is_unsafe": False,
            }
        ],
        unsafe_sequence_overlap_pct=0.0,
        family_somatic_max_tpm=0.0,
    )
    assert proc_res["closed_repeat_family_cleared"] is False
    assert proc_res["disposition"] == "retained_processed_pseudogene_or_host_partner"

    # Negative control 3: NaN / None ceilings must fail closed (never clear via NaN >= 5.0 == False)
    nan_res = repeat_family_equivalence.partition_sequence_partners(
        partners=[
            {
                "gene_id": "ENSG00000260832",
                "symbol": "DUX4L2",
                "biotype": "unprocessed_pseudogene",
                "cds_31mer_recall": 0.96,
                "is_unsafe": False,
            }
        ],
        unsafe_sequence_overlap_pct=float("nan"),
        family_somatic_max_tpm=0.0,
    )
    assert nan_res["closed_repeat_family_cleared"] is False
    nan_som = repeat_family_equivalence.partition_sequence_partners(
        partners=[
            {
                "gene_id": "ENSG00000260832",
                "symbol": "DUX4L2",
                "biotype": "unprocessed_pseudogene",
                "cds_31mer_recall": 0.96,
                "is_unsafe": False,
            }
        ],
        unsafe_sequence_overlap_pct=0.5,
        family_somatic_max_tpm=float("nan"),
    )
    assert nan_som["closed_repeat_family_cleared"] is False


def test_dux4_joint_765_and_767_dependency_and_symmetric_fshd_muscle_silence():
    adj = repeat_family_equivalence.repeat_family_adjudications()
    by_lbl = adj.set_index("protein_group_label")

    assert bool(by_lbl.loc["DUX4", "closed_repeat_family_cleared"]) is True
    assert by_lbl.loc["DUX4", "disposition"] == "closed_repeat_family_cleared"
    # FSHD skeletal muscle locus is verified silent (<0.1 TPM) across all DUX4 + DUX4L* copies
    assert float(by_lbl.loc["DUX4", "gtex_skeletal_muscle_p95_tpm"]) < 0.1
    assert float(by_lbl.loc["DUX4", "family_somatic_max_tpm"]) < 5.0

    # Negative controls with Unsafe paralogs remain blocked
    for neg_lbl in ("PSG8", "POTEC", "TUBB8B", "GAGE2E", "GAGE10"):
        assert bool(by_lbl.loc[neg_lbl, "closed_repeat_family_cleared"]) is False

    # 1. #765 alone: DUX4 remains Other:Review due to paralog_ambiguous = True
    res_765_only = human.evaluate_integrated_catalog(
        include_nuclear_spillover_qc=True,
        include_repeat_family_equivalence=False,
    ).catalog.set_index("protein_group_label")
    assert res_765_only.loc["DUX4", "classification"] == "Other:Review"

    # 2. #767 alone: DUX4 remains Other:Review due to unresolved brain snRNA-seq flag
    res_767_only = human.evaluate_integrated_catalog(
        include_nuclear_spillover_qc=False,
        include_repeat_family_equivalence=True,
    ).catalog.set_index("protein_group_label")
    assert res_767_only.loc["DUX4", "classification"] == "Other:Review"

    # 3. Joint #765 + #767: both holds clear and DUX4 advances to CTA
    res_joint = human.evaluate_integrated_catalog(
        include_nuclear_spillover_qc=True,
        include_repeat_family_equivalence=True,
    ).catalog.set_index("protein_group_label")
    assert res_joint.loc["DUX4", "classification"] == "CTA"
    assert res_joint.loc["DUX4", "category"] == "CTA"
    assert res_joint.loc["DUX4", "sex_restriction_code"] == "CTA-P"
    assert res_joint.loc["DUX4", "rare_cell_status"] == "nuclear_only"
    assert res_joint.loc["DUX4", "repeat_family_767_disposition"] == "closed_repeat_family_cleared"

    # VCX3B has whole-cell scRNA-seq signal in eye (not nucleus-only), so even with #765 + #767 it stays Other:Review
    assert res_joint.loc["VCX3B", "classification"] == "Other:Review"


def test_tsarina_closed_coding_family_sibling_contract_on_catalog_and_bundle():
    ts = target_selection.results().set_index("protein_group_label")
    disp = sequence_bundle.target_selection_dispositions().set_index("protein_group_label")
    cta_proteins = antigen_sets.proteins("CTA").set_index("protein_group_label")

    for col in repeat_family_equivalence.REPEAT_FAMILY_COLUMNS:
        assert col in ts.columns
        assert col in disp.columns
        assert col in cta_proteins.columns

    # CT47A1..12 and CT47B1 list somatic-safe RA-only sibling CT47C1 (ENSG00000277535)
    # and never include untranslated/Other:Review ENSG00000278646
    ct47_key = "CT47A1/10/11/12/2/3/4/5/6/7/8/9"
    assert ts.loc[ct47_key, "closed_homolog_family_id"] == "CT47"
    ct47_sibs = set(json.loads(ts.loc[ct47_key, "co_silenced_coding_sibling_ids"]))
    assert "ENSG00000277535" in ct47_sibs  # CT47C1 (RA only)
    assert "ENSG00000236446" in ct47_sibs  # CT47B1 (CTA)
    assert "ENSG00000278646" not in ct47_sibs

    # NY-ESO-1 and CTAG2 share closed_homolog_family_id == "CTAG", and CTAG2 lists NY-ESO-1 genes
    assert ts.loc["NY-ESO-1", "closed_homolog_family_id"] == "CTAG"
    assert ts.loc["CTAG2", "closed_homolog_family_id"] == "CTAG"
    ctag2_sibs = set(json.loads(ts.loc["CTAG2", "co_silenced_coding_sibling_ids"]))
    assert {"ENSG00000184033", "ENSG00000268651"} <= ctag2_sibs

    # TSPY1 lists somatic-safe RA-only TSPY9 (ENSG00000238074) and NEVER Unsafe TSPYL1..6
    assert ts.loc["TSPY1", "closed_homolog_family_id"] == "TSPY"
    tspy_sibs = set(json.loads(ts.loc["TSPY1", "co_silenced_coding_sibling_ids"]))
    assert "ENSG00000238074" in tspy_sibs  # TSPY9 (RA only)
    unsafe_tspyl_ids = {
        "ENSG00000178021",
        "ENSG00000180543",
        "ENSG00000184205",
        "ENSG00000187189",
        "ENSG00000189241",
    }
    assert tspy_sibs.isdisjoint(unsafe_tspyl_ids)

    # SPANXN1 lists RA-only SPANXN2 (ENSG00000268988) and SPANXN5 (ENSG00000204363)
    spanxn_sibs = set(json.loads(ts.loc["SPANXN1", "co_silenced_coding_sibling_ids"]))
    assert {"ENSG00000268988", "ENSG00000204363"} <= spanxn_sibs

    # CSAG1 lists CA-only CSAG2 (ENSG00000268902) and CSAG3 (ENSG00000268916)
    csag_sibs = set(json.loads(ts.loc["CSAG1", "co_silenced_coding_sibling_ids"]))
    assert {"ENSG00000268902", "ENSG00000268916"} <= csag_sibs

    # DUX4 lists all 20 primary-assembly DUX4L* D4Z4 repeat paralogs in closed_self_paralog_ids
    dux_pars = json.loads(ts.loc["DUX4", "closed_self_paralog_ids"])
    assert len(dux_pars) >= 20

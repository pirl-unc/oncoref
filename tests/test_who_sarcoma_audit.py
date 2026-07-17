# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""WHO soft-tissue/bone parity and reference-backing contracts (#378)."""

from oncoref import (
    WHO_AUDIT_STATUS_VALUES,
    cancer_type_descendants,
    cancer_type_records,
    cancer_type_registry,
    cancer_who_soft_tissue_bone_audit,
    cohort_aggregate_members,
    resolve_cancer_type,
)
from oncoref.expression_registry import expression_source_candidates
from oncoref.fusions import cancer_fusions
from oncoref.load_dataset import get_data

NEW_ENTITY_CODES = {
    "SARC_ADAMANTINOMA",
    "SARC_AFS",
    "SARC_CHON_CLEAR_CELL",
    "SARC_CHON_DEDIFF",
    "SARC_CHON_MESENCHYMAL",
    "SARC_CHON_PERIOSTEAL",
    "SARC_CHOR_CONVENTIONAL",
    "SARC_CHOR_DEDIFF",
    "SARC_CHOR_POORLY_DIFF",
    "SARC_EBV_SMT",
    "SARC_ECTOMES",
    "SARC_EIMS",
    "SARC_EWSR1_NONETS",
    "SARC_FIBROSARCOMA_BONE",
    "SARC_GCT_MALIGNANT",
    "SARC_GLOMUS_MALIGNANT",
    "SARC_ICMT",
    "SARC_ILMS",
    "SARC_INTIMAL",
    "SARC_MIFS",
    "SARC_MYOFIBROBLASTIC",
    "SARC_NTRK_SPINDLE",
    "SARC_OFMT_MALIGNANT",
    "SARC_OS_CONVENTIONAL",
    "SARC_OS_EXTRASKELETAL",
    "SARC_OS_HIGH_GRADE_SURFACE",
    "SARC_OS_LOW_GRADE_CENTRAL",
    "SARC_OS_PAROSTEAL",
    "SARC_OS_PERIOSTEAL",
    "SARC_OS_SECONDARY",
    "SARC_PERINEURIOMA_MALIGNANT",
    "SARC_PMT_MALIGNANT",
    "SARC_SCD34FT",
    "SARC_SEF",
    "SARC_TGCT_MALIGNANT",
}


def test_who_audit_is_checked_and_machine_readable():
    audit = cancer_who_soft_tissue_bone_audit()
    registry_codes = set(cancer_type_registry()["code"])

    assert list(audit.columns) == [
        "who_entity",
        "who_category",
        "who_behavior",
        "registry_status",
        "registry_code",
        "source_url",
        "notes",
    ]
    assert set(audit["registry_status"]) == set(WHO_AUDIT_STATUS_VALUES)
    assert audit["source_url"].str.startswith("https://").all()
    assert audit["notes"].str.strip().ne("").all()

    linked = audit[audit["registry_status"].isin({"represented", "alias", "axis"})]
    unlinked = audit[audit["registry_status"].isin({"missing", "out_of_scope"})]
    assert linked["registry_code"].notna().all()
    assert set(linked["registry_code"]) <= registry_codes
    assert unlinked["registry_code"].isna().all()

    represented = audit[audit["registry_status"] == "represented"]
    assert not represented["registry_code"].duplicated().any()
    assert audit[
        (audit["who_behavior"] == "malignant") & (audit["registry_status"] == "missing")
    ].empty


def test_who_audit_tracks_known_intermediate_gaps():
    missing = set(cancer_who_soft_tissue_bone_audit(registry_status="missing")["who_entity"])
    assert {
        "Desmoid fibromatosis",
        "Plexiform fibrohistiocytic tumour",
        "Kaposiform haemangioendothelioma",
        "Atypical fibroxanthoma",
        "Desmoplastic fibroma of bone",
    } <= missing


def test_who_metadata_is_public_and_queryable_without_changing_the_tree():
    registry = cancer_type_registry().set_index("code")
    assert registry.loc["SARC_EWS", "who_category"] == (
        "undifferentiated_small_round_cell_sarcomas"
    )
    assert registry.loc["SARC_EWS", "who_behavior"] == "malignant"
    assert registry.loc["SARC_KS", "who_behavior"] == "intermediate_rarely_metastasizing"

    records = cancer_type_records(who_category="vascular_tumours")
    assert {"SARC_KS", "SARC_EHE", "SARC_ANGIO"} <= set(records["code"])
    malignant = cancer_type_records(under="SARC_CHON", who_behavior="malignant", include_self=True)
    assert set(malignant["code"]) == {
        "SARC_CHON_PERIOSTEAL",
        "SARC_CHON_CLEAR_CELL",
        "SARC_CHON_MESENCHYMAL",
        "SARC_CHON_DEDIFF",
    }


def test_new_entities_exist_but_are_not_claimed_as_exact_reference_targets():
    records = cancer_type_records(sorted(NEW_ENTITY_CODES)).set_index("code")
    assert set(records.index) == NEW_ENTITY_CODES
    assert set(records["reference_source"]) == {"parent"}
    assert not records["is_classification_target"].any()
    assert not records["has_expression_matrix"].any()
    assert records["classification_reference_code"].notna().all()

    for code in NEW_ENTITY_CODES:
        assert resolve_cancer_type(code) == code


def test_round_cell_tree_is_molecular_and_not_a_sample_aggregate():
    assert set(cancer_type_descendants("SARC_ROUND_CELL")) == {
        "SARC_EWS",
        "SARC_CIC",
        "SARC_BCOR",
        "SARC_EWSR1_NONETS",
    }
    assert "SARC_DSRCT" not in cancer_type_descendants("SARC_ROUND_CELL")

    grouping = cancer_type_records("SARC_ROUND_CELL").iloc[0]
    assert grouping["ontology_level"] == "grouping"
    assert grouping["reference_source"] == "parent"
    assert grouping["classification_reference_code"] == "SARC"
    assert not bool(grouping["is_classification_target"])
    assert "SARC_ROUND_CELL" not in cohort_aggregate_members("SARC")


def test_chondrosarcoma_tree_and_fusion_ownership_are_entity_specific():
    assert set(cancer_type_descendants("SARC_CHON")) == {
        "SARC_CHON_PERIOSTEAL",
        "SARC_CHON_CLEAR_CELL",
        "SARC_CHON_MESENCHYMAL",
        "SARC_CHON_DEDIFF",
    }

    broad_fusions = set(cancer_fusions("SARC_CHON")["fusion_family"])
    mesenchymal_fusions = set(cancer_fusions("SARC_CHON_MESENCHYMAL")["fusion_family"])
    assert "HEY1-NCOA2" not in broad_fusions
    assert {"HEY1-NCOA2", "IRF2BP2-CDX1"} <= mesenchymal_fusions

    key_genes = get_data("cancer-key-genes")
    hey1 = key_genes[
        (key_genes["symbol"] == "HEY1") & key_genes["cancer_code"].astype(str).str.contains("CHON")
    ]
    assert hey1["cancer_code"].tolist() == ["SARC_CHON_MESENCHYMAL"]


def test_osteosarcoma_and_chordoma_trees_are_explicit():
    assert set(cancer_type_descendants("SARC_OS")) == {
        "SARC_OS_LOW_GRADE_CENTRAL",
        "SARC_OS_CONVENTIONAL",
        "SARC_OS_PAROSTEAL",
        "SARC_OS_PERIOSTEAL",
        "SARC_OS_HIGH_GRADE_SURFACE",
        "SARC_OS_SECONDARY",
        "SARC_OS_EXTRASKELETAL",
    }
    assert set(cancer_type_descendants("SARC_CHOR")) == {
        "SARC_CHOR_CONVENTIONAL",
        "SARC_CHOR_DEDIFF",
        "SARC_CHOR_POORLY_DIFF",
    }


def test_new_entities_have_one_explicit_expression_source_follow_up():
    candidates = expression_source_candidates()
    rows = candidates[candidates["cancer_code"].isin(NEW_ENTITY_CODES)]
    assert set(rows["cancer_code"]) == NEW_ENTITY_CODES
    assert not rows["cancer_code"].duplicated().any()
    assert set(rows["source_status"]) == {"source_needed"}
    assert rows["reference_code"].notna().all()
    assert rows["notes"].str.contains("not exact entity evidence").all()


def test_common_who_name_variants_resolve():
    assert resolve_cancer_type("mesenchymal chondrosarcoma") == "SARC_CHON_MESENCHYMAL"
    assert resolve_cancer_type("extrarenal rhabdoid tumour") == "RT"
    assert resolve_cancer_type("intracranial mesenchymal tumor") == "SARC_ICMT"
    assert resolve_cancer_type("EBV-associated smooth muscle tumor") == "SARC_EBV_SMT"


def test_registry_still_has_no_dangling_parents():
    registry = cancer_type_registry()
    codes = set(registry["code"])
    parents = set(registry["parent_code"].dropna().astype(str)) - {""}
    assert parents <= codes

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Per-type gene biology + fusion-rule tables (#35, R-onto)."""

import pandas as pd

from oncoref import cancer_genes as cg
from oncoref import fusions as f


def _split_refs(value):
    return {x.strip() for x in str(value).split(";") if x.strip()}


def test_viral_antigens():
    assert cg.cancer_viral_antigens("HPV") == ["E6", "E7"]
    assert cg.cancer_viral_antigens("not_a_virus") == []
    all_map = cg.cancer_viral_antigens()
    assert "HPV" in all_map
    # reverse lookup: a cervical/HPV cohort should surface HPV antigens
    pairs = cg.viral_antigens_for_cancer("CESC")
    assert any(v == "HPV" for v, _ in pairs)


def test_key_genes_roles():
    df = cg.cancer_key_genes_df()
    assert {"cancer_code", "symbol", "role"} <= set(df.columns)
    # PRAD biomarkers should include AR (androgen receptor)
    biomarkers = cg.cancer_biomarker_genes("PRAD")
    assert "AR" in biomarkers
    targets = cg.cancer_therapy_targets("PRAD")
    assert (targets["role"].astype(str) == "target").all()


def test_key_gene_citation_audit_covers_every_citation():
    key_columns = ["cancer_code", "subtype", "symbol", "role", "agent", "pmid"]
    key_genes = cg.cancer_key_genes_df().fillna("")
    cited = key_genes.assign(pmid=key_genes["source"].str.split("; ")).explode("pmid")
    cited = cited[key_columns].sort_values(key_columns).reset_index(drop=True)

    audit = cg.cancer_key_gene_citation_audit().fillna("")
    audited = audit[key_columns].sort_values(key_columns).reset_index(drop=True)
    pd.testing.assert_frame_equal(audited, cited)

    assert audit["pubmed_title"].str.strip().ne("").all()
    assert audit["publication_year"].between(1900, 2100).all()
    assert audit["supports_gene_or_agent"].all()
    assert audit["supports_disease_context"].all()
    assert audit["supports_role_and_phase"].all()
    reviewed_on = pd.to_datetime(audit["reviewed_on"], format="%Y-%m-%d", errors="raise")
    assert reviewed_on.notna().all()


def test_key_gene_citation_audit_returns_a_defensive_copy():
    first = cg.cancer_key_gene_citation_audit()
    first.loc[0, "pubmed_title"] = "changed"
    assert cg.cancer_key_gene_citation_audit().loc[0, "pubmed_title"] != "changed"


def test_key_gene_known_wrong_pmids_are_replaced():
    key = cg.cancer_key_genes_df()

    expected_sources = {
        ("PRAD", "AR", "biomarker", ""): "PMID:22894553",
        ("PRAD", "KLK3", "biomarker", ""): "PMID:25153393",
        ("PRAD", "KLK2", "biomarker", ""): "PMID:25153393",
        ("PRAD", "NKX3-1", "biomarker", ""): "PMID:35265947",
        ("PRAD", "DLL3", "target", "tarlatamab"): "PMID:40689871",
        ("PRAD", "KLK2", "target", "pasritamig"): "PMID:40450573",
        ("PRAD", "PSCA", "biomarker", ""): "PMID:15342669",
        ("BRCA", "ERBB2", "biomarker", ""): "PMID:11248153",
        ("LUAD", "RET", "biomarker", ""): "PMID:32846060",
        ("LUAD", "MET", "biomarker", ""): "PMID:32877583",
        ("LUAD", "EGFR", "target", "osimertinib"): "PMID:29151359",
        ("COAD", "BRAF", "biomarker", ""): "PMID:31566309",
        ("COAD", "CEACAM5", "biomarker", ""): "PMID:17060676",
        ("DLBC", "MS4A1", "biomarker", ""): "PMID:39017945",
        ("DLBC", "MS4A1", "target", "rituximab"): "PMID:16702182",
        ("DLBC", "CD19", "target", "axicabtagene ciloleucel"): "PMID:29226797",
        ("LAML", "CD33", "biomarker", ""): "PMID:22482940",
        ("LAML", "CD33", "target", "gemtuzumab ozogamicin"): "PMID:30076173",
        ("LAML", "FLT3", "target", "gilteritinib"): "PMID:31665578",
        ("THCA", "RET", "biomarker", ""): "PMID:32846061",
    }
    for (code, symbol, role, agent), expected in expected_sources.items():
        mask = (
            (key["cancer_code"] == code)
            & (key["symbol"] == symbol)
            & (key["role"] == role)
            & (key["agent"].fillna("") == agent)
        )
        rows = key[mask]
        assert len(rows) == 1, (code, symbol, role, agent)
        assert rows.iloc[0]["source"] == expected

    bad_pmids = {
        "PMID:23332746",
        "PMID:17538631",
        "PMID:36417474",
        "PMID:38110199",
        "PMID:9653118",
        "PMID:17544441",
        "PMID:32273263",
        "PMID:31950082",
        "PMID:31733398",
        "PMID:24637364",
        "PMID:15767599",
        "PMID:9562152",
        "PMID:29091450",
        "PMID:21233305",
        "PMID:31634902",
        "PMID:12075054",
        "PMID:15178638",
        "PMID:20181787",
        "PMID:32427717",
        "PMID:33667670",
        "PMID:36535566",
    }
    all_refs = set().union(*(_split_refs(v) for v in key["source"]))
    assert bad_pmids.isdisjoint(all_refs)


def test_key_gene_nonexistent_pmids_are_replaced():
    key = cg.cancer_key_genes_df()

    bad_pmids = {"PMID:34428009", "PMID:35379757", "PMID:37173835"}
    all_refs = set().union(*(_split_refs(v) for v in key["source"]))
    assert bad_pmids.isdisjoint(all_refs)

    expected_refs = {
        ("SKCM", "", "PRAME", "biomarker", ""): {"PMID:38338862"},
        ("SKCM", "", "PRAME", "target", "IMA203"): {"PMID:40205198"},
        ("SARC", "myxoid_liposarcoma", "PRAME", "biomarker", ""): {"PMID:27499900"},
        ("SARC", "synovial_sarcoma", "PRAME", "biomarker", ""): {"PMID:30524904"},
        ("SARC", "ewing_sarcoma", "PRAME", "biomarker", ""): {"PMID:24973179"},
        ("SARC", "MPNST", "NF1", "biomarker", ""): {"PMID:36598417"},
        ("SARC", "MPNST", "CDKN2A", "biomarker", ""): {"PMID:36598417"},
        ("SARC", "MPNST", "CDKN2B", "biomarker", ""): {"PMID:14519636"},
        ("SARC", "MPNST", "TP53", "biomarker", ""): {"PMID:36598417"},
        ("SARC", "MPNST", "SOX10", "biomarker", ""): {"PMID:28551330"},
        ("SARC", "MPNST", "S100B", "biomarker", ""): {"PMID:28551330"},
        ("SARC", "MPNST", "MAP2K1", "target", "trametinib"): {"PMID:32975370"},
        ("SARC", "MPNST", "PTPN11", "target", "SHP2 inhibitors"): {"PMID:33032988"},
    }
    for (code, subtype, symbol, role, agent), expected in expected_refs.items():
        mask = (
            (key["cancer_code"] == code)
            & (key["subtype"].fillna("") == subtype)
            & (key["symbol"] == symbol)
            & (key["role"] == role)
            & (key["agent"].fillna("") == agent)
        )
        rows = key[mask]
        assert len(rows) == 1, (code, subtype, symbol, role, agent)
        assert _split_refs(rows.iloc[0]["source"]) == expected



def test_key_gene_claims_are_structurally_unambiguous():
    key = cg.cancer_key_genes_df().fillna("")
    claim_key = ["cancer_code", "subtype", "symbol", "role", "agent"]

    assert not key.duplicated(claim_key).any()
    assert key["symbol"].str.strip().ne("").all()
    assert key["source"].str.fullmatch(r"PMID:\d+(; PMID:\d+)*").all()
    assert key.loc[key["role"] == "target", "agent"].str.strip().ne("").all()

    noncanonical_symbols = {
        "CD20",
        "CD123",
        "HIF2A",
        "MAGE-A4",
        "MEK1",
        "T",
        "TROP2",
        "VEGFR2",
    }
    assert noncanonical_symbols.isdisjoint(set(key["symbol"]))


def test_issue_161_status_and_context_corrections():
    key = cg.cancer_key_genes_df().fillna("")

    corrected_phases = {
        ("BLCA", "sacituzumab govitecan + pembrolizumab"): "phase_2",
        ("BLCA", "atezolizumab + platinum chemotherapy"): "phase_3",
        ("KIRP", "pembrolizumab"): "phase_2",
        ("MESO", "bevacizumab + pemetrexed/cisplatin"): "phase_3",
        ("MTC", "pralsetinib"): "phase_2",
        ("THYM", "pembrolizumab"): "phase_2",
        ("THYM", "sunitinib"): "phase_2",
        ("FL", "tazemetostat"): "phase_2",
    }
    for (cancer_code, agent), expected_phase in corrected_phases.items():
        rows = key[(key["cancer_code"] == cancer_code) & (key["agent"] == agent)]
        assert not rows.empty, (cancer_code, agent)
        assert set(rows["phase"]) == {expected_phase}

    removed_claims = {
        ("LGG", "temozolomide"),
        ("NPC", "pembrolizumab"),
    }
    current_claims = set(zip(key["cancer_code"], key["agent"]))
    assert removed_claims.isdisjoint(current_claims)

    thca_ret = key[
        (key["cancer_code"] == "THCA")
        & (key["symbol"] == "RET")
        & (key["role"] == "target")
    ]
    assert thca_ret["indication"].str.contains("RET-fusion-positive").all()
    assert "PMID:37870969" not in set(thca_ret["source"])

    assert not (
        (key["symbol"] == "B4GALNT1") & key["rationale"].str.contains("GD2", case=False)
    ).any()


def test_type_gene_sets():
    sets = cg.cancer_type_gene_sets("PRAD")
    assert sets  # non-empty role->{ensembl:symbol}
    assert all(isinstance(v, dict) for v in sets.values())

    mmnst = cg.cancer_biomarker_genes("SARC", subtype="MMNST")
    assert {
        "TYR",
        "PMEL",
        "MLANA",
        "DCT",
        "MITF",
        "SOX10",
        "S100B",
        "PMP22",
        "PMP2",
        "MPZ",
        "PRKAR1A",
    } <= set(mmnst)


def test_narrative_and_rule_loaders():
    assert {"set_name", "members"} <= set(cg.narrative_gene_sets_df().columns)
    assert "rule_id" in cg.disease_state_rules_df().columns
    assert "pair_id" in cg.degenerate_subtype_pairs_df().columns
    assert {"Symbol", "Ensembl_Gene_ID"} <= set(cg.cancer_driver_genes_df().columns)


def test_fusion_rule_tables():
    assert {"cancer_code", "gene_a", "gene_b"} <= set(f.rare_cancer_fusion_rules_df().columns)
    assert {"fusion_class", "surrogate_gene"} <= set(f.fusion_surrogate_expression_df().columns)
    assert {"gene_a", "gene_b"} <= set(f.fusion_expression_effect_rules_df().columns)
    assert isinstance(f.fusion_surrogate_genes_for_cancer("SARC_EWS"), list)


def test_mcl_diagnostic_references_are_source_anchored():
    bad_pmids = {
        "PMID:9500537",
        "PMID:8049438",
        "PMID:20554603",
        "PMID:14983945",
        "PMID:12438234",
        "PMID:18832546",
    }

    key = cg.cancer_key_genes_df()
    small_b_cell_key = key[key["cancer_code"].isin(["CLL", "MCL"])]
    all_key_refs = set().union(*(_split_refs(v) for v in small_b_cell_key["source"]))
    assert bad_pmids.isdisjoint(all_key_refs)

    mcl_key = key[key["cancer_code"] == "MCL"]
    ccnd1 = mcl_key[mcl_key["symbol"] == "CCND1"]
    assert len(ccnd1) == 1
    ccnd1 = ccnd1.iloc[0]
    assert {"PMID:34114641", "PMID:40381701"} <= _split_refs(ccnd1["source"])
    assert "pathognomonic" not in str(ccnd1["rationale"]).lower()
    assert "cyclin D1 overexpression" in str(ccnd1["rationale"])

    sox11 = mcl_key[mcl_key["symbol"] == "SOX11"]
    assert len(sox11) == 1
    sox11 = sox11.iloc[0]
    assert sox11["source"] == "PMID:19801969"
    assert "~95%" in str(sox11["rationale"])

    cll_key = key[key["cancer_code"] == "CLL"]
    for symbol in ("CD5", "CD23"):
        row = cll_key[cll_key["symbol"] == symbol]
        assert len(row) == 1
        assert row.iloc[0]["source"] == "PMID:32249238"

    fusion = f.cancer_fusions("MCL")
    ccnd1_igh = fusion[fusion["fusion_family"] == "CCND1-IGH"]
    assert len(ccnd1_igh) == 1
    ccnd1_igh = ccnd1_igh.iloc[0]
    assert ccnd1_igh["pmid"] == "PMID:40381701"
    assert ccnd1_igh["frequency"] == "~95%"
    assert "unverified" not in str(ccnd1_igh["notes"]).lower()

    pair = cg.degenerate_subtype_pairs_df()
    pair = pair[pair["pair_id"] == "CLL_vs_MCL_vs_FL"]
    assert len(pair) == 1
    pair = pair.iloc[0]
    assert bad_pmids.isdisjoint(_split_refs(pair["refs"]))
    assert {"PMID:34114641", "PMID:40381701", "PMID:23897248", "PMID:32249238"} <= _split_refs(
        pair["refs"]
    )

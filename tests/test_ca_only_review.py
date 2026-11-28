"""Tests for the 66-protein CA-only functional biology & adult micro-niche review."""

from __future__ import annotations

import oncoref
from oncoref import ca_only_review, human, target_selection


def test_ca_only_functional_review_integrity_and_counts():
    meta = ca_only_review.ca_only_review_metadata()
    df = ca_only_review.ca_only_functional_review()
    assert len(df) == meta["total_reviewed_proteins"] == 66
    assert df["entity_id"].is_unique
    assert df["protein_group_label"].is_unique
    assert set(df["adjudication_disposition"]) == set(ca_only_review.CA_ONLY_DISPOSITIONS)
    assert set(df["functional_class"]) == set(ca_only_review.CA_ONLY_FUNCTIONAL_CLASSES)
    assert df["primary_citations"].str.len().gt(5).all()
    assert df["primary_urls"].str.startswith("https://").all()
    assert df["molecular_function_summary"].str.len().gt(25).all()

    counts = df["adjudication_disposition"].value_counts().to_dict()
    assert counts == {
        "retained_ca_only": 53,
        "adjudicated_unsafe_adult_human": 7,
        "adjudicated_review_niche_or_paralog": 6,
    }

    adj = human.ca_only_niche_adjudications()
    assert len(adj) == 13
    summary = human.ca_only_class_summary()
    assert int(summary["n_proteins"].sum()) == 66


def test_ca_only_adjudications_match_active_human_catalog():
    cat = target_selection.results().set_index("entity_id")
    df = oncoref.ca_only_review.ca_only_functional_review()

    for r in df.itertuples(index=False):
        got_cls = cat.loc[r.entity_id, "classification"]
        assert got_cls == r.post_adjudication_classification, r.protein_group_label
        assert cat.loc[r.entity_id, "ca_only_review_disposition"] == r.adjudication_disposition
        assert cat.loc[r.entity_id, "ca_only_functional_class"] == r.functional_class

    galp = ca_only_review.ca_only_functional_review("GALP").iloc[0]
    assert galp["post_adjudication_classification"] == "Unsafe"
    assert "39910309" in galp["primary_citations"]

    ocm2 = ca_only_review.ca_only_functional_review("OCM2").iloc[0]
    assert ocm2["post_adjudication_classification"] == "Other:Review"
    assert int(ocm2["safe_2mm_9mers"]) == 0

    csag3 = ca_only_review.ca_only_functional_review("CSAG3").iloc[0]
    assert csag3["post_adjudication_classification"] == "Cancer Antigen only"
    assert csag3["functional_class"] == "sub_10x_ctx_cancer_germline_paralog"

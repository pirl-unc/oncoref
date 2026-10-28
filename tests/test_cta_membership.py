"""Biological identity survives normal-expression concerns; evidence does not."""

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scripts.build_cta_membership import expression_facts, tumor_facts

from oncoref import cta_membership as policy
from oncoref import cta_research_policy

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/audits/cta-membership-20261026"


def facts(**changes):
    row = {
        "Gene ID": "g1",
        "Gene": "UNFAMILIAR",
        "Taxonomy ID": 9606,
        "Reference": "reference",
        "Core RNA (nTPM)": 4.0,
        "Testis RNA (nTPM)": 4.0,
        "Ovary RNA (nTPM)": 0.0,
        "Placenta RNA (nTPM)": 0.0,
        "Core comparator maximum (nTPM)": 1.0,
        "Core tissues measured": 3,
        "Core comparator tissues measured": 47,
        "Other reproductive RNA (nTPM)": 0.0,
        "Other reproductive tissues measured": 8,
        "Nonreproductive maximum (nTPM)": 1.0,
        "Nonreproductive tissues measured": 39,
        "Cancer-germline papers": 0,
        "Cancer-placental papers": 0,
        "Publication classes": "",
        "Tumor RNA positive": True,
        "Tumor IHC positive": False,
        "Tumor minimum samples": 3,
        "Tumor minimum prevalence": 0.05,
        "Tumor IHC minimum cases": 2,
        "Expression forbidden": True,
        "Antibody concern": True,
        "Normal IHC positive": True,
    }
    row.update(changes)
    return pd.DataFrame([row])


def review():
    records = policy.evidence()
    d = records.loc[records["Admission role"].eq("Membership support")].iloc[[0]].copy()
    d["Gene ID"], d["Gene"], d["Reference"] = "g1", "UNFAMILIAR", "reference"
    return d


def classify(f=None, reviewed=False, thresholds=None):
    return policy.classify(
        facts() if f is None else f,
        review() if reviewed else review().iloc[:0],
        thresholds=thresholds,
    ).iloc[0]


def test_established_identity_survives_all_normal_flags_and_no_bulk_expression():
    d = classify(
        facts(
            **{
                "Core RNA (nTPM)": 0,
                "Testis RNA (nTPM)": 0,
                "Core comparator maximum (nTPM)": 100000,
                "Tumor RNA positive": False,
                "Tumor IHC positive": False,
            }
        ),
        reviewed=True,
    )
    assert d.Member and d["Evidence tier"] == "Established"
    assert d["Expression forbidden"] and d["Normal IHC positive"] and d["Antibody concern"]


def test_familiar_gene_symbol_is_not_evidence_and_control_labels_are_not_inputs():
    assert not classify(facts(Gene="PRAME", **{"Tumor RNA positive": False})).Member
    assert classify(reviewed=True).Member


def test_context_only_primary_evidence_never_forces_admission():
    context = review()
    context["Admission role"] = "Context only"
    f = facts(**{"Tumor RNA positive": False})
    result = policy.classify(f, context).iloc[0]
    assert not result.Member and not result["Established route"]
    assert result["Membership sources"] == context.Sources.iloc[0]
    f["Cancer-germline papers"] = 2
    result = policy.classify(f, context).iloc[0]
    assert result.Member and result["Evidence tier"] == "Publication consensus"


def test_magea4_is_the_only_mage_control_and_other_mages_qualify_organically():
    controls = policy.controls()
    assert set(controls.loc[controls.Gene.str.startswith("MAGE"), "Gene"]) == {"MAGEA4"}
    assert {"CT83", "SSX2", "TBXT"} <= set(
        controls.loc[controls["Control role"].eq("User-required"), "Gene"]
    )
    reviews = policy.evidence()
    other = reviews.loc[reviews.Gene.str.startswith("MAGE") & reviews.Gene.ne("MAGEA4")]
    assert len(other) == 6 and other["Admission role"].eq("Context only").all()
    results = policy.results().set_index("Gene")
    assert not results.loc[other.Gene, "Established route"].any()
    assert results.loc[other.Gene, "Consensus route"].all()
    assert results.loc["MAGEA3", "Expression category"] == "Unsafe"


def test_source_consensus_can_survive_a_modern_atlas_non_detection():
    row = classify(
        facts(**{"Cancer-germline papers": 2, "Core RNA (nTPM)": 0, "Tumor RNA positive": False})
    )
    assert row.Member and row["Evidence tier"] == "Publication consensus"
    assert not row["Expression route"]


def test_expression_route_accepts_unlisted_genes_at_inclusive_boundary():
    assert classify().Member
    assert not classify(facts(**{"Core comparator maximum (nTPM)": 1.001})).Member
    assert not classify(facts(**{"Tumor RNA positive": False})).Member


@pytest.mark.parametrize(
    "field,value",
    [
        ("Core tissues measured", 2),
        ("Core comparator tissues measured", 46),
        ("Core comparator maximum (nTPM)", np.nan),
        ("Core RNA (nTPM)", 0.9),
    ],
)
def test_missing_or_weak_expression_never_means_absence(field, value):
    assert not classify(facts(**{field: value})).Member


def test_other_reproductive_tissues_need_explicit_cta_publication():
    d = facts(**{"Core RNA (nTPM)": 0.0, "Other reproductive RNA (nTPM)": 8.0})
    assert not classify(d).Member
    d["Cancer-germline papers"] = 1
    result = classify(d)
    assert result.Member and result.Classes == "Cancer-reproductive"
    d["Nonreproductive tissues measured"] = 38
    assert not classify(d).Member


def test_zero_normal_values_are_observed_and_positive_floor_still_applies():
    assert classify(facts(**{"Core comparator maximum (nTPM)": 0.0})).Member
    assert not classify(
        facts(**{"Core comparator maximum (nTPM)": 0.0, "Core RNA (nTPM)": 0.0})
    ).Member


@pytest.mark.parametrize(
    "field,value",
    [
        ("Class", "Unknown"),
        ("Reference", "different-reference"),
        ("Gene", "wrong-symbol"),
        ("Identity scope", "Unresolved family"),
        ("Biological evidence", ""),
        ("Tumor evidence", ""),
        ("Sources", "uncited"),
        ("Taxonomy ID", 10090),
    ],
)
def test_primary_evidence_must_be_reference_bound_and_explicit(field, value):
    r = review()
    r[field] = value
    with pytest.raises(ValueError):
        policy.classify(facts(), r)


@pytest.mark.parametrize(
    "field,value",
    [
        ("Core RNA (nTPM)", -1),
        ("Core RNA (nTPM)", np.inf),
        ("Cancer-germline papers", 0.5),
        ("Tumor RNA positive", "True"),
    ],
)
def test_invalid_facts_fail_closed(field, value):
    with pytest.raises(ValueError):
        classify(facts(**{field: value}))


def test_changed_tumor_cutoffs_require_reaggregation():
    with pytest.raises(ValueError, match="Reaggregate"):
        classify(thresholds=replace(policy.Thresholds(), tumor_samples=5))


def test_no_pooling_of_subthreshold_tumor_cohorts_or_repeated_papers():
    rna = pd.DataFrame(
        [
            {
                "gene_id": "g1",
                "cancer_code": code,
                "cohort": "TCGA",
                "samples": 50,
                "expressed_samples_ptpm_ge_1": 2,
                "prevalence_ptpm_ge_1": 0.04,
            }
            for code in ["A", "B"]
        ]
    )
    ihc = pd.DataFrame(
        [
            {
                "gene_id": "g1",
                "cancer": "A",
                "medium_high": 0,
                "measurement_status": "measured",
                "total": 10,
            }
        ]
    )
    result, _, _ = tumor_facts(facts(), rna, ihc, policy.Thresholds())
    assert not result["Tumor RNA positive"].iloc[0]
    rna.loc[0, "expressed_samples_ptpm_ge_1"] = 3
    rna.loc[0, "prevalence_ptpm_ge_1"] = 0.06
    result, _, _ = tumor_facts(facts(), rna, ihc, policy.Thresholds())
    assert result["Tumor RNA positive"].iloc[0]


def test_multiple_tables_of_one_doi_do_not_count_as_multiple_papers():
    base = cta_research_policy.results().iloc[[0]]
    gid = base["Gene ID"].iloc[0]
    names = [
        "testis",
        "ovary",
        "placenta",
        "breast",
        "cervix",
        "endometrium",
        "epididymis",
        "fallopian tube",
        "prostate",
        "seminal vesicle",
        "vagina",
    ] + [f"normal-{i}" for i in range(39)]
    normal = pd.DataFrame({"Gene": gid, "Tissue": names, "nTPM": 1.0})
    membership = pd.DataFrame(
        [
            {
                "Ensembl_Gene_ID": gid,
                "doi": "one-paper",
                "Class": "Cancer-germline",
                "Counts toward consensus": True,
                "source_tag": tag,
            }
            for tag in ["table1", "table2"]
        ]
    )
    result = expression_facts(base, normal, membership)
    assert result["Cancer-germline papers"].iloc[0] == 1


def test_published_reprocessing_and_benchmark_contract():
    manifest = json.loads((AUDIT / "run-manifest.json").read_text())
    for name, digest in manifest["inputs"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    for name, digest in manifest["outputs"].items():
        assert hashlib.sha256((AUDIT / name).read_bytes()).hexdigest() == digest, name
    shipped = policy.results()
    assert len(shipped) == 20097 and shipped["Gene ID"].is_unique
    benchmark = policy.controls().merge(shipped, on=["Gene ID", "Gene"], validate="one_to_one")
    assert benchmark.Member.eq(benchmark["Expected member"]).all()
    assert {"PRAME", "MAGEA4", "CTAG1A", "CTAG1B"} <= set(shipped.loc[shipped.Member, "Gene"])
    assert set(shipped.loc[shipped.Member & shipped.Classes.str.contains("Oncofetal"), "Gene"]) >= {
        "AFP",
        "GPC3",
        "SALL4",
        "TBXT",
        "CLDN6",
    }
    assert len(policy.gene_ids(established_only=True)) == 15
    assert policy.gene_ids(category="Oncofetal") <= policy.gene_ids()
    sensitivity = pd.read_csv(AUDIT / "threshold-sensitivity.csv")
    assert len(sensitivity) == 17
    assert sensitivity["Positive controls selected"].eq(15).all()
    assert sensitivity["Negative controls selected"].eq(0).all()
    facts_table = pd.read_csv(AUDIT / "decision-facts.csv.gz", low_memory=False)
    rebuilt = policy.classify(facts_table, policy.evidence())
    for col in [
        "Member",
        "Classes",
        "Evidence tier",
        "Established route",
        "Consensus route",
        "Expression route",
    ]:
        assert rebuilt[col].equals(shipped[col]), col
    old = cta_research_policy.results().set_index("Gene ID")
    assert shipped.set_index("Gene ID")["Expression forbidden"].equals(old["Expression forbidden"])


def test_ctag_duplicates_and_developmental_gaps_remain_explicit():
    d = policy.results().set_index("Gene")
    assert (
        d.loc[["CTAG1A", "CTAG1B"], "Identity scope"].eq("Shared antigen; locus unresolved").all()
    )
    assert "nucleus pulposus" in d.loc["TBXT", "Membership limitations"]
    assert d["Developmental atlas coverage"].str.startswith("Not imported").all()
    assert d["Normal breadth vital upper"].ge(d["Normal breadth vital lower"]).all()
    assert d["Normal breadth vital missing_weight"].gt(0).all()


def test_nonselection_and_expression_risk_counts_reconcile_without_collapsing_axes():
    d = policy.results()
    counts = pd.read_csv(AUDIT / "expression-status-counts.csv").set_index("Scope")
    assert counts[["CTA", "Unsafe", "Neither"]].sum(axis=1).equals(counts.Genes)
    assert counts.loc["Selected biological members", ["CTA", "Unsafe", "Neither"]].tolist() == [
        211,
        726,
        479,
    ]
    exclusions = pd.read_csv(AUDIT / "nonselection-reasons.csv.gz")
    assert exclusions["Gene ID"].is_unique
    assert set(exclusions["Gene ID"]) == set(d.loc[~d.Member, "Gene ID"])
    reasons = pd.read_csv(AUDIT / "nonselection-counts.csv")
    nominees = exclusions.loc[exclusions.Nominated]
    assert len(nominees) == 1229
    assert (
        nominees["Nonselection reason"].value_counts().to_dict()
        == reasons.loc[reasons.Scope.eq("Previously nominated nonmembers")]
        .set_index("Reason")
        .Genes.to_dict()
    )

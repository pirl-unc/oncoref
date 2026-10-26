"""Decision precedence, evidence identity, faithful migration and current artifacts."""

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from oncoref import cta_policy as policy
from oncoref import target_sets

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/audits/cta-policy-20261025"


def facts(**overrides):
    row = dict.fromkeys(policy.FLAGS, False)
    row.update(
        {
            "Gene": "A",
            "Gene ID": "g1",
            "Taxonomy ID": 9606,
            "Reference": "test-reference",
            "Reproductive RNA (nTPM)": 10.0,
            "Nominated": True,
            "Restriction passed": True,
            "Coverage adequate": True,
        }
    )
    row.update(overrides)
    return pd.DataFrame([row])


def review(kind="Normal protein", status="Accepted", **overrides):
    row = dict.fromkeys(policy.REVIEW_COLUMNS, "")
    row.update(
        {
            "Review ID": "study-1",
            "Gene ID": "g1",
            "Gene": "A",
            "Taxonomy ID": 9606,
            "Reference": "test-reference",
            "Finding": kind,
            "Status": status,
            "Tissue category": "vital",
            "Context": "adult normal human neurons",
            "Identity scope": "gene",
            "Assay": "endogenous protein",
            "Specificity": "orthogonal identification",
            "Source": "https://example.org/study",
            "Study": "primary-study-1",
            "Reviewed by": "test",
            "Review date": "2026-10-25",
            "Review depth": "Full text",
            "Rationale": "Reported endogenous normal protein",
        }
    )
    if kind.startswith("Reproductive"):
        row["Tissue category"] = "reproductive"
    row.update(overrides)
    return row


def classify(f=None, rs=()):
    return policy.classify(
        facts() if f is None else f, pd.DataFrame(rs, columns=policy.REVIEW_COLUMNS)
    ).iloc[0]


def test_normal_protein_overrides_reproductive_abundance_and_support():
    r = classify(
        facts(**{"Reproductive RNA (nTPM)": 10000, "Reproductive protein": True}), [review()]
    )
    assert r.Category == "Unsafe" and r.Forbidden and not r.Target
    assert r.Detail == "Normal"


@pytest.mark.parametrize(
    "status,expected", [("Accepted", "Unsafe"), ("Unresolved", "Neither"), ("Rejected", "CTA")]
)
def test_finding_acceptance_is_explicit(status, expected):
    r = classify(rs=[review(status=status)])
    assert r.Category == expected
    assert bool(r["Accepted normal protein"]) == (status == "Accepted")


def test_unresolved_does_not_hide_decisive_exclusion_and_all_reasons_survive():
    r = classify(
        facts(**{"High normal RNA": True, "Antibody concern": True}), [review(status="Unresolved")]
    )
    assert r.Category == "Unsafe" and r.Detail == "High RNA"
    assert "30 TPM" in r["All reasons"]
    assert "unresolved" in r["All reasons"]
    assert "antibody" in r["All reasons"]


def test_protein_and_reviewed_biology_bypass_floor_but_not_normal_gates():
    low = facts(**{"Reproductive RNA (nTPM)": 1.2})
    assert classify(low).Category == "Neither"
    assert classify(low, [review("Reproductive protein")]).Category == "CTA"
    assert classify(low, [review("Reproductive expression")]).Category == "Neither"
    assert classify(low, [review("Tumor recognition")]).Category == "Neither"
    support = [
        review("Reproductive expression"),
        review("Tumor recognition", **{"Review ID": "recognition"}),
    ]
    assert classify(low, support).Category == "CTA"
    low["Normal IHC positive"] = True
    assert classify(low, support).Category == "Unsafe"


@pytest.mark.parametrize("veto", ["High normal RNA", "Normal IHC positive", "Safety veto"])
def test_variable_is_narrow_and_cannot_override_veto(veto):
    f = facts(**{"Variable": True, "Normal RNA failure": True, "Restriction passed": False})
    assert classify(f).Category == "Neither"
    f[veto] = True
    assert classify(f).Category == "Unsafe"
    assert classify(facts(**{"Variable": True}), [review()]).Category == "Unsafe"


@pytest.mark.parametrize("flag", ["Coverage adequate", "Restriction passed", "Nominated"])
def test_inadequacy_is_neither_not_forbidden(flag):
    r = classify(facts(**{flag: False}))
    assert r.Category == "Neither" and not r.Forbidden


@pytest.mark.parametrize("value", [None, "False", 0])
def test_no_missing_or_string_boolean_inference(value):
    with pytest.raises(ValueError, match="explicit booleans"):
        classify(facts(**{"Normal IHC positive": value}))


@pytest.mark.parametrize(
    "overrides",
    [
        {"Taxonomy ID": 10090},
        {"Reference": "another"},
        {"Identity scope": "family"},
        {"Identity scope": "exact protein"},
        {"Tissue category": "disposable"},
        {"Finding": "Trial registration"},
        {"Specificity": ""},
    ],
)
def test_review_identity_and_scope_are_not_assumed(overrides):
    with pytest.raises(ValueError):
        classify(rs=[review(**overrides)])


def test_no_gene_name_branches_or_paper_voting():
    a = classify(rs=[review()])
    f = facts(Gene="FAMOUS_CTA")
    b = classify(
        f,
        [
            review(Gene="FAMOUS_CTA"),
            review(Gene="FAMOUS_CTA", **{"Review ID": "duplicate-study-observation"}),
        ],
    )
    assert (a.Category, a.Detail) == (b.Category, b.Detail)
    with pytest.raises(ValueError, match="Unique review"):
        classify(rs=[review(), review()])


@pytest.mark.parametrize("value", [0, float("nan")])
def test_inconsistent_measurement_facts_cannot_make_cta(value):
    with pytest.raises(ValueError, match="positive reproductive RNA"):
        classify(facts(**{"Reproductive RNA (nTPM)": value, "Reproductive protein": True}))


def test_frozen_measurement_replay_and_adjudicated_changes():
    f = pd.read_csv(REPORT / "decision-facts.csv.gz")
    rs = policy.findings()
    original = pd.read_csv(
        ROOT / "docs/audits/cta-neither-review-20261024/gene-results.csv.gz", low_memory=False
    )
    replay = policy.classify(f, rs.loc[rs["Review ID"].str.startswith("legacy-")])
    assert replay.Category.tolist() == original.Category.tolist()
    current = policy.classify(f, rs)
    assert current.Category.tolist() == target_sets.results().Category.tolist()
    changes = current.loc[current.Category.ne(original.Category)]
    initial = {"FTMT", "ABCC11", "YBX2", "CCDC158", "ACTRT1", "MAGEC2", "LGALS13", "PLCZ1"}
    cases = policy.adjudications()
    assert set(changes.Gene) == initial | set(cases.Gene)
    assert len(changes) == 186
    assert changes.Category.value_counts().to_dict() == {"Neither": 180, "Unsafe": 6}
    assert set(original.loc[changes.index, "Category"]) == {"CTA"}


def test_current_policy_package_and_report_hashes():
    run = json.loads((REPORT / "run-manifest.json").read_text())
    for section, folder in [("implementation_and_inputs", ROOT), ("outputs", REPORT)]:
        for name, expected in run[section].items():
            assert hashlib.sha256((folder / name).read_bytes()).hexdigest() == expected, name
    for packaged, artifact in [
        ("target-gene-sets.csv.gz", "gene-results.csv.gz"),
        ("cta-followup-reviews.csv", "follow-up.csv"),
        ("sequence-gene-rna-reviews.csv.gz", "sequence-gene-reviews.csv.gz"),
        ("cta-policy.json", "policy.json"),
    ]:
        assert (ROOT / "oncoref/data" / packaged).read_bytes() == (REPORT / artifact).read_bytes()
    m = policy.manifest()
    assert m["policy"] == policy.POLICY_ID
    for name, expected in m["data_sha256"].items():
        assert hashlib.sha256((ROOT / "oncoref/data" / name).read_bytes()).hexdigest() == expected
    rows = target_sets.results()
    assert set(rows.Policy) == {policy.POLICY_ID}
    assert rows.Category.value_counts().to_dict() == m["counts"]["genome"]
    assert rows.loc[rows.Nominated].Category.value_counts().to_dict() == m["counts"]["nominated"]
    assert m["counts"]["proteins"] == {"CTA": 359, "Neither": 3230, "Unsafe": 126396}
    recipe = json.loads((REPORT / "exclusion-basis.json").read_text())
    assert recipe["policy"] == policy.POLICY_ID
    assert recipe["forbidden"] == ["unsafe.fasta.gz", "additional-forbidden.fasta.gz"]
    assert not recipe["neither_forbidden"]


@pytest.mark.parametrize("status,expected", [("Unresolved", "Neither"), ("Rejected", "CTA")])
def test_qualitative_normal_rna_is_not_protein(status, expected):
    r = classify(rs=[review("Normal RNA", status=status)])
    assert r.Category == expected
    assert not r["Accepted normal protein"]
    with pytest.raises(ValueError, match="Qualitative Normal RNA"):
        classify(rs=[review("Normal RNA")])


def test_all_remaining_raw_antibody_cases_have_complete_evidence_derived_dispositions():
    cases = policy.adjudications()
    assert len(cases) == 178 and cases["Review complete"].all()
    assert not cases["Gene ID"].duplicated().any()
    assert cases.Decision.value_counts().to_dict() == {"Neither": 177, "Unsafe": 1}
    current = target_sets.results()
    merged = cases.merge(current, on="Gene ID", validate="one_to_one")
    assert merged.Decision.eq(merged.Category).all()
    assert not merged.Target.any()
    raw = pd.read_csv(
        ROOT / "docs/audits/cta-literature-20261025/retained-primary-normal-signals.csv"
    )
    assert not set(raw.gene_id) & set(current.loc[current.Target, "Gene ID"])
    assert policy.manifest()["counts"]["remaining_primary_ihc_queue"] == 0
    rs = policy.findings()
    plcz1 = rs.loc[rs.Gene.eq("PLCZ1")]
    assert plcz1.loc[plcz1.Finding.eq("Normal protein"), "Status"].eq("Rejected").all()
    assert plcz1.loc[plcz1.Finding.eq("Normal RNA"), "Status"].eq("Unresolved").all()
    run = json.loads((ROOT / "docs/audits/cta-adjudication-20261025/run-manifest.json").read_text())
    for section, root in [
        ("inputs", ROOT),
        ("outputs", ROOT / "docs/audits/cta-adjudication-20261025"),
    ]:
        for name, expected in run[section].items():
            assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected

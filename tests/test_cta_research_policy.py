import pandas as pd
import pytest

from oncoref import cta_policy
from oncoref import cta_research_policy as policy
from tests.test_cta_policy import facts, review


def clinical(**changes):
    row = policy.clinical_evidence().iloc[0].to_dict()
    row.update(
        {
            "Gene": "A",
            "Gene ID": "g1",
            "Reference": "test-reference",
            "Reviewed finding IDs": "study-1",
        }
    )
    row.update(changes)
    return pd.DataFrame([row])


def classify(f=None, c=None, rs=None):
    reviewed = pd.DataFrame([review()] if rs is None else rs, columns=cta_policy.REVIEW_COLUMNS)
    evidence = clinical() if c is None else c.copy()
    evidence["Findings SHA256"] = policy.finding_digest(reviewed, "g1")
    return policy.classify(
        facts() if f is None else f,
        reviewed,
        evidence,
    ).iloc[0]


def test_clinical_route_preserves_true_normal_protein_and_background_scope():
    d = classify()
    assert d.Category == "CTA" and d.Target
    assert d["Accepted normal protein"] and d["Expression category"] == "Unsafe"
    assert d["Expression forbidden"]
    assert "Forbidden" not in d.index
    assert d["Selection basis"] == "Clinical qualification with normal-expression caveats"


@pytest.mark.parametrize(
    "veto",
    [
        "Safety veto",
        "High normal RNA",
        "Normal IHC positive",
        "Normal RNA failure",
        "Variable",
        "Cell RNA concern",
        "Antibody concern",
    ],
)
def test_clinical_qualification_cannot_erase_measured_vetoes(veto):
    assert not classify(facts(**{veto: True})).Target


@pytest.mark.parametrize("field", policy.BOOLEAN_FIELDS[:-1])
def test_registration_or_incomplete_clinical_evidence_cannot_qualify(field):
    assert not classify(c=clinical(**{field: False})).Target


def test_substantial_normal_protein_and_unresolved_organ_toxicity_prevent_route():
    assert not classify(c=clinical(**{"Normal finding scope": "Substantial"})).Target
    assert not classify(c=clinical(**{"Unresolved target-mediated organ toxicity": True})).Target


def test_clinical_organ_toxicity_also_holds_an_expression_qualified_gene():
    d = classify(
        c=clinical(
            **{
                "Reviewed finding IDs": "None",
                "Normal finding scope": "None",
                "Unresolved target-mediated organ toxicity": True,
            }
        ),
        rs=[],
    )
    assert d["Expression category"] == "CTA"
    assert not d.Target and d.Category == "Neither"
    assert not d["Expression forbidden"]


def test_unresolved_normal_findings_cannot_be_labelled_absent():
    with pytest.raises(ValueError, match="Active normal findings"):
        classify(c=clinical(**{"Normal finding scope": "None"}), rs=[review(status="Unresolved")])


def test_clinical_protein_reconsideration_does_not_relax_a_normal_rna_hold():
    assert not classify(rs=[review(kind="Normal RNA", status="Unresolved")]).Target


def test_no_famous_gene_exception_or_reproductive_maximum_rescue():
    a = classify()
    b = classify(
        facts(Gene="UNFAMILIAR", **{"Reproductive RNA (nTPM)": 50000}),
        clinical(Gene="UNFAMILIAR"),
        [review(Gene="UNFAMILIAR")],
    )
    assert a.Category == b.Category
    assert not classify(facts(**{"High normal RNA": True, "Reproductive RNA (nTPM)": 50000})).Target


def test_new_findings_require_explicit_re_review():
    rs = [review(), review(status="Unresolved", **{"Review ID": "new-conflict"})]
    with pytest.raises(ValueError, match="every active normal finding"):
        classify(rs=rs)
    with pytest.raises(ValueError, match="Accepted normal protein"):
        classify(c=clinical(**{"Normal finding scope": "None"}))


def test_changed_finding_contents_invalidate_clinical_review():
    rs = pd.DataFrame([review()], columns=cta_policy.REVIEW_COLUMNS)
    evidence = clinical()
    evidence["Findings SHA256"] = policy.finding_digest(rs, "g1")
    rs.loc[0, "Rationale"] = "New independent evidence of substantial normal protein"
    with pytest.raises(ValueError, match="contents changed"):
        policy.classify(facts(), rs, evidence)


@pytest.mark.parametrize(
    "field,value",
    [
        ("Taxonomy ID", 10090),
        ("Reference", "wrong"),
        ("Normal selectivity assessed", "False"),
        ("Evidence sources", "uncited"),
    ],
)
def test_reference_binding_and_no_string_truthiness(field, value):
    with pytest.raises(ValueError):
        classify(c=clinical(**{field: value}))

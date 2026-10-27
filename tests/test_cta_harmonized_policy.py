"""Unusable antibodies lower confidence without adding a second gene-wide veto."""

import pandas as pd
from cta_harmonized_report import ASSAYS, decision_layers

from oncoref.cta_funnel import classify_layers


def test_antibody_flag_alone_does_not_override_strict_rna_but_positive_evidence_does():
    facts = pd.DataFrame(
        {
            "gene_id": ["strict_rna_pass", "rna_failure", "protein_positive"],
            "positive_context": True,
            "tumor_union_ge2": True,
            "required_coverage": True,
            "ihc_detected": [False, False, True],
            "ihc_reliable_detected": [False, False, True],
            "Antibody concern": [True, True, False],
            "Accepted normal protein": False,
            "Safety veto": False,
            "external_eligible_positive": False,
        }
    )
    for assay in ASSAYS:
        facts["hold:" + assay] = False
        facts["unsafe:" + assay] = False
    facts.loc[1, ["hold:HPA tissue clusters", "unsafe:HPA tissue clusters"]] = True
    result = classify_layers(facts, decision_layers(facts)).classifications
    assert result.category.tolist() == ["CTA", "Unsafe", "Unsafe"]
    assert result["Antibody concern"].tolist() == [True, True, False]
    assert "hold:Unresolved antibody review" not in result

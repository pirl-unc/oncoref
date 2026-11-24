"""The comparison must preserve assay context and complete protein identities."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))
from audit_reproductive_evidence import group_values, scope


def test_group_requires_every_locus_in_the_same_cell_context():
    raw = pd.DataFrame(
        {
            "Gene": ["a", "b", "a", "b"],
            "Tissue": ["testis"] * 4,
            "Cluster": ["x", "x", "y", "y"],
            "nTPM": [3, 5, 100, np.nan],
        }
    )
    result = group_values(
        raw,
        {"a": "g1", "b": "g2"},
        {"g1": "protein", "g2": "protein"},
        {"protein": 2},
        ["Tissue", "Cluster"],
    )
    assert result[["Cluster", "value"]].to_dict("records") == [{"Cluster": "x", "value": 8}]
    with pytest.raises(ValueError, match="Duplicate"):
        group_values(
            pd.concat([raw, raw.iloc[:1]]),
            {"a": "g1", "b": "g2"},
            {"g1": "protein", "g2": "protein"},
            {"protein": 2},
            ["Tissue", "Cluster"],
        )


@pytest.mark.parametrize("value", [-1, float("inf")])
def test_invalid_member_cannot_be_hidden_by_aggregation(value):
    raw = pd.DataFrame({"Gene": ["a", "b"], "Tissue": ["testis"] * 2, "nTPM": [value, 10]})
    with pytest.raises(ValueError, match="Invalid source expression"):
        group_values(raw, {"a": "a", "b": "b"}, {"a": "p", "b": "p"}, {"p": 2}, ["Tissue"])


@pytest.mark.parametrize("tissue", ["breast", "prostate", "testis", "ovary"])
def test_reproductive_scope_keeps_accessory_tissues(tissue):
    assert scope(tissue) == "reproductive"


@pytest.mark.parametrize("tissue", ["retina", "pbmc", "cerebral cortex"])
def test_normal_scope_keeps_vital_cells_and_blood(tissue):
    assert scope(tissue) == "normal"


def test_ocular_paper_is_counterevidence_not_antigen_support():
    from oncoref import cta_literature

    data = cta_literature.publication_registry()
    paper = next(p for p in data["sources"] if p["source_id"] == "10.1167/iovs.09-3622")
    (claim,) = paper["gene_claims"]
    assert claim["gene_ids"] == ["ENSG00000124467"]
    assert claim["allowed_uses"] == ["normal_tissue_review"]
    assert "Probe specificity" in claim["limitations"]


def test_review_holds_remain_in_default_subtraction_background():
    from oncoref import antigen_sets, target_selection

    root = Path(__file__).parents[1]
    cohort = pd.read_csv(
        root / "docs/audits/antigen-evidence-review-20261118/requested-cohorts.csv"
    )
    holds = set(cohort.loc[cohort.baseline_classification.eq("Other:Review"), "entity_id"])
    assert len(holds) == 31
    live = target_selection.results().set_index("entity_id")
    # CTAG2 (#766) and DUX4 (#765 + #767) resolved to CTA; the remaining 29 stay in Unsafe|Other
    resolved_to_cta = set(
        live.loc[live.index.isin(holds) & live["category"].eq("CTA"), "protein_group_label"]
    )
    assert resolved_to_cta <= {"CTAG2", "DUX4", "TENT5D"}
    still_held = {eid for eid in holds if live.at[eid, "category"] in {"Other", "Unsafe"}}
    selection = antigen_sets.select("CA|RA")
    assert not still_held & set(selection.proteins().entity_id)
    assert still_held <= set(selection.background_proteins().entity_id)

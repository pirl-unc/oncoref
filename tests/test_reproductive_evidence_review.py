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

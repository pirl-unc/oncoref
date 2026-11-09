# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Rare vital-tissue cell expression: scope, 1% boundary and the cascade layer."""

import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

from oncoref import cta_rare_cell
from oncoref import cta_selection_policy as policy
from oncoref.cta_funnel import EvidenceLayer

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/build_cta_rare_cell_screen.py"
spec = importlib.util.spec_from_file_location("rare_cell_screen_script", SCRIPT)
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)


def raw(rows):
    return pd.DataFrame(
        rows,
        columns=["gene_id", "tissue", "cell_type", "n_expressing", "n_cells", "tissue_cells"],
    ).assign(mean_expr_log=1.6)


def kept(rows):
    return screen.eligible(raw(rows)).cell_type.tolist()


def test_vital_tissues_are_screened_and_other_tissues_are_allowed():
    rows = [
        ("g", "brain", "midbrain dopaminergic neuron", 293, 19760, 1_000_000),
        ("g", "heart", "cardiac muscle myocyte", 100, 8000, 500_000),
        ("g", "spinal cord", "motor neuron", 100, 8000, 500_000),
        ("g", "spleen", "follicular B cell", 400, 9000, 500_000),
        ("g", "adipose tissue", "subcutaneous adipocyte", 400, 9000, 500_000),
        ("g", "liver", "hepatocyte", 400, 9000, 500_000),
    ]
    assert kept(rows) == [
        "midbrain dopaminergic neuron",
        "cardiac muscle myocyte",
        "motor neuron",
    ]


def test_lung_tiers_follow_the_shared_cell_type_table():
    rows = [
        ("g", "lung", "type II pneumocyte", 200, 9000, 800_000),
        ("g", "lung", "pulmonary capillary endothelial cell", 200, 9000, 800_000),
        ("g", "lung", "basal cell", 200, 9000, 800_000),
        ("g", "lung", "respiratory ciliated cell", 200, 9000, 800_000),
        ("g", "lung", "club cell", 200, 9000, 800_000),
        ("g", "lung", "alveolar macrophage", 200, 9000, 800_000),
        ("g", "lung", "lung fibroblast", 200, 9000, 800_000),
        ("g", "lung", "CD8-positive T cell", 200, 9000, 800_000),
    ]
    table = screen.eligible(raw(rows)).set_index("cell_type").tier
    assert table.to_dict() == {
        "type II pneumocyte": "vital",
        "pulmonary capillary endothelial cell": "vital",
        "basal cell": "vital",
        "respiratory ciliated cell": "protected",
        "club cell": "vital",
    }


def test_airway_tier_is_shared_with_the_rna_cascade():
    contexts = pd.DataFrame(
        {
            "tissue": ["lung", "bronchus", "bronchus", "lung", "endometrium"],
            "cell_type": [
                "ciliated cells",
                "club cells",
                "basal respiratory cells",
                "alveolar cells type 2",
                "ciliated cells",
            ],
            "category": ["vital"] * 4 + ["reproductive"],
            "mapping_status": ["mapped"] * 5,
            "source": ["hpa_single_cell_tissue"] * 5,
        }
    )
    assert policy.compartment_tiers(contexts).tolist() == [
        "protected",
        "vital",
        "vital",
        "vital",
        "reproductive",
    ]


def test_immune_stromal_and_developmental_cell_types_are_ignored_in_vital_tissues():
    rows = [
        ("g", "brain", "microglial cell", 300, 9000, 1_000_000),
        ("g", "brain", "neuroblast (sensu Vertebrata)", 300, 9000, 1_000_000),
        ("g", "brain", "cerebellar granule cell precursor", 300, 9000, 1_000_000),
        ("g", "brain", "radial glial cell", 300, 9000, 1_000_000),
        ("g", "heart", "fibroblast of cardiac tissue", 300, 9000, 500_000),
        ("g", "heart", "cardiac macrophage", 300, 9000, 500_000),
        ("g", "brain", "malignant cell", 300, 9000, 1_000_000),
        ("g", "brain", "Purkinje cell", 300, 9000, 1_000_000),
    ]
    assert kept(rows) == ["Purkinje cell"]


def test_size_floors_and_broad_parent_terms_do_not_count():
    rows = [
        ("g", "brain", "small population", 50, 2999, 1_000_000),  # under 3,000 cells
        ("g", "brain", "few detections", 9, 9000, 1_000_000),  # under 10 expressing
        ("g", "brain", "broad parent population", 5000, 200_000, 1_000_000),  # over 10%
        ("g", "brain", "neuron", 500, 9000, 1_000_000),  # generic parent term
        ("g", "brain", "L5 pyramidal neuron", 500, 9000, 1_000_000),
    ]
    assert kept(rows) == ["L5 pyramidal neuron"]


def test_threshold_is_inclusive_at_one_percent(tmp_path, monkeypatch):
    rows = [
        ("exact", "brain", "astrocyte of the cortex", 100, 10_000, 1_000_000),  # 1.00%
        ("under", "brain", "astrocyte of the cortex", 99, 10_000, 1_000_000),  # 0.99%
    ]
    eligible = screen.eligible(raw(rows))
    assert dict(zip(eligible.gene_id, eligible.pct_expressing)) == pytest.approx(
        {"exact": 1.0, "under": 0.99}
    )
    assert screen.FLAG_PCT == cta_rare_cell.FLAG_PCT == 1.0


@pytest.fixture
def fake_screen(monkeypatch):
    table = pd.DataFrame(
        {
            "gene_id": ["A", "B", "B", "C"],
            "tissue": ["brain", "heart", "brain", "brain"],
            "cell_type": ["neuron x", "cardiac y", "neuron z", "neuron q"],
            "n_expressing": [150, 120, 30, 99],
            "n_cells": [10_000, 10_000, 10_000, 10_000],
            "pct_expressing": [1.5, 1.2, 0.3, 0.99],
            "tier": ["vital"] * 4,
        }
    )
    monkeypatch.setattr(cta_rare_cell, "screen_table", lambda: table)
    monkeypatch.setattr(
        cta_rare_cell, "provenance", lambda: {"screened_gene_ids": ["A", "B", "C", "D"]}
    )


def test_a_group_is_flagged_when_any_member_gene_reaches_the_threshold(fake_screen):
    facts = pd.DataFrame(
        {
            "member_gene_ids": [
                json.dumps(ids) for ids in (["A"], ["C"], ["C", "B"], ["D"], ["NOT_IN_SERVICE"])
            ]
        }
    )
    out = cta_rare_cell.annotate(facts)
    assert out.rare_cell_flag.tolist() == [True, False, True, False, False]
    # A multi-member group reports its highest member.
    assert out.rare_cell_type.tolist()[2] == "heart: cardiac y"
    assert out.rare_cell_max_pct.tolist()[:3] == pytest.approx([1.5, 0.99, 1.2])
    # Screened but clean is distinct from absent from the service.
    assert out.rare_cell_screened.tolist() == [True, True, True, True, False]
    assert out.rare_cell_type.tolist()[3] == out.rare_cell_type.tolist()[4] == ""
    assert out.rare_cell_max_pct.iloc[3:].isna().all()


def test_a_flag_makes_an_otherwise_passing_protein_unsafe_not_other():
    facts = pd.concat([_passing()] * 2, ignore_index=True)
    facts["entity"] = ["clean", "rare"]
    flag = pd.Series([False, True])
    no = pd.Series([False, False])
    layers = [
        EvidenceLayer("Reproductive / embryonic association", no, no),
        EvidenceLayer("Rare vital-tissue cell expression", flag, flag),
    ]
    result = policy.classify(facts, layers, entity_column="entity")
    assert result.classifications.classification.tolist() == ["CTA", "Unsafe"]
    effects = result.independent_effects.set_index("layer")
    assert effects.loc["Rare vital-tissue cell expression", "Unsafe_evidence_anywhere"] == 1


def test_the_layer_can_be_a_review_hold_instead_of_unsafe():
    facts = pd.concat([_passing()] * 2, ignore_index=True)
    facts["entity"] = ["clean", "rare"]
    flag, no = pd.Series([False, True]), pd.Series([False, False])
    layers = [
        EvidenceLayer("Reproductive / embryonic association", no, no),
        EvidenceLayer("Rare vital-tissue cell expression", flag, no),
    ]
    result = policy.classify(facts, layers, entity_column="entity").classifications
    assert result.classification.tolist() == ["CTA", "Other:Review"]


def test_protected_airway_expression_creates_hold_without_vital_flag(monkeypatch):
    table = pd.DataFrame(
        {
            "gene_id": ["A", "B"],
            "tissue": ["lung", "lung"],
            "cell_type": ["ciliated cells", "goblet cells"],
            "n_expressing": [100, 99],
            "n_cells": [10000, 10000],
            "pct_expressing": [1.0, 0.99],
            "tier": ["protected", "protected"],
        }
    )
    monkeypatch.setattr(cta_rare_cell, "screen_table", lambda: table)
    monkeypatch.setattr(cta_rare_cell, "provenance", lambda: {"screened_gene_ids": ["A", "B"]})
    out = cta_rare_cell.annotate(pd.DataFrame({"member_gene_ids": ['["A"]', '["B"]']}))
    assert out.rare_cell_hold.tolist() == [True, False]
    assert not out.rare_cell_flag.any()


def _passing():
    from tests.test_cta_selection_policy import profiles

    facts = policy.other_profiles(profiles().iloc[:1])
    facts["positive_context"] = True
    return facts

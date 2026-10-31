"""Source- and identity-bound positive RNA routes cannot borrow normal clearance."""

import pandas as pd
import pytest

from oncoref import cta_selection_policy as policy
from oncoref import reproductive_evidence as evidence


def test_full_reproductive_reference_replaces_stale_input_without_imputing_missing():
    facts = pd.DataFrame(
        {"protein_equivalence_group_id": ["a", "b"], "reproductive_ntpm": [4.3, 9]}
    )
    scores = pd.DataFrame(
        {
            "entity_id": ["a", "b"],
            "reference_value": [20.5, float("nan")],
            "source": "hpa_rna_consensus",
            "version": "v23",
            "unit": "nTPM",
            "taxonomy_id": 9606,
            "summary": "consensus",
        }
    )
    result = evidence.bulk_reference(facts, scores)
    assert result.reproductive_ntpm.iloc[0] == 20.5
    assert pd.isna(result.reproductive_ntpm.iloc[1])
    assert facts.reproductive_ntpm.tolist() == [4.3, 9]
    with pytest.raises(ValueError, match="HPA v23"):
        evidence.bulk_reference(facts, scores.assign(source="gtex_tissue_rna"))
    with pytest.raises(ValueError, match="per protein"):
        evidence.bulk_reference(facts, scores.iloc[:1])


def test_cell_presence_requires_support_identity_and_species_but_generic_sites_are_recorded():
    base = dict(
        source="hpa_single_cell_tissue",
        version="v23",
        taxonomy_id=9606,
        unit="nTPM",
        mapping_status="mapped",
        category="reproductive",
        tissue="testis",
        cell_type="spermatogonia",
        complete_members=True,
        member_count=1,
        value=5,
        **{"Read count": 20},
    )
    cases = [
        base,
        base | {"value": 4.9},
        base | {"Read count": 19},
        base | {"complete_members": False, "member_count": 2},
        base | {"taxonomy_id": 10090},
        base | {"version": "v24"},
        base | {"cell_type": "endothelial cells"},
        base
        | {
            "tissue": "breast",
            "cell_type": "macrophages",
            "category": "neutral_reproductive_immune",
        },
        base | {"tissue": "brain", "cell_type": "macrophages", "category": "vital"},
    ]
    result = evidence.cell_observations(pd.DataFrame(cases))
    assert result.reproductive_cell_positive.tolist() == [True] + [False] * 8
    assert result.reproductive_cell_indirect_candidate.tolist() == [False] * 6 + [True, True, False]
    assert result.reproductive_cell_reason.iloc[3] == "incomplete_members"


def test_supported_cell_route_is_independent_of_bulk_rna_floor():
    facts = pd.DataFrame(
        {
            "reproductive_ntpm": [1.2, 1.2],
            "developmental_mean_tpm": [0, 0],
            "reproductive_cell_association": [True, False],
        }
    )
    result = policy.association_evidence(facts)
    assert result.positive_context.tolist() == [True, False]
    assert not result.reproductive_rna_association.any()


def test_indirect_site_requires_complete_low_outside_panel_and_no_induction():
    base = {
        "reproductive_cell_indirect_ntpm": 5,
        "outside_cluster_complete": True,
        "outside_cluster_max_ntpm": 4.9,
        "normal_induction_evidence": False,
    }
    cases = [
        base,
        base | {"outside_cluster_complete": False},
        base | {"outside_cluster_max_ntpm": 5},
        base | {"outside_cluster_max_ntpm": float("nan")},
        base | {"normal_induction_evidence": True},
        base | {"normal_induction_evidence": None},
        base | {"reproductive_cell_indirect_ntpm": 4.9},
    ]
    assert evidence.indirect_association(pd.DataFrame(cases)).tolist() == [True] + [False] * 6

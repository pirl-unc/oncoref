"""Expression identities, same-context sums and unpaired-specimen bounds."""

import itertools
import json

import numpy as np
import pandas as pd
import pytest

from oncoref.protein_equivalence import (
    equivalence_groups,
    group_identifier,
    sum_context_observations,
    unpaired_rank_bounds,
)

CTAG = ["ENSG00000268651", "ENSG00000184033"]


def test_group_identity_is_membership_based_and_keeps_member_loci():
    groups = equivalence_groups(CTAG)
    assert groups.protein_group_label.tolist() == ["NY-ESO-1"]
    assert group_identifier(CTAG) == group_identifier(reversed(CTAG))
    assert set(json.loads(groups.member_gene_ids.iloc[0])) == set(CTAG)
    assert len(json.loads(groups.member_identities.iloc[0])) == 2
    assert groups.protein_equivalence_group_id.iloc[0].startswith("PEG:")


def test_context_sums_do_not_add_different_samples_or_turn_absence_into_zero():
    data = pd.DataFrame(
        {
            "gene_id": CTAG * 2,
            "sample": ["a", "a", "b", "b"],
            "value": [1.2, 1.3, np.nan, 0.5],
            "reads": [10, 11, np.nan, 2],
        }
    )
    result = sum_context_observations(
        data,
        equivalence_groups(CTAG),
        context_columns=["sample"],
        value_columns=["value", "reads"],
        unit="TPM",
    )
    assert result.value.tolist() == [2.5, 0.5]
    assert result.reads.tolist() == [21, 2]
    assert result.complete_members.tolist() == [True, False]
    with pytest.raises(ValueError, match="Repeated"):
        sum_context_observations(
            pd.concat([data, data.iloc[[0]]]),
            equivalence_groups(CTAG),
            context_columns=["sample"],
            value_columns=["value"],
            unit="TPM",
        )


def test_unpaired_rank_bounds_contain_every_possible_pairing():
    a, b = [10.0, 4.0, 1.0, 0.1], [6.0, 3.0, 2.0, 0.0]
    lower, upper = unpaired_rank_bounds([a, b], rank=3)
    for permutation in itertools.permutations(b):
        actual = sorted(np.array(a) + permutation, reverse=True)[2]
        assert lower <= actual <= upper
    assert upper < max(a) + max(b)
    assert unpaired_rank_bounds([a, b[:-1]])[1] == np.inf
    with pytest.raises(ValueError, match="nonnegative"):
        unpaired_rank_bounds([[1, 2, -1]])


def test_unknown_registry_members_are_visible_instead_of_silently_dropped():
    groups = equivalence_groups(["ENSG00000183628"])
    # Ordinary singletons are still full entities, never inferred missing groups.
    assert len(groups) == 1
    assert json.loads(groups.member_gene_ids.iloc[0])


def test_unnamed_locus_uses_its_stable_id_as_display_label():
    group = equivalence_groups(["ENSG00000198211"])
    assert group.protein_group_label.iloc[0] == "ENSG00000198211"

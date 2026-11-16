import random

import pandas as pd
import pytest

from oncoref.peptide_subtraction import MatchPolicy, PeptideIndex, edit_distance, kmers


def test_boundaries_occurrences_equivalence_and_subtraction():
    rows = pd.DataFrame(
        [
            {"gene_id": "unsafe1", "sequence": "ACDEFGHIKXACDEFGHLK"},
            {"gene_id": "unsafe2", "sequence": "ACDEFGHIKXACDEFGHLK"},
        ]
    )
    index = PeptideIndex(rows, lengths=[9], policy=MatchPolicy(equivalence_groups=("IL",)))
    hits = index.matches("ACDEFGHIK")
    assert len(hits) == 4
    assert set(hits.start) == {0, 10}
    assert set(hits.gene_id) == {"unsafe1", "unsafe2"}
    assert index.subtract("ACDEFGHLK", lengths=[9]).empty
    assert not PeptideIndex(rows, lengths=[9]).subtract("ACDEFGHVK", lengths=[9]).empty
    assert not list(kmers("ACDEXFGHIK", [9]))
    assert index.matches("HIKACDEFG").empty


def test_indels_and_directed_weighted_substitutions():
    policy = MatchPolicy(
        max_distance=0.5, substitution_costs={("I", "V"): 0.25}, insertion_cost=0.5
    )
    assert edit_distance("ACDI", "ACDV", policy=policy) == 0.25
    assert edit_distance("ACDV", "ACDI", policy=policy) == 1
    index = PeptideIndex(pd.DataFrame({"sequence": ["ACDVEFGHIK"]}), policy=policy)
    assert (index.matches("ACDIEFGHIK").distance == 0.25).any()
    assert (index.matches("ACDEFGHIK").distance == 0.5).any()


def test_no_matches_preserve_occurrence_schema_and_provenance():
    records = pd.DataFrame({"gene_id": ["unsafe1"], "sequence": ["ACDEFGHIK"]})
    records.attrs["reference"] = "test-reference"
    for background in (records, records.iloc[:0]):
        index = PeptideIndex(background)
        hits = index.matches("WWWWWWWWW")
        assert hits.empty
        assert list(hits.columns) == [
            "gene_id",
            "query",
            "matched_peptide",
            "start",
            "end",
            "distance",
        ]
        assert hits.attrs["reference"] == "test-reference"


@pytest.mark.parametrize(
    "policy",
    [
        MatchPolicy(max_distance=1),
        MatchPolicy(
            max_distance=0.6,
            substitution_costs={("A", "C"): 0.3},
            insertion_cost=0.6,
            deletion_cost=0.4,
            equivalence_groups=("IL",),
        ),
        MatchPolicy(max_distance=4),
    ],
)
def test_seed_search_equals_exhaustive_windows(policy):
    rng = random.Random(421)
    sequences = ["".join(rng.choices("ACDIL", k=14)) for _ in range(12)]
    records = pd.DataFrame({"sequence": sequences, "record_id": range(len(sequences))})
    index = PeptideIndex(records, policy=policy, lengths=range(3, 7))
    for query in [s[3:8] for s in sequences] + ["AAAA"]:
        expected = {
            (r, start, end)
            for r, s in enumerate(sequences)
            for kmer, start, end in kmers(s, range(3, 7))
            if edit_distance(query, kmer, policy=policy) <= policy.max_distance + 1e-12
        }
        actual = index.matches(query)
        observed = set(zip(actual.record_id, actual.start, actual.end)) if len(actual) else set()
        assert observed == expected


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_distance": -1},
        {"max_distance": float("nan")},
        {"equivalence_groups": ("IL", "LV")},
        {"substitution_costs": {("A", "C"): 0}},
        {"insertion_cost": 0},
        {"equivalence_groups": ("IX",)},
    ],
)
def test_invalid_policies(kwargs):
    with pytest.raises(ValueError):
        MatchPolicy(**kwargs)

import pandas as pd
import pytest

from oncoref.cta_funnel import classify_layers


def test_later_evidence_can_classify_a_gene_held_at_an_earlier_layer():
    genes = pd.DataFrame({"gene_id": ["missing", "later_unsafe", "candidate", "overlap"]})
    layers = [
        ("coverage", pd.Series([True, True, False, True]), pd.Series([False] * 4)),
        ("RNA", pd.Series([False, True, False, True]), pd.Series([False, True, False, True])),
        ("protein", pd.Series([False, False, False, True]), pd.Series([False, False, False, True])),
    ]
    result, funnel, independent = classify_layers(genes, layers)
    assert result.category.tolist() == ["Neither", "Unsafe", "CTA", "Unsafe"]
    assert funnel.iloc[1].new_Unsafe == 2
    assert funnel.iloc[1].new_candidate_removals == 0
    assert independent.iloc[2].Unsafe_evidence_anywhere == 1
    assert independent.iloc[2].Unsafe_reclassified_if_layer_removed == 0
    reverse, _, _ = classify_layers(genes, list(reversed(layers)))
    assert reverse.category.equals(result.category)
    assert (funnel.CTA_eligible_so_far + funnel.Unsafe + funnel.Neither).eq(4).all()


def test_no_missing_or_misaligned_decisions():
    genes = pd.DataFrame({"gene_id": ["a"]})
    with pytest.raises(ValueError, match="boolean"):
        classify_layers(genes, [("RNA", pd.Series([None]), pd.Series([False]))])
    with pytest.raises(ValueError, match="boolean"):
        classify_layers(genes, [("RNA", pd.Series([False], index=[2]), pd.Series([False]))])

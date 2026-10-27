"""Anatomical aggregation must follow the source's consensus convention."""

import pandas as pd
from cta_protein_group_normal import consensus_bulk_regions


def test_brain_subregions_are_mapped_after_locus_sums_and_never_added_together():
    components = pd.DataFrame(
        {
            "gene_id": ["group-representative"] * 4,
            "Tissue": ["caudate", "putamen", "nucleus accumbens", "hippocampus"],
            "cohort": ["gtex"] * 4,
            "nTPM": [1.0, 2.0, 3.0, 4.0],
        }
    )
    result = consensus_bulk_regions(components).set_index("Tissue")
    assert set(result.index) == {"basal ganglia", "hippocampal formation"}
    assert result.loc["basal ganglia", "nTPM"] == 3
    assert result.loc["basal ganglia", "source_tissue"] == "nucleus accumbens"

"""Cross-check the published study and rank tables through the catalog API."""

from pathlib import Path

import pandas as pd

from oncoref import antigen_data, target_selection

DOCS = Path(__file__).resolve().parents[1] / "docs"
ASSETS = DOCS / "antigen-whitepaper"


def test_tables_agree_with_unfiltered_membership_and_preserve_unmeasured_peaks():
    catalog = target_selection.results().set_index("entity_id")
    ranks = pd.read_csv(ASSETS / "tumor-ranks.csv").set_index("entity_id")
    assert ranks.index.is_unique and set(ranks.index) == set(
        catalog.index[catalog.is_target_candidate]
    )
    for name, member in [
        ("ca", "is_cancer_antigen"),
        ("ra", "is_reproductive_antigen"),
        ("cta", "is_cta"),
    ]:
        expected = (
            catalog[catalog[member] & catalog.tumor_max_tpm.notna()]
            .reset_index()
            .sort_values(
                ["tumor_max_tpm", "protein_group_label", "entity_id"], ascending=[False, True, True]
            )
        )
        actual = ranks[f"rank_{name}"].dropna().sort_values()
        assert actual.index.tolist() == expected.entity_id.tolist()
        assert actual.tolist() == list(range(1, len(actual) + 1))
    studied = pd.read_csv(ASSETS / "studied-proteins.csv")
    mage = studied.set_index("protein_group_label").loc["MAGEA3"]
    assert mage.classification == "Unsafe" and pd.isna(mage.tumor_max_tpm)
    assert "colon" in mage.matched_normal_assay and "nTPM" in mage.matched_normal_assay
    sources = pd.read_csv(ASSETS / "sources.csv", keep_default_na=False)
    pd.testing.assert_frame_equal(sources[antigen_data.sources().columns], antigen_data.sources())

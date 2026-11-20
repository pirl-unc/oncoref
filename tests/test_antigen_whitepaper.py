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


def test_whitepaper_separates_measured_peaks_by_category_and_explains_skipped_comparisons():
    paper = (DOCS / "antigen-whitepaper.md").read_text()
    catalog = target_selection.results()
    for title, category in [
        ("CTA", "CTA"),
        ("CA only", "Cancer Antigen only"),
        ("RA only", "Reproductive Antigen only"),
        ("Other", "Other"),
    ]:
        section = paper.split(f"### {title}\n", 1)[1].split("\n### ", 1)[0]
        rows = [
            line.split(" | ")
            for line in section.splitlines()
            if line.startswith("| ") and line[2].isdigit()
        ]
        measured = catalog[catalog.category.eq(category)].dropna(subset=["tumor_max_tpm"])
        expected = measured.sort_values(
            ["tumor_max_tpm", "protein_group_label", "entity_id"], ascending=[False, True, True]
        ).head(10)
        assert [row[1] for row in rows] == expected.protein_group_label.tolist()
    assert "| MAGEA3 | 8 | Unsafe | NR | NR | NR |" in paper
    assert "human antigen" not in paper.lower()
    for figure in [
        "summary/classification-overview.png",
        "summary/decision-tree.png",
        "summary/normal-evidence-layers.png",
        "summary/funnel-all.png",
        "summary/funnel-published.png",
        "summary/cancer-enrichment.png",
        "summary/reproductive-enrichment.png",
        "summary/host-profiles-matrix.png",
        "summary/host-profiles-outcomes.png",
        "summary/sex-restricted-ra-cta.png",
        "summary/sex-restricted-landscape.png",
        "summary/venn.png",
        "summary/literature-sources-counts.png",
        "summary/literature-sources-retention.png",
        "normal-rna/tpm-coverage-somatic.png",
        "normal-rna/tpm-coverage-neural-repro.png",
    ]:
        path = f"audits/target-selection-20261110/rebuild/{figure}"
        assert f"]({path})" in paper and (DOCS / path).is_file()
    lit_csv = DOCS / "audits/target-selection-20261110/rebuild/summary/literature-sources.csv"
    assert lit_csv.is_file()
    lit = pd.read_csv(lit_csv)
    assert {
        "CTA / cancer-germline",
        "Oncofetal / DUX4 / embryo",
        "Placental / reproductive",
    } <= set(lit["antigen_class"])
    excluded_section = paper.split("### Studied proteins in Other or Unsafe\n", 1)[1].split(
        "\n## ", 1
    )[0]
    excluded_rows = [
        line
        for line in excluded_section.splitlines()
        if line.startswith("| ")
        and not line.startswith("| ---")
        and not line.startswith("| Studied protein")
    ]
    assert len(excluded_rows) == 11
    assert all(not row.endswith("| — |") for row in excluded_rows)

"""Full paper unions, exact covers and citation coverage are one contract."""

import pytest

from oncoref import cta
from oncoref.cta_provenance import (
    candidate_evidence,
    candidate_provenance,
    legacy_only_candidates,
    minimum_covers,
    paper_membership,
    selected_membership,
    source_cover,
)


def test_exact_cover_beats_greedy_and_preserves_ties():
    sets = {"a": {1, 2, 3, 4}, "b": {1, 2, 5}, "c": {3, 4, 6}}
    assert minimum_covers(range(1, 7), sets) == [("b", "c")]
    assert minimum_covers({1}, {"b": {1}, "a": {1}}) == [("a",), ("b",)]
    assert minimum_covers(set(), {}) == [()]
    with pytest.raises(ValueError, match="lack paper provenance"):
        minimum_covers({1, 2}, {"a": {1}})


def test_minimum_counts_papers_once_and_has_indispensable_witnesses():
    cover = source_cover()
    assert len(paper_membership().source_tag.unique()) > len(cover["papers_considered"])
    assert cover["minimum_papers"] == 10
    assert len(cover["all_minimum_covers"]) == 1
    assert all(cover["indispensable_witness_gene_ids"].values())
    assert "10.1038/s41467-020-19141-w" not in cover["selected_dois"]


def test_complete_lists_are_not_truncated_to_retained_genes():
    refs = selected_membership()
    candidates = candidate_evidence()
    ids = set(candidates.Ensembl_Gene_ID)
    assert ids == set(refs.loc[refs.protein_candidate_eligible, "Ensembl_Gene_ID"])
    assert len(ids) == 2537
    assert len(cta.cta_gene_ids()) == 624
    assert cta.cta_gene_ids() < ids
    coverage = candidate_provenance()
    assert coverage.Ensembl_Gene_ID.is_unique
    assert set(coverage.Ensembl_Gene_ID) == ids
    assert coverage.paper_dois.ne("").all() and coverage.source_locations.ne("").all()
    assert set(coverage.loc[coverage.default_panel, "Ensembl_Gene_ID"]) == cta.cta_gene_ids()
    legacy = legacy_only_candidates()
    assert len(legacy) == 9 and "TRIM64" in set(legacy.Symbol)
    assert set(legacy.Ensembl_Gene_ID).isdisjoint(cta.cta_gene_ids())


def test_normal_only_multi_source_nominations_cannot_bypass_holdback():
    from oncoref.load_dataset import get_data

    raw = get_data("cancer-testis-antigens")
    normal = {
        "daSilva2017_testis_biased",
        "Gong2021_reproductive_PC",
        "Gong2021_reproductive_ncRNA",
    }
    mask = (
        raw.source_databases.fillna("")
        .str.split(";")
        .map(lambda tags: bool(set(tags) - {""}) and set(tags) - {""} <= normal)
    )
    assert mask.sum() > 100
    assert cta.passes_filters_mask(raw[mask]).any()
    assert set(raw.loc[mask, "Ensembl_Gene_ID"]).isdisjoint(cta.cta_gene_ids())


def test_targeted_citations_preserve_gene_resolution_and_normal_tissue_caveats():
    from oncoref.cta_sources import gene_publication_evidence

    evidence = gene_publication_evidence()
    sun5 = evidence[evidence.Symbol.eq("SUN5")].iloc[0]
    assert sun5.source_anchor == "Figure 1C-E" and "Adjacent" in sun5.finding
    ssx4b = evidence[evidence.Symbol.eq("SSX4B")]
    assert set(ssx4b.source_tag) == {"Caballero2014_targeted"}
    assert "unresolved" in ssx4b.iloc[0].limitations
    assert not evidence.hla_peptide_validation.any()
    assert not evidence.antigen_specific_t_cell_validation.any()


def test_source_expansion_preserves_the_entire_previously_retained_panel():
    from pathlib import Path

    import pandas as pd

    old = pd.read_csv(Path(__file__).parent / "fixtures/cta-default-1.8.205.csv")
    assert len(old) == 297
    assert set(old.gene_id) <= cta.cta_gene_ids()
    assert "TRIM64" not in cta.cta_gene_names()


def test_historical_protein_tags_are_not_confused_with_paper_membership():
    from oncoref.cta_provenance import historical_tag_audit

    audit = historical_tag_audit()
    garin4 = audit[audit.Symbol.eq("GARIN4") & audit.historical_tag.eq("daSilva2017_protein")]
    assert garin4.status.tolist() == ["not_reproduced_by_complete_import"]
    refs = paper_membership()
    s3 = refs[refs.source_tag.eq("daSilva2017_tumor_proteomics")]
    assert {"GARIN5B", "MEIOC"} <= set(s3.Symbol)
    assert "GARIN4" not in set(s3.Symbol)


def test_per_gene_evidence_modalities_do_not_transfer_antigen_validation():
    from oncoref.cta_provenance import gene_citation_evidence_report

    report = gene_citation_evidence_report().set_index("Symbol")
    assert len(report) == 2537
    assert report.published_list_nomination_sources.ne("").all()
    assert set(report.index[report.hla_peptide_evidence_sources.ne("")]) == {"VGLL1"}
    assert set(report.index[report.antigen_specific_t_cell_evidence_sources.ne("")]) == {"VGLL1"}
    assert "10.3390/cancers14215368" in report.loc["SUN5", "protein_evidence_sources"]
    assert report.loc["CT47A8", "legacy_retention_evidence_caveat"]
    assert report.loc["CT47A8", "nomination_scope"] == "normal_reproductive_expression_only"


def test_paper_identity_counts_reconcile_with_the_pooled_membership():
    from oncoref.cta_curation_plots import stage_membership
    from oncoref.cta_provenance import paper_intake_counts

    pooled = stage_membership()
    counts = paper_intake_counts()
    tags = pooled.source_databases.str.split(";").map(set)
    for doi, group in selected_membership().groupby("doi"):
        source_tags = set(group.source_tag)
        expected = tags.map(
            lambda values, source_tags=source_tags: bool(values & source_tags)
        ).sum()
        actual = counts[counts.source_tag.eq(doi) & counts.stage.eq("published")].remaining.item()
        assert actual == expected, doi

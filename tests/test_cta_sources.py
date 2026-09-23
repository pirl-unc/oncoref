"""Published nominations retain source identity and do not bypass CTA gates."""

import pandas as pd

from oncoref import cta
from oncoref.cta_sources import (
    BRADLEY,
    BRADLEY_GENES,
    GONG_NC,
    GONG_PC,
    add_publication_candidates,
    intake_counts,
    intake_membership,
    publication_membership,
    publication_sources,
    resolve_membership,
)
from oncoref.load_dataset import get_data


def test_complete_published_lists_and_citations():
    refs = publication_membership()
    assert refs.groupby("source_tag").size().to_dict() == {GONG_PC: 71, GONG_NC: 74, BRADLEY: 10}
    assert not refs.duplicated(["source_tag", "source_symbol"]).any()
    assert set(refs.loc[refs.source_tag.eq(BRADLEY), "source_symbol"]) == set(BRADLEY_GENES)
    assert set(publication_sources().source_tag) == set(refs.source_tag)
    assert publication_sources().input_sha256.str.fullmatch("[a-f0-9]{64}").all()


def test_import_rejects_an_edited_workbook_before_parsing(tmp_path):
    import pytest

    from oncoref.cta_sources import extract_publications

    path = tmp_path / "edited.xlsx"
    path.write_bytes(b"not the published workbook")
    with pytest.raises(ValueError, match="checksum"):
        extract_publications(path)


def test_historical_annotation_is_preserved_without_promoting_noncoding():
    refs = publication_membership()
    erv = refs[refs.source_symbol.eq("ERVH48-1")].iloc[0]
    assert erv.source_tag == GONG_NC and erv.source_biotype == "noncoding"
    assert erv.biotype == "protein_coding" and erv.protein_candidate_eligible
    dscr = refs[refs.source_symbol.eq("DSCR4")].iloc[0]
    assert dscr.source_tag == GONG_PC and dscr.source_biotype == "protein_coding"
    assert dscr.biotype == "lncRNA" and not dscr.protein_candidate_eligible
    unmapped = refs[refs.mapping_status.eq("unmapped")]
    assert len(unmapped) == 10
    assert unmapped.Ensembl_Gene_ID.isna().all()
    assert not unmapped.protein_candidate_eligible.any()


def test_source_id_mapping_does_not_fall_back_to_potentially_wrong_symbol():
    rows = pd.DataFrame(
        [{"source_tag": GONG_PC, "source_gene_id": "ENSG99999999999", "source_symbol": "VGLL1"}]
    )
    resolved = resolve_membership(rows)
    assert resolved.iloc[0].mapping_status == "unmapped"
    assert not resolved.iloc[0].protein_candidate_eligible


def test_import_is_idempotent_and_preserves_existing_evidence():
    table = get_data("cancer-testis-antigens")
    updated = add_publication_candidates(table, publication_membership())
    pd.testing.assert_frame_equal(table, updated)


def test_only_vgll1_is_claimed_as_newly_validated_by_bradley():
    refs = publication_membership()
    validated = refs[refs.validation_scope.str.contains("HLA peptide", na=False)]
    assert validated.source_symbol.tolist() == ["VGLL1"]
    assert validated.source_tag.tolist() == [BRADLEY]


def test_source_funnels_are_nested_and_reconcile_to_public_default():
    refs = intake_membership()
    counts = intake_counts()
    for _, group in counts.groupby("source_tag", sort=False):
        assert (group.remaining.diff().dropna() <= 0).all()
        assert group.dropped.sum() + group.iloc[-1].remaining == group.iloc[0].remaining
    coding_ids = set(refs.loc[refs.protein_coding, "Ensembl_Gene_ID"])
    assert set(refs.loc[refs.default_panel, "Ensembl_Gene_ID"]) == coding_ids & cta.cta_gene_ids()
    assert not (refs.default_panel & ~refs.hpa_restriction).any()
    assert not (refs.hpa_restriction & ~refs.family_eligible).any()
    assert not (refs.family_eligible & ~refs.protein_coding).any()


def test_all_prior_placental_nominations_have_publication_provenance():
    from oncoref.cta_sources import placental_source_coverage

    coverage = placental_source_coverage()
    assert len(coverage) == 19 and coverage.covered_by_publication.all()
    assert coverage.default_panel.sum() == 9
    assert set(
        coverage.loc[
            coverage.Symbol.isin(["CGB1", "CGB2", "CGB7"]) & coverage.default_panel, "Symbol"
        ]
    ) == {"CGB2"}


def test_cgb_assay_resolution_does_not_imply_protein_validation():
    from oncoref.cta_sources import add_gene_evidence_tags, gene_publication_evidence

    rows = gene_publication_evidence()
    assert len(rows) == 10
    direct = rows[(rows.Symbol == "CGB2") & (rows.source_tag == "Rull2005_CGB_placenta")].iloc[0]
    assert direct.assay_resolution == "gene_resolved_restriction_digest"
    assert direct.assay_gene_scope == "CGB2"
    combined = rows[rows.source_tag == "Kubiczak2013_CGB_ovarian"]
    assert set(combined.assay_gene_scope) == {"CGB1;CGB2"}
    assert set(combined.assay_resolution) == {"combined_paralog_assay"}
    for column in [
        "gene_specific_protein_validation",
        "hla_peptide_validation",
        "antigen_specific_t_cell_validation",
    ]:
        assert not rows[column].any()
    assert set(rows.loc[rows.source_tag.str.startswith("McKellar"), "publication_status"]) == {
        "preprint"
    }
    table = get_data("cancer-testis-antigens")
    pd.testing.assert_frame_equal(add_gene_evidence_tags(table), table)

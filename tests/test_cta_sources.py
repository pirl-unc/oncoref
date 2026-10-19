"""Published nominations retain source identity and do not bypass CTA gates."""

import io

import pandas as pd
import pytest

from oncoref import cta
from oncoref.cta_sources import (
    BRADLEY,
    BRADLEY_GENES,
    GONG_NC,
    GONG_PC,
    add_gene_evidence_tags,
    add_publication_candidates,
    extract_publications,
    intake_counts,
    intake_membership,
    publication_membership,
    publication_sources,
    resolve_membership,
)
from oncoref.load_dataset import get_data


def test_shipped_candidate_import_is_idempotent():
    original = get_data("cancer-testis-antigens")
    imported = add_gene_evidence_tags(
        add_publication_candidates(original, publication_membership())
    )
    pd.testing.assert_frame_equal(original, imported)


def test_complete_published_lists_and_citations():
    refs = publication_membership()
    counts = refs.groupby("source_tag").size().to_dict()
    assert {tag: counts[tag] for tag in (GONG_PC, GONG_NC, BRADLEY)} == {
        GONG_PC: 71,
        GONG_NC: 74,
        BRADLEY: 10,
    }
    assert not refs.duplicated(["source_tag", "source_symbol"]).any()
    assert set(refs.loc[refs.source_tag.eq(BRADLEY), "source_symbol"]) == set(BRADLEY_GENES)
    assert set(publication_sources().source_tag) == set(refs.source_tag)
    assert publication_sources().input_sha256.str.fullmatch("[a-f0-9]{64}").all()


def test_historical_annotation_is_preserved_without_promoting_noncoding():
    refs = publication_membership()
    erv = refs[refs.source_symbol.eq("ERVH48-1")].iloc[0]
    assert erv.source_tag == GONG_NC and erv.source_biotype == "noncoding"
    assert erv.biotype == "protein_coding" and erv.protein_candidate_eligible
    dscr = refs[refs.source_symbol.eq("DSCR4")].iloc[0]
    assert dscr.source_tag == GONG_PC and dscr.source_biotype == "protein_coding"
    assert dscr.biotype == "lncRNA" and not dscr.protein_candidate_eligible
    unmapped = refs[refs.source_tag.isin([GONG_PC, GONG_NC]) & refs.mapping_status.eq("unmapped")]
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


@pytest.mark.parametrize("add_new_rows", [False, True])
def test_legacy_core_import_is_complete_idempotent_and_auditable(monkeypatch, add_new_rows):
    from oncoref import cta_regen

    shipped = get_data("cancer-testis-antigens")
    legacy = shipped.drop(columns="rna_reproductive_tissue_scope").copy()
    if add_new_rows:
        # The original 397 candidate rows remain at the start of the expanded
        # snapshot. Exercise importing the full current intake into that legacy
        # schema, without requiring repository history in installed tests.
        legacy = legacy.iloc[:397].copy()
    original = legacy.copy(deep=True)
    new_ids = set(shipped.Ensembl_Gene_ID) - set(legacy.Ensembl_Gene_ID)
    regenerated_ids = []

    def regenerate_new_rows(seed):
        regenerated_ids.extend(seed.Ensembl_Gene_ID)
        # HPA regeneration itself has real-data parity tests. This regression
        # isolates the old/new schema merge and must never regenerate old rows.
        assert set(seed.Ensembl_Gene_ID) == new_ids
        return (
            shipped.set_index("Ensembl_Gene_ID", drop=False)
            .loc[seed.Ensembl_Gene_ID]
            .reset_index(drop=True)
        )

    monkeypatch.setattr(cta_regen, "regenerate_cta_columns", regenerate_new_rows)
    membership = publication_membership()
    migrated = add_publication_candidates(legacy, membership)
    assert len(migrated) == len(shipped)
    assert list(migrated.columns) == list(shipped.columns)
    assert migrated.rna_reproductive_tissue_scope.eq("core").all()
    assert not migrated.rna_reproductive_tissue_scope.isna().any()
    pd.testing.assert_frame_equal(legacy, original)
    pd.testing.assert_frame_equal(migrated, add_publication_candidates(migrated, membership))
    # CSV output is the downstream contract, including re-import after reload.
    reloaded = pd.read_csv(io.StringIO(migrated.to_csv(index=False)))
    pd.testing.assert_frame_equal(
        reloaded.sort_values("Ensembl_Gene_ID").reset_index(drop=True),
        shipped.sort_values("Ensembl_Gene_ID").reset_index(drop=True),
    )
    pd.testing.assert_frame_equal(reloaded, add_publication_candidates(reloaded, membership))
    pd.testing.assert_frame_equal(intake_membership(table=reloaded), intake_membership())
    pd.testing.assert_frame_equal(intake_counts(table=reloaded), intake_counts())
    assert len(regenerated_ids) == len(new_ids) == (len(shipped) - 397 if add_new_rows else 0)


@pytest.mark.parametrize("scope", [pd.NA, "", "extended", "unknown"])
def test_import_and_intake_reject_incomplete_or_mixed_core_scope(scope):
    table = get_data("cancer-testis-antigens").iloc[:2].copy()
    table["rna_reproductive_tissue_scope"] = pd.Series(["core", scope], dtype="string")
    with pytest.raises(ValueError, match="core"):
        add_publication_candidates(table, publication_membership())
    with pytest.raises(ValueError, match="tissue_scope"):
        intake_membership(table=table)


def test_legacy_core_evidence_cannot_be_audited_as_extended():
    legacy = get_data("cancer-testis-antigens").drop(columns="rna_reproductive_tissue_scope")
    pd.testing.assert_frame_equal(intake_membership(table=legacy), intake_membership())
    with pytest.raises(ValueError, match="tissue_scope"):
        intake_membership(table=legacy, tissue_scope="extended")


def test_only_vgll1_is_claimed_as_newly_validated_by_bradley():
    refs = publication_membership()
    validated = refs[refs.validation_scope.str.contains("HLA peptide", na=False)]
    assert validated.source_symbol.tolist() == ["VGLL1"]
    assert validated.source_tag.tolist() == [BRADLEY]


@pytest.mark.parametrize("tissue_scope", ["core", "extended"])
def test_source_funnels_are_nested_and_reconcile_to_public_default(tissue_scope):
    refs = intake_membership(tissue_scope=tissue_scope)
    counts = intake_counts(tissue_scope=tissue_scope)
    for _, group in counts.groupby("source_tag", sort=False):
        assert (group.remaining.diff().dropna() <= 0).all()
        assert group.dropped.sum() + group.iloc[-1].remaining == group.iloc[0].remaining
    coding_ids = set(refs.loc[refs.protein_coding, "Ensembl_Gene_ID"])
    default = cta.cta_gene_ids() if tissue_scope == "core" else cta.cta_extended_gene_ids()
    assert set(refs.loc[refs.default_panel, "Ensembl_Gene_ID"]) == coding_ids & default
    assert not (refs.default_panel & ~refs.hpa_restriction).any()
    assert not (refs.hpa_restriction & ~refs.family_eligible).any()
    assert not (refs.family_eligible & ~refs.protein_coding).any()


def test_distinct_source_symbols_cannot_collapse_within_one_source():
    rows = pd.DataFrame(
        [
            {"source_tag": GONG_PC, "source_gene_id": "ENSG00000101951", "source_symbol": "PAGE4"},
            {"source_tag": GONG_PC, "source_gene_id": "ENSG00000101951.2", "source_symbol": "CT16"},
        ]
    )
    with pytest.raises(ValueError, match=r"PAGE4.*ENSG00000101951"):
        resolve_membership(rows)
    rows.loc[1, "source_tag"] = BRADLEY
    assert len(resolve_membership(rows)) == 2  # Cross-source agreement is legitimate.


def test_versioned_candidate_is_annotated_without_duplicate_addition():
    table = get_data("cancer-testis-antigens")
    table["Ensembl_Gene_ID"] += ".1"
    updated = add_publication_candidates(table, publication_membership())
    pd.testing.assert_frame_equal(table, updated)
    duplicate = pd.concat(
        [table, table.iloc[[0]].assign(Ensembl_Gene_ID=table.iloc[0].Ensembl_Gene_ID.split(".")[0])]
    )
    with pytest.raises(ValueError, match="Duplicate candidate Ensembl IDs"):
        add_publication_candidates(duplicate, publication_membership())


def test_intake_uses_supplied_table_and_normalizes_versions():
    table = get_data("cancer-testis-antigens")
    table = table[table.Symbol.eq("PAGE4")].copy()
    table["Ensembl_Gene_ID"] += ".2"
    table["passes_filters"] = True
    refs = intake_membership(table=table).set_index("Symbol")
    assert refs.loc[["PAGE4"], "family_eligible"].all()
    assert refs.loc[["PAGE4"], "hpa_restriction"].all()
    assert not refs.loc[["PAGE4"], "default_panel"].any()  # Reviewed exclusion still applies.
    assert not refs.loc[["INSL4"], "family_eligible"].any()  # Absent from the supplied table.


def test_new_candidates_without_hpa_fail_closed(monkeypatch):
    from oncoref import hpa

    monkeypatch.setattr(
        hpa, "hpa_rna_consensus", lambda: pd.DataFrame(columns=["Gene", "Tissue", "nTPM"])
    )
    monkeypatch.setattr(
        hpa,
        "hpa_normal_tissue",
        lambda: pd.DataFrame(columns=["Gene", "Tissue", "Level", "Reliability"]),
    )
    refs = publication_membership()
    refs = refs[refs.source_symbol.isin(["PAGE4", "DSCR4", "TXNRD3NB"])]
    with pytest.warns(UserWarning, match="absent from HPA"):
        result = add_publication_candidates(get_data("cancer-testis-antigens").iloc[:0], refs)
    assert result.Symbol.tolist() == ["PAGE4"]
    row = result.iloc[0]
    assert pd.isna(row.rna_deflated_reproductive_frac)
    assert pd.isna(row.rna_max_ntpm)
    assert not row.passes_filters
    assert row.rna_restriction_level == "NO_DATA"


@pytest.mark.parametrize(
    "fault, message",
    [
        ("row_count", "70 placenta rows, expected 71"),
        ("duplicate_id", "only 70 distinct source Ensembl IDs"),
        ("missing_id", "Missing source Ensembl ID"),
        (None, None),
    ],
)
def test_workbook_extraction_guards_source_shape(tmp_path, monkeypatch, fault, message):
    import openpyxl

    workbook = openpyxl.Workbook()
    workbook.remove(workbook.active)
    refs = publication_membership()
    for name, tag in [
        ("Data 5 - tissue-enriched PC", GONG_PC),
        ("Data 6 tissue-enriched lncR", GONG_NC),
    ]:
        sheet = workbook.create_sheet(name)
        for row in refs.loc[refs.source_tag.eq(tag)].itertuples():
            sheet.append(["placenta", row.source_gene_id, row.source_symbol])
    sheet = workbook.worksheets[0]
    if fault == "row_count":
        sheet.delete_rows(1)
    elif fault == "duplicate_id":
        sheet["B2"] = sheet["B1"].value
    elif fault == "missing_id":
        sheet["B1"] = None
    path = tmp_path / "source.xlsx"
    workbook.save(path)
    workbook.close()
    # Isolate shape validation; a separate regression checks the real input pin.
    import hashlib

    from oncoref import cta_sources

    sources = publication_sources()
    sources.loc[sources.source_tag.isin([GONG_PC, GONG_NC]), "input_sha256"] = hashlib.sha256(
        path.read_bytes()
    ).hexdigest()
    monkeypatch.setattr(cta_sources, "publication_sources", lambda: sources)
    if message:
        with pytest.raises(ValueError, match=message):
            extract_publications(path)
    else:
        extracted = extract_publications(path)
        pd.testing.assert_frame_equal(
            extracted.fillna(""),
            refs.loc[refs.source_tag.isin([GONG_PC, GONG_NC, BRADLEY]), extracted.columns]
            .reset_index(drop=True)
            .fillna(""),
            check_dtype=False,
        )


@pytest.mark.parametrize("apply", [False, True])
def test_import_cli_requires_explicit_date_and_defaults_to_sidecars(tmp_path, monkeypatch, apply):
    import sys

    import import_cta_publications as importer

    tags = {GONG_PC, GONG_NC, BRADLEY}
    workbook = tmp_path / "source.xlsx"
    workbook.write_bytes(b"source bytes used for hashing")
    original_membership = publication_membership()
    monkeypatch.setattr(
        importer,
        "extract_publications",
        lambda _: original_membership[original_membership.source_tag.isin(tags)],
    )
    monkeypatch.setattr(
        importer,
        "regenerate_cta_columns",
        lambda table, **kw: get_data("cancer-testis-antigens-extended"),
    )
    names = [
        "cancer-testis-antigens.csv",
        "cancer-testis-antigens-extended.csv",
        "cta-publication-sources.csv",
        "cta-publication-membership.csv",
    ]
    for name in names:
        (tmp_path / name).write_text("original\n")
    args = ["import", "--gong-xlsx", str(workbook), "--out-dir", str(tmp_path)]
    monkeypatch.setattr(sys, "argv", args)
    with pytest.raises(SystemExit) as exc:
        importer.main()
    assert exc.value.code == 2
    args += ["--retrieved-date", "2026-09-23"] + (["--apply"] if apply else [])
    monkeypatch.setattr(sys, "argv", args)
    importer.main()
    for name in names:
        dest = tmp_path / (name if apply else name.replace(".csv", ".import.csv"))
        assert dest.exists()
        if not apply:
            assert (tmp_path / name).read_text() == "original\n"
    sources = pd.read_csv(
        tmp_path
        / ("cta-publication-sources.csv" if apply else "cta-publication-sources.import.csv")
    )
    assert set(sources.loc[sources.source_tag.isin(tags), "retrieved_date"]) == {"2026-09-23"}
    memberships = pd.read_csv(
        tmp_path
        / ("cta-publication-membership.csv" if apply else "cta-publication-membership.import.csv")
    )
    pd.testing.assert_frame_equal(memberships, original_membership)
    pd.testing.assert_frame_equal(
        sources.loc[~sources.source_tag.isin(tags)].reset_index(drop=True),
        publication_sources()
        .loc[~publication_sources().source_tag.isin(tags)]
        .reset_index(drop=True),
    )


def test_import_rejects_an_edited_workbook_before_parsing(tmp_path):
    import pytest

    from oncoref.cta_sources import extract_publications

    path = tmp_path / "edited.xlsx"
    path.write_bytes(b"not the published workbook")
    with pytest.raises(ValueError, match="checksum"):
        extract_publications(path)


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
    assert len(rows[rows.Symbol.isin(["CGB1", "CGB2", "CGB7"])]) == 10
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
    annotated = add_gene_evidence_tags(table)
    pd.testing.assert_frame_equal(
        annotated.drop(columns="source_databases"), table.drop(columns="source_databases")
    )
    assert "Song2022_targeted" in annotated.set_index("Symbol").loc["SUN5", "source_databases"]

"""Complete literature intake precedes curation and preserves assay limitations."""

import pandas as pd
import pytest

from oncoref import cta, cta_curation_plots
from oncoref.cta_landscape import INPUTS, extract_landscapes, resolve_landscape_membership
from oncoref.cta_sources import publication_membership, publication_sources
from oncoref.load_dataset import get_data


def test_full_paper_sets_are_imported_not_just_the_twelve_audit_genes():
    refs = publication_membership()
    counts = refs.groupby("source_tag").size()
    expected = {
        "Wang2016_CT": 1019,
        "Bruggeman2018_GC": 756,
        "daSilva2017_testis_biased": 1103,
        "daSilva2017_CT": 745,
        "Jamin2021_CT": 607,
        "Jamin2021_core": 125,
        "Chang2019_TGCT": 1036,
        "Carter2023_CT": 103,
        "Seager2024_CTA": 17,
        "Loriot2025_S1": 146,
        "Loriot2025_S2": 134,
        "Gong2021_reproductive_PC": 744,
        "Gong2021_reproductive_ncRNA": 1067,
        "daSilva2017_tumor_proteomics": 136,
        "Bai2016_EGFL6": 1,
    }
    assert counts.loc[list(expected)].to_dict() == expected
    source_counts = publication_sources().set_index("source_tag").published_rows
    assert (source_counts.loc[counts.index] == counts).all()
    loci = set(refs.loc[refs.protein_candidate_eligible, "Ensembl_Gene_ID"])
    assert loci <= set(get_data("cancer-testis-antigens").Ensembl_Gene_ID)
    assert len(loci) > 2000


def test_exact_locations_for_previously_unsourced_genes_and_scope():
    refs = publication_membership()
    required = {
        "CT47A8",
        "CT47A9",
        "CT47A10",
        "CT45A8",
        "GAGE12D",
        "GAGE10",
        "MAGEA2B",
        "SSX4B",
        "SSX7",
        "NLRP9",
        "ZFP42",
        "SUN5",
    }
    assert required <= set(refs.Symbol)
    wang = refs[refs.source_tag.eq("Wang2016_CT")].set_index("source_symbol")
    assert str(wang.loc["SSX7", "source_rows"]) == "814"
    assert wang.loc["SSX7", "source_gene_id"] == "ENSG00000187754"
    brug = refs[refs.source_tag.eq("Bruggeman2018_GC")].set_index("source_symbol")
    assert set(brug.loc[["NLRP9", "ZFP42"], "source_subset"]) == {"TGCT_dependent"}
    carter = set(refs.loc[refs.source_tag.eq("Carter2023_CT"), "Symbol"])
    assert {"SUN5", "ZFP42"}.isdisjoint(carter)
    jamin = refs[refs.source_tag.eq("Jamin2021_CT")].set_index("source_symbol")
    assert "SSX4B;SSX4" in jamin.loc["SSX4B", "assay_resolution"]
    assert "shared probe" in jamin.loc["SSX4B", "assay_resolution"]


def test_unknown_probes_and_historical_ids_are_not_guessed():
    rows = pd.DataFrame(
        [
            {
                "source_gene_id": "ENSG99999999999",
                "source_symbol": "MAGEA4",
                "source_id_type": "ensembl",
            },
            {
                "source_gene_id": "unknown_probe",
                "source_symbol": "MAGEA4",
                "source_id_type": "probe",
            },
            {"source_gene_id": "", "source_symbol": "MAGEA4", "source_id_type": "symbol"},
        ]
    )
    resolved = resolve_landscape_membership(rows)
    assert list(resolved.mapping_status) == ["unmapped", "unmapped", "mapped"]
    assert not resolved.iloc[:2].protein_candidate_eligible.any()
    assert resolved.iloc[2].Ensembl_Gene_ID == "ENSG00000147381"
    refs = publication_membership()
    probes = refs[refs.source_id_type.eq("probe")]
    assert len(probes) == 56
    assert probes.Ensembl_Gene_ID.isna().all()


def test_input_bytes_are_verified_before_selection(tmp_path):
    first = next(iter(INPUTS))
    file = tmp_path / first
    file.parent.mkdir(parents=True)
    file.write_bytes(b"changed supplement")
    with pytest.raises(ValueError, match="checksum"):
        extract_landscapes(tmp_path)


def test_testis_only_nominations_cannot_become_default_by_hpa_alone():
    raw = get_data("cancer-testis-antigens")
    broad = raw[raw.source_databases.eq("daSilva2017_testis_biased")]
    assert len(broad) > 100 and cta.passes_filters_mask(broad).any()
    assert set(broad.Ensembl_Gene_ID).isdisjoint(cta.cta_gene_ids())
    evidence = cta.cta_evidence().set_index("Ensembl_Gene_ID")
    present = broad.Ensembl_Gene_ID[broad.Ensembl_Gene_ID.isin(evidence.index)]
    assert set(evidence.loc[present, "specificity_action"]) == {"candidate_only"}
    assert {"CT45A8", "CT47A8"} <= cta.cta_gene_names()


def test_missing_hpa_and_existing_review_decisions_remain_excluded():
    raw = get_data("cancer-testis-antigens")
    missing = raw[raw.rna_max_ntpm.isna()]
    assert len(missing) == 6
    assert not cta.passes_filters_mask(missing).any()
    assert set(missing.Ensembl_Gene_ID).isdisjoint(cta.cta_gene_ids())
    assert {"TRIM64", "CTAG2", "CSAG1", "SPAG4"}.isdisjoint(cta.cta_gene_names())
    assert raw.Symbol.notna().all()


def test_funnel_and_overlap_count_distinct_loci_and_reconcile():
    raw = get_data("cancer-testis-antigens")
    membership = cta_curation_plots.stage_membership()
    assert membership.identity.is_unique
    assert len(membership) > len(raw)
    assert membership.protein_coding.sum() == len(cta_curation_plots._evidence())
    assert set(membership.loc[membership.default_panel, "Ensembl_Gene_ID"]) == cta.cta_gene_ids()
    assert not membership.loc[~membership.mapped, "protein_coding"].any()
    overlaps = cta_curation_plots.source_overlap_counts()
    a = overlaps.pivot(index="source_a", columns="source_b", values="candidate_overlap")
    pd.testing.assert_frame_equal(a, a.T.rename_axis(index="source_a", columns="source_b"))
    sets = cta_curation_plots._tag_sets(raw)
    for source, genes in sets.items():
        assert a.loc[source, source] == len(genes)
    assert a.loc["Wang 2016", "Wang 2016"] > 900

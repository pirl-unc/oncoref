"""Paper-level topics cannot become gene-level evidence without a claim."""

from copy import deepcopy

import pandas as pd
import pytest

from oncoref import cta_literature as literature
from oncoref import target_filters, target_selection


def test_reviewed_gene_uses_override_source_roles_and_old_group_ids(monkeypatch):
    catalog = pd.DataFrame(
        [
            {"entity_id": "A", "protein_group_label": "A", "member_gene_ids": '["a"]'},
            {"entity_id": "B", "protein_group_label": "B", "member_gene_ids": '["b"]'},
        ]
    )
    record = {
        "source_id": "paper",
        "source_roles": ["tumor_expression"],
        "entity_id": "old-A/B",
        "protein": "A/B",
        "gene_ids": ["a", "b"],
        "nomination_gene_ids": [],
        "published_nomination": False,
        "eligible_for_coverage": True,
        "claim_scope": "Old paper-wide tumor scope",
        "membership_record": {},
        "source_record": {},
        "publication_claims": [],
    }
    claims = [
        {"gene_ids": [gene], "scope": scope, "allowed_uses": [use], "eligible_for_coverage": True}
        for gene, scope, use in [
            ("a", "Reproductive RNA only", "reproductive_support"),
            ("b", "Normal-tissue concern", "normal_tissue_review"),
        ]
    ]
    monkeypatch.setattr(literature, "attestations", lambda: pd.DataFrame([record]))
    monkeypatch.setattr(
        literature,
        "publication_registry",
        lambda: {"sources": [{"source_id": "paper", "gene_claims": claims}]},
    )
    monkeypatch.setattr(target_selection, "results", lambda gene=None: catalog)
    mapped = literature.current_attestations().set_index("entity_id")
    assert mapped.loc["A", "allowed_uses"] == ["reproductive_support"]
    assert mapped.loc["B", "allowed_uses"] == ["normal_tissue_review"]
    assert mapped.loc["B", "recorded_claim_scope"] == record["claim_scope"]
    # Avoid a cached production result while testing the actual filter consumer.
    table = target_filters.paper_table.__wrapped__()
    assert table.to_dict("records") == [
        {"entity_id": "A", "family": "reproductive", "source_id": "paper"}
    ]
    # A rebuild can supply new identities without reading the previous catalog.
    monkeypatch.setattr(target_selection, "results", lambda *args: pytest.fail("Stale catalog"))
    assert len(literature.current_attestations(catalog=catalog)) == 2


@pytest.mark.parametrize(
    "gene,old_name,row",
    [
        ("FOXL3", "WI2-2373I1.1", 1274),
        ("SSU72L6", "SSU72P8", 1242),
    ],
)
def test_native_testis_rows_override_placenta_paper_title(gene, old_name, row):
    claims = literature.paper_gene_claims(gene)
    assert len(claims) == 1
    claim = claims.iloc[0]
    assert claim.reported_name == old_name
    assert claim.provenance["source_row"] == row
    assert claim.provenance["source_tissue"] == "Testis"
    assert claim.allowed_uses == ["reproductive_support"]
    assert claim.curation_status == "reviewed_gene_claim"
    paper = literature.papers(gene).iloc[0]
    assert paper.gene_ids == [claim.gene_id]
    assert len(paper.gene_claims) == 1


def test_normal_counterevidence_is_not_ra_or_tumor_support():
    claims = literature.paper_gene_claims("SCGB1D4")
    assert len(claims) == 2
    assert all("normal_tissue_review" in uses for uses in claims.allowed_uses)
    assert not {"reproductive_support", "tumor_support", "cta_nomination"} & {
        use for uses in claims.allowed_uses for use in uses
    }
    assert literature.paper_gene_claims("SCGB1D4", use="reproductive_support").empty


def test_group_coverage_does_not_transfer_eligibility_between_member_loci(monkeypatch):
    catalog = pd.DataFrame(
        [
            {
                "entity_id": "AB",
                "protein_group_label": "A/B",
                "member_gene_ids": '["a", "b"]',
                "category": "CTA",
            }
        ]
    )
    monkeypatch.setattr(target_selection, "results", lambda **kwargs: catalog)
    monkeypatch.setattr(
        literature,
        "current_attestations",
        lambda: pd.DataFrame(
            [
                {
                    "source_id": "paper",
                    "gene_ids": ["a", "b"],
                    "eligible_for_coverage": True,
                    "published_nomination": False,
                    "gene_claims": [
                        {"gene_id": "a", "eligible_for_coverage": True},
                        {"gene_id": "b", "eligible_for_coverage": False},
                    ],
                }
            ]
        ),
    )
    coverage = literature.current_coverage().iloc[0]
    assert coverage.paper_coverage == "some_members"
    assert coverage.covered_gene_ids == ["a"]
    assert coverage.uncovered_gene_ids == ["b"]
    assert coverage.per_gene_sources == {"a": ["paper"], "b": []}


def test_functional_citation_is_not_a_second_perturbation_experiment():
    claims = literature.paper_gene_claims("TPRX2").set_index("source_id")
    primary = claims.loc["10.1126/science.abo7923"]
    assert "combined_gene_perturbation" in primary.evidence_types
    assert "translation_support" in primary.allowed_uses
    later = claims.loc["https://pmc.ncbi.nlm.nih.gov/articles/PMC9978639/"]
    assert later.allowed_uses == ["reproductive_support"]
    assert "not a second independent" in later.limitations
    assert literature.paper_gene_claims("unmapped-test-gene").empty
    assert literature.papers("unmapped-test-gene").empty


def test_curation_schema_rejects_duplicate_genes_and_unsupported_uses():
    data = literature.publication_registry()
    duplicate = deepcopy(data)
    duplicate["sources"][0]["gene_claims"] *= 2
    with pytest.raises(ValueError, match="Consolidate"):
        literature.validate_publications(duplicate)
    invalid = deepcopy(data)
    invalid["sources"][0]["gene_claims"][0]["allowed_uses"] = ["clinical_safety"]
    with pytest.raises(ValueError, match="Unknown publication use"):
        literature.validate_publications(invalid)
    invalid["sources"][0]["gene_claims"][0]["allowed_uses"] = ["cta_nomination"]
    with pytest.raises(ValueError, match="explicit nomination"):
        literature.validate_publications(invalid)


def test_cell_and_isoform_study_claims_resolve_pinned_sources():
    from oncoref import reference_data

    sources = {s["source_id"]: s for s in literature.publication_registry()["sources"]}
    studies = {
        "10.1038/s41593-020-0587-9": ("GFY", "normal_tissue_review"),
        "10.1016/j.devcel.2023.07.013": ("AMELX", "normal_tissue_review"),
        "10.1038/s41586-024-08189-z": ("SP7", "normal_tissue_review"),
        "10.18632/aging.203452": ("SP7", "normal_tissue_review"),
        "10.1038/s41467-021-23758-w": ("CER1", "reproductive_support"),
        "10.1186/1471-2199-11-33": ("RPS4Y2", "reproductive_support"),
    }
    for doi, (gene, use) in studies.items():
        claim = next(c for c in sources[doi]["gene_claims"] if c["protein"] == gene)
        assert use in claim["allowed_uses"]
        assert "cta_nomination" not in claim["allowed_uses"]
        assert claim["provenance"]["reference_sources"]
        for key in claim["provenance"]["reference_sources"]:
            version = reference_data.resolve_version(key)
            assert reference_data.REFERENCE_SOURCES[key]["pins"][version]["sha256"]


def test_alarin_has_one_current_product_specific_claim():
    claims = literature.paper_gene_claims("GALP")
    alarin = claims[claims.source_id.eq("10.1385/jmn:29:2:145")]
    assert len(alarin) == 1
    assert "identity_review" in alarin.iloc[0].allowed_uses
    assert not alarin.iloc[0].eligible_for_coverage


def test_unscreened_adult_home_requires_a_scoped_human_claim():
    ledger = deepcopy(literature.publication_registry())
    claim = next(g for s in ledger["sources"] for g in s["gene_claims"] if g["protein"] == "GFY")
    assert claim["unscreened_adult_home"] == "olfactory epithelium"
    assert "olfactory epithelium" in set(
        literature.paper_gene_claims("GFY").unscreened_adult_home.dropna()
    )
    claim["unscreened_adult_home"] = "typo"
    with pytest.raises(ValueError, match="canonical unscreened tissue"):
        literature.validate_publications(ledger)

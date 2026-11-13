"""Coverage compression must never erase identity limits or adverse evidence."""

import gzip
import importlib.util
import json
from collections import Counter, defaultdict
from pathlib import Path

import pytest

from oncoref import cta_literature
from oncoref._cta_literature_followup import location_facet

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-literature-20261107"


@pytest.fixture(scope="module")
def data():
    return cta_literature.registry()


def test_attestations_preserve_every_link_and_its_exact_scope(data):
    from oncoref import cta_selection

    result = cta_literature.attestations()
    followup = cta_literature.followup_findings()
    assert result.groupby("origin").size().to_dict() == {
        "current_nomination": len(cta_selection.audit().memberships),
        "literature_registry": len(data["memberships"]),
        "literature_followup": len(followup["memberships"]),
    }
    assert result.source_url.str.startswith(("https://", "http://")).all()
    assert result.citation.str.len().gt(0).all()
    assert result.claim_scope.str.len().gt(0).all()
    assert (
        result[result.origin.eq("literature_registry")].membership_record.tolist()
        == data["memberships"]
    )
    assert (
        result[result.origin.eq("literature_followup")].membership_record.tolist()
        == followup["memberships"]
    )
    # A broad-program or exact biological paper is not a catalog nomination.
    program = result.membership_record.map(lambda r: r.get("kind") == "program_membership")
    assert program.any() and not result.loc[program, "published_nomination"].any()
    reproductive_only = result.source_roles.map(
        lambda roles: roles == ["reproductive_developmental"]
    )
    assert (
        reproductive_only.any() and not result.loc[reproductive_only, "published_nomination"].any()
    )
    assert not result.loc[result.origin.eq("literature_followup"), "published_nomination"].any()
    assert "supplementary" in set(result.collection)
    # Preserve historical Unsafe nominations and counterevidence as well as CTA links.
    unsafe = result.membership_record.map(lambda r: r.get("category") == "Unsafe")
    assert result.loc[unsafe, "published_nomination"].any()


def test_attestations_gene_filter_does_not_expand_related_copies():
    rows = cta_literature.attestations("ENSG00000183128")
    assert not rows.empty
    assert rows.gene_ids.map(lambda ids: "ENSG00000183128" in ids).all()
    correction = rows[rows.origin.eq("literature_followup")].iloc[0]
    assert "taste" in correction.claim_scope
    assert "limitations" in correction.source_record
    assert not correction.published_nomination
    assert cta_literature.attestations("unmapped-test-gene").empty


def test_background_table_is_not_a_published_tumor_target_nomination():
    background = cta_literature.attestations("KRTAP20-2")
    background = background[background.source_id.eq("10.18632/oncotarget.21715")]
    assert not background.empty and not background.published_nomination.any()
    assert all(not genes for genes in background.nomination_gene_ids)
    assert all(
        {c["source_tag"] for c in claims} == {"daSilva2017_testis_biased"}
        for claims in background.publication_claims
    )
    positive = cta_literature.attestations("OR1K1")
    positive = positive[positive.source_id.eq("10.18632/oncotarget.21715")]
    assert positive.published_nomination.all()
    assert all("ENSG00000165204" in genes for genes in positive.nomination_gene_ids)


def builder():
    spec = importlib.util.spec_from_file_location(
        "literature_builder", ROOT / "scripts/build_cta_literature_review.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_every_requested_protein_has_an_individual_review_and_search_receipts(data):
    population = json.loads((OUT / "review-population.json").read_text())
    reviews = {r["protein"]: r for r in data["reviews"]}
    assert set(reviews) == {r["label"] for r in population}
    assert len(reviews) == 96
    assert Counter(r["baseline_gap"] for r in reviews.values()) == {
        "no_match": 41,
        "broad_table_only": 55,
    }
    receipts = json.loads((OUT / "search-ledger.json").read_text())
    assert len(receipts) == 192
    for protein, review in reviews.items():
        searches = [r for r in receipts if r["protein"] == protein]
        assert {s["kind"] for s in searches} == {"gene_and_aliases", "fulltext_symbol_or_locus"}
        assert all(
            s["status"] == "retrieved" and s["sha256"] and s["url"].startswith("https://")
            for s in searches
        )
        assert review["summary"] and review["identity"] and review["review_disposition"]
        assert review["released_category"] == "CTA"


def test_core_preserves_all_eligible_coverage_and_is_inclusion_minimal(data):
    sets = defaultdict(set)
    for edge in data["memberships"]:
        if edge["eligible_for_coverage"]:
            sets[edge["source_id"]].add(edge["entity_id"])
    core = {s["source_id"] for s in data["sources"] if s["collection"] == "core"}
    full = set().union(*sets.values())
    assert set().union(*(sets[k] for k in core)) == full
    for key in core:
        assert set().union(*(sets[k] for k in core - {key})) != full
    for source in data["sources"]:
        if source["collection"] == "supplementary" and sets[source["source_id"]]:
            assert source["covered_by"]
            assert sets[source["source_id"]] <= set().union(
                *(sets[k] for k in source["covered_by"])
            )
        assert source["retention_reason"]


def test_original_catalog_and_cta_coverage_are_preserved(data):
    sets = defaultdict(set)
    cta_ids = {r["entity_id"] for r in data["cta_coverage"]}
    for edge in data["memberships"]:
        if edge["kind"] == "historical_catalog":
            sets[edge["source_id"]].add(edge["entity_id"])
    assert len(sets) == 15
    assert len(set().union(*sets.values())) == 2514
    catalog_core = data["counts"]["historical_catalog_core"]
    cta_core = data["counts"]["historical_cta_core"]
    assert len(catalog_core) == 12
    assert len(cta_core) == 8
    assert set().union(*(sets[k] for k in catalog_core)) == set().union(*sets.values())
    assert len(set().union(*(sets[k] for k in cta_core)) & cta_ids) == 480


def test_identity_species_and_proteoforms_are_not_silently_collapsed(data):
    reviews = {r["protein"]: r for r in data["reviews"]}
    for label in [
        "CIMIP3",
        "H2AL1Q",
        "DEFB113",
        "LCN9",
        "TEDDM1",
        "UBTFL1",
        "TRIM75",
        "PRAMEF18",
        "TAF11L11",
        "TAF11L13",
        "TAF11L14",
    ]:
        assert not reviews[label]["exact_human_paper_coverage"]
        assert reviews[label]["citations"]
    assert reviews["H3Y1"]["exact_human_paper_coverage"]
    assert "both variants" in reviews["H3Y1"]["summary"]
    assert reviews["FRG2"]["outcome"] == "identity_review"
    assert reviews["FRG2B"]["outcome"] == "identity_review"
    assert reviews["SLCO1B3-SLCO1B7"]["outcome"] == "normal_review"
    for label in ["BECN2", "SRY", "OBP2A", "ZACN"]:
        assert reviews[label]["outcome"] == "normal_review"
    assert data["classification_effect"].startswith("Evidence annotation only")


def test_complete_program_tables_keep_unmapped_rows_and_upstream_provenance(data):
    rows = data["program_rows"]
    assert len(rows) == 1505
    assert any(not r["matches"] and r["type"] != "gene" for r in rows)
    assert any(not r["matches"] and r["type"] == "gene" for r in rows)
    h3y = [r for r in rows if r["symbol"] == "H3.Y"]
    assert h3y and all({m["gene_id"] for m in r["matches"]} == {"ENSG00000269466"} for r in h3y)
    assert all(
        r["upstream_source_id"] == "10.1016/j.stem.2022.01.014" for r in rows if r["table"] == "S3A"
    )
    shared = [
        e
        for e in data["memberships"]
        if e["kind"] == "program_membership" and e["protein"] == "PRAMEF25/6"
    ]
    assert shared and all(e["gene_ids"] == ["ENSG00000229571"] for e in shared)
    assert data["background_screen"]["gene_rows"] == 12697
    assert not any(e.get("table") == "S1" for e in data["memberships"])


def test_human_paper_coverage_is_separate_from_normal_review_and_mouse_mechanism(data):
    coverage = {r["protein"]: r for r in data["cta_coverage"]}
    assert (
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC280673/" in coverage["SP7"]["exact_paper_sources"]
    )
    assert "https://pubmed.ncbi.nlm.nih.gov/11880645/" in coverage["PATE1"]["exact_paper_sources"]
    assert not coverage["PABPN1L"]["exact_paper_sources"]
    assert not coverage["CALHM3"]["exact_paper_sources"]
    pate = next(
        e
        for e in data["memberships"]
        if e["protein"] == "PATE1" and e["source_id"].endswith("11880645/")
    )
    assert "AF462605" in pate["coverage_rationale"]
    assert "https://www.ncbi.nlm.nih.gov/gene/160065" in pate["references"]
    sp7 = next(s for s in data["sources"] if s["source_id"].endswith("PMC280673/"))
    assert "normal_tissue" in sp7["roles"]
    assert "not detected" in sp7["retention_reason"]


def test_role_accessors_retain_noncore_normal_evidence_and_validate_filters(data):
    sources = {s["source_id"]: s for s in data["sources"]}
    assert len(sources) == len(data["sources"])
    for review in data["reviews"]:
        for citation in review["citations"]:
            assert citation["source_id"] in sources
            assert set(citation["roles"]) <= set(data["roles"])
            assert citation["pmid"].isdigit()
    assert not cta_literature.sources(collection="supplementary", role="normal_tissue").empty
    assert len(cta_literature.reviews("SRY")) == 1
    sry = cta_literature.reviews("SRY").iloc[0]
    assert len(cta_literature.reviews(sry.gene_ids[0])) == 1
    assert cta_literature.reviews("not-a-gene").empty
    with pytest.raises(ValueError):
        cta_literature.sources(collection="trusted")
    with pytest.raises(ValueError):
        cta_literature.sources(role="safe")


def test_registry_and_report_reproduce_offline(data):
    assert builder().build_registry() == data
    assert json.loads(gzip.decompress((OUT / "registry.json.gz").read_bytes())) == data
    assert len(data["cta_coverage"]) == 602
    assert (
        sum(bool(r["exact_paper_sources"]) for r in data["cta_coverage"])
        == data["counts"]["covered_ctas"]
    )
    assert len(cta_literature.memberships()) == len(data["memberships"])
    assert len(cta_literature.cta_coverage()) == 602


def test_gap_assessments_keep_coding_qualification_separate_from_missing_papers(data):
    gaps = {r["protein"]: r for r in data["gap_assessments"]}
    assert len(gaps) == 33
    assert all(r["translation_qualified"] for r in gaps.values())
    assert sum(r["normal_question"] for r in gaps.values()) == 5
    assert sum(r["exact_product_question"] for r in gaps.values()) == 2
    assert data["gap_comparison"]["paper_gap"]["routes"]["Developmental only"] == 18
    assert len(cta_literature.gap_assessments("CALHM3")) == 1
    concerns = {r["protein"]: r for r in data["concerns"]}
    assert len(concerns) == 29
    assert Counter(r["outcome"] for r in concerns.values()) == {
        "normal_review": 18,
        "identity_review": 11,
    }
    assert sum(r["concern_group"] == "little_characterization" for r in concerns.values()) == 5
    assert concerns["FRG2"]["concern_group"] == "coding_conflict"
    assert concerns["ENSG00000273155"]["concern_group"] == "shared_identity"
    assert len(cta_literature.concerns("SRY")) == 1


@pytest.mark.parametrize(
    "name,topology,orientation,expected",
    [
        ("Secreted", "", "", "secreted"),
        ("Cytoplasmic vesicle, secretory vesicle, acrosome", "", "", "intracellular"),
        (
            "Endoplasmic reticulum membrane",
            "Multi-pass membrane protein",
            "",
            "intracellular_membrane",
        ),
        ("Membrane", "Multi-pass membrane protein", "", "membrane_unspecified"),
        ("Cell membrane", "Peripheral membrane protein", "", "cell_membrane_unresolved"),
        ("Cell membrane", "Lipid-anchor", "Cytoplasmic side", "cell_membrane_unresolved"),
        ("Cell membrane", "GPI-anchor", "", "cell_surface"),
        ("Cell membrane", "Multi-pass membrane protein", "", "cell_surface"),
        ("Cell surface", "", "", "cell_surface"),
    ],
)
def test_localization_does_not_turn_membranes_or_secretory_vesicles_into_surface_evidence(
    name, topology, orientation, expected
):
    assert (
        location_facet(
            {
                "location": {"value": name},
                "topology": {"value": topology},
                "orientation": {"value": orientation},
            }
        )
        == expected
    )


def test_localization_keeps_shared_identity_model_and_evidence_scope(data):
    locations = {r["protein"]: r for r in data["localization"]}
    assert len(locations) == 602
    for protein in ["OOSP4B", "EPPIN-WFDC6"]:
        assert "secreted" in locations[protein]["shared_facets"]
        assert "secreted" not in locations[protein]["facets"]
    assert "cell_surface" not in locations["PRAME"]["facets"]
    assert "cell_membrane_unresolved" in locations["PRAME"]["shared_facets"]
    assert locations["SPADH"]["facets"] == ["unknown"]
    assert locations["CALHM3"]["annotations"][0]["location_evidence"] == ["By similarity"]
    assert "cell_surface" in locations["ZP4"]["facets"]
    assert "extracellular_detection" in locations["CRIPTO3"]["facets"]
    assert "cell_surface" not in locations["CRIPTO3"]["facets"]
    readthrough = locations["SLCO1B3-SLCO1B7"]
    assert "intracellular_membrane" in readthrough["facets"]
    assert any("Native human hepatocytes" in n["scope"] for n in readthrough["paper_notes"])
    assert all(not r["sequence_localization_verified"] for r in locations.values())
    assert len(cta_literature.localization("ENSG00000255393")) == 1
    assert cta_literature.localization("not-a-gene").empty

"""Coverage compression must never erase identity limits or adverse evidence."""

import gzip
import importlib.util
import json
from collections import Counter, defaultdict
from pathlib import Path

import pytest

from oncoref import cta_literature

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-literature-20261107"


@pytest.fixture(scope="module")
def data():
    return cta_literature.registry()


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

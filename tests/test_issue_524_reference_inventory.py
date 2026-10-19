"""Regressions for searchable identifiers and source-document provenance."""

import csv
import io
import json
from collections import Counter
from html.parser import HTMLParser

import build_issue_524_reference_inventory as inventory
import pytest


class BibliographyParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.cards = {}
        self.card = None
        self.link = None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "article":
            self.card = {"text": "", "links": [], "group": attrs["data-group"]}
            self.cards[attrs["id"]] = self.card
        elif tag == "a" and self.card is not None:
            self.link = {"url": attrs["href"], "label": ""}
            self.card["links"].append(self.link)

    def handle_endtag(self, tag):
        if tag == "article":
            self.card = None
        elif tag == "a":
            self.link = None

    def handle_data(self, data):
        if self.card is not None:
            self.card["text"] += data
        if self.link is not None:
            self.link["label"] += data


@pytest.fixture(scope="module")
def report():
    artifacts = inventory.build_inventory()
    parser = BibliographyParser()
    parser.feed(artifacts["bibliography.html"])
    return {
        "artifacts": artifacts,
        "cards": parser.cards,
        "references": list(csv.DictReader(io.StringIO(artifacts["references.csv"]))),
        "records": list(csv.DictReader(io.StringIO(artifacts["flagged-records.csv"]))),
        "summary": json.loads(artifacts["lookup-audit.json"]),
    }


@pytest.mark.parametrize("field", ["pmid", "doi", "pmcid"])
def test_all_known_identifiers_are_visible_and_searchable(report, field):
    # The browser filters textContent, so identifiers in href alone do not qualify.
    for reference in report["references"]:
        if reference[field]:
            text = report["cards"][reference["reference_id"]]["text"]
            assert f"{field.upper()}:{reference[field]}" in text, reference["reference"]
            assert reference[field].lower() in text.lower(), reference["reference"]


def test_all_source_document_urls_are_rendered_once(report):
    for reference in report["references"]:
        urls = [link["url"] for link in report["cards"][reference["reference_id"]]["links"]]
        normalized = [url.lower() if url.startswith("https://doi.org/") else url for url in urls]
        for field in ("full_text_or_existing_document_urls", "followup_source_urls"):
            for url in reference[field].split(" | "):
                if url:
                    key = url.lower() if url.startswith("https://doi.org/") else url
                    assert key in normalized, reference["reference"]
        assert len(normalized) == len(set(normalized)), reference["reference"]


def test_dart_pmc_document_keeps_existing_locator_provenance(report):
    identifier = "DOI:10.1158/1078-0432.CCR-19-3356"
    url = "https://pmc.ncbi.nlm.nih.gov/articles/PMC7231627/"
    reference = next(row for row in report["references"] if row["reference"] == identifier)
    card = report["cards"][reference["reference_id"]]
    assert {"url": url, "label": "Existing source document"} in card["links"]
    assert reference["pmcid"] == ""
    assert reference["metadata_provider"] == "Crossref DOI registry"
    assert "PMCID:" not in card["text"]
    assert reference["numeric_claim_review_status"] == "unresolved_audit_flags"
    records = [row for row in report["records"] if row["reference"] == identifier]
    assert records
    for row in records:
        assert row["source_document_url"] == url
        assert row["source_document_kind"] == "pmc_full_text"
        assert "numeric_source_location_not_verified" in row["findings"].split(";")


def test_pmc_counts_distinguish_snapshot_identifiers_from_existing_urls(report):
    assert report["summary"]["references_with_pmc"] == 61
    assert report["summary"]["references_with_pmc_document_links"] == 62
    assert "61 references with snapshot PMCIDs" in report["artifacts"]["bibliography.html"]
    assert "62 references with PMC document links" in report["artifacts"]["bibliography.html"]


def test_scientific_review_annotations_and_bibliographic_snapshot_are_preserved(report):
    ledger = inventory.read_csv(inventory.ROOT / inventory.SOURCE_FILES[0])
    selected = {
        (row["dataset"], row["row_key"]): row
        for row in ledger
        if set(row["findings"].split(";")) & inventory.FLAGS
    }
    assert len(report["records"]) == len(selected) == 435
    for row in report["records"]:
        original = selected[(row["dataset"], row["row_key"])]
        assert {key: row[key] for key in original} == original
    flags = Counter(
        flag
        for row in report["records"]
        for flag in row["findings"].split(";")
        if flag in inventory.FLAGS
    )
    assert (
        flags
        == report["summary"]["flags"]
        == {
            "tmb_source_not_revalidated": 67,
            "numeric_source_location_not_verified": 338,
            "response_denominator_missing": 33,
            "small_denominator": 47,
        }
    )
    snapshot = json.loads((inventory.DEFAULT_DIRECTORY / "bibliographic-metadata.json").read_text())
    metadata = {row["reference"]: row for row in snapshot}
    for row in report["references"]:
        for field in ("pmid", "doi", "pmcid", "metadata_provider", "metadata_status"):
            assert row[field] == metadata[row["reference"]][field]
        assert row["numeric_claim_review_status"] == "unresolved_audit_flags"


def test_checked_in_inventory_reproduces(report):
    for name, content in report["artifacts"].items():
        assert (inventory.DEFAULT_DIRECTORY / name).read_bytes() == content.encode("utf-8")

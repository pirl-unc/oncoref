"""Protect the boundary between retrieved literature and accepted findings."""

import concurrent.futures
import gzip
import hashlib
import json
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import pytest
from scripts.audit_cta_literature import SearchClient
from scripts.build_cta_literature_audit import (
    build_rows,
    load_adjudications,
    search_evidence,
)

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/audits/cta-literature-20261025"


def test_duplicate_queries_share_one_immutable_response(tmp_path):
    response = Mock()
    response.url = "https://example.test/public-search"
    response.json.return_value = {"hitCount": 0, "resultList": {"result": []}}
    client = SearchClient(tmp_path, interval=0)
    with (
        patch("scripts.audit_cta_literature.requests.get", return_value=response) as get,
        concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor,
    ):
        results = list(executor.map(lambda _: client.fetch("same family query", 20), range(8)))
    assert get.call_count == 1
    assert len({record["retrieved_utc"] for record, _ in results}) == 1
    assert len({path for _, path in results}) == 1


def test_search_hits_remain_unadjudicated_and_cache_tampering_fails(tmp_path):
    cache = tmp_path / "paper.json.gz"
    with gzip.open(cache, "wt") as handle:
        json.dump(
            {
                "response": {
                    "resultList": {"result": [{"id": "123", "title": "Normal neural protein"}]}
                }
            },
            handle,
        )
    row = {
        "Gene": "FTMT",
        "Gene ID": "id1",
        "Category": "CTA",
        "Status": "Retrieved",
        "Search": "normal",
        "Cache file": cache.name,
        "SHA256": hashlib.sha256(cache.read_bytes()).hexdigest(),
    }
    evidence = search_evidence(pd.DataFrame([row]), tmp_path)
    assert evidence.iloc[0]["Category"] == "CTA"
    assert evidence.iloc[0]["Review"].startswith("Retrieved lead")
    cache.write_bytes(b"changed")
    with pytest.raises(ValueError, match="checksum mismatch"):
        search_evidence(pd.DataFrame([row]), tmp_path)


def test_every_gene_has_a_provisional_row_without_silent_category_changes():
    atlas = pd.read_csv(
        ROOT / "docs/audits/cta-atlas-20261024/all-genes.csv", keep_default_na=False
    )
    decisions = load_adjudications(AUDIT / "adjudications.yml", set(atlas.Gene))
    columns = ["Gene ID", "Status", "Search", "Hits", "Truncated"]
    result = build_rows(atlas, pd.DataFrame(columns=columns), decisions, {})
    assert len(result) == 2540
    assert result["Gene ID"].is_unique
    assert result.Category.equals(atlas.Category)
    assert result["Search coverage"].eq("Incomplete").all()
    ftmt = result.loc[result.Gene.eq("FTMT")].iloc[0]
    assert ftmt.Category == "CTA" and ftmt["Proposed category"] == "Unsafe"
    assert ftmt["Change proposed"]
    # Unknown literature is not converted into a safety endorsement.
    unreviewed = result.loc[~result.Gene.isin(decisions)]
    assert unreviewed["Assessment basis"].str.contains("not individually adjudicated").all()


def test_neither_cannot_hide_a_hard_normal_rna_failure():
    atlas = pd.read_csv(
        ROOT / "docs/audits/cta-atlas-20261024/all-genes.csv", keep_default_na=False
    )
    atlas.loc[atlas.Category.eq("Neither"), "Normal RNA maximum (TPM)"] = "30.1"
    columns = ["Gene ID", "Status", "Search", "Hits", "Truncated"]
    with pytest.raises(ValueError, match="hard normal-RNA limit"):
        build_rows(atlas, pd.DataFrame(columns=columns), {}, {})


def test_adjudications_reject_duplicate_or_unmapped_gene_names(tmp_path):
    path = tmp_path / "adjudications.yml"
    path.write_text("records:\n  - genes: [FTMT, FTMT]\n    proposed: Unsafe\n")
    with pytest.raises(ValueError, match="duplicate"):
        load_adjudications(path, {"FTMT"})
    with pytest.raises(ValueError, match="Unknown"):
        load_adjudications(path, {"ABCC11"})

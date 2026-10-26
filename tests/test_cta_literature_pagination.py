import json
from unittest.mock import Mock

import pytest
from expand_cta_literature import qualify, retrieve


def client(tmp_path, responses):
    c = Mock()
    records = []
    for i, data in enumerate(responses):
        record = {"response": data, "url": "https://example.test", "retrieved_utc": "test"}
        p = tmp_path / f"{i}.json"
        p.write_text(json.dumps(record))
        records.append((record, p))
    c.fetch.side_effect = records
    return c


def page(ids, total, cursor=None):
    return {
        "hitCount": total,
        "resultList": {"result": [{"id": str(i)} for i in ids]},
        "nextCursorMark": cursor,
    }


def test_cursor_pagination_and_raw_receipts(tmp_path):
    c = client(tmp_path, [page([1, 2], 3, "next"), page([2, 3], 3)])
    papers, pages, total, complete = retrieve(c, "query")
    assert len(papers) == total == 3 and complete and len(pages) == 2
    assert all(len(p["sha256"]) == 64 for p in pages)


@pytest.mark.parametrize(
    "responses",
    [
        [page([1], 3, "next"), page([2], 3, "next")],
        [page([1], 3, "next"), page([2], 4)],
        [page([1], 3)],
    ],
)
def test_broken_cursor_or_changed_index_never_reports_completion(tmp_path, responses):
    with pytest.raises(ValueError):
        retrieve(client(tmp_path, responses), "query")


def test_page_cap_is_incomplete(tmp_path):
    _, _, _, complete = retrieve(client(tmp_path, [page([1], 3, "next")]), "q", max_pages=1)
    assert not complete


def test_ambiguous_aliases_are_contextualized_without_dropping_canonical_identity():
    result = qualify('TITLE_ABS:"FATE" OR TITLE_ABS:"FATE1" OR TITLE_ABS:"MRP8"')
    assert 'TITLE_ABS:"FATE"' not in result
    assert 'TITLE_ABS:"FATE1"' in result
    assert "transporter" in result
    assert qualify('TITLE_ABS:"MAGEA4"') == 'TITLE_ABS:"MAGEA4"'

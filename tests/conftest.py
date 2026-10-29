import hashlib
import io
import json
from pathlib import Path

import pytest
import requests


@pytest.fixture
def historical_cta_input():
    """Resolve exact frozen input bytes without requiring Git history or network."""
    root = Path(__file__).resolve().parents[1]
    archive = root / "tests/fixtures/cta-1.21-inputs"
    receipt = json.loads((archive / "manifest.json").read_text())

    def resolve(name, expected):
        current = root / name
        if hashlib.sha256(current.read_bytes()).hexdigest() == expected:
            return current
        assert receipt["sha256"].get(name) == expected, name
        frozen = archive / name
        assert hashlib.sha256(frozen.read_bytes()).hexdigest() == expected, name
        return frozen

    return resolve


@pytest.fixture
def legacy_cta_membership():
    """Nomination inputs matching the deliberately frozen legacy candidate table."""
    import pandas as pd

    root = Path(__file__).parent / "fixtures/cta-1.21-inputs"
    name = "oncoref/data/cta-publication-membership.csv"
    receipt = json.loads((root / "manifest.json").read_text())
    assert hashlib.sha256((root / name).read_bytes()).hexdigest() == receipt["sha256"][name]
    return pd.read_csv(root / name)


@pytest.fixture
def http_downloads(monkeypatch):
    """Supply HTTP response bodies while exercising the real datacache transfer."""

    def install(body_for_url):
        def get(url, **kwargs):
            response = requests.Response()
            response.status_code = 200
            response.url = url
            body = body_for_url(url)
            response.raw = io.BytesIO(body)
            response.headers["Content-Length"] = str(len(body))
            return response

        monkeypatch.setattr(requests, "get", get)

    return install

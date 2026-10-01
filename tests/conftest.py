import io

import pytest
import requests


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

"""Retries of HTTP requests after a transient failure (gwtc_analysis.net)."""
import io

import pytest
import requests

from gwtc_analysis import net

URL = "https://zenodo.org/api/records/22685054/versions"


def response(status, headers=None):
    r = requests.Response()
    r.status_code, r.url = status, URL
    r.headers.update(headers or {})
    r._content, r.raw = b"{}", io.BytesIO(b"")
    return r


@pytest.fixture
def waits(monkeypatch):
    w = []
    monkeypatch.setattr(net.time, "sleep", w.append)
    return w


@pytest.fixture
def server(monkeypatch):
    """requests.request answering the statuses queued in the returned list, in turn."""
    queue, calls = [], []

    def request(method, url, **kwargs):
        calls.append((method, url))
        s = queue.pop(0)
        if isinstance(s, Exception):
            raise s
        return response(*s) if isinstance(s, tuple) else response(s)

    monkeypatch.setattr(requests, "request", request)
    return queue, calls


def test_503_is_retried_then_succeeds(server, waits):
    """Zenodo's API answered 503 in the tools-iuc CI (zenodo_releases test)."""
    queue, calls = server
    queue += [503, 502, 200]
    logs = []
    assert net.get(URL, log=logs.append).status_code == 200
    assert len(calls) == 3 and waits == [5, 15]
    assert "HTTP 503" in logs[0] and "zenodo.org" in logs[0]


def test_retry_after_is_followed_up_to_a_limit(server, waits):
    queue, _ = server
    queue += [(429, {"Retry-After": "2"}), (503, {"Retry-After": "3600"}), 200]
    net.get(URL, log=lambda m: None)
    assert waits == [2, 60]


def test_dns_failure_is_retried(server, waits):
    queue, _ = server
    queue += [requests.ConnectionError("Failed to resolve 'dcc.ligo.org'"), 200]
    assert net.get(URL, log=lambda m: None).status_code == 200 and waits == [5]


def test_client_error_is_left_to_the_caller_at_once(server, waits):
    queue, calls = server
    queue += [404]
    r = net.get(URL, log=lambda m: None)
    assert r.status_code == 404 and len(calls) == 1 and waits == []
    with pytest.raises(requests.HTTPError):
        r.raise_for_status()


def test_persistent_503_is_raised_after_the_last_try(server, waits):
    queue, calls = server
    queue += [503, 503, 503]
    with pytest.raises(requests.HTTPError) as e:
        net.get(URL, log=lambda m: None)
    assert e.value.response.status_code == 503 and len(calls) == 3


def test_post_is_retried(server, waits):
    queue, calls = server
    queue += [504, 200]
    assert net.post(URL, data={"q": 1}, log=lambda m: None).status_code == 200
    assert calls == [("POST", URL)] * 2


def test_retry_of_a_whole_operation(waits):
    """A download interrupted mid-file (ChunkedEncodingError) is started again; a bug is not retried."""
    tries = []

    def download():
        tries.append(1)
        if len(tries) == 1:
            raise requests.exceptions.ChunkedEncodingError("connection broken")
        return "done"

    assert net.retry(download, URL, log=lambda m: None) == "done" and len(tries) == 2
    with pytest.raises(KeyError):
        net.retry(lambda: {}["x"], URL, log=lambda m: None)
    assert waits == [5]

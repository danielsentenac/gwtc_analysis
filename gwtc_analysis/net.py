"""HTTP requests retried after a transient failure.

A connection error or a timeout (a DNS failure, a reset), or a status the server sends while it is overloaded or
throttling (429, 500, 502, 503, 504), is retried twice: after the delay the server asks for (Retry-After, at most
60 s), else after 5 s then 15 s. Any other HTTP error (404, 403...) is raised at once: retrying would not change it.

`retry` wraps a whole operation (a download written to a file, so an interruption mid-file is retried too);
`get` and `post` wrap a single request.
"""
from __future__ import annotations

import time
from typing import Callable, Optional, TypeVar
from urllib.parse import urlparse

import requests

RETRY_STATUS = frozenset({429, 500, 502, 503, 504})
TRIES = 3
PAUSE_S = 5.0
MAX_WAIT_S = 60.0

T = TypeVar("T")


def _status(e: BaseException) -> Optional[int]:
    r = getattr(e, "response", None)
    return getattr(r, "status_code", None)


def transient(e: BaseException) -> bool:
    """Whether retrying after `e` may succeed."""
    if isinstance(e, requests.HTTPError):
        return _status(e) in RETRY_STATUS
    return isinstance(e, (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError))


def _retry_after(e: BaseException) -> Optional[float]:
    """The delay of the Retry-After header of an HTTP error, in seconds (the HTTP-date form is ignored)."""
    r = getattr(e, "response", None)
    v = (getattr(r, "headers", None) or {}).get("Retry-After")
    try:
        return max(0.0, float(v)) if v is not None else None
    except ValueError:
        return None


def _why(e: BaseException) -> str:
    s = _status(e)
    return f"HTTP {s}" if isinstance(e, requests.HTTPError) and s else type(e).__name__


def retry(fn: Callable[[], T], url: str, *, log: Callable[[str], None] = print, tries: int = TRIES,
          pause: float = PAUSE_S, max_wait: float = MAX_WAIT_S) -> T:
    """fn() retried after a transient failure (see the module docstring)."""
    for i in range(tries):
        try:
            return fn()
        except Exception as e:
            if i == tries - 1 or not transient(e):
                raise
            asked = _retry_after(e)
            wait = min(asked if asked is not None else pause * 3 ** i, max_wait)
            log(f"⚠️ [RETRY] {urlparse(url).netloc}: {_why(e)}, retrying in {wait:g} s")
            time.sleep(wait)
    raise AssertionError("unreachable")


def _send(method: str, url: str, session: Optional[requests.Session], **kwargs) -> requests.Response:
    r = (session or requests).request(method, url, **kwargs)
    if r.status_code in RETRY_STATUS:
        r.close()
        r.raise_for_status()
    return r


def get(url: str, *, log: Callable[[str], None] = print, session: Optional[requests.Session] = None,
        **kwargs) -> requests.Response:
    """requests.get(url, **kwargs), retried after a transient failure. Other HTTP errors are left to the caller
    (raise_for_status), as with requests.get."""
    return retry(lambda: _send("GET", url, session, **kwargs), url, log=log)


def post(url: str, *, log: Callable[[str], None] = print, session: Optional[requests.Session] = None,
         **kwargs) -> requests.Response:
    """requests.post(url, **kwargs), retried after a transient failure (for idempotent queries only)."""
    return retry(lambda: _send("POST", url, session, **kwargs), url, log=log)

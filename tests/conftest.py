"""Shared pytest helpers that make the smoke tests self-explaining.

Each test prints (1) its purpose, taken from the docstring, before it runs, and
(2) step-by-step lines via the ``explain`` fixture. Both are written straight to
the terminal reporter so they show up regardless of stdout capture.
"""
from __future__ import annotations

import pytest


def _reporter(config):
    return config.pluginmanager.getplugin("terminalreporter")


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_call(item):
    """Print the test's purpose (its docstring) just before it executes."""
    tr = _reporter(item.config)
    func = getattr(item, "function", None)
    doc = (func.__doc__ or "").strip() if func is not None else ""
    if tr is not None and doc:
        tr.ensure_newline()
        tr.write_line(f"\n▶ {item.nodeid}", bold=True)
        for line in doc.splitlines():
            tr.write_line(f"  {line.strip()}")
    yield


@pytest.fixture
def explain(request):
    """Return ``explain(msg)`` to narrate what a test step is doing."""
    tr = _reporter(request.config)

    def _explain(message: str) -> None:
        if tr is not None:
            tr.ensure_newline()
            tr.write_line(f"    · {message}", cyan=True)
        else:  # pragma: no cover - fallback when no terminal reporter
            print(f"    · {message}")

    return _explain

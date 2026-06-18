"""Smoke tests for the ``zenodo`` data-repo (public Zenodo PE records).

Offline test: the CLI accepts ``--data-repo zenodo``.
Network test: build the Zenodo PE index and HEAD-check a resolved download link
(no multi-GB download); skips cleanly when zenodo.org is unreachable.
"""
from __future__ import annotations

import pytest
import requests

from gwtc_analysis import parameters_estimation as pe
from gwtc_analysis.cli import build_parser


def test_cli_pe_accepts_zenodo_data_repo(explain):
    """CLI wiring: parameters_estimation accepts --data-repo zenodo (offline)."""
    explain("Parsing CLI args: parameters_estimation --src-name GW240420_175625 --data-repo zenodo")
    args = build_parser().parse_args(
        ["parameters_estimation", "--src-name", "GW240420_175625", "--data-repo", "zenodo"]
    )
    explain(f"Parsed mode={args.mode!r}, data_repo={args.data_repo!r}")
    assert args.mode == "parameters_estimation"
    assert args.data_repo == "zenodo"


@pytest.mark.network
def test_zenodo_pe_index_resolves(tmp_path, explain):
    """Resolve a real PE link from public Zenodo records and HEAD-check it (no download)."""
    explain("Building the Zenodo PE index from the configured PE records")
    try:
        index = pe.build_zenodo_pe_index(cache_dir=str(tmp_path), force_refresh=True)
    except requests.RequestException as e:
        pytest.skip(f"zenodo.org unreachable: {e}")
    explain(f"Index built: {len(index)} events carry PE files on Zenodo")
    assert index, "expected a non-empty Zenodo PE index"

    event = "GW240420_175625" if "GW240420_175625" in index else next(iter(index))
    chosen = pe.choose_best_pe_file(index[event])
    explain(f"Best PE file for {event}: {chosen['filename'] if chosen else None}")
    assert chosen is not None
    assert chosen["url"].startswith("https://zenodo.org/")

    explain(f"HEAD-checking the download URL (no body): {chosen['url']}")
    try:
        resp = requests.head(chosen["url"], allow_redirects=True, timeout=60)
    except requests.RequestException as e:
        pytest.skip(f"zenodo.org unreachable: {e}")
    explain(f"Server replied HTTP {resp.status_code}, content-length={resp.headers.get('content-length')}")
    assert resp.status_code == 200

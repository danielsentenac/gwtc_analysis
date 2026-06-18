"""Smoke tests for the ``galaxy`` data-repo (usegalaxy.org published history).

Offline tests: the CLI accepts ``--data-repo galaxy`` and the URL helpers build
correct REST/download endpoints.
Network test: resolve real PE download links from the public Galaxy history and
HEAD-check one (no multi-GB download); skips cleanly when usegalaxy.org is down.
"""
from __future__ import annotations

import pytest
import requests

from gwtc_analysis import data_repo
from gwtc_analysis import parameters_estimation as pe
from gwtc_analysis.cli import build_parser


def test_cli_pe_accepts_galaxy_data_repo(explain):
    """CLI wiring: parameters_estimation accepts --data-repo galaxy (offline)."""
    explain("Parsing CLI args: parameters_estimation --src-name GW240420_175625 --data-repo galaxy")
    args = build_parser().parse_args(
        ["parameters_estimation", "--src-name", "GW240420_175625", "--data-repo", "galaxy"]
    )
    explain(f"Parsed mode={args.mode!r}, data_repo={args.data_repo!r}")
    assert args.mode == "parameters_estimation"
    assert args.data_repo == "galaxy"


def test_galaxy_url_helpers(explain):
    """URL helpers build the expected REST/download endpoints (offline)."""
    explain("Checking galaxy_server_url / galaxy_history_id / galaxy_contents_url")
    assert data_repo.galaxy_server_url().startswith("http")
    assert data_repo.galaxy_history_id()
    assert data_repo.galaxy_contents_url().endswith("/contents")
    assert data_repo.galaxy_history_id() in data_repo.galaxy_contents_url()

    url = data_repo.galaxy_dataset_download_url("ABC123", "h5")
    explain(f"galaxy_dataset_download_url('ABC123','h5') -> {url}")
    assert url == f"{data_repo.galaxy_server_url()}/api/datasets/ABC123/display?to_ext=h5"


@pytest.mark.network
def test_galaxy_history_pe_links_resolve(tmp_path, explain):
    """Resolve PE links from the public history and HEAD-check one (no download)."""
    explain("Building the Galaxy PE index from the public history's REST contents")
    try:
        index = pe.build_galaxy_pe_index(cache_dir=str(tmp_path), force_refresh=True)
    except requests.RequestException as e:
        pytest.skip(f"usegalaxy.org unreachable: {e}")
    explain(f"Index built: {len(index)} events with PE files in the history")
    assert index, "expected a non-empty PE index from the public Galaxy history"

    event = next(iter(index))
    chosen = pe.choose_best_pe_file(index[event])
    explain(f"Best PE file for {event}: {chosen['name'] if chosen else None}")
    assert chosen is not None
    assert chosen["name"].endswith((".h5", ".hdf5"))
    assert chosen["url"].startswith(data_repo.galaxy_server_url())

    explain(f"HEAD-checking the download URL (no body): {chosen['url']}")
    try:
        resp = requests.head(chosen["url"], allow_redirects=True, timeout=60)
    except requests.RequestException as e:
        pytest.skip(f"usegalaxy.org unreachable: {e}")
    explain(f"Server replied HTTP {resp.status_code}, content-length={resp.headers.get('content-length')}")
    assert resp.status_code == 200

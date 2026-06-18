"""Smoke tests for the ``galaxy`` data-repo (usegalaxy.org published history).

The offline tests check CLI wiring and URL construction. The network test
resolves real PE download links from the public Galaxy history and HEAD-checks
one of them (no multi-GB download); it skips cleanly when usegalaxy.org is
unreachable.
"""
from __future__ import annotations

import pytest
import requests

from gwtc_analysis import data_repo
from gwtc_analysis import parameters_estimation as pe
from gwtc_analysis.cli import build_parser


def test_cli_pe_accepts_galaxy_data_repo():
    """The parameters_estimation CLI must accept ``--data-repo galaxy``."""
    args = build_parser().parse_args(
        ["parameters_estimation", "--src-name", "GW240420_175625", "--data-repo", "galaxy"]
    )
    assert args.mode == "parameters_estimation"
    assert args.data_repo == "galaxy"


def test_galaxy_url_helpers():
    """URL helpers build the expected REST/download endpoints."""
    assert data_repo.galaxy_server_url().startswith("http")
    assert data_repo.galaxy_history_id()
    assert data_repo.galaxy_contents_url().endswith("/contents")
    assert data_repo.galaxy_history_id() in data_repo.galaxy_contents_url()

    url = data_repo.galaxy_dataset_download_url("ABC123", "h5")
    assert url == f"{data_repo.galaxy_server_url()}/api/datasets/ABC123/display?to_ext=h5"


@pytest.mark.network
def test_galaxy_history_pe_links_resolve(tmp_path):
    """Resolve PE links from the public history and HEAD-check one (no download)."""
    try:
        index = pe.build_galaxy_pe_index(cache_dir=str(tmp_path), force_refresh=True)
    except requests.RequestException as e:  # offline / server down
        pytest.skip(f"usegalaxy.org unreachable: {e}")

    assert index, "expected a non-empty PE index from the public Galaxy history"

    event = next(iter(index))
    chosen = pe.choose_best_pe_file(index[event])
    assert chosen is not None
    assert chosen["name"].endswith((".h5", ".hdf5"))
    assert chosen["url"].startswith(data_repo.galaxy_server_url())

    try:
        resp = requests.head(chosen["url"], allow_redirects=True, timeout=60)
    except requests.RequestException as e:
        pytest.skip(f"usegalaxy.org unreachable: {e}")
    assert resp.status_code == 200

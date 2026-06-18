"""Smoke tests for the ``s3`` data-repo (MinIO bucket 'gwtc').

Offline test: the CLI accepts ``--data-repo s3``.
Network test: list the 'gwtc' bucket and confirm a PEDataRelease file is present.
The bucket is private, so the network test requires S3 credentials in the
``S3_CREDENTIALS`` env var and skips cleanly when they are absent or listing fails.
"""
from __future__ import annotations

import json
import os

import pytest

from gwtc_analysis.cli import build_parser

DEFAULT_ENDPOINT = "minio-dev.odahub.fr"


def test_cli_pe_accepts_s3_data_repo(explain):
    """CLI wiring: parameters_estimation accepts --data-repo s3 (offline)."""
    explain("Parsing CLI args: parameters_estimation --src-name GW240420_175625 --data-repo s3")
    args = build_parser().parse_args(
        ["parameters_estimation", "--src-name", "GW240420_175625", "--data-repo", "s3"]
    )
    explain(f"Parsed mode={args.mode!r}, data_repo={args.data_repo!r}")
    assert args.mode == "parameters_estimation"
    assert args.data_repo == "s3"


def _make_minio_client():
    """Build a MinIO client the same way run_parameters_estimation does."""
    from minio import Minio

    env = os.environ.get("S3_CREDENTIALS")
    creds = json.loads(env) if env else {"endpoint": DEFAULT_ENDPOINT, "secure": True}
    return Minio(
        endpoint=creds["endpoint"],
        secure=creds.get("secure", True),
        access_key=creds.get("access_key"),
        secret_key=creds.get("secret_key"),
    )


@pytest.mark.network
def test_s3_bucket_lists_pe_files(explain):
    """List the MinIO 'gwtc' bucket and confirm a PEDataRelease .h5/.hdf5 object exists.

    Requires S3_CREDENTIALS (private bucket); skips if absent or listing is denied.
    """
    if not os.environ.get("S3_CREDENTIALS"):
        explain("S3_CREDENTIALS not set; the 'gwtc' bucket is private, so skipping")
        pytest.skip("set S3_CREDENTIALS to run the S3 network smoke test")

    explain(f"Connecting to MinIO endpoint and listing bucket 'gwtc' recursively")
    client = _make_minio_client()
    try:
        found = None
        for obj in client.list_objects("gwtc", recursive=True):
            name = obj.object_name
            if "PEDataRelease" in name and name.endswith((".h5", ".hdf5")):
                found = name
                break
    except Exception as e:  # network / auth / S3Error
        pytest.skip(f"S3 listing failed: {type(e).__name__}: {e}")

    explain(f"First PEDataRelease object found: {found}")
    assert found is not None, "no PEDataRelease .h5/.hdf5 object found in bucket 'gwtc'"

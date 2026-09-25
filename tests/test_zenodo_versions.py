"""Offline tests for Zenodo release-version resolution (``data_repo``).

The Zenodo API is replaced by a fake version listing, so these tests need no
network: latest-by-default, explicit vN, errors, and the offline fallbacks.
"""
from __future__ import annotations

import json

import pytest
import requests

from gwtc_analysis import data_repo as dr
from gwtc_analysis.cli import build_parser, _parse_zenodo_versions

FAKE_GWTC3 = [  # oldest first, as returned by data_repo._fetch_versions
    {"record_id": "5546663", "publication_date": "2021-11-08",
     "files": [{"key": "skymaps.tar.gz"}, {"key": "contour_data.tar.gz"}]},
    {"record_id": "8177023", "publication_date": "2023-10-23",
     "files": [{"key": "IGWN-GWTC3p0-v2-PESkyLocalizations.tar.gz"},
               {"key": "IGWN-GWTC3p0-v2-PEContours.tar.gz"}]},
    {"record_id": "22685054", "publication_date": "2026-09-21",
     "files": [{"key": "IGWN-GWTC3p0-v3-PEContours.tar.gz"},
               {"key": "IGWN-GWTC3p0-v3-PESkyLocalizations.tar.gz"}]},
]


@pytest.fixture
def fake_zenodo(monkeypatch, tmp_path):
    """Serve FAKE_GWTC3 for GWTC-3 and cache listings under tmp_path."""
    calls = []

    def fetch(record_id):
        calls.append(record_id)
        return [dict(v) for v in FAKE_GWTC3]

    monkeypatch.setattr(dr, "_fetch_versions", fetch)
    monkeypatch.setattr(dr, "zenodo_cache_dir", lambda: tmp_path)
    return calls


@pytest.mark.parametrize("raw, expected", [
    (None, None), ("latest", None), ("", None), ("v2", 2), ("2", 2), (3, 3), ("V1", 1),
])
def test_parse_zenodo_version(raw, expected):
    """Version strings: 'latest' means None, vN / N mean N."""
    assert dr.parse_zenodo_version(raw) == expected


@pytest.mark.parametrize("raw", ["v0", "beta", "v2.1", "-1"])
def test_parse_zenodo_version_rejects(raw):
    """Malformed versions raise ValueError."""
    with pytest.raises(ValueError):
        dr.parse_zenodo_version(raw)


def test_latest_is_default(fake_zenodo):
    """Without a version, the newest record is used."""
    (rec,) = dr.resolve_zenodo_records("GWTC-3")
    assert (rec.record_id, rec.version, rec.n_versions) == ("22685054", 3, 3)
    assert rec.label == "v3 (latest), record 22685054"


def test_explicit_older_version(fake_zenodo):
    """vN counts from the oldest version."""
    (rec,) = dr.resolve_zenodo_records("GWTC-3", "v1")
    assert rec.record_id == "5546663"


def test_unknown_version_lists_available(fake_zenodo):
    """An out-of-range version names the available ones."""
    with pytest.raises(ValueError, match=r"v1 \(record 5546663.*v3 \(record 22685054"):
        dr.resolve_zenodo_records("GWTC-3", "v4")


def test_unknown_catalog():
    """Catalogs without a Zenodo release are rejected."""
    with pytest.raises(ValueError, match="No Zenodo release configured"):
        dr.resolve_zenodo_records("GWTC-1")


@pytest.mark.parametrize("version, filename", [
    (None, "IGWN-GWTC3p0-v3-PESkyLocalizations.tar.gz"),
    ("v2", "IGWN-GWTC3p0-v2-PESkyLocalizations.tar.gz"),
    ("v1", "skymaps.tar.gz"),  # GWTC-3 v1 used a plain name
])
def test_skymap_tarball_found_by_name(fake_zenodo, version, filename):
    """The skymap tarball is picked by name pattern, never the contour tarball."""
    rec, fname = dr.zenodo_skymap_tarball("GWTC-3", version)
    assert fname == filename
    assert dr.zenodo_skymap_url("GWTC-3", version) == (
        f"https://zenodo.org/records/{rec.record_id}/files/{filename}?download=1"
    )


def test_listing_is_cached(fake_zenodo):
    """A fresh cached listing avoids a second API call."""
    dr.resolve_zenodo_records("GWTC-3")
    dr.resolve_zenodo_records("GWTC-3", "v2")
    assert len(fake_zenodo) == 1


def test_stale_cache_used_when_offline(monkeypatch, tmp_path):
    """If Zenodo is unreachable, a stale cached listing still resolves versions."""
    cache = tmp_path / "versions_22685054.json"
    cache.write_text(json.dumps({"fetched": 0, "versions": FAKE_GWTC3}))
    monkeypatch.setattr(dr, "zenodo_cache_dir", lambda: tmp_path)

    def offline(record_id):
        raise requests.ConnectionError("offline")

    monkeypatch.setattr(dr, "_fetch_versions", offline)
    (rec,) = dr.resolve_zenodo_records("GWTC-3", "v2")
    assert rec.record_id == "8177023"


def test_offline_without_cache_falls_back_to_configured_record(monkeypatch, tmp_path):
    """Offline and uncached: latest falls back to the configured record; vN fails."""
    monkeypatch.setattr(dr, "zenodo_cache_dir", lambda: tmp_path)

    def offline(record_id):
        raise requests.ConnectionError("offline")

    monkeypatch.setattr(dr, "_fetch_versions", offline)
    rec, fname = dr.zenodo_skymap_tarball("GWTC-3")
    assert (rec.record_id, rec.version) == ("22685054", None)
    assert fname == "IGWN-GWTC3p0-v3-PESkyLocalizations.tar.gz"
    with pytest.raises(ValueError, match="Cannot list the Zenodo versions"):
        dr.resolve_zenodo_records("GWTC-3", "v2")


def test_cli_zenodo_version_parsing():
    """--zenodo-version takes CATALOG=VERSION items and needs --data-repo zenodo."""
    args = build_parser().parse_args([
        "search_skymaps", "--catalogs", "GWTC-3", "--ra-deg", "1", "--dec-deg", "2",
        "--zenodo-version", "GWTC-3=v2", "GWTC-4=latest",
    ])
    assert _parse_zenodo_versions(args.zenodo_version, args.data_repo) == {"GWTC-3": "v2", "GWTC-4": "latest"}
    assert _parse_zenodo_versions(None, "zenodo") is None
    with pytest.raises(ValueError, match="only applies with --data-repo zenodo"):
        _parse_zenodo_versions(["GWTC-3=v2"], "s3")
    with pytest.raises(ValueError, match="expected CATALOG=VERSION"):
        _parse_zenodo_versions(["v2"], "zenodo")
    with pytest.raises(ValueError, match="No Zenodo release for catalog"):
        _parse_zenodo_versions(["GWTC-1=v1"], "zenodo")

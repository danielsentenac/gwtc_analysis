"""Offline tests for the unofficial PE bundle sources and cache.

Downloads are served by a fake ``requests.get``; the bundle writer is replaced
by a stub, so no network access or PESummary build is needed.
"""
from __future__ import annotations

import io
import tarfile
from dataclasses import replace

import pytest

from gwtc_analysis import unofficial_pe as up


def _tar_gz(members: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
    return buf.getvalue()


class _Resp:
    def __init__(self, content: bytes):
        self.content = content

    def raise_for_status(self):
        pass


@pytest.fixture
def spec(tmp_path):
    """A GW170817-like spec whose public sources live under tmp_path."""
    cal = tmp_path / "CalEnv"
    return replace(
        up.GW170817_SPEC,
        raw_samples_path=tmp_path / "samples.hdf5",
        psd_path=tmp_path / "psds.dat",
        skymap_path=tmp_path / "skymap.fits.gz",
        calibration_paths=(("H1", cal / "H.txt"), ("L1", cal / "L.txt")),
        downloads=(
            up.PublicSource("https://dcc.example/samples.hdf5", tmp_path / "samples.hdf5"),
            up.PublicSource("https://dcc.example/psds.dat", tmp_path / "psds.dat"),
            up.PublicSource("https://dcc.example/skymap.fits.gz", tmp_path / "skymap.fits.gz"),
            up.PublicSource("https://dcc.example/cal.tar.gz", cal / "H.txt", member="CalEnv/H.txt"),
            up.PublicSource("https://dcc.example/cal.tar.gz", cal / "L.txt", member="CalEnv/L.txt"),
        ),
    )


@pytest.fixture
def fake_dcc(monkeypatch):
    """Serve fixed payloads for the fake DCC urls and record the requests."""
    payloads = {
        "https://dcc.example/samples.hdf5": b"samples",
        "https://dcc.example/psds.dat": b"psds",
        "https://dcc.example/skymap.fits.gz": b"skymap",
        "https://dcc.example/cal.tar.gz": _tar_gz({"./CalEnv/H.txt": b"H", "./CalEnv/L.txt": b"L"}),
    }
    calls = []

    def get(url, timeout=None):
        calls.append(url)
        return _Resp(payloads[url])

    import requests
    monkeypatch.setattr(requests, "get", get)
    return calls


def test_gw170817_spec_uses_public_dcc_sources_only():
    """Every GW170817 input is downloaded from a public DCC release."""
    spec = up.GW170817_SPEC
    assert spec.fit_extrinsic
    assert not spec.asd_paths
    assert all(a.lalinference_samples_path is None for a in spec.analyses)
    assert all(s.url.startswith("https://dcc.ligo.org/public/") for s in spec.downloads)
    downloaded = {s.path for s in spec.downloads}
    needed = {spec.raw_samples_path, spec.psd_path, spec.skymap_path, *(p for _, p in spec.calibration_paths)}
    assert needed <= downloaded


def test_download_fetches_missing_files_and_extracts_tar_members(spec, fake_dcc):
    """Missing sources are downloaded; a shared archive is fetched once."""
    up._download_public_sources(spec, log_cb=None)
    assert spec.raw_samples_path.read_bytes() == b"samples"
    assert dict(spec.calibration_paths)["H1"].read_bytes() == b"H"
    assert dict(spec.calibration_paths)["L1"].read_bytes() == b"L"
    assert fake_dcc.count("https://dcc.example/cal.tar.gz") == 1

    fake_dcc.clear()
    up._download_public_sources(spec, log_cb=None)
    assert fake_dcc == []  # nothing left to download


def test_failed_download_is_reported_not_raised(spec, monkeypatch):
    """A download error is logged; the build then reports the source as missing."""
    import requests

    def offline(url, timeout=None):
        raise requests.ConnectionError("offline")

    monkeypatch.setattr(requests, "get", offline)
    logs = []
    up._download_public_sources(spec, log_cb=logs.append)
    assert any("Could not download" in m for m in logs)
    assert not spec.raw_samples_path.exists()


def test_bundle_cache_is_rebuilt_when_the_recipe_changes(spec, fake_dcc, monkeypatch, tmp_path):
    """A cached bundle is reused only if it was built by the same recipe."""
    builds = []

    def fake_write(s, out_path, *, log_cb=None):
        builds.append(s)
        out_path.write_bytes(b"bundle")

    monkeypatch.setattr(up, "_write_unofficial_pesummary_bundle", fake_write)
    monkeypatch.setitem(up.UNOFFICIAL_PE_BUNDLES, "GW170817", spec)
    cache = tmp_path / "cache"

    target = up.build_unofficial_pe_bundle("GW170817", cache_dir=cache)
    assert target is not None and len(builds) == 1
    up.build_unofficial_pe_bundle("GW170817", cache_dir=cache)
    assert len(builds) == 1  # same recipe: cached bundle reused

    monkeypatch.setitem(up.UNOFFICIAL_PE_BUNDLES, "GW170817", replace(spec, fit_extrinsic=False))
    up.build_unofficial_pe_bundle("GW170817", cache_dir=cache)
    assert len(builds) == 2  # recipe changed: rebuilt

    target.with_name(target.name + ".recipe.json").unlink()
    up.build_unofficial_pe_bundle("GW170817", cache_dir=cache)
    assert len(builds) == 3  # bundle without a recipe file (older builder): rebuilt

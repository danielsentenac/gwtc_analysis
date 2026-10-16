"""Public files: chosen by their URL, downloaded from the S3 mirror with --data-repo s3 when it holds them."""
from pathlib import Path

import pytest

from gwtc_analysis import mirror, net, report
from gwtc_analysis.cli import build_parser


def test_mirror_keys():
    assert mirror.mirror_key("https://zenodo.org/records/19500052/files/mixture-real_o3.hdf?download=1") == \
        "public/zenodo/19500052/mixture-real_o3.hdf"
    assert mirror.mirror_key("https://dcc.ligo.org/public/0157/P1800370/005/GW170817_GWTC-1.hdf5") == \
        "public/dcc/P1800370/005/GW170817_GWTC-1.hdf5"
    assert mirror.mirror_key("https://arxiv.org/e-print/2501.01495") == "public/arxiv/2501.01495.tar.gz"
    assert mirror.mirror_key("https://gwosc.org/eventapi/jsonfull/GWTC/") is None


@pytest.fixture
def clean(monkeypatch):
    monkeypatch.setattr(report, "_DATA_FILES", {})
    monkeypatch.setattr(report, "_PE_SOURCES", {})
    yield
    mirror.set_data_repo("zenodo")


URL = "https://zenodo.org/records/16877102/files/GW250114_data_release.tar.gz?download=1"


def test_fetch_from_mirror(tmp_path, monkeypatch, clean):
    mirror.set_data_repo("s3")
    keys = []
    monkeypatch.setattr(mirror, "_from_s3", lambda url, dest, log: (keys.append(mirror.mirror_key(url)),
                                                                     Path(dest).write_bytes(b"s3"))[0] is None)
    out = mirror.fetch_to(URL, tmp_path / "a.tgz", lambda u, p: pytest.fail("not from the public source"))
    assert out.read_bytes() == b"s3" and keys == ["public/zenodo/16877102/GW250114_data_release.tar.gz"]
    assert report._DATA_FILES == {URL.removesuffix("?download=1"): "S3 mirror"}
    assert "Public data files of this run" in report.data_files_html()


def test_fetch_falls_back_to_the_source(tmp_path, monkeypatch, clean):
    mirror.set_data_repo("s3")
    monkeypatch.setattr(mirror, "_from_s3", lambda url, dest, log: False)
    out = mirror.fetch_to(URL, tmp_path / "a.tgz", lambda u, p: Path(p).write_bytes(b"zenodo"))
    assert out.read_bytes() == b"zenodo" and not (tmp_path / "a.tgz.part").exists()
    assert report._DATA_FILES[URL.removesuffix("?download=1")] == "zenodo.org (not on the S3 mirror)"
    assert mirror.fetch_bytes(URL, lambda u: b"zenodo") == b"zenodo"


def test_zenodo_repo_never_asks_the_mirror(tmp_path, monkeypatch, clean):
    mirror.set_data_repo("zenodo")
    monkeypatch.setattr("gwtc_analysis.gw_stat._s3_client_from_env", lambda: pytest.fail("no S3 with zenodo"))
    mirror.fetch_to(URL, tmp_path / "a.tgz", lambda u, p: Path(p).write_bytes(b"z"))
    assert report._DATA_FILES[URL.removesuffix("?download=1")] == "zenodo.org"
    assert not mirror.on_mirror(URL)


def test_every_downloading_mode_has_data_repo():
    modes = build_parser()._mode_parsers
    for m in ("rates", "hubble_constant", "counterpart", "continuous_waves", "area_law", "scalar_cloud",
              "neutron_star_eos", "build_unofficial_pe", "spin_population", "lens_pairs", "parameters_estimation"):
        assert any("--data-repo" in a.option_strings for a in modes[m]._actions), m


def test_failed_s3_transfer_falls_back(tmp_path, monkeypatch, clean):
    """A mirror that holds the file but fails the transfer: the public source."""
    mirror.set_data_repo("s3")

    class Client:
        def stat_object(self, bucket, key):
            return object()

    monkeypatch.setattr("gwtc_analysis.gw_stat._s3_client_from_env", lambda: Client())

    def broken(*a, **k):
        raise OSError("connection reset")

    monkeypatch.setattr("gwtc_analysis.parameters_estimation._download_s3_with_progress", broken)
    out = mirror.fetch_to(URL, tmp_path / "a.tgz", lambda u, p: Path(p).write_bytes(b"zenodo"), log=lambda m: None)
    assert out.read_bytes() == b"zenodo"


def test_unreachable_mirror_pe_file_from_zenodo(tmp_path, monkeypatch, clean):
    from gwtc_analysis import parameters_estimation as pe
    from gwtc_analysis import pe_files

    monkeypatch.setattr(pe, "_shared_pe_file", lambda f: None)

    def unreachable(cache):
        raise ConnectionError("minio-dev.odahub.fr unreachable")

    monkeypatch.setattr(pe_files, "s3_objects", unreachable)
    monkeypatch.setattr(pe, "_download_http_with_progress", lambda url, out, desc=None: Path(out).write_bytes(b"z"))
    entry = {"filename": "IGWN-x-GW150914_095045-combined_PEDataRelease.hdf5", "url": "https://zenodo.org/x"}
    path, downloaded = pe_files.fetch(entry, tmp_path, "s3", log=lambda m: None)
    assert path.read_bytes() == b"z" and downloaded
    assert report._PE_SOURCES[entry["filename"]] == "Zenodo (S3 mirror failed)"


def test_s3_skymaps_of_the_selected_release_only(monkeypatch):
    from gwtc_analysis import data_repo, gw_stat

    monkeypatch.setattr(data_repo, "zenodo_skymap_tarball",
                        lambda cat, version=None: (None, "IGWN-GWTC4p0-38214bd95_724-Archived_Skymaps.tar.gz"))
    tag = data_repo.skymap_release_tag("GWTC-4")
    assert tag == "IGWN-GWTC4p0-38214bd95_724"

    class Obj:
        def __init__(self, n):
            self.object_name = n

    class Client:
        def list_objects(self, bucket, prefix, recursive):
            return [Obj(f"GWTC-4/IGWN-GWTC4p0-1a206db3d_721-GW230529_181500-Mixed_Skymap_PEDataRelease.fits.gz"),
                    Obj(f"GWTC-4/{tag}-GW230529_181500-Mixed_Skymap_PEDataRelease.fits.gz"),
                    Obj(f"GWTC-4/IGWN-GWTC4p0-1a206db3d_721-GW230601_224134-Mixed_Skymap_PEDataRelease.fits.gz")]

    idx = gw_stat._build_s3_skymap_index(Client(), "gwtc", "GWTC-4/", release_tag=tag)
    assert idx == {"GW230529_181500": f"GWTC-4/{tag}-GW230529_181500-Mixed_Skymap_PEDataRelease.fits.gz"}


def test_bright_siren_takes_data_repo():
    from gwtc_analysis.cli import _apply_settings_file, _resolve_h0_method

    p = build_parser()
    a = _apply_settings_file(p, ["hubble_constant", "--method", "bright", "--data-repo", "galaxy"])
    _resolve_h0_method(a, p._mode_parsers["hubble_constant"])
    assert a.data_repo == "galaxy"


def test_network_error_is_retried(tmp_path, monkeypatch, clean):
    """A DNS failure or a reset (tools-iuc CI, dcc.ligo.org): retried, then the file is written whole."""
    import requests

    waits, tries = [], []
    monkeypatch.setattr(net.time, "sleep", waits.append)

    def flaky(u, p):
        tries.append(u)
        if len(tries) < 3:
            raise requests.ConnectionError("Failed to resolve 'dcc.ligo.org'")
        Path(p).write_bytes(b"ok")

    out = mirror.fetch_to(URL, tmp_path / "a.tgz", flaky, log=lambda m: None)
    assert out.read_bytes() == b"ok" and len(tries) == 3 and waits == [5, 15]
    tries.clear()
    assert mirror.fetch_bytes(URL, lambda u: flaky(u, tmp_path / "b") or b"bytes", log=lambda m: None) == b"bytes"


def test_persistent_network_error_and_http_error_are_raised(tmp_path, monkeypatch, clean):
    import requests

    monkeypatch.setattr(net.time, "sleep", lambda s: None)
    calls = []

    def down(u, p):
        calls.append(u)
        raise requests.Timeout("read timed out")

    with pytest.raises(requests.Timeout):
        mirror.fetch_to(URL, tmp_path / "a.tgz", down, log=lambda m: None)
    assert len(calls) == 3 and not (tmp_path / "a.tgz").exists()

    calls.clear()

    def missing(u, p):
        calls.append(u)
        raise requests.HTTPError("404 Client Error")

    with pytest.raises(requests.HTTPError):
        mirror.fetch_to(URL, tmp_path / "a.tgz", missing, log=lambda m: None)
    assert len(calls) == 1


def test_server_overload_is_retried(tmp_path, monkeypatch, clean):
    """A 503 raised by the download function (raise_for_status) is transient: retried, unlike a 404."""
    import requests

    monkeypatch.setattr(net.time, "sleep", lambda s: None)
    tries = []

    def overloaded(u, p):
        tries.append(u)
        if len(tries) == 1:
            r = requests.Response()
            r.status_code = 503
            raise requests.HTTPError("503 Server Error: Service Unavailable", response=r)
        Path(p).write_bytes(b"ok")

    assert mirror.fetch_to(URL, tmp_path / "a.tgz", overloaded, log=lambda m: None).read_bytes() == b"ok"
    assert len(tries) == 2


@pytest.fixture
def backup(monkeypatch):
    """The Galaxy history as a fake mirror holding every file; records the repositories asked."""
    asked = []

    def from_repo(repo, url, dest, log, *, otherwise):
        asked.append(repo)
        Path(dest).write_bytes(b"galaxy copy")
        return True

    monkeypatch.setattr(mirror, "_from_repo", from_repo)
    monkeypatch.setattr(net.time, "sleep", lambda s: None)
    return asked


def _down(u, p):
    import requests

    raise requests.ConnectionError("Failed to resolve 'dcc.ligo.org'")


DCC_URL = "https://dcc.ligo.org/public/0201/P2500416/006/uls_soap_noheader.txt"


def test_failing_source_falls_back_on_the_galaxy_copy(tmp_path, backup, clean):
    """--data-repo zenodo (the default): the source is retried, then the Galaxy copy is read."""
    mirror.set_data_repo("zenodo")
    logs = []
    out = mirror.fetch_to(DCC_URL, tmp_path / "a.txt", _down, log=logs.append)
    assert out.read_bytes() == b"galaxy copy" and backup == ["galaxy"]
    assert report._DATA_FILES[DCC_URL] == "Galaxy history (dcc.ligo.org failed)"
    assert any("trying the copy on the Galaxy history" in m for m in logs)
    assert mirror.fetch_bytes(DCC_URL, lambda u: _down(u, None), log=lambda m: None) == b"galaxy copy"


def test_no_fallback_for_a_missing_file_or_a_run_already_on_galaxy(tmp_path, monkeypatch, backup, clean):
    import requests

    mirror.set_data_repo("zenodo")

    def missing(u, p):
        r = requests.Response()
        r.status_code = 404
        raise requests.HTTPError("404 Client Error", response=r)

    with pytest.raises(requests.HTTPError):
        mirror.fetch_to(DCC_URL, tmp_path / "a.txt", missing, log=lambda m: None)
    assert backup == []

    mirror.set_data_repo("galaxy")        # the history was asked first (and lacks the file): not again
    asked = []
    monkeypatch.setattr(mirror, "_from_repo", lambda repo, url, dest, log, *, otherwise: asked.append(repo) and False)
    with pytest.raises(requests.ConnectionError):
        mirror.fetch_to(DCC_URL, tmp_path / "b.txt", _down, log=lambda m: None)
    assert asked == ["galaxy"] and not (tmp_path / "b.txt").exists()

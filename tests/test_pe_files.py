"""The PE file read depends only on the options: the Zenodo listing and the rule of the mode, not --data-repo."""
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from gwtc_analysis import hubble_constant as hc
from gwtc_analysis import parameters_estimation as pe
from gwtc_analysis import pe_files

COSMO = "IGWN-GWTC3p0-v3-GW191103_012549_PEDataRelease_mixed_cosmo.h5"
NOCOSMO = "IGWN-GWTC3p0-v3-GW191103_012549_PEDataRelease_mixed_nocosmo.h5"
COMBINED = "IGWN-GWTC4p0-38214bd95_724-GW230529_181500-combined_PEDataRelease.hdf5"


def _c(*names):
    return [{"filename": n, "url": f"https://zenodo.org/records/1/files/{n}?download=1"} for n in names]


def test_rules():
    assert pe_files.pick(_c(COSMO, NOCOSMO), pe_files.SAMPLING)["filename"] == NOCOSMO
    assert pe_files.pick(_c(NOCOSMO, COSMO), pe_files.CATALOG)["filename"] == COSMO
    assert pe_files.pick(_c(COMBINED), pe_files.SAMPLING)["filename"] == COMBINED
    assert pe_files.pick(_c(COMBINED), pe_files.CATALOG)["filename"] == COMBINED
    assert hc._pick_pe_file(_c(COSMO, NOCOSMO))["filename"] == NOCOSMO
    assert pe.choose_best_pe_file(_c(NOCOSMO, COSMO))["filename"] == COSMO
    assert pe_files.pick([], pe_files.SAMPLING) is None


def test_cosmo_reweighted_prior():
    recorded = "PowerLaw(alpha=2, minimum=100, maximum=10000, name='luminosity_distance')"
    assert pe_files.is_cosmo_reweighted(COSMO) and not pe_files.is_cosmo_reweighted(NOCOSMO)
    assert not pe_files.is_cosmo_reweighted(COMBINED)
    assert pe_files.distance_prior_desc(NOCOSMO, recorded) == recorded
    desc = pe_files.distance_prior_desc(COSMO, recorded)
    dl = np.array([200.0, 1000.0, 3000.0])
    p, _ = hc.pe_distance_prior(dl, desc, "O3b")
    q, _ = hc.pe_distance_prior(dl, "", "O3b", kind="source-frame")
    # uniform in comoving volume and source-frame time (astropy Planck15, close to Planck15_LAL), not D_L^2
    assert np.allclose(p / p[0], q / q[0], rtol=0.02)
    assert not np.allclose(p / p[0], (dl / dl[0]) ** 2, rtol=0.1)


def _no_shared(monkeypatch):
    monkeypatch.setattr(pe, "_shared_pe_file", lambda fname: None)


def test_fetch_s3_has_the_file(tmp_path, monkeypatch):
    _no_shared(monkeypatch)
    got = []
    monkeypatch.setattr(pe_files, "s3_objects", lambda cache: {NOCOSMO: f"GWTC-3/PE/{NOCOSMO}"})
    monkeypatch.setattr("gwtc_analysis.gw_stat._s3_client_from_env", lambda: object())
    monkeypatch.setattr(pe, "_download_s3_with_progress",
                        lambda client, bucket, object_name, out_path, desc: (got.append(object_name),
                                                                             Path(out_path).write_bytes(b"x")))
    monkeypatch.setattr(pe, "_download_http_with_progress", lambda *a, **k: (_ for _ in ()).throw(AssertionError))
    path, downloaded = pe_files.fetch(_c(NOCOSMO)[0], tmp_path, "s3")
    assert got == [f"GWTC-3/PE/{NOCOSMO}"] and downloaded and path == tmp_path / "files" / NOCOSMO


def test_fetch_s3_without_the_file_uses_zenodo(tmp_path, monkeypatch):
    _no_shared(monkeypatch)
    urls = []
    # the mirror holds only another version and another kind: neither is read
    monkeypatch.setattr(pe_files, "s3_objects", lambda cache: {COSMO.replace("v3", "v2"): "GWTC-3/PE/old.h5",
                                                               COSMO: f"GWTC-3/PE/{COSMO}"})
    monkeypatch.setattr(pe, "_download_s3_with_progress", lambda *a, **k: (_ for _ in ()).throw(AssertionError))
    monkeypatch.setattr(pe, "_download_http_with_progress",
                        lambda url, out, desc=None: (urls.append(url), Path(out).write_bytes(b"x")))
    path, downloaded = pe_files.fetch(_c(NOCOSMO)[0], tmp_path, "s3", log=lambda m: None)
    assert urls == [_c(NOCOSMO)[0]["url"]] and downloaded and path.name == NOCOSMO
    # already there: not downloaded again, and not to be deleted by the caller
    assert pe_files.fetch(_c(NOCOSMO)[0], tmp_path, "s3") == (path, False)


def test_extract_is_current(tmp_path):
    p = tmp_path / "x.h5"
    with h5py.File(p, "w") as h:
        h.attrs["source"] = NOCOSMO
    assert pe_files.extract_is_current(p, NOCOSMO)
    assert not pe_files.extract_is_current(p, COSMO)
    with h5py.File(p, "w") as h:
        h.attrs["source_file"] = COMBINED          # hubble_constant extracts
    assert pe_files.extract_is_current(p, COMBINED)
    with h5py.File(p, "w"):
        pass
    assert not pe_files.extract_is_current(p, NOCOSMO)
    assert not pe_files.extract_is_current(tmp_path / "missing.h5", NOCOSMO)


def test_hubble_extract_from_another_file_is_redone(tmp_path, monkeypatch):
    """An extract of an earlier release version (or of the cosmo file) is made again from the selected file."""
    name = "GW191103_012549"
    (tmp_path / "samples").mkdir()
    with h5py.File(tmp_path / "samples" / f"{name}.h5", "w") as h:
        h.attrs["source_file"] = COSMO.replace("v3", "v2")
    monkeypatch.setattr(pe, "build_zenodo_pe_index", lambda **k: {name: _c(COSMO, NOCOSMO)})
    fetched = []

    def fake_fetch(entry, cache, data_repo="zenodo", desc="", log=print):
        fetched.append(entry["filename"])
        return Path(cache) / "files" / entry["filename"], False

    monkeypatch.setattr(pe_files, "fetch", fake_fetch)
    monkeypatch.setattr(hc, "_extract_pe", lambda src, out: h5py.File(out, "w").attrs.update(source_file=src.name))
    ev = pd.DataFrame([dict(event=name, common_name="", gps=1256779567.0)])
    out = hc.fetch_pe_samples(ev, tmp_path)
    assert fetched == [NOCOSMO] and pe_files.extract_is_current(out[name], NOCOSMO)
    assert hc.fetch_pe_samples(ev, tmp_path) == out and fetched == [NOCOSMO]      # current: kept


def test_hubble_extract_of_a_cosmo_file_records_its_prior(tmp_path):
    src = tmp_path / COSMO
    with h5py.File(src, "w") as f:
        g = f.create_group("C01:IMRPhenomXPHM")
        g.create_dataset("posterior_samples", data=np.zeros(3, dtype=[(c, "f8") for c in hc._PE_COLUMNS]))
        g.create_group("priors/analytic").create_dataset("luminosity_distance", data=np.array([b"PowerLaw(alpha=2)"]))
    out = tmp_path / "x.h5"
    hc._extract_pe(src, out)
    with h5py.File(out, "r") as h:
        assert h["C01:IMRPhenomXPHM"].attrs["prior:luminosity_distance"] == pe_files.COSMO_REWEIGHTED_PRIOR


def test_staged_galaxy_files_follow_the_rule(tmp_path):
    d = tmp_path / "GWTC-3-PE"
    d.mkdir()
    for n in (NOCOSMO, COSMO):
        (d / n).write_bytes(b"x")
    assert pe._find_local_pe_file("GW191103_012549", [d]).name == COSMO


def test_report_lists_the_pe_files(monkeypatch):
    from gwtc_analysis import report

    monkeypatch.setattr(report, "_PE_FILES", {})
    monkeypatch.setattr(report, "_PE_SOURCES", {})
    monkeypatch.setattr(report, "_DATA_FILES", {})
    monkeypatch.setattr(report, "_RUN_OPTIONS", None)
    assert report.options_html() == ""
    report.record_pe_file("GW191103_012549", NOCOSMO)
    report.note_pe_source(NOCOSMO, "Zenodo (not on the S3 mirror)")
    report.record_pe_file("GW230529_181500", f"/cache/files/{COMBINED}")
    report.record_pe_file("GW170817", "Unofficial-GWTC1-GW170817_PEDataRelease.h5", "build_unofficial_pe")
    html = report.options_html()
    assert "PE files of this run" in html and "3 event(s)" in html
    assert "1 nocosmo" in html and "1 combined" in html and "1 unofficial bundle" in html
    assert "Zenodo (not on the S3 mirror)" in html and "build_unofficial_pe" in html
    assert COMBINED in html and "/cache/files" not in html
    assert "samples already extracted from it" in html          # no source noted: an extract of the PE cache
    assert pe_files.describe(COSMO)[0] == "cosmo" and pe_files.describe("my.h5")[0] == "other"

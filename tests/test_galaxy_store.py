"""--data-repo galaxy works as s3: the public Galaxy history seen as the bucket, Zenodo as fallback."""
from pathlib import Path

import pytest

from gwtc_analysis import galaxy_store as gs
from gwtc_analysis import mirror, pe_files, report

PE3 = "IGWN-GWTC3p0-v3-GW191103_012549_PEDataRelease_mixed_nocosmo.h5"
PE4 = "IGWN-GWTC4p0-38214bd95_724-GW230529_181500-combined_PEDataRelease.hdf5"
MAP41 = "IGWN-GWTC4p1-18965dda8_5-GW230529_181500-Mixed_Skymap_PEDataRelease.fits.gz"
PUB = "public/zenodo/16877102/GW250114_data_release.tar.gz"
FILES = {n: {"id": f"id{i}", "ext": "h5", "size": 10 + i} for i, n in enumerate((PE3, PE4, MAP41, PUB, "notes.txt"))}


@pytest.fixture
def history(monkeypatch):
    monkeypatch.setattr(gs, "history_files", lambda refresh=False: FILES)
    monkeypatch.setattr(report, "_PE_SOURCES", {})
    monkeypatch.setattr(report, "_DATA_FILES", {})
    yield
    mirror.set_data_repo("zenodo")


def test_object_names():
    assert gs.object_name(PE3) == f"GWTC-3/PE/{PE3}"
    assert gs.object_name(PE4) == f"GWTC-4/PE/{PE4}"
    assert gs.object_name(MAP41) == f"GWTC-4.1/{MAP41}"          # not inside GWTC-4/
    assert gs.object_name(PUB) == PUB
    assert gs.object_name("notes.txt") is None


def test_store_calls(history):
    st = gs.GalaxyHistoryStore()
    assert [o.object_name for o in st.list_objects("gwtc", prefix="GWTC-4/", recursive=True)] == [f"GWTC-4/PE/{PE4}"]
    assert st.stat_object("gwtc", PUB).size == FILES[PUB]["size"]
    with pytest.raises(KeyError):
        st.stat_object("gwtc", "public/zenodo/1/missing")
    assert pe_files.s3_objects("unused", "galaxy") == {PE3: f"GWTC-3/PE/{PE3}", PE4: f"GWTC-4/PE/{PE4}"}


def test_pe_file_from_the_history_else_zenodo(tmp_path, history, monkeypatch):
    from gwtc_analysis import parameters_estimation as pe

    monkeypatch.setattr(pe, "_shared_pe_file", lambda f: None)
    got = []
    monkeypatch.setattr(pe, "_download_s3_with_progress",
                        lambda client, bucket, object_name, out_path, desc: (got.append(object_name),
                                                                             Path(out_path).write_bytes(b"g")))
    monkeypatch.setattr(pe, "_download_http_with_progress", lambda url, out, desc=None: Path(out).write_bytes(b"z"))
    path, _ = pe_files.fetch({"filename": PE3, "url": "https://zenodo.org/x"}, tmp_path, "galaxy")
    assert path.read_bytes() == b"g" and got == [f"GWTC-3/PE/{PE3}"]
    assert report._PE_SOURCES[PE3] == "Galaxy history"
    other = "IGWN-GWTC3p0-v3-GW191105_143521_PEDataRelease_mixed_nocosmo.h5"
    path, _ = pe_files.fetch({"filename": other, "url": "https://zenodo.org/y"}, tmp_path, "galaxy",
                             log=lambda m: None)
    assert path.read_bytes() == b"z" and report._PE_SOURCES[other] == "Zenodo (not on the Galaxy history)"


def test_public_file_from_the_history(tmp_path, history, monkeypatch):
    mirror.set_data_repo("galaxy")
    monkeypatch.setattr("gwtc_analysis.parameters_estimation._download_s3_with_progress",
                        lambda client, bucket, object_name, out_path, desc: Path(out_path).write_bytes(b"g"))
    url = "https://zenodo.org/records/16877102/files/GW250114_data_release.tar.gz?download=1"
    out = mirror.fetch_to(url, tmp_path / "a.tgz", lambda u, p: pytest.fail("not from Zenodo"), log=lambda m: None)
    assert out.read_bytes() == b"g" and report._DATA_FILES[url.removesuffix("?download=1")] == "Galaxy history"

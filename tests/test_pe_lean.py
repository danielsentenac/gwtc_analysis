"""The PE file read one analysis at a time: what the copies hold, and the samples taken from the file's table."""
import h5py
import numpy as np
import pytest

from gwtc_analysis import pe_lean


@pytest.fixture
def pe_file(tmp_path):
    f = tmp_path / "pe.h5"
    rng = np.random.default_rng(0)
    with h5py.File(f, "w") as h:
        h.create_group("version")["pesummary"] = [b"1.0"]
        for lab, n, psd in (("C01:A", 40, True), ("C01:B", 30, False), ("C01:Mixed", 70, None)):
            g = h.create_group(lab)
            t = np.zeros(n, dtype=[("mass_1", float), ("mass_2", float), ("log_likelihood", float)])
            t["mass_1"], t["mass_2"], t["log_likelihood"] = rng.normal(30, 2, n), rng.normal(20, 2, n), rng.normal(0, 1, n)
            g["posterior_samples"] = t
            g.create_group("skymap")["data"] = np.ones(48)
            g.create_group("priors").create_group("samples")["mass_1"] = np.ones(1000)
            g.create_group("config_file")["x"] = [b"1"]
            if psd is not None:
                g.create_group("psds")["H1"] = np.ones((50, 2)) if psd else np.ones((0, 2))
    return f


def test_first_copy_drops_what_fills_the_memory(pe_file, tmp_path):
    with h5py.File(pe_file) as src, h5py.File(tmp_path / "stub.h5", "w") as dst:
        pe_lean._copy(src, dst)
        assert set(dst) == set(src)
        for lab in ("C01:A", "C01:B", "C01:Mixed"):
            assert len(dst[lab]["posterior_samples"]) == pe_lean.STUB_ROWS
            assert set(dst[lab]) == {"posterior_samples", "config_file"}


def test_copy_of_one_label(pe_file, tmp_path):
    with h5py.File(pe_file) as src, h5py.File(tmp_path / "one.h5", "w") as dst:
        pe_lean._copy(src, dst, label="C01:A", psds=True)
        assert set(dst) == {"version", "C01:A"}
        assert set(dst["C01:A"]) == {"posterior_samples", "config_file", "psds"}
        assert dst["C01:A"]["psds"]["H1"].shape == (50, 2)


def test_samples_are_the_table_of_the_file(pe_file, tmp_path, monkeypatch):
    pytest.importorskip("pesummary")
    monkeypatch.setattr(pe_lean, "KEEP", 2)
    s = pe_lean._LazySamples(pe_lean._OneLabel(pe_file, tmp_path), ["C01:A", "C01:B", "C01:Mixed"])
    assert sorted(s.keys()) == ["C01:A", "C01:B", "C01:Mixed"] and "C01:A" in s and len(s) == 3
    with h5py.File(pe_file) as h:
        for lab in s:
            t = h[lab]["posterior_samples"][()]
            assert list(s[lab].keys()) == list(t.dtype.names)
            assert np.array_equal(np.asarray(s[lab]["mass_1"]), t["mass_1"])
            assert np.array_equal(s.light(lab)["mass_2"], t["mass_2"])
    assert list(s._loaded) == ["C01:B", "C01:Mixed"]          # the first label read was dropped
    with pytest.raises(KeyError):
        s["C01:none"]


def test_table_pesummary_would_edit_is_not_taken_directly(pe_file, tmp_path):
    with h5py.File(pe_file, "a") as h:
        t = h["C01:A"]["posterior_samples"][()]
        t["log_likelihood"][3] = np.nan
        del h["C01:A"]["posterior_samples"]
        h["C01:A"]["posterior_samples"] = t
    s = pe_lean._LazySamples(pe_lean._OneLabel(pe_file, tmp_path), ["C01:A"])
    assert s._from_table("C01:A") is None


def test_psd_availability_without_reading(pe_file, tmp_path):
    from gwtc_analysis.gwpe_utils import label_has_psd

    class Reader:
        def read(self, label, **what):
            raise AssertionError("the spectra must not be read to know whether they exist")

    class Data:
        psd = pe_lean._LazyPsds(Reader(), {"C01:A": True, "C01:B": False})

    assert label_has_psd(Data, "C01:A") and not label_has_psd(Data, "C01:B") and not label_has_psd(Data, "C01:Mixed")
    Data.psd["C01:B"] = {"H1": np.ones((5, 2))}               # a supplementary spectrum attached
    assert label_has_psd(Data, "C01:B")

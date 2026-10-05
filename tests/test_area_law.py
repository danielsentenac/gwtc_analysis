"""Offline tests of the area_law mode, on a synthetic copy of the GW250114 release layout."""
from __future__ import annotations

import h5py
import numpy as np
import pytest

from gwtc_analysis import area_law as al
from gwtc_analysis.cli import build_parser


def test_kerr_area():
    """16π M² without spin, 8π M² at maximal spin; equal non-spinning masses: the remnant area grows."""
    assert al.kerr_area(1.0, 0.0) == pytest.approx(16 * np.pi)
    assert al.kerr_area(1.0, 1.0) == pytest.approx(8 * np.pi)
    m = 30.0
    assert al.kerr_area(0.952 * 2 * m, 0.686) / (2 * al.kerr_area(m, 0.0)) == pytest.approx(1.55, abs=0.02)
    assert al.MSUN_KM == pytest.approx(1.47663, abs=1e-5)


def test_significances():
    rng = np.random.default_rng(1)
    f, i = rng.normal(5, 0.3, 50_000), rng.normal(3, 0.4, 50_000)
    assert al.gaussian_significance(f, i) == pytest.approx(2 / 0.5, rel=0.01)
    assert al.pair_probability(f, i) == pytest.approx(3.2e-5, abs=4e-5)       # Φ(-4)


def _release(tmp_path, rng):
    """The files the mode reads: inspiral areas by truncation, ringdown areas, pyRing, prior."""
    d = tmp_path / "data"
    (d / "ringdown_areas").mkdir(parents=True)
    times = [-250, -40, 0]
    with h5py.File(d / "area_law_inspiral_data.hdf5", "w") as h:
        h["times"] = np.array(times)
        for k, (t, s) in enumerate(zip(times, (0.20, 0.10, 0.08))):
            h[f"area_insp_{k}"] = rng.normal(1.0, s, 6000) * 1.3e5
    for model, shift in (("220", 0.0), ("220+221", -0.2)):
        for t in (6.0, 10.5, 15.0):
            with h5py.File(d / "ringdown_areas" / f"{model}_{t:g}M_final_mass_spin_area.hdf5", "w") as h:
                h["Area_f"] = rng.normal(1.7 + shift - 0.02 * t, 0.06 + 0.005 * t, 4000) * 1.3e5
    np.save(d / "remnant_area_pyring_reweighted.npy", rng.normal(1.7, 0.12, 5000) * 1.3e5)
    x = np.linspace(-0.95, 9.5, 50)
    np.savetxt(d / "area_change_prior.dat", np.vstack([x, np.exp(-x)]))
    return d


def test_analyse_and_report(tmp_path):
    rng = np.random.default_rng(2)
    _release(tmp_path, rng)
    d = al.load(tmp_path / "data")
    assert sorted(d["inspiral"]) == [-250, -40, 0] and sorted(d["ringdown"]["220"]) == [6.0, 10.5, 15.0]
    r = al.analyse(d)
    assert r["main"] > 3 and r["min_all_cuts"] < r["main"] and r["overtone"] < r["main"]
    assert r["ratio_q"][1] == pytest.approx(0.67, abs=0.1)
    assert list(r["cut_scan"]["t_cut"]) == [-250, -40, 0] and len(r["start_scan"]) == 6
    out = al.run_area_law(cache_dir=tmp_path, out_report_html=tmp_path / "r.html", out_summary_tsv=tmp_path / "s.tsv",
                          plots_dir=tmp_path / "p")
    assert out["table"]["published"].tolist() == [4.4, 3.4, -10, 3.6]
    assert (tmp_path / "s.truncation.tsv").exists() and (tmp_path / "s.ringdown.tsv").exists()
    assert (tmp_path / "p" / "area_law_GW250114.png").exists() and (tmp_path / "p" / "area_law_scans_GW250114.png").exists()
    assert "different parts of the signal" in (tmp_path / "r.html").read_text()


def test_only_gw250114():
    with pytest.raises(ValueError, match="only GW250114"):
        al.run_area_law(src_name="GW150914")
    a = build_parser().parse_args(["area_law", "--with-imr"])
    assert a.src_name == "GW250114" and a.with_imr

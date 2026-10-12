"""continuous_waves: frequency of a boson cloud, parsing of the paper's table, the boson-cloud limits in a band."""
from __future__ import annotations

import io
import tarfile

import numpy as np
import pandas as pd
import pytest

from gwtc_analysis import continuous_waves as cw


def test_boson_frequency():
    """f = 2 m c^2 / h: 483.6 Hz for 1e-12 eV; the GW190728 band at 339-566 Hz."""
    assert cw.boson_frequency_hz(1e-12) == pytest.approx(483.6, rel=1e-3)
    lo, _, hi = cw.GW190728_BOSON_MASS_EV
    assert cw.boson_frequency_hz([lo, hi]) == pytest.approx([338.5, 565.8], rel=1e-3)


def test_latex_numbers_and_names():
    assert cw._scinum(r"\scinum{6.4(8.9)}{-27}") == pytest.approx(6.4e-27)
    assert cw._scinum(r"\scinum{1.9}{-8}") == pytest.approx(1.9e-8)
    assert cw._scinum("0.0078(0.0067)") == pytest.approx(0.0078)
    assert cw._clean_name(r"J0537\textminus6910\tablenotemark{$\beta$}") == "J0537-6910"


TABLE = r"""\startdata
\multirow{3}{*}{J0437\textminus4715\tablenotemark{$\delta$}} & \multirow{3}{*}{173.69} & \multirow{3}{*}{\scinum{-4.15}{-16}} & \multirow{3}{*}{0.16\tablenotemark{e}} & \multirow{3}{*}{\scinum{7.95}{-27}} & Bayesian & \scinum{7.2}{-27} & \scinum{8.8}{-9} & \scinum{6.8}{29} & 0.90 &  \scinum{1.1}{-26} & \scinum{3.4}{-27} & -5.4 & -10 \\
&  &  &  &  & $\mathcal{F}$-statistic & \scinum{9.0}{-27} & \scinum{1.1}{-8} & \scinum{8.5}{29} & 1.1 & x & x & 0.31 & 0.26 \\
"""


def test_known_pulsars_table(tmp_path):
    """The Bayesian row of each pulsar: distance, spin-down limit, h0, ellipticity, ratio."""
    p = tmp_path / "src.tar.gz"
    with tarfile.open(p, "w:gz") as tf:
        data = TABLE.encode()
        info = tarfile.TarInfo("Bayes-table.tex")
        info.size = len(data)
        tf.addfile(info, io.BytesIO(data))
    t = cw._known_pulsars_table(p)
    assert len(t) == 1
    r = t.iloc[0]
    assert r["pulsar"] == "J0437-4715" and r["distance_kpc"] == pytest.approx(0.16)
    assert r["h0_spindown"] == pytest.approx(7.95e-27) and r["ellipticity_ul"] == pytest.approx(8.8e-9)
    assert r["ratio_to_spindown"] == pytest.approx(0.90)


def test_boson_limits_in_band():
    """The largest excluded distance and the excluded black-hole masses, restricted to the boson-mass band."""
    m = np.array([1e-13, 5e-13, 8e-13, 1e-12, 2e-12])
    boson = {"distance": {(50.0, 1e3): pd.DataFrame(dict(m_b_ev=m, d_kpc=[0.1, 50, 30, 20, 1]))},
             "mass": {(1.0, 0.9, 1e3): pd.DataFrame(dict(m_bh=[8, 15, 30, 60], m_b_ev=[1.1e-12, 9e-13, 4e-13, 2e-13]))}}
    r = cw.distance_reach(boson, 7e-13, 1.2e-12)
    assert r.loc[0, "d_max_kpc_in_band"] == 30 and r.loc[0, "d_max_kpc_min_in_band"] == 20
    e = cw.mass_exclusion_in_band(boson, 7e-13, 1.2e-12)
    assert (e.loc[0, "m_bh_min"], e.loc[0, "m_bh_max"]) == (8, 15)


def test_narrowband(tmp_path, monkeypatch):
    """Post-glitch analyses flagged; ratio to the spin-down limit from the targeted table."""
    p = tmp_path / "nb.txt"
    p.write_text("J0534+2200 59.14 4.16e-26\nJ0537-6910 123.78 3.52e-26\nJ0537-6910pg 123.78 3.75e-26\n")
    monkeypatch.setattr(cw, "fetch", lambda name, refresh=False: p)
    kp = pd.DataFrame(dict(pulsar=["J0534+2200", "J0537-6910"], h0_spindown=[1.43e-24, 2.91e-26]))
    d = cw.narrowband(kp)
    assert list(d["post_glitch"]) == [False, False, True] and list(d["pulsar"]) == ["J0534+2200", "J0537-6910", "J0537-6910"]
    assert d.loc[0, "ratio_to_spindown"] == pytest.approx(0.0291, rel=1e-2) and d.loc[1, "ratio_to_spindown"] == pytest.approx(1.21, rel=1e-2)

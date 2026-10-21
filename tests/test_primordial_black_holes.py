"""Offline tests of the primordial_black_holes mode: merger rates, suppression, limits and the mass windows."""
from __future__ import annotations

import io
import zipfile

import numpy as np
import pandas as pd
import pytest

from gwtc_analysis import primordial_black_holes as pbh


def test_suppression_models():
    """The LVK fit equals S1 S2 at f = 1 and falls below it at small f; S2 saturates at 1 and S1 < 1."""
    f = np.array([1e-4, 1e-3, 1e-2, 0.1, 1.0])
    s_h, s_l = pbh.suppression_hutsi(f), pbh.suppression_lvk(f)
    assert s_h[-1] == pytest.approx(s_l[-1], rel=0.05)
    assert s_h[2] / s_l[2] == pytest.approx(2.9, rel=0.05)
    assert np.all(s_h < 1)
    assert s_l[0] < 1 < pbh.suppression_lvk(1e-5)          # the fit exceeds 1 below f ~ 1e-4


def test_limit_inverts_the_rate():
    for sup in pbh.SUPPRESSIONS:
        for m in (0.35, 1.0, 30.0):
            r = float(pbh.rate_early(m, 3e-3, sup))
            assert pbh.f_pbh_limit(m, r, "early", sup) == pytest.approx(3e-3, rel=1e-4)
    assert pbh.f_pbh_limit(10.0, float(pbh.rate_late(10.0, 0.02)), "late") == pytest.approx(0.02)
    assert np.isnan(pbh.f_pbh_limit(1.0, 0.0))


def test_rate_upper_limit():
    assert pbh.rate_upper_limit(0, 1.0) == pytest.approx(-np.log(0.1))     # 2.30 / VT, as in the LVK papers
    assert pbh.rate_upper_limit(3, 2.0) == pytest.approx(6.681 / 2.0, rel=1e-3)


def _ssm_zip(tmp_path):
    """A results-PBH.zip with the first rows of the O4b release (PyCBC; the same rows for the others)."""
    txt = ("mc_min mc_max VT_O4a efficiency_O4a variance_O4a neff_O4a VT_O4b efficiency_O4b variance_O4b neff_O4b "
           "VT_O4ab efficiency_O4ab variance_O4ab neff_O4ab\n"
           "0.16 0.44444 0 0 0 0 0 0 0 0 0.00079975 0 0 508.52\n"
           "0.72889 1.0133 0 0 0 0 0 0 0 0 0.0072779 0 0 3000\n"
           "2.4400 2.7200 0 0 0 0 0 0 0 0 0.005 0 0 300\n")
    p = tmp_path / "results-PBH.zip"
    with zipfile.ZipFile(p, "w") as z:
        for key in pbh.SSM_PIPELINES.values():
            z.writestr(f"PBH/vt-results-final/{key}_pbh_vt_inputs_mc.txt", txt)
    return p


def test_ssm_reproduces_the_lvk(tmp_path):
    """With the LVK fit, the PyCBC VT of the first bin gives their 0.096 at 0.35 Msun; bins above 1 Msun are left out."""
    vt = pbh.read_ssm_vt(_ssm_zip(tmp_path))
    assert set(vt["pipeline"]) == set(pbh.SSM_PIPELINES)
    t = pbh.ssm_limits(vt[vt["pipeline"] == "PyCBC"], "lvk")
    assert t["m_pbh"].max() <= 1.01 and len(t) == 2
    assert t["m_pbh"].iloc[0] == pytest.approx(0.347, abs=1e-3)
    assert t["f_early"].iloc[0] == pytest.approx(0.0958, rel=0.01)
    assert t["f_late"].iloc[0] == pytest.approx(2.55, rel=0.01)
    assert pbh.ssm_limits(vt[vt["pipeline"] == "PyCBC"], "hutsi")["f_early"].iloc[0] < t["f_early"].iloc[0]


def _injections(n=200000, seed=2):
    """Found injections drawn log-uniform in both masses on [1, 100] (m1 >= m2) and uniform in z, unit weights."""
    rng = np.random.default_rng(seed)
    a, b = np.exp(rng.uniform(0, np.log(100), (2, n)))
    m1, m2 = np.maximum(a, b), np.minimum(a, b)
    z = rng.uniform(0.01, 1.0, n)
    ln_draw = np.log(2.0) - np.log(m1) - np.log(m2) - 2 * np.log(np.log(100)) - np.log(0.99)
    spin = np.zeros(n)
    return dict(mass1_source=m1, mass2_source=m2, redshift=z, weights=np.ones(n), lnpdraw=ln_draw,
                a1=np.full(n, 0.1), a2=np.full(n, 0.1), dvdz=np.ones(n), T_yr=1.0, N_gen=float(n),
                spin1x=spin, spin1y=spin, spin1z=spin, spin2x=spin, spin2y=spin, spin2z=spin)


def test_window_vt_selects_the_window():
    inj = _injections()
    vt, neff = pbh.window_vt(inj, 10.0, 0.2, evolve=False)
    assert neff > 500
    # dV/dz = 1 Gpc^3, all found: VT = T * <1 / (1 + z)> * (spin density ratio), the same in every window
    vt2, _ = pbh.window_vt(inj, 30.0, 0.2, evolve=False)
    assert vt == pytest.approx(vt2, rel=0.15)
    vt_e, _ = pbh.window_vt(inj, 10.0, 0.2, evolve=True)
    assert vt_e > vt                                   # (t/t0)^(-34/37) > 1 at z > 0


def test_events_in_window_use_the_chirp_mass():
    ev = pd.DataFrame(dict(event=["A", "B", "C"], mass_1_source=[10.0, 14.0, 30.0], mass_2_source=[10.0, 7.0, 30.0],
                           chirp_mass_source=[8.71, None, 26.1]))
    sel = pbh.events_in_window(ev, 10.0, 0.2)          # chirp mass of 10 + 10: 8.71; B: 8.61 from its masses
    assert list(sel["event"]) == ["A", "B"]


def test_catalog_limits_attribution():
    inj = _injections()
    ev = pd.DataFrame(dict(event=["A"], mass_1_source=[10.0], mass_2_source=[10.0], chirp_mass_source=[8.71]))
    w = pbh.catalog_limits(inj, ev, np.array([10.0, 40.0]), 0.2, "window")
    n = pbh.catalog_limits(inj, ev, np.array([10.0, 40.0]), 0.2, "none")
    assert list(w["n_events"]) == [1, 0] and list(n["n_events"]) == [1, 0]
    assert w["rate_90"].iloc[0] > n["rate_90"].iloc[0]              # the event loosens the limit
    assert w["f_early"].iloc[1] == pytest.approx(n["f_early"].iloc[1])

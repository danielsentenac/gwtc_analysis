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


def test_suppression_values():
    """S1 S2 at hand-checked points (monochromatic, Hütsi et al. 2021 eqs. 2.2, 2.3, 2.7, A.5)."""
    assert float(pbh.suppression_hutsi(0.1)) == pytest.approx(0.0131, rel=0.03)
    assert float(pbh.suppression_hutsi(1e-3)) == pytest.approx(0.419, rel=0.03)
    assert float(pbh.suppression_hutsi(1e-6)) < float(pbh.suppression_hutsi(1e-3))   # S1 falls at small f


def test_late_rate_independent_of_mass_for_equal_masses():
    r = pbh.rate_late(np.array([0.5, 5.0, 50.0]), 1.0)
    assert np.allclose(r, 2 ** (10 / 7) * pbh.R_LATE_CLUSTERS)


def test_early_rate_falls_with_mass():
    r = pbh.rate_early(np.array([1.0, 10.0, 100.0]), 1e-3)
    assert np.all(np.diff(r) < 0)
    assert r[0] / r[1] == pytest.approx(10 ** (32 / 37), rel=1e-6)


def test_read_ssm_vt_needs_every_pipeline(tmp_path):
    p = tmp_path / "results-PBH.zip"
    with zipfile.ZipFile(p, "w") as z:
        z.writestr("PBH/vt-results-final/pycbc_pbh_vt_inputs_mc.txt", "mc_min mc_max VT_O4ab neff_O4ab\n0.16 0.44 1e-3 100\n")
    with pytest.raises(ValueError, match="GstLAL"):
        pbh.read_ssm_vt(p)


def test_cli_choices_match_the_module():
    from gwtc_analysis import cli

    sp = cli.build_parser()._mode_parsers["primordial_black_holes"]
    choices = {a.dest: a.choices for a in sp._actions if a.choices}
    assert tuple(choices["suppression"]) == pbh.SUPPRESSIONS
    assert sp.get_default("suppression") == pbh.SUPPRESSIONS[0]
    assert tuple(choices["attribute"]) == ("window", "none")


def test_cli_dispatch(monkeypatch):
    from gwtc_analysis import cli

    seen = {}
    monkeypatch.setattr(pbh, "run_primordial_black_holes", lambda **kw: seen.update(kw))
    assert cli.main(["primordial_black_holes", "--masses", "3", "30", "--mass-window", "0.1",
                     "--suppression", "lvk", "--attribute", "none", "--no-sub-solar"]) == 0
    assert seen["masses"] == [3.0, 30.0] and seen["mass_window"] == 0.1
    assert seen["suppression"] == "lvk" and seen["attribute"] == "none" and seen["sub_solar"] is False


@pytest.mark.parametrize("kw, msg", [(dict(mass_window=1.5), "mass-window"), (dict(masses=[1.05]), "e\\^-w"),
                                     (dict(suppression="x"), "suppression"), (dict(attribute="x"), "attribute")])
def test_run_rejects_bad_options(kw, msg):
    with pytest.raises(ValueError, match=msg):
        pbh.run_primordial_black_holes(**kw)


def test_run_end_to_end_offline(tmp_path, monkeypatch):
    """The whole mode on synthetic injections and events and the sub-solar zip: TSV, report and plot."""
    from gwtc_analysis import catalogs

    inj = _injections(seed=4)
    inj["segments"] = [(0.0, 1.0)]
    ev = pd.DataFrame(dict(event=["GWA", "GWB"], mass_1_source=[10.0, 31.0], mass_2_source=[9.0, 29.0],
                           chirp_mass_source=[8.27, 26.1], far_per_yr=[1e-3, 1e-2]))
    zip_path = _ssm_zip(tmp_path)
    monkeypatch.setattr(pbh, "fetch_ssm", lambda refresh=False: zip_path)
    monkeypatch.setattr(catalogs, "_rates_sensitivity_path", lambda f, r: tmp_path / "inj.hdf")
    monkeypatch.setattr(catalogs, "_load_found_injections", lambda p, far: inj)
    monkeypatch.setattr(catalogs, "_rates_events", lambda seg, far, ns, lists: ev)
    out = tmp_path / "out"
    t = pbh.run_primordial_black_holes(masses=[2.0, 10.0, 30.0], out_limits_tsv=out / "l.tsv",
                                       out_report_html=out / "r.html", plots_dir=out / "plots")
    assert (out / "l.tsv").exists() and (out / "plots" / "f_pbh.png").exists()
    html = (out / "r.html").read_text()
    assert "Hütsi" in html and "arXiv:2610.20276" in html
    cat = t[t["source"] == "GWTC injections and events"].set_index("m_pbh")
    assert list(cat["n_events"]) == [0, 1, 1]
    assert set(t[t["source"] != "GWTC injections and events"]["pipeline"]) == set(pbh.SSM_PIPELINES)
    assert (cat["f_early"] > 0).all() and (cat["f_late"] > 0).all()

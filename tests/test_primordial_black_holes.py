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
            z.writestr(f"PBH/vt-results-final/{key}_pbh_vt_inputs_m1m2.txt", M1M2)
    return p


# first rows of the GstLAL (m1, m2) table of the O4b release, and its row without O4a values (no injections)
M1M2 = ("m1_min m1_max m2_min m2_max VT_O4a efficiency_O4a variance_O4a neff_O4a VT_O4b efficiency_O4b variance_O4b "
        "neff_O4b VT_O4ab efficiency_O4ab variance_O4ab neff_O4ab\n"
        "0.19 2.372 0.19 0.372 0.0005936 0.000359 1.7e-10 743.7 0.00050628 0.00015156 5.5e-11 414.23 0.00091707 "
        "0.00026304 8.4412e-11 819.65\n"
        "0.19 2.372 0.918 1.1 0.005 0.003 1e-8 700 0.005 0.001 1e-9 500 0.010943 0.003 1e-9 1100\n"
        "8.918 11.1 0.19 0.372     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n")


def test_ssm_map_is_the_lvk_formula(tmp_path):
    """F = 2.30 / VT / (dln m1 dln m2 R(m1, m2; F = 1)) at the arithmetic bin centre; the empty bin is NaN."""
    vt2 = pbh.read_ssm_vt_m1m2(_ssm_zip(tmp_path))
    g = pbh.ssm_map(vt2[vt2["pipeline"] == "GstLAL"]).reset_index(drop=True)
    assert len(g) == 3 and np.isnan(g["vt_gpc3_yr"].iloc[2]) and np.isnan(g["F_early"].iloc[2])
    m1, m2 = 1.281, 0.281
    mt = m1 + m2
    dln = np.log(2.372 / 0.19) * np.log(0.372 / 0.19)
    r90 = -np.log(0.1) / 0.00091707
    eb = 1.6e6 * 2.3e-3 * mt ** (-32 / 37) * (m1 * m2 / mt ** 2) ** (-34 / 37)
    lb = 420.0 * mt ** (10 / 7) / (m1 * m2) ** (5 / 7)
    assert g["F_early"].iloc[0] == pytest.approx(r90 / (dln * eb), rel=1e-9)
    assert g["F_late"].iloc[0] == pytest.approx(r90 / (dln * lb), rel=1e-9)


def test_ordered_ln_area():
    d = np.log(2.0)
    assert pbh.ordered_ln_area((1, 2), (1, 2)) == pytest.approx(d * d / 2, rel=1e-4)       # triangle
    assert pbh.ordered_ln_area((4, 8), (1, 2)) == pytest.approx(d * d, rel=1e-6)           # whole rectangle
    assert pbh.ordered_ln_area((1, 2), (4, 8)) == pytest.approx(0.0, abs=1e-9)             # all m1 < m2
    rng = np.random.default_rng(0)                                                          # straddling bin
    u1, u2 = rng.uniform(np.log(0.19), np.log(2.372), 400000), rng.uniform(np.log(0.918), np.log(1.1), 400000)
    mc = np.mean(u1 >= u2) * np.log(2.372 / 0.19) * np.log(1.1 / 0.918)
    assert pbh.ordered_ln_area((0.19, 2.372), (0.918, 1.1)) == pytest.approx(mc, rel=0.01)


def test_region_vt_is_window_vt_on_the_diagonal():
    inj = _injections()
    x = pbh._injection_terms(inj, evolve=False)
    w = 0.2
    b = (10 * np.exp(-w), 10 * np.exp(w))
    assert pbh.region_vt(inj, b, b, x) == pytest.approx(pbh.window_vt(inj, 10.0, w, evolve=False), rel=1e-3)


def _dense_points(vt, lo=0.2, hi=2.0, step=0.004):
    """Exposure points of a constant VT on a fine (ln m1, ln m2) grid, m1 >= m2."""
    u = np.arange(np.log(lo), np.log(hi), step)
    m1, m2, c = pbh.grid_exposure(u, u, vt)
    return pbh._ExposurePoints(m1, m2, c, c)


def test_grid_exposure_counts_half_the_diagonal():
    u = np.linspace(0.0, 1.0, 11)
    _, _, c = pbh.grid_exposure(u, u, 1.0)
    assert c.sum() == pytest.approx(2.0 * 0.5, rel=1e-3)          # 2 VT x the triangle of area 1/2


@pytest.mark.parametrize("sup", pbh.SUPPRESSIONS)
def test_narrow_lognormal_is_monochromatic(sup):
    """As sigma -> 0 the log-normal limit tends to the monochromatic one at m_c, for both channels."""
    vt, mc = 2e-3, 0.6
    pts = _dense_points(vt)
    f_e, f_l = pbh.lognormal_limit(pts, mc, 0.03, {}, sup)
    r90 = pbh.rate_upper_limit(0, vt)
    assert f_e == pytest.approx(pbh.f_pbh_limit(mc, r90, "early", sup), rel=0.02)
    assert f_l == pytest.approx(pbh.f_pbh_limit(mc, r90, "late"), rel=0.02)


def test_broad_lognormal_and_events():
    """A wider log-normal at the same m_c spreads over more masses; an event loosens the limit."""
    pts = _dense_points(2e-3, 0.05, 20.0, 0.01)
    f_narrow, _ = pbh.lognormal_limit(pts, 0.6, 0.1, {}, "lvk")
    f_wide, _ = pbh.lognormal_limit(pts, 0.6, 1.0, {}, "lvk")
    assert np.isfinite(f_wide) and f_wide != pytest.approx(f_narrow, rel=0.01)
    f_one, _ = pbh.lognormal_limit(pts, 0.6, 1.0, {pbh.SUB_SOLAR_GROUP: 1}, "lvk")
    assert f_one > f_wide


def test_lognormal_windows():
    """Split in two windows, the limit is set by the window that holds no event; pooled, by all the exposure."""
    pts = _dense_points(2e-3)
    pts.group = np.where(pts.m1 * pts.m2 > 0.36, 1, 2)          # split around m1 m2 = 0.6^2
    f_split, _ = pbh.lognormal_limit(pts, 0.6, 0.3, {1: 50}, "lvk")
    f_only2, _ = pbh.lognormal_limit(pbh._ExposurePoints(pts.m1[pts.group == 2], pts.m2[pts.group == 2],
                                                         pts.c_early[pts.group == 2], pts.c_late[pts.group == 2]),
                                     0.6, 0.3, {}, "lvk")
    assert f_split == pytest.approx(f_only2, rel=1e-9)
    f_pool, _ = pbh.lognormal_limit(pts, 0.6, 0.3, {1: 50}, "lvk", pooled=True)
    assert f_pool < f_split


def test_chirp_window_matches_the_monochromatic_window():
    """The window of an equal-mass m + m binary is centred on it: m = 1 is window 0, m = e^(2w) the next."""
    mc = 2 ** -0.2 * np.array([1.0, np.exp(0.39), np.exp(0.41)])
    assert list(pbh.chirp_window(mc, 0.2)) == [0, 1, 1]
    assert list(pbh.chirp_window(2 ** -0.2 * np.array([np.exp(0.19), np.exp(0.21)]), 0.2)) == [0, 1]


def test_hutsi_suppression_mass_function():
    """S1 falls with M/<m> (more PBHs near the pair) and with the spread <m^2>/<m>^2; the defaults are the
    monochromatic values."""
    f = 1e-2
    assert float(pbh.suppression_hutsi(f, 2.0, 1.0)) == pytest.approx(float(pbh.suppression_hutsi(f)))
    assert float(pbh.suppression_hutsi(f, 4.0)) < float(pbh.suppression_hutsi(f, 2.0))
    assert float(pbh.suppression_hutsi(f, 2.0, np.e)) < float(pbh.suppression_hutsi(f, 2.0, 1.0))


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
                     "--suppression", "lvk", "--attribute", "none", "--no-sub-solar",
                     "--lognormal-widths", "0.3", "--lognormal-masses", "1", "10", "--no-extended"]) == 0
    assert seen["masses"] == [3.0, 30.0] and seen["mass_window"] == 0.1
    assert seen["lognormal_widths"] == [0.3] and seen["lognormal_masses"] == [1.0, 10.0] and seen["extended"] is False
    assert seen["suppression"] == "lvk" and seen["attribute"] == "none" and seen["sub_solar"] is False


@pytest.mark.parametrize("kw, msg", [(dict(mass_window=1.5), "mass-window"), (dict(masses=[1.05]), "e\\^-w"),
                                     (dict(suppression="x"), "suppression"), (dict(attribute="x"), "attribute"),
                                     (dict(lognormal_widths=[0.0]), "widths"),
                                     (dict(lognormal_masses=[-1.0]), "lognormal-masses")])
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
    t = pbh.run_primordial_black_holes(masses=[2.0, 10.0, 30.0], lognormal_masses=[0.5, 5.0, 30.0],
                                       map_edges=np.geomspace(1.1, 100.0, 5), out_limits_tsv=out / "l.tsv",
                                       out_report_html=out / "r.html", plots_dir=out / "plots")
    assert (out / "l.tsv").exists() and (out / "plots" / "f_pbh.png").exists()
    for f in ("l.map.tsv", "l.lognormal.tsv", "plots/f_pbh_map.png", "plots/f_pbh_lognormal.png"):
        assert (out / f).exists(), f
    mp = pd.read_csv(out / "l.map.tsv", sep="\t")
    assert set(mp["source"]) == {"sub-solar searches O4a+O4b", "GWTC injections and events"}
    ln = pd.read_csv(out / "l.lognormal.tsv", sep="\t")
    assert set(ln["pipeline"]) == set(pbh.SSM_PIPELINES) and set(ln["sigma"]) == {0.5, 1.0, 2.0}
    html = (out / "r.html").read_text()
    assert "Hütsi" in html and "arXiv:2610.20276" in html
    cat = t[t["source"] == "GWTC injections and events"].set_index("m_pbh")
    assert list(cat["n_events"]) == [0, 1, 1]
    assert set(t[t["source"] != "GWTC injections and events"]["pipeline"]) == set(pbh.SSM_PIPELINES)
    assert (cat["f_early"] > 0).all() and (cat["f_late"] > 0).all()


def test_run_without_sub_solar(tmp_path, monkeypatch):
    """--no-sub-solar: the catalog part alone still writes the TSV, the report and the plot."""
    from gwtc_analysis import catalogs

    inj = _injections(seed=5)
    inj["segments"] = [(0.0, 1.0)]
    ev = pd.DataFrame(dict(event=["GWA"], mass_1_source=[10.0], mass_2_source=[9.0],
                           chirp_mass_source=[8.27], far_per_yr=[1e-3]))
    monkeypatch.setattr(pbh, "fetch_ssm", lambda refresh=False: pytest.fail("sub-solar release fetched"))
    monkeypatch.setattr(catalogs, "_rates_sensitivity_path", lambda f, r: tmp_path / "inj.hdf")
    monkeypatch.setattr(catalogs, "_load_found_injections", lambda p, far: inj)
    monkeypatch.setattr(catalogs, "_rates_events", lambda seg, far, ns, lists: ev)
    out = tmp_path / "out"
    t = pbh.run_primordial_black_holes(masses=[2.0, 10.0], sub_solar=False, lognormal_masses=[5.0],
                                       lognormal_widths=[0.5], map_edges=np.geomspace(1.1, 100.0, 4),
                                       out_limits_tsv=out / "l.tsv",
                                       out_report_html=out / "r.html", plots_dir=out / "plots")
    assert (out / "l.tsv").exists() and (out / "r.html").exists() and (out / "plots" / "f_pbh.png").exists()
    assert set(t["source"]) == {"GWTC injections and events"}
    ln = pd.read_csv(out / "l.lognormal.tsv", sep="\t")
    assert list(ln["pipeline"]) == ["all (FAR threshold)"] and np.isfinite(ln["f_early"]).all()


def test_first_excluded_finds_the_first_crossing():
    """A count that rises, falls and rises again: the limit is the first crossing, and the dip is reported."""
    count = lambda d, f: d * (f / 1e-4 if f < 1e-3 else 10 * (1e-3 / f) ** 3 if f < 1e-1 else 1e-5 * (f / 0.1) ** 9)
    f_lim, f_re = pbh._first_excluded([(1.0, 2.0)], count)
    assert f_lim == pytest.approx(2e-4, rel=1e-4)
    assert f_re is not None and 1e-3 < f_re < 1e-1
    f_lim2, f_re2 = pbh._first_excluded([(1.0, 2.0), (1.0, 0.5)], count)        # the second window crosses first
    assert f_lim2 == pytest.approx(5e-5, rel=1e-4)
    assert pbh._first_excluded([(1e-20, 2.0)], count)[0] == np.inf

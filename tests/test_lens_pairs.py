"""Offline tests of the lens_pairs mode: overlap integrals, prior Jacobian, lens time delays, ranking."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gwtc_analysis import lens_pairs as lp


def _gauss_ln_overlap(mu1, c1, mu2, c2, ln_pi):
    """ln ∫ N(mu1, c1) N(mu2, c2) / pi with a constant prior density exp(ln_pi)."""
    c = c1 + c2
    d = mu1 - mu2
    sign, logdet = np.linalg.slogdet(2 * np.pi * c)
    return -0.5 * logdet - 0.5 * d @ np.linalg.solve(c, d) - ln_pi


def test_intrinsic_overlap_matches_gaussian_integral():
    rng = np.random.default_rng(3)
    dim, n = 6, 3000
    c1 = np.diag(rng.uniform(0.05, 0.3, dim) ** 2)
    c2 = np.diag(rng.uniform(0.05, 0.3, dim) ** 2)
    mu1 = np.zeros(dim)
    mu2 = mu1 + 0.1
    ya = rng.multivariate_normal(mu1, c1, n)
    yb = rng.multivariate_normal(mu2, c2, n)
    ln_pi = -2.0
    got = lp.log_overlap_intrinsic(ya, np.full(n, ln_pi), yb, np.full(n, ln_pi))
    # the kernel density widens each posterior by (1 + h^2), h the Scott factor
    h2 = n ** (-2 / (dim + 4))
    want = _gauss_ln_overlap(mu1, c1 * (1 + h2), mu2, c2 * (1 + h2), ln_pi)
    assert got == pytest.approx(want, abs=0.25)
    # far apart: much smaller
    yc = rng.multivariate_normal(mu1 + 2.0, c2, n)
    assert lp.log_overlap_intrinsic(ya, np.full(n, ln_pi), yc, np.full(n, ln_pi)) < got - 20


def test_mass_prior_jacobian():
    """ln_pi of intrinsic_coords (mass part) equals ln of the uniform (m1, m2) density times |d(m1,m2)/d(y1,y2)|."""
    m1, m2 = np.array([36.0, 80.0, 12.0]), np.array([29.0, 20.0, 11.0])
    d = dict(mass_1=m1, mass_2=m2, a_1=np.full(3, A := 0.495), a_2=np.full(3, A), cos_tilt_1=np.zeros(3),
             cos_tilt_2=np.zeros(3))
    y, ln_pi = lp.intrinsic_coords(d)
    # spin part at u = 0.5, c = 0: 2 ln(0.25) + 2 ln(0.5)
    ln_mass = ln_pi - 2 * np.log(0.25) - 2 * np.log(0.5)

    def to_m(y1, y2):
        mc, q = np.exp(y1), 1 / (1 + np.exp(-y2))
        m1_ = mc * (1 + q) ** 0.2 / q ** 0.6
        return m1_, q * m1_

    eps = 1e-6
    for i in range(3):
        j = np.empty((2, 2))
        for k in range(2):
            dy = np.zeros(2); dy[k] = eps
            a = np.array(to_m(*(y[i, :2] + dy))); b = np.array(to_m(*(y[i, :2] - dy)))
            j[:, k] = (a - b) / (2 * eps)
        lo, hi = lp.MASS_BOX
        want = np.log(2 / (hi - lo) ** 2) + np.log(abs(np.linalg.det(j)))
        assert ln_mass[i] == pytest.approx(want, abs=1e-5)


def test_sky_overlap():
    import healpy as hp

    npix = hp.nside2npix(8)
    uniform = np.full(npix, 1 / npix)
    assert lp.log_overlap_sky(uniform, uniform) == pytest.approx(0.0, abs=1e-12)
    peaked = np.full(npix, 1e-12); peaked[10] = 1.0; peaked /= peaked.sum()
    assert lp.log_overlap_sky(peaked, peaked) == pytest.approx(np.log(npix), abs=1e-6)


def test_sky_density_concentrates():
    rng = np.random.default_rng(1)
    d = dict(ra=rng.normal(1.0, 0.01, 4000), dec=rng.normal(0.3, 0.01, 4000))
    m = lp.sky_density(d, nside=32, fwhm_deg=1.0)
    assert m.sum() == pytest.approx(1.0)
    assert lp.log_overlap_sky(m, m) > 5            # a few-degree blob against the whole sky


def test_sis_delay():
    # sigma 200 km/s, z_l 0.3, z_s 1, y 0.5: about a month (Planck18)
    assert lp.sis_time_delay_days(200.0, 0.3, 1.0, 0.5) == pytest.approx(29.3, rel=0.02)
    # scales as sigma^4 and y
    assert lp.sis_time_delay_days(400.0, 0.3, 1.0, 0.5) == pytest.approx(16 * 29.3, rel=0.02)
    assert lp.sis_time_delay_days(200.0, 0.3, 1.0, 0.25) == pytest.approx(29.3 / 2, rel=0.02)


def test_lensed_delays_typical_range():
    dt, w = lp.lensed_delays(np.array([0.3, 0.5, 1.0]), n=50_000)
    med = np.exp(np.average(np.log(dt), weights=w))
    assert 1 < med < 100                            # days to weeks for galaxy lenses


def test_delay_ratio_favours_short_delays():
    rng = np.random.default_rng(2)
    t = np.sort(rng.uniform(0, 240, 60))            # 60 events over an 8-month run
    pairs = np.abs(t[:, None] - t[None, :])[np.triu_indices(60, 1)]
    r = lp.DelayRatio(pairs, np.full(60, 0.4))
    ln_r = r.ln_r(np.array([3.0, 200.0]))
    assert ln_r[0] > 0 > ln_r[1]


def _mock_samples(rng, mc, q, ra, dec, n=600, t=0.0):
    m1 = mc * (1 + q) ** 0.2 / q ** 0.6 * rng.normal(1, 0.02, n)
    return dict(mass_1=m1, mass_2=q * m1, a_1=rng.uniform(0, 0.99, n), a_2=rng.uniform(0, 0.99, n),
                cos_tilt_1=rng.uniform(-1, 1, n), cos_tilt_2=rng.uniform(-1, 1, n),
                ra=rng.normal(ra, 0.02, n), dec=rng.normal(dec, 0.02, n), geocent_time=np.full(n, t))


def test_search_ranks_the_lensed_pair_first():
    rng = np.random.default_rng(5)
    names, samples, z = [], {}, []
    for i in range(12):
        nm = f"GW{i:02d}"
        samples[nm] = _mock_samples(rng, rng.uniform(15, 40), rng.uniform(0.4, 0.95), rng.uniform(0, 6),
                                    rng.uniform(-1, 1), t=rng.uniform(0, 2e7))
        names.append(nm); z.append(0.4)
    # a lensed copy of GW03, 10 days later
    s = samples["GW03"]
    samples["GW99"] = _mock_samples(rng, float(np.median(lp._chirp_mass_det(s))),
                                    float(np.median(s["mass_2"] / s["mass_1"])), float(np.median(s["ra"])),
                                    float(np.median(s["dec"])), t=float(s["geocent_time"][0]) + 10 * 86400)
    names.append("GW99"); z.append(0.4)
    events = pd.DataFrame(dict(event=names, redshift=z))
    df, _ = lp.search_pairs(samples, events, progress_every=0)
    assert len(df) == 13 * 12 // 2
    assert {df.iloc[0]["event_1"], df.iloc[0]["event_2"]} == {"GW03", "GW99"}
    assert df.iloc[0]["fpp"] == pytest.approx(1 / len(df))
    assert (df["log10_B_overlap"] >= lp.FLOOR_LOG10).all()


def test_delay_ratio_bounded_in_the_tails():
    rng = np.random.default_rng(4)
    t = np.sort(rng.uniform(0, 240, 60))
    pairs = np.abs(t[:, None] - t[None, :])[np.triu_indices(60, 1)]
    r = lp.DelayRatio(pairs, np.full(60, 0.4))
    ln_r = r.ln_r(np.array([0.01, 0.02, 2000.0, 9000.0]))
    assert np.all(np.abs(ln_r) / np.log(10) < 3)      # |log10 R| < 3 where both densities vanish


def test_window_overlap_and_delay_cutoff():
    seg = [(0.0, 100 * 86400.0), (200 * 86400.0, 300 * 86400.0)]       # two 100-day runs with a 100-day gap
    w = lp.window_overlap_days(np.array([0.0, 50.0, 150.0, 250.0, 400.0]), seg)
    assert w == pytest.approx([200.0, 100.0, 50.0, 50.0, 0.0])
    r = lp.DelayRatio(np.abs(np.random.default_rng(1).uniform(0, 300, (200, 2)) @ [1, -1]), np.full(20, 0.5),
                      segments=[(0.0, 30 * 86400.0)])
    # no lensed pair can be longer than a 30-day run: R there is the floor of the uniform part
    assert r.ln_r(np.array([5.0]))[0] > r.ln_r(np.array([100.0]))[0]

"""Offline tests of the spin_population mode: the isotropic chi_eff density and recovery of a mock population."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gwtc_analysis import spin_population as sp


def _iso_draws(rng, q, n):
    s1 = rng.uniform(0, 1, n) * rng.uniform(-1, 1, n)
    s2 = rng.uniform(0, 1, n) * rng.uniform(-1, 1, n)
    return (s1 + q * s2) / (1 + q)


def test_iso_density_matches_monte_carlo():
    rng = np.random.default_rng(1)
    chi = np.linspace(-1, 1, 401)
    for q in (1.0, 0.5, 0.1):
        p = sp.iso_chi_eff_density(q, chi)
        assert sp._trapz(p, chi) == pytest.approx(1.0, abs=2e-3)
        h, e = np.histogram(_iso_draws(rng, q, 1_000_000), bins=40, range=(-1, 1), density=True)
        assert np.interp(0.5 * (e[1:] + e[:-1]), chi, p) == pytest.approx(h, abs=0.03)
    t = sp.IsoTable(nq=20, nchi=201)
    assert t(np.array([0.1]), np.array([0.7]))[0] == pytest.approx(sp.iso_chi_eff_density(0.7, [0.1])[0], rel=0.02)


def test_trunc_normal_normalized():
    x = np.linspace(-1, 1, 20001)
    for mu, s in ((0.05, 0.1), (0.9, 0.3), (-0.5, 0.6)):
        assert sp._trapz(sp.trunc_normal(x, mu, s), x) == pytest.approx(1.0, abs=1e-3)


def _mock(rng, mu0, sigma, alpha, n_events=120, n_pe=600, meas=0.08):
    """Events whose chi_eff follows the population; PE samples from an isotropic prior times a Gaussian
    measurement; injections from the isotropic prior, all found (no selection)."""
    iso = sp.IsoTable(nq=25, nchi=401)
    events = []
    for _ in range(n_events):
        q = rng.uniform(0.3, 1.0)
        mu = mu0 + alpha * (q - 0.5)
        while True:
            chi = rng.normal(mu, sigma)
            if -1 < chi < 1:
                break
        obs = chi + rng.normal(0, meas)
        draw = _iso_draws(rng, q, 40000)
        w = np.exp(-0.5 * ((draw - obs) / meas) ** 2)
        c = draw[rng.choice(len(draw), n_pe, p=w / w.sum())]
        qq = np.full(n_pe, q)
        events.append(dict(event="E", w=np.ones(n_pe), chi=c, q=qq, iso=np.maximum(iso(c, qq), 1e-12)))
    qi = rng.uniform(0.3, 1.0, 60000)
    ci = (rng.uniform(0, 1, 60000) * rng.uniform(-1, 1, 60000)
          + qi * rng.uniform(0, 1, 60000) * rng.uniform(-1, 1, 60000)) / (1 + qi)
    inj = dict(w=np.ones(len(qi)), chi=ci, q=qi, iso=np.maximum(iso(ci, qi), 1e-12), ntotal=float(len(qi)))
    return dict(events=events, inj=inj)


def test_recovers_mock_population(tmp_path):
    rng = np.random.default_rng(3)
    truth = dict(mu0=0.06, sigma=0.1, alpha=-0.5)
    data = _mock(rng, **truth)
    table = sp.run_spin_population(data=data, n_walkers=16, n_steps=600, out_report_html=tmp_path / "r.html",
                                   out_summary_tsv=tmp_path / "s.tsv", plots_dir=tmp_path / "p").set_index("parameter")
    for k, v in truth.items():
        assert table.loc[k, "low_90"] < v < table.loc[k, "high_90"], k
    post = pd.read_csv(tmp_path / "s.posterior.tsv", sep="\t")
    assert np.mean(post["alpha"] < 0) > 0.9
    assert (tmp_path / "p" / "spin_population.png").exists() and "P(α &lt; 0)" in (tmp_path / "r.html").read_text()


def test_uncorrelated_model():
    rng = np.random.default_rng(4)
    data = _mock(rng, 0.05, 0.12, 0.0, n_events=60)
    lik = sp.Likelihood(data, correlated=False)
    assert np.isfinite(lik.ln_prob([0.05, 0.12])) and lik.ln_prob([0.05, 2.0]) == -np.inf
    post = sp.sample(lik, n_walkers=12, n_steps=400, burn=150)
    assert list(post.columns) == ["mu0", "sigma"] and post["mu0"].quantile(0.05) < 0.05 < post["mu0"].quantile(0.95)

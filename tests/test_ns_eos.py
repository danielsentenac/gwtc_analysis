"""Offline tests of the joint neutron-star EOS mode: relations, the Lambda_1.4 estimator on mock events, report."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gwtc_analysis import ns_eos as ne


def test_relations():
    """Lambda~ = Lambda for equal masses and deformabilities; Lambda(1.4) = Lambda_1.4; radius relation."""
    a, b = ne.tilde_coefficients(np.array([1.4]), np.array([1.4]))
    assert a[0] + b[0] == pytest.approx(1.0)
    assert ne.lambda_of_m([1.4, 2.8], 300)[0] == pytest.approx([300, 300 / 64])
    assert ne.radius_14(2.88e-6 * 12.0 ** 7.5) == pytest.approx(12.0)
    assert 10.5 < ne.radius_14(190) < 11.5


def test_trapezoid_prior_is_normalized():
    a, b = 0.6, 0.45
    y = np.linspace(-10, (a + b) * ne.LAMBDA_MAX + 10, 200001)
    p = ne.prior_tilde(y, a, b)
    assert ne._trapz(p, y) == pytest.approx(1.0, abs=1e-3)
    rng = np.random.default_rng(1)
    s = a * rng.uniform(0, ne.LAMBDA_MAX, 400000) + b * rng.uniform(0, ne.LAMBDA_MAX, 400000)
    h, e = np.histogram(s, bins=40, density=True)
    assert h == pytest.approx(ne.prior_tilde(0.5 * (e[1:] + e[:-1]), a, b), rel=0.08, abs=2e-6)


def _mock(rng, truth, m1c, m2c, sigma, n=150_000, keep=5000):
    """PE-like samples: Lambda_1, Lambda_2 uniform on [0, 5000], a Gaussian likelihood of Lambda~ around the
    value of the common-EOS model at Lambda_1.4 = truth, with noise."""
    m1, m2 = rng.normal(m1c, 0.05, n), rng.normal(m2c, 0.05, n)
    l1, l2 = rng.uniform(0, ne.LAMBDA_MAX, n), rng.uniform(0, ne.LAMBDA_MAX, n)
    a, b = ne.tilde_coefficients(m1, m2)
    lt = a * l1 + b * l2
    true_lt = float(np.mean(a * ne.lambda_of_m(m1, truth)[0] + b * ne.lambda_of_m(m2, truth)[0]))
    w = np.exp(-0.5 * ((lt - (true_lt + rng.normal(0, sigma))) / sigma) ** 2)
    i = rng.choice(n, keep, p=w / w.sum())
    return m1[i], m2[i], lt[i]


def test_estimator_recovers_lambda_14():
    """Over 20 mock GW170817-like events the truth is inside the 90% interval about 90% of the time and at the
    posterior median on average (the posterior median itself is pulled up by the bound at 0 for small truths);
    dividing by the prior at the samples keeps small Lambda_1.4 finite."""
    rng = np.random.default_rng(7)
    grid = np.linspace(0, 3000, 601)
    for truth in (150, 400):
        u = []
        for _ in range(20):
            p = ne.normalize(grid, ne.event_ln_likelihood(grid, *_mock(rng, truth, 1.46, 1.27, 150)))
            c = np.cumsum(p); c /= c[-1]
            u.append(np.interp(truth, grid, c))
        u = np.array(u)
        assert np.mean((u > 0.05) & (u < 0.95)) >= 0.8
        assert 0.35 < u.mean() < 0.65             # the truth sits at the posterior median on average


def test_joint_tighter_than_each(tmp_path, monkeypatch):
    rng = np.random.default_rng(3)
    fake = {"GW170817": _mock(rng, 250, 1.46, 1.27, 150), "GW190425": _mock(rng, 250, 1.75, 1.56, 250)}

    def load(name, spin, cache, pe_cache):
        m1, m2, lt = fake[name]
        a, b = ne.tilde_coefficients(m1, m2)
        return dict(name=name, label=f"C0:{name}", file="x", m1=m1, m2=m2, l1=lt, l2=lt, lam_tilde=lt)

    monkeypatch.setattr(ne, "load_event", load)
    t = ne.run_ns_eos(out_report_html=tmp_path / "r.html", out_summary_tsv=tmp_path / "s.tsv",
                      plots_dir=tmp_path / "p").set_index("analysis")
    width = lambda k: t.loc[k, "high_90"] - t.loc[k, "low_90"]
    assert width("joint") < width("GW170817") and width("joint") < width("GW190425")
    assert 100 < t.loc["joint", "median"] < 500 and 9 < t.loc["joint", "r_median"] < 13
    assert (tmp_path / "s.posterior.tsv").exists() and (tmp_path / "p" / "ns_eos.png").exists()
    assert "common equation of state" in (tmp_path / "r.html").read_text()
    with pytest.raises(ValueError, match="no tidal analysis"):
        ne.run_ns_eos(events=["GW200115"], out_report_html=None, out_summary_tsv=None)


def test_cli():
    from gwtc_analysis.cli import build_parser

    a = build_parser().parse_args(["neutron_star_eos", "--spin-prior", "high", "--events", "GW170817"])
    assert a.spin_prior == "high" and a.events == ["GW170817"] and a.lambda_max == 5000.0

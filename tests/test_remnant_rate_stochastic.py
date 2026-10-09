"""Offline tests of the remnant section, the rate evolution and the stochastic-background prediction."""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
import pytest

from gwtc_analysis import rate_evolution as rev
from gwtc_analysis import remnants as rm
from gwtc_analysis import stochastic as st


# ---------------------------------------------------------------------------
# remnants
# ---------------------------------------------------------------------------
def test_remnant_summary(tmp_path):
    rng = np.random.default_rng(1)
    n = 3000
    s = dict(final_mass_source=rng.normal(61.5, 1.8, n), final_spin=rng.normal(0.685, 0.03, n),
             radiated_energy=rng.normal(3.0, 0.25, n), peak_luminosity=rng.normal(3.6, 0.2, n),
             mass_1_source=rng.normal(35.0, 1.0, n), mass_2_source=rng.normal(29.5, 1.0, n))
    table, files, html = rm.analyse({"C01:Mixed": s, "C01:Other": {"mass_1": np.ones(3)}}, "C01:Mixed", "GWX", tmp_path)
    assert set(table["quantity"]) == {"final mass", "final spin", "radiated energy", "peak luminosity",
                                      "radiated fraction"}
    frac = table[table["quantity"] == "radiated fraction"]["median"].iloc[0]
    assert frac == pytest.approx(3.0 / 64.5, rel=0.03)
    assert [p.name for p in files] == ["GWX_remnant.tsv", "remnant_GWX.png"]
    assert "5.4e+54 erg" in html and "numerical-relativity" in html
    assert "stores no remnant" in rm.analyse({"A": {"mass_1": np.ones(3)}}, "A", "GWY", tmp_path)[2]


# ---------------------------------------------------------------------------
# rate evolution
# ---------------------------------------------------------------------------
def test_md_shape_matches_icarogw_and_md14():
    """psi(0) = 1; the icarogw formula; the star-formation peak at z ≈ 1.86 and R(1)/R(0) ≈ 5.8."""
    z = np.array([0.0, 0.5, 1.0, 3.0])
    g, k, p = 3.3, 2.9, 2.7
    ref = np.exp(np.log1p((1 + p) ** (-g - k)) + g * np.log1p(z) - np.log1p(((1 + z) / (1 + p)) ** (g + k)))
    assert rev.md_shape(z, g, k, p)[0] == pytest.approx(ref)
    assert rev.md_shape(0.0, g, k, p)[0, 0] == pytest.approx(1.0)
    zz = np.linspace(0, 6, 60001)
    assert zz[np.argmax(rev.md_shape(zz, **rev.MD14)[0])] == pytest.approx(rev.z_peak(2.7, 2.9, 1.9), abs=1e-3)
    assert rev.z_peak(2.7, 2.9, 1.9) == pytest.approx(1.86, abs=0.01)
    assert rev.md_shape(1.0, **rev.MD14)[0, 0] == pytest.approx(5.79, abs=0.02)
    assert rev.md_shape(z, [1, 2], [1, 1], [1, 1]).shape == (2, 4)


def test_rate_report(tmp_path):
    rng = np.random.default_rng(2)
    post = pd.DataFrame(dict(gamma=rng.normal(3.3, 0.5, 500), kappa=rng.uniform(0, 6, 500),
                             zp=rng.uniform(1, 4, 500)))
    s = rev.summary(post)
    assert s["gamma"][1] == pytest.approx(3.3, abs=0.1) and s["R(1)/R(0)"][1] > 5
    assert (tmp_path / "rz.png") == rev.plot(post, tmp_path / "rz.png", z_max_events=1.0)
    assert "follows the prior" in rev.report_paragraph(s, 1.0)


# ---------------------------------------------------------------------------
# stochastic background
# ---------------------------------------------------------------------------
def test_energy_spectra():
    """Inspiral energy to the ISCO is the Newtonian binding energy there; the IMR frequencies of GW150914."""
    m1 = m2 = 1.4
    m, mu = m1 + m2, m1 * m2 / (m1 + m2)
    v2 = (1 / 6)                                   # (pi G M f_isco / c^3)^(2/3) = 1/6
    assert st.radiated_energy(m1, m2, st.dE_df_inspiral) == pytest.approx(0.5 * mu * v2, rel=0.01)
    fr = st._ajith_freqs(np.array([65 * st.MSUN]), np.array([36 * 29 / 65 ** 2]))
    assert fr["merg"][0] == pytest.approx(124, abs=2) and fr["ring"][0] == pytest.approx(248, abs=3)
    assert 3.0 < st.radiated_energy(36, 29) < 5.0     # Ajith et al. 2008 overestimates NR's ~3.0
    f = np.array([10.0, 20.0])
    e = st.dE_df_imr(f, 30, 30)[0]
    assert e[1] / e[0] == pytest.approx(2 ** (-1 / 3), rel=1e-6)


def test_omega_scales_as_f_two_thirds_and_rate():
    m1, m2 = np.full(10, 1.4), np.full(10, 1.4)
    spec = st.mean_spectrum(m1, m2, st.dE_df_inspiral)
    f = np.array([10.0, 20.0, 40.0])
    om = st.omega_gw(f, lambda z: np.ones_like(z), spec)
    assert om[1] / om[0] == pytest.approx(2 ** (2 / 3), rel=0.01) and om[2] / om[1] == pytest.approx(2 ** (2 / 3), rel=0.01)
    assert st.omega_gw(f, lambda z: 3 * np.ones_like(z), spec) == pytest.approx(3 * om)


def test_bbh_rate_high_z():
    """Normalized at z_ref; with high_z = sfr, continuous at z_h and with the star-formation shape beyond."""
    z = np.array([0.2, 0.99, 1.0, 1.01, 2.0, 3.0])
    a = st.bbh_rate(z, 3.3, 2.9, 2.7, 25.0, 0.2, "posterior", 1.0)
    b = st.bbh_rate(z, 3.3, 2.9, 2.7, 25.0, 0.2, "sfr", 1.0)
    assert a[0] == pytest.approx(25.0) and b[0] == pytest.approx(25.0)
    assert b[1:3] == pytest.approx(a[1:3]) and b[3] == pytest.approx(b[2], rel=0.02)
    sfr = rev.md_shape(z, **rev.MD14)[0]
    assert b[5] / b[4] == pytest.approx(sfr[5] / sfr[4])
    with pytest.raises(ValueError, match="high_z"):
        st.bbh_rate(z, 3.3, 2.9, 2.7, 25.0, 0.2, "flat", 1.0)


def test_plp_sampler():
    rng = np.random.default_rng(3)
    m1, m2 = st.sample_plp(rng, 5000, 3.4, 1.1, 5.1, 87.0, 4.8, 34.0, 3.6, 0.04)
    assert (m2 <= m1 + 1e-9).all() and m1.min() >= 5.1 and m1.max() <= 87.0 + 1e-6
    assert 0.05 < np.mean(np.abs(m1 - 34) < 7) < 0.25


def test_run_stochastic(tmp_path):
    rng = np.random.default_rng(4)
    n = 40
    post = pd.DataFrame(dict(H0=70.0, alpha=rng.normal(3.4, 0.2, n), beta=1.1, mmin=5.0, mmax=87.0, delta_m=4.8,
                             mu_g=34.0, sigma_g=3.6, lambda_peak=0.04, gamma=rng.normal(2.7, 0.3, n), kappa=2.9, zp=1.9))
    post.to_csv(tmp_path / "posterior.tsv", sep="\t", index=False)
    pd.DataFrame([dict(population="BNS", model="", z_ref=0.0, rate_median=100.0, rate_05=20.0, rate_95=300.0),
                  dict(population="NSBH", model="", z_ref=0.0, rate_median=30.0, rate_05=10.0, rate_95=70.0),
                  dict(population="BBH", model="", z_ref=0.2, rate_median=25.0, rate_05=22.0, rate_95=28.0)]
                 ).to_csv(tmp_path / "rates.tsv", sep="\t", index=False)
    table = st.run_stochastic(tmp_path, tmp_path / "rates.tsv", n_draws=20, out_report_html=tmp_path / "r.html",
                              out_summary_tsv=tmp_path / "s.tsv", plots_dir=tmp_path / "p")
    t = table.set_index("population")
    assert 1e-10 < t.loc["BBH", "omega_25Hz_median"] < 3e-9
    assert t.loc["total", "omega_25Hz_median"] == pytest.approx(
        t.loc[["BBH", "BNS", "NSBH"], "omega_25Hz_median"].sum(), rel=0.3)
    assert (tmp_path / "s.spectrum.tsv").exists() and (tmp_path / "p" / "omega_gw.png").exists()
    assert "Systematics" in (tmp_path / "r.html").read_text()
    with pytest.raises(ValueError, match="Power Law"):
        pd.DataFrame(dict(H0=[70.0])).to_csv(tmp_path / "bad.tsv", sep="\t", index=False)
        st.run_stochastic(tmp_path / "bad.tsv", tmp_path / "rates.tsv", out_report_html=None, out_summary_tsv=None)


MLTP = dict(alpha=3.5, beta=1.1, mmin=5.0, mmax=90.0, delta_m=4.0, mu_g_low=10.0, sigma_g_low=1.0, lambda_g_low=0.6,
            mu_g_high=33.0, sigma_g_high=4.0, lambda_g=0.3)


def test_mltp_sampler_and_model_detection():
    """The Multi Peak sampler puts the expected weight in its two peaks; the model is read from the columns."""
    rng = np.random.default_rng(3)
    m1, m2 = st.sample_mltp(rng, 20000, *(MLTP[k] for k in ("alpha", "beta", "mmin", "mmax", "delta_m", "mu_g_low",
                                                              "sigma_g_low", "lambda_g_low", "mu_g_high",
                                                              "sigma_g_high", "lambda_g")))
    assert (m2 <= m1 + 1e-9).all() and m1.min() >= 5.0
    assert 0.45 < np.mean(np.abs(m1 - 10) < 2) < 0.6          # icarogw: 0.50 between 8 and 12
    assert 0.15 < np.mean((m1 > 25) & (m1 < 40)) < 0.25      # icarogw: 0.20
    common = dict(alpha=3.4, beta=1.1, mmin=5.0, mmax=87.0, delta_m=4.8, gamma=2.7, kappa=2.9, zp=1.9)
    assert st.mass_model_of(pd.DataFrame([dict(common, mu_g=34.0, sigma_g=3.6, lambda_peak=0.04)])) == "plp"
    assert st.mass_model_of(pd.DataFrame([dict(common, **{k: MLTP[k] for k in st.MASS_MODELS["mltp"]["params"]})])) \
        == "mltp"
    with pytest.raises(ValueError, match="neither"):
        st.mass_model_of(pd.DataFrame([common]))


def test_run_stochastic_mltp(tmp_path):
    rng = np.random.default_rng(5)
    n = 30
    post = pd.DataFrame(dict(H0=70.0, gamma=rng.normal(2.7, 0.3, n), kappa=2.9, zp=1.9, **MLTP))
    post.to_csv(tmp_path / "posterior.tsv", sep="\t", index=False)
    pd.DataFrame([dict(population="BNS", model="", z_ref=0.0, rate_median=100.0, rate_05=20.0, rate_95=300.0),
                  dict(population="NSBH", model="", z_ref=0.0, rate_median=30.0, rate_05=10.0, rate_95=70.0),
                  dict(population="BBH", model="", z_ref=0.2, rate_median=25.0, rate_05=22.0, rate_95=28.0)]
                 ).to_csv(tmp_path / "rates.tsv", sep="\t", index=False)
    t = st.run_stochastic(tmp_path, tmp_path / "rates.tsv", n_draws=15, out_report_html=tmp_path / "r.html",
                          out_summary_tsv=None, plots_dir=tmp_path / "p").set_index("population")
    assert 1e-10 < t.loc["BBH", "omega_25Hz_median"] < 3e-9
    assert "Multi Peak masses" in (tmp_path / "r.html").read_text()


@pytest.mark.skipif(not os.environ.get("GWTC_ICAROGW_PYTHON"), reason="set GWTC_ICAROGW_PYTHON to the icarogw interpreter")
def test_mltp_sampler_matches_icarogw(tmp_path):
    """The m1 distribution of sample_mltp is icarogw's massprior_MultiPeak with the low-mass smoothing."""
    import subprocess

    code = (
        "import sys, json, numpy as np\n"
        f"sys.path.insert(0, {str(tmp_path)!r})\n"
        "import icarogw\n"
        f"p = json.loads({json.dumps(json.dumps(MLTP))})\n"
        "w = icarogw.wrappers.m1m2_conditioned_lowpass(icarogw.wrappers.massprior_MultiPeak()); w.update(**p)\n"
        "g = np.linspace(5, 120, 20000); pg = np.exp(w.prior.pdf1.log_pdf(g))\n"
        "out = [float(np.trapz(pg[(g > a) & (g < b)], g[(g > a) & (g < b)])) for a, b in ((8, 12), (25, 40), (40, 90))]\n"
        "print('OUT', json.dumps(out))\n")
    (tmp_path / "config.py").write_text("CUPY=False\n")
    r = subprocess.run([os.environ["GWTC_ICAROGW_PYTHON"], "-c", code], capture_output=True, text=True, cwd=tmp_path)
    line = [ln for ln in r.stdout.splitlines() if ln.startswith("OUT")]
    assert r.returncode == 0 and line, r.stderr[-1500:]
    ref = json.loads(line[0][4:])
    rng = np.random.default_rng(1)
    m1, _ = st.sample_mltp(rng, 200000, *(MLTP[k] for k in ("alpha", "beta", "mmin", "mmax", "delta_m", "mu_g_low",
                                                              "sigma_g_low", "lambda_g_low", "mu_g_high",
                                                              "sigma_g_high", "lambda_g")))
    for (a, b), r_ in zip(((8, 12), (25, 40), (40, 90)), ref):
        assert np.mean((m1 > a) & (m1 < b)) == pytest.approx(r_, abs=0.005)

"""Offline tests for the `bright_siren` mode: cosmology, likelihood, selection, report, CLI, and mock bright
sirens (simulated detections with known H0, which the analysis must recover)."""
from __future__ import annotations

import json
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy.stats import kstest

from gwtc_analysis import bright_siren as bs
from gwtc_analysis.cli import build_parser

CP = bs.COUNTERPARTS["GW170817"]


# ---------------------------------------------------------------------------
# building blocks
# ---------------------------------------------------------------------------
def test_gw170817_hubble_flow_redshift():
    """v_H = 3327 - 310 = 3017 ± 166 km/s, as in LVK 2017."""
    z, s, _ = bs.hubble_flow_redshift(CP)
    assert z * bs.C_KMS == pytest.approx(3017) and s * bs.C_KMS == pytest.approx(166, abs=0.5)
    assert bs.hubble_flow_redshift(CP, redshift=(0.1, 0.01))[:2] == (0.1, 0.01)


def test_cosmology():
    """Linear Hubble law at low z, astropy's flat ΛCDM at high z, and the inverse."""
    from astropy.cosmology import FlatLambdaCDM

    assert bs.luminosity_distance(1e-4, 70) == pytest.approx(1e-4 * bs.C_KMS / 70, rel=1e-4)
    ref = FlatLambdaCDM(H0=70, Om0=bs.OM0, Tcmb0=0).luminosity_distance([0.01, 0.438, 2.0]).value
    assert bs.luminosity_distance(np.array([0.01, 0.438, 2.0]), 70) == pytest.approx(ref, rel=1e-5)
    assert bs.redshift_at(ref, 70) == pytest.approx([0.01, 0.438, 2.0], rel=1e-5)


def test_low_redshift_limit_is_the_lvk_formula():
    """With a d_L^2 PE prior and the euclidean selection, the posterior is Σ N(v_H; H0 d_i, σ) (LVK 2017),
    up to the ΛCDM correction of d_L (0.8% at z = 0.01)."""
    rng = np.random.default_rng(1)
    d = rng.normal(43.0, 4.0, 4000)
    h0 = np.linspace(40, 120, 801)
    p = bs.posterior_from_ln(h0, bs.event_ln_likelihood(h0, d, d ** 2, 3017 / bs.C_KMS, 166 / bs.C_KMS)
                             - bs.ln_selection_euclidean(h0))
    lvk = np.exp(-0.5 * ((3017 - h0[:, None] * d[None, :]) / 166) ** 2).mean(axis=1)
    assert bs.summarize(h0, p)["median"] == pytest.approx(bs.summarize(h0, lvk)["median"] * 1.0077, rel=0.003)


def test_summarize_gaussian():
    """MAP, 68.3% HPD and 90% interval of a Gaussian."""
    h0 = np.linspace(10, 200, 19001)
    s = bs.summarize(h0, np.exp(-0.5 * ((h0 - 70) / 5) ** 2))
    assert s["map"] == pytest.approx(70, abs=0.02)
    assert s["hpd68_low"] == pytest.approx(65, abs=0.1) and s["hpd68_high"] == pytest.approx(75, abs=0.1)
    assert s["low_90"] == pytest.approx(70 - 1.645 * 5, abs=0.05)
    with pytest.raises(ValueError, match="zero everywhere"):
        bs.posterior_from_ln(h0, np.full(len(h0), -np.inf))


def test_sky_conditioning():
    """Samples fixed to the counterpart are all kept; otherwise only those near it, or an error."""
    n = 1000
    ra0, dec0 = np.radians(CP.ra_deg), np.radians(CP.dec_deg)
    d = np.linspace(20, 50, n)
    keep, how = bs.sky_conditioned(dict(luminosity_distance=d, ra=np.full(n, ra0), dec=np.full(n, dec0)), CP, 3.0)
    assert keep.all() and "fixed" in how
    spread = dict(luminosity_distance=d, ra=ra0 + np.radians(np.linspace(-10, 10, n)), dec=np.full(n, dec0))
    keep, how = bs.sky_conditioned(spread, CP, 3.0)
    assert 250 < keep.sum() < 330 and "within 3.0°" in how
    with pytest.raises(ValueError, match="increase --sky-radius"):
        bs.sky_conditioned(spread, CP, 0.5)


def test_spectral_density_and_reader(tmp_path):
    """A spectral posterior is read from a TSV or a work directory, and its density integrates to 1."""
    rng = np.random.default_rng(2)
    pd.DataFrame({"H0": np.clip(rng.normal(90, 30, 4000), 10.5, 199.5)}).to_csv(tmp_path / "posterior.tsv", sep="\t",
                                                                                  index=False)
    assert bs.read_spectral_posterior(tmp_path).shape == (4000,)
    h0 = np.linspace(10, 200, 3801)
    assert bs._trapz(bs.spectral_density(h0, bs.read_spectral_posterior(tmp_path / "posterior.tsv")), h0) == \
        pytest.approx(1, abs=0.02)
    (tmp_path / "empty").mkdir()
    with pytest.raises(ValueError, match="no posterior"):
        bs.read_spectral_posterior(tmp_path / "empty")
    # a dark-siren work directory is named as such (from its selection.json)
    assert bs.siren_kind(tmp_path) == "spectral siren"
    (tmp_path / "selection.json").write_text(json.dumps({"galaxy_catalog": {"band": "K-glade+"}}))
    assert bs.siren_kind(tmp_path).startswith("dark siren (K-glade+")
    assert bs.siren_kind(tmp_path / "posterior.tsv").startswith("dark siren")


def _pe_file(path, labels=("C02:Test-HighSpin", "C02:Test-LowSpin")):
    rng = np.random.default_rng(3)
    with h5py.File(path, "w") as h:
        for lab in labels:
            n = 3000
            arr = np.zeros(n, dtype=[("luminosity_distance", "f8"), ("ra", "f8"), ("dec", "f8"), ("theta_jn", "f8")])
            arr["luminosity_distance"] = rng.normal(43.0, 3.0, n)
            arr["theta_jn"] = np.radians(rng.uniform(0, 180, n))
            arr["ra"], arr["dec"] = np.radians(CP.ra_deg), np.radians(CP.dec_deg)
            h.create_group(lab).create_dataset("posterior_samples", data=arr)
    return path


def test_run_bright_siren_end_to_end(tmp_path):
    """Both labels (LowSpin first), the combination, the TSVs and the report."""
    pe = _pe_file(tmp_path / "pe.h5")
    rng = np.random.default_rng(4)
    pd.DataFrame({"H0": np.clip(rng.normal(90, 30, 4000), 10.5, 199.5)}).to_csv(tmp_path / "spec.tsv", sep="\t",
                                                                                  index=False)
    table = bs.run_bright_siren(pe_file=pe, spectral_posterior=tmp_path / "spec.tsv",
                                out_report_html=tmp_path / "r.html", out_summary_tsv=tmp_path / "s.tsv",
                                plots_dir=tmp_path / "plots")
    assert list(table["analysis"])[:2] == ["bright siren, C02:Test-LowSpin", "bright siren, C02:Test-HighSpin"]
    assert table["map"].iloc[0] == pytest.approx(3017 / 43.0 * 1.0077, rel=0.03)
    assert "spectral" in table["analysis"].iloc[-1]
    grid = pd.read_csv(tmp_path / "s.posterior.tsv", sep="\t")
    assert {"H0", "p_spectral", "p_combined"} <= set(grid.columns)
    assert (tmp_path / "plots" / "h0_bright_siren_GW170817.png").exists()
    assert (tmp_path / "plots" / "distance_inclination_GW170817.png").exists()
    assert "By viewing angle" in (tmp_path / "r.html").read_text()
    assert (tmp_path / "plots" / "h0_combined.png").exists()
    assert "NGC 4993" in (tmp_path / "r.html").read_text()
    with pytest.raises(ValueError, match="not in pe.h5"):
        bs.run_bright_siren(pe_file=pe, pe_labels=["C02:Nope"], out_report_html=None, out_summary_tsv=None)
    with pytest.raises(ValueError, match="unknown selection"):
        bs.run_bright_siren(pe_file=pe, selection="none", out_report_html=None, out_summary_tsv=None)


def test_viewing_angle_and_degeneracy_table():
    """theta_jn folded to 0-90 degrees; H0 implied by the host redshift; the table by viewing angle."""
    v = bs.viewing_angle(dict(theta_jn=np.radians([10.0, 170.0, 90.0])))
    assert v == pytest.approx([10.0, 10.0, 90.0])
    assert bs.viewing_angle(dict(iota=np.radians([30.0]))) == pytest.approx([30.0])
    assert bs.viewing_angle(dict(luminosity_distance=[1.0])) is None
    z = 3017 / bs.C_KMS
    assert bs.implied_h0(np.array([43.0]), z)[0] == pytest.approx(3017 / 43.0 * 1.0077, rel=0.002)
    d = np.array([45.0, 44.0, 36.0, 23.0])
    tab = bs.degeneracy_table(d, np.array([10.0, 20.0, 45.0, 70.0]), z)
    assert tab["fraction"].tolist() == [0.5, 0.25, 0.25] and tab["h0_median"].is_monotonic_increasing


def test_cli_parses_bright_siren():
    args = build_parser().parse_args(["bright_siren", "--v-peculiar", "300", "100", "--pe-label", "A", "B",
                                      "--selection", "injections"])
    assert args.src_name == "GW170817" and args.v_peculiar == [300.0, 100.0] and args.pe_label == ["A", "B"]
    assert args.h0_range == [10.0, 200.0] and args.selection == "injections"
    assert build_parser().parse_args(["bright_siren", "--src-name", "GW190521"]).src_name == "GW190521"


# ---------------------------------------------------------------------------
# mock bright sirens
#
# Sources uniform in comoving volume, masses uniform in [20, 40] M☉, isotropic orientations. The network SNR
# is ρ = A(M_c,det) w(ι) / d_L plus unit Gaussian noise, detected above 12. Each detected event gets the
# posterior of d_L given its observed SNR (d_L^2 prior, inclination marginalized: the real degeneracy) and
# its host redshift with an error of 0.001. Injections are drawn with broader masses than the population, as
# the LVK ones, and carry their draw density in the detector frame.
# ---------------------------------------------------------------------------
ZMAX, MLO, MHI, ILO, IHI, RHO_TH, SIGZ = 4.0, 20.0, 40.0, 5.0, 120.0, 12.0, 0.001
_ZG = np.linspace(0, ZMAX, 40001)
_CDF = np.cumsum(bs.pop_z(_ZG))
_CDF /= _CDF[-1]
_ZNORM = bs._trapz(bs.pop_z(_ZG), _ZG)


def _ln_pop_mass(m1, m2):
    ok = (m1 >= MLO) & (m1 <= MHI) & (m2 >= MLO) & (m2 <= MHI)
    return np.where(ok, -2 * np.log(MHI - MLO), -np.inf)


def _w(ci):
    return np.sqrt(((1 + ci ** 2) / 2) ** 2 + ci ** 2) / np.sqrt(2)


def _amp(m1d, m2d, rho_ref):
    mc = (m1d * m2d) ** 0.6 / (m1d + m2d) ** 0.2
    return rho_ref * (mc / 26.0) ** (5 / 6) * 1000.0


def _sources(rng, n, h0, rho_ref, mlo=MLO, mhi=MHI, zmax=ZMAX):
    cdf = _CDF / np.interp(zmax, _ZG, _CDF)
    z = np.interp(rng.uniform(size=n), cdf, _ZG)
    m1, m2 = rng.uniform(mlo, mhi, n), rng.uniform(mlo, mhi, n)
    d = bs.luminosity_distance(z, h0)
    a = _amp(m1 * (1 + z), m2 * (1 + z), rho_ref)
    rho = a * _w(rng.uniform(-1, 1, n)) / d + rng.normal(size=n)
    return z, m1, m2, d, a, rho


def _injections(rng, n, h0, rho_ref, zmax=ZMAX):
    z, m1, m2, d, _, rho = _sources(rng, n, h0, rho_ref, ILO, IHI, zmax)
    f = rho > RHO_TH
    pdraw = bs.pop_z(z) / (_ZNORM * np.interp(zmax, _ZG, _CDF)) / (IHI - ILO) ** 2
    prior = pdraw / ((1 + z) ** 2 * bs.ddl_dz(z, h0))
    return dict(mass_1=(m1 * (1 + z))[f], mass_2=(m2 * (1 + z))[f], luminosity_distance=d[f], prior=prior[f],
                ntotal=float(n))


def _distance_samples(rng, a, rho_obs, ns=2000, m=30000):
    """Posterior of d_L given rho_obs: d_L^2 prior, uniform cos ι, proposal d = A w / r with r ~ N(rho_obs, 1)."""
    ci = rng.uniform(-1, 1, m)
    r = rng.normal(rho_obs, 1.0, m)
    ci, r = ci[r > 1], r[r > 1]
    d = a * _w(ci) / r
    wt = d ** 4 / _w(ci)
    return d[rng.choice(len(d), ns, p=wt / wt.sum())]


def _events(rng, n_ev, h0, rho_ref):
    out = []
    while len(out) < n_ev:
        z, _, _, _, a, rho = _sources(rng, 20000, h0, rho_ref)
        for i in np.flatnonzero(rho > RHO_TH)[:n_ev - len(out)]:
            out.append((z[i] + rng.normal(0, SIGZ), _distance_samples(rng, a[i], rho[i])))
    return out


def _median_sd(h0, ln_p):
    p = bs.posterior_from_ln(h0, ln_p)
    med = bs.summarize(h0, p)["median"]
    return med, float(np.sqrt(bs._trapz(p * (h0 - med) ** 2, h0)))


def test_mock_selection_from_injections_matches_direct_simulation():
    """β(H0) from injections made at H0 = 70 and reweighted equals the detected fraction simulated at each H0,
    up to the constant normalization of p_pop(z)."""
    rng = np.random.default_rng(2)
    h0 = np.linspace(40, 120, 5)
    ln_b, neff = bs.ln_selection_injections(h0, _injections(rng, 600_000, 70.0, 43.0), _ln_pop_mass, n_grid=5)
    direct = np.array([np.mean(_sources(rng, 300_000, h, 43.0)[-1] > RHO_TH) for h in h0])
    ratio = np.exp(ln_b) / direct / _ZNORM
    assert ratio == pytest.approx(1.0, abs=0.03)
    assert direct[-1] / direct[0] > 10            # the selection effect is strong


def test_mock_catalog_recovers_h0_and_needs_the_selection():
    """100 detections at z ~ 0.5 with H0 = 70: unbiased with the selection term, biased without it."""
    rng = np.random.default_rng(1)
    h0 = np.linspace(40, 120, 201)
    events = _events(rng, 100, 70.0, 43.0)
    ln_num = np.sum([bs.event_ln_likelihood(h0, d, d ** 2, z, SIGZ) for z, d in events], axis=0)
    ln_b, _ = bs.ln_selection_injections(h0, _injections(rng, 600_000, 70.0, 43.0), _ln_pop_mass)
    med, sd = _median_sd(h0, ln_num - len(events) * ln_b)
    assert abs(med - 70) < 2.5 * sd and sd < 2.0
    med0, sd0 = _median_sd(h0, ln_num)
    assert med0 - 70 > 3 * sd0                    # without selection: biased high


def test_mock_posteriors_are_calibrated():
    """H0 drawn from the prior for each event: the posterior CDF at the truth is uniform (P-P test)."""
    rng = np.random.default_rng(5)
    h0 = np.linspace(50, 100, 251)
    ln_b, _ = bs.ln_selection_injections(h0, _injections(rng, 600_000, 70.0, 43.0), _ln_pop_mass)
    u = []
    for _ in range(100):
        truth = rng.uniform(50, 100)
        (z, d), = _events(rng, 1, truth, 43.0)
        p = bs.posterior_from_ln(h0, bs.event_ln_likelihood(h0, d, d ** 2, z, SIGZ) - ln_b)
        cdf = np.cumsum(p)
        u.append(np.interp(truth, h0, cdf / cdf[-1]))
    assert kstest(u, "uniform").pvalue > 0.01


def test_mock_nearby_sources_euclidean_selection():
    """For nearby sources (GW170817-like) β from injections scales as H0^3, the euclidean selection."""
    rng = np.random.default_rng(8)
    h0 = np.linspace(40, 120, 41)
    ln_b, _ = bs.ln_selection_injections(h0, _injections(rng, 400_000, 70.0, 1.1, zmax=0.1), _ln_pop_mass)
    assert np.polyfit(np.log(h0), ln_b, 1)[0] == pytest.approx(3.0, abs=0.15)

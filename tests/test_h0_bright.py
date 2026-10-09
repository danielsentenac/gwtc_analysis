"""Offline tests for `hubble_constant --method bright` and `--method joint`: cosmology, likelihood, selection,
work directory, joint posterior, CLI, and mock bright sirens (simulated detections with known H0, which the
analysis must recover)."""
from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy.stats import kstest

from gwtc_analysis import h0_bright as bs
from gwtc_analysis import h0_joint as hj
from gwtc_analysis.cli import _resolve_h0_method, build_parser

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


def test_sample_density_and_kind(tmp_path):
    """Spectral samples give a density that integrates to 1; a dark-siren work directory is named as such."""
    rng = np.random.default_rng(2)
    h0 = np.linspace(10, 200, 3801)
    assert bs._trapz(hj.sample_density(h0, np.clip(rng.normal(90, 30, 4000), 10.5, 199.5)), h0) == \
        pytest.approx(1, abs=0.02)
    assert hj.siren_kind(tmp_path) == "spectral siren"
    (tmp_path / "selection.json").write_text(json.dumps({"galaxy_catalog": {"band": "K-glade+"}}))
    assert hj.siren_kind(tmp_path).startswith("dark siren (K-glade+")
    assert hj.siren_kind(tmp_path / "posterior.tsv").startswith("dark siren")


def _pe_file(path, labels=("C02:Test-HighSpin", "C02:Test-LowSpin"), masses=False):
    rng = np.random.default_rng(3)
    cols = [("luminosity_distance", "f8"), ("ra", "f8"), ("dec", "f8"), ("theta_jn", "f8")]
    with h5py.File(path, "w") as h:
        for lab in labels:
            n = 3000
            arr = np.zeros(n, dtype=cols + ([("mass_1", "f8"), ("mass_2", "f8")] if masses else []))
            if masses:
                arr["mass_1"], arr["mass_2"] = rng.normal(1.48, 0.05, n), rng.normal(1.27, 0.05, n)
            arr["luminosity_distance"] = rng.normal(43.0, 3.0, n)
            arr["theta_jn"] = np.radians(rng.uniform(0, 180, n))
            arr["ra"], arr["dec"] = np.radians(CP.ra_deg), np.radians(CP.dec_deg)
            h.create_group(lab).create_dataset("posterior_samples", data=arr)
    return path


def _spectral_workdir(path, events=("GW150914_095045",), common=("GW150914",), mean=90.0):
    rng = np.random.default_rng(4)
    path.mkdir()
    pd.DataFrame({"H0": np.clip(rng.normal(mean, 30, 4000), 10.5, 199.5)}).to_csv(path / "posterior_reweighted.tsv",
                                                                                   sep="\t", index=False)
    pd.DataFrame({"event": list(events), "common_name": list(common)}).to_csv(path / "events.tsv", sep="\t",
                                                                               index=False)
    (path / "summary.json").write_text(json.dumps({"mass_model": "plp"}))
    return path


def test_run_h0_bright_end_to_end(tmp_path):
    """Both labels (LowSpin first), the work directory, the TSV and the report."""
    pe = _pe_file(tmp_path / "pe.h5")
    wd = tmp_path / "bright"
    table = bs.run_h0_bright(pe_file=pe, workdir=wd, out_report_html=tmp_path / "r.html",
                             out_summary_tsv=tmp_path / "s.tsv")
    assert list(table["analysis"]) == ["bright siren GW170817, C02:Test-LowSpin",
                                       "bright siren GW170817, C02:Test-HighSpin"]
    assert table["map"].iloc[0] == pytest.approx(3017 / 43.0 * 1.0077, rel=0.03)
    grid = pd.read_csv(wd / bs.GRID_FILE, sep="\t")
    assert {"H0", "p", "p_C02:Test-LowSpin"} <= set(grid.columns) and len(grid) == 3801
    info = json.loads((wd / bs.BRIGHT_FILE).read_text())
    assert info["method"] == "bright" and info["label"] == "C02:Test-LowSpin" and info["h0_range"] == [10.0, 200.0]
    assert (wd / "plots" / "h0_posterior.png").exists()
    html = (tmp_path / "r.html").read_text()
    assert "NGC 4993" in html and "counterpart" in html
    with pytest.raises(ValueError, match="not in pe.h5"):
        bs.run_h0_bright(pe_file=pe, pe_labels=["C02:Nope"], workdir=wd, out_report_html=None, out_summary_tsv=None)
    with pytest.raises(ValueError, match="unknown selection"):
        bs.run_h0_bright(pe_file=pe, selection="none", workdir=wd, out_report_html=None, out_summary_tsv=None)
    with pytest.raises(ValueError, match="unknown population"):
        bs.run_h0_bright(pe_file=pe, population="bbh", workdir=wd, out_report_html=None, out_summary_tsv=None)
    with pytest.raises(ValueError, match="no detector-frame masses"):
        bs.run_h0_bright(pe_file=pe, population="fullpop4", selection="euclidean", workdir=wd, out_report_html=None,
                         out_summary_tsv=None)


def test_fullpop4_matches_icarogw():
    """FullPop-4.0 at the GWTC-4.0 medians against icarogw's m1m2_paired_massratio_bplmulti_dip (values relative to
    the first point, computed with icarogw): the BNS region, the dip, the low peak, the BBH range, near mmax."""
    m1 = np.array([1.4, 1.03, 4.5, 9.0, 30.0, 93.79])
    m2 = np.array([1.3, 1.02, 1.4, 8.0, 25.0, 20.0])
    lp = bs.ln_mass_fullpop4(m1, m2)
    assert lp - lp[0] == pytest.approx([0.0, -0.23554, -5.282925, -4.824762, -10.50374, -19.719131], abs=1e-5)
    assert np.isneginf(bs.ln_mass_fullpop4(np.array([1.3, 0.9, 95.0]), np.array([1.4, 0.8, 10.0]))).all()
    z = np.array([0.0, 1.0, 2.59, 5.0])
    psi = (1 + z) ** 3.57 / (1 + ((1 + z) / 3.59) ** 6.54)
    assert bs.ln_rate_madau(z, **bs.MADAU_GWTC4) == pytest.approx(np.log(psi / psi[0]))


def test_population_and_distance_prior_options(tmp_path):
    """--population fullpop4 adds the mass and rate terms (a small change at 43 Mpc) and is recorded; an explicit
    distance prior equal to the default leaves the result unchanged."""
    pe = _pe_file(tmp_path / "pe.h5", labels=("C02:Test-LowSpin",), masses=True)
    kw = dict(pe_file=pe, out_report_html=None, out_summary_tsv=None, selection="euclidean")
    base = bs.run_h0_bright(workdir=tmp_path / "a", **kw)
    pop = bs.run_h0_bright(workdir=tmp_path / "b", population="fullpop4", **kw)
    same = bs.run_h0_bright(workdir=tmp_path / "c", pe_distance_prior="dl2", **kw)
    assert pop["median"].iloc[0] == pytest.approx(base["median"].iloc[0], rel=0.02)
    assert pop["median"].iloc[0] != base["median"].iloc[0]
    assert same["median"].iloc[0] == pytest.approx(base["median"].iloc[0], rel=1e-9)
    info = json.loads((tmp_path / "b" / bs.BRIGHT_FILE).read_text())
    assert info["population"] == "fullpop4" and info["selection"] == "euclidean"


def test_viewing_angle_constraint_narrows_h0(tmp_path):
    """A viewing-angle constraint weights the samples: with distance and inclination correlated, it narrows H0
    and moves it to the distance of the constrained angles."""
    rng = np.random.default_rng(5)
    n = 6000
    view = rng.uniform(0, 80, n)
    d = 47.0 - 0.3 * view + rng.normal(0, 1.0, n)          # inclined orbits look closer
    arr = np.zeros(n, dtype=[("luminosity_distance", "f8"), ("ra", "f8"), ("dec", "f8"), ("theta_jn", "f8")])
    arr["luminosity_distance"], arr["theta_jn"] = d, np.radians(view)
    arr["ra"], arr["dec"] = np.radians(CP.ra_deg), np.radians(CP.dec_deg)
    with h5py.File(tmp_path / "pe.h5", "w") as h:
        h.create_group("C02:Test-LowSpin").create_dataset("posterior_samples", data=arr)
    free = bs.run_h0_bright(pe_file=tmp_path / "pe.h5", workdir=tmp_path / "a", out_report_html=None,
                            out_summary_tsv=None).iloc[0]
    con = bs.run_h0_bright(pe_file=tmp_path / "pe.h5", workdir=tmp_path / "b", out_report_html=None,
                           out_summary_tsv=None, viewing_angle_constraint=(20.0, 3.0)).iloc[0]
    assert con["high_90"] - con["low_90"] < free["high_90"] - free["low_90"]
    assert con["median"] == pytest.approx(3017 / (47.0 - 6.0) * 1.0077, rel=0.06)
    assert json.loads((tmp_path / "b" / bs.BRIGHT_FILE).read_text())["viewing_angle"] == [20.0, 3.0]
    del arr
    with h5py.File(tmp_path / "noview.h5", "w") as h:
        a = np.zeros(300, dtype=[("luminosity_distance", "f8")])
        a["luminosity_distance"] = 43.0
        h.create_group("L").create_dataset("posterior_samples", data=a)
    with pytest.raises(ValueError, match="no viewing angle"):
        bs.run_h0_bright(pe_file=tmp_path / "noview.h5", workdir=tmp_path / "c", out_report_html=None,
                         out_summary_tsv=None, viewing_angle_constraint=(20.0, 3.0))


def test_joint_posterior(tmp_path):
    """The joint posterior is the normalized product; independence is checked."""
    pe = _pe_file(tmp_path / "pe.h5")
    bs.run_h0_bright(pe_file=pe, workdir=tmp_path / "bright", out_report_html=None, out_summary_tsv=None)
    spec = _spectral_workdir(tmp_path / "spec")
    table = hj.run_h0_joint([tmp_path / "bright", spec], workdir=tmp_path / "joint",
                            out_report_html=tmp_path / "j.html", out_summary_tsv=tmp_path / "j.tsv")
    assert table["analysis"].tolist()[1:] == ["spectral siren, plp", "joint"]
    b, sp, j = (table.iloc[k] for k in range(3))
    assert b["map"] < j["map"] < sp["map"] and j["high_90"] - j["low_90"] < b["high_90"] - b["low_90"]
    h0 = np.linspace(10, 200, 3801)
    ref = bs._normalize(h0, hj.read_input(tmp_path / "bright", h0)["density"] * hj.read_input(spec, h0)["density"])
    assert j["median"] == pytest.approx(bs.summarize(h0, ref)["median"], abs=0.1)
    grid = pd.read_csv(tmp_path / "joint" / "posterior_joint.tsv", sep="\t")
    assert bs._trapz(grid["p"], grid["H0"]) == pytest.approx(1, abs=1e-3)
    assert (tmp_path / "joint" / "plots" / "h0_joint.png").exists() and "joint" in (tmp_path / "j.html").read_text()
    # a posterior TSV: samples, or a grid with a p column
    grid[["H0", "p"]].to_csv(tmp_path / "g.tsv", sep="\t", index=False)
    assert hj.read_input(tmp_path / "g.tsv", h0)["method"] == "file"
    with pytest.raises(ValueError, match="at least two"):
        hj.run_h0_joint([spec])
    with pytest.raises(ValueError, match="at most one spectral or dark"):
        hj.run_h0_joint([spec, _spectral_workdir(tmp_path / "spec2", events=("GW190412_053044",))])
    with pytest.raises(ValueError, match="GW170817 in both"):
        hj.run_h0_joint([tmp_path / "bright", _spectral_workdir(tmp_path / "spec3", events=("GW170817",))])
    with pytest.raises(ValueError, match="not a hubble_constant work directory"):
        (tmp_path / "empty").mkdir()
        hj.run_h0_joint([tmp_path / "bright", tmp_path / "empty"])


def test_cli_hubble_constant_methods():
    """--method selects the options; those of another method are refused; defaults are named after the method."""
    p = build_parser()
    sp = p._mode_parsers["hubble_constant"]

    def parse(*argv):
        args = p.parse_args(["hubble_constant", *argv])
        _resolve_h0_method(args, sp)
        return args

    args = parse("--method", "bright", "--v-peculiar", "300", "100", "--pe-label", "A", "B", "--selection",
                 "injections", "--viewing-angle", "20", "5")
    assert args.event == "GW170817" and args.v_peculiar == [300.0, 100.0] and args.pe_label == ["A", "B"]
    assert args.h0_range == [10.0, 200.0] and args.viewing_angle == [20.0, 5.0]
    assert (args.workdir, args.out_report, args.out_summary) == ("hubble_constant_bright", "hubble_constant_bright.html",
                                                                 "hubble_constant_bright.tsv")
    assert parse("--method", "bright", "--event", "GW190521").event == "GW190521"
    assert parse().workdir == "hubble_constant_spectral"
    assert parse("--method", "joint", "--inputs", "a", "b").inputs == ["a", "b"]
    assert parse("--method", "dark", "--galaxy-catalog", "c.hdf5", "--workdir", "w").workdir == "w"
    # later stages take the method from the work directory
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "selection.json").write_text(json.dumps({"galaxy_catalog": {"band": "K"}}))
        assert parse("--method", "dark", "--stages", "sample", "--workdir", d).method == "dark"
        with pytest.raises(ValueError, match="use --method dark"):
            parse("--stages", "report", "--workdir", d)
    for argv, msg in ((("--method", "bright", "--nlive", "50"), "--nlive: not an option of hubble_constant --method bright"),
                      (("--event", "GW190521"), "--event: not an option"),
                      (("--method", "joint", "--inputs", "a", "--pe-cache", "x"), "--pe-cache: not an option"),
                      (("--galaxy-catalog", "c.hdf5"), "use --method dark"),
                      (("--method", "dark",), "needs --galaxy-catalog"),
                      (("--method", "dark", "--stages", "sample", "--workdir", "/nonexistent"), "not been prepared"),
                      (("--method", "joint",), "needs --inputs")):
        with pytest.raises(ValueError, match=msg):
            parse(*argv)


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

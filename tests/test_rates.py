"""Offline tests for the `rates` mode (catalogs.run_merger_rates and helpers).

A small synthetic injection file stands in for the LVK release, and the GWOSC
event lists are replaced by a stub, so no network access is needed.
"""
from __future__ import annotations

import csv

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy.stats import gamma

from gwtc_analysis import catalogs as cat
from gwtc_analysis.cli import build_parser

O3 = (1_238_166_018.0, 1_269_363_618.0)
O4A = (1_368_975_618.0, 1_389_456_018.0)
LNPDRAW = "lnpdraw_mass1_source_mass2_source_redshift_spin1x_spin1y_spin1z_spin2x_spin2y_spin2z"


def _iso_spins(rng, n, amax=0.99):
    a = rng.uniform(0, amax, n)
    cos_t = rng.uniform(-1, 1, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    sin_t = np.sqrt(1 - cos_t**2)
    return a * sin_t * np.cos(phi), a * sin_t * np.sin(phi), a * cos_t, a


def _write_injections(path, n=6000, seed=1):
    """Injections with a known draw density.

    Masses (m1 >= m2): half uniform in [1, 3]^2, half uniform in [1, 60]^2; z uniform in [0, 1];
    isotropic spins with magnitudes uniform in [0, 0.99]. Found = z < 0.3.
    """
    rng = np.random.default_rng(seed)
    low = rng.random(n) < 0.5
    a = np.where(low, rng.uniform(1, 3, n), rng.uniform(1, 60, n))
    b = np.where(low, rng.uniform(1, 3, n), rng.uniform(1, 60, n))
    m1, m2 = np.maximum(a, b), np.minimum(a, b)
    z = rng.uniform(0, 1, n)
    s1x, s1y, s1z, a1 = _iso_spins(rng, n)
    s2x, s2y, s2z, a2 = _iso_spins(rng, n)
    p_mass = 0.5 * (2.0 / 2.0**2) * (m1 <= 3) + 0.5 * (2.0 / 59.0**2)
    t = np.where(rng.random(n) < 0.5, rng.uniform(*O3, n), rng.uniform(*O4A, n))
    names = ["mass1_source", "mass2_source", "redshift", "spin1x", "spin1y", "spin1z", "spin2x", "spin2y",
             "spin2z", "weights", "time_geocenter", LNPDRAW, "o3_gstlal_far", "o4a_gstlal_far"]
    arr = np.zeros(n, dtype=[(k, "f8") for k in names])
    arr["mass1_source"], arr["mass2_source"], arr["redshift"], arr["time_geocenter"] = m1, m2, z, t
    arr["spin1x"], arr["spin1y"], arr["spin1z"] = s1x, s1y, s1z
    arr["spin2x"], arr["spin2y"], arr["spin2z"] = s2x, s2y, s2z
    arr["weights"] = 1.0
    arr[LNPDRAW] = (np.log(p_mass) + np.log(1.0)                  # masses x redshift U[0, 1]
                    + cat._ln_iso_spin(a1, 0.99) + cat._ln_iso_spin(a2, 0.99))
    arr["o3_gstlal_far"] = np.where(z < 0.3, 1e-3, np.inf)        # "found" = nearby
    arr["o4a_gstlal_far"] = np.inf
    with h5py.File(path, "w") as h:
        h.create_dataset("events", data=arr)
        h.attrs["total_analysis_time"] = 2.0 * cat._YEAR_S
        h.attrs["total_generated"] = float(n)
        h.attrs["searches"] = ["o3_gstlal", "o4a_gstlal"]
    return path


def test_injection_segments_split_at_gaps():
    """Observing periods come from the injection times, split at gaps longer than a week."""
    t = np.array([0.0, 1e5, 2e5, 2e5 + 30 * 86400, 2e5 + 31 * 86400])
    assert cat._injection_segments(t) == [(0.0, 2e5), (2e5 + 30 * 86400, 2e5 + 31 * 86400)]


def test_sensitive_vt_reweighting(tmp_path):
    """With p_pop = p_draw over the found injections, <VT> = T * N_found / N_gen."""
    inj = cat._load_found_injections(_write_injections(tmp_path / "inj.hdf"), far_threshold=1.0)
    n_found = len(inj["redshift"])
    assert 0 < n_found < 6000
    # a "population" equal to the draw: remove the spin and volume factors that _sensitive_vt adds
    ln_mass = (inj["lnpdraw"] - cat._ln_iso_spin(inj["a1"], 0.99) - cat._ln_iso_spin(inj["a2"], 0.99)
               - np.log(inj["dvdz"]) + np.log1p(inj["redshift"]))
    vt, neff = cat._sensitive_vt(inj, ln_mass, 0.99, 0.99, kappa=0.0)
    assert vt == pytest.approx(2.0 * n_found / 6000, rel=1e-9)
    assert neff == pytest.approx(n_found)


def test_rate_quantiles_jeffreys():
    """Rate posterior quantiles are Gamma(N + 1/2) / VT."""
    q = cat._rate_quantiles(4, 2.0)
    assert q == pytest.approx(gamma.ppf([0.05, 0.5, 0.95], 4.5) / 2.0)
    assert cat._rate_quantiles(0, 1.0)[2] == pytest.approx(gamma.ppf(0.95, 0.5))


FAKE_LISTS = {
    "GWTC-2.1-confident": {
        "GW190425-v3": dict(commonName="GW190425", GPS=1240215503.0, far=0.034, mass_1_source=2.1, mass_2_source=1.3),
        "GW190521-v3": dict(commonName="GW190521", GPS=1242442967.4, far=2e-4, mass_1_source=98.4, mass_2_source=57.2),
        "GW170817-v3": dict(commonName="GW170817", GPS=1187008882.4, far=1e-7, mass_1_source=1.46, mass_2_source=1.27),
    },
    "GWTC-3-marginal": {
        "GW200105_162426-v2": dict(commonName="GW200105_162426", GPS=1262276684.0, far=0.2,
                                   mass_1_source=9.1, mass_2_source=1.91),
        "GW200322_091133-v1": dict(commonName="GW200322_091133", GPS=1268903511.0, far=140.0,
                                   mass_1_source=38.0, mass_2_source=11.3),
    },
    "GWTC-4.0": {
        "GW230518_125908-v1": dict(commonName="GW230518_125908", GPS=1368449966.2, far=1e-5,
                                   mass_1_source=8.17, mass_2_source=1.45),
        "GW230529_181500-v2": dict(commonName="GW230529_181500", GPS=1369419318.7, far=2.2e-4,
                                   mass_1_source=3.66, mass_2_source=1.42),
        "GW230630_070659-v1": dict(commonName="GW230630_070659", GPS=1372143837.0, far=0.5,
                                   mass_1_source=None, mass_2_source=None),
    },
}


@pytest.fixture
def fake_gwosc(monkeypatch):
    monkeypatch.setattr(cat.gw, "fetch_gwtc_events", lambda c: {"events": FAKE_LISTS.get(c, {})})


def test_rates_events_selection(fake_gwosc):
    """Only events inside the injection periods and below the FAR threshold count, classified by mass."""
    df = cat._rates_events([O3, O4A], far_threshold=1.0, ns_max_mass=2.5)
    got = dict(zip(df["event"], df["class"]))
    assert got == {"GW190425": "BNS", "GW190521": "BBH", "GW200105_162426": "NSBH",
                   "GW230529_181500": "NSBH", "GW230630_070659": "unknown"}
    # GW170817 (O2), GW230518 (ER15, before O4a) and GW200322 (FAR 140/yr) are out


def test_run_merger_rates_end_to_end(tmp_path, fake_gwosc):
    """The rates mode writes the rates TSV, the events TSV and the HTML report."""
    inj = _write_injections(tmp_path / "inj.hdf")
    rates = cat.run_merger_rates(
        out_rates_tsv=tmp_path / "rates.tsv", out_events_tsv=tmp_path / "events.tsv",
        out_report_html=tmp_path / "rates.html", plots_dir=tmp_path / "plots", sensitivity_file=inj,
    )
    assert list(rates["population"]) == ["BNS", "NSBH", "BBH", "BBH", "BBH"]
    assert list(rates["n_detected"]) == [1, 2, 1, 1, 1]
    # the evolving BBH rate is also given at z = 0 (same model, divided by (1 + 0.2)^2.9)
    bbh = rates[rates["population"] == "BBH"].reset_index(drop=True)
    assert list(bbh["z_ref"]) == [0.2, 0.0, 0.0] and bbh["model"][1].endswith("at z = 0")
    assert bbh["rate_median"][1] == pytest.approx(bbh["rate_median"][0] / 1.2 ** 2.9)
    from gwtc_analysis.stochastic import read_rates
    assert read_rates(tmp_path / "rates.tsv")["BBH"][1] == pytest.approx(bbh["rate_median"][0])
    assert (rates["vt_gpc3_yr"] > 0).all()
    assert (rates["rate_05"] < rates["rate_median"]).all() and (rates["rate_median"] < rates["rate_95"]).all()
    with open(tmp_path / "rates.tsv") as f:
        assert len(list(csv.DictReader(f, delimiter="\t"))) == 5
    assert (tmp_path / "events.tsv").exists()
    assert "data:image/png;base64" in (tmp_path / "rates.html").read_text()


def test_cli_rates_arguments():
    """CLI wiring for the rates mode."""
    args = build_parser().parse_args(["rates", "--far-threshold", "0.5", "--bbh-kappa", "3.2"])
    assert args.mode == "rates" and args.far_threshold == 0.5 and args.bbh_kappa == 3.2
    assert args.ns_max_mass == 2.5 and args.bbh_z_ref == 0.2


def test_sensitivity_release_is_retrieved_automatically(tmp_path, monkeypatch):
    """--sensitivity-release resolves the latest Zenodo version and downloads its mixture file once."""
    from gwtc_analysis import data_repo

    listing = [{"record_id": "999", "publication_date": "2026-01-01", "files": [
        {"key": "mixture-semi_o1_o2-real_o3_o4a_o4b-cartesian_spins_X.hdf"},
        {"key": "mixture-real_o3_o4a_o4b-polar_spins_X.hdf"},
        {"key": "mixture-real_o3_o4a_o4b-cartesian_spins_X-clipped.hdf"},
        {"key": "psds-o1234ab.tar.gz"},
    ]}]
    monkeypatch.setattr(data_repo, "zenodo_release_versions", lambda release, **kw: listing)
    monkeypatch.setattr(data_repo, "zenodo_cache_dir", lambda: tmp_path)
    downloads = []

    def fake_download(url, dest, **kw):
        downloads.append(url)
        dest.write_bytes(b"hdf")
        return dest

    monkeypatch.setattr(cat.gw, "_download_with_byte_progress", fake_download)
    p = cat._rates_sensitivity_path(None, "gwtc5")
    assert p.name == "zenodo_999_mixture-real_o3_o4a_o4b-cartesian_spins_X-clipped.hdf"
    assert downloads == ["https://zenodo.org/records/999/files/mixture-real_o3_o4a_o4b-cartesian_spins_X-clipped.hdf?download=1"]
    assert cat._rates_sensitivity_path(None, "gwtc5") == p and len(downloads) == 1   # cached
    with pytest.raises(ValueError, match="Unknown sensitivity release"):
        cat._rates_sensitivity_path(None, "gwtc9")


def test_cli_sensitivity_release_choices():
    """The CLI offers the registered releases, defaulting to the latest catalog."""
    assert build_parser().parse_args(["rates"]).sensitivity_release == cat.RATES_DEFAULT_RELEASE == "gwtc5"
    assert build_parser().parse_args(["rates", "--sensitivity-release", "gwtc4"]).sensitivity_release == "gwtc4"
    with pytest.raises(SystemExit):
        build_parser().parse_args(["rates", "--sensitivity-release", "gwtc3"])


def test_event_selection_plot(tmp_path, fake_gwosc, monkeypatch):
    """event_selection --out-plot writes the PNG of the selected events among all the events."""
    from gwtc_analysis import event_selection as es

    monkeypatch.setattr(es.gw, "fetch_gwtc_events", lambda catalog: {"events": FAKE_LISTS.get(catalog, {})})
    monkeypatch.setattr(es.gw, "events_to_dataframe", lambda ev: __import__("pandas").DataFrame(
        [dict(event_id=k, mass_1_source=v["mass_1_source"], mass_2_source=v["mass_2_source"],
              luminosity_distance=500.0) for k, v in ev.items()]))
    es.run_event_selection(catalogs=["GWTC-4.0"], out_tsv=tmp_path / "sel.tsv", m2_max=3.0,
                           out_plot=tmp_path / "sel.png")
    assert (tmp_path / "sel.png").stat().st_size > 1000
    assert len(open(tmp_path / "sel.tsv").read().splitlines()) == 3      # header + GW230518, GW230529


def test_multi_peak_reduces_to_power_law_peak():
    """Multi Peak with all the peak weight in its high Gaussian is Power Law + Peak; both are normalized."""
    rng = np.random.default_rng(3)
    m1 = rng.uniform(3, 100, 2000)
    m2 = m1 * rng.uniform(0.05, 1, 2000)
    plp = cat._ln_power_law_peak(m1, m2, alpha=3.0, beta=1.5, mmin=4.0, mmax=90.0, delta=3.0,
                                 lam=0.1, mu=33.0, sig=4.0)
    mltp = cat._ln_multi_peak(m1, m2, alpha=3.0, beta=1.5, mmin=4.0, mmax=90.0, delta_m=3.0, lambda_g=0.1,
                              lambda_g_low=0.0, mu_g_low=10.0, sigma_g_low=2.0, mu_g_high=33.0, sigma_g_high=4.0)
    assert np.allclose(mltp, plp)
    # m2 integrated on [2, m1] for each m1 (the edge m2 = m1 on the grid), then m1
    trapz = getattr(np, "trapezoid", None) or np.trapz
    g1, u = np.linspace(2, 150, 3000), np.linspace(0, 1, 1500)
    G1, U = np.meshgrid(g1, u, indexing="ij")
    G2 = 2 + (G1 - 2) * U
    p = np.exp(cat._ln_multi_peak(G1.ravel(), G2.ravel(), alpha=3.0, beta=1.5, mmin=4.0, mmax=90.0, delta_m=3.0,
                                  lambda_g=0.3, lambda_g_low=0.6, mu_g_low=10.0, sigma_g_low=2.0,
                                  mu_g_high=33.0, sigma_g_high=4.0)).reshape(G1.shape)
    assert trapz(trapz(p, u, axis=1) * (g1 - 2), g1) == pytest.approx(1.0, abs=5e-3)


def _write_population_posterior(path, mass_model, n=40, seed=2):
    rng = np.random.default_rng(seed)
    cols = dict(alpha=rng.normal(3.2, 0.3, n), beta=rng.normal(1.2, 0.3, n), mmin=rng.uniform(4, 5, n),
                mmax=rng.uniform(80, 90, n), delta_m=rng.uniform(2, 4, n),
                gamma=rng.normal(3.0, 0.5, n), kappa=rng.uniform(2, 4, n), zp=rng.uniform(1.5, 3, n))
    if mass_model == "plp":
        cols.update(mu_g=rng.normal(34, 1, n), sigma_g=rng.uniform(2, 5, n), lambda_peak=rng.uniform(0.02, 0.1, n))
    else:
        cols.update(mu_g_low=rng.normal(10, 0.5, n), sigma_g_low=rng.uniform(1, 3, n), mu_g_high=rng.normal(34, 1, n),
                    sigma_g_high=rng.uniform(2, 5, n), lambda_g=rng.uniform(0.2, 0.5, n),
                    lambda_g_low=rng.uniform(0.5, 0.9, n))
    pd.DataFrame(cols).to_csv(path, sep="\t", index=False)
    return path


@pytest.mark.parametrize("mass_model", ["plp", "mltp"])
def test_run_merger_rates_population_posterior(tmp_path, fake_gwosc, mass_model):
    """With a population posterior, the BBH rate is reported at z_ref and z = 0 over the posterior draws."""
    inj = _write_injections(tmp_path / "inj.hdf")
    post = _write_population_posterior(tmp_path / f"{mass_model}.tsv", mass_model)
    rates = cat.run_merger_rates(
        out_rates_tsv=tmp_path / "rates.tsv", out_events_tsv=None, out_report_html=tmp_path / "rates.html",
        plots_dir=tmp_path / "plots", sensitivity_file=inj, bbh_mass_model=mass_model,
        population_posterior=post, n_draws=20,
    )
    assert list(rates["population"]) == ["BNS", "NSBH", "BBH", "BBH"]
    bbh = rates[rates["population"] == "BBH"]
    assert list(bbh["z_ref"]) == [0.2, 0.0]
    assert bbh["model"].str.contains("20 draws").all()
    assert (bbh["rate_05"] < bbh["rate_median"]).all() and (bbh["rate_median"] < bbh["rate_95"]).all()
    # the rate grows with redshift (gamma > 0 over the draws)
    assert bbh["rate_median"].iloc[0] > bbh["rate_median"].iloc[1]
    # the stochastic mode still finds the evolving BBH rate in the TSV
    from gwtc_analysis.stochastic import read_rates
    assert read_rates(tmp_path / "rates.tsv")["BBH_z_ref"] == 0.2
    assert "population posterior" in (tmp_path / "rates.html").read_text()


def test_run_merger_rates_mltp_needs_posterior(tmp_path):
    with pytest.raises(ValueError, match="needs --population-posterior"):
        cat.run_merger_rates(out_rates_tsv=tmp_path / "r.tsv", bbh_mass_model="mltp")


def test_population_posterior_columns_checked(tmp_path):
    post = _write_population_posterior(tmp_path / "plp.tsv", "plp")
    with pytest.raises(ValueError, match="Multi Peak population posterior"):
        cat._read_population_posterior(post, "mltp")
    df, path = cat._read_population_posterior(post, "plp")
    assert len(df) == 40 and path == post
    work = tmp_path / "work"
    work.mkdir()
    _write_population_posterior(work / "posterior_reweighted.tsv", "plp")
    assert cat._read_population_posterior(work, "plp")[1] == work / "posterior_reweighted.tsv"


def test_cli_rates_mass_model_arguments():
    args = build_parser().parse_args(["rates", "--mass-model", "mltp", "--population-posterior", "w", "--n-draws", "50"])
    assert args.mass_model == "mltp" and args.population_posterior == "w" and args.n_draws == 50
    assert build_parser().parse_args(["rates"]).mass_model == "plp"

"""Offline tests for the `hubble_constant` mode (preparation, report and CLI; the icarogw sampling itself
needs its own environment and is not run here)."""
from __future__ import annotations

import json

import h5py
import numpy as np
import pandas as pd
import pytest
from astropy.cosmology import FlatLambdaCDM

from gwtc_analysis import hubble_constant as hc
from gwtc_analysis.cli import build_parser

LNPDRAW = hc._LNPDRAW


def test_pe_distance_prior_from_recorded_description():
    """The PE distance prior density is read from its bilby description, with catalog defaults when missing."""
    dl = np.array([500.0, 1000.0, 3000.0])
    p, _ = hc.pe_distance_prior(dl, "PowerLaw(alpha=2, minimum=10, maximum=10000, name='luminosity_distance')", "O3a")
    assert p == pytest.approx(dl ** 2)

    desc = "bilby.gw.prior.UniformSourceFrame(minimum=600.0, maximum=9000.0, cosmology='Planck15_LAL', name='luminosity_distance')"
    p, _ = hc.pe_distance_prior(dl, desc, "O4a")
    c = FlatLambdaCDM(H0=67.9, Om0=0.3065)
    z = np.linspace(0, 5, 20001)
    d = c.luminosity_distance(z).value
    ref = np.interp(dl, d, c.differential_comoving_volume(z).value / (1 + z) / np.gradient(d, z))
    assert p / p[0] == pytest.approx(ref / ref[0], rel=1e-6)

    p_default, how = hc.pe_distance_prior(dl, "", "O4a")          # not recorded: O4 default
    assert p_default == pytest.approx(p) and "not recorded" in how
    assert hc.pe_distance_prior(dl, "", "O2")[0] == pytest.approx(dl ** 2)
    with pytest.raises(ValueError, match="Unsupported PE distance prior"):
        hc.pe_distance_prior(dl, "Uniform(minimum=10, maximum=5000)", "O3a")


def test_choose_label_prefers_the_cosmology_paper_waveforms():
    """O4a SpinTaylor first, then the O1-O3 IMRPhenomXPHM run; mixed-waveform labels are not used."""
    assert hc._choose_label(["C00:Mixed", "C00:IMRPhenomXPHM-SpinTaylor", "C00:SEOBNRv5PHM"]) == "C00:IMRPhenomXPHM-SpinTaylor"
    assert hc._choose_label(["C01:Mixed", "C01:IMRPhenomXPHM", "C01:SEOBNRv4PHM"]) == "C01:IMRPhenomXPHM"
    assert hc._choose_label(["C01:IMRPhenomXPHM:HighSpin", "C01:Mixed"]) == "C01:IMRPhenomXPHM:HighSpin"
    assert hc._choose_label(["C01:IMRPhenomPv2_NRTidal:HighSpin"]) is None


FAKE_LISTS = {
    "GWTC-1-confident": {
        "GW150914-v3": dict(commonName="GW150914", GPS=1126259462.4, far=1e-7, mass_1_source=35.6, mass_2_source=30.6),
        "GW170817-v3": dict(commonName="GW170817", GPS=1187008882.4, far=1e-7, mass_1_source=1.46, mass_2_source=1.27),
    },
    "GWTC-3-confident": {
        "GW190814-v2": dict(commonName="GW190814", GPS=1249852257.0, far=1e-5, mass_1_source=23.2, mass_2_source=2.59),
        "GW200105_162426-v2": dict(commonName="GW200105_162426", GPS=1262276684.1, far=0.2, mass_1_source=9.1, mass_2_source=1.91),
    },
    "GWTC-3-marginal": {
        "GW200322_091133-v1": dict(commonName="GW200322_091133", GPS=1268903511.0, far=140.0, mass_1_source=38.0, mass_2_source=11.3),
    },
    "GWTC-4.0": {
        "GW230518_125908-v1": dict(commonName="GW230518_125908", GPS=1368449966.2, far=1e-5, mass_1_source=8.17, mass_2_source=1.45),
        "GW230601_224134-v1": dict(commonName="GW230601_224134", GPS=1369694512.0, far=1e-4, mass_1_source=71.0, mass_2_source=49.0),
        "GW231123_135430-v1": dict(commonName="GW231123_135430", GPS=1384782888.6, far=1e-5, mass_1_source=137.0, mass_2_source=103.0),
    },
}


@pytest.fixture
def fake_gwosc(monkeypatch):
    monkeypatch.setattr(hc.gw, "fetch_gwtc_events", lambda c: {"events": FAKE_LISTS.get(c, {})})


def test_event_selection(fake_gwosc):
    """BBHs with FAR < 0.25/yr inside the runs; NS candidates, ER15, loud-noise and excluded events are out."""
    df = hc.select_h0_events("gwtc4", 0.25, 3.0, hc.H0_DEFAULT_EXCLUDE)
    assert list(df["common_name"]) == ["GW150914", "GW230601_224134"]
    assert list(df["run"]) == ["O1", "O4a"]
    assert df["event"].iloc[0] == "GW150914_095045"        # full name from the GPS time
    # GW170817 and GW190814 (NS masses), GW200105 and GW231123 (excluded), GW200322 (FAR), GW230518 (ER15)
    assert len(hc.select_h0_events("gwtc4", 0.25, 3.0, [])) == 3   # GW231123 back in


def _write_mixture(path, n=4000, seed=3):
    """Semi-analytic (O1+O2) and real injections, with known draw densities."""
    rng = np.random.default_rng(seed)
    names = ["mass1_source", "mass2_source", "redshift", "luminosity_distance", "dluminosity_distance_dredshift",
             "spin1x", "spin1y", "spin1z", "spin2x", "spin2y", "spin2z", "weights", "time_geocenter", LNPDRAW,
             "semianalytic_observed_phase_maximized_snr_net", "o3_gstlal_far"]
    a = np.zeros(n, dtype=[(k, "f8") for k in names])
    a["mass1_source"], a["mass2_source"] = rng.uniform(20, 60, n), rng.uniform(5, 20, n)
    a["redshift"] = rng.uniform(0.05, 1.0, n)
    a["luminosity_distance"] = 5000 * a["redshift"]
    a["dluminosity_distance_dredshift"] = 5000.0
    for k in ("spin1x", "spin1y", "spin1z", "spin2x", "spin2y", "spin2z"):
        a[k] = rng.uniform(-0.3, 0.3, n)
    a["weights"] = rng.uniform(0.5, 2.0, n)
    a[LNPDRAW] = rng.normal(-10, 1, n)
    semi = np.arange(n) < n // 2
    a["time_geocenter"] = np.where(semi, rng.uniform(1.17e9, 1.18e9, n), rng.uniform(1.24e9, 1.25e9, n))  # O2 / O3a
    a["semianalytic_observed_phase_maximized_snr_net"] = np.where(semi, rng.uniform(0, 20, n), np.nan)
    a["o3_gstlal_far"] = np.where(semi, np.inf, 10 ** rng.uniform(-3, 1, n))
    with h5py.File(path, "w") as h:
        h.create_dataset("events", data=a)
        h.attrs["total_generated"] = 1e6
        h.attrs["total_analysis_time"] = 2 * hc._YEAR_S
        h.attrs["searches"] = ["o3_gstlal"]
    return a


def test_detector_frame_injections(tmp_path):
    """Found = SNR > 10 (semi-analytic) or FAR < threshold (real); density in (m1_det, m2_det, D_L) without spins."""
    a = _write_mixture(tmp_path / "mix.hdf")
    inj = hc.detector_frame_injections(tmp_path / "mix.hdf", 0.25, 10.0)
    found = (a["semianalytic_observed_phase_maximized_snr_net"] > 10) | (a["o3_gstlal_far"] < 0.25)
    f = a[found]
    assert len(inj["prior"]) == found.sum() and inj["n_recorded"] == len(a)
    assert inj["mass_1"] == pytest.approx(f["mass1_source"] * (1 + f["redshift"]))
    a1 = np.sqrt(f["spin1x"] ** 2 + f["spin1y"] ** 2 + f["spin1z"] ** 2)
    a2 = np.sqrt(f["spin2x"] ** 2 + f["spin2y"] ** 2 + f["spin2z"] ** 2)
    expect = (np.exp(f[LNPDRAW]) * (4 * np.pi * a1 ** 2) * (4 * np.pi * a2 ** 2)
              / ((1 + f["redshift"]) ** 2 * 5000.0) / f["weights"])
    assert inj["prior"] == pytest.approx(expect, rel=1e-9)
    assert inj["Tobs"] == pytest.approx(2.0) and inj["ntotal"] == 1e6


def test_prepare_inputs_from_cached_extracts(tmp_path, fake_gwosc):
    """prepare writes inputs.h5 and events.tsv from cached PE extracts, without any download."""
    a = _write_mixture(tmp_path / "mix.hdf")
    cache = tmp_path / "pe"
    (cache / "samples").mkdir(parents=True)
    rng = np.random.default_rng(0)
    for name, lab, prior in (
        ("GW150914_095045", "C01:IMRPhenomXPHM", "PowerLaw(alpha=2, minimum=10, maximum=10000)"),
        ("GW230601_224134", "C00:IMRPhenomXPHM-SpinTaylor",
         "bilby.gw.prior.UniformSourceFrame(minimum=600.0, maximum=9000.0, cosmology='Planck15_LAL')"),
    ):
        with h5py.File(cache / "samples" / f"{name}.h5", "w") as h:
            for l in (lab, "C01:Mixed"):
                g = h.create_group(l)
                g.create_dataset("mass_1", data=rng.uniform(30, 40, 800))
                g.create_dataset("mass_2", data=rng.uniform(20, 30, 800))
                g.create_dataset("luminosity_distance", data=rng.uniform(300, 2000, 800))
                g.attrs["prior:luminosity_distance"] = prior
    ev = hc.prepare_inputs(tmp_path / "work", "gwtc4", tmp_path / "mix.hdf", 0.25, 10.0, 3.0, hc.H0_DEFAULT_EXCLUDE,
                           cache, keep_pe_files=False)
    assert list(ev["pe_label"]) == ["C01:IMRPhenomXPHM", "C00:IMRPhenomXPHM-SpinTaylor"]
    assert list(ev["pe_distance_prior"]) == ["PowerLaw", "UniformSourceFrame"]
    with h5py.File(tmp_path / "work" / "inputs.h5") as h:
        assert sorted(k for k in h if not k.startswith("_")) == ["GW150914_095045", "GW230601_224134"]
        g = h["GW150914_095045"]
        assert g["prior"][:] == pytest.approx(g["dl"][:] ** 2)
        assert h["_injections"].attrs["n_found"] == len(h["_injections"]["prior"]) > 0
    assert (tmp_path / "work" / "events.tsv").exists()


def test_report_from_combined_run(tmp_path):
    """The report stage writes the HTML report and the quantile TSV from summary.json and posterior.tsv."""
    w = tmp_path / "work"
    w.mkdir()
    rng = np.random.default_rng(1)
    post = pd.DataFrame({k: rng.uniform(10, 200, 500) if k == "H0" else rng.normal(size=500) for k in ("H0", "mu_g", "alpha")})
    post.to_csv(w / "posterior.tsv", sep="\t", index=False)
    q = {k: list(np.quantile(post[k], [0.05, 0.16, 0.5, 0.84, 0.95])) for k in post}
    (w / "summary.json").write_text(json.dumps(dict(
        n_runs=2, n_samples=500, log_evidence=-3824.1, log_evidence_err=0.4, run_log_evidences=[-3824.5, -3823.9],
        quantiles=q, settings=dict(nlive=100, pe_samples=1500, inj_fraction=0.1),
        diagnostics=dict(n_points=200, neff_inj_threshold=544, neff_pe_threshold=10, neff_inj_min=3800.0,
                         neff_inj_median=8500.0, neff_pe_min=12.0, neff_pe_median_of_min=27.0,
                         lowest_neff_pe_events=["GW190924_021846"]))))
    table = hc.write_h0_report(w, tmp_path / "h0.html", tmp_path / "h0.tsv", "gwtc4")
    assert list(table["parameter"]) == ["H0", "mu_g", "alpha"]
    html = (tmp_path / "h0.html").read_text()
    assert "data:image/png;base64" in html and "arXiv:2509.04348" in html and "OK." in html
    assert len(pd.read_csv(tmp_path / "h0.tsv", sep="\t")) == 3


def test_runner_command_uses_the_environment_libraries(tmp_path, monkeypatch):
    """The icarogw interpreter runs the standalone sampler with its environment's lib/ on LD_LIBRARY_PATH."""
    (tmp_path / "env" / "bin").mkdir(parents=True)
    (tmp_path / "env" / "lib").mkdir()
    monkeypatch.setenv("LD_LIBRARY_PATH", "/opt/x")
    cmd, env = hc._runner_command(str(tmp_path / "env" / "bin" / "python"))
    assert cmd[1].endswith("h0_icarogw.py")
    assert env["LD_LIBRARY_PATH"] == f"{tmp_path / 'env' / 'lib'}:/opt/x"


def test_cli_hubble_constant_arguments():
    """CLI wiring: defaults reproduce the GWTC-4.0 cosmology setup."""
    a = build_parser().parse_args(["hubble_constant"])
    assert a.mode == "hubble_constant" and a.stages == list(hc.STAGES)
    assert a.sensitivity_release == "gwtc4" and a.far_threshold == 0.25 and a.snr_threshold == 10.0
    assert a.min_mass == 3.0 and a.exclude == list(hc.H0_DEFAULT_EXCLUDE)
    a = build_parser().parse_args(["hubble_constant", "--stages", "combine", "report", "--seeds", "1", "2", "3",
                                   "--icarogw-python", "/env/bin/python"])
    assert a.stages == ["combine", "report"] and a.seeds == [1, 2, 3] and a.icarogw_python == "/env/bin/python"
    assert a.parallel == 1 and build_parser().parse_args(["hubble_constant", "--parallel", "4"]).parallel == 4
    assert a.mass_model == "plp" and build_parser().parse_args(["hubble_constant", "--mass-model", "mltp"]).mass_model == "mltp"


def test_seed_lock(tmp_path, monkeypatch):
    """A seed cannot run twice at once; the lock of a process that died is taken over."""
    import os
    import socket

    from gwtc_analysis import h0_icarogw as runner

    monkeypatch.chdir(tmp_path)
    runner._lock_seed(1)
    lock = tmp_path / "result" / "plp_seed1.lock"
    assert lock.read_text().split() == [socket.gethostname(), str(os.getpid())]
    with pytest.raises(SystemExit, match="already running"):
        runner._lock_seed(1)
    (tmp_path / "result" / "plp_seed2.lock").write_text(f"{socket.gethostname()} 999999999\n")   # stale
    runner._lock_seed(2)
    assert (tmp_path / "result" / "plp_seed2.lock").read_text().split()[1] == str(os.getpid())
    (tmp_path / "result" / "plp_seed3.lock").write_text("farmn7 1234\n")                        # other host
    with pytest.raises(SystemExit, match="farmn7"):
        runner._lock_seed(3)


def test_parallel_runs_report_failures(tmp_path, monkeypatch):
    """Seeds run concurrently with one log each; a failed seed is reported after the others finish."""
    import sys
    import time

    fake = tmp_path / "fake_runner.py"
    fake.write_text(
        "import sys, time\n"
        "seed = int(sys.argv[sys.argv.index('--seed') + 1])\n"
        "time.sleep(1)\n"
        "print(f'[h0] seed {seed}: H0 = 100')\n"
        "sys.exit(1 if seed == 3 else 0)\n")
    monkeypatch.setattr(hc, "_runner_command", lambda python: ([sys.executable, str(fake)], None))
    monkeypatch.setattr(time, "sleep", lambda s: None)      # no 5-second polling in the test
    t0 = time.monotonic()
    with pytest.raises(ValueError, match=r"seed\(s\) \[3\]"):
        hc.run_seeds_parallel(sys.executable, tmp_path, [1, 2, 3, 4], 4, [])
    logs = sorted(p.name for p in (tmp_path / "logs").iterdir())
    assert logs == [f"run_seed{i}.log" for i in (1, 2, 3, 4)]
    assert "H0 = 100" in (tmp_path / "logs" / "run_seed2.log").read_text()
    assert time.monotonic() - t0 < 3.5                      # four 1-second runs at once, not in sequence


def test_mass_model_priors_match_the_paper():
    """PLP and MLTP priors (Tables 3 and 4 of the paper) cover exactly the parameters of each icarogw model."""
    pytest.importorskip("bilby")
    from gwtc_analysis import h0_icarogw as runner

    for model, cfg in runner.MASS_MODELS.items():
        P = runner.priors(model)
        assert set(P) == set(cfg["params"]) | {"Om0"}
    P = runner.priors("mltp")
    assert (P["mu_g_low"].minimum, P["mu_g_low"].maximum) == (5, 100)
    assert (P["sigma_g_low"].minimum, P["sigma_g_low"].maximum) == (0.4, 5)
    assert (P["sigma_g_high"].minimum, P["sigma_g_high"].maximum) == (0.4, 10)
    with pytest.raises(SystemExit, match="unknown mass model"):
        runner.priors("bpl")
    P5 = runner.priors("mltp", "gwtc5")                  # GWTC-5.0 cosmology, Table 5: wider peaks
    assert P5["sigma_g_low"].maximum == 10 and P5["sigma_g_high"].maximum == 15
    assert runner.priors("plp", "gwtc5")["sigma_g"].maximum == 10


def test_report_compares_with_the_published_value_of_its_model(tmp_path):
    """The report quotes the published H0 of the mass model the runs used."""
    w = tmp_path / "work"
    w.mkdir()
    post = pd.DataFrame({"H0": np.random.default_rng(2).uniform(40, 120, 300)})
    post.to_csv(w / "posterior.tsv", sep="\t", index=False)
    (w / "summary.json").write_text(json.dumps(dict(
        mass_model="mltp", n_runs=1, n_samples=300, log_evidence=-1.0, log_evidence_err=0.1, run_log_evidences=[-1.0],
        quantiles={"H0": list(np.quantile(post["H0"], [0.05, 0.16, 0.5, 0.84, 0.95]))}, settings={})))
    hc.write_h0_report(w, tmp_path / "h0.html", None, "gwtc4")
    html = (tmp_path / "h0.html").read_text()
    assert "Multi Peak" in html and "72.3" in html and "(MLTP)" in html


def _probe(t_all=1.0, **fracs):
    """Probe result with fractions given as name=(seconds, predicted ESS fraction, rejected fraction)."""
    return {"seconds_per_eval": {"1.0": t_all, **{k: v[0] for k, v in fracs.items()}},
            "fractions": {k: dict(predicted_ess_fraction=v[1], rejected_fraction=v[2], dlnl_sd=0.5) for k, v in fracs.items()}}


def test_choose_injection_fraction():
    """The smallest subset that is faster, accurate enough and rejects few points; otherwise all the injections."""
    f, why = hc.choose_injection_fraction(_probe(**{"0.1": (0.3, 0.6, 0.0), "0.2": (0.4, 0.8, 0.0), "0.5": (0.7, 0.95, 0.0)}))
    assert f == 0.1 and "reweighted" in why
    # 0.1 too inaccurate, 0.2 rejects too many points: 0.5
    f, _ = hc.choose_injection_fraction(_probe(**{"0.1": (0.3, 0.2, 0.0), "0.2": (0.4, 0.8, 0.2), "0.5": (0.7, 0.9, 0.0)}))
    assert f == 0.5
    # no real speed-up (the PE part dominates): all the injections
    f, why = hc.choose_injection_fraction(_probe(**{"0.1": (0.9, 0.9, 0.0), "0.5": (0.95, 0.99, 0.0)}))
    assert f == 1.0 and "faster" in why
    assert hc.choose_injection_fraction(_probe(**{"0.1": (0.3, 0.6, 0.0)}), min_ess_fraction=0.7)[0] == 1.0


def test_reweighting_math():
    """Importance weights, effective sample size, rejected points and weighted quantiles."""
    from gwtc_analysis import h0_icarogw as runner

    rng = np.random.default_rng(3)
    ll_runs = rng.normal(-3800, 2, 2000)
    r = runner.reweight_weights(ll_runs, ll_runs)                 # same likelihood: uniform weights
    assert r["ess"] == pytest.approx(2000) and r["rejected"] == 0
    target = ll_runs + rng.normal(0, 0.5, 2000)
    target[:10] = -1.797e308                                       # rejected by the target likelihood
    r = runner.reweight_weights(ll_runs, target)
    assert r["rejected"] == 10 and (r["weights"][:10] == 0).all()
    assert r["ess_fraction"] == pytest.approx(np.exp(-0.25), abs=0.08)
    assert runner.predicted_ess_fraction(rng.normal(0, 0.5, 5000)) == pytest.approx(np.exp(-0.25), abs=0.02)
    # weighted quantiles: weights exp(x) on a uniform sample tilt the median upwards
    x = np.linspace(0, 1, 20001)
    assert runner.weighted_quantiles(x, np.ones_like(x))[2] == pytest.approx(0.5, abs=1e-3)
    w = np.exp(x)
    assert runner.weighted_quantiles(x, w)[2] == pytest.approx(np.log((1 + np.e) / 2), abs=1e-3)
    with pytest.raises(SystemExit, match="rejects every"):
        runner.reweight_weights(ll_runs[:3], np.full(3, -1.797e308))


def test_cli_injection_fraction_values():
    """--inj-fraction takes 'auto' (the default) or a fraction in (0, 1]."""
    assert build_parser().parse_args(["hubble_constant"]).inj_fraction == "auto"
    assert build_parser().parse_args(["hubble_constant", "--inj-fraction", "0.2"]).inj_fraction == 0.2
    assert build_parser().parse_args(["hubble_constant", "--inj-fraction", "1"]).inj_fraction == 1.0
    for bad in ("0", "1.5", "all"):
        with pytest.raises(SystemExit):
            build_parser().parse_args(["hubble_constant", "--inj-fraction", bad])


def test_report_leads_with_the_reweighted_result(tmp_path):
    """With a reweighting, the report quotes the reweighted posterior, the runs' result and the ESS."""
    w = tmp_path / "work"
    w.mkdir()
    rng = np.random.default_rng(4)
    pd.DataFrame({"H0": rng.uniform(60, 180, 400)}).to_csv(w / "posterior.tsv", sep="\t", index=False)
    pd.DataFrame({"H0": rng.uniform(50, 170, 400)}).to_csv(w / "posterior_reweighted.tsv", sep="\t", index=False)
    (w / "summary.json").write_text(json.dumps(dict(
        mass_model="plp", n_runs=10, n_samples=400, log_evidence=-3824.1, log_evidence_err=0.4, run_log_evidences=[-3824.1],
        quantiles={"H0": [62.9, 84.4, 119.3, 165.4, 186.1]}, settings=dict(nlive=100, pe_samples=1500, inj_fraction=0.1),
        reweighted=dict(target_inj_fraction=1.0, target_pe_samples=1500, ess=286.0, ess_fraction=0.716, rejected=0,
                        dlnl_mean=-1.38, dlnl_sd=0.75, runs_lnl_match=True,
                        quantiles={"H0": [54.3, 72.5, 106.5, 151.5, 176.9]}))))
    (w / "plan.json").write_text(json.dumps(dict(inj_fraction=0.1, reason="fraction 0.1: 3.5x faster")))
    table = hc.write_h0_report(w, tmp_path / "h0.html", tmp_path / "h0.tsv", "gwtc4")
    assert table.loc[0, "median"] == pytest.approx(106.5)
    html = (tmp_path / "h0.html").read_text()
    assert "106.5" in html and "119.3" in html and "72%" in html and "Reliable reweighting" in html and "3.5x faster" in html


def test_injections_and_events_restricted_to_runs(tmp_path, fake_gwosc):
    """--catalogs restricts the found injections and the events to the observing runs of the catalogs."""
    a = _write_mixture(tmp_path / "mix.hdf")
    all_ = hc.detector_frame_injections(tmp_path / "mix.hdf", 0.25, 10.0)
    o2 = hc.detector_frame_injections(tmp_path / "mix.hdf", 0.25, 10.0, runs=hc.catalog_runs(["GWTC-1"]))
    o3a = hc.detector_frame_injections(tmp_path / "mix.hdf", 0.25, 10.0, runs=hc.catalog_runs(["GWTC-2.1"]))
    assert len(o2["prior"]) + len(o3a["prior"]) == len(all_["prior"]) and len(o2["prior"]) > 0 and len(o3a["prior"]) > 0
    semi = a["time_geocenter"] < 1.2e9
    assert len(o2["prior"]) == int((a["semianalytic_observed_phase_maximized_snr_net"][semi] > 10).sum())
    assert list(hc.select_h0_events("gwtc4", 0.25, 3.0, [], runs=["O4a"])["run"]) == ["O4a", "O4a"]
    assert hc.catalog_runs(["GWTC-5", "GWTC-1"]) == ("O1", "O2", "O4b") and len(hc.catalog_runs(["ALL"])) == 6
    with pytest.raises(ValueError, match="Unknown catalog"):
        hc.catalog_runs(["GWTC-6"])


def test_event_far_threshold_is_inclusive(monkeypatch):
    """Published FARs are rounded: an event listed at exactly the threshold is kept (GW191127, FAR 0.25)."""
    lists = {"GWTC-3-confident": {"GW191127_050227-v1": dict(commonName="GW191127_050227", GPS=1258866165.5, far=0.25,
                                                             mass_1_source=53.0, mass_2_source=24.0),
                                  "GW200216_220804-v1": dict(commonName="GW200216_220804", GPS=1265926102.9, far=0.35,
                                                             mass_1_source=51.0, mass_2_source=30.0)}}
    monkeypatch.setattr(hc.gw, "fetch_gwtc_events", lambda c: {"events": lists.get(c, {})})
    assert list(hc.select_h0_events("gwtc4", 0.25, 3.0, [])["common_name"]) == ["GW191127_050227"]


def test_slurm_chain_auto_fraction(tmp_path):
    """--executor slurm: probe and plan, one array task per seed reading the planned fraction, combine (clearing an
    earlier reweighting), reweighting array and merge that skip when nothing is to be done."""
    (tmp_path / "inputs.h5").write_bytes(b"")
    s = hc.write_h0_slurm_chain(tmp_path, "/env/bin/python", ["sample", "combine", "reweight"], seeds=[1, 2, 3, 2],
                                mass_model="plp", prior_set="gwtc4", nlive=100, npool=8, naccept=60, pe_samples=1500,
                                inj_fraction="auto", probe_points=30, min_ess_fraction=0.5, reweight_pe_samples=None,
                                reweight_jobs=12, slurm_options=["--partition=htc", "--mem=16G"])
    d = tmp_path / "slurm"
    names = sorted(p.name for p in d.glob("[0-9]*.sh"))
    assert names == ["00_probe.sh", "01_plan.sh", "02_run.sh", "03_combine.sh", "04_reweight.sh", "05_reweight_merge.sh"]
    run = (d / "02_run.sh").read_text()
    assert "#SBATCH --array=0-2" in run and "SEEDS=(1 2 3)" in run and "--seed ${SEEDS[$SLURM_ARRAY_TASK_ID]}" in run
    assert "--inj-fraction plan" in run and "#SBATCH --cpus-per-task=8" in run and "#SBATCH --mem=16G" in run
    assert (d / "h0_icarogw.py").exists()
    assert "rm -rf" in (d / "03_combine.sh").read_text()
    rw = (d / "04_reweight.sh").read_text()
    assert "#SBATCH --array=0-11" in rw and "--skip-if-done" in rw and "--chunk $SLURM_ARRAY_TASK_ID" in rw
    assert s.read_text().count("sbatch --parsable") == 6


def test_slurm_chain_fixed_fraction_and_existing_runs(tmp_path):
    (tmp_path / "run_settings.json").write_text(json.dumps(dict(inj_fraction=0.2, pe_samples=1500)))
    hc.write_h0_slurm_chain(tmp_path, "/env/bin/python", ["sample"], seeds=[5], mass_model="mltp", prior_set="gwtc4",
                            nlive=100, npool=4, naccept=60, pe_samples=1500, inj_fraction="auto", probe_points=30,
                            min_ess_fraction=0.5, reweight_pe_samples=None, reweight_jobs=4)
    names = sorted(p.name for p in (tmp_path / "slurm").glob("[0-9]*.sh"))
    assert names == ["00_run.sh"]                                       # the runs' fraction, no new probe
    assert "--inj-fraction 0.2" in (tmp_path / "slurm" / "00_run.sh").read_text()


def test_runner_plan_and_skip(tmp_path):
    from gwtc_analysis import h0_icarogw as runner

    probe = {"seconds_per_eval": {"1.0": 2.0, "0.1": 0.4},
             "fractions": {"0.1": {"predicted_ess_fraction": 0.8, "rejected_fraction": 0.0}}}
    (tmp_path / "probe.json").write_text(json.dumps(probe))
    assert runner.plan(tmp_path) == 0.1
    assert runner._inj_fraction("plan", tmp_path) == 0.1 and runner._inj_fraction("0.3", tmp_path) == 0.3
    assert runner._reweight_needed(tmp_path, 1.0, None)                  # no runs yet
    (tmp_path / "run_settings.json").write_text(json.dumps(dict(inj_fraction=1.0, pe_samples=1500)))
    assert not runner._reweight_needed(tmp_path, 1.0, None) and runner._reweight_needed(tmp_path, 1.0, 3000)
    assert runner.plan(tmp_path) == 1.0                                  # the runs' fraction wins
    assert hc.choose_injection_fraction is runner.choose_injection_fraction


def test_cli_hubble_constant_slurm_arguments():
    a = build_parser().parse_args(["hubble_constant", "--executor", "slurm", "--slurm-option=--mem=16G", "--submit",
                                   "--reweight-jobs", "20", "--galaxy-catalog", "c.hdf5"])
    assert a.executor == "slurm" and a.slurm_option == ["--mem=16G"] and a.submit and a.reweight_jobs == 20

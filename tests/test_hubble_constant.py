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

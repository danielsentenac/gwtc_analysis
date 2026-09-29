"""icarogw sampler of the `hubble_constant` mode (spectral siren; Power Law + Peak or Multi Peak BBH mass model).

This file is self-contained on purpose: icarogw needs its own Python environment (Python >= 3.12,
CPU mode through a `config.py` holding CUPY=False), which usually does not have gwtc_analysis
installed. It only needs numpy, h5py, icarogw and bilby, and reads the `inputs.h5` written by the
`prepare` stage of `gwtc_analysis hubble_constant`.

    python h0_icarogw.py run --workdir DIR --seed 1 [--mass-model plp|mltp --nlive 100 --npool 4 --naccept 60
                                                    --pe-samples 1500 --inj-fraction 0.1]
    python h0_icarogw.py combine --workdir DIR

Independent runs (different seeds, possibly on different machines sharing DIR) are merged by
`combine`, which also writes the posterior as TSV, a corner plot, and a summary with the
numerical-stability diagnostics (effective numbers of injections and PE samples).
"""
from __future__ import annotations

import argparse
import atexit
import json
import os
import socket
import sys
from pathlib import Path

import numpy as np

# BBH primary-mass models of the GWTC-4.0 cosmology paper (arXiv:2509.04348): the icarogw mass prior and
# the population parameters, in the order of the posterior tables
MASS_MODELS = {
    "plp": dict(name="Power Law + Peak", prior="massprior_PowerLawPeak",
                params=("H0", "alpha", "beta", "mmin", "mmax", "delta_m", "mu_g", "sigma_g", "lambda_peak",
                        "gamma", "kappa", "zp")),
    "mltp": dict(name="Multi Peak", prior="massprior_MultiPeak",
                 params=("H0", "alpha", "beta", "mmin", "mmax", "delta_m", "mu_g_low", "sigma_g_low", "mu_g_high",
                         "sigma_g_high", "lambda_g", "lambda_g_low", "gamma", "kappa", "zp")),
}
OM0 = 0.3065
NEFF_PE = 10


def _enter(workdir: Path) -> None:
    """icarogw reads `config` from the import path: CPU mode (CUPY=False) from the work directory."""
    workdir.mkdir(parents=True, exist_ok=True)
    os.chdir(workdir)
    Path("config.py").write_text("CUPY=False\n")
    sys.path.insert(0, str(workdir))


def build_likelihood(pe_samples: int, inj_fraction: float, mass_model: str = "plp"):
    """Hierarchical likelihood on the events and found injections of inputs.h5 (cwd)."""
    import h5py
    import icarogw

    with h5py.File("inputs.h5", "r") as h:
        pes = {}
        for name in h:
            if name.startswith("_"):
                continue
            g = h[name]
            n = min(pe_samples, g["dl"].shape[0])
            pes[name] = icarogw.posterior_samples.posterior_samples(
                {"mass_1": g["m1"][:n], "mass_2": g["m2"][:n], "luminosity_distance": g["dl"][:n]},
                prior=g["prior"][:n])
        gi = h["_injections"]
        # random subset of the found injections: unbiased when ntotal is scaled by the kept fraction
        keep = np.random.default_rng(2024).random(gi["prior"].shape[0]) < inj_fraction
        inj = icarogw.injections.injections(
            {k: gi[k][:][keep] for k in ("mass_1", "mass_2", "luminosity_distance")},
            prior=gi["prior"][:][keep], ntotal=float(gi.attrs["ntotal"]) * keep.mean(), Tobs=float(gi.attrs["Tobs"]))
    print(f"[h0] {len(pes)} events, {pe_samples} PE samples each at most; {int(keep.sum())} injections "
          f"(fraction {keep.mean():.3f})", flush=True)
    cat = icarogw.posterior_samples.posterior_samples_catalog(pes)
    rate = icarogw.rates.CBC_vanilla_rate(
        icarogw.wrappers.FlatLambdaCDM_wrap(zmax=20.0),
        icarogw.wrappers.m1m2_conditioned_lowpass(getattr(icarogw.wrappers, MASS_MODELS[mass_model]["prior"])()),
        icarogw.wrappers.rateevolution_Madau(), scale_free=True)
    like = icarogw.likelihood.hierarchical_likelihood(cat, inj, rate, nparallel=pe_samples, neffPE=NEFF_PE,
                                                      neffINJ=None)
    return like, rate, cat, inj


def priors(mass_model: str = "plp"):
    """GWTC-4.0 cosmology paper (arXiv:2509.04348): Tables 3 (PLP), 4 (MLTP) and 6 (Madau-Dickinson)."""
    import bilby

    U = bilby.core.prior.Uniform
    P = bilby.core.prior.PriorDict()
    P["H0"] = U(10, 200, "H0")
    P["Om0"] = OM0
    P["alpha"], P["beta"] = U(1.5, 12, "alpha"), U(-4, 12, "beta")
    P["mmin"], P["mmax"], P["delta_m"] = U(2, 10, "mmin"), U(50, 200, "mmax"), U(1e-3, 10, "delta_m")
    if mass_model == "plp":
        P["mu_g"], P["sigma_g"], P["lambda_peak"] = U(20, 50, "mu_g"), U(0.4, 10, "sigma_g"), U(0, 1, "lambda_peak")
    elif mass_model == "mltp":
        P["mu_g_low"], P["sigma_g_low"] = U(5, 100, "mu_g_low"), U(0.4, 5, "sigma_g_low")
        P["mu_g_high"], P["sigma_g_high"] = U(5, 100, "mu_g_high"), U(0.4, 10, "sigma_g_high")
        P["lambda_g"], P["lambda_g_low"] = U(0, 1, "lambda_g"), U(0, 1, "lambda_g_low")
    else:
        raise SystemExit(f"[h0] unknown mass model {mass_model!r}; choose from {', '.join(MASS_MODELS)}")
    P["gamma"], P["kappa"], P["zp"] = U(0, 12, "gamma"), U(0, 6, "kappa"), U(0, 4, "zp")
    return P


def _settings(workdir: Path) -> dict:
    p = workdir / "run_settings.json"
    s = json.loads(p.read_text()) if p.exists() else {}
    if s:
        s.setdefault("mass_model", "plp")        # work directories made before the model choice existed
    return s


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _lock_seed(seed: int, mass_model: str = "plp") -> None:
    """One process per seed: result/<model>_seed<N>.lock holds 'host pid' while the run is going (cwd = workdir)."""
    lock = Path("result") / f"{mass_model}_seed{seed}.lock"
    lock.parent.mkdir(exist_ok=True)
    host = socket.gethostname()
    for _ in range(2):
        try:
            fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            break
        except FileExistsError:
            parts = lock.read_text().split()
            h, pid = (parts + ["?", "-1"])[:2]
            if h == host and not _pid_alive(int(pid)):
                lock.unlink(missing_ok=True)       # stale: the process is gone
                continue
            raise SystemExit(f"[h0] seed {seed} is already running on {h} (pid {pid}); "
                             f"if it is not, remove {lock.resolve()}")
    else:
        raise SystemExit(f"[h0] cannot take the lock {lock.resolve()}")
    os.write(fd, f"{host} {os.getpid()}\n".encode())
    os.close(fd)
    atexit.register(lambda: lock.unlink(missing_ok=True))


def run(workdir: Path, seed: int, nlive: int, npool: int, naccept: int, pe_samples: int, inj_fraction: float,
        mass_model: str = "plp") -> Path:
    """One dynesty run; resumable from its checkpoint."""
    if mass_model not in MASS_MODELS:
        raise SystemExit(f"[h0] unknown mass model {mass_model!r}; choose from {', '.join(MASS_MODELS)}")
    _enter(workdir)
    _lock_seed(seed, mass_model)
    import bilby

    s = dict(nlive=nlive, pe_samples=pe_samples, inj_fraction=inj_fraction, mass_model=mass_model)
    old = _settings(workdir)
    if old and old != s:
        raise SystemExit(f"[h0] {workdir} holds runs made with {old}; runs with {s} cannot be combined with them. "
                         "Use another work directory.")
    (workdir / "run_settings.json").write_text(json.dumps(s))
    like, _, _, _ = build_likelihood(pe_samples, inj_fraction, mass_model)
    res = bilby.run_sampler(like, priors(mass_model), sampler="dynesty", nlive=nlive, npool=npool, outdir="result",
                            label=f"{mass_model}_seed{seed}", sample="acceptance-walk", naccept=naccept, seed=seed,
                            resume=True, check_point_delta_t=600)
    q = np.quantile(res.posterior["H0"], [0.05, 0.16, 0.5, 0.84, 0.95])
    print(f"[h0] seed {seed}: H0 = {q[2]:.1f} (+{q[3] - q[2]:.1f} / -{q[2] - q[1]:.1f}) km/s/Mpc [68%], "
          f"90%: {q[0]:.1f}-{q[4]:.1f}; ln Z = {res.log_evidence:.2f}", flush=True)
    return workdir / "result" / f"{mass_model}_seed{seed}_result.json"


def _diagnostics(post, n_points: int = 200) -> dict:
    """Effective numbers of injections and PE samples over posterior draws (icarogw's stability criteria)."""
    s = _settings(Path.cwd())
    like, rate, cat, inj = build_likelihood(int(s.get("pe_samples", 1500)), float(s.get("inj_fraction", 1.0)),
                                            s.get("mass_model", "plp"))
    rows = post.sample(min(n_points, len(post)), random_state=1)
    neff_inj, neff_pe, worst = [], [], {}
    names = list(cat.posterior_samples_dict.keys()) if hasattr(cat, "posterior_samples_dict") else None
    for _, r in rows.iterrows():
        rate.update(**{k: float(r[k]) if k in r else OM0 for k in rate.population_parameters})
        inj.update_weights(rate)
        neff_inj.append(float(inj.effective_injections_number()))
        cat.update_weights(rate)
        e = np.asarray(cat.get_effective_number_of_PE(), dtype=float)
        neff_pe.append(float(e.min()))
        if names is not None:
            k = names[int(e.argmin())]
            worst[k] = worst.get(k, 0) + 1
    ni, npe = np.array(neff_inj), np.array(neff_pe)
    return dict(n_points=len(rows), neff_inj_threshold=int(like.neffINJ), neff_pe_threshold=NEFF_PE,
                neff_inj_min=float(ni.min()), neff_inj_median=float(np.median(ni)),
                neff_pe_min=float(npe.min()), neff_pe_median_of_min=float(np.median(npe)),
                lowest_neff_pe_events=sorted(worst, key=worst.get, reverse=True)[:5])


def combine(workdir: Path, diagnostics: bool = True) -> dict:
    """Merge the per-seed runs, write posterior.tsv, corner.png and summary.json."""
    _enter(workdir)
    import bilby

    model = _settings(Path.cwd()).get("mass_model", "plp")
    params = MASS_MODELS[model]["params"]
    files = sorted(Path("result").glob(f"{model}_seed*_result.json"))
    if not files:
        raise SystemExit(f"[h0] no finished run in {workdir / 'result'}")
    results = [bilby.core.result.read_in_result(str(f)) for f in files]
    res = bilby.core.result.ResultList(results).combine() if len(results) > 1 else results[0]
    res.label, res.outdir = f"{model}_combined", "result"
    res.save_to_file(overwrite=True, extension="json")
    post = res.posterior
    post[[k for k in params if k in post]].to_csv("posterior.tsv", sep="\t", index=False, float_format="%.6g")
    res.plot_corner(parameters=[k for k in params if k in post], filename="corner.png", quantiles=[0.05, 0.95])
    quant = {}
    for k in params:
        if k in post:
            quant[k] = [float(x) for x in np.quantile(post[k], [0.05, 0.16, 0.5, 0.84, 0.95])]
    summary = dict(mass_model=model, n_runs=len(results), runs=[f.name for f in files], n_samples=len(post),
                   log_evidence=float(res.log_evidence), log_evidence_err=float(res.log_evidence_err),
                   run_log_evidences=[float(r.log_evidence) for r in results], quantiles=quant,
                   settings=_settings(Path.cwd()))
    if diagnostics:
        summary["diagnostics"] = _diagnostics(post)
    Path("summary.json").write_text(json.dumps(summary, indent=1))
    q = quant["H0"]
    print(f"[h0] {len(results)} run(s), {len(post)} samples: H0 = {q[2]:.1f} (+{q[3] - q[2]:.1f} / -{q[2] - q[1]:.1f}) "
          f"km/s/Mpc [68%], 90%: {q[0]:.1f}-{q[4]:.1f}", flush=True)
    return summary


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="icarogw spectral-siren sampler of gwtc_analysis hubble_constant.")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--workdir", required=True)
    r.add_argument("--seed", type=int, default=1)
    r.add_argument("--mass-model", choices=list(MASS_MODELS), default="plp")
    r.add_argument("--nlive", type=int, default=100)
    r.add_argument("--npool", type=int, default=4)
    r.add_argument("--naccept", type=int, default=60)
    r.add_argument("--pe-samples", type=int, default=1500)
    r.add_argument("--inj-fraction", type=float, default=0.1)
    c = sub.add_parser("combine")
    c.add_argument("--workdir", required=True)
    c.add_argument("--no-diagnostics", action="store_true")
    a = p.parse_args(argv)
    wd = Path(a.workdir).expanduser().resolve()
    if a.cmd == "run":
        run(wd, a.seed, a.nlive, a.npool, a.naccept, a.pe_samples, a.inj_fraction, a.mass_model)
    else:
        combine(wd, diagnostics=not a.no_diagnostics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""icarogw sampler of the `hubble_constant` mode (spectral siren; Power Law + Peak or Multi Peak BBH mass model).

This file is self-contained on purpose: icarogw needs its own Python environment (Python >= 3.12,
CPU mode through a `config.py` holding CUPY=False), which usually does not have gwtc_analysis
installed. It only needs numpy, h5py, icarogw and bilby, and reads the `inputs.h5` written by the
`prepare` stage of `gwtc_analysis hubble_constant`.

    python h0_icarogw.py run --workdir DIR --seed 1 [--mass-model plp|mltp --nlive 100 --npool 4 --naccept 60
                                                    --pe-samples 1500 --inj-fraction 0.1]
    python h0_icarogw.py combine --workdir DIR
    python h0_icarogw.py probe --workdir DIR [--mass-model plp --pe-samples 1500 --fractions 0.1 0.2 0.5]
    python h0_icarogw.py reweight --workdir DIR --chunk I --nchunks N [--target-inj-fraction 1]
    python h0_icarogw.py reweight-merge --workdir DIR

Independent runs (different seeds, possibly on different machines sharing DIR) are merged by
`combine`, which also writes the posterior as TSV, a corner plot, and a summary with the
numerical-stability diagnostics (effective numbers of injections and PE samples).

`probe` measures, before sampling, the speed and the accuracy of likelihoods built on subsets of the
found injections, to choose the fastest strategy. `reweight` and `reweight-merge` turn a posterior
sampled with a subset into the posterior with all the injections, by importance reweighting of its
samples (weights exp(ln L_target - ln L_runs)), with the effective sample size as the validity check.
"""
from __future__ import annotations

import argparse
import atexit
import json
import os
import socket
import sys
import time
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


def build_likelihood(pe_samples: int, inj_fraction: float, mass_model: str = "plp", inputs: str = "inputs.h5"):
    """Hierarchical likelihood on the events and found injections of `inputs` (default: inputs.h5 of the cwd)."""
    import h5py
    import icarogw

    with h5py.File(inputs, "r") as h:
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


PRIOR_SETS = ("gwtc4", "gwtc5")


def priors(mass_model: str = "plp", prior_set: str = "gwtc4"):
    """Priors of the LVK cosmology papers: GWTC-4.0 (arXiv:2509.04348, Tables 3, 4 and 6) or GWTC-5.0
    (arXiv:2605.27227, Tables 5 and 8), which widens the MLTP peak widths."""
    if prior_set not in PRIOR_SETS:
        raise SystemExit(f"[h0] unknown prior set {prior_set!r}; choose from {', '.join(PRIOR_SETS)}")
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
        wide = prior_set == "gwtc5"
        P["mu_g_low"], P["sigma_g_low"] = U(5, 100, "mu_g_low"), U(0.4, 10 if wide else 5, "sigma_g_low")
        P["mu_g_high"], P["sigma_g_high"] = U(5, 100, "mu_g_high"), U(0.4, 15 if wide else 10, "sigma_g_high")
        P["lambda_g"], P["lambda_g_low"] = U(0, 1, "lambda_g"), U(0, 1, "lambda_g_low")
    else:
        raise SystemExit(f"[h0] unknown mass model {mass_model!r}; choose from {', '.join(MASS_MODELS)}")
    P["gamma"], P["kappa"], P["zp"] = U(0, 12, "gamma"), U(0, 6, "kappa"), U(0, 4, "zp")
    return P


def _settings(workdir: Path) -> dict:
    p = workdir / "run_settings.json"
    s = json.loads(p.read_text()) if p.exists() else {}
    if s:
        s.setdefault("mass_model", "plp")        # work directories made before these choices existed
        s.setdefault("prior_set", "gwtc4")
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
        mass_model: str = "plp", prior_set: str = "gwtc4") -> Path:
    """One dynesty run; resumable from its checkpoint."""
    if mass_model not in MASS_MODELS:
        raise SystemExit(f"[h0] unknown mass model {mass_model!r}; choose from {', '.join(MASS_MODELS)}")
    _enter(workdir)
    _lock_seed(seed, mass_model)
    import bilby

    s = dict(nlive=nlive, pe_samples=pe_samples, inj_fraction=inj_fraction, mass_model=mass_model, prior_set=prior_set)
    old = _settings(workdir)
    if old and old != s:
        raise SystemExit(f"[h0] {workdir} holds runs made with {old}; runs with {s} cannot be combined with them. "
                         "Use another work directory.")
    (workdir / "run_settings.json").write_text(json.dumps(s))
    like, _, _, _ = build_likelihood(pe_samples, inj_fraction, mass_model)
    res = bilby.run_sampler(like, priors(mass_model, prior_set), sampler="dynesty", nlive=nlive, npool=npool, outdir="result",
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


# ---------------------------------------------------------------------------
# importance reweighting (numpy only, so that it can be tested without icarogw)
# ---------------------------------------------------------------------------

REJECTED = -1e300     # icarogw returns nan_to_num(-inf) for the points it rejects


def reweight_weights(ll_runs: np.ndarray, ll_target: np.ndarray) -> dict:
    """Weights exp(ln L_target - ln L_runs), normalized, with their effective sample size.

    Points rejected by the target likelihood get zero weight; points rejected by the runs cannot occur
    (they are not in the posterior)."""
    ll_runs, ll_target = np.asarray(ll_runs, float), np.asarray(ll_target, float)
    ok = (ll_target > REJECTED) & np.isfinite(ll_target)
    d = np.where(ok, ll_target - ll_runs, -np.inf)
    if not ok.any():
        raise SystemExit("[h0] the target likelihood rejects every posterior sample: reweighting impossible")
    w = np.exp(d - d[ok].max())
    w /= w.sum()
    ess = float(1.0 / np.sum(w ** 2))
    return dict(weights=w, ess=ess, ess_fraction=ess / len(w), rejected=int((~ok).sum()),
                dlnl_mean=float(d[ok].mean()), dlnl_sd=float(d[ok].std()))


def weighted_quantiles(x: np.ndarray, w: np.ndarray, q=(0.05, 0.16, 0.5, 0.84, 0.95)) -> list:
    o = np.argsort(x)
    c = np.cumsum(w[o]) - 0.5 * w[o]
    c /= w.sum()
    return [float(v) for v in np.interp(q, c, np.asarray(x)[o])]


def predicted_ess_fraction(dlnl: np.ndarray) -> float:
    """ESS fraction of log-normal weights with the scatter of dlnl: exp(-sigma^2)."""
    return float(np.exp(-np.var(np.asarray(dlnl, float))))


# ---------------------------------------------------------------------------
# probe, reweight
# ---------------------------------------------------------------------------

_PILOT = {}   # likelihood and prior of the pilot chain, shared with the forked worker processes


def _pilot_log_prob(x):
    keys, P, like = _PILOT["keys"], _PILOT["prior"], _PILOT["like"]
    point = dict(zip(keys, (float(v) for v in x)))
    lp = P.ln_prob(point)
    if not np.isfinite(lp):
        return -np.inf
    like.parameters.update(point | {"Om0": OM0})
    ll = like.log_likelihood()
    return ll + lp if ll > REJECTED else -np.inf


def probe(workdir: Path, mass_model: str, pe_samples: int, fractions, npoints: int = 30,
          max_draws: int = 4000, seed: int = 0, npool: int = 1, pilot_steps: int = 100,
          measure_points: int = 60, prior_set: str = "gwtc4") -> dict:
    """Speed and accuracy of likelihoods on injection subsets, before sampling.

    1. prior points with a finite likelihood (all the injections): the fraction of them that each subset
       rejects;
    2. a short ensemble MCMC (emcee), started from those points, with the likelihood of the smallest
       subset: its walkers move to the region the sampler runs will explore;
    3. at the final walker positions: seconds per evaluation of each likelihood, and the scatter of
       ln L_f - ln L_all, which predicts the ESS fraction of a reweighting to all the injections."""
    import multiprocessing

    import emcee

    _enter(workdir)
    fractions = sorted({float(f) for f in fractions if 0 < float(f) < 1})
    if not fractions:
        raise SystemExit("[h0] probe: no subset fraction in (0, 1) to test")
    np.random.seed(seed)
    P = priors(mass_model, prior_set)
    keys = [k for k in MASS_MODELS[mass_model]["params"]]
    full, _, _, _ = build_likelihood(pe_samples, 1.0, mass_model)
    subs = {f: build_likelihood(pe_samples, f, mass_model)[0] for f in fractions}
    # 1. finite prior points
    pts, draws = [], 0
    while len(pts) < npoints and draws < max_draws:
        draws += 1
        smp = P.sample()
        full.parameters.update(smp)
        if full.log_likelihood() > REJECTED:
            pts.append({k: float(smp[k]) for k in keys})
    if len(pts) < 4:
        raise SystemExit(f"[h0] probe: only {len(pts)} finite points in {draws} prior draws")
    rejected = {}
    for f, like in subs.items():
        n = 0
        for pt in pts:
            like.parameters.update(pt | {"Om0": OM0})
            n += like.log_likelihood() <= REJECTED
        rejected[f] = n / len(pts)
    # 2. pilot ensemble MCMC with the smallest subset
    ndim = len(keys)
    nwalkers = max(32, 2 * ndim + 2)
    rng = np.random.default_rng(seed)
    lo = np.array([P[k].minimum for k in keys]); hi = np.array([P[k].maximum for k in keys])
    start = np.array([[pt[k] for k in keys] for pt in (pts[i] for i in rng.integers(0, len(pts), nwalkers))])
    start = np.clip(start + 1e-4 * (hi - lo) * rng.standard_normal(start.shape), lo + 1e-9 * (hi - lo), hi - 1e-9 * (hi - lo))
    _PILOT.update(keys=keys, prior=P, like=subs[fractions[0]])
    t0 = time.time()
    if npool > 1:
        with multiprocessing.get_context("fork").Pool(npool) as pool:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, _pilot_log_prob, pool=pool)
            sampler.run_mcmc(start, pilot_steps, progress=False)
    else:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, _pilot_log_prob)
        sampler.run_mcmc(start, pilot_steps, progress=False)
    chain, lnp = sampler.get_chain(), sampler.get_log_prob()
    last = max(1, int(np.ceil(measure_points / nwalkers)))
    cand = chain[-last:].reshape(-1, ndim)[np.isfinite(lnp[-last:].reshape(-1))][:measure_points]
    print(f"[h0] probe: pilot MCMC of {nwalkers} walkers x {pilot_steps} steps with fraction {fractions[0]:g} "
          f"({time.time() - t0:.0f} s); ln posterior median {np.median(lnp[0][np.isfinite(lnp[0])]):.1f} -> "
          f"{np.median(lnp[-1][np.isfinite(lnp[-1])]):.1f}", flush=True)
    # 3. speed and accuracy at the pilot positions
    def evaluate(like):
        ll, tt = [], []
        for x in cand:
            like.parameters.update(dict(zip(keys, map(float, x))) | {"Om0": OM0})
            t1 = time.perf_counter(); ll.append(like.log_likelihood()); tt.append(time.perf_counter() - t1)
        return np.array(ll), float(np.median(tt))

    ll_full, t_full = evaluate(full)
    out = dict(mass_model=mass_model, prior_set=prior_set, pe_samples=pe_samples, draws=draws, points=len(pts),
               finite_fraction=len(pts) / draws, pilot=dict(walkers=nwalkers, steps=pilot_steps, fraction=fractions[0],
                                                             measured_points=int(len(cand))),
               seconds_per_eval={"1.0": t_full}, fractions={})
    for f, like in subs.items():
        ll, tf = evaluate(like)
        ok = (ll > REJECTED) & (ll_full > REJECTED)
        d = (ll - ll_full)[ok]
        out["seconds_per_eval"][str(f)] = tf
        out["fractions"][str(f)] = dict(
            dlnl_sd=float(d.std()) if len(d) > 1 else float("inf"),
            predicted_ess_fraction=predicted_ess_fraction(d) if len(d) > 1 else 0.0,
            rejected_fraction=float(max(rejected[f], 1 - ok.mean())))
        print(f"[h0] probe f={f:g}: {tf:.3f} s/eval, sd(dlnL) = {out['fractions'][str(f)]['dlnl_sd']:.2f}, predicted "
              f"ESS fraction {out['fractions'][str(f)]['predicted_ess_fraction']:.2f}, rejected "
              f"{out['fractions'][str(f)]['rejected_fraction']:.0%}", flush=True)
    print(f"[h0] probe f=1: {t_full:.3f} s/eval; {len(pts)} finite prior points in {draws} draws", flush=True)
    Path("probe.json").write_text(json.dumps(out, indent=1))
    return out


def _combined_posterior(workdir: Path):
    import pandas as pd

    model = _settings(workdir).get("mass_model", "plp")
    res = json.loads((workdir / "result" / f"{model}_combined_result.json").read_text())
    return model, pd.DataFrame({k: np.asarray(v) for k, v in res["posterior"]["content"].items()
                                if isinstance(v, list) and len(v) and not isinstance(v[0], (dict, str))})


def reweight(workdir: Path, chunk: int, nchunks: int, target_inj_fraction: float = 1.0,
             target_pe_samples: int | None = None, target_inputs: str | None = None) -> Path:
    """ln L of the runs and of the target settings at a chunk of the combined posterior samples.

    `target_inputs`: another inputs.h5 for the target likelihood (e.g. with one more event)."""
    _enter(workdir)
    s = _settings(workdir)
    model, post = _combined_posterior(workdir)
    params = [k for k in MASS_MODELS[model]["params"]]
    idx = np.array_split(np.arange(len(post)), nchunks)[chunk]
    npe = int(target_pe_samples or s["pe_samples"])
    out = {"idx": idx, "stored": post["log_likelihood"].to_numpy()[idx]}
    target_file = str(Path(target_inputs).expanduser().resolve()) if target_inputs else "inputs.h5"
    for key, (frac, pe, inputs) in (("runs", (float(s["inj_fraction"]), int(s["pe_samples"]), "inputs.h5")),
                                    ("target", (float(target_inj_fraction), npe, target_file))):
        like, _, _, _ = build_likelihood(pe, frac, model, inputs)
        ll = np.empty(len(idx))
        for j, i in enumerate(idx):
            like.parameters.update({k: float(post[k].iloc[i]) for k in params} | {"Om0": OM0})
            ll[j] = like.log_likelihood()
        out[key] = ll
        del like
    Path("reweight").mkdir(exist_ok=True)
    dest = Path("reweight") / f"chunk{chunk:03d}_of_{nchunks:03d}.npz"
    np.savez(dest, target_settings=np.array([target_inj_fraction, npe]), **out)
    print(f"[h0] reweight chunk {chunk + 1}/{nchunks}: {len(idx)} samples", flush=True)
    return dest


def reweight_merge(workdir: Path, seed: int = 1) -> dict:
    """Merge the chunks: weights, ESS, reweighted quantiles and a resampled posterior_reweighted.tsv."""
    _enter(workdir)
    model, post = _combined_posterior(workdir)
    files = sorted(Path("reweight").glob("chunk*_of_*.npz"))
    if not files:
        raise SystemExit("[h0] no reweighting chunk found")
    nchunks = int(files[0].name.split("_of_")[1][:3])
    if len(files) != nchunks:
        raise SystemExit(f"[h0] {len(files)} of {nchunks} reweighting chunks present")
    d = [np.load(f) for f in files]
    idx = np.concatenate([x["idx"] for x in d])
    o = np.argsort(idx)
    ll_runs, ll_target = (np.concatenate([x[k] for x in d])[o] for k in ("runs", "target"))
    stored = np.concatenate([x["stored"] for x in d])[o]
    target = [float(v) for v in d[0]["target_settings"]]
    r = reweight_weights(ll_runs, ll_target)
    params = [k for k in MASS_MODELS[model]["params"] if k in post]
    quant = {k: weighted_quantiles(post[k].to_numpy(), r["weights"]) for k in params}
    rng = np.random.default_rng(seed)
    pick = rng.choice(len(post), size=len(post), replace=True, p=r["weights"])
    post[params].iloc[pick].to_csv("posterior_reweighted.tsv", sep="\t", index=False, float_format="%.6g")
    summary = json.loads(Path("summary.json").read_text()) if Path("summary.json").exists() else {}
    summary["reweighted"] = dict(target_inj_fraction=target[0], target_pe_samples=int(target[1]), ess=r["ess"],
                                 ess_fraction=r["ess_fraction"], rejected=r["rejected"], dlnl_mean=r["dlnl_mean"],
                                 dlnl_sd=r["dlnl_sd"], runs_lnl_match=bool(np.allclose(ll_runs, stored)),
                                 quantiles=quant)
    Path("summary.json").write_text(json.dumps(summary, indent=1))
    q = quant["H0"]
    print(f"[h0] reweighted to injection fraction {target[0]:g}, {int(target[1])} PE samples: H0 = {q[2]:.1f} "
          f"(+{q[3] - q[2]:.1f} / -{q[2] - q[1]:.1f}), 90%: {q[0]:.1f}-{q[4]:.1f}; ESS {r['ess']:.0f} of {len(post)} "
          f"({r['ess_fraction']:.0%}), {r['rejected']} rejected; runs ln L reproduced: {summary['reweighted']['runs_lnl_match']}",
          flush=True)
    if not summary["reweighted"]["runs_lnl_match"]:
        print("[h0] WARN: ln L of the runs not reproduced; the inputs or settings changed since sampling", flush=True)
    return summary["reweighted"]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="icarogw spectral-siren sampler of gwtc_analysis hubble_constant.")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--workdir", required=True)
    r.add_argument("--seed", type=int, default=1)
    r.add_argument("--mass-model", choices=list(MASS_MODELS), default="plp")
    r.add_argument("--prior-set", choices=list(PRIOR_SETS), default="gwtc4")
    r.add_argument("--nlive", type=int, default=100)
    r.add_argument("--npool", type=int, default=4)
    r.add_argument("--naccept", type=int, default=60)
    r.add_argument("--pe-samples", type=int, default=1500)
    r.add_argument("--inj-fraction", type=float, default=0.1)
    c = sub.add_parser("combine")
    c.add_argument("--workdir", required=True)
    c.add_argument("--no-diagnostics", action="store_true")
    pr = sub.add_parser("probe")
    pr.add_argument("--workdir", required=True)
    pr.add_argument("--mass-model", choices=list(MASS_MODELS), default="plp")
    pr.add_argument("--prior-set", choices=list(PRIOR_SETS), default="gwtc4")
    pr.add_argument("--pe-samples", type=int, default=1500)
    pr.add_argument("--fractions", type=float, nargs="+", default=[0.1, 0.2, 0.5])
    pr.add_argument("--npoints", type=int, default=30)
    pr.add_argument("--npool", type=int, default=1)
    pr.add_argument("--pilot-steps", type=int, default=100)
    rw = sub.add_parser("reweight")
    rw.add_argument("--workdir", required=True)
    rw.add_argument("--chunk", type=int, required=True)
    rw.add_argument("--nchunks", type=int, required=True)
    rw.add_argument("--target-inj-fraction", type=float, default=1.0)
    rw.add_argument("--target-pe-samples", type=int, default=None)
    rw.add_argument("--target-inputs", default=None, help="inputs.h5 of the target likelihood (default: the runs')")
    rm = sub.add_parser("reweight-merge")
    rm.add_argument("--workdir", required=True)
    a = p.parse_args(argv)
    wd = Path(a.workdir).expanduser().resolve()
    if a.cmd == "run":
        run(wd, a.seed, a.nlive, a.npool, a.naccept, a.pe_samples, a.inj_fraction, a.mass_model, a.prior_set)
    elif a.cmd == "combine":
        combine(wd, diagnostics=not a.no_diagnostics)
    elif a.cmd == "probe":
        probe(wd, a.mass_model, a.pe_samples, a.fractions, a.npoints, npool=a.npool, pilot_steps=a.pilot_steps,
              prior_set=a.prior_set)
    elif a.cmd == "reweight":
        reweight(wd, a.chunk, a.nchunks, a.target_inj_fraction, a.target_pe_samples, a.target_inputs)
    else:
        reweight_merge(wd)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

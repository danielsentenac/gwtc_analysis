"""Hubble constant from the BBH mass spectrum (spectral siren), with icarogw.

Stages of `run_hubble_constant`:

- ``prepare``: BBH events from the GWOSC lists (lowest FAR below a threshold, inside the observing
  runs of the sensitivity release, both masses above a minimum); their PE samples downloaded from
  Zenodo and reduced to (m1, m2, D_L) in the detector frame, with the density of the PE prior they
  were drawn from; the found LVK injections with their draw density in the same variables. All of it
  goes to ``<workdir>/inputs.h5``.
- ``sample``: independent icarogw + bilby/dynesty runs, one per seed (`h0_icarogw.py`, run by the
  icarogw interpreter).
- ``combine``: the runs merged, with numerical-stability diagnostics.
- ``report``: HTML report.

The default setup reproduces the Power Law + Peak measurement of the GWTC-4.0 cosmology paper
(arXiv:2509.04348): H0 = 112.7 (+51.0 / -35.9) km/s/Mpc.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from . import gw_stat as gw
from .report import write_simple_html_report

# GWOSC observing-run boundaries (GPS)
OBSERVING_RUNS = {
    "O1": (1126051217, 1137254417), "O2": (1164556817, 1187733618),
    "O3a": (1238166018, 1253977218), "O3b": (1256655618, 1269363618),
    "O4a": (1368975618, 1389456018), "O4b": (1396796418, 1422118818),
}

# Cumulative search-sensitivity releases with semi-analytic O1+O2 injections, so that the whole
# catalog since O1 can be used.
H0_SENSITIVITY_RELEASES = {
    "gwtc4": dict(
        record="16740128", label="GWTC-4.0 cumulative, semi-analytic O1+O2 + real O3+O4a injections",
        file_re=r"^mixture-semi_o1_o2-real_o3_o4a-cartesian_spins.*\.hdf5?$",
        runs=("O1", "O2", "O3a", "O3b", "O4a"),
        catalogs=("GWTC-1-confident", "GWTC-2.1-confident", "GWTC-2.1-marginal", "GWTC-3-confident",
                  "GWTC-3-marginal", "GWTC-4.0"),
        published=dict(ref="GWTC-4.0 cosmology, arXiv:2509.04348 (PLP)", median=112.7, plus=51.0, minus=35.9,
                       lo90=57.6, hi90=186.7),
    ),
    "gwtc5": dict(
        record="19500052", label="GWTC-5.0 cumulative, semi-analytic O1+O2 + real O3+O4a+O4b injections",
        file_re=r"^mixture-semi_o1_o2-real_o3_o4a_o4b-cartesian_spins.*\.hdf5?$",
        runs=("O1", "O2", "O3a", "O3b", "O4a", "O4b"),
        catalogs=("GWTC-1-confident", "GWTC-2.1-confident", "GWTC-2.1-marginal", "GWTC-3-confident",
                  "GWTC-3-marginal", "GWTC-4.0", "GWTC-5.0"),
        published=None,
    ),
}
H0_DEFAULT_RELEASE = "gwtc4"
# GWTC-4.0 cosmology: GW231123 (PE systematics) and GW200105 (NSBH) are not used
H0_DEFAULT_EXCLUDE = ("GW231123_135430", "GW200105_162426")
# PE labels, in order of preference (the GWTC-4.0 cosmology choice first)
PE_LABELS = ("C00:IMRPhenomXPHM-SpinTaylor", "C01:IMRPhenomXPHM", "C00:IMRPhenomXPHM")
_PE_COLUMNS = ("mass_1", "mass_2", "luminosity_distance")
_LNPDRAW = "lnpdraw_mass1_source_mass2_source_redshift_spin1x_spin1y_spin1z_spin2x_spin2y_spin2z"
_YEAR_S = 3.15576e7
STAGES = ("prepare", "sample", "combine", "report")


def default_pe_cache() -> Path:
    return Path(os.environ.get("GWTC_PE_CACHE", Path.home() / ".cache_gwtc_analysis" / "pe_catalog"))


def _log(msg: str) -> None:
    print(f"[hubble_constant] {msg}", flush=True)


# ---------------------------------------------------------------------------
# events
# ---------------------------------------------------------------------------

def _full_name(gps: float) -> str:
    from astropy.time import Time

    return Time(float(gps), format="gps", scale="utc").datetime.strftime("GW%y%m%d_%H%M%S")


def select_h0_events(release: str, far_threshold: float, min_mass: float, exclude: Iterable[str]) -> pd.DataFrame:
    """BBH events of the release's catalogs: lowest FAR below threshold, inside its runs, both masses >= min_mass."""
    cfg = H0_SENSITIVITY_RELEASES[release]
    exclude = set(exclude)
    best: dict[str, dict] = {}
    for cat in cfg["catalogs"]:
        try:
            raw = gw.fetch_gwtc_events(cat)
        except Exception as e:
            _log(f"WARN: could not fetch {cat}: {type(e).__name__}: {e}")
            continue
        for v in (raw.get("events") or {}).values():
            g, far = v.get("GPS"), v.get("far")
            m1, m2 = v.get("mass_1_source"), v.get("mass_2_source")
            if g is None or far is None or m1 is None or m2 is None or float(far) >= far_threshold:
                continue
            run = next((r for r in cfg["runs"] if OBSERVING_RUNS[r][0] <= float(g) <= OBSERVING_RUNS[r][1]), None)
            if run is None or min(float(m1), float(m2)) < min_mass:
                continue
            name, common = _full_name(g), v.get("commonName") or ""
            if name in exclude or common in exclude or any(e.split("_")[0] == common for e in exclude):
                continue
            if name not in best or best[name]["far_per_yr"] > float(far):
                best[name] = dict(event=name, common_name=common, catalog=cat, run=run, gps=float(g),
                                  far_per_yr=float(far), mass_1_source=float(m1), mass_2_source=float(m2))
    cols = ["event", "common_name", "catalog", "run", "gps", "far_per_yr", "mass_1_source", "mass_2_source"]
    return pd.DataFrame(list(best.values()), columns=cols).sort_values("gps").reset_index(drop=True)


# ---------------------------------------------------------------------------
# PE samples
# ---------------------------------------------------------------------------

def _pick_pe_file(cands: list[dict]) -> dict:
    """Prefer the files with the original PE distance prior (nocosmo); GWTC-4/5 have one combined file."""
    def score(c):
        f = c["filename"].lower()
        return 0 if "nocosmo" in f else 1 if "combined" in f else 2 if "mixed_cosmo" in f else 3
    return sorted(cands, key=score)[0]


def _extract_pe(path: Path, out: Path) -> None:
    """Keep (m1, m2, D_L) of every analysis in a PESummary file, with its recorded priors."""
    import h5py

    part = out.with_suffix(".part")
    with h5py.File(path, "r") as f, h5py.File(part, "w") as o:
        o.attrs["source_file"] = path.name
        for lab in f:
            g = f[lab]
            if not isinstance(g, h5py.Group) or "posterior_samples" not in g:
                continue
            ps = g["posterior_samples"]
            if not all(c in (ps.dtype.names or ()) for c in _PE_COLUMNS):
                continue
            og = o.create_group(lab)
            for c in _PE_COLUMNS:
                og.create_dataset(c, data=np.asarray(ps[c], dtype=float), compression="gzip")
            pr = g.get("priors")
            if pr is not None and "analytic" in pr:
                for k in pr["analytic"]:
                    try:
                        v = pr["analytic"][k][()]
                        v = v[0] if hasattr(v, "__len__") and not isinstance(v, (bytes, str)) and len(v) else v
                        og.attrs[f"prior:{k}"] = v.decode() if isinstance(v, bytes) else str(v)
                    except Exception:
                        pass
    part.replace(out)


def _pe_extract_path(cache: Path, name: str) -> Optional[Path]:
    """Extract of an event, allowing for a +-2 s GPS rounding in the name."""
    from astropy.time import Time

    p = cache / "samples" / f"{name}.h5"
    if p.exists():
        return p
    t0 = Time.strptime(name, "GW%y%m%d_%H%M%S").gps
    for d in (-2, -1, 1, 2):
        q = cache / "samples" / f"{_full_name(t0 + d)}.h5"
        if q.exists():
            return q
    return None


def fetch_pe_samples(events: pd.DataFrame, cache: Path, keep_files: bool = False, workers: int = 3) -> dict[str, Path]:
    """Download (once) the PE file of each event from Zenodo and extract its samples; restartable."""
    from . import parameters_estimation as pe

    (cache / "samples").mkdir(parents=True, exist_ok=True)
    (cache / "files").mkdir(parents=True, exist_ok=True)
    out = {n: p for n in events["event"] if (p := _pe_extract_path(cache, n)) is not None}
    todo = [r for _, r in events.iterrows() if r["event"] not in out]
    if not todo:
        return out
    index = pe.build_zenodo_pe_index(cache_dir=str(cache / "index"), force_refresh=False)
    jobs, missing = [], []
    for r in todo:
        cands = index.get(r["event"]) or index.get(r["common_name"])
        for d in (-2, -1, 1, 2):
            if cands:
                break
            cands = index.get(_full_name(r["gps"] + d))
        if cands:
            jobs.append((r["event"], _pick_pe_file(cands)))
        else:
            missing.append(r["event"])
    if missing:
        raise ValueError(f"No Zenodo PE release found for: {', '.join(missing)}")
    _log(f"downloading the PE files of {len(jobs)} events to {cache / 'files'}")

    def one(name, entry):
        dest = cache / "files" / entry["filename"]
        tmp = dest.with_suffix(dest.suffix + ".part")
        if not (dest.exists() and dest.stat().st_size > 0):
            pe._download_http_with_progress(entry["url"], tmp, desc=name)
            tmp.replace(dest)
        target = cache / "samples" / f"{name}.h5"
        _extract_pe(dest, target)
        if not keep_files:
            dest.unlink(missing_ok=True)
        return name, target

    errors = []
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(one, n, e): n for n, e in jobs}
        for i, fu in enumerate(as_completed(futs), 1):
            try:
                n, p = fu.result()
                out[n] = p
                _log(f"[{i}/{len(jobs)}] {n}: extracted")
            except Exception as e:
                errors.append(f"{futs[fu]}: {type(e).__name__}: {e}")
    if errors:
        raise ValueError("PE download failed (run again to resume): " + "; ".join(errors))
    return out


def _choose_label(labels: Iterable[str]) -> Optional[str]:
    labels = list(labels)
    for want in PE_LABELS:
        if want in labels:
            return want
    alt = [l for l in labels if "IMRPhenomXPHM" in l and "Mixed" not in l]
    return sorted(alt)[0] if alt else None


def _pe_cosmology(name: str):
    from astropy.cosmology import FlatLambdaCDM, Planck15

    return FlatLambdaCDM(H0=67.90, Om0=0.3065) if name.lower() == "planck15_lal" else Planck15


def pe_distance_prior(dl: np.ndarray, desc: str, run: str) -> tuple[np.ndarray, str]:
    """Density (up to a constant) of the PE luminosity-distance prior, from its recorded bilby description.

    The mass priors of the LVK analyses are uniform in the detector-frame component masses, so only the
    distance prior matters. Without a recorded prior, the catalog default is used: D_L^2 up to O3, uniform in
    source-frame comoving volume (Planck15_LAL) from O4.
    """
    desc = desc or ""
    kind = None
    if re.match(r"(bilby\.core\.prior\.)?PowerLaw\(", desc) and re.search(r"alpha=2(\.0)?\b", desc):
        kind = "D_L^2"
    elif "UniformSourceFrame" in desc or "UniformComovingVolume" in desc:
        kind = "UniformSourceFrame" if "UniformSourceFrame" in desc else "UniformComovingVolume"
    elif not desc:
        kind = "UniformSourceFrame" if run.startswith("O4") else "D_L^2"
        desc = f"(not recorded: {kind} assumed for {run})"
    else:
        raise ValueError(f"Unsupported PE distance prior: {desc[:120]}")
    if kind == "D_L^2":
        return dl ** 2, desc
    m = re.search(r"cosmology='?([A-Za-z0-9_]+)", desc)
    cosmo = _pe_cosmology(m.group(1) if m else "Planck15_LAL")
    z = np.linspace(0, 5.0, 20001)
    d = cosmo.luminosity_distance(z).value
    dens = cosmo.differential_comoving_volume(z).value / np.gradient(d, z)
    if kind == "UniformSourceFrame":
        dens = dens / (1 + z)
    return np.interp(dl, d, dens), desc


# ---------------------------------------------------------------------------
# injections
# ---------------------------------------------------------------------------

def h0_sensitivity_path(sensitivity_file: str | Path | None, release: str) -> Path:
    from .catalogs import _zenodo_sensitivity_file

    if sensitivity_file:
        p = Path(sensitivity_file).expanduser()
        if not p.exists():
            raise ValueError(f"Sensitivity file not found: {p}")
        return p
    if release not in H0_SENSITIVITY_RELEASES:
        raise ValueError(f"Unknown sensitivity release {release!r}; choose from {', '.join(H0_SENSITIVITY_RELEASES)}")
    cfg = H0_SENSITIVITY_RELEASES[release]
    return _zenodo_sensitivity_file(cfg["record"], cfg["label"], cfg["file_re"], tag="hubble_constant")


def detector_frame_injections(path: Path, far_threshold: float, snr_threshold: float) -> dict:
    """Found injections in (m1_det, m2_det, D_L), with their draw density in these variables.

    Found: semi-analytic network SNR above `snr_threshold` (O1+O2) or lowest search FAR below `far_threshold`.
    The spin part of the draw is divided out (the population spins are then the injected, isotropic ones), the
    density is carried to the detector frame by 1/[(1+z)^2 dD_L/dz], and the mixture weights enter as p/w.
    """
    import h5py

    with h5py.File(path, "r") as fi:
        e = fi["events"]
        searches = [s.decode() if isinstance(s, bytes) else str(s) for s in fi.attrs["searches"]]
        far = np.min([e[f"{s}_far"][:] for s in searches], axis=0)
        sel = far < far_threshold
        if "semianalytic_observed_phase_maximized_snr_net" in e.dtype.names:
            sel |= e["semianalytic_observed_phase_maximized_snr_net"][:] > snr_threshold
        get = lambda k: e[k][:][sel]
        z = get("redshift")
        a1 = np.sqrt(get("spin1x") ** 2 + get("spin1y") ** 2 + get("spin1z") ** 2)
        a2 = np.sqrt(get("spin2x") ** 2 + get("spin2y") ** 2 + get("spin2z") ** 2)
        ln_spin = -np.log(4 * np.pi * np.maximum(a1, 1e-12) ** 2) - np.log(4 * np.pi * np.maximum(a2, 1e-12) ** 2)
        ln_p = get(_LNPDRAW) - ln_spin
        p_det = np.exp(ln_p) / ((1 + z) ** 2 * get("dluminosity_distance_dredshift")) / get("weights")
        return dict(mass_1=get("mass1_source") * (1 + z), mass_2=get("mass2_source") * (1 + z),
                    luminosity_distance=get("luminosity_distance"), prior=p_det,
                    ntotal=float(fi.attrs["total_generated"]), Tobs=float(fi.attrs["total_analysis_time"]) / _YEAR_S,
                    n_recorded=int(len(sel)))


# ---------------------------------------------------------------------------
# stages
# ---------------------------------------------------------------------------

def prepare_inputs(workdir: Path, release: str, sensitivity_file, far_threshold: float, snr_threshold: float,
                   min_mass: float, exclude: Iterable[str], pe_cache: Path, keep_pe_files: bool,
                   max_pe_samples: int = 5000) -> pd.DataFrame:
    """Write <workdir>/inputs.h5 (events + injections) and <workdir>/events.tsv."""
    import h5py

    workdir.mkdir(parents=True, exist_ok=True)
    inj_path = h0_sensitivity_path(sensitivity_file, release)
    events = select_h0_events(release, far_threshold, min_mass, exclude)
    if events.empty:
        raise ValueError("No event passes the selection")
    _log(f"{len(events)} BBH events: " + ", ".join(f"{r} {n}" for r, n in events["run"].value_counts().sort_index().items()))
    extracts = fetch_pe_samples(events, pe_cache, keep_files=keep_pe_files)
    rng = np.random.default_rng(12345)
    labels, priors, nsamp = [], [], []
    with h5py.File(workdir / "inputs.h5.part", "w") as h:
        for _, r in events.iterrows():
            with h5py.File(extracts[r["event"]], "r") as f:
                lab = _choose_label(f)
                if lab is None:
                    raise ValueError(f"{r['event']}: no IMRPhenomXPHM analysis in {extracts[r['event']].name} "
                                     f"(labels {list(f)})")
                g = f[lab]
                m1, m2, dl = g["mass_1"][:], g["mass_2"][:], g["luminosity_distance"][:]
                desc = str(g.attrs.get("prior:luminosity_distance", ""))
            idx = rng.permutation(len(dl))[:max_pe_samples]   # shuffled: the sampler takes the first N
            prior, desc = pe_distance_prior(dl[idx], desc, r["run"])
            gr = h.create_group(r["event"])
            for k, v in (("m1", m1[idx]), ("m2", m2[idx]), ("dl", dl[idx]), ("prior", prior)):
                gr.create_dataset(k, data=v)
            gr.attrs.update(run=r["run"], label=lab, prior_desc=desc[:200])
            labels.append(lab); priors.append(re.sub(r"\(.*", "", desc).replace("bilby.gw.prior.", "")); nsamp.append(len(idx))
        inj = detector_frame_injections(inj_path, far_threshold, snr_threshold)
        gi = h.create_group("_injections")
        for k in ("mass_1", "mass_2", "luminosity_distance", "prior"):
            gi.create_dataset(k, data=inj[k])
        gi.attrs.update(ntotal=inj["ntotal"], Tobs=inj["Tobs"], n_found=len(inj["prior"]), source=inj_path.name)
    (workdir / "inputs.h5.part").replace(workdir / "inputs.h5")
    events["pe_label"], events["pe_distance_prior"], events["pe_samples"] = labels, priors, nsamp
    events.to_csv(workdir / "events.tsv", sep="\t", index=False)
    _log(f"{len(inj['prior'])} found injections of {inj['n_recorded']} recorded ({inj_path.name}); "
         f"inputs written to {workdir / 'inputs.h5'}")
    return events


def _runner_command(python: str) -> tuple[list[str], dict]:
    env = dict(os.environ)
    lib = Path(python).expanduser().parent.parent / "lib"
    if lib.is_dir():   # icarogw's compiled dependencies need the environment's libstdc++
        env["LD_LIBRARY_PATH"] = f"{lib}:{env['LD_LIBRARY_PATH']}" if env.get("LD_LIBRARY_PATH") else str(lib)
    return [str(Path(python).expanduser()), str(Path(__file__).with_name("h0_icarogw.py"))], env


def _check_runner(python: str) -> None:
    cmd, env = _runner_command(python)
    r = subprocess.run([cmd[0], "-c", "import icarogw, bilby"], env=env, capture_output=True, text=True,
                       cwd=str(Path.home()))
    if r.returncode != 0:
        raise ValueError(f"icarogw/bilby cannot be imported by {cmd[0]} (install icarogw there, or pass "
                         f"--icarogw-python): {r.stderr.strip().splitlines()[-1] if r.stderr.strip() else ''}")


def _run_runner(python: str, args: list[str]) -> None:
    cmd, env = _runner_command(python)
    _log("running " + " ".join(cmd[1:2] + args))
    r = subprocess.run(cmd + args, env=env)
    if r.returncode != 0:
        raise ValueError(f"icarogw runner failed ({' '.join(args[:1])}, exit code {r.returncode})")


def _last_line(log: Path, prefix: str = "[h0]") -> str:
    try:
        lines = [l for l in log.read_text(errors="replace").splitlines() if l.startswith(prefix)]
    except OSError:
        return ""
    return lines[-1] if lines else ""


def run_seeds_parallel(python: str, workdir: Path, seeds: list[int], parallel: int, run_args: list[str]) -> None:
    """Up to `parallel` sampler runs at a time on this machine, each logging to <workdir>/logs/run_seed<N>.log.

    Interrupting (Ctrl-C) stops the runs, which write their checkpoint and resume on the next launch.
    """
    import signal
    import time

    def _stop(signum, frame):
        raise KeyboardInterrupt

    previous = signal.signal(signal.SIGTERM, _stop)   # a killed launcher also stops (and checkpoints) its runs
    cmd, env = _runner_command(python)
    logs = workdir / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    todo, running, failed = list(seeds), {}, []
    _log(f"{len(seeds)} run(s), {parallel} at a time; logs in {logs}")
    try:
        while todo or running:
            while todo and len(running) < parallel:
                seed = todo.pop(0)
                log = logs / f"run_seed{seed}.log"
                fh = open(log, "a")
                running[seed] = (subprocess.Popen(cmd + ["run", "--workdir", str(workdir), "--seed", str(seed)] + run_args,
                                                  env=env, stdout=fh, stderr=subprocess.STDOUT,
                                                  start_new_session=True), fh, log)   # signals reach the launcher only
                _log(f"seed {seed}: started (pid {running[seed][0].pid})")
            time.sleep(5)
            for seed, (proc, fh, log) in list(running.items()):
                if proc.poll() is None:
                    continue
                fh.close()
                del running[seed]
                if proc.returncode == 0:
                    _log(f"seed {seed}: finished. {_last_line(log)}")
                else:
                    failed.append(seed)
                    _log(f"seed {seed}: FAILED (exit code {proc.returncode}): {_last_line(log) or 'see ' + str(log)}")
    except KeyboardInterrupt:
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        _log(f"stopping seed(s) {sorted(running)}: they write their checkpoint first")
        for proc, fh, _ in running.values():
            proc.terminate()
        for proc, fh, _ in running.values():
            proc.wait()
            fh.close()
        raise ValueError("interrupted: the running seeds wrote their checkpoint and resume on the next launch")
    finally:
        signal.signal(signal.SIGTERM, previous)
        signal.signal(signal.SIGINT, signal.default_int_handler)
    if failed:
        raise ValueError(f"sampler run(s) failed for seed(s) {failed}; see {logs}")


def _plot_h0(post: pd.DataFrame, published: Optional[dict], out_png: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid, blue, orange = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6", "#eb6834"
    h0 = post["H0"].to_numpy()
    q = np.quantile(h0, [0.05, 0.5, 0.95])
    fig, ax = plt.subplots(figsize=(8.2, 4.2), dpi=150)
    fig.patch.set_facecolor(surface); ax.set_facecolor(surface)
    ax.hist(h0, bins=np.linspace(10, 200, 39), density=True, color=blue, alpha=0.85, linewidth=0, zorder=2,
            label=f"This run: {q[1]:.0f}, 90% {q[0]:.0f}–{q[2]:.0f}")
    if published:
        ax.axvspan(published["lo90"], published["hi90"], color=orange, alpha=0.12, zorder=1, linewidth=0)
        ax.axvline(published["median"], color=orange, lw=2, zorder=3,
                   label=f"{published['ref']}: {published['median']:.0f}, 90% {published['lo90']:.0f}–{published['hi90']:.0f}")
    for x, t, ha, dx in ((67.4, "Planck", "right", -3), (73.0, "SH0ES", "left", 3)):
        ax.axvline(x, color=ink2, lw=0.8, ls=(0, (3, 3)), zorder=1)
        ax.annotate(t, xy=(x, 1), xycoords=("data", "axes fraction"), xytext=(dx, -12), textcoords="offset points",
                    fontsize=8, color=ink2, ha=ha)
    ax.set_xlim(10, 200)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.08)
    ax.set_xlabel("H$_0$ (km s$^{-1}$ Mpc$^{-1}$)", color=ink2, fontsize=10)
    ax.set_ylabel("Posterior density", color=ink2, fontsize=10)
    ax.grid(False)
    ax.grid(axis="y", color=grid, lw=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(grid)
    ax.tick_params(colors=ink2, labelsize=9)
    ax.set_title("Hubble constant from the BBH mass spectrum (Power Law + Peak)", color=ink, fontsize=10.5, loc="left",
                 pad=34)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=ink2, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2,
              borderaxespad=0.3, handlelength=1.2)
    fig.tight_layout()
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def write_h0_report(workdir: Path, out_report_html: Path, out_summary_tsv: Optional[Path], release: str) -> pd.DataFrame:
    summary = json.loads((workdir / "summary.json").read_text())
    post = pd.read_csv(workdir / "posterior.tsv", sep="\t")
    events = pd.read_csv(workdir / "events.tsv", sep="\t") if (workdir / "events.tsv").exists() else pd.DataFrame()
    published = H0_SENSITIVITY_RELEASES.get(release, {}).get("published")
    rows = [dict(parameter=k, median=q[2], minus_68=q[2] - q[1], plus_68=q[3] - q[2], low_90=q[0], high_90=q[4])
            for k, q in summary["quantiles"].items()]
    table = pd.DataFrame(rows)
    if out_summary_tsv:
        Path(out_summary_tsv).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_summary_tsv, sep="\t", index=False, float_format="%.6g")
    plots = workdir / "plots"
    plots.mkdir(exist_ok=True)
    images = [_plot_h0(post, published, plots / "h0_posterior.png")]
    if (workdir / "corner.png").exists():
        images.append(workdir / "corner.png")
    q = summary["quantiles"]["H0"]
    s = summary.get("settings", {})
    paras = [
        f"H<sub>0</sub> = <b>{q[2]:.1f} (+{q[3] - q[2]:.1f} / −{q[2] - q[1]:.1f}) km/s/Mpc</b> (median, 68%); "
        f"90%: {q[0]:.1f}–{q[4]:.1f}. Spectral siren with {len(events) or '?'} BBH events: the redshift comes from the "
        "source-frame mass distribution (Power Law + Peak, fitted together with H<sub>0</sub>) and the "
        "Madau–Dickinson rate evolution, with flat ΛCDM (Ω<sub>m</sub> = 0.3065).",
        f"{summary['n_runs']} dynesty run(s) of {s.get('nlive', '?')} live points, {summary['n_samples']} posterior samples, "
        f"ln Z = {summary['log_evidence']:.2f} ± {summary['log_evidence_err']:.2f} (runs: "
        + ", ".join(f"{x:.1f}" for x in summary["run_log_evidences"]) + f"). PE samples per event: {s.get('pe_samples', '?')}; "
        f"fraction of the found injections used: {s.get('inj_fraction', '?')}.",
    ]
    if published:
        paras.append(f"Published ({published['ref']}): {published['median']} (+{published['plus']} / −{published['minus']}), "
                     f"90%: {published['lo90']}–{published['hi90']} km/s/Mpc.")
    d = summary.get("diagnostics")
    if d:
        ok = d["neff_inj_min"] >= d["neff_inj_threshold"] and d["neff_pe_min"] >= d["neff_pe_threshold"]
        paras.append(
            f"Numerical stability over {d['n_points']} posterior draws: effective injections min {d['neff_inj_min']:.0f}, "
            f"median {d['neff_inj_median']:.0f} (threshold {d['neff_inj_threshold']}); smallest per-event effective PE "
            f"samples min {d['neff_pe_min']:.1f}, median {d['neff_pe_median_of_min']:.1f} (threshold {d['neff_pe_threshold']}); "
            f"lowest events: {', '.join(d['lowest_neff_pe_events'])}. " + ("OK." if ok else "<b>Below threshold: increase "
            "--inj-fraction or --pe-samples.</b>"))
    paras.append("The upper part of the H<sub>0</sub> interval depends on the prior bound (uniform 10–200 km/s/Mpc).")
    tables = [("Posterior (median, 68% and 90% intervals)",
               table.to_html(index=False, float_format=lambda x: f"{x:.3g}"))]
    if not events.empty:
        tables.append(("Events", events.to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}")))
    Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
    write_simple_html_report(out_report_html, title="Hubble constant (spectral siren)", paragraphs=paras,
                             images=images, tables=tables)
    _log(f"report written to {out_report_html}")
    return table


def run_hubble_constant(
    stages: Iterable[str] = STAGES,
    workdir: str | Path = "hubble_constant_run",
    out_report_html: Optional[str | Path] = "hubble_constant.html",
    out_summary_tsv: Optional[str | Path] = "hubble_constant.tsv",
    sensitivity_release: str = H0_DEFAULT_RELEASE,
    sensitivity_file: Optional[str | Path] = None,
    far_threshold: float = 0.25,
    snr_threshold: float = 10.0,
    min_mass: float = 3.0,
    exclude: Iterable[str] = H0_DEFAULT_EXCLUDE,
    pe_cache: Optional[str | Path] = None,
    keep_pe_files: bool = False,
    seeds: Iterable[int] = (1,),
    parallel: int = 1,
    nlive: int = 100,
    npool: int = 4,
    naccept: int = 60,
    pe_samples: int = 1500,
    inj_fraction: float = 0.1,
    icarogw_python: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """Spectral-siren H0 with the Power Law + Peak model; see the module docstring for the stages."""
    stages = [s for s in STAGES if s in set(stages)]
    workdir = Path(workdir).expanduser().resolve()
    python = icarogw_python or sys.executable
    if "prepare" in stages:
        prepare_inputs(workdir, sensitivity_release, sensitivity_file, far_threshold, snr_threshold, min_mass,
                       exclude, Path(pe_cache).expanduser() if pe_cache else default_pe_cache(), keep_pe_files)
    if ("sample" in stages or "combine" in stages):
        _check_runner(python)
    if "sample" in stages:
        if not (workdir / "inputs.h5").exists():
            raise ValueError(f"{workdir / 'inputs.h5'} not found: run the prepare stage first")
        seeds = list(dict.fromkeys(int(x) for x in seeds))
        run_args = ["--nlive", str(nlive), "--npool", str(npool), "--naccept", str(naccept),
                    "--pe-samples", str(pe_samples), "--inj-fraction", str(inj_fraction)]
        parallel = max(1, min(int(parallel), len(seeds)))
        ncpu = os.cpu_count() or 1
        if parallel * npool > ncpu:
            _log(f"WARN: {parallel} parallel runs x {npool} processes = {parallel * npool} > {ncpu} CPUs; "
                 "reduce --parallel or --npool")
        if parallel == 1:
            for seed in seeds:
                _run_runner(python, ["run", "--workdir", str(workdir), "--seed", str(seed)] + run_args)
        else:
            run_seeds_parallel(python, workdir, seeds, parallel, run_args)
    if "combine" in stages:
        _run_runner(python, ["combine", "--workdir", str(workdir)])
    if "report" in stages and out_report_html:
        return write_h0_report(workdir, Path(out_report_html), Path(out_summary_tsv) if out_summary_tsv else None,
                               sensitivity_release)
    return None

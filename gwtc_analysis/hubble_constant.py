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
- ``reweight``: when the runs used a subset of the injections, their posterior reweighted to all of
  them (importance weights exp(ln L_all - ln L_runs), checked by their effective sample size).
- ``report``: HTML report.

With ``galaxy_catalog`` (an icarogw line-of-sight catalog built by the ``galaxy_catalog`` mode), the analysis is
a dark siren with a galaxy catalog: ``prepare`` also keeps the sky position of every PE sample and records the
catalog in ``inputs.h5``, and the likelihood (``h0_icarogw.build_likelihood``) uses icarogw's
``CBC_catalog_vanilla_rate``: the redshift prior of each event along its line of sight is the in-catalog galaxy
density plus the out-of-catalog completeness term. The LVK injection files do not record sky positions; the
injections enter through the sky-averaged galaxy density, for which their position does not matter, and are given
isotropic positions.

With ``inj_fraction="auto"`` (the default), a probe measures before sampling the speed and the
accuracy of likelihoods built on subsets of the injections, and the fastest reliable strategy is
chosen (`choose_injection_fraction`): sample with a subset and reweight, or sample with all.

The default setup reproduces the spectral-siren measurements of the GWTC-4.0 cosmology paper
(arXiv:2509.04348, published version v3): H0 = 105.5 (+46.4 / -35.8) km/s/Mpc with the Power Law + Peak
mass model (`plp`), 72.3 (+42.5 / -25.6) with the Multi Peak model (`mltp`).
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
from .h0_icarogw import choose_injection_fraction  # noqa: E402,F401  (one implementation, in the runner)

from .catalogs import CATALOG_RUNS, OBSERVING_RUNS, SEMI_ANALYTIC_RUNS, _in_runs, catalog_runs  # noqa: E402

# Cumulative search-sensitivity releases with semi-analytic O1+O2 injections, so that the whole
# catalog since O1 can be used (catalog_registry)
from . import catalog_registry as _reg  # noqa: E402

H0_SENSITIVITY_RELEASES = {
    k: dict(record=r.record, label=r.semi_label, file_re=r.semi_file_re, runs=r.runs,
            catalogs=_reg.gwosc_lists(_reg.release_catalogs(k)), published=r.published or None)
    for k, r in sorted(_reg.SENSITIVITY_RELEASES.items())
}
H0_DEFAULT_RELEASE = _reg.DEFAULT_H0_RELEASE
# GWTC-4.0 cosmology: GW231123 (PE systematics) and GW200105 (NSBH) are not used
H0_DEFAULT_EXCLUDE = ("GW231123_135430", "GW200105_162426")
# PE labels, in order of preference (the GWTC-4.0 cosmology choice first)
PE_LABELS = ("C00:IMRPhenomXPHM-SpinTaylor", "C01:IMRPhenomXPHM", "C00:IMRPhenomXPHM")
_PE_COLUMNS = ("mass_1", "mass_2", "luminosity_distance")
_PE_SKY = ("ra", "dec")                    # kept when present, for the dark siren with a galaxy catalog
_LNPDRAW = "lnpdraw_mass1_source_mass2_source_redshift_spin1x_spin1y_spin1z_spin2x_spin2y_spin2z"
_YEAR_S = 3.15576e7
STAGES = ("prepare", "sample", "combine", "reweight", "report")
PROBE_FRACTIONS = (0.1, 0.2, 0.5)
MASS_MODELS = {"plp": "Power Law + Peak", "mltp": "Multi Peak"}
H0_DEFAULT_MASS_MODEL = "plp"


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


def select_h0_events(release: str, far_threshold: float, min_mass: float, exclude: Iterable[str],
                     runs: Optional[Iterable[str]] = None, extra_lists: Iterable[str] = ()) -> pd.DataFrame:
    """BBH events of the release's catalogs: lowest FAR at most the threshold, inside its runs (or `runs`), both
    masses >= min_mass."""
    cfg = dict(H0_SENSITIVITY_RELEASES[release])
    if runs:
        cfg["runs"] = tuple(r for r in cfg["runs"] if r in set(runs))
    exclude = set(exclude)
    best: dict[str, dict] = {}
    for cat in tuple(cfg["catalogs"]) + tuple(extra_lists):     # extra_lists: of update catalogs (GWTC-4.1)
        try:
            raw = gw.fetch_gwtc_events(cat)
        except Exception as e:
            _log(f"WARN: could not fetch {cat}: {type(e).__name__}: {e}")
            continue
        for v in (raw.get("events") or {}).values():
            g, far = v.get("GPS"), v.get("far")
            m1, m2 = v.get("mass_1_source"), v.get("mass_2_source")
            # the published FARs are rounded to two decimals (0.25 stands for 0.245-0.255), while the LVK
            # analyses cut at full precision: the rounded value is compared inclusively
            if g is None or far is None or m1 is None or m2 is None or float(far) > far_threshold:
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
    """Keep (m1, m2, D_L) and the sky position (ra, dec) of every analysis in a PESummary file, with its recorded
    priors."""
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
            for c in _PE_COLUMNS + tuple(c for c in _PE_SKY if c in ps.dtype.names):
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


def _extract_has_sky(path: Path) -> bool:
    """Whether an extract holds the sky position (extracts made before it was kept do not)."""
    import h5py

    with h5py.File(path, "r") as f:
        lab = _choose_label(f)
        return lab is not None and all(c in f[lab] for c in _PE_SKY)


def fetch_pe_samples(events: pd.DataFrame, cache: Path, keep_files: bool = False, workers: int = 3,
                     prefer_catalogs: Iterable[str] = (), need_sky: bool = False) -> dict[str, Path]:
    """Download (once) the PE file of each event from Zenodo and extract its samples; restartable.

    With `need_sky`, an extract without the sky position is made again (from the kept PE file, else downloaded)."""
    from . import parameters_estimation as pe

    (cache / "samples").mkdir(parents=True, exist_ok=True)
    (cache / "files").mkdir(parents=True, exist_ok=True)
    out = {n: p for n in events["event"] if (p := _pe_extract_path(cache, n)) is not None
           and (not need_sky or _extract_has_sky(p))}
    todo = [r for _, r in events.iterrows() if r["event"] not in out]
    if not todo:
        return out
    record_ids = None
    if prefer_catalogs:       # update catalogs first (their PE files are then chosen), then the default ones
        from .data_repo import resolve_zenodo_records

        record_ids = [r.record_id for c in prefer_catalogs for r in resolve_zenodo_records(c, None)]
        record_ids += pe.zenodo_pe_record_ids()
    index = pe.build_zenodo_pe_index(cache_dir=str(cache / "index"), force_refresh=False, record_ids=record_ids)
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


def detector_frame_injections(path: Path, far_threshold: float, snr_threshold: float,
                              runs: Optional[Iterable[str]] = None) -> dict:
    """Found injections in (m1_det, m2_det, D_L), with their draw density in these variables.

    Found: semi-analytic network SNR above `snr_threshold` for the O1+O2 injections, lowest search FAR below
    `far_threshold` for the others; restricted to the observing `runs` if given.
    The spin part of the draw is divided out (the population spins are then the injected, isotropic ones), the
    density is carried to the detector frame by 1/[(1+z)^2 dD_L/dz], and the mixture weights enter as p/w.
    """
    import h5py

    with h5py.File(path, "r") as fi:
        e = fi["events"]
        searches = [s.decode() if isinstance(s, bytes) else str(s) for s in fi.attrs["searches"]]
        from .catalogs import _found_mask

        sel = _found_mask(e, searches, far_threshold, snr_threshold)
        if runs:
            sel &= _in_runs(e["time_geocenter"][:], runs)
        get = lambda k: e[k][:][sel]
        z = get("redshift")
        a1 = np.sqrt(get("spin1x") ** 2 + get("spin1y") ** 2 + get("spin1z") ** 2)
        a2 = np.sqrt(get("spin2x") ** 2 + get("spin2y") ** 2 + get("spin2z") ** 2)
        ln_spin = -np.log(4 * np.pi * np.maximum(a1, 1e-12) ** 2) - np.log(4 * np.pi * np.maximum(a2, 1e-12) ** 2)
        ln_p = get(_LNPDRAW) - ln_spin
        p_det = np.exp(ln_p) / ((1 + z) ** 2 * get("dluminosity_distance_dredshift")) / get("weights")
        m1s, m2s = get("mass1_source"), get("mass2_source")
        return dict(mass_1=m1s * (1 + z), mass_2=m2s * (1 + z),
                    luminosity_distance=get("luminosity_distance"), prior=p_det,
                    chi_eff=(m1s * get("spin1z") + m2s * get("spin2z")) / (m1s + m2s),
                    mass_ratio=np.minimum(m1s, m2s) / np.maximum(m1s, m2s),
                    ntotal=float(fi.attrs["total_generated"]), Tobs=float(fi.attrs["total_analysis_time"]) / _YEAR_S,
                    n_recorded=int(len(sel)))


def galaxy_catalog_info(path: str | Path) -> dict:
    """Path, group names and construction settings of an icarogw catalog file made by the galaxy_catalog mode."""
    import h5py

    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise ValueError(f"Galaxy catalog not found: {p}")
    with h5py.File(p, "r") as h:
        raw = h.attrs.get("gwtc_analysis_settings")
        if raw is None:
            raise ValueError(f"{p} is not a finished galaxy_catalog file (no gwtc_analysis_settings attribute)")
        st = json.loads(raw)
        if st["grouping"] not in h or st["subgrouping"] not in h[st["grouping"]]:
            raise ValueError(f"{p}: group {st['grouping']}/{st['subgrouping']} missing; run the finish stage")
    return dict(path=p, grouping=st["grouping"], subgrouping=st["subgrouping"], settings=st)


# ---------------------------------------------------------------------------
# stages
# ---------------------------------------------------------------------------

def prepare_inputs(workdir: Path, release: str, sensitivity_file, far_threshold: float, snr_threshold: float,
                   min_mass: float, exclude: Iterable[str], pe_cache: Path, keep_pe_files: bool,
                   max_pe_samples: int = 5000, runs: Optional[Iterable[str]] = None,
                   updates: tuple[str, ...] = (), galaxy_catalog: Optional[Path] = None) -> pd.DataFrame:
    """Write <workdir>/inputs.h5 (events + injections), <workdir>/events.tsv and <workdir>/selection.json.

    With `galaxy_catalog`, the sky positions of the PE samples are kept and the catalog is recorded in inputs.h5."""
    import h5py

    workdir.mkdir(parents=True, exist_ok=True)
    release_runs = H0_SENSITIVITY_RELEASES[release]["runs"]
    runs = tuple(release_runs if not runs else [r for r in release_runs if r in set(runs)])
    missing = [r for r in (runs or ()) if r not in release_runs]
    if missing:
        raise ValueError(f"The {release} sensitivity release does not cover {', '.join(missing)}")
    inj_path = h0_sensitivity_path(sensitivity_file, release)
    events = select_h0_events(release, far_threshold, min_mass, exclude, runs,
                              extra_lists=_reg.gwosc_lists(updates) if updates else ())
    if events.empty:
        raise ValueError("No event passes the selection")
    _log(f"{len(events)} BBH events: " + ", ".join(f"{r} {n}" for r, n in events["run"].value_counts().sort_index().items()))
    catalog = galaxy_catalog_info(galaxy_catalog) if galaxy_catalog else None
    extracts = fetch_pe_samples(events, pe_cache if not updates else pe_cache / "_".join(updates),
                                keep_files=keep_pe_files, prefer_catalogs=updates, need_sky=catalog is not None)
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
                sky = {c: g[c][:] for c in _PE_SKY} if catalog else {}
                desc = str(g.attrs.get("prior:luminosity_distance", ""))
            idx = rng.permutation(len(dl))[:max_pe_samples]   # shuffled: the sampler takes the first N
            prior, desc = pe_distance_prior(dl[idx], desc, r["run"])
            gr = h.create_group(r["event"])
            for k, v in (("m1", m1[idx]), ("m2", m2[idx]), ("dl", dl[idx]), ("prior", prior)):
                gr.create_dataset(k, data=v)
            for k, v in sky.items():
                gr.create_dataset(k, data=v[idx])
            gr.attrs.update(run=r["run"], label=lab, prior_desc=desc[:200])
            labels.append(lab); priors.append(re.sub(r"\(.*", "", desc).replace("bilby.gw.prior.", "")); nsamp.append(len(idx))
        inj = detector_frame_injections(inj_path, far_threshold, snr_threshold, runs)
        gi = h.create_group("_injections")
        for k in ("mass_1", "mass_2", "luminosity_distance", "prior"):
            gi.create_dataset(k, data=inj[k])
        gi.attrs.update(ntotal=inj["ntotal"], Tobs=inj["Tobs"], n_found=len(inj["prior"]), source=inj_path.name,
                        runs=",".join(runs))
        if catalog:
            # isotropic positions: the injections enter through the sky-averaged galaxy density (see module doc)
            sky_rng = np.random.default_rng(777)
            n_inj = len(inj["prior"])
            gi.create_dataset("ra", data=sky_rng.uniform(0, 2 * np.pi, n_inj))
            gi.create_dataset("dec", data=np.arcsin(sky_rng.uniform(-1, 1, n_inj)))
            h.attrs.update(galaxy_catalog=str(catalog["path"]), catalog_grouping=catalog["grouping"],
                           catalog_subgrouping=catalog["subgrouping"], catalog_settings=json.dumps(catalog["settings"]))
    (workdir / "inputs.h5.part").replace(workdir / "inputs.h5")
    events["pe_label"], events["pe_distance_prior"], events["pe_samples"] = labels, priors, nsamp
    events.to_csv(workdir / "events.tsv", sep="\t", index=False)
    (workdir / "selection.json").write_text(json.dumps(dict(release=release, runs=list(runs), updates=list(updates),
                                                            default_runs=runs == tuple(release_runs) and not updates,
                                                            galaxy_catalog=catalog["settings"] if catalog else None)))
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
    jobs = [(f"run_seed{seed}", ["run", "--workdir", str(workdir), "--seed", str(seed)] + run_args) for seed in seeds]
    try:
        run_jobs_parallel(python, workdir, jobs, parallel)
    except ValueError as e:
        if str(e).startswith("job(s) failed"):
            bad = [int(n[len("run_seed"):]) for n in jobs_failed(e)]
            raise ValueError(f"sampler run(s) failed for seed(s) {bad}; see {workdir / 'logs'}") from None
        raise


def jobs_failed(err: ValueError) -> list[str]:
    import ast

    txt = str(err)
    return ast.literal_eval(txt[txt.index("["):txt.index("]") + 1])


def run_jobs_parallel(python: str, workdir: Path, jobs: list[tuple[str, list[str]]], parallel: int) -> None:
    """Runner jobs (name, arguments), up to `parallel` at a time, each logging to <workdir>/logs/<name>.log.

    Interrupting (Ctrl-C) stops the jobs; sampler runs write their checkpoint first."""
    import signal
    import time

    def _stop(signum, frame):
        raise KeyboardInterrupt

    previous = signal.signal(signal.SIGTERM, _stop)   # a killed launcher also stops (and checkpoints) its jobs
    cmd, env = _runner_command(python)
    logs = workdir / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    todo, running, failed = list(jobs), {}, []
    _log(f"{len(jobs)} job(s), {parallel} at a time; logs in {logs}")
    try:
        while todo or running:
            while todo and len(running) < parallel:
                name, args = todo.pop(0)
                log = logs / f"{name}.log"
                fh = open(log, "a")
                running[name] = (subprocess.Popen(cmd + args, env=env, stdout=fh, stderr=subprocess.STDOUT,
                                                  start_new_session=True), fh, log)   # signals reach the launcher only
                _log(f"{name}: started (pid {running[name][0].pid})")
            time.sleep(5)
            for name, (proc, fh, log) in list(running.items()):
                if proc.poll() is None:
                    continue
                fh.close()
                del running[name]
                if proc.returncode == 0:
                    _log(f"{name}: finished. {_last_line(log)}")
                else:
                    failed.append(name)
                    _log(f"{name}: FAILED (exit code {proc.returncode}): {_last_line(log) or 'see ' + str(log)}")
    except KeyboardInterrupt:
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        _log(f"stopping {sorted(running)}")
        for proc, fh, _ in running.values():
            proc.terminate()
        for proc, fh, _ in running.values():
            proc.wait()
            fh.close()
        raise ValueError("interrupted: the running jobs were stopped (sampler runs resume from their checkpoint)")
    finally:
        signal.signal(signal.SIGTERM, previous)
        signal.signal(signal.SIGINT, signal.default_int_handler)
    if failed:
        raise ValueError(f"job(s) failed: {failed}; see {logs}")


def _plot_h0(post: pd.DataFrame, published: Optional[dict], out_png: Path, model_name: str = "Power Law + Peak") -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid, blue, orange = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6", "#eb6834"
    h0 = post["H0"].to_numpy()
    q = np.quantile(h0, [0.05, 0.5, 0.95])
    fig, ax = plt.subplots(figsize=(8.2, 4.2), dpi=150)
    fig.patch.set_facecolor(surface); ax.set_facecolor(surface)
    ax.hist(h0, bins=np.linspace(10, 200, 39), density=True, color=blue, alpha=0.85, linewidth=0, zorder=2,
            label=f"gwtc_analysis: {q[1]:.0f}, 90% {q[0]:.0f}–{q[2]:.0f}")
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
    ax.set_title(f"Hubble constant from the BBH mass spectrum ({model_name})", color=ink, fontsize=10.5, loc="left",
                 pad=34)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=ink2, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2,
              borderaxespad=0.3, handlelength=1.2)
    fig.tight_layout()
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def write_h0_report(workdir: Path, out_report_html: Path, out_summary_tsv: Optional[Path], release: str) -> pd.DataFrame:
    summary = json.loads((workdir / "summary.json").read_text())
    rw = summary.get("reweighted")
    use_rw = bool(rw) and (workdir / "posterior_reweighted.tsv").exists()
    post = pd.read_csv(workdir / ("posterior_reweighted.tsv" if use_rw else "posterior.tsv"), sep="\t")
    events = pd.read_csv(workdir / "events.tsv", sep="\t") if (workdir / "events.tsv").exists() else pd.DataFrame()
    model = summary.get("mass_model", "plp")
    model_name = MASS_MODELS.get(model, model)
    sel_file = workdir / "selection.json"
    selection = json.loads(sel_file.read_text()) if sel_file.exists() else {}
    cat_settings = selection.get("galaxy_catalog")
    published = (H0_SENSITIVITY_RELEASES.get(release, {}).get("published") or {}).get(
        f"dark_{model}" if cat_settings else model)
    if cat_settings and published and not _matches_published_catalog(cat_settings, published):
        published = None          # the published dark siren is for its own catalog settings
    if selection and not selection.get("default_runs", True):
        published = None          # the published value is for all the runs of the release
    quantiles = rw["quantiles"] if use_rw else summary["quantiles"]
    rows = [dict(parameter=k, median=q[2], minus_68=q[2] - q[1], plus_68=q[3] - q[2], low_90=q[0], high_90=q[4])
            for k, q in quantiles.items()]
    table = pd.DataFrame(rows)
    if out_summary_tsv:
        Path(out_summary_tsv).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_summary_tsv, sep="\t", index=False, float_format="%.6g")
    plots = workdir / "plots"
    plots.mkdir(exist_ok=True)
    images = [_plot_h0(post, published, plots / "h0_posterior.png", model_name)]
    rate_para = None
    if {"gamma", "kappa", "zp"} <= set(post.columns):
        from . import rate_evolution as rev

        z_ev = None
        if (workdir / "inputs.h5").exists():          # the farthest event: largest median distance of its samples
            import h5py
            from astropy import units as u
            from astropy.cosmology import Planck15, z_at_value

            with h5py.File(workdir / "inputs.h5", "r") as h:
                dmax = max((float(np.median(h[k]["dl"][:])) for k in h if not k.startswith("_")), default=None)
            if dmax:
                z_ev = float(z_at_value(Planck15.luminosity_distance, dmax * u.Mpc))
        images.append(rev.plot(post, plots / "rate_evolution.png", z_max_events=z_ev))
        rate_para = rev.report_paragraph(rev.summary(post), z_ev)
    if (workdir / "corner.png").exists():
        images.append(workdir / "corner.png")
    q = quantiles["H0"]
    s = summary.get("settings", {})
    fmt = lambda q: f"{q[2]:.1f} (+{q[3] - q[2]:.1f} / −{q[2] - q[1]:.1f}) km/s/Mpc, 90%: {q[0]:.1f}–{q[4]:.1f}"
    paras = [
        f"H<sub>0</sub> = <b>{q[2]:.1f} (+{q[3] - q[2]:.1f} / −{q[2] - q[1]:.1f}) km/s/Mpc</b> (median, 68%); "
        f"90%: {q[0]:.1f}–{q[4]:.1f}" + (" (reweighted to all the injections)" if use_rw and rw["target_inj_fraction"] >= 1 else "")
        + (f". Dark siren with a galaxy catalog, {len(events) or '?'} BBH events" if cat_settings else
           f". Spectral siren with {len(events) or '?'} BBH events")
        + (f" (runs {', '.join(selection['runs'])})" if selection.get("runs") else "") + ": the redshift comes from the "
        f"source-frame mass distribution ({model_name}, fitted together with H<sub>0</sub>)"
        + (" and from the galaxies along each line of sight" if cat_settings else "") + ", with the "
        "Madau–Dickinson rate evolution and flat ΛCDM (Ω<sub>m</sub> = 0.3065).",
        f"{summary['n_runs']} dynesty run(s) of {s.get('nlive', '?')} live points, {summary['n_samples']} posterior samples, "
        f"ln Z = {summary['log_evidence']:.2f} ± {summary['log_evidence_err']:.2f} (runs: "
        + ", ".join(f"{x:.1f}" for x in summary["run_log_evidences"]) + f"). PE samples per event: {s.get('pe_samples', '?')}; "
        f"fraction of the found injections used by the runs: {s.get('inj_fraction', '?')}.",
    ]
    if use_rw:
        qr = summary["quantiles"]["H0"]
        low = rw["ess_fraction"] < 0.1
        paras.append(
            f"The runs give H<sub>0</sub> = {fmt(qr)}. Their samples were reweighted to injection fraction "
            f"{rw['target_inj_fraction']:g} and {rw['target_pe_samples']} PE samples per event (importance weights "
            f"exp(ln L<sub>target</sub> − ln L<sub>runs</sub>), mean Δln L = {rw['dlnl_mean']:.2f}, scatter {rw['dlnl_sd']:.2f}): "
            f"effective sample size {rw['ess']:.0f} of {summary['n_samples']} ({rw['ess_fraction']:.0%}), "
            f"{rw['rejected']} sample(s) rejected by the target likelihood"
            + ("" if rw.get("runs_lnl_match", True) else "; <b>the ln L of the runs was not reproduced</b>") + ". "
            + ("<b>Low effective sample size: sample directly with --inj-fraction 1.</b>" if low else "Reliable reweighting."))
    if cat_settings:
        paras.append(_catalog_paragraph(cat_settings))
    if published:
        paras.append(f"Published ({published['ref']}): {published['median']} (+{published['plus']} / −{published['minus']}), "
                     f"90%: {published['lo90']}–{published['hi90']} km/s/Mpc.")
    plan = workdir / "plan.json"
    if plan.exists():
        paras.append("Strategy: " + json.loads(plan.read_text())["reason"] + ".")
    d = summary.get("diagnostics")
    if d:
        ok = d["neff_inj_min"] >= d["neff_inj_threshold"] and d["neff_pe_min"] >= d["neff_pe_threshold"]
        paras.append(
            f"Numerical stability of the runs over {d['n_points']} posterior draws: effective injections min {d['neff_inj_min']:.0f}, "
            f"median {d['neff_inj_median']:.0f} (threshold {d['neff_inj_threshold']}); smallest per-event effective PE "
            f"samples min {d['neff_pe_min']:.1f}, median {d['neff_pe_median_of_min']:.1f} (threshold {d['neff_pe_threshold']}); "
            f"lowest events: {', '.join(d['lowest_neff_pe_events'])}. " + ("OK." if ok else "<b>Below threshold: increase "
            "--inj-fraction or --pe-samples.</b>"))
    paras.append("The upper part of the H<sub>0</sub> interval depends on the prior bound (uniform 10–200 km/s/Mpc).")
    if rate_para:
        paras.append(rate_para)
    tables = [("Posterior (median, 68% and 90% intervals)" + (", reweighted" if use_rw else ""),
               table.to_html(index=False, float_format=lambda x: f"{x:.3g}"))]
    if not events.empty:
        tables.append(("Events", events.to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}")))
    Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
    write_simple_html_report(out_report_html, title="Hubble constant (" + ("dark siren, galaxy catalog" if cat_settings
                                                                     else "spectral siren") + ")", paragraphs=paras,
                             images=images, tables=tables)
    _log(f"report written to {out_report_html}")
    return table


def _matches_published_catalog(st: dict, published: dict) -> bool:
    """Whether the catalog has the settings of the published dark siren (those the paper gives: band, epsilon,
    nside, threshold map, redshift range, galaxy selection). Settings missing from a catalog built by an older version
    take their default; an unrecorded galaxy selection is the default one."""
    from .galaxy_catalog import DEFAULT_SETTINGS

    want = published.get("catalog") or {}
    for k, v in want.items():
        have = st.get(k, DEFAULT_SETTINGS.get(k, v if k == "galaxy_selection" else None))
        if isinstance(v, float) and have is not None:
            try:
                if float(have) != v:
                    return False
                continue
            except (TypeError, ValueError):
                return False
        if str(have) != str(v):
            return False
    return True


def _catalog_paragraph(st: dict) -> str:
    return (f"Galaxy catalog: {st.get('source', st.get('galaxies', '?'))}, band {st.get('band')}, luminosity weight "
            f"ε = {st.get('epsilon')}; HEALPix nside {st.get('nside')} (apparent-magnitude threshold: percentile "
            f"{st.get('mthr_percentile')} in nside-{st.get('nside_mthr')} pixels); galaxy redshift likelihood "
            f"{st.get('ptype')}, ±{st.get('numsigma')}σ, redshift grid {st.get('nintegration')}, in-catalog part from "
            f"z = {st.get('zmin', 0.0):g} to {st.get('zcut')}"
            + (f"; galaxy selection: {st['galaxy_selection']}" if st.get("galaxy_selection") else "")
            + ". Built with icarogw's pixelated catalog pipeline (galaxy_catalog mode).")


def _plan_injection_fraction(python: str, workdir: Path, mass_model: str, pe_samples: int, probe_points: int,
                             min_ess_fraction: float, npool: int = 1, prior_set: str = "gwtc4") -> float:
    """The injection fraction of the runs: the one already used in the work directory, else chosen by a probe."""
    settings = workdir / "run_settings.json"
    if settings.exists():
        f = float(json.loads(settings.read_text())["inj_fraction"])
        _log(f"injection fraction {f:g}, as in the runs already in {workdir}")
        return f
    probe_file = workdir / "probe.json"
    probe = json.loads(probe_file.read_text()) if probe_file.exists() else None
    if (not probe or probe.get("mass_model") != mass_model or probe.get("pe_samples") != pe_samples
            or probe.get("prior_set", "gwtc4") != prior_set):
        _log("probing the likelihood speed and accuracy on injection subsets")
        _run_runner(python, ["probe", "--workdir", str(workdir), "--mass-model", mass_model, "--prior-set", prior_set,
                             "--pe-samples",
                             str(pe_samples), "--npoints", str(probe_points), "--npool", str(npool), "--fractions"]
                    + [str(f) for f in PROBE_FRACTIONS])
        probe = json.loads(probe_file.read_text())
    f, reason = choose_injection_fraction(probe, min_ess_fraction=min_ess_fraction)
    (workdir / "plan.json").write_text(json.dumps(dict(inj_fraction=f, reason=reason, min_ess_fraction=min_ess_fraction),
                                                  indent=1))
    _log(f"strategy: {reason}")
    return f


def write_h0_slurm_chain(workdir: Path, python: str, stages: list[str], *, seeds: Iterable[int], mass_model: str,
                         prior_set: str, nlive: int, npool: int, naccept: int, pe_samples: int,
                         inj_fraction: float | str, probe_points: int, min_ess_fraction: float,
                         reweight_pe_samples: Optional[int], reweight_jobs: int, slurm_options: Iterable[str] = (),
                         env_setup: str = "") -> Path:
    """Slurm chain of the sampling stages in <workdir>/slurm (see `slurm.write_chain`): probe and plan (with
    inj_fraction "auto" and no runs yet), the runs as an array job (one task per seed, `npool` CPUs each), combine,
    and the reweighting to all the injections (an array of `reweight_jobs` chunks, then the merge; both do nothing
    when the runs already used all the injections)."""
    from .slurm import ARRAY_INDEX, Step, write_chain

    wd = str(workdir)
    seeds = list(dict.fromkeys(int(x) for x in seeds))
    settings = json.loads((workdir / "run_settings.json").read_text()) if (workdir / "run_settings.json").exists() else {}
    steps = []
    if "sample" in stages:
        if inj_fraction == "auto" and not settings:
            steps.append(Step("probe", ["probe", "--workdir", wd, "--mass-model", mass_model, "--prior-set", prior_set,
                                        "--pe-samples", str(pe_samples), "--npoints", str(probe_points), "--npool",
                                        str(npool), "--fractions"] + [str(f) for f in PROBE_FRACTIONS],
                              options=[f"--cpus-per-task={npool}"]))
            steps.append(Step("plan", ["plan", "--workdir", wd, "--min-ess-fraction", str(min_ess_fraction)],
                              options=["--cpus-per-task=1"]))
            frac = "plan"
        else:
            frac = str(float(settings["inj_fraction"]) if settings and inj_fraction == "auto" else float(inj_fraction))
        steps.append(Step("run", ["run", "--workdir", wd, "--seed", "${SEEDS[$SLURM_ARRAY_TASK_ID]}", "--mass-model",
                                  mass_model, "--prior-set", prior_set, "--nlive", str(nlive), "--npool", str(npool),
                                  "--naccept", str(naccept), "--pe-samples", str(pe_samples), "--inj-fraction", frac],
                          array=len(seeds), options=[f"--cpus-per-task={npool}"],
                          pre=["SEEDS=(" + " ".join(str(x) for x in seeds) + ")"]))
    if "combine" in stages:
        steps.append(Step("combine", ["combine", "--workdir", wd], options=["--cpus-per-task=1"],
                          pre=[f"rm -rf {wd}/reweight    # chunks of an earlier reweighting"]))
    if "reweight" in stages:
        target_pe = str(int(reweight_pe_samples or settings.get("pe_samples", pe_samples)))
        steps.append(Step("reweight", ["reweight", "--workdir", wd, "--chunk", ARRAY_INDEX, "--nchunks",
                                       str(reweight_jobs), "--target-inj-fraction", "1.0", "--target-pe-samples",
                                       target_pe, "--skip-if-done"], array=reweight_jobs,
                          options=["--cpus-per-task=1"]))
        steps.append(Step("reweight_merge", ["reweight-merge", "--workdir", wd, "--skip-if-done",
                                             "--target-pe-samples", target_pe], options=["--cpus-per-task=1"]))
    if not steps:
        raise ValueError("no sampling stage to run on Slurm (sample, combine, reweight)")
    cmd, env = _runner_command(python)
    return write_chain(workdir / "slurm", cmd[0], Path(cmd[1]), steps, logs=workdir / "logs", prefix="h0",
                       options=slurm_options, env_setup=env_setup, env=env)


def set_likelihood_settings(workdir: Path, neff_pe: Optional[float] = None, neff_inj: Optional[int] = None) -> dict:
    """Write the likelihood thresholds of a work directory (<workdir>/likelihood.json, read by every runner stage).
    Unset values keep those already recorded, else the runner defaults; changing them once runs exist is refused,
    since the runs, their combination and the reweighting must share one likelihood."""
    from .h0_icarogw import LIKELIHOOD_FILE, likelihood_settings

    workdir.mkdir(parents=True, exist_ok=True)
    f = workdir / LIKELIHOOD_FILE
    old = likelihood_settings(workdir)
    new = dict(old)
    if neff_pe is not None:
        if neff_pe <= 0:
            raise ValueError("--neff-pe must be positive")
        new["neff_pe"] = float(neff_pe)
    if neff_inj is not None:
        if neff_inj <= 0:
            raise ValueError("--neff-inj must be positive")
        new["neff_inj"] = int(neff_inj)
    runs = list((workdir / "result").glob("*_result.json")) if (workdir / "result").exists() else []
    if new != old and runs:
        raise ValueError(f"{workdir} has runs made with neff_pe={old['neff_pe']}, neff_inj={old['neff_inj']}; use "
                         "another work directory to change them")
    if new != old or not f.exists():
        f.write_text(json.dumps(new, indent=1))
    return new


def run_hubble_constant(
    stages: Iterable[str] = STAGES,
    workdir: str | Path = "hubble_constant_run",
    out_report_html: Optional[str | Path] = "hubble_constant.html",
    out_summary_tsv: Optional[str | Path] = "hubble_constant.tsv",
    sensitivity_release: str = H0_DEFAULT_RELEASE,
    sensitivity_file: Optional[str | Path] = None,
    catalogs: Optional[Iterable[str]] = None,
    far_threshold: float = 0.25,
    snr_threshold: float = 10.0,
    min_mass: float = 3.0,
    exclude: Iterable[str] = H0_DEFAULT_EXCLUDE,
    pe_cache: Optional[str | Path] = None,
    keep_pe_files: bool = False,
    seeds: Iterable[int] = (1,),
    parallel: int = 1,
    mass_model: str = H0_DEFAULT_MASS_MODEL,
    nlive: int = 100,
    npool: int = 4,
    naccept: int = 60,
    pe_samples: int = 1500,
    inj_fraction: float | str = "auto",
    min_ess_fraction: float = 0.5,
    probe_points: int = 30,
    reweight_pe_samples: Optional[int] = None,
    icarogw_python: Optional[str] = None,
    galaxy_catalog: Optional[str | Path] = None,
    executor: str = "local",
    slurm_options: Iterable[str] = (),
    slurm_env_setup: str = "",
    submit: bool = False,
    reweight_jobs: Optional[int] = None,
    neff_pe: Optional[float] = None,
    neff_inj: Optional[int] = None,
) -> Optional[pd.DataFrame]:
    """Spectral-siren H0 with the `mass_model` BBH mass distribution, or dark siren with `galaxy_catalog`; see the
    module docstring for the stages."""
    if mass_model not in MASS_MODELS:
        raise ValueError(f"Unknown mass model {mass_model!r}; choose from {', '.join(MASS_MODELS)}")
    stages = [s for s in STAGES if s in set(stages)]
    workdir = Path(workdir).expanduser().resolve()
    python = icarogw_python or sys.executable
    if "prepare" in stages:
        prepare_inputs(workdir, sensitivity_release, sensitivity_file, far_threshold, snr_threshold, min_mass,
                       exclude, Path(pe_cache).expanduser() if pe_cache else default_pe_cache(), keep_pe_files,
                       runs=catalog_runs(catalogs) if catalogs else None,
                       updates=_reg.update_catalogs(catalogs or ()),
                       galaxy_catalog=Path(galaxy_catalog) if galaxy_catalog else None)
    if any(st in stages for st in ("sample", "combine", "reweight")) or neff_pe is not None or neff_inj is not None:
        ls = set_likelihood_settings(workdir, neff_pe, neff_inj)
        _log(f"likelihood thresholds: {ls['neff_pe']:g} effective PE samples per event, "
             + (f"{ls['neff_inj']} effective injections" if ls["neff_inj"] else "4 x the events effective injections"))
    if executor == "slurm" and any(st in stages for st in ("sample", "combine", "reweight")):
        if not (workdir / "inputs.h5").exists():
            raise ValueError(f"{workdir / 'inputs.h5'} not found: run the prepare stage first")
        script = write_h0_slurm_chain(workdir, python, stages, seeds=seeds, mass_model=mass_model,
                                      prior_set=sensitivity_release, nlive=nlive, npool=npool, naccept=naccept,
                                      pe_samples=pe_samples, inj_fraction=inj_fraction, probe_points=probe_points,
                                      min_ess_fraction=min_ess_fraction, reweight_pe_samples=reweight_pe_samples,
                                      reweight_jobs=reweight_jobs or 16, slurm_options=slurm_options,
                                      env_setup=slurm_env_setup)
        _log(f"Slurm scripts written to {script.parent}; submit with {script}; then run the report stage")
        if submit:
            from .slurm import submit as _submit

            print(_submit(script), end="")
        return None
    elif executor != "local":
        raise ValueError("executor must be local or slurm")
    if any(st in stages for st in ("sample", "combine", "reweight")):
        _check_runner(python)
    ncpu = os.cpu_count() or 1
    if "sample" in stages:
        if not (workdir / "inputs.h5").exists():
            raise ValueError(f"{workdir / 'inputs.h5'} not found: run the prepare stage first")
        if inj_fraction == "auto":
            inj_fraction = _plan_injection_fraction(python, workdir, mass_model, pe_samples, probe_points,
                                                    min_ess_fraction, npool=max(1, min(ncpu, int(parallel) * int(npool))),
                                                    prior_set=sensitivity_release)
        seeds = list(dict.fromkeys(int(x) for x in seeds))
        run_args = ["--mass-model", mass_model, "--prior-set", sensitivity_release, "--nlive", str(nlive), "--npool", str(npool), "--naccept", str(naccept),
                    "--pe-samples", str(pe_samples), "--inj-fraction", str(float(inj_fraction))]
        parallel = max(1, min(int(parallel), len(seeds)))
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
    if "reweight" in stages:
        settings = json.loads((workdir / "run_settings.json").read_text()) if (workdir / "run_settings.json").exists() else {}
        target_pe = int(reweight_pe_samples or settings.get("pe_samples", pe_samples))
        if not settings:
            _log("reweight: no runs in the work directory, skipped")
        elif float(settings["inj_fraction"]) >= 1 and target_pe == int(settings["pe_samples"]):
            _log("reweight: the runs already use all the injections, nothing to do")
        else:
            import shutil

            shutil.rmtree(workdir / "reweight", ignore_errors=True)
            nchunks = max(1, min(ncpu, int(parallel) * int(npool)))
            jobs = [(f"reweight_{c:03d}", ["reweight", "--workdir", str(workdir), "--chunk", str(c), "--nchunks",
                                           str(nchunks), "--target-inj-fraction", "1.0",
                                           "--target-pe-samples", str(target_pe)]) for c in range(nchunks)]
            _log(f"reweighting the posterior to all the injections and {target_pe} PE samples per event, "
                 f"{nchunks} chunk(s)")
            run_jobs_parallel(python, workdir, jobs, nchunks)
            _run_runner(python, ["reweight-merge", "--workdir", str(workdir)])
    if "report" in stages and out_report_html:
        return write_h0_report(workdir, Path(out_report_html), Path(out_summary_tsv) if out_summary_tsv else None,
                               sensitivity_release)
    return None

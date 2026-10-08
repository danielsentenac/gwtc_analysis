"""Galaxy catalog for the dark-siren analysis: from a galaxy survey to icarogw's line-of-sight catalog.

Stages of `run_galaxy_catalog` (each restartable; the icarogw part runs in the icarogw environment through
`dark_catalog_icarogw.py`):

- ``galaxies``: the standard galaxy file (see `galaxies.py`): GLADE+ Ks band from VizieR (``source="glade-kband"``),
  or any catalog converted in chunks with a column mapping;
- ``shard``, ``pixels``, ``gather``, ``prepare``, ``init``, ``interpolate``, ``finish``: icarogw's pixelated catalog
  pipeline (HEALPix pixels, apparent-magnitude threshold, redshift grids, in-catalog line-of-sight interpolants,
  assembly of the catalog file);
- ``summary``: threshold map and completeness (icarogw environment);
- ``report``: HTML report.

The chunked stages (``pixels``, ``prepare``, ``interpolate``) run ``jobs`` chunks: as parallel processes with
``executor="local"``, or as Slurm array jobs with ``executor="slurm"`` (scripts written to <workdir>/slurm, chained
with ``afterok`` dependencies and submitted with ``submit=True``): the way to build catalogs deeper than GLADE+
(DES, Rubin) on a cluster such as CC-IN2P3.

The defaults are those of the GWTC-4.0 cosmology analysis (arXiv:2509.04348, Sections 2.2 and 3.2): GLADE+ K band,
HEALPix nside 64 for the catalog, apparent-magnitude threshold at the median in nside-32 pixels, luminosity
weighting epsilon = 1. The redshift-grid settings (`nintegration`, `numsigma`, `zmin`, `zcut`) and the form of the
galaxy redshift probability (`ptype`) are not given in the paper: the grid is one logarithmic grid fine enough for the
spectroscopic redshifts of GLADE+ (see DEFAULT_SETTINGS), and `ptype` is icarogw's Gaussian likelihood times a
uniform-in-comoving-volume prior (the paper tests it against the Gaussian taken as the posterior, gaussian_nocom, and
finds negligible differences).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, Optional

from .report import write_simple_html_report

GC_STAGES = ("galaxies", "shard", "pixels", "gather", "prepare", "init", "interpolate", "finish", "summary", "report")
ICAROGW_STAGES = GC_STAGES[1:-1]
CHUNKED = ("pixels", "prepare", "interpolate")
ASSEMBLY = ("init", "finish", "summary")       # single jobs holding the whole catalog in memory
RUNNER = Path(__file__).with_name("dark_catalog_icarogw.py")

# the redshift grid: one logarithmic grid for every pixel (icarogw's fixed-array path, as its LVK job templates
# write it). icarogw's adaptive grid (an integer) refines around every galaxy, and with the spectroscopic redshifts of
# GLADE+ (sigma_z of a few 1e-4) the merged grid grows to tens of thousands of redshifts, too slow to build and too
# large to load. 5000 points from 1e-4 to zcut = 0.5: a relative step of 0.17%, 8.5e-5 at z = 0.05.
DEFAULT_SETTINGS = dict(nside=64, nside_mthr=32, mthr_percentile=50.0, epsilon=1.0, nintegration="logspace:0.0001:5000",
                        numsigma=3, zmin=0.0, zcut=0.5, ptype="gaussian", cosmo_zmax=10.0, nshards=1024)
# icarogw's galaxy redshift probabilities (icarogw.catalog.EM_likelihood_prior_differential_volume): gaussian =
# Gaussian likelihood x uniform-in-comoving-volume prior, normalized; gaussian_nocom = the Gaussian itself;
# uniform = uniform in comoving volume within +-numsigma sigma
PTYPES = ("gaussian", "gaussian_nocom", "uniform")


def _log(msg: str) -> None:
    print(f"[galaxy_catalog] {msg}", flush=True)


def catalog_settings(band: str, **kw) -> dict:
    """Settings of the icarogw catalog construction (DEFAULT_SETTINGS overridden by `kw`)."""
    s = dict(DEFAULT_SETTINGS)
    s.update({k: v for k, v in kw.items() if v is not None})
    s["band"] = band
    s["grouping"] = band
    s["subgrouping"] = f"eps_{float(s['epsilon']):g}"
    s["outfile"] = f"catalog_{band}_nside{int(s['nside'])}_eps{float(s['epsilon']):g}.hdf5"
    if int(s["nside_mthr"]) > int(s["nside"]):
        raise ValueError("nside_mthr must not exceed nside")
    if not 0 < float(s["mthr_percentile"]) <= 100:
        raise ValueError("mthr_percentile must be in (0, 100]")
    if s["ptype"] not in PTYPES:
        raise ValueError(f"ptype must be one of {', '.join(PTYPES)}")
    s["zmin"] = float(s["zmin"])
    if not 0 <= s["zmin"] < float(s["zcut"]):
        raise ValueError("zmin must be in [0, zcut)")
    if s["zmin"] > 0 and not str(s["nintegration"]).startswith("logspace:"):
        raise ValueError("zmin needs a logarithmic redshift grid (nintegration logspace:ZMIN:N): with the adaptive "
                         "grid, redshifts below zmin would count as an empty catalog instead of out of its range")
    return s


def _runner_command(python: str) -> tuple[list[str], dict]:
    from .hubble_constant import _runner_command as h0_cmd

    cmd, env = h0_cmd(python)
    env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    return [cmd[0], str(RUNNER)], env


def _check_runner(python: str) -> None:
    cmd, env = _runner_command(python)
    r = subprocess.run([cmd[0], "-c", "import icarogw, healpy, h5py"], env=env, capture_output=True, text=True,
                       cwd=str(Path.home()))
    if r.returncode != 0:
        raise ValueError(f"icarogw cannot be imported by {cmd[0]} (pass --icarogw-python): "
                         f"{r.stderr.strip().splitlines()[-1] if r.stderr.strip() else ''}")


def _run(python: str, args: list[str], log: Optional[Path] = None) -> None:
    cmd, env = _runner_command(python)
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        with open(log, "a") as fh:
            r = subprocess.run(cmd + args, env=env, stdout=fh, stderr=subprocess.STDOUT)
    else:
        r = subprocess.run(cmd + args, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"catalog stage {args[0]} failed (exit code {r.returncode})"
                           + (f"; see {log}" if log else ""))


def _stage_args(stage: str, workdir: Path, chunk: Optional[int] = None, nchunks: Optional[int] = None,
                galaxies: Optional[Path] = None) -> list[str]:
    a = [stage, "--workdir", str(workdir)]
    if stage == "shard":
        a += ["--galaxies", str(galaxies), "--settings", str(workdir / "settings_input.json")]
    if chunk is not None:
        a += ["--chunk", str(chunk), "--nchunks", str(nchunks)]
    elif stage == "gather":
        a += ["--nchunks", str(nchunks)]
    return a


# ---------------------------------------------------------------------------
# Slurm
# ---------------------------------------------------------------------------

def write_slurm_scripts(workdir: Path, python: str, stages: list[str], jobs: int, galaxies: Path,
                        slurm_options: Iterable[str] = (), env_setup: str = "",
                        assembly_options: Iterable[str] = ()) -> Path:
    """One sbatch script per stage in <workdir>/slurm (array jobs for the chunked stages) and submit.sh, which
    submits them chained by afterok dependencies (see `slurm.write_chain`). `assembly_options` are added (after
    `slurm_options`, so they win) to the init, finish and summary jobs, which hold the whole catalog in memory."""
    from .slurm import ARRAY_INDEX, Step, write_chain

    steps = []
    for st in stages:
        chunked = st in CHUNKED
        args = _stage_args(st, workdir, chunk=0 if chunked else None, nchunks=jobs, galaxies=galaxies)
        if chunked:
            args[args.index("--chunk") + 1] = ARRAY_INDEX
        steps.append(Step(st, args, array=jobs if chunked else None,
                          options=list(assembly_options) if st in ASSEMBLY else []))
    _, env = _runner_command(python)
    return write_chain(workdir / "slurm", python, RUNNER, steps, logs=workdir / "logs", prefix="gwcat",
                       options=slurm_options, env_setup=env_setup, env=env)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

def write_catalog_report(workdir: Path, out_report_html: Path) -> None:
    summ = json.loads((workdir / "summary.json").read_text())
    s = summ["settings"]
    gal_attrs = {}
    gfile = Path(s.get("galaxies", ""))
    if gfile.exists():
        from .galaxies import galaxy_summary

        gal_attrs = galaxy_summary(gfile)
    comp = ", ".join(f"z = {z:.2f}: {m:.0%} ({lo:.0%}–{hi:.0%})" for z, m, lo, hi in
                     zip(summ["completeness_z"], summ["completeness_median"], summ["completeness_p20"],
                         summ["completeness_p80"]))
    paras = [
        f"Galaxy catalog for the dark-siren analysis: <b>{summ['n_galaxies']:,}</b> galaxies"
        + (f" from {gal_attrs.get('source')}" if gal_attrs.get("source") else "")
        + (f" (selection: {s['galaxy_selection']})" if s.get("galaxy_selection") else "")
        + f", band {s['band']} (median redshift {gal_attrs.get('z_median', float('nan')):.3f}, median redshift "
        f"uncertainty {gal_attrs.get('sigmaz_median', float('nan')):.3f}); {summ['n_galaxies_in_catalog_term']:,} are "
        "brighter than the apparent-magnitude threshold of their pixel and enter the in-catalog term.",
        f"HEALPix nside {s['nside']}: {summ['n_filled_pixels']:,} pixels hold galaxies ({summ['sky_fraction']:.0%} of the "
        f"sky). Apparent-magnitude threshold: percentile {s['mthr_percentile']:g} of the magnitudes in nside-"
        f"{s['nside_mthr']} pixels, median {summ['mthr_median']:.2f} (10–90%: {summ['mthr_p10_p90'][0]:.2f}–"
        f"{summ['mthr_p10_p90'][1]:.2f}).",
        f"Line-of-sight interpolants: luminosity weight ε = {s['epsilon']:g}, galaxy redshift probability {s['ptype']} "
        f"(±{s['numsigma']}σ), redshift grid {s['nintegration']}, in-catalog term from z = {s.get('zmin', 0):g} to "
        f"{s['zcut']}: {summ['n_redshifts']} "
        f"redshifts × {summ['n_moc_pixels']:,} sky pixels, {summ['file_mb']:.0f} MB "
        f"(<code>{workdir / s['outfile']}</code>).",
        "Completeness (fraction of the luminosity-weighted galaxy density in the catalog, median and 20–80% over the "
        f"pixels): {comp}.",
        "Built with icarogw's pixelated catalog pipeline (the functions reviewed for the LVK analyses); use it with "
        f"<code>gwtc_analysis hubble_constant --galaxy-catalog {workdir / s['outfile']}</code>.",
    ]
    images = [p for p in (workdir / "plots" / "mthr_map.png", workdir / "plots" / "completeness.png") if p.exists()]
    Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
    write_simple_html_report(out_report_html, title=f"Galaxy catalog ({s['band']})", paragraphs=paras, images=images)
    _log(f"report written to {out_report_html}")


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def run_galaxy_catalog(
    stages: Iterable[str] = GC_STAGES,
    workdir: str | Path = "galaxy_catalog_run",
    source: Optional[str] = "glade-kband",
    galaxies: Optional[str | Path] = None,
    input_catalog: Optional[str | Path] = None,
    columns: Optional[dict] = None,
    band: Optional[str] = None,
    angle_unit: str = "deg",
    sigmaz: Optional[float] = None,
    sigmaz_relative: bool = False,
    where: Optional[str] = None,
    glade_types: Optional[str | Iterable[str]] = None,
    glade_redshift: Optional[str] = None,
    glade_sigmaz: Optional[str] = None,
    input_format: Optional[str] = None,
    nside: Optional[int] = None,
    nside_mthr: Optional[int] = None,
    mthr_percentile: Optional[float] = None,
    epsilon: Optional[float] = None,
    nintegration: Optional[int | str] = None,
    numsigma: Optional[int] = None,
    zmin: Optional[float] = None,
    zcut: Optional[float] = None,
    ptype: Optional[str] = None,
    nshards: Optional[int] = None,
    jobs: int = 4,
    executor: str = "local",
    slurm_options: Iterable[str] = (),
    slurm_assembly_options: Iterable[str] = (),
    slurm_env_setup: str = "",
    submit: bool = False,
    icarogw_python: Optional[str] = None,
    out_report_html: Optional[str | Path] = "galaxy_catalog.html",
) -> Optional[Path]:
    """Build the icarogw line-of-sight galaxy catalog; returns the catalog file (when built locally)."""
    from . import galaxies as gx

    stages = [s for s in GC_STAGES if s in set(stages)]
    workdir = Path(workdir).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    python = str(Path(icarogw_python).expanduser()) if icarogw_python else sys.executable
    if executor not in ("local", "slurm"):
        raise ValueError("executor must be local or slurm")
    jobs = max(1, int(jobs))

    # the standard galaxy file
    gal = Path(galaxies).expanduser().resolve() if galaxies else workdir / "galaxies.h5"
    if "galaxies" in stages and not galaxies:
        if input_catalog:
            if not (columns and band):
                raise ValueError("A converted catalog needs --columns and --band")
            gx.convert_catalog(input_catalog, gal, band=band, columns=columns, angle_unit=angle_unit, sigmaz=sigmaz,
                               sigmaz_relative=sigmaz_relative, where=where, fmt=input_format)
        elif source == "glade-kband":
            sel = gx.glade_selection(types=glade_types, redshift=glade_redshift, sigmaz=glade_sigmaz,
                                     sigmaz_const=sigmaz, sigmaz_relative=sigmaz_relative or None, where=where)
            if gal.exists():
                import h5py

                with h5py.File(gal, "r") as h:
                    prev = json.loads(h.attrs["glade_selection"]) if "glade_selection" in h.attrs else \
                        dict(gx.GLADE_DEFAULT_SELECTION)       # files written before the selection was recorded
                prev["types"] = tuple(prev["types"])
                if prev != sel:
                    raise ValueError(f"{gal} was made with another GLADE+ selection ({gx.glade_selection_text(prev)}); "
                                     "use another work directory, or remove it")
                _log(f"{gal} exists, kept")
            else:
                gx.fetch_glade_kband(gal, selection=sel)
        else:
            raise ValueError(f"Unknown galaxy source {source!r} (glade-kband, or --input-catalog)")
    icarogw_stages = [s for s in stages if s in ICAROGW_STAGES]
    if icarogw_stages:
        if not gal.exists():
            raise ValueError(f"{gal} not found: run the galaxies stage, or pass --galaxies")
        import h5py

        with h5py.File(gal, "r") as h:
            file_band = str(h.attrs["band"])
        st = catalog_settings(band or file_band, nside=nside, nside_mthr=nside_mthr, mthr_percentile=mthr_percentile,
                              epsilon=epsilon, nintegration=nintegration, numsigma=numsigma, zmin=zmin, zcut=zcut,
                              ptype=ptype, nshards=nshards)
        summ = gx.galaxy_summary(gal)
        st["source"] = summ["source"]
        st["galaxy_selection"] = summ.get("selection", "")
        st["galaxies"] = str(gal)
        old = workdir / "catalog_settings.json"
        if old.exists():
            prev = json.loads(old.read_text())
            # settings added in later versions take their default in an older catalog_settings.json
            diff = {k: (prev.get(k, DEFAULT_SETTINGS.get(k)), v) for k, v in st.items()
                    if k not in ("galaxies", "source", "nshards") and not (k == "galaxy_selection" and k not in prev)
                    and str(prev.get(k, DEFAULT_SETTINGS.get(k))) != str(v)}
            if diff:
                raise ValueError(f"{workdir} holds a catalog built with other settings {diff}; use another work "
                                 "directory")
        (workdir / "settings_input.json").write_text(json.dumps(st, indent=1))
        if executor == "slurm":
            script = write_slurm_scripts(workdir, python, icarogw_stages, jobs, gal, slurm_options, slurm_env_setup,
                                         slurm_assembly_options)
            _log(f"Slurm scripts written to {script.parent}; submit with {script}")
            if submit:
                from .slurm import submit as _submit

                print(_submit(script), end="")
            return None
        _check_runner(python)
        logs = workdir / "logs"
        for st_name in icarogw_stages:
            if st_name in CHUNKED:
                _log(f"{st_name}: {jobs} chunk(s) in parallel")
                with ThreadPoolExecutor(jobs) as ex:
                    futs = {ex.submit(_run, python, _stage_args(st_name, workdir, c, jobs, gal),
                                      logs / f"{st_name}_{c:03d}.log"): c for c in range(jobs)}
                    errors = []
                    for fu in as_completed(futs):
                        try:
                            fu.result()
                        except Exception as e:
                            errors.append(str(e))
                    if errors:
                        raise RuntimeError(f"{st_name}: {len(errors)} chunk(s) failed: {errors[0]}")
            else:
                _log(st_name)
                _run(python, _stage_args(st_name, workdir, nchunks=jobs, galaxies=gal), logs / f"{st_name}.log")
    if "report" in stages and out_report_html:
        if not (workdir / "summary.json").exists():
            raise ValueError(f"{workdir / 'summary.json'} not found: run the summary stage first")
        write_catalog_report(workdir, Path(out_report_html))
    st_file = workdir / "catalog_settings.json"
    if st_file.exists():
        out = workdir / json.loads(st_file.read_text())["outfile"]
        return out if out.exists() else None
    return None

"""``scalar_cloud`` mode: is a binary black hole better fitted with a scalar-field cloud around it?

The physics, the fits and the cluster jobs are in the self-contained runner ``scalar_cloud.py`` (the waveform
IMRPhenomXPHM_Scalar of Roy et al., arXiv:2510.17967). This module prepares the work directory, runs the steps
locally or writes them as Slurm jobs, and writes the report:

- **best fit**: the vacuum model and the cloud model fitted to the event's strain from the catalog posterior;
  the statistic Lambda = 2 (ln L_cloud - ln L_vacuum) and the whitened strain against the catalog, vacuum and
  cloud waveforms;
- **background**: the same procedure on the vacuum best fit injected into off-source noise around the event: the
  distribution of Lambda without a cloud, hence the p-value of the event's Lambda;
- **sampling** (optional, on a cluster): nested sampling of the vacuum and cloud models with the paper's settings,
  hence the Bayes factor.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import numpy as np

from . import scalar_cloud as sc

INK, INK2, GRID = "#1f2328", "#57606a", "#d0d7de"
COLORS = dict(data="#8c959f", catalog="#7a5195", vacuum="#2a78d6", cloud="#d1495b")


def _log(msg: str) -> None:
    print(f"[scalar_cloud] {msg}", flush=True)


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, color=GRID, lw=0.5, alpha=0.6)


# ---------------------------------------------------------------------------
# whitened strain
# ---------------------------------------------------------------------------
def _whiten(ifo, h_fd: np.ndarray, f_lo: float = 20.0, f_hi: float = 350.0) -> np.ndarray:
    """Whitened, band-passed time series of a frequency-domain series of `ifo`."""
    f = ifo.frequency_array
    asd = ifo.amplitude_spectral_density_array
    w = np.zeros_like(h_fd)
    band = (f >= f_lo) & (f <= f_hi) & np.isfinite(asd) & (asd > 0)
    w[band] = h_fd[band] / asd[band]
    # a smooth edge so that the band limits do not ring
    taper = np.ones_like(f)
    for lo, hi in ((f_lo, f_lo + 5), (f_hi - 20, f_hi)):
        k = (f >= lo) & (f <= hi)
        x = (f[k] - lo) / (hi - lo)
        taper[k] = 0.5 * (1 - np.cos(np.pi * x)) if lo == f_lo else 0.5 * (1 + np.cos(np.pi * x))
    return np.fft.irfft(w * taper, n=len(ifo.time_array))


def plot_whitened(an, curves: list[tuple[str, dict, str]], out_png: Path,
                  windows: tuple[tuple[float, float, float], ...] = ((-5.0, -3.0, 60.0), (-0.5, 0.06, 350.0))
                  ) -> Path:
    """Whitened strain of H1 and L1 with the waveforms of `curves` (label, parameters, colour), in an early
    inspiral window and around the merger (times from the merger of the first curve; each window band-passed from
    20 Hz to its own upper frequency, so that the early inspiral is not buried in high-frequency noise), and a last
    row with the difference between the last and the first waveform (where a cloud acts)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dets = [ifo for ifo in an.ifos if ifo.name != "V1"] or list(an.ifos)
    nrow, ncol = len(dets) + 1, len(windows)
    fig, axes = plt.subplots(nrow, ncol, figsize=(12, 2.4 * nrow), squeeze=False,
                             gridspec_kw=dict(width_ratios=[w[1] - w[0] for w in windows]))
    p_ref = curves[0][1]
    for r, ifo in enumerate(dets):
        t = ifo.time_array
        t_merger = p_ref["geocent_time"] + ifo.time_delay_from_geocenter(p_ref["ra"], p_ref["dec"],
                                                                         p_ref["geocent_time"])
        resp = [ifo.get_detector_response(an.waveform_generator.frequency_domain_strain(p), p) for _, p, _ in curves]
        hw = {}
        for c, (a, b, f_hi) in enumerate(windows):
            ax = axes[r, c]
            k = (t >= t_merger + a) & (t <= t_merger + b)
            d = _whiten(ifo, ifo.frequency_domain_strain, f_hi=f_hi)
            hw[c] = [_whiten(ifo, h, f_hi=f_hi) for h in resp]
            ax.plot(t[k] - t_merger, d[k], color=COLORS["data"], lw=0.7,
                    label=f"{ifo.name} data (whitened, 20–{f_hi:g} Hz)")
            for (label, p, color), h in zip(curves, hw[c]):
                ax.plot(t[k] - t_merger, h[k], color=color, lw=1.3, label=label)
            _style(ax)
            if c == 0:
                ax.set_ylabel(f"{ifo.name}", color=INK2, fontsize=9)
            if r == 0:
                ax.legend(frameon=False, fontsize=7.5, loc="upper left", labelcolor=INK2)
        if r == len(dets) - 1:
            for c, (a, b, _) in enumerate(windows):
                ax = axes[-1, c]
                k = (t >= t_merger + a) & (t <= t_merger + b)
                ax.plot(t[k] - t_merger, hw[c][-1][k] - hw[c][0][k], color=curves[-1][2], lw=1.2,
                        label=f"{ifo.name}: {curves[-1][0]} − {curves[0][0]}")
                _style(ax)
                ax.set_xlabel("time from merger (s)", color=INK2, fontsize=9)
                if c == 0:
                    ax.set_ylabel("difference", color=INK2, fontsize=9)
                    ax.legend(frameon=False, fontsize=7.5, loc="upper left", labelcolor=INK2)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
# dephasing, profile, background
# ---------------------------------------------------------------------------
def _masses(p: dict) -> tuple[float, float, float]:
    from bilby.gw.conversion import chirp_mass_and_mass_ratio_to_total_mass

    m_det = chirp_mass_and_mass_ratio_to_total_mass(p["chirp_mass"], p["mass_ratio"])
    return m_det, sc._source_frame_total_mass(m_det, p["luminosity_distance"]), p["mass_ratio"]


def plot_dephasing(an, p_cloud: dict, p_vac: dict, out_png: Path) -> Path:
    """Delta Psi_22(f) of the best cloud fit, against the cumulative SNR^2 of the vacuum fit."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    f = np.geomspace(an.s["f_low"], 400, 400)
    m_det, m_src, q = _masses(p_cloud)
    d = sc.dephasing_22(f, m_det, m_src, q, p_cloud["rho_phi"], p_cloud["alpha_phi"])
    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    ax.plot(f, -d, color=COLORS["cloud"], lw=2, label="|ΔΨ₂₂(f)| of the best cloud fit (rad)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-3, max(1.0, -d.min() * 2))
    _style(ax)
    ax.set_xlabel("GW frequency of the (2, 2) mode (Hz)", color=INK2, fontsize=9)
    ax.set_ylabel("dephasing (rad)", color=INK2, fontsize=9)
    ax2 = ax.twinx()
    pol = an.waveform_generator.frequency_domain_strain(p_vac)
    cum = 0
    fr = an.ifos[0].frequency_array
    snr2 = np.zeros_like(fr)
    for ifo in an.ifos:
        h = ifo.get_detector_response(pol, p_vac)
        psd = ifo.power_spectral_density_array
        ok = ifo.frequency_mask & np.isfinite(psd) & (psd > 0)
        snr2[ok] += 4 * np.abs(h[ok]) ** 2 / psd[ok] / an.s["duration"]
    cum = np.cumsum(snr2) / snr2.sum()
    ax2.plot(fr, cum, color=COLORS["vacuum"], lw=1.2, ls="--", label="cumulative SNR² of the vacuum fit")
    ax2.set_ylim(0, 1.05)
    ax.set_xlim(15, 500)
    ax2.tick_params(colors=INK2, labelsize=9)
    ax2.set_ylabel("cumulative fraction of SNR²", color=INK2, fontsize=9)
    for s in ("top",):
        ax2.spines[s].set_visible(False)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8, loc="lower left", labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png


def plot_profile(points: list[dict], ln_vac: float, out_png: Path) -> Path:
    """ln L - ln L_vacuum over the (alpha, dephasing at f_low) grid, the binary re-fitted at each point."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    alphas = sorted({p["alpha"] for p in points})
    targets = sorted({p["target_dephasing"] for p in points})
    z = np.full((len(targets), len(alphas)), np.nan)
    for p in points:
        z[targets.index(p["target_dephasing"]), alphas.index(p["alpha"])] = p["log_likelihood_ratio"] - ln_vac
    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    vmax = np.nanmax(np.abs(z))
    im = ax.imshow(z, origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    for i in range(len(targets)):
        for j in range(len(alphas)):
            if np.isfinite(z[i, j]):
                ax.text(j, i, f"{z[i, j]:+.1f}", ha="center", va="center", fontsize=8, color=INK)
    ax.set_xticks(range(len(alphas)), [f"{a:g}" for a in alphas])
    ax.set_yticks(range(len(targets)), [f"{t:g}" for t in targets])
    ax.set_xlabel("α = M m_φ", color=INK2, fontsize=9)
    ax.set_ylabel("|ΔΨ₂₂| at f_low (rad): sets ρ̄", color=INK2, fontsize=9)
    ax.tick_params(colors=INK2, labelsize=9)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("ln L − ln L_vacuum (binary re-fitted)", color=INK2, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png


def plot_background(lams: np.ndarray, lam_event: float, out_png: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    bins = np.linspace(0, max(lams.max(), lam_event) * 1.1 + 1e-9, 30)
    ax.hist(lams, bins=bins, color=COLORS["vacuum"], alpha=0.75, label=f"background: {len(lams)} vacuum injections")
    ax.axvline(lam_event, color=COLORS["cloud"], lw=2, label=f"GW event: Λ = {lam_event:.2f}")
    _style(ax)
    ax.set_xlabel("Λ = 2 (ln L_cloud − ln L_vacuum)", color=INK2, fontsize=9)
    ax.set_ylabel("number", color=INK2, fontsize=9)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png


def background_p_value(lams: np.ndarray, lam_event: float) -> tuple[float, int]:
    """Fraction of the background with Lambda >= the event's ((k + 1) / (n + 1), conservative), and k."""
    k = int(np.sum(lams >= lam_event))
    return (k + 1) / (len(lams) + 1), k


def lnl_by_frequency(an, p: dict) -> tuple[np.ndarray, np.ndarray]:
    """(f, cumulative ln L ratio up to f) of the waveform `p`: ln L = sum_f Re<d, h> - <h, h>/2 over the bins."""
    pol = an.waveform_generator.frequency_domain_strain(p)
    f = an.ifos[0].frequency_array
    contrib = np.zeros_like(f)
    for ifo in an.ifos:
        h = ifo.get_detector_response(pol, p)
        d = ifo.frequency_domain_strain
        psd = ifo.power_spectral_density_array
        ok = ifo.frequency_mask & np.isfinite(psd) & (psd > 0)
        w = np.zeros_like(f)
        w[ok] = 4.0 / an.s["duration"] / psd[ok]
        contrib += w * (np.real(np.conj(d) * h) - 0.5 * np.abs(h) ** 2)
    return f, np.cumsum(contrib)


LNL_BANDS = (20.0, 30.0, 35.0, 60.0, None)


def lnl_gain_by_band(an, p_cloud: dict, p_vac: dict, edges=LNL_BANDS):
    """ln L_cloud - ln L_vacuum per detector and frequency band (a table: where the improvement comes from)."""
    import pandas as pd

    f = an.ifos[0].frequency_array
    gain = {ifo.name: np.zeros_like(f) for ifo in an.ifos}
    for sign, p in ((1.0, p_cloud), (-1.0, p_vac)):
        pol = an.waveform_generator.frequency_domain_strain(p)
        for ifo in an.ifos:
            h = ifo.get_detector_response(pol, p)
            psd = ifo.power_spectral_density_array
            ok = ifo.frequency_mask & np.isfinite(psd) & (psd > 0)
            w = np.zeros_like(f)
            w[ok] = 4.0 / an.s["duration"] / psd[ok]
            gain[ifo.name] += sign * w * (np.real(np.conj(ifo.frequency_domain_strain) * h) - 0.5 * np.abs(h) ** 2)
    hi = [e if e is not None else np.inf for e in edges[1:]]
    cols = [f"{lo:g}–{h:g} Hz" if np.isfinite(h) else f"> {lo:g} Hz" for lo, h in zip(edges[:-1], hi)]
    rows = []
    for det, g in gain.items():
        rows.append(dict(detector=det, **{c: g[(f >= lo) & (f < h)].sum() for c, lo, h in zip(cols, edges[:-1], hi)},
                         total=g.sum()))
    tab = pd.DataFrame(rows)
    tab.loc[len(tab)] = dict(detector="all", **{c: tab[c].sum() for c in cols + ["total"]})
    return tab


def plot_lnl_by_frequency(an, p_cloud: dict, p_vac: dict, out_png: Path) -> Path:
    """Cumulative ln L_cloud - ln L_vacuum against frequency: which part of the signal prefers the cloud."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    f, lc = lnl_by_frequency(an, p_cloud)
    _, lv = lnl_by_frequency(an, p_vac)
    k = (f >= an.s["f_low"] * 0.95) & (f <= 500)
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    ax.plot(f[k], (lc - lv)[k], color=COLORS["cloud"], lw=2)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xscale("log")
    _style(ax)
    ax.set_xlabel("frequency (Hz)", color=INK2, fontsize=9)
    ax.set_ylabel("cumulative ln L_cloud − ln L_vacuum", color=INK2, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
# the mode
# ---------------------------------------------------------------------------
STAGES = ("prepare", "fit", "background", "sample", "seeded", "report")
SAMPLE_MODELS = ("vacuum", "agnostic", "superradiance")


def _sample_label(model: str, tau_d: float) -> str:
    return model if model != "superradiance" else f"superradiance_tau{tau_d:.0e}".replace("+", "")


def _profile_worker(args):
    wd, i = args
    an = sc.Analysis(wd)
    vac = json.loads((Path(wd) / "vacuum.json").read_text())
    sc._dump(Path(wd) / "profile" / f"point_{i:03d}.json", sc.profile_point(an, vac, i))
    return i


def _background_worker(args):
    wd, k = args
    sc.main(["background", "--workdir", str(wd), "--task", str(k)])
    return k


def _pool_map(fn, items, jobs: int):
    if jobs <= 1:
        return [fn(x) for x in items]
    from concurrent.futures import ProcessPoolExecutor

    with ProcessPoolExecutor(max_workers=jobs) as ex:
        return list(ex.map(fn, items))


def write_slurm_chains(workdir: Path, python: str, stages: list[str], n_background: int, sample_models, tau_d: float,
                       npool: int, slurm_options=(), env_setup: str = "") -> list[Path]:
    """Slurm chains in <workdir>/slurm (fit, then background) and <workdir>/slurm_sample (one array task per
    model), each with its submit.sh. The jobs run the self-contained runner scalar_cloud.py."""
    from .slurm import ARRAY_INDEX, Step, write_chain

    wd = str(workdir)
    runner = Path(sc.__file__)
    lib = Path(python).expanduser().parent.parent / "lib"
    env = {"LD_LIBRARY_PATH": str(lib)} if lib.is_dir() else None
    pre = ["export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1"]
    out = []
    steps = []
    if "fit" in stages:
        n_grid = len(sc.ALPHAS) * len(sc.DEPHASING_TARGETS)
        steps += [Step("vacuum", ["vacuum", "--workdir", wd], options=["--cpus-per-task=1", "--mem=3G",
                                                                        "--time=06:00:00"], pre=pre),
                  Step("profile", ["profile", "--workdir", wd, "--task", ARRAY_INDEX], array=n_grid,
                       options=["--cpus-per-task=1", "--mem=3G", "--time=08:00:00"], pre=pre),
                  Step("polish", ["polish", "--workdir", wd], options=["--cpus-per-task=1", "--mem=3G",
                                                                      "--time=12:00:00"], pre=pre)]
    if "background" in stages:
        steps.append(Step("background", ["background", "--workdir", wd, "--task", ARRAY_INDEX], array=n_background,
                          options=["--cpus-per-task=1", "--mem=3G", "--time=1-00:00:00"], pre=pre))
    if steps:
        out.append(write_chain(workdir / "slurm", python, runner, steps, logs=workdir / "logs", prefix="sc",
                               options=slurm_options, env_setup=env_setup, env=env))
    if "sample" in stages:
        models = list(sample_models)
        step = Step("sample", ["sample", "--workdir", wd, "--model", "${MODELS[$SLURM_ARRAY_TASK_ID]}", "--tau-d",
                               f"{tau_d:g}", "--npool", str(npool)], array=len(models),
                    options=[f"--cpus-per-task={npool}", f"--mem={max(8, int(1.5 * npool))}G", "--time=4-00:00:00"],
                    pre=pre + ["MODELS=(" + " ".join(models) + ")"])
        out.append(write_chain(workdir / "slurm_sample", python, runner, [step], logs=workdir / "logs",
                               prefix="sc_sample", options=slurm_options, env_setup=env_setup, env=env))
    if "seeded" in stages:
        steps = [Step("seeds", ["seeds", "--workdir", wd, "--result", str(_seed_result(workdir, tau_d))],
                      options=["--cpus-per-task=1", "--mem=4G", "--time=06:00:00"], pre=pre),
                 Step("seeded", ["seeded", "--workdir", wd, "--task", "$((SLURM_ARRAY_TASK_ID-1))"],
                      array=n_background + 1, options=["--cpus-per-task=1", "--mem=3G", "--time=1-00:00:00"],
                      pre=pre)]
        out.append(write_chain(workdir / "slurm_seeded", python, runner, steps, logs=workdir / "logs",
                               prefix="sc_seeded", options=slurm_options, env_setup=env_setup, env=env))
    return out


def _seed_result(workdir: Path, tau_d: float) -> Path:
    """The sampling result whose best samples seed the cloud fits: the superradiance prior's."""
    return Path(workdir) / "sampling" / f"{_sample_label('superradiance', tau_d)}_result.hdf5"


def _seeded_worker(args):
    wd, k = args
    sc.main(["seeded", "--workdir", wd, "--task", str(k)])


def run_scalar_cloud(event: str = "GW190728_064510", workdir: str | Path = "scalar_cloud",
                     stages=("prepare", "fit", "report"), pe_file: Optional[str] = None,
                     pe_label: str = "C01:IMRPhenomXPHM", sampling_frequency: float = 2048.0,
                     n_background: int = 150, jobs: int = 1, executor: str = "local", python: Optional[str] = None,
                     slurm_options=(), slurm_env_setup: str = "", submit: bool = False,
                     sample_models=SAMPLE_MODELS, tau_d: float = 1e6, npool: int = 32,
                     out_report_html: Optional[str] = "scalar_cloud.html",
                     out_summary_tsv: Optional[str] = "scalar_cloud.tsv") -> dict:
    """Run the stages of the scalar_cloud mode in `workdir` (see the module docstring)."""
    import sys

    wd = Path(workdir).expanduser()
    stages = [s for s in STAGES if s in stages]
    bad = [m for m in sample_models if m not in SAMPLE_MODELS]
    if bad:
        raise ValueError(f"unknown sampling model(s) {bad}; choose among {SAMPLE_MODELS}")
    if "prepare" in stages:
        if pe_file is None:
            from . import hubble_constant as hc
            from .counterpart import _event_pe_file

            pe_file = str(_event_pe_file(event, hc.default_pe_cache()))
        s = sc.extract_pe(pe_file, pe_label, wd)
        s["event"] = event
        (wd / "event.json").write_text(json.dumps(s, indent=1))
        sc.fetch_strain(wd, n_offsource=n_background, sampling_frequency=sampling_frequency)
    heavy = [s for s in stages if s in ("fit", "background", "sample", "seeded")]
    if heavy and executor == "slurm":
        scripts = write_slurm_chains(wd, python or sys.executable, heavy, n_background, sample_models, tau_d, npool,
                                     slurm_options, slurm_env_setup)
        for s in scripts:
            _log(f"Slurm chain written: {s}")
            if submit:
                from .slurm import submit as sbatch

                _log(sbatch(s).strip())
        heavy = []
    if "fit" in heavy:
        an = sc.Analysis(wd)
        sc._dump(wd / "vacuum.json", sc.fit_vacuum(an, an.catalog_starts()))
        n_grid = len(sc.ALPHAS) * len(sc.DEPHASING_TARGETS)
        _pool_map(_profile_worker, [(str(wd), i) for i in range(n_grid)], jobs)
        sc.main(["polish", "--workdir", str(wd)])
    if "background" in heavy:
        starts = json.loads((wd / "offsource" / "starts.json").read_text())
        _pool_map(_background_worker, [(str(wd), k) for k in range(min(n_background, len(starts)))], jobs)
    if "sample" in heavy:
        for m in sample_models:
            sc.sample(wd, m, tau_d if m == "superradiance" else None, npool=max(1, jobs))
    if "seeded" in heavy:
        sc.main(["seeds", "--workdir", str(wd), "--result", str(_seed_result(wd, tau_d))])
        starts = json.loads((wd / "offsource" / "starts.json").read_text())
        _pool_map(_seeded_worker, [(str(wd), k) for k in range(-1, min(n_background, len(starts)))], jobs)
    result = {}
    if "report" in stages:
        result = write_report(wd, out_report_html, out_summary_tsv, tau_d=tau_d)
    return result


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def _sampling_results(wd: Path, tau_d: float) -> dict:
    """log evidences and posteriors of the sampling runs present in <workdir>/sampling."""
    out = {}
    for m in SAMPLE_MODELS:
        f = wd / "sampling" / f"{_sample_label(m, tau_d)}_result.hdf5"
        if f.exists():
            import bilby

            r = bilby.core.result.read_in_result(str(f))
            out[m] = dict(log_evidence=float(r.log_evidence), log_evidence_err=float(r.log_evidence_err),
                          posterior=r.posterior, log_bayes_factor=float(r.log_bayes_factor),
                          nlive=(r.sampler_kwargs or {}).get("nlive"))
    return out


def plot_cloud_posteriors(samp: dict, out_png: Path, m_src_ref: float) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.3))
    for m, color in (("agnostic", COLORS["vacuum"]), ("superradiance", COLORS["cloud"])):
        if m not in samp:
            continue
        post = samp[m]["posterior"]
        axes[0].hist(post["rho_phi"] / 1e8, bins=40, density=True, histtype="step", lw=1.8, color=color, label=m)
        axes[1].hist(sc.scalar_mass_ev(post["alpha_phi"], m_src_ref) * 1e12, bins=40, density=True,
                     histtype="step", lw=1.8, color=color, label=m)
    axes[0].set_xlabel("ρ̄ (10⁸ g/cm³)", color=INK2, fontsize=9)
    axes[1].set_xlabel(f"m_φ (10⁻¹² eV, for M = {m_src_ref:.1f} M☉)", color=INK2, fontsize=9)
    for ax in axes:
        _style(ax)
        ax.set_ylabel("posterior density", color=INK2, fontsize=9)
        ax.legend(frameon=False, fontsize=8, labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png


def write_report(wd: Path, out_report_html: Optional[str], out_summary_tsv: Optional[str],
                 tau_d: float = 1e6) -> dict:
    import pandas as pd

    from .report import write_simple_html_report

    s = json.loads((wd / "event.json").read_text())
    event = s.get("event", "the event")
    plots = wd / "plots"
    an = sc.Analysis(wd)
    cat = dict(an.catalog_starts()[0], rho_phi=0.0, alpha_phi=0.1)
    rows = [dict(fit="catalog maximum-likelihood sample (IMRPhenomXPHM)", log_likelihood_ratio=an.lnl(cat))]
    paras, images, result, bands = [], [], {}, None
    paras.append(
        f"<b>{event}</b> with the IMRPhenomXPHM_Scalar waveform of Roy et al. (arXiv:2510.17967): IMRPhenomXPHM "
        "plus the angular momentum the binary loses to a scalar cloud of mean density ρ̄ and α = M m<sub>φ</sub> "
        "(Eq. 2 and 6 of the paper), as a phase correction of the inspiral. Data: GWOSC strain, "
        f"{s['duration']:g} s at {s.get('sampling_frequency', 2048):g} Hz from {s['f_low']:g} Hz, noise spectra "
        f"and catalog samples of {s['label']}.")
    if "pe_max_log_likelihood_ratio" in s:
        paras.append(f"The PE file records a maximum ln L ratio of {s['pe_max_log_likelihood_ratio']:.1f}; this "
                     f"likelihood gives {rows[0]['log_likelihood_ratio']:.1f} at the catalog's best sample (the "
                     "phases agree; the difference comes from the noise-spectrum, window and calibration handling, "
                     "and cancels in the comparison of the two models, which use the same likelihood).")
    vac_f, sc_f = wd / "vacuum.json", wd / "scalar.json"
    seeded_f = wd / "seeded" / "event.json"
    if seeded_f.exists():
        # the statistic with the cloud fits started in the cloud mode of the event's nested sampling (it reaches the
        # mode that fits started from vacuum miss)
        prof = json.loads(sc_f.read_text())["profile"] if sc_f.exists() else []
        sc_f = seeded_f
        bkg_dir = wd / "seeded"
    else:
        prof, bkg_dir = None, wd / "background"
    if vac_f.exists():
        vac = json.loads(vac_f.read_text())
        rows.append(dict(fit="vacuum best fit (IMRPhenomXPHM)", log_likelihood_ratio=vac["log_likelihood_ratio"]))
    if sc_f.exists():
        res = json.loads(sc_f.read_text())
        p = res["params"]
        m_det, m_src, q = _masses(p)
        d20 = sc.dephasing_22(np.array([s["f_low"]]), m_det, m_src, q, p["rho_phi"], p["alpha_phi"])[0]
        rows.append(dict(fit="cloud best fit (IMRPhenomXPHM_Scalar)", log_likelihood_ratio=res["log_likelihood_ratio"],
                         rho_g_cm3=p["rho_phi"], alpha=p["alpha_phi"],
                         m_phi_eV=sc.scalar_mass_ev(p["alpha_phi"], m_src),
                         cloud_mass_fraction=sc.cloud_mass_fraction(p["rho_phi"], p["alpha_phi"], m_src),
                         dephasing_at_f_low_rad=abs(d20)))
        result["Lambda"] = res["Lambda"]
        paras.append(f"<b>Best fit.</b> Λ = 2 (ln L<sub>cloud</sub> − ln L<sub>vacuum</sub>) = <b>{res['Lambda']:.2f}</b>"
                     f" (ln L {res['vacuum_log_likelihood_ratio']:.2f} → {res['log_likelihood_ratio']:.2f}). "
                     "The cloud model contains vacuum (ρ̄ = 0), so Λ ≥ 0 even without a cloud: it is read against "
                     "the background below.")
        images.append(plot_whitened(an, [("vacuum best fit", vac["params"], COLORS["vacuum"]),
                                         ("catalog max-likelihood sample", cat, COLORS["catalog"]),
                                         ("cloud best fit", p, COLORS["cloud"])], plots / "whitened.png"))
        images.append(plot_lnl_by_frequency(an, p, vac["params"], plots / "lnl_by_frequency.png"))
        bands = lnl_gain_by_band(an, p, res.get("vacuum", vac)["params"])
        images.append(plot_dephasing(an, p, vac["params"], plots / "dephasing.png"))
        points = res["profile"] if prof is None else prof
        extra = [json.loads(f.read_text()) for f in sorted((wd / "profile_extra").glob("*.json"))] \
            if (wd / "profile_extra").is_dir() else []
        if points:
            images.append(plot_profile(points + extra, vac["log_likelihood_ratio"], plots / "profile.png"))
        if seeded_f.exists():
            paras.append(f"The cloud fits start from the {len(res['fits'])} best samples (at fixed merger time) of "
                         "the nested sampling of the event with the superradiance prior: fits started from the vacuum "
                         "best fit stay in its mode, a few ln L units below the cloud mode at a lower chirp mass. The "
                         "same starting points, moved to the time of each injection, are used on the background, which "
                         "therefore does not include the event's search of the parameter space: its p-value is "
                         "optimistic.")
    bkg = sorted(bkg_dir.glob("inj_*.json")) if bkg_dir.is_dir() else []
    if bkg and "Lambda" in result:
        lams = np.array([json.loads(f.read_text())["Lambda"] for f in bkg])
        pv, k = background_p_value(lams, result["Lambda"])
        result.update(background_n=len(lams), p_value=pv)
        images.append(plot_background(lams, result["Lambda"], plots / "background.png"))
        paras.append(f"<b>Background.</b> {len(lams)} injections of the vacuum best fit into off-source GWOSC noise "
                     f"within an hour of the event, analysed with the same procedure: Λ ≥ {result['Lambda']:.2f} in "
                     f"{k} of them, p = {pv:.3f}. Median background Λ = {np.median(lams):.2f}, 90th percentile "
                     f"{np.percentile(lams, 90):.2f}.")
    samp = _sampling_results(wd, tau_d)
    if samp:
        tab = []
        for m, r in samp.items():
            row = dict(model=m, log_evidence=r["log_evidence"], log_evidence_err=r["log_evidence_err"])
            if m != "vacuum" and "vacuum" in samp:
                row["ln_B_cloud_vs_vacuum"] = r["log_evidence"] - samp["vacuum"]["log_evidence"]
                result[f"ln_B_{m}"] = row["ln_B_cloud_vs_vacuum"]
            tab.append(row)
        rows_s = pd.DataFrame(tab)
        paras.append("<b>Nested sampling</b> (bilby, dynesty, acceptance-walk; the paper uses 1000 live points): "
            + "; ".join(
            f"{r['model']} ({samp[r['model']]['nlive']} live points): ln Z = {r['log_evidence']:.2f} ± {r['log_evidence_err']:.2f}"
            + (f", ln B = {r['ln_B_cloud_vs_vacuum']:.2f}" if "ln_B_cloud_vs_vacuum" in r and
               np.isfinite(r["ln_B_cloud_vs_vacuum"]) else "") for r in tab)
            + f". The paper finds ln B ≈ 0.4 (agnostic) and 3.5 (superradiance, τ_d = 10⁶ yr).")
        if vac_f.exists():
            m_ref = _masses(json.loads(vac_f.read_text())["params"])[1]
            images.append(plot_cloud_posteriors(samp, plots / "cloud_posteriors.png", m_src_ref=m_ref))
    table = pd.DataFrame(rows)
    if out_summary_tsv:
        table.to_csv(out_summary_tsv, sep="\t", index=False, float_format="%.6g")
    if out_report_html:
        paras.append("<b>Model limits.</b> Newtonian torque, inspiral only, valid for β = R / r<sub>φ</sub> ≪ 1; the "
                     "correction is applied per J-frame mode (exact for aligned spins); 1/(1 + δ) is kept "
                     "unlinearised, so that a torque larger than the GW loss stays finite. A large Λ means that the "
                     "vacuum model fits worse than usual in the inspiral: other effects (eccentricity, waveform "
                     "systematics, noise) can do the same.")
        write_simple_html_report(out_report_html, title=f"Scalar cloud: {event}", paragraphs=paras,
                                 images=[str(i) for i in images],
                                 tables=[("Fits", table.to_html(index=False, float_format=lambda x: f"{x:.4g}",
                                                               na_rep=""))]
                                 + ([("ln L cloud − ln L vacuum by detector and frequency band",
                                      bands.to_html(index=False, float_format=lambda x: f"{x:+.2f}"))]
                                    if bands is not None else [])
                                 + ([("Evidences", rows_s.to_html(index=False, float_format=lambda x: f"{x:.3f}",
                                                                   na_rep=""))] if samp else []))
        _log(f"report written to {out_report_html}")
    return result

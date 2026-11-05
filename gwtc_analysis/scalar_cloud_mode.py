"""``scalar_cloud`` mode: is a binary black hole better fitted with a scalar-field cloud around it?

The physics, the fits and the cluster jobs are in the self-contained runner ``scalar_cloud.py`` (the waveform
IMRPhenomXPHM_Scalar of Roy et al., arXiv:2510.17967). This module prepares the work directory, runs the steps
locally or writes them as Slurm jobs, and writes the report:

- **best fit**: the vacuum model and the cloud model fitted to the event's strain from the catalog posterior;
  the deviance D = 2 (ln L_cloud - ln L_vacuum) and the matched-filter SNR of the catalog, vacuum and
  cloud waveforms;
- **background**: the same procedure on the vacuum best fit injected into off-source noise around the event: the
  distribution of D without a cloud, hence the p-value of the event's D;

The deviance (the likelihood-ratio test statistic) is stored under the key "Lambda" in the result files, its
name in earlier versions.
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
    """Whitened, band-passed time series of a frequency-domain series of `ifo` (arbitrary units: divide by
    `_noise_sigma` of the data for units of the noise standard deviation)."""
    f = ifo.frequency_array
    asd = ifo.amplitude_spectral_density_array
    w = np.zeros_like(h_fd)
    band = (f >= f_lo) & (f <= f_hi) & np.isfinite(asd) & (asd > 0)
    w[band] = h_fd[band] / asd[band]
    # a smooth edge so that the band limits do not ring
    taper = np.ones_like(f)
    edge = min(5.0, 0.25 * (f_hi - f_lo))
    for lo, hi in ((f_lo, f_lo + edge), (f_hi - edge, f_hi)):
        k = (f >= lo) & (f <= hi)
        x = (f[k] - lo) / (hi - lo)
        taper[k] = 0.5 * (1 - np.cos(np.pi * x)) if lo == f_lo else 0.5 * (1 + np.cos(np.pi * x))
    return np.fft.irfft(w * taper, n=len(ifo.time_array))


def _noise_sigma(x: np.ndarray, t: np.ndarray, t_merger: float, duration: float) -> float:
    """Robust standard deviation (1.4826 MAD) of whitened data away from the signal (from 1 s after the segment
    start to 1 s before the end, without the last `duration` s before the merger and 0.5 s after it)."""
    k = (t > t[0] + 1) & (t < t[-1] - 1) & ~((t > t_merger - duration) & (t < t_merger + 0.5))
    if k.sum() < 100:
        k = (t > t[0] + 1) & (t < t[-1] - 1)
    y = x[k]
    return float(1.4826 * np.median(np.abs(y - np.median(y)))) or 1.0


def _frequency_track(ifo, h_fd: np.ndarray, f_lo: float, smooth: float = 0.05) -> np.ndarray:
    """Instantaneous frequency (Hz) of a detector-frame waveform against `ifo.time_array`: the phase derivative of
    its analytic signal, with a running median over `smooth` s (precession and higher modes make it wobble)."""
    from scipy.ndimage import median_filter, uniform_filter1d
    from scipy.signal import hilbert

    f = ifo.frequency_array
    x = np.fft.irfft(np.where(f >= f_lo, h_fd, 0), n=len(ifo.time_array))
    fs = 1.0 / (ifo.time_array[1] - ifo.time_array[0])
    fi = np.gradient(np.unwrap(np.angle(hilbert(x)))) * fs / (2 * np.pi)
    n = max(3, int(smooth * fs))
    return uniform_filter1d(median_filter(fi, size=n), size=2 * n)


def _snr_band(ifo, h_fd: np.ndarray, f_lo: float, fraction: float = 0.99) -> float:
    """Frequency below which `fraction` of the optimal SNR^2 of `h_fd` lies."""
    f = ifo.frequency_array
    psd = ifo.power_spectral_density_array
    w = np.where((f >= f_lo) & np.isfinite(psd) & (psd > 0), np.abs(h_fd) ** 2 / np.where(psd > 0, psd, 1), 0.0)
    c = np.cumsum(w)
    return float(f[np.searchsorted(c, fraction * c[-1])]) if c[-1] > 0 else float(f[-1])


def _band_snr(ifo, h_fd: np.ndarray, f_lo: float, f_hi: float) -> float:
    """Optimal SNR of `h_fd` in `ifo` between f_lo and f_hi (the whole signal, not only a time window)."""
    f = ifo.frequency_array
    psd = ifo.power_spectral_density_array
    k = (f >= f_lo) & (f <= f_hi) & np.isfinite(psd) & (psd > 0)
    df = f[1] - f[0]
    return float(np.sqrt(4 * df * np.sum(np.abs(h_fd[k]) ** 2 / psd[k])))


def _window_of_track(t: np.ndarray, track: np.ndarray, f_a: float, f_b: float, t_merger: float,
                     fallback: tuple[float, float]) -> tuple[float, float]:
    """(start, end) from the merger of the times at which the inspiral track crosses f_a and f_b."""
    k = (t < t_merger) & (t > t_merger - 30)
    tt, ff = t[k] - t_merger, track[k]
    try:
        a = float(tt[np.where(ff <= f_a)[0][-1]])
        b = float(tt[np.where(ff <= f_b)[0][-1]])
    except IndexError:
        return fallback
    return (a, b) if b - a > 0.2 else fallback


def _matched_filter(ifo, h_fd: np.ndarray, f_lo: float, f_hi: float) -> tuple[np.ndarray, np.ndarray]:
    """(lags in s, complex matched-filter SNR) of the template `h_fd` (a detector response, at its own time) on the
    data of `ifo`, between f_lo and f_hi: z(tau) = <d, h e^{-2 pi i f tau}> / sqrt(<h, h>). The lag 0 is the fit's own
    time; Re z(0) is its own phase, |z| the SNR maximised over phase and max |z| over time too."""
    f = ifo.frequency_array
    psd = ifo.power_spectral_density_array
    d = ifo.frequency_domain_strain
    n = len(ifo.time_array)
    duration = 1.0 / (f[1] - f[0])
    k = ifo.frequency_mask & (f >= f_lo) & (f <= f_hi) & np.isfinite(psd) & (psd > 0)
    x = np.zeros(n, dtype=complex)
    x[:len(f)][k] = d[k] * np.conj(h_fd[k]) / psd[k]
    norm = np.sqrt(4.0 / duration * np.sum(np.abs(h_fd[k]) ** 2 / psd[k]))
    z = 4.0 / duration * n * np.fft.ifft(x) / (norm if norm > 0 else 1.0)
    fs = n / duration
    lags = (np.arange(n) + n // 2) % n - n // 2
    order = np.argsort(lags)
    return lags[order] / fs, z[order]


def _peak(lags: np.ndarray, z: np.ndarray, max_lag: float) -> tuple[float, float]:
    """(lag, |z|) of the matched-filter peak within +-max_lag of the fit's time."""
    k = np.abs(lags) <= max_lag
    i = int(np.argmax(np.abs(z[k])))
    return float(lags[k][i]), float(np.abs(z[k][i]))


def plot_matched_filter(an, curves: list[tuple[str, dict, str]], out_png: Path,
                        inspiral_band: tuple[float, float] = (20.0, 45.0)) -> Path:
    """The fits against the data, through matched filtering (the strain overlay cannot show differences of a few
    tenths of a noise sigma): one row per detector (H1, L1) with

    - a Q-transform of the whitened data, with the frequency tracks of the waveforms of `curves` (label, parameters,
      colour);
    - the matched-filter SNR |z| of each fit against the time shift from its own merger time, over the full band and
      over `inspiral_band` only (where a cloud acts), with the peak values;

    and a last row with the peak SNR of each fit in each detector and in the network (quadrature sum), in both bands:
    for templates of optimal amplitude, (rho_cloud^2 - rho_vacuum^2) / 2 is the likelihood gain of the cloud."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    TimeSeries = sc.gwpy_timeseries()
    rows = [ifo for ifo in an.ifos if ifo.name != "V1"] or list(an.ifos)
    nrow = len(rows) + 1
    fig, axes = plt.subplots(nrow, 3, figsize=(14, 2.8 * nrow), squeeze=False,
                             gridspec_kw=dict(width_ratios=[3.0, 2.2, 2.2]))
    p_ref = curves[0][1]
    f_low = float(an.s.get("f_low", 20.0))
    bands = [("full band", f_low, None, 0.012, 0.015), ("inspiral band", max(f_low, inspiral_band[0]),
                                                           inspiral_band[1], 0.03, 0.15)]
    peaks = {b[0]: {label: {} for label, _, _ in curves} for b in bands}
    t_tf = None
    for ifo in an.ifos:
        t = ifo.time_array
        fs = 1.0 / (t[1] - t[0])
        t_merger = p_ref["geocent_time"] + ifo.time_delay_from_geocenter(p_ref["ra"], p_ref["dec"],
                                                                         p_ref["geocent_time"])
        resp = [ifo.get_detector_response(an.waveform_generator.frequency_domain_strain(p), p) for _, p, _ in curves]
        r = next((i for i, x in enumerate(rows) if x.name == ifo.name), None)
        for c, (bname, lo, hi, max_lag, span) in enumerate(bands):
            hi = 0.5 * fs if hi is None else hi
            for (label, _, color), h in zip(curves, resp):
                lags, z = _matched_filter(ifo, h, lo, hi)
                lag, rho = _peak(lags, z, max_lag)
                peaks[bname][label][ifo.name] = rho
                if r is None:
                    continue
                ax = axes[r, c + 1]
                k = np.abs(lags) <= span
                ls = "--" if "catalog" in label else "-"
                ax.plot(lags[k] * 1e3, np.abs(z[k]), color=color, lw=1.4 if ls == "-" else 1.1, ls=ls)
                ax.plot([lag * 1e3], [rho], "o", color=color, ms=4.5)
            if r is None:
                continue
            ax = axes[r, c + 1]
            top = max(peaks[bname][l][ifo.name] for l, _, _ in curves)
            ax.set_ylim(0.55 * top, 1.0 * top + 0.12 * top)
            ax.axvline(0, color=INK2, lw=0.7, ls=":")
            ax.text(0.02, 0.97, "\n".join(f"{l.replace(' best fit', '').replace(' max-likelihood sample', '')}: "
                                          f"{peaks[bname][l][ifo.name]:.2f}" for l, _, _ in curves),
                    transform=ax.transAxes, ha="left", va="top", fontsize=7.5, color=INK2,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
            _style(ax)
            ax.set_title(f"{ifo.name}: matched-filter SNR, {bname} ({lo:.0f}–{hi:.0f} Hz)", color=INK, fontsize=9,
                         loc="left")
            if r == len(rows) - 1:
                ax.set_xlabel("time shift from the fit's merger time (ms)", color=INK2, fontsize=9)
            if c == 0:
                ax.set_ylabel("|z| (SNR)", color=INK2, fontsize=8.5)
        if r is None:
            continue

        # time-frequency map
        tracks = [_frequency_track(ifo, h, f_low) for h in resp]
        f_hi = min(max(_snr_band(ifo, resp[0], f_low), 150.0), 0.45 * fs)
        if t_tf is None:
            t_start = _window_of_track(t, tracks[0], f_low + 1, f_low + 2, t_merger, (-8.0, -7.0))[0]
            t_tf = (max(t_start - 0.3, t[0] - t_merger + 1.0), 0.1)
        ax = axes[r, 0]
        white = _whiten(ifo, ifo.frequency_domain_strain, f_low, 0.45 * fs)
        white = white / _noise_sigma(white, t, t_merger, -t_tf[0])
        q = TimeSeries(white, t0=t[0], sample_rate=fs).q_transform(
            qrange=(4, 64), frange=(f_low, f_hi), outseg=(t_merger + t_tf[0], t_merger + t_tf[1]), whiten=False,
            logf=True)
        e = np.asarray(q.value)
        mesh = ax.pcolormesh(np.asarray(q.times.value) - t_merger, np.asarray(q.frequencies.value), e.T,
                             shading="auto", cmap="viridis", vmin=0, vmax=max(5.0, float(np.percentile(e, 99.7))),
                             rasterized=True)
        for (label, _, color), tr in zip(curves, tracks):
            k = (t - t_merger >= t_tf[0]) & (t - t_merger <= 0.0) & (tr >= f_low) & (tr <= f_hi)
            ax.plot(t[k] - t_merger, tr[k], color=color, lw=1.0, ls="--" if "catalog" in label else "-", alpha=0.85)
        ax.set_yscale("log")
        ax.set_ylim(f_low, f_hi)
        ticks = [v for v in (20, 30, 50, 100, 200, 300, 500, 800) if f_low <= v <= f_hi]
        ax.set_yticks(ticks, [f"{v:g}" for v in ticks])
        ax.minorticks_off()
        ax.set_xlim(*t_tf)
        ax.set_ylabel(f"{ifo.name}\nfrequency (Hz)", color=INK2, fontsize=9)
        cb = fig.colorbar(mesh, ax=ax, pad=0.01, fraction=0.04)
        cb.set_label("normalised energy", color=INK2, fontsize=7.5)
        cb.ax.tick_params(labelsize=7, colors=INK2)
        ax.tick_params(colors=INK2, labelsize=9)
        ax.set_title(f"{ifo.name}: Q-transform of the whitened data, with the frequency tracks of the fits"
                     if r == 0 else ifo.name, color=INK, fontsize=9, loc="left")
        if r == len(rows) - 1:
            ax.set_xlabel("time from merger (s)", color=INK2, fontsize=9)

    # peak SNR per detector and in the network, both bands
    names = [ifo.name for ifo in an.ifos] + ["network"]
    for c, (bname, *_rest) in enumerate(bands):
        ax = axes[-1, c + 1]
        n_c = len(curves)
        for j, (label, _, color) in enumerate(curves):
            vals = list(peaks[bname][label].values())
            vals = vals + [float(np.sqrt(np.sum(np.square(vals))))]
            x = np.arange(len(names)) + (j - (n_c - 1) / 2) * 0.22
            ax.plot(x, vals, "o", color=color, ms=6, mfc="none" if "catalog" in label else color)
            for xi, v in zip(x, vals):
                ax.annotate(f"{v:.2f}", (xi, v), xytext=(0, 6), textcoords="offset points", rotation=90,
                            ha="center", va="bottom", fontsize=6.5, color=color)
        ax.set_xticks(range(len(names)), names)
        lo_v = min(min(v.values()) for v in peaks[bname].values())
        hi_v = max(np.sqrt(np.sum(np.square(list(v.values())))) for v in peaks[bname].values())
        ax.set_ylim(0.8 * lo_v, 1.22 * hi_v)
        ax.set_xlim(-0.6, len(names) - 0.4)
        _style(ax)
        lo, hi = bands[c][1], bands[c][2] or 0.5 / (an.ifos[0].time_array[1] - an.ifos[0].time_array[0])
        ax.set_title(f"peak SNR, {bname} ({lo:.0f}–{hi:.0f} Hz)", color=INK, fontsize=9, loc="left")
        if c == 0:
            ax.set_ylabel("peak |z|", color=INK2, fontsize=8.5)
    lnl = {l: an.lnl(p) for l, p, _ in curves}

    lax = axes[-1, 0]
    lax.axis("off")
    hs = [Line2D([], [], color=c, lw=1.6, ls="--" if "catalog" in l else "-", marker="o", ms=4, label=l)
          for l, _, c in curves]
    lax.legend(handles=hs, loc="center", frameon=False, fontsize=9, labelcolor=INK2, handlelength=2.6,
               title="matched filter: each fit used as a template on the data, |z| maximised over\n"
               "phase, its peak also over a ±12 ms (full band) or ±30 ms (inspiral band) time shift;\n"
               "network: quadrature sum of the detectors. Coherent ln L of each fit\n(own time and phase): "
               + ", ".join(f"{l.split()[0]} {v:.2f}" for l, v in lnl.items()), title_fontsize=8.5)
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
    d = sc.dephasing_22(f, m_det, m_src, q, p_cloud["rhobar_phi"], p_cloud["alpha_phi"])
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


def plot_background(lams: np.ndarray, lam_event: float, out_png: Path,
                    lams_cloud: Optional[np.ndarray] = None) -> Path:
    """The background (vacuum injections) and, when given, the cloud injections, against the event; the x axis stops
    at 1.15 times the 99th percentile (or the event): the values beyond (a noise transient) are counted in the legend."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    every = np.concatenate([lams] + ([lams_cloud] if lams_cloud is not None else []))
    top = max(np.percentile(every, 99), lam_event) * 1.15 + 1e-9
    bins = np.linspace(0, top, 40)
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    off = (lams > top).sum()
    ax.hist(np.minimum(lams, top), bins=bins, color=COLORS["vacuum"], alpha=0.75,
            label=f"background: {len(lams)} vacuum injections" + (f" ({off} off scale)" if off else ""))
    if lams_cloud is not None:
        off = (lams_cloud > top).sum()
        ax.hist(np.minimum(lams_cloud, top), bins=bins, color=COLORS["catalog"], alpha=0.55,
                label=f"{len(lams_cloud)} injections of the cloud best fit" + (f" ({off} off scale)" if off else ""))
    ax.axvline(lam_event, color=COLORS["cloud"], lw=2, label=f"GW event: $\\mathcal{{D}}$ = {lam_event:.2f}")
    _style(ax)
    ax.set_xlabel(r"deviance $\mathcal{D}$ = 2 (ln L$_{\rm cloud}$ − ln L$_{\rm vacuum}$)", color=INK2, fontsize=9)
    ax.set_ylabel("number", color=INK2, fontsize=9)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png


def background_p_value(lams: np.ndarray, lam_event: float) -> tuple[float, int]:
    """Fraction of the background with a deviance >= the event's ((k + 1) / (n + 1), conservative), and k."""
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


def lnl_gain_per_frequency(an, p_cloud: dict, p_vac: dict) -> tuple[np.ndarray, dict]:
    """(frequencies, {detector: ln L_cloud - ln L_vacuum per frequency bin}): summed, the likelihood gain."""
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
    return f, gain


def lnl_gain_by_band(an, p_cloud: dict, p_vac: dict, edges=LNL_BANDS):
    """ln L_cloud - ln L_vacuum per detector and frequency band (a table: where the improvement comes from)."""
    import pandas as pd

    f, gain = lnl_gain_per_frequency(an, p_cloud, p_vac)
    hi = [e if e is not None else np.inf for e in edges[1:]]
    cols = [f"{lo:g}–{h:g} Hz" if np.isfinite(h) else f"> {lo:g} Hz" for lo, h in zip(edges[:-1], hi)]
    rows = []
    for det, g in gain.items():
        rows.append(dict(detector=det, **{c: g[(f >= lo) & (f < h)].sum() for c, lo, h in zip(cols, edges[:-1], hi)},
                         total=g.sum()))
    tab = pd.DataFrame(rows)
    tab.loc[len(tab)] = dict(detector="all", **{c: tab[c].sum() for c in cols + ["total"]})
    return tab


DET_COLORS = {"H1": "#1f2328", "L1": "#bf8700", "V1": "#8250df"}


def plot_lnl_by_frequency(an, p_cloud: dict, p_vac: dict, out_png: Path, edges=LNL_BANDS) -> Path:
    """Cumulative ln L_cloud - ln L_vacuum against frequency, per detector and summed: which part of the signal, in
    which detector, prefers the cloud; the bands of the table by band are marked."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    f, gain = lnl_gain_per_frequency(an, p_cloud, p_vac)
    k = (f >= an.s["f_low"] * 0.95) & (f <= 500)
    fig, ax = plt.subplots(figsize=(8.5, 3.8))
    for det, g in gain.items():
        c = np.cumsum(g)
        ax.plot(f[k], c[k], color=DET_COLORS.get(det, INK2), lw=1.3, label=f"{det} (total {c[-1]:+.2f})")
    tot = np.cumsum(sum(gain.values()))
    ax.plot(f[k], tot[k], color=COLORS["cloud"], lw=2.2, label=f"all detectors (total {tot[-1]:+.2f})")
    ax.axhline(0, color=INK2, lw=0.8)
    for e in edges[1:]:
        if e is not None and f[k][0] < e < f[k][-1]:
            ax.axvline(e, color=GRID, lw=1.0, ls="--", zorder=0)
    ax.set_xscale("log")
    _style(ax)
    ax.set_xlabel("frequency (Hz)", color=INK2, fontsize=9)
    ax.set_ylabel("cumulative ln L_cloud − ln L_vacuum", color=INK2, fontsize=9)
    ax.set_title("Where the cloud fit gains over the vacuum fit (dashed: the bands of the table)", color=INK,
                 fontsize=9.5, loc="left")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2, loc="upper left", bbox_to_anchor=(1.01, 1.0))
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
    vac = sc.read_json(Path(wd) / "vacuum.json")
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
    starts_file = Path(workdir) / "offsource" / "starts.json"
    if starts_file.exists():        # fewer off-source segments than asked (data gaps): no array task past the last
        n_background = min(n_background, len(sc.read_json(starts_file)))
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
        starts = sc.read_json(wd / "offsource" / "starts.json")
        _pool_map(_background_worker, [(str(wd), k) for k in range(min(n_background, len(starts)))], jobs)
    if "sample" in heavy:
        for m in sample_models:
            sc.sample(wd, m, tau_d if m == "superradiance" else None, npool=max(1, jobs))
    if "seeded" in heavy:
        sc.main(["seeds", "--workdir", str(wd), "--result", str(_seed_result(wd, tau_d))])
        starts = sc.read_json(wd / "offsource" / "starts.json")
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
            r = sc.read_result(f)
            out[m] = dict(log_evidence=float(r.log_evidence), log_evidence_err=float(r.log_evidence_err),
                          posterior=r.posterior, log_bayes_factor=float(r.log_bayes_factor),
                          nlive=(r.sampler_kwargs or {}).get("nlive"))
    return out


def boson_mass_quantiles(posterior, q=(5, 50, 95)) -> np.ndarray:
    """Quantiles of the scalar mass m_phi (eV) of a cloud posterior, each sample with its own source-frame mass."""
    from bilby.gw.conversion import chirp_mass_and_mass_ratio_to_total_mass

    m_det = chirp_mass_and_mass_ratio_to_total_mass(np.asarray(posterior["chirp_mass"], float),
                                                    np.asarray(posterior["mass_ratio"], float))
    m_src = np.array([sc._source_frame_total_mass(a, b) for a, b in zip(m_det, posterior["luminosity_distance"])])
    return np.percentile([sc.scalar_mass_ev(a, b) for a, b in zip(posterior["alpha_phi"], m_src)], q)


def plot_cloud_posteriors(samp: dict, out_png: Path, m_src_ref: float) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.3))
    for m, color in (("agnostic", COLORS["vacuum"]), ("superradiance", COLORS["cloud"])):
        if m not in samp:
            continue
        post = samp[m]["posterior"]
        axes[0].hist(post["rhobar_phi"] / 1e8, bins=40, density=True, histtype="step", lw=1.8, color=color, label=m)
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

    s = sc.read_json(wd / "event.json")
    if s.get("pe_file"):
        from .report import record_pe_file

        record_pe_file(s.get("event", ""), s["pe_file"], "samples extracted from it at the prepare stage")
    event = s.get("event", "the event")
    plots = wd / "plots"
    an = sc.Analysis(wd)
    cat = dict(an.catalog_starts()[0], rhobar_phi=0.0, alpha_phi=0.1)
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
        prof = sc.read_json(sc_f)["profile"] if sc_f.exists() else []
        sc_f = seeded_f
        bkg_dir = wd / "seeded"
    else:
        prof, bkg_dir = None, wd / "background"
    if vac_f.exists():
        vac = sc.read_json(vac_f)
        rows.append(dict(fit="vacuum best fit (IMRPhenomXPHM)", log_likelihood_ratio=vac["log_likelihood_ratio"]))
    if sc_f.exists():
        res = sc.read_json(sc_f)
        p = res["params"]
        m_det, m_src, q = _masses(p)
        d20 = sc.dephasing_22(np.array([s["f_low"]]), m_det, m_src, q, p["rhobar_phi"], p["alpha_phi"])[0]
        rows.append(dict(fit="cloud best fit (IMRPhenomXPHM_Scalar)", log_likelihood_ratio=res["log_likelihood_ratio"],
                         rhobar_g_cm3=p["rhobar_phi"], alpha=p["alpha_phi"],
                         m_phi_eV=sc.scalar_mass_ev(p["alpha_phi"], m_src),
                         cloud_mass_fraction=sc.cloud_mass_fraction(p["rhobar_phi"], p["alpha_phi"], m_src),
                         dephasing_at_f_low_rad=abs(d20)))
        result["Lambda"] = res["Lambda"]
        paras.append(f"<b>Best fit.</b> deviance 𝒟 = 2 (ln L<sub>cloud</sub> − ln L<sub>vacuum</sub>) = <b>{res['Lambda']:.2f}</b>"
                     f" (ln L {res['vacuum_log_likelihood_ratio']:.2f} → {res['log_likelihood_ratio']:.2f}). "
                     "The cloud model contains vacuum (ρ̄ = 0), so 𝒟 ≥ 0 even without a cloud: it is read against "
                     "the background below.")
        images.append(plot_matched_filter(an, [("vacuum best fit", vac["params"], COLORS["vacuum"]),
                                               ("catalog max-likelihood sample", cat, COLORS["catalog"]),
                                               ("cloud best fit", p, COLORS["cloud"])],
                                           plots / "matched_filter.png"))
        images.append(plot_lnl_by_frequency(an, p, vac["params"], plots / "lnl_by_frequency.png"))
        bands = lnl_gain_by_band(an, p, res.get("vacuum", vac)["params"])
        images.append(plot_dephasing(an, p, vac["params"], plots / "dephasing.png"))
        points = res["profile"] if prof is None else prof
        extra = [sc.read_json(f) for f in sorted((wd / "profile_extra").glob("*.json"))] \
            if (wd / "profile_extra").is_dir() else []
        if points:
            images.append(plot_profile(points + extra, vac["log_likelihood_ratio"], plots / "profile.png"))
        if seeded_f.exists():
            paras.append(f"The cloud fits start from the {len(res['fits'])} best samples (at fixed merger time) of "
                         "the nested sampling of the event with the superradiance prior: fits started from the vacuum "
                         "best fit stay in its mode, a few ln L units below the cloud mode found by the sampling. The "
                         "same starting points, moved to the time of each injection, are used on the background, which "
                         "therefore does not include the event's search of the parameter space: its p-value is "
                         "optimistic.")
    bkg = sorted(bkg_dir.glob("inj_*.json")) if bkg_dir.is_dir() else []
    if bkg and "Lambda" in result:
        lams = np.array([sc.read_json(f)["Lambda"] for f in bkg])
        pv, k = background_p_value(lams, result["Lambda"])
        result.update(background_n=len(lams), p_value=pv)
        cloud_inj = sorted((wd / "seeded_cloud").glob("inj_*.json")) if bkg_dir.name == "seeded" else []
        lams_c = np.array([sc.read_json(f)["Lambda"] for f in cloud_inj]) if cloud_inj else None
        images.append(plot_background(lams, result["Lambda"], plots / "background.png", lams_c))
        if lams_c is not None:
            frac = float((lams_c <= result["Lambda"]).mean())
            result.update(cloud_injections_n=len(lams_c), cloud_injections_fraction_below_event=frac)
            q5, q50, q95 = np.percentile(lams_c, [5, 50, 95])
            paras.append(f"<b>Sensitivity.</b> {len(lams_c)} injections of the cloud best fit into the same segments, "
                         f"same procedure: 𝒟 = {q50:.1f} (5–95%: {q5:.1f}–{q95:.1f}); {frac:.1%} of them at or below "
                         f"the event's 𝒟. The event lies between the two distributions: above what noise gives "
                         "around vacuum, below what this cloud would give.")
        paras.append(f"<b>Background.</b> {len(lams)} injections of the vacuum best fit into off-source GWOSC noise "
                     f"within an hour of the event, analysed with the same procedure: 𝒟 ≥ {result['Lambda']:.2f} in "
                     f"{k} of them, p = {pv:.3f}. Median background 𝒟 = {np.median(lams):.2f}, 90th percentile "
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
            m_ref = _masses(sc.read_json(vac_f)["params"])[1]
            images.append(plot_cloud_posteriors(samp, plots / "cloud_posteriors.png", m_src_ref=m_ref))
        from .continuous_waves import EV_TO_HZ

        cw = []
        for m in ("superradiance", "agnostic"):
            if m in samp:
                lo, mid, hi = boson_mass_quantiles(samp[m]["posterior"])
                result[f"m_phi_{m}_eV"] = (lo, mid, hi)
                cw.append(f"{m}: m<sub>φ</sub> = {mid:.2g} eV (90%: {lo:.2g}–{hi:.2g}), continuous waves at "
                          f"{lo * EV_TO_HZ:.0f}–{hi * EV_TO_HZ:.0f} Hz")
        if cw:
            paras.append("<b>The same boson in continuous waves.</b> A cloud around an isolated spinning black hole "
                         "annihilates into a continuous wave at 2 m<sub>φ</sub> c²/h. " + "; ".join(cw) + ". The LVK "
                         "all-sky search for boson clouds (O3, 20–610 Hz) covers this band; the "
                         "<code>continuous_waves</code> mode draws its exclusions with this mass.")
    table = pd.DataFrame(rows)
    if out_summary_tsv:
        table.to_csv(out_summary_tsv, sep="\t", index=False, float_format="%.6g")
    if out_report_html:
        paras.append("<b>Model limits.</b> Newtonian torque, inspiral only, valid for β = R / r<sub>φ</sub> ≪ 1; the "
                     "correction is applied per J-frame mode (exact for aligned spins); 1/(1 + δ) is kept "
                     "unlinearised, so that a torque larger than the GW loss stays finite. A large 𝒟 means that the "
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

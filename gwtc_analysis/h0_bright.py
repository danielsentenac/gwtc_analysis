"""Hubble constant from bright sirens (``hubble_constant --method bright``): the GW luminosity distance and the
redshift of an identified host.

For an event whose host is identified, the GW signal gives the luminosity distance d_L and the host gives
the Hubble-flow redshift z (for a nearby galaxy, from its recession velocity minus its peculiar velocity).
In flat ΛCDM, d_L = (c / H0) D(z; Ωm), with Ωm = 0.3065 as in the spectral siren.

Following Chen, Fishbach & Holz 2018 (arXiv:1712.06531) and Mandel, Farr & Gair 2019 (arXiv:1809.02063),
for sources uniform in comoving volume and source-frame time, p_pop(z) ∝ dV_c/dz / (1 + z),

    p(H0 | data) ∝ p(H0) / β(H0) ∫ N(z_obs; z, σ_z) p_pop(z) L_GW(d_L(z, H0)) dz,

with L_GW = (posterior of d_L) / (PE distance prior). Over the PE samples d_i, with z_i the redshift of
d_i at this H0, the integral is a kernel estimate

    Σ_i N(z_obs; z_i, σ) p_pop(z_i) / [π_PE(d_i) dd_L/dz(z_i)],   σ = max(σ_z, h),

where h is the kernel width of the z_i (Silverman's rule): it only matters when σ_z is narrower than the
sample spacing (a distant host with a spectroscopic redshift).

β(H0) is the fraction of the population that is detected:

- ``euclidean``: GW-limited detection of nearby sources, β ∝ H0^3; it then cancels the H0^3 of the
  numerator and, for a d_L^2 PE prior, p(H0) ∝ Σ_i N(v_H; H0 d_i, σ_v) (LVK 2017, arXiv:1710.05835).
- ``injections``: the LVK sensitivity injections of the event's observing run, reweighted to the
  population (the mass model of its class; the injected isotropic spins) at each trial H0.
- ``auto`` (default): ``euclidean`` below z = 0.05, ``injections`` above; ``injections`` with a population model.

The population (``--population``) is, by default, sources uniform in comoving volume and source-frame time, with no
mass term in the numerator (the mass of the event's class enters the selection only). ``fullpop4`` is the model of the
GWTC-4.0 reanalysis of GW170817 (arXiv:2509.04348, Appendix E): the FullPop-4.0 mass distribution (BNS, NSBH and BBH
in one, icarogw's ``m1m2_paired_massratio_bplmulti_dip``) and the Madau-Dickinson rate, fixed to the medians of the
paper's spectral siren. The numerator then carries p_m(m1/(1+z), m2/(1+z)) ψ(z) / (1+z)^2 for each sample (the PE
mass prior is uniform in the detector-frame masses), and the selection the same density.

The samples must be conditioned on the sky position of the counterpart: the GWTC-1 GW170817 samples are
(the sky was fixed to AT2017gfo); for samples that are not, only those near the counterpart are kept.
An independent constraint on the viewing angle (e.g. from the jet) multiplies the GW likelihood: the samples are
weighted by it. The counterparts and the PE samples come from the `counterpart` module.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Optional

import numpy as np
import pandas as pd

from .counterpart import (COUNTERPARTS, Counterpart, default_style, get_counterpart, hubble_flow_redshift,
                          load_event_samples, sky_conditioned, viewing_angle, viewing_angle_weights)
from .report import write_simple_html_report

C_KMS = 299792.458
OM0 = 0.3065                      # as the spectral siren (hubble_constant)
H0_PRIOR = (10.0, 200.0)          # the prior of the spectral siren, so that they combine directly
EUCLIDEAN_ZMAX = 0.05             # selection "auto": euclidean below, injections above
MIN_NEFF_INJECTIONS = 50
_trapz = getattr(np, "trapezoid", None) or np.trapz

# ---------------------------------------------------------------------------
# cosmology: dimensionless flat ΛCDM, distances in units of c / H0
# ---------------------------------------------------------------------------
_ZG = np.concatenate([np.linspace(0.0, 0.1, 20001)[:-1], np.geomspace(0.1, 10.0, 40001)])
_EZ = np.sqrt(OM0 * (1 + _ZG) ** 3 + 1 - OM0)
_DC = np.concatenate([[0.0], np.cumsum(0.5 * (1 / _EZ[1:] + 1 / _EZ[:-1]) * np.diff(_ZG))])
_DL = (1 + _ZG) * _DC
_DDL = _DC + (1 + _ZG) / _EZ
_POPZ = _DC ** 2 / _EZ / (1 + _ZG)          # dV_c/dz / (1 + z), up to (c / H0)^3 (normalized over z: H0-free)


def _at_dimensionless_distance(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(z, p_pop(z), dD_L/dz) at dimensionless luminosity distances x = d_L H0 / c, with one search."""
    i = np.clip(np.searchsorted(_DL, x) - 1, 0, len(_DL) - 2)
    f = np.clip((x - _DL[i]) / (_DL[i + 1] - _DL[i]), 0.0, 1.0)
    lin = lambda tab: tab[i] + f * (tab[i + 1] - tab[i])
    return lin(_ZG), lin(_POPZ), lin(_DDL)


def luminosity_distance(z, h0):
    """d_L (Mpc) of redshift z for H0 (km/s/Mpc)."""
    return C_KMS / np.asarray(h0, float) * np.interp(z, _ZG, _DL)


def redshift_at(dl, h0):
    """Redshift of luminosity distance dl (Mpc) for H0."""
    return np.interp(np.asarray(dl, float) * np.asarray(h0, float) / C_KMS, _DL, _ZG)


def ddl_dz(z, h0):
    return C_KMS / np.asarray(h0, float) * np.interp(z, _ZG, _DDL)


def pop_z(z):
    """Redshift density of sources uniform in comoving volume and source-frame time (unnormalized)."""
    return np.interp(z, _ZG, _POPZ)


# ---------------------------------------------------------------------------
# population models (--population): FullPop-4.0 and the Madau-Dickinson rate, as in icarogw
# ---------------------------------------------------------------------------
# medians of the GWTC-4.0 spectral siren, used by the paper for GW170817 (arXiv:2509.04348, Appendix E); names of
# icarogw's m1m2_paired_massratio_bplmulti_dip (bottomsmooth, topsmooth: δm min and max; leftdip...deep: the dip)
FULLPOP4_GWTC4 = dict(alpha_1=2.52, alpha_2=1.79, beta_bottom=1.30, beta_top=2.62, mmin=0.99, mmax=93.8,
                      bottomsmooth=0.074, topsmooth=0.035, mu_g_low=8.98, sigma_g_low=0.77, mu_g_high=26.7,
                      sigma_g_high=7.8, lambda_g=0.18, lambda_g_low=0.72, leftdip=2.33, rightdip=7.23,
                      leftdipsmooth=0.11, rightdipsmooth=0.14, deep=0.52)
MADAU_GWTC4 = dict(gamma=3.57, kappa=2.97, zp=2.59)


def _highpass(m: np.ndarray, mmin: float, delta: float) -> np.ndarray:
    """icarogw's _highpass_filter: 0 below mmin, a smooth rise over delta, 1 above."""
    out = np.where(m >= mmin + delta, 1.0, 0.0)
    w = (m > mmin) & (m < mmin + delta)
    x = m[w] - mmin
    with np.errstate(over="ignore"):
        out[w] = 1.0 / (np.exp(delta / x + delta / (x - delta)) + 1.0)
    return out


def _lowpass(m: np.ndarray, mmax: float, delta: float) -> np.ndarray:
    """icarogw's _lowpass_filter: 1 below mmax - delta, a smooth fall over delta, 0 above mmax."""
    out = np.where(m <= mmax - delta, 1.0, 0.0)
    w = (m < mmax) & (m > mmax - delta)
    x = mmax - m[w]
    with np.errstate(over="ignore"):
        out[w] = 1.0 / (np.exp(delta / x + delta / (x - delta)) + 1.0)
    return out


def _ln_power_law(m: np.ndarray, lo: float, hi: float, slope: float) -> np.ndarray:
    """ln of the normalized m^slope on [lo, hi]."""
    norm = np.log(hi / lo) if slope == -1 else (hi ** (slope + 1) - lo ** (slope + 1)) / (slope + 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where((m >= lo) & (m <= hi), slope * np.log(m) - np.log(norm), -np.inf)


def _ln_truncated_gaussian(m: np.ndarray, mu: float, sigma: float, lo: float, hi: float) -> np.ndarray:
    from math import erf, sqrt

    norm = 0.5 * (erf((hi - mu) / (sigma * sqrt(2))) - erf((lo - mu) / (sigma * sqrt(2))))
    return np.where((m >= lo) & (m <= hi),
                    -0.5 * ((m - mu) / sigma) ** 2 - np.log(sigma * np.sqrt(2 * np.pi) * norm), -np.inf)


def _fullpop4_mbreak(p: dict) -> float:
    return 0.5 * (p["leftdip"] + p["leftdipsmooth"] + p["rightdip"] - p["rightdipsmooth"])


def ln_mass_1d_fullpop4(m: np.ndarray, p: dict = FULLPOP4_GWTC4) -> np.ndarray:
    """ln p_S(m) of FullPop-4.0 up to a constant: a broken power law (break at m_break, between the dip edges) and
    two Gaussians, times the low- and high-mass windows and the dip (icarogw's SmoothedPlusDipProb)."""
    m = np.asarray(m, float)
    lo, hi, b = p["mmin"], p["mmax"], _fullpop4_mbreak(p)
    with np.errstate(divide="ignore"):
        ln_pl1 = _ln_power_law(m, lo, b, -p["alpha_1"])
        ln_pl2 = _ln_power_law(m, b, hi, -p["alpha_2"])
        bb = np.array([b])
        ln_join = (_ln_power_law(bb, lo, b, -p["alpha_1"]) - _ln_power_law(bb, b, hi, -p["alpha_2"]))[0]
        ln_bpl = np.logaddexp(ln_pl1, ln_pl2 + ln_join) - np.log1p(np.exp(ln_join))
        lam, lam_low = p["lambda_g"], p["lambda_g_low"]
        ln_mix = np.logaddexp.reduce([
            np.log1p(-lam) + ln_bpl,
            np.log(lam * lam_low) + _ln_truncated_gaussian(m, p["mu_g_low"], p["sigma_g_low"], lo,
                                                           p["mu_g_low"] + 5 * p["sigma_g_low"]),
            np.log(lam * (1 - lam_low)) + _ln_truncated_gaussian(m, p["mu_g_high"], p["sigma_g_high"], lo,
                                                                 p["mu_g_high"] + 5 * p["sigma_g_high"]),
        ])
        notch = 1.0 - p["deep"] * _highpass(m, p["leftdip"], p["leftdipsmooth"]) * _lowpass(m, p["rightdip"],
                                                                                             p["rightdipsmooth"])
        window = _highpass(m, lo, p["bottomsmooth"]) * _lowpass(m, hi, p["topsmooth"]) * notch
        return np.where((m >= lo) & (m <= hi), ln_mix + np.log(window), -np.inf)


def ln_mass_fullpop4(m1: np.ndarray, m2: np.ndarray, p: dict = FULLPOP4_GWTC4) -> np.ndarray:
    """ln p(m1, m2) of FullPop-4.0 up to a constant: p_S(m1) p_S(m2) times the pairing function q^β, with
    β = beta_bottom when m2 <= m_break, beta_top above, and m2 <= m1."""
    m1, m2 = np.broadcast_arrays(np.asarray(m1, float), np.asarray(m2, float))
    q = m2 / m1
    beta = np.where(m2 <= _fullpop4_mbreak(p), p["beta_bottom"], p["beta_top"])
    with np.errstate(divide="ignore", invalid="ignore"):
        ln_pair = np.where(q <= 1, beta * np.log(q), -np.inf)
    return ln_mass_1d_fullpop4(m1, p) + ln_mass_1d_fullpop4(m2, p) + ln_pair


def ln_rate_madau(z: np.ndarray, gamma: float, kappa: float, zp: float) -> np.ndarray:
    """ln ψ(z) of the Madau-Dickinson rate, normalized to 1 at z = 0 (icarogw's md_rate)."""
    z = np.asarray(z, float)
    return (np.log1p((1 + zp) ** (-gamma - kappa)) + gamma * np.log1p(z)
            - np.log1p(((1 + z) / (1 + zp)) ** (gamma + kappa)))


@dataclass(frozen=True)
class Population:
    """Mass and redshift density of the sources, beyond uniform in comoving volume and source-frame time."""
    ln_mass: Callable          # ln p(m1, m2) in the source frame
    ln_rate: Callable          # ln ψ(z)
    description: str


POPULATIONS: dict[str, Optional[Population]] = {
    "volume": None,
    "fullpop4": Population(
        ln_mass_fullpop4, lambda z: ln_rate_madau(z, **MADAU_GWTC4),
        "FullPop-4.0 masses and the Madau–Dickinson rate (γ = 3.57, κ = 2.97, z<sub>p</sub> = 2.59), fixed to "
        "the medians of the GWTC-4.0 spectral siren, as the GWTC-4.0 reanalysis of GW170817"),
}


def _log(msg: str) -> None:
    print(f"[hubble_constant] {msg}", flush=True)


# ---------------------------------------------------------------------------
# likelihood
# ---------------------------------------------------------------------------
def _logsumexp(a: np.ndarray, axis: int) -> np.ndarray:
    m = np.max(a, axis=axis, keepdims=True)
    m = np.where(np.isfinite(m), m, 0.0)
    with np.errstate(divide="ignore"):
        return np.squeeze(m, axis) + np.log(np.sum(np.exp(a - m), axis=axis))


def event_ln_likelihood(h0: np.ndarray, distances: np.ndarray, pe_prior: np.ndarray, z_obs: float, sigma_z: float,
                        weights: Optional[np.ndarray] = None, rows: int = 64,
                        population: Optional[Population] = None,
                        masses: Optional[tuple[np.ndarray, np.ndarray]] = None) -> np.ndarray:
    """ln ∫ N(z_obs; z, σ_z) p_pop(z) L_GW(d_L(z, H0)) dz on the H0 grid, up to a constant (no selection).
    `weights`: an independent likelihood of each sample (e.g. a viewing-angle constraint). With a `population`,
    `masses` are the detector-frame (mass_1, mass_2) of the samples, and each sample also carries
    p_m(m1/(1+z), m2/(1+z)) ψ(z) / (1+z)^2."""
    if population is not None:
        m1d, m2d = (np.asarray(m, float)[None, :] for m in masses)
    d = np.asarray(distances, float)[None, :]
    ln_w0 = -np.log(np.asarray(pe_prior, float))[None, :]
    if weights is not None:
        with np.errstate(divide="ignore"):
            ln_w0 = ln_w0 + np.log(np.asarray(weights, float))[None, :]
    n = d.shape[1]
    out = np.empty(len(h0))
    for i in range(0, len(h0), rows):
        h = np.asarray(h0[i:i + rows], float)[:, None]
        z, pz, dd = _at_dimensionless_distance(d * h / C_KMS)
        with np.errstate(divide="ignore"):
            ln_w = ln_w0 + np.log(pz) - np.log(dd) + np.log(h / C_KMS)
            if population is not None:
                ln_w = ln_w + (population.ln_rate(z) + population.ln_mass(m1d / (1 + z), m2d / (1 + z))
                               - 2 * np.log1p(z))
        bw = 1.06 * np.std(z, axis=1, keepdims=True) * n ** -0.2
        s = np.maximum(sigma_z, bw)
        out[i:i + rows] = _logsumexp(ln_w - 0.5 * ((z_obs - z) / s) ** 2 - np.log(s), axis=1) - np.log(n)
    return out


# ---------------------------------------------------------------------------
# selection
# ---------------------------------------------------------------------------
def ln_mass_bns(m1: np.ndarray, m2: np.ndarray, lo: float = 1.0, hi: float = 2.5) -> np.ndarray:
    """Neutron-star masses uniform in [lo, hi], m2 <= m1."""
    ok = (m2 >= lo) & (m2 <= m1) & (m1 <= hi)
    return np.where(ok, np.log(2.0 / (hi - lo) ** 2), -np.inf)


def ln_mass_bbh(m1: np.ndarray, m2: np.ndarray) -> np.ndarray:
    """GWTC-3 Power Law + Peak (the mass model of the rates mode)."""
    from .catalogs import _ln_power_law_peak

    return _ln_power_law_peak(m1, m2)


MASS_POPULATIONS: dict[str, Callable] = {"BNS": ln_mass_bns, "BBH": ln_mass_bbh}


def ln_selection_euclidean(h0: np.ndarray) -> np.ndarray:
    """ln β for a GW-limited detection of nearby sources: the detectable redshifts scale as H0."""
    return 3 * np.log(np.asarray(h0, float))


def ln_selection_injections(h0: np.ndarray, inj: dict, ln_pop_mass: Callable, n_grid: int = 77,
                            ln_rate: Optional[Callable] = None) -> tuple[np.ndarray, np.ndarray]:
    """(ln β, effective injections) on the H0 grid, from found injections in the detector frame.

    `inj` has detector-frame mass_1, mass_2, luminosity_distance, their draw density `prior` in these
    variables, and `ntotal` (as `hubble_constant.detector_frame_injections` returns). At each trial H0 the
    injections are carried to the source frame and weighted by the population density in the detector frame,
    p_m(m1_s, m2_s) p_pop(z) / [(1+z)^2 dd_L/dz], times the rate ψ(z) when `ln_rate` is given.
    """
    hg = np.linspace(h0[0], h0[-1], n_grid)
    ln_b, neff = np.empty(n_grid), np.empty(n_grid)
    d, m1, m2, prior = (np.asarray(inj[k], float) for k in ("luminosity_distance", "mass_1", "mass_2", "prior"))
    for i, h in enumerate(hg):
        z, pz, dd = _at_dimensionless_distance(d * h / C_KMS)
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            lnp = (ln_pop_mass(m1 / (1 + z), m2 / (1 + z)) + np.log(pz) - 2 * np.log1p(z)
                   - np.log(dd) + np.log(h / C_KMS))
            if ln_rate is not None:
                lnp = lnp + ln_rate(z)
            w = np.exp(lnp) / prior
        w = np.where(np.isfinite(w), w, 0.0)
        s = w.sum()
        ln_b[i] = np.log(s / inj["ntotal"]) if s > 0 else -np.inf
        neff[i] = s * s / (w * w).sum() if s > 0 else 0.0
    return np.interp(h0, hg, ln_b), np.interp(h0, hg, neff)


def load_injections(cp: Counterpart, sensitivity_release: Optional[str], sensitivity_file, far_threshold: float,
                    snr_threshold: float) -> dict:
    from . import hubble_constant as hc

    release = sensitivity_release or hc.H0_DEFAULT_RELEASE
    path = hc.h0_sensitivity_path(sensitivity_file, release)
    inj = hc.detector_frame_injections(path, far_threshold, snr_threshold, runs=[cp.run])
    _log(f"{len(inj['prior'])} found injections in {cp.run} ({path.name})")
    return inj


# ---------------------------------------------------------------------------
# summaries
# ---------------------------------------------------------------------------
def _normalize(h0: np.ndarray, p: np.ndarray) -> np.ndarray:
    return p / _trapz(p, h0)


def posterior_from_ln(h0: np.ndarray, ln_p: np.ndarray) -> np.ndarray:
    ln_p = np.where(np.isfinite(ln_p), ln_p, -np.inf)
    if not np.isfinite(ln_p).any():
        raise ValueError("the H0 posterior is zero everywhere on the grid (no effective injections?)")
    return _normalize(h0, np.exp(ln_p - ln_p.max()))


def summarize(h0: np.ndarray, p: np.ndarray) -> dict:
    """Maximum a posteriori, 68.3% highest-density interval, median and 90% equal-tailed interval."""
    p = _normalize(h0, p)
    dx = np.gradient(h0)
    order = np.argsort(-p)
    n = np.searchsorted(np.cumsum(p[order] * dx[order]), 0.683) + 1
    hpd = h0[order[:n]]
    cdf = np.cumsum(p * dx)
    cdf /= cdf[-1]
    q = np.interp([0.05, 0.5, 0.95], cdf, h0)
    return dict(map=float(h0[np.argmax(p)]), hpd68_low=float(hpd.min()), hpd68_high=float(hpd.max()),
                median=float(q[1]), low_90=float(q[0]), high_90=float(q[2]))


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
@default_style
def _plot(h0: np.ndarray, curves: list[tuple[str, np.ndarray, str, str]], published: Optional[dict], ref: str,
          title: str, out_png: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid, orange = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df", "#eb6834"
    fig, ax = plt.subplots(figsize=(8.2, 4.2), dpi=150)
    fig.patch.set_facecolor(surface); ax.set_facecolor(surface)
    if published:
        ax.axvspan(published["lo68"], published["hi68"], color=orange, alpha=0.12, zorder=1, linewidth=0)
        ax.axvline(published["map"], color=orange, lw=2, zorder=3,
                   label=f"{ref}: {published['map']:.0f}, 68% {published['lo68']:.0f}–{published['hi68']:.0f}")
    for label, p, color, ls in curves:
        ax.plot(h0, _normalize(h0, p), color=color, lw=2, ls=ls, zorder=4, label=label)
    for x, t, ha, dx in ((67.4, "Planck", "right", -3), (73.0, "SH0ES", "left", 3)):
        ax.axvline(x, color=ink2, lw=0.8, ls=(0, (3, 3)), zorder=1)
        ax.annotate(t, xy=(x, 1), xycoords=("data", "axes fraction"), xytext=(dx, -12), textcoords="offset points",
                    fontsize=8, color=ink2, ha=ha)
    ax.set_xlim(h0[0], h0[-1])
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
    labels = [c[0] for c in curves] + ([ref] if published else [])
    ncol = 2 if max(len(t) for t in labels) < 48 else 1
    ax.set_title(title, color=ink, fontsize=10.5, loc="left", pad=20 + 13 * -(-len(labels) // ncol))
    leg = ax.legend(frameon=False, fontsize=8.5, labelcolor=ink2, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=ncol,
                    borderaxespad=0.3, handlelength=2.6)
    for h in leg.legend_handles:
        h.set_linewidth(2)                 # gwpy, once imported, widens the legend lines
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def _fmt(s: dict) -> str:
    return (f"{s['map']:.1f} (+{s['hpd68_high'] - s['map']:.1f} / −{s['map'] - s['hpd68_low']:.1f}) km/s/Mpc "
            f"(maximum a posteriori, 68%); median {s['median']:.1f}, 90%: {s['low_90']:.1f}–{s['high_90']:.1f}")


# ---------------------------------------------------------------------------
# method
# ---------------------------------------------------------------------------
BRIGHT_FILE = "bright.json"
GRID_FILE = "posterior_grid.tsv"


def run_h0_bright(
    event: str = "GW170817",
    workdir: str | Path = "hubble_constant_bright",
    out_report_html: Optional[str | Path] = "hubble_constant_bright.html",
    out_summary_tsv: Optional[str | Path] = "hubble_constant_bright.tsv",
    pe_labels: Optional[Iterable[str]] = None,
    pe_file: Optional[str | Path] = None,
    cache_dir: str | Path = ".cache_gwosc",
    pe_cache: Optional[str | Path] = None,
    v_recession: Optional[tuple[float, float]] = None,
    v_peculiar: Optional[tuple[float, float]] = None,
    redshift: Optional[tuple[float, float]] = None,
    sky_radius_deg: float = 3.0,
    viewing_angle_constraint: Optional[tuple[float, float]] = None,
    selection: str = "auto",
    sensitivity_release: Optional[str] = None,
    sensitivity_file: Optional[str | Path] = None,
    far_threshold: float = 0.25,
    snr_threshold: float = 10.0,
    h0_range: tuple[float, float] = H0_PRIOR,
    population: str = "volume",
    pe_distance_prior: Optional[str] = None,
) -> pd.DataFrame:
    """Bright-siren H0 of `event` for each of its PE labels, in `workdir` (the posterior grid, ``bright.json``, plots).
    The first label is the result; `hubble_constant --method joint` combines it with other methods."""
    if selection not in ("auto", "euclidean", "injections"):
        raise ValueError(f"unknown selection {selection!r}; choose auto, euclidean or injections")
    if population not in POPULATIONS:
        raise ValueError(f"unknown population {population!r}; choose {', '.join(POPULATIONS)}")
    pop = POPULATIONS[population]
    lo, hi = h0_range
    if not 0 < lo < hi or sky_radius_deg <= 0:
        raise ValueError("the H0 range must be 0 < MIN < MAX and the sky radius > 0")
    if viewing_angle_constraint is not None and not (0 <= viewing_angle_constraint[0] <= 90
                                                     and viewing_angle_constraint[1] > 0):
        raise ValueError("the viewing-angle constraint is MEAN SIGMA in degrees, 0 <= MEAN <= 90, SIGMA > 0")
    from . import hubble_constant as hc

    cp = get_counterpart(event)
    z_obs, sigma_z, z_desc = hubble_flow_redshift(cp, v_recession, v_peculiar, redshift)
    defaults = (v_recession is None and v_peculiar is None and redshift is None and viewing_angle_constraint is None
                and pop is None and pe_distance_prior is None)
    _log(f"bright siren {event}, {cp.host}: {z_desc.replace('<sub>', '').replace('</sub>', '')}")
    if cp.candidate:
        _log(f"WARN: {cp.transient} is a candidate counterpart; the association is not established")
    pe_path, samples = load_event_samples(event, cp, pe_file, pe_labels, cache_dir, pe_cache, log_cb=_log)

    h0 = np.linspace(lo, hi, int(round((hi - lo) / 0.05)) + 1)
    sel = selection if selection != "auto" else (
        "euclidean" if z_obs < EUCLIDEAN_ZMAX and pop is None else "injections")
    if sel == "euclidean":
        ln_beta, sel_note = ln_selection_euclidean(h0), (
            "Selection: GW-limited detection of nearby sources, β(H<sub>0</sub>) ∝ H<sub>0</sub><sup>3</sup>, which "
            "cancels the volume factor of the population.")
    else:
        inj = load_injections(cp, sensitivity_release, sensitivity_file, far_threshold, snr_threshold)
        if pop is None:
            ln_beta, neff = ln_selection_injections(h0, inj, MASS_POPULATIONS[cp.population])
            mass_desc = "Power Law + Peak" if cp.population == "BBH" else "uniform 1–2.5 M☉"
            pop_desc = f"the {cp.population} population ({mass_desc}, uniform in comoving volume)"
        else:
            ln_beta, neff = ln_selection_injections(h0, inj, pop.ln_mass, ln_rate=pop.ln_rate)
            pop_desc = "the population"
        sel_note = (f"Selection: {len(inj['prior'])} found injections of {cp.run} (FAR &lt; {far_threshold:g}/yr), "
                    f"reweighted at each H<sub>0</sub> to {pop_desc}; effective injections "
                    f"{neff.min():.0f}–{neff.max():.0f}.")
        if neff.min() < MIN_NEFF_INJECTIONS:
            sel_note += f" <b>Fewer than {MIN_NEFF_INJECTIONS} effective injections: β is noisy.</b>"
            _log(f"WARN: only {neff.min():.0f} effective injections")
    _log(f"selection: {sel}; population: {population}")

    rows, posts, notes = [], {}, []
    for label, s in samples.items():
        keep, how = sky_conditioned(s, cp, sky_radius_deg)
        d = np.asarray(s["luminosity_distance"], float)[keep]
        prior, prior_desc = hc.pe_distance_prior(d, s["prior_desc"], cp.run, kind=pe_distance_prior)
        masses = None
        if pop is not None:
            if "mass_1" not in s or "mass_2" not in s:
                raise ValueError(f"{label}: no detector-frame masses (mass_1, mass_2) in the samples for "
                                 f"--population {population}")
            masses = (np.asarray(s["mass_1"], float)[keep], np.asarray(s["mass_2"], float)[keep])
        w, w_note = None, ""
        if viewing_angle_constraint is not None:
            view = viewing_angle(s)
            if view is None:
                raise ValueError(f"{label}: no viewing angle (theta_jn or iota) in the samples for the constraint")
            w = viewing_angle_weights(view[keep], viewing_angle_constraint)
            w_note = f"; {w.sum() ** 2 / (w * w).sum():.0f} effective samples with the viewing-angle constraint"
        p = posterior_from_ln(h0, event_ln_likelihood(h0, d, prior, z_obs, sigma_z, weights=w, population=pop,
                                                          masses=masses) - ln_beta)
        posts[label] = p
        dq = np.percentile(d, [5, 50, 95])
        rows.append(dict(analysis=f"bright siren {event}, {label}", n_samples=len(d), distance_median=dq[1],
                         distance_low_90=dq[0], distance_high_90=dq[2], **summarize(h0, p)))
        notes.append(f"{label}: {how}{w_note}; d<sub>L</sub> = {dq[1]:.1f} Mpc (90%: {dq[0]:.1f}–{dq[2]:.1f}); "
                     f"PE distance prior {prior_desc[:80]}.")
        _log(f"{label}: H0 = {_fmt(rows[-1])}")

    workdir = Path(workdir).expanduser()
    workdir.mkdir(parents=True, exist_ok=True)
    main = next(iter(posts))
    grid = pd.DataFrame({"H0": h0, "p": posts[main], **{f"p_{k}": v for k, v in posts.items()}})
    grid.to_csv(workdir / GRID_FILE, sep="\t", index=False, float_format="%.6g")
    (workdir / BRIGHT_FILE).write_text(json.dumps(dict(
        method="bright", event=event, host=cp.host, transient=cp.transient, candidate=cp.candidate,
        events=[e for e in (event, cp.pe_event) if e], label=main, labels=list(posts), pe_file=str(pe_path), z=z_obs, sigma_z=sigma_z, selection=sel,
        population=population, pe_distance_prior=pe_distance_prior, h0_range=[lo, hi], viewing_angle=list(viewing_angle_constraint) if viewing_angle_constraint else None,
        summary={r["analysis"]: {k: r[k] for k in ("map", "hpd68_low", "hpd68_high", "median", "low_90", "high_90")}
                 for r in rows}), indent=1))
    table = pd.DataFrame(rows)
    if out_summary_tsv:
        Path(out_summary_tsv).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_summary_tsv, sep="\t", index=False, float_format="%.6g")
    if out_report_html:
        published = cp.published if defaults else None
        colors = ("#2a78d6", "#1a9e77", "#7a5195")
        curves = [(f"{k}: {r['map']:.0f}, 68% {r['hpd68_low']:.0f}–{r['hpd68_high']:.0f}", posts[k],
                   colors[i % len(colors)], "-" if i == 0 else "--") for i, (k, r) in enumerate(zip(posts, rows))]
        images = [_plot(h0, curves, published, cp.ref, f"Hubble constant from {event} and {cp.host} (bright siren)",
                        workdir / "plots" / "h0_posterior.png")]
        paras = [
            f"H<sub>0</sub> = <b>{_fmt(rows[0])}</b> from {event} ({main}) and its "
            f"{'candidate ' if cp.candidate else ''}host {cp.host} ({cp.transient}): the luminosity distance comes "
            "from the GW signal alone, the redshift from the host.",
        ]
        if cp.caveat:
            paras.append(f"<b>Caveat.</b> {cp.caveat}")
        paras += [
            f"Host: {z_desc}. Flat ΛCDM (Ω<sub>m</sub> = {OM0}), flat H<sub>0</sub> prior {lo:g}–{hi:g} km/s/Mpc, "
            + ("sources uniform in comoving volume. " if pop is None else f"Population: {pop.description}. ")
            + sel_note,
            *notes,
        ]
        if viewing_angle_constraint is not None:
            paras.append(f"The viewing angle is constrained to {viewing_angle_constraint[0]:g} ± "
                         f"{viewing_angle_constraint[1]:g}° (Gaussian, independent of the GW signal): it breaks part "
                         "of the distance–inclination degeneracy, which dominates the uncertainty.")
        else:
            paras.append("The distance is degenerate with the inclination of the orbit, which dominates the "
                         "uncertainty: the low-distance tail, hence the high-H<sub>0</sub> tail, comes from inclined "
                         "orbits. The <code>counterpart</code> mode shows this degeneracy; an independent constraint "
                         "on the viewing angle (<code>--viewing-angle</code>) narrows it.")
        if published:
            paras.append(f"Published ({cp.ref}): {published['map']:.1f} (+{published['hi68'] - published['map']:.1f} / "
                         f"−{published['map'] - published['lo68']:.1f}) km/s/Mpc. It used the distance posterior of the "
                         "2017 analysis and a linear Hubble law; the public GWTC-1 samples are a later reanalysis, "
                         "whose longer low-distance tail widens the upper side of the H<sub>0</sub> interval.")
        paras.append(f"Work directory: {workdir} (posterior grid {GRID_FILE}). Combine with the spectral or dark siren: "
                     f"<code>hubble_constant --method joint --inputs {workdir.name} DIR</code>.")
        Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
        write_simple_html_report(out_report_html, title=f"Hubble constant (bright siren, {event})", paragraphs=paras,
                                 images=images,
                                 tables=[("H0 posterior", table.to_html(index=False, float_format=lambda x: f"{x:.3g}",
                                                                        na_rep=""))])
        _log(f"report written to {out_report_html}")
    return table

"""Redshift evolution of the merger rate: the Madau-Dickinson shape fitted by the spectral siren.

icarogw's `rateevolution_Madau` (the `hubble_constant` mode) uses, normalized to 1 at z = 0,

    psi(z) = [1 + (1 + z_p)^(-gamma-kappa)] (1 + z)^gamma / [1 + ((1 + z) / (1 + z_p))^(gamma + kappa)],

the merger rate per comoving volume and source-frame time being R(z) = R_0 psi(z). The spectral-siren runs
marginalize over R_0 (scale-free likelihood), so their posterior gives the shape R(z)/R_0 only; the
absolute rate comes from the `rates` mode. The cosmic star-formation history of Madau & Dickinson 2014
(arXiv:1403.0007) has the same form with gamma = 2.7, kappa = 2.9, z_p = 1.9. The rate peaks at
z_peak = (1 + z_p) (gamma / kappa)^(1 / (gamma + kappa)) - 1.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

MD14 = dict(gamma=2.7, kappa=2.9, zp=1.9)
PRIOR = dict(gamma=(0.0, 12.0), kappa=(0.0, 6.0), zp=(0.0, 4.0))     # those of h0_icarogw.py


def md_shape(z, gamma, kappa, zp):
    """psi(z), normalized to 1 at z = 0, as an array (n_draws, n_z) for n_draws values of gamma, kappa, zp."""
    z = np.atleast_1d(np.asarray(z, float))
    g, k, p = (np.atleast_1d(np.asarray(x, float))[:, None] for x in (gamma, kappa, zp))
    return (1 + (1 + p) ** (-g - k)) * (1 + z) ** g / (1 + ((1 + z) / (1 + p)) ** (g + k))


def z_peak(gamma, kappa, zp):
    gamma, kappa, zp = (np.asarray(x, float) for x in (gamma, kappa, zp))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(kappa > 0, (1 + zp) * (gamma / kappa) ** (1 / (gamma + kappa)) - 1, np.inf)


def summary(post: pd.DataFrame) -> dict:
    """Quantiles of gamma, kappa, z_p, z_peak and R(1)/R(0), R(2)/R(0)."""
    g, k, p = (post[c].to_numpy(float) for c in ("gamma", "kappa", "zp"))
    out = {c: np.percentile(post[c], [5, 50, 95]) for c in ("gamma", "kappa", "zp")}
    out["z_peak"] = np.percentile(z_peak(g, k, p), [5, 50, 95])
    for zz in (1.0, 2.0):
        out[f"R({zz:g})/R(0)"] = np.percentile(md_shape(zz, g, k, p)[:, 0], [5, 50, 95])
    return out


def plot(post: pd.DataFrame, out_png: Path, z_max_events: Optional[float] = None, seed: int = 1) -> Path:
    """R(z)/R(0): posterior median and 90% band, the prior band, and the star-formation history."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid, blue, orange, grey = ("#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6", "#eb6834",
                                                   "#b8b6b0")
    z = np.linspace(0, 4, 201)
    rng = np.random.default_rng(seed)
    n = min(len(post), 4000)
    d = post.sample(n, random_state=seed) if len(post) > n else post
    psi = md_shape(z, d["gamma"], d["kappa"], d["zp"])
    pri = md_shape(z, *(rng.uniform(*PRIOR[c], 4000) for c in ("gamma", "kappa", "zp")))
    fig, ax = plt.subplots(figsize=(8.2, 4.4), dpi=150)
    fig.patch.set_facecolor(surface); ax.set_facecolor(surface)
    lo, hi = np.percentile(pri, [5, 95], axis=0)
    ax.fill_between(z, lo, hi, color=grey, alpha=0.25, lw=0, zorder=1, label="prior, 90%")
    lo, med, hi = np.percentile(psi, [5, 50, 95], axis=0)
    ax.fill_between(z, lo, hi, color=blue, alpha=0.3, lw=0, zorder=2, label="spectral-siren posterior, 90%")
    ax.plot(z, med, color=blue, lw=2, zorder=3, label="median")
    ax.plot(z, md_shape(z, **MD14)[0], color=orange, lw=1.8, ls=(0, (4, 2)), zorder=4,
            label="star formation (Madau & Dickinson 2014)")
    if z_max_events:
        ax.axvspan(z_max_events, z[-1], color=grey, alpha=0.12, lw=0, zorder=0)
        ax.annotate("beyond the detected events", xy=(z_max_events, 1), xycoords=("data", "axes fraction"),
                    xytext=(4, -12), textcoords="offset points", fontsize=8, color=ink2)
    ax.set_yscale("log")
    ax.set_ylim(0.5, 300)
    ax.set_xlim(0, z[-1])
    ax.set_xlabel("Redshift z", color=ink2)
    ax.set_ylabel("R(z) / R(0)", color=ink2)
    ax.grid(color=grid, lw=0.7, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(grid)
    ax.tick_params(colors=ink2, labelsize=9)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=ink2, loc="upper left")
    ax.set_title("BBH merger-rate evolution (Madau–Dickinson shape, fitted with H$_0$)", color=ink, fontsize=10.5,
                 loc="left")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def report_paragraph(s: dict, z_max_events: Optional[float] = None) -> str:
    f = lambda q, d=1: f"{q[1]:.{d}f} (90%: {q[0]:.{d}f}–{q[2]:.{d}f})"
    return ("<b>Merger-rate evolution.</b> The Madau–Dickinson shape fitted together with H<sub>0</sub>: "
            f"γ = {f(s['gamma'])}, κ = {f(s['kappa'])}, z<sub>p</sub> = {f(s['zp'])}; the rate grows by "
            f"R(1)/R(0) = {f(s['R(1)/R(0)'])}. The star-formation history has γ = 2.7, κ = 2.9, z<sub>p</sub> = 1.9 "
            f"(R(1)/R(0) = {md_shape(1.0, **MD14)[0, 0]:.1f}). "
            + (f"The detected events reach z ≈ {z_max_events:.1f}: beyond it the shape (and κ, z<sub>p</sub>) follows "
               "the prior. " if z_max_events else "")
            + "The likelihood is scale-free: the local rate R(0) itself comes from the <code>rates</code> mode.")

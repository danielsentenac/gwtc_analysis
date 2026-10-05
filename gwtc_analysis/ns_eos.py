"""Joint neutron-star equation of state from GW170817 and GW190425: one tidal deformability curve for all.

All neutron stars follow the same equation of state (EOS), so the four stars of the two binary neutron
stars lie on one curve Lambda(m). With stars of nearly equal radius, Lambda ∝ (R / m)^5 k_2 gives
Lambda(m) = Lambda_1.4 (m / 1.4 M_sun)^-6 (De et al. 2018, arXiv:1804.08583, the Lambda_1 = q^6 Lambda_2
prescription), so the EOS is summarized by Lambda_1.4, the deformability of a 1.4 M_sun star.

The data measure mostly the mass-weighted combination

    Lambda~ = (16/13) [(m1 + 12 m2) m1^4 Lambda_1 + (m2 + 12 m1) m2^4 Lambda_2] / (m1 + m2)^5 = a Lambda_1 + b Lambda_2.

For each event the likelihood of Lambda_1.4 is estimated from its PE samples (m1_i, m2_i, Lambda~_i),

    L(Lambda_1.4) ∝ Σ_i K_h(Lambda~_i - Lambda~_model,i) / pi(Lambda~_i | m_i),

with Lambda~_model,i the combination predicted at the masses of sample i, K_h a Gaussian kernel (reflected
at 0) and pi the prior of Lambda~ given the masses implied by the PE priors on Lambda_1 and Lambda_2
(uniform on [0, L], a trapezoid). The mass population is taken as the PE mass prior. The events multiply; the prior on
Lambda_1.4 is flat. The radius follows from the empirical relation Lambda_1.4 = 2.88e-6 (R_1.4 / km)^7.5
of Annala et al. 2018 (arXiv:1711.02644), approximate (~0.5 km).

Published (GW170817, common EOS, LVK 2018, arXiv:1805.11581): Lambda_1.4 = 190 (+390 / -120), 90%.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from .report import write_simple_html_report

LAMBDA_MAX = 5000.0               # uniform PE priors on Lambda_1, Lambda_2 (GWTC-1, GWTC-2.1)
GRID_MAX = 3000.0
PUBLISHED = dict(median=190.0, lo90=70.0, hi90=580.0, ref="LVK 2018, arXiv:1805.11581 (GW170817)")
EVENTS = {
    "GW170817": dict(pe_event=None, labels=dict(low="C02:IMRPhenomPv2_NRTidal-LowSpin",
                                                high="C02:IMRPhenomPv2_NRTidal-HighSpin")),
    "GW190425": dict(pe_event="GW190425_081805", labels=dict(low="C01:IMRPhenomPv2_NRTidal:LowSpin",
                                                              high="C01:IMRPhenomPv2_NRTidal:HighSpin")),
}
_trapz = getattr(np, "trapezoid", None) or np.trapz


def _log(msg: str) -> None:
    print(f"[ns_eos] {msg}", flush=True)


# ---------------------------------------------------------------------------
# tidal relations
# ---------------------------------------------------------------------------
def tilde_coefficients(m1, m2) -> tuple[np.ndarray, np.ndarray]:
    """(a, b) with Lambda~ = a Lambda_1 + b Lambda_2."""
    m1, m2 = np.asarray(m1, float), np.asarray(m2, float)
    m = m1 + m2
    return 16 / 13 * (m1 + 12 * m2) * m1 ** 4 / m ** 5, 16 / 13 * (m2 + 12 * m1) * m2 ** 4 / m ** 5


def lambda_of_m(m, lambda_14):
    """Lambda(m) = Lambda_1.4 (m / 1.4)^-6, broadcast as (n_lambda, n_m)."""
    return np.atleast_1d(np.asarray(lambda_14, float))[:, None] * (np.asarray(m, float)[None, :] / 1.4) ** -6


def radius_14(lambda_14):
    """R_1.4 (km) from Lambda_1.4 = 2.88e-6 (R / km)^7.5 (Annala et al. 2018)."""
    return (np.asarray(lambda_14, float) / 2.88e-6) ** (1 / 7.5)


def prior_tilde(y, a, b, lam_max: float = LAMBDA_MAX):
    """Density of y = a U1 + b U2, U1, U2 uniform on [0, lam_max]: a trapezoid."""
    y = np.asarray(y, float)
    A, B = np.minimum(a, b) * lam_max, np.maximum(a, b) * lam_max
    return np.where(y < 0, 0.0, np.where(y < A, y / (A * B), np.where(y < B, 1 / B,
                    np.where(y <= A + B, (A + B - y) / (A * B), 0.0))))


def event_ln_likelihood(grid: np.ndarray, m1: np.ndarray, m2: np.ndarray, lam_tilde: np.ndarray,
                        lam_max: float = LAMBDA_MAX, bw_factor: float = 0.5, rows: int = 50) -> np.ndarray:
    """ln L(Lambda_1.4) on the grid for one event, from its PE samples (source-frame masses).

    Each sample is first reweighted to a flat prior on Lambda~ (weight 1 / pi(Lambda~_i | m_i), finite at the
    sample), then smoothed with a Gaussian kernel reflected at Lambda~ = 0 and evaluated at the prediction
    of the model at its masses. Dividing a smoothed density by the prior at the model point instead diverges
    at small Lambda_1.4, where the prior of Lambda~ vanishes.
    """
    m1, m2, lt = (np.asarray(x, float) for x in (m1, m2, lam_tilde))
    a, b = tilde_coefficients(m1, m2)
    pri = prior_tilde(lt, a, b, lam_max)
    ok = pri > 0
    m1, m2, lt, a, b, w = m1[ok], m2[ok], lt[ok], a[ok], b[ok], 1.0 / pri[ok]
    w = np.minimum(w, np.percentile(w, 99.5))          # a few samples at Lambda~ ≈ 0 carry huge weights
    h = bw_factor * 1.06 * np.std(lt) * len(lt) ** -0.2
    out = np.empty(len(grid))
    for i in range(0, len(grid), rows):
        g = grid[i:i + rows]
        l1, l2 = lambda_of_m(m1, g), lambda_of_m(m2, g)
        model = a[None, :] * l1 + b[None, :] * l2
        inside = (l1 <= lam_max) & (l2 <= lam_max)
        k = np.exp(-0.5 * ((lt[None, :] - model) / h) ** 2) + np.exp(-0.5 * ((lt[None, :] + model) / h) ** 2)
        s = np.where(inside, w[None, :] * k, 0.0).sum(axis=1)
        with np.errstate(divide="ignore"):
            out[i:i + rows] = np.log(s / w.sum())
    return out


def normalize(grid, ln_l):
    ln_l = np.where(np.isfinite(ln_l), ln_l, -np.inf)
    p = np.exp(ln_l - ln_l.max())
    return p / _trapz(p, grid)


def summary(grid, p) -> dict:
    cdf = np.cumsum(p * np.gradient(grid)); cdf /= cdf[-1]
    q = np.interp([0.05, 0.5, 0.95], cdf, grid)
    r = radius_14(q)
    return dict(median=q[1], low_90=q[0], high_90=q[2], upper_90=float(np.interp(0.9, cdf, grid)),
                r_median=r[1], r_low_90=r[0], r_high_90=r[2])


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def load_event(name: str, spin_prior: str, cache_dir, pe_cache) -> dict:
    import h5py

    from . import hubble_constant as hc

    cfg = EVENTS[name]
    if cfg["pe_event"] is None:
        from .unofficial_pe import build_unofficial_pe_bundle

        path = build_unofficial_pe_bundle(name, cache_dir=cache_dir, log_cb=_log)
        if path is None:
            raise ValueError(f"could not build the PE bundle of {name}")
    else:
        from .bright_siren import _event_pe_file

        path = _event_pe_file(cfg["pe_event"], Path(pe_cache).expanduser() if pe_cache else hc.default_pe_cache())
    label = cfg["labels"][spin_prior]
    with h5py.File(path, "r") as h:
        ps = h[label]["posterior_samples"][()]
    m1, m2 = ps["mass_1_source"], ps["mass_2_source"]
    swap = m2 > m1
    l1, l2 = np.where(swap, ps["lambda_2"], ps["lambda_1"]), np.where(swap, ps["lambda_1"], ps["lambda_2"])
    m1, m2 = np.where(swap, m2, m1), np.where(swap, m1, m2)
    a, b = tilde_coefficients(m1, m2)
    return dict(name=name, label=label, file=Path(path).name, m1=m1, m2=m2, l1=l1, l2=l2, lam_tilde=a * l1 + b * l2)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
_C = dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", grid="#e4e3df", joint="#0b0b0b",
          GW170817="#2a78d6", GW190425="#eb6834", pub="#7a5195", band="#2a78d6")


def _style(ax):
    ax.set_facecolor(_C["surface"])
    ax.grid(color=_C["grid"], lw=0.7, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(_C["grid"])
    ax.tick_params(colors=_C["ink2"], labelsize=9)


def plot(grid, posts: dict, joint, events: list[dict], out_png: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), dpi=150)
    fig.patch.set_facecolor(_C["surface"])
    ax = axes[0]
    for name, p in posts.items():
        ax.plot(grid, p, color=_C.get(name, "#1a9e77"), lw=1.6, ls="--", label=f"{name}")
    s = summary(grid, joint)
    ax.plot(grid, joint, color=_C["joint"], lw=2.2,
            label=f"joint: {s['median']:.0f} (90%: {s['low_90']:.0f}–{s['high_90']:.0f})")
    ax.axvspan(PUBLISHED["lo90"], PUBLISHED["hi90"], color=_C["pub"], alpha=0.10, lw=0, zorder=0)
    ax.axvline(PUBLISHED["median"], color=_C["pub"], lw=1.4, zorder=1,
               label=f"GW170817, LVK 2018: {PUBLISHED['median']:.0f} (90%: {PUBLISHED['lo90']:.0f}–{PUBLISHED['hi90']:.0f})")
    ax.set_xlim(0, 1500)
    ax.set_xlabel("Λ$_{1.4}$, tidal deformability of a 1.4 M$_\\odot$ star", color=_C["ink2"])
    ax.set_yticks([])
    top = ax.secondary_xaxis("top", functions=(radius_14, lambda r: 2.88e-6 * np.asarray(r) ** 7.5))
    top.set_xticks([9, 10, 11, 12, 13, 14])
    top.set_xlabel("R$_{1.4}$ (km, Annala et al. 2018)", color=_C["ink2"], fontsize=9)
    top.tick_params(colors=_C["ink2"], labelsize=8)
    ax.legend(frameon=False, fontsize=8, labelcolor=_C["ink2"], loc="upper right")
    ax = axes[1]
    m = np.linspace(1.0, 2.3, 120)
    cdf = np.cumsum(joint * np.gradient(grid)); cdf /= cdf[-1]
    draws = np.interp(np.random.default_rng(1).uniform(size=4000), cdf, grid)
    lam = lambda_of_m(m, draws)
    lo, med, hi = np.percentile(lam, [5, 50, 95], axis=0)
    ax.fill_between(m, lo, hi, color=_C["band"], alpha=0.25, lw=0, label="joint, 90%")
    ax.plot(m, med, color=_C["joint"], lw=2, label="joint, median")
    for e in events:
        for mm, ll in ((e["m1"], e["l1"]), (e["m2"], e["l2"])):
            ax.errorbar([np.median(mm)], [np.median(ll)], xerr=[[np.median(mm) - np.percentile(mm, 5)],
                        [np.percentile(mm, 95) - np.median(mm)]], fmt="o", ms=4, color=_C.get(e["name"]), capsize=2,
                        alpha=0.8)
        ax.plot([], [], "o", color=_C.get(e["name"]), label=f"{e['name']} stars (PE median Λ, uncorrelated)")
    ax.set_yscale("log")
    ax.set_ylim(5, 5000)
    ax.set_xlabel("Neutron-star mass (M$_\\odot$, source frame)", color=_C["ink2"])
    ax.set_ylabel("Tidal deformability Λ", color=_C["ink2"])
    ax.legend(frameon=False, fontsize=8, labelcolor=_C["ink2"], loc="upper right")
    for a in axes:
        _style(a)
    fig.suptitle("Joint neutron-star equation of state: Λ(m) = Λ$_{1.4}$ (m / 1.4 M$_\\odot$)$^{-6}$",
                 color=_C["ink"], fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=_C["surface"])
    plt.close(fig)
    return out_png


def run_ns_eos(
    events: Iterable[str] = ("GW170817", "GW190425"),
    spin_prior: str = "low",
    lambda_max: float = LAMBDA_MAX,
    cache_dir: str | Path = ".cache_gwosc",
    pe_cache: Optional[str | Path] = None,
    out_report_html: Optional[str | Path] = "ns_eos.html",
    out_summary_tsv: Optional[str | Path] = "ns_eos.tsv",
    plots_dir: str | Path = "ns_eos_plots",
) -> pd.DataFrame:
    """Lambda_1.4 and R_1.4 of each event and of all of them, under a common equation of state."""
    if spin_prior not in ("low", "high"):
        raise ValueError("spin_prior must be 'low' or 'high'")
    events = list(dict.fromkeys(events))
    unknown = [e for e in events if e not in EVENTS]
    if unknown:
        raise ValueError(f"no tidal analysis registered for {', '.join(unknown)}; supported: {', '.join(EVENTS)}")
    grid = np.linspace(0.0, GRID_MAX, 1201)
    data = [load_event(e, spin_prior, cache_dir, pe_cache) for e in events]
    lns, posts, rows = {}, {}, []
    for e in data:
        lns[e["name"]] = event_ln_likelihood(grid, e["m1"], e["m2"], e["lam_tilde"], lambda_max)
        posts[e["name"]] = normalize(grid, lns[e["name"]])
        s = summary(grid, posts[e["name"]])
        rows.append(dict(analysis=e["name"], label=e["label"], n_samples=len(e["m1"]),
                         m1_median=float(np.median(e["m1"])), m2_median=float(np.median(e["m2"])), **s))
        _log(f"{e['name']} ({e['label']}): Λ1.4 = {s['median']:.0f} [{s['low_90']:.0f}, {s['high_90']:.0f}], "
             f"R1.4 = {s['r_median']:.1f} km")
    joint = normalize(grid, np.sum(list(lns.values()), axis=0))
    s = summary(grid, joint)
    rows.append(dict(analysis="joint", label=f"{spin_prior} spin", n_samples=None, m1_median=None, m2_median=None, **s))
    _log(f"joint: Λ1.4 = {s['median']:.0f} [{s['low_90']:.0f}, {s['high_90']:.0f}], R1.4 = {s['r_median']:.1f} "
         f"[{s['r_low_90']:.1f}, {s['r_high_90']:.1f}] km")
    table = pd.DataFrame(rows)
    if out_summary_tsv:
        Path(out_summary_tsv).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_summary_tsv, sep="\t", index=False, float_format="%.4g")
        pd.DataFrame({"lambda_1.4": grid, **{f"p_{k}": v for k, v in posts.items()}, "p_joint": joint}).iloc[::4].to_csv(
            Path(out_summary_tsv).with_suffix(".posterior.tsv"), sep="\t", index=False, float_format="%.5g")
    if out_report_html:
        img = plot(grid, posts, joint, data, Path(plots_dir) / "ns_eos.png")
        t = table.set_index("analysis")
        per = "; ".join(f"{n}: {t.loc[n, 'median']:.0f} (90%: {t.loc[n, 'low_90']:.0f}–{t.loc[n, 'high_90']:.0f})"
                        for n in posts)
        paras = [
            f"<b>Λ<sub>1.4</sub> = {s['median']:.0f}</b> (90%: {s['low_90']:.0f}–{s['high_90']:.0f}), i.e. "
            f"<b>R<sub>1.4</sub> ≈ {s['r_median']:.1f} km</b> ({s['r_low_90']:.1f}–{s['r_high_90']:.1f}), from "
            f"{' and '.join(posts)} under a common equation of state ({spin_prior}-spin PE priors). Each event: {per}. "
            f"Published for GW170817 ({PUBLISHED['ref']}): {PUBLISHED['median']:.0f} (+{PUBLISHED['hi90'] - PUBLISHED['median']:.0f} / "
            f"−{PUBLISHED['median'] - PUBLISHED['lo90']:.0f}).",
            "All neutron stars follow one equation of state, so the four stars lie on one curve Λ(m) = Λ<sub>1.4</sub> "
            "(m / 1.4 M☉)<sup>−6</sup> (stars of nearly equal radius; De et al. 2018). Heavier stars are less "
            "deformable: the 1.5–1.9 M☉ stars of GW190425 have Λ ≈ 0.2–0.7 Λ<sub>1.4</sub>, which is why that event, "
            "louder at low frequency but seen essentially by one detector, adds little.",
            "Each event's likelihood is a kernel estimate on its samples of the measured combination Λ̃, evaluated at "
            "the value predicted at each sample's masses and divided by the prior on Λ̃ that the PE priors (Λ<sub>1</sub>, "
            f"Λ<sub>2</sub> uniform on [0, {lambda_max:.0f}]) imply at those masses. The mass population is the PE mass "
            "prior. Flat prior on Λ<sub>1.4</sub>.",
            "R<sub>1.4</sub> uses the empirical relation Λ<sub>1.4</sub> = 2.88 × 10⁻⁶ (R<sub>1.4</sub>/km)<sup>7.5</sup> "
            "of Annala et al. 2018 (scatter ~0.5 km across equations of state). The right panel shows Λ(m) with the "
            "stars' own PE medians, which the common-EOS fit pulls together.",
        ]
        Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
        write_simple_html_report(out_report_html, title="Joint neutron-star equation of state", paragraphs=paras,
                                 images=[img], tables=[("Λ1.4 and R1.4", table.to_html(
                                     index=False, na_rep="", float_format=lambda x: f"{x:.3g}"))])
        _log(f"report written to {out_report_html}")
    return table

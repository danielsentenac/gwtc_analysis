"""Hubble constant from a bright siren: the GW luminosity distance and the Hubble-flow velocity of the host.

For an event with an identified host galaxy, the GW signal gives the distance d and the host gives the
recession velocity v_r. Removing the host's peculiar velocity <v_p> leaves the Hubble-flow velocity
v_H = v_r - <v_p>, and at these distances (z ~ 0.01) the Hubble law v_H = H0 d is linear.

With sources uniform in volume and a detection limited by the GW signal, the selection term and the
volume prior on d cancel (Chen, Fishbach & Holz 2018, arXiv:1712.06531; Mandel, Farr & Gair 2019,
arXiv:1809.02063). With Gaussian velocity measurements and posterior samples d_i drawn with a d^2
distance prior, the H0 posterior is then, for a flat H0 prior,

    p(H0 | data) ∝ Σ_i N(v_H; H0 d_i, σ),   σ² = σ_r² + σ_p².

The samples must be conditioned on the sky position of the counterpart: the GWTC-1 GW170817 samples
are (the sky was fixed to AT2017gfo); for samples that are not, only those near the counterpart are
kept.

The default GW170817 setup follows LVK 2017 (arXiv:1710.05835), which gives
H0 = 70.0 (+12.0 / -8.0) km/s/Mpc (maximum a posteriori and 68% interval). The public GWTC-1 samples
used here are a later reanalysis, with a longer low-distance tail.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from .report import write_simple_html_report

C_KMS = 299792.458
H0_PRIOR = (10.0, 200.0)          # the prior of the spectral siren (hubble_constant), so they combine directly
MIN_SKY_SAMPLES = 200
_trapz = getattr(np, "trapezoid", None) or np.trapz


@dataclass(frozen=True)
class Counterpart:
    host: str
    transient: str
    ra_deg: float
    dec_deg: float
    v_recession: float          # km/s, CMB frame
    sigma_recession: float
    v_peculiar: float           # km/s
    sigma_peculiar: float
    ref: str
    published: dict             # maximum a posteriori and 68% interval


COUNTERPARTS = {
    "GW170817": Counterpart(
        host="NGC 4993", transient="AT2017gfo", ra_deg=197.450374, dec_deg=-23.381495,
        v_recession=3327.0, sigma_recession=72.0, v_peculiar=310.0, sigma_peculiar=150.0,
        ref="LVK 2017, arXiv:1710.05835",
        published=dict(map=70.0, lo68=62.0, hi68=82.0),
    ),
}


def _log(msg: str) -> None:
    print(f"[bright_siren] {msg}", flush=True)


def angular_separation_deg(ra: np.ndarray, dec: np.ndarray, ra0_deg: float, dec0_deg: float) -> np.ndarray:
    """Angle between (ra, dec) in radians and a position in degrees, in degrees."""
    ra0, dec0 = np.radians(ra0_deg), np.radians(dec0_deg)
    c = np.sin(dec) * np.sin(dec0) + np.cos(dec) * np.cos(dec0) * np.cos(ra - ra0)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def sky_conditioned_distances(samples: dict, cp: Counterpart, sky_radius_deg: float) -> tuple[np.ndarray, str]:
    """Luminosity distances of the samples at the counterpart's position, with a note on how they were chosen."""
    d = np.asarray(samples["luminosity_distance"], dtype=float)
    if "ra" not in samples or "dec" not in samples:
        return d, "no sky position in the samples: all of them used"
    sep = angular_separation_deg(np.asarray(samples["ra"], float), np.asarray(samples["dec"], float),
                                 cp.ra_deg, cp.dec_deg)
    if np.median(sep) < 0.1:
        return d, f"sky position fixed to {cp.transient} in the samples"
    keep = sep < sky_radius_deg
    if keep.sum() < MIN_SKY_SAMPLES:
        raise ValueError(f"only {keep.sum()} samples within {sky_radius_deg}° of {cp.transient}; "
                         "increase --sky-radius")
    return d[keep], f"{keep.sum()} of {len(d)} samples within {sky_radius_deg}° of {cp.transient}"


def h0_likelihood(h0: np.ndarray, distances: np.ndarray, v_hubble: float, sigma_v: float,
                  chunk: int = 2000) -> np.ndarray:
    """Σ_i N(v_H; H0 d_i, σ) on the H0 grid (unnormalized)."""
    out = np.zeros(len(h0))
    for i in range(0, len(distances), chunk):
        d = distances[i:i + chunk][None, :]
        out += np.exp(-0.5 * ((v_hubble - h0[:, None] * d) / sigma_v) ** 2).sum(axis=1)
    return out / len(distances)


def _normalize(h0: np.ndarray, p: np.ndarray) -> np.ndarray:
    return p / _trapz(p, h0)


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


def read_spectral_posterior(path: str | Path) -> np.ndarray:
    """H0 samples of a spectral-siren run: a posterior TSV, or a hubble_constant work directory."""
    path = Path(path).expanduser()
    if path.is_dir():
        for name in ("posterior_reweighted.tsv", "posterior.tsv"):
            if (path / name).exists():
                path = path / name
                break
        else:
            raise ValueError(f"no posterior_reweighted.tsv or posterior.tsv in {path}")
    df = pd.read_csv(path, sep="\t")
    if "H0" not in df:
        raise ValueError(f"{path} has no H0 column")
    return df["H0"].to_numpy(float)


def spectral_density(h0: np.ndarray, samples: np.ndarray) -> np.ndarray:
    """Density of spectral-siren H0 samples on the grid (Gaussian KDE, reflected at the prior bounds)."""
    from scipy.stats import gaussian_kde

    lo, hi = H0_PRIOR
    kde = gaussian_kde(samples)
    p = kde(h0) + kde(2 * lo - h0) + kde(2 * hi - h0)
    return np.where((h0 >= lo) & (h0 <= hi), p, 0.0)


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
    ax.set_title(title, color=ink, fontsize=10.5, loc="left", pad=34 if len(curves) < 3 else 46)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=ink2, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2,
              borderaxespad=0.3, handlelength=1.6)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def _read_pe_labels(pe_file: Path, labels: Optional[Iterable[str]]) -> dict[str, dict]:
    import h5py

    out = {}
    with h5py.File(pe_file, "r") as h:
        available = [k for k in h if isinstance(h[k], h5py.Group) and "posterior_samples" in h[k]]
        # LowSpin first: the usual primary analysis of binary neutron stars
        wanted = list(labels) if labels else sorted(available, key=lambda k: "lowspin" not in k.lower())
        missing = [k for k in wanted if k not in available]
        if missing:
            raise ValueError(f"label(s) {', '.join(missing)} not in {pe_file.name}; available: {', '.join(available)}")
        for k in wanted:
            ps = h[k]["posterior_samples"][()]
            out[k] = {c: ps[c] for c in ("luminosity_distance", "ra", "dec") if c in ps.dtype.names}
    return out


def _fmt(s: dict) -> str:
    return (f"{s['map']:.1f} (+{s['hpd68_high'] - s['map']:.1f} / −{s['map'] - s['hpd68_low']:.1f}) km/s/Mpc "
            f"(maximum a posteriori, 68%); median {s['median']:.1f}, 90%: {s['low_90']:.1f}–{s['high_90']:.1f}")


def run_bright_siren(
    src_name: str = "GW170817",
    pe_labels: Optional[Iterable[str]] = None,
    pe_file: Optional[str | Path] = None,
    cache_dir: str | Path = ".cache_gwosc",
    v_recession: Optional[tuple[float, float]] = None,
    v_peculiar: Optional[tuple[float, float]] = None,
    sky_radius_deg: float = 3.0,
    spectral_posterior: Optional[str | Path] = None,
    h0_range: tuple[float, float] = H0_PRIOR,
    out_report_html: Optional[str | Path] = "bright_siren.html",
    out_summary_tsv: Optional[str | Path] = "bright_siren.tsv",
    plots_dir: str | Path = "bright_siren_plots",
) -> pd.DataFrame:
    """Bright-siren H0 for each PE label of `src_name`, optionally combined with a spectral-siren posterior."""
    if src_name not in COUNTERPARTS:
        raise ValueError(f"no counterpart registered for {src_name}; supported: {', '.join(COUNTERPARTS)}")
    cp = COUNTERPARTS[src_name]
    vr, sr = v_recession or (cp.v_recession, cp.sigma_recession)
    vp, sp = v_peculiar or (cp.v_peculiar, cp.sigma_peculiar)
    v_h, sigma_v = vr - vp, float(np.hypot(sr, sp))
    default_velocities = v_recession is None and v_peculiar is None
    _log(f"{cp.host}: v_r = {vr:.0f} ± {sr:.0f}, <v_p> = {vp:.0f} ± {sp:.0f} → v_H = {v_h:.0f} ± {sigma_v:.0f} km/s")

    if pe_file is None:
        from .unofficial_pe import build_unofficial_pe_bundle

        pe_file = build_unofficial_pe_bundle(src_name, cache_dir=cache_dir, log_cb=_log)
        if pe_file is None:
            raise ValueError(f"could not build the PE bundle of {src_name}")
    pe_file = Path(pe_file)
    samples = _read_pe_labels(pe_file, pe_labels)

    lo, hi = h0_range
    h0 = np.linspace(lo, hi, int(round((hi - lo) / 0.05)) + 1)
    rows, posts, notes = [], {}, []
    for label, s in samples.items():
        d, how = sky_conditioned_distances(s, cp, sky_radius_deg)
        p = _normalize(h0, h0_likelihood(h0, d, v_h, sigma_v))
        posts[label] = p
        dq = np.percentile(d, [5, 50, 95])
        rows.append(dict(analysis=f"bright siren, {label}", n_samples=len(d), distance_median=dq[1],
                         distance_low_90=dq[0], distance_high_90=dq[2], **summarize(h0, p)))
        notes.append(f"{label}: {how}; d<sub>L</sub> = {dq[1]:.1f} Mpc (90%: {dq[0]:.1f}–{dq[2]:.1f}).")
        _log(f"{label}: H0 = {_fmt(rows[-1])}")

    combined = None
    if spectral_posterior:
        if (lo, hi) != H0_PRIOR:
            _log(f"WARN: H0 range {lo:g}–{hi:g} differs from the spectral-siren prior {H0_PRIOR[0]:g}–{H0_PRIOR[1]:g}")
        spec = spectral_density(h0, read_spectral_posterior(spectral_posterior))
        rows.append(dict(analysis="spectral siren", n_samples=None, **summarize(h0, spec)))
        main = next(iter(posts))
        combined = (main, _normalize(h0, posts[main] * spec), spec)
        rows.append(dict(analysis=f"bright ({main}) × spectral", n_samples=None, **summarize(h0, combined[1])))
        _log(f"combined with the spectral siren: H0 = {_fmt(rows[-1])}")

    table = pd.DataFrame(rows)
    if out_summary_tsv:
        Path(out_summary_tsv).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_summary_tsv, sep="\t", index=False, float_format="%.6g")
        grid = pd.DataFrame({"H0": h0, **{f"p_{k}": v for k, v in posts.items()}})
        if combined:
            grid["p_spectral"] = _normalize(h0, combined[2])
            grid["p_combined"] = combined[1]
        grid.iloc[::10].to_csv(Path(out_summary_tsv).with_suffix(".posterior.tsv"), sep="\t", index=False,
                               float_format="%.6g")

    if out_report_html:
        published = cp.published if default_velocities else None
        colors = ("#2a78d6", "#1a9e77", "#7a5195")
        curves = [(f"{k}: {r['map']:.0f}, 68% {r['hpd68_low']:.0f}–{r['hpd68_high']:.0f}", posts[k],
                   colors[i % len(colors)], "-" if i == 0 else "--") for i, (k, r) in enumerate(zip(posts, rows))]
        images = [_plot(h0, curves, published, cp.ref, f"Hubble constant from {src_name} and {cp.host} (bright siren)",
                        Path(plots_dir) / f"h0_bright_siren_{src_name}.png")]
        if combined:
            r_spec, r_comb = rows[-2], rows[-1]
            images.append(_plot(h0, [
                (f"bright siren: {rows[0]['map']:.0f}", posts[combined[0]], "#2a78d6", "--"),
                (f"spectral siren: {r_spec['map']:.0f}", combined[2], "#52514e", ":"),
                (f"combined: {r_comb['map']:.0f}, 68% {r_comb['hpd68_low']:.0f}–{r_comb['hpd68_high']:.0f}",
                 combined[1], "#eb6834", "-")], None, cp.ref,
                f"Bright siren ({src_name}) combined with the spectral siren", Path(plots_dir) / "h0_combined.png"))
        paras = [
            f"H<sub>0</sub> = <b>{_fmt(rows[0])}</b> from {src_name} ({next(iter(posts))}) and its host {cp.host}: "
            "the luminosity distance comes from the GW signal alone, the velocity from the redshift of the host.",
            f"Hubble-flow velocity v<sub>H</sub> = v<sub>r</sub> − ⟨v<sub>p</sub>⟩ = {vr:.0f} − {vp:.0f} = {v_h:.0f} ± "
            f"{sigma_v:.0f} km/s (recession velocity ± {sr:.0f}, peculiar velocity ± {sp:.0f}). Linear Hubble law, flat "
            f"H<sub>0</sub> prior {lo:g}–{hi:g} km/s/Mpc; for sources uniform in volume and a GW-limited detection, the "
            "selection term cancels the volume prior on the distance.",
            *notes,
        ]
        if published:
            paras.append(f"Published ({cp.ref}): {published['map']:.1f} (+{published['hi68'] - published['map']:.1f} / "
                         f"−{published['map'] - published['lo68']:.1f}) km/s/Mpc. It used the distance posterior of the "
                         "2017 analysis; the public GWTC-1 samples are a later reanalysis, whose longer low-distance "
                         "tail widens the upper side of the H<sub>0</sub> interval.")
        paras.append("The distance is degenerate with the inclination of the orbit, which dominates the uncertainty: "
                     "the low-distance tail corresponds to inclined orbits (smaller amplitude at given distance).")
        if combined:
            paras.append(f"Combined with the spectral siren ({Path(spectral_posterior).name}, same flat prior): "
                         f"H<sub>0</sub> = <b>{_fmt(rows[-1])}</b>. The two measurements are independent (different "
                         "events), so their posteriors multiply.")
        Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
        write_simple_html_report(out_report_html, title=f"Hubble constant (bright siren, {src_name})", paragraphs=paras,
                                 images=images,
                                 tables=[("H0 posterior", table.to_html(index=False, float_format=lambda x: f"{x:.3g}",
                                                                        na_rep=""))])
        _log(f"report written to {out_report_html}")
    return table

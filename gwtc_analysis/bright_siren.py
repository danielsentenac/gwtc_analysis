"""Hubble constant from bright sirens: the GW luminosity distance and the redshift of an identified host.

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
- ``auto`` (default): ``euclidean`` below z = 0.05, ``injections`` above.

The samples must be conditioned on the sky position of the counterpart: the GWTC-1 GW170817 samples are
(the sky was fixed to AT2017gfo); for samples that are not, only those near the counterpart are kept.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Optional

import numpy as np
import pandas as pd

from .report import write_simple_html_report

C_KMS = 299792.458
OM0 = 0.3065                      # as the spectral siren (hubble_constant)
H0_PRIOR = (10.0, 200.0)          # the prior of the spectral siren, so that they combine directly
EUCLIDEAN_ZMAX = 0.05             # selection "auto": euclidean below, injections above
MIN_SKY_SAMPLES = 200
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
# counterparts
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Counterpart:
    host: str
    transient: str
    ra_deg: float
    dec_deg: float
    run: str                                  # observing run: injections and PE prior default
    population: str                           # "BNS" or "BBH": mass model of the selection term
    ref: str
    z: Optional[float] = None                 # Hubble-flow redshift, for a distant host
    sigma_z: Optional[float] = None
    v_recession: Optional[float] = None       # km/s, CMB frame, for a nearby host
    sigma_recession: Optional[float] = None
    v_peculiar: Optional[float] = None        # km/s
    sigma_peculiar: Optional[float] = None
    pe_event: Optional[str] = None            # full name of the event's Zenodo PE file (None: unofficial bundle)
    pe_labels: tuple = ()                     # default PE labels (empty: all of them)
    published: Optional[dict] = None          # maximum a posteriori and 68% interval
    candidate: bool = False                   # association not established
    caveat: str = ""


COUNTERPARTS = {
    "GW170817": Counterpart(
        host="NGC 4993", transient="AT2017gfo", ra_deg=197.450374, dec_deg=-23.381495, run="O2",
        population="BNS", ref="LVK 2017, arXiv:1710.05835",
        v_recession=3327.0, sigma_recession=72.0, v_peculiar=310.0, sigma_peculiar=150.0,
        published=dict(map=70.0, lo68=62.0, hi68=82.0),
    ),
    "GW190521": Counterpart(
        host="AGN J124942.3+344929", transient="ZTF19abanrhr", ra_deg=192.42625, dec_deg=34.82472, run="O3a",
        population="BBH", ref="Graham et al. 2020, arXiv:2006.14122", z=0.438, sigma_z=0.0015,
        pe_event="GW190521_030229", pe_labels=("C01:IMRPhenomXPHM",), candidate=True,
        caveat="The association of ZTF19abanrhr with GW190521 is not established: the odds of a common source "
               "are 1 to 12 depending on the waveform model (Ashton et al. 2021, arXiv:2009.12346). The result "
               "is H<sub>0</sub> <i>if</i> the flare is the counterpart.",
    ),
}


def hubble_flow_redshift(cp: Counterpart, v_recession=None, v_peculiar=None, redshift=None) -> tuple[float, float, str]:
    """(z, σ_z, description) of the host, from overrides or the registry."""
    if redshift is not None:
        z, s = redshift
        return float(z), float(s), f"z = {z:g} ± {s:g} (given)"
    if v_recession is not None or v_peculiar is not None or cp.v_recession is not None:
        vr, sr = v_recession or (cp.v_recession, cp.sigma_recession)
        vp, sp = v_peculiar or (cp.v_peculiar or 0.0, cp.sigma_peculiar or 0.0)
        vh, sv = vr - vp, float(np.hypot(sr, sp))
        return vh / C_KMS, sv / C_KMS, (f"v<sub>H</sub> = v<sub>r</sub> − ⟨v<sub>p</sub>⟩ = {vr:.0f} − {vp:.0f} = "
                                        f"{vh:.0f} ± {sv:.0f} km/s (recession ± {sr:.0f}, peculiar ± {sp:.0f}), "
                                        f"z = {vh / C_KMS:.5f}")
    return float(cp.z), float(cp.sigma_z), f"z = {cp.z:g} ± {cp.sigma_z:g}"


def _log(msg: str) -> None:
    print(f"[bright_siren] {msg}", flush=True)


# ---------------------------------------------------------------------------
# likelihood
# ---------------------------------------------------------------------------
def angular_separation_deg(ra: np.ndarray, dec: np.ndarray, ra0_deg: float, dec0_deg: float) -> np.ndarray:
    """Angle between (ra, dec) in radians and a position in degrees, in degrees."""
    ra0, dec0 = np.radians(ra0_deg), np.radians(dec0_deg)
    c = np.sin(dec) * np.sin(dec0) + np.cos(dec) * np.cos(dec0) * np.cos(ra - ra0)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def sky_conditioned(samples: dict, cp: Counterpart, sky_radius_deg: float) -> tuple[np.ndarray, str]:
    """Mask of the samples at the counterpart's position, with a note on how they were chosen."""
    n = len(samples["luminosity_distance"])
    if "ra" not in samples or "dec" not in samples:
        return np.ones(n, bool), "no sky position in the samples: all of them used"
    sep = angular_separation_deg(np.asarray(samples["ra"], float), np.asarray(samples["dec"], float),
                                 cp.ra_deg, cp.dec_deg)
    if np.median(sep) < 0.1:
        return np.ones(n, bool), f"sky position fixed to {cp.transient} in the samples"
    keep = sep < sky_radius_deg
    if keep.sum() < MIN_SKY_SAMPLES:
        raise ValueError(f"only {keep.sum()} samples within {sky_radius_deg}° of {cp.transient}; "
                         "increase --sky-radius")
    return keep, f"{keep.sum()} of {n} samples within {sky_radius_deg}° of {cp.transient}"


def _logsumexp(a: np.ndarray, axis: int) -> np.ndarray:
    m = np.max(a, axis=axis, keepdims=True)
    m = np.where(np.isfinite(m), m, 0.0)
    with np.errstate(divide="ignore"):
        return np.squeeze(m, axis) + np.log(np.sum(np.exp(a - m), axis=axis))


def event_ln_likelihood(h0: np.ndarray, distances: np.ndarray, pe_prior: np.ndarray, z_obs: float, sigma_z: float,
                        rows: int = 64) -> np.ndarray:
    """ln ∫ N(z_obs; z, σ_z) p_pop(z) L_GW(d_L(z, H0)) dz on the H0 grid, up to a constant (no selection)."""
    d = np.asarray(distances, float)[None, :]
    ln_w0 = -np.log(np.asarray(pe_prior, float))[None, :]
    n = d.shape[1]
    out = np.empty(len(h0))
    for i in range(0, len(h0), rows):
        h = np.asarray(h0[i:i + rows], float)[:, None]
        z, pz, dd = _at_dimensionless_distance(d * h / C_KMS)
        with np.errstate(divide="ignore"):
            ln_w = ln_w0 + np.log(pz) - np.log(dd) + np.log(h / C_KMS)
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


def ln_selection_injections(h0: np.ndarray, inj: dict, ln_pop_mass: Callable, n_grid: int = 77
                            ) -> tuple[np.ndarray, np.ndarray]:
    """(ln β, effective injections) on the H0 grid, from found injections in the detector frame.

    `inj` has detector-frame mass_1, mass_2, luminosity_distance, their draw density `prior` in these
    variables, and `ntotal` (as `hubble_constant.detector_frame_injections` returns). At each trial H0 the
    injections are carried to the source frame and weighted by the population density in the detector frame,
    p_m(m1_s, m2_s) p_pop(z) / [(1+z)^2 dd_L/dz].
    """
    hg = np.linspace(h0[0], h0[-1], n_grid)
    ln_b, neff = np.empty(n_grid), np.empty(n_grid)
    d, m1, m2, prior = (np.asarray(inj[k], float) for k in ("luminosity_distance", "mass_1", "mass_2", "prior"))
    for i, h in enumerate(hg):
        z, pz, dd = _at_dimensionless_distance(d * h / C_KMS)
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            lnp = (ln_pop_mass(m1 / (1 + z), m2 / (1 + z)) + np.log(pz) - 2 * np.log1p(z)
                   - np.log(dd) + np.log(h / C_KMS))
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
# summaries and combination
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


# ---------------------------------------------------------------------------
# PE samples
# ---------------------------------------------------------------------------
def _event_pe_file(event: str, cache: Path) -> Path:
    """The event's Zenodo PE file, downloaded once into `cache`/files."""
    from . import hubble_constant as hc
    from . import parameters_estimation as pe

    files = cache / "files"
    local = sorted(files.glob(f"*{event}*PEDataRelease*.h*5"))
    if local:
        return Path(hc._pick_pe_file([dict(filename=p.name, path=p) for p in local])["path"])
    index = pe.build_zenodo_pe_index(cache_dir=str(cache / "index"), force_refresh=False)
    cands = index.get(event)
    if not cands:
        raise ValueError(f"no Zenodo PE release found for {event}")
    entry = hc._pick_pe_file(cands)
    files.mkdir(parents=True, exist_ok=True)
    dest = files / entry["filename"]
    tmp = dest.with_suffix(dest.suffix + ".part")
    _log(f"downloading the PE file of {event} to {files}")
    pe._download_http_with_progress(entry["url"], tmp, desc=event)
    tmp.replace(dest)
    return dest


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
            s = {c: ps[c] for c in ("luminosity_distance", "ra", "dec", "theta_jn", "iota") if c in ps.dtype.names}
            desc = ""
            try:
                v = h[k]["priors"]["analytic"]["luminosity_distance"][()]
                v = v[0] if hasattr(v, "__len__") and not isinstance(v, (bytes, str)) and len(v) else v
                desc = v.decode() if isinstance(v, bytes) else str(v)
            except (KeyError, TypeError, ValueError):
                pass
            s["prior_desc"] = desc
            out[k] = s
    return out


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
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


def viewing_angle(samples: dict) -> Optional[np.ndarray]:
    """Angle between the line of sight and the orbital (or total) angular momentum, folded to 0-90 degrees."""
    for k in ("theta_jn", "iota"):
        if k in samples:
            th = np.degrees(np.asarray(samples[k], float))
            return np.minimum(th, 180 - th)
    return None


def implied_h0(distances: np.ndarray, z_obs: float) -> np.ndarray:
    """H0 that puts each distance sample at the host redshift: d_L(z_obs; H0) = d."""
    return C_KMS / np.asarray(distances, float) * float(np.interp(z_obs, _ZG, _DL))


def degeneracy_table(distances: np.ndarray, view: np.ndarray, z_obs: float) -> pd.DataFrame:
    h = implied_h0(distances, z_obs)
    rows = []
    for lo, hi in ((0, 30), (30, 60), (60, 90)):
        m = (view >= lo) & (view <= hi)
        rows.append(dict(viewing_angle=f"{lo}–{hi}°", fraction=float(m.mean()),
                         distance_median=float(np.median(distances[m])) if m.any() else np.nan,
                         h0_median=float(np.median(h[m])) if m.any() else np.nan))
    return pd.DataFrame(rows)


def plot_degeneracy(distances: np.ndarray, view: np.ndarray, z_obs: float, title: str, out_png: Path) -> Path:
    """Distance against viewing angle, colored by the H0 each sample implies."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
    h = implied_h0(distances, z_obs)
    lo, hi = np.percentile(h, [2, 98])
    fig, ax = plt.subplots(figsize=(7.6, 4.6), dpi=150)
    fig.patch.set_facecolor(surface); ax.set_facecolor(surface)
    sc = ax.scatter(view, distances, c=np.clip(h, lo, hi), s=4, cmap="viridis", lw=0, alpha=0.7, zorder=2)
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("H$_0$ implied by the host redshift (km s$^{-1}$ Mpc$^{-1}$)", color=ink2, fontsize=9)
    cb.ax.tick_params(colors=ink2, labelsize=8)
    ax.set_xlim(0, 90)
    ax.set_xlabel("Viewing angle (deg): 0 = face-on, 90 = edge-on", color=ink2)
    ax.set_ylabel("Luminosity distance d$_L$ (Mpc)", color=ink2)
    ax.set_title(title, color=ink, fontsize=10.5, loc="left")
    ax.grid(color=grid, lw=0.7, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(grid)
    ax.tick_params(colors=ink2, labelsize=9)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def _fmt(s: dict) -> str:
    return (f"{s['map']:.1f} (+{s['hpd68_high'] - s['map']:.1f} / −{s['map'] - s['hpd68_low']:.1f}) km/s/Mpc "
            f"(maximum a posteriori, 68%); median {s['median']:.1f}, 90%: {s['low_90']:.1f}–{s['high_90']:.1f}")


# ---------------------------------------------------------------------------
# mode
# ---------------------------------------------------------------------------
def run_bright_siren(
    src_name: str = "GW170817",
    pe_labels: Optional[Iterable[str]] = None,
    pe_file: Optional[str | Path] = None,
    cache_dir: str | Path = ".cache_gwosc",
    pe_cache: Optional[str | Path] = None,
    v_recession: Optional[tuple[float, float]] = None,
    v_peculiar: Optional[tuple[float, float]] = None,
    redshift: Optional[tuple[float, float]] = None,
    sky_radius_deg: float = 3.0,
    selection: str = "auto",
    sensitivity_release: Optional[str] = None,
    sensitivity_file: Optional[str | Path] = None,
    far_threshold: float = 0.25,
    snr_threshold: float = 10.0,
    spectral_posterior: Optional[str | Path] = None,
    h0_range: tuple[float, float] = H0_PRIOR,
    out_report_html: Optional[str | Path] = "bright_siren.html",
    out_summary_tsv: Optional[str | Path] = "bright_siren.tsv",
    plots_dir: str | Path = "bright_siren_plots",
) -> pd.DataFrame:
    """Bright-siren H0 for each PE label of `src_name`, optionally combined with a spectral-siren posterior."""
    from . import hubble_constant as hc

    if src_name not in COUNTERPARTS:
        raise ValueError(f"no counterpart registered for {src_name}; supported: {', '.join(COUNTERPARTS)}")
    if selection not in ("auto", "euclidean", "injections"):
        raise ValueError(f"unknown selection {selection!r}; choose auto, euclidean or injections")
    cp = COUNTERPARTS[src_name]
    z_obs, sigma_z, z_desc = hubble_flow_redshift(cp, v_recession, v_peculiar, redshift)
    defaults = v_recession is None and v_peculiar is None and redshift is None
    _log(f"{cp.host}: {z_desc.replace('<sub>', '').replace('</sub>', '')}")
    if cp.candidate:
        _log(f"WARN: {cp.transient} is a candidate counterpart; the association is not established")

    if pe_file is None:
        if cp.pe_event:
            pe_file = _event_pe_file(cp.pe_event, Path(pe_cache).expanduser() if pe_cache else hc.default_pe_cache())
        else:
            from .unofficial_pe import build_unofficial_pe_bundle

            pe_file = build_unofficial_pe_bundle(src_name, cache_dir=cache_dir, log_cb=_log)
            if pe_file is None:
                raise ValueError(f"could not build the PE bundle of {src_name}")
    pe_file = Path(pe_file)
    samples = _read_pe_labels(pe_file, pe_labels or cp.pe_labels or None)

    lo, hi = h0_range
    h0 = np.linspace(lo, hi, int(round((hi - lo) / 0.05)) + 1)
    sel = selection if selection != "auto" else ("euclidean" if z_obs < EUCLIDEAN_ZMAX else "injections")
    if sel == "euclidean":
        ln_beta, sel_note = ln_selection_euclidean(h0), (
            "Selection: GW-limited detection of nearby sources, β(H<sub>0</sub>) ∝ H<sub>0</sub><sup>3</sup>, which "
            "cancels the volume factor of the population.")
    else:
        inj = load_injections(cp, sensitivity_release, sensitivity_file, far_threshold, snr_threshold)
        ln_beta, neff = ln_selection_injections(h0, inj, MASS_POPULATIONS[cp.population])
        sel_note = (f"Selection: {len(inj['prior'])} found injections of {cp.run} (FAR &lt; {far_threshold:g}/yr), "
                    f"reweighted at each H<sub>0</sub> to the {cp.population} population "
                    f"({'Power Law + Peak' if cp.population == 'BBH' else 'uniform 1–2.5 M☉'}, uniform in comoving "
                    f"volume); effective injections {neff.min():.0f}–{neff.max():.0f}.")
        if neff.min() < MIN_NEFF_INJECTIONS:
            sel_note += f" <b>Fewer than {MIN_NEFF_INJECTIONS} effective injections: β is noisy.</b>"
            _log(f"WARN: only {neff.min():.0f} effective injections")
    _log(f"selection: {sel}")

    rows, posts, notes, degen = [], {}, [], None
    for label, s in samples.items():
        keep, how = sky_conditioned(s, cp, sky_radius_deg)
        d = np.asarray(s["luminosity_distance"], float)[keep]
        prior, prior_desc = hc.pe_distance_prior(d, s["prior_desc"], cp.run)
        p = posterior_from_ln(h0, event_ln_likelihood(h0, d, prior, z_obs, sigma_z) - ln_beta)
        posts[label] = p
        if not degen and (view := viewing_angle(s)) is not None:
            degen = (label, d, view[keep])
        dq = np.percentile(d, [5, 50, 95])
        rows.append(dict(analysis=f"bright siren, {label}", n_samples=len(d), distance_median=dq[1],
                         distance_low_90=dq[0], distance_high_90=dq[2], **summarize(h0, p)))
        notes.append(f"{label}: {how}; d<sub>L</sub> = {dq[1]:.1f} Mpc (90%: {dq[0]:.1f}–{dq[2]:.1f}); PE distance "
                     f"prior {prior_desc[:80]}.")
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
        published = cp.published if defaults else None
        colors = ("#2a78d6", "#1a9e77", "#7a5195")
        curves = [(f"{k}: {r['map']:.0f}, 68% {r['hpd68_low']:.0f}–{r['hpd68_high']:.0f}", posts[k],
                   colors[i % len(colors)], "-" if i == 0 else "--") for i, (k, r) in enumerate(zip(posts, rows))]
        images = [_plot(h0, curves, published, cp.ref, f"Hubble constant from {src_name} and {cp.host} (bright siren)",
                        Path(plots_dir) / f"h0_bright_siren_{src_name}.png")]
        dtab = None
        if degen:
            dlab, dd, dview = degen
            images.append(plot_degeneracy(dd, dview, z_obs, f"{src_name} ({dlab}): the distance–inclination degeneracy",
                                          Path(plots_dir) / f"distance_inclination_{src_name}.png"))
            dtab = degeneracy_table(dd, dview, z_obs)
        if combined:
            r_spec, r_comb = rows[-2], rows[-1]
            images.append(_plot(h0, [
                (f"bright siren: {rows[0]['map']:.0f}", posts[combined[0]], "#2a78d6", "--"),
                (f"spectral siren: {r_spec['map']:.0f}", combined[2], "#52514e", ":"),
                (f"combined: {r_comb['map']:.0f}, 68% {r_comb['hpd68_low']:.0f}–{r_comb['hpd68_high']:.0f}",
                 combined[1], "#eb6834", "-")], None, cp.ref,
                f"Bright siren ({src_name}) combined with the spectral siren", Path(plots_dir) / "h0_combined.png"))
        paras = [
            f"H<sub>0</sub> = <b>{_fmt(rows[0])}</b> from {src_name} ({next(iter(posts))}) and its "
            f"{'candidate ' if cp.candidate else ''}host {cp.host} ({cp.transient}): the luminosity distance comes "
            "from the GW signal alone, the redshift from the host.",
        ]
        if cp.caveat:
            paras.append(f"<b>Caveat.</b> {cp.caveat}")
        paras += [
            f"Host: {z_desc}. Flat ΛCDM (Ω<sub>m</sub> = {OM0}), flat H<sub>0</sub> prior {lo:g}–{hi:g} km/s/Mpc, "
            "sources uniform in comoving volume. " + sel_note,
            *notes,
        ]
        if published:
            paras.append(f"Published ({cp.ref}): {published['map']:.1f} (+{published['hi68'] - published['map']:.1f} / "
                         f"−{published['map'] - published['lo68']:.1f}) km/s/Mpc. It used the distance posterior of the "
                         "2017 analysis and a linear Hubble law; the public GWTC-1 samples are a later reanalysis, "
                         "whose longer low-distance tail widens the upper side of the H<sub>0</sub> interval.")
        if dtab is not None:
            parts = "; ".join(f"{r.viewing_angle}: {100 * r.fraction:.0f}% of the samples, d<sub>L</sub> ≈ "
                              f"{r.distance_median:.3g} Mpc, H<sub>0</sub> ≈ {r.h0_median:.0f}"
                              for r in dtab.itertuples() if r.fraction > 0)
            paras.append("<b>The distance is degenerate with the inclination of the orbit</b>, which dominates the "
                         "uncertainty: an inclined binary is fainter than a face-on one at the same distance, so the "
                         "samples follow a band where the distance falls as the viewing angle grows (plot below the "
                         f"posterior). By viewing angle ({degen[0]}): {parts}. The upper tail of the H<sub>0</sub> "
                         "posterior comes from the inclined orbits; an independent constraint on the viewing angle "
                         "(e.g. from the jet of GW170817) narrows it.")
        else:
            paras.append("The distance is degenerate with the inclination of the orbit, which dominates the "
                         "uncertainty: the low-distance tail corresponds to inclined orbits (smaller amplitude at given "
                         "distance).")
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

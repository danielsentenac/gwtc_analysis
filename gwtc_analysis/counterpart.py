"""Counterpart of a GW event: where an electromagnetic counterpart and its host sit in the GW posterior.

For an event with a registered counterpart (or another position given by the user):

- association: the searched probability of the counterpart's position in the sky posterior of the PE samples;
- distance along its line of sight, against the distance of the host redshift for reference H0 values;
- viewing angle and the distance-inclination degeneracy;
- an independent constraint on the viewing angle (e.g. from the jet of GW170817), applied as weights on the samples.

The registry of counterparts (`COUNTERPARTS`) is shared with the bright-siren H0 (``hubble_constant --method bright``).
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable, Iterable, Optional

import numpy as np
import pandas as pd

from .report import write_simple_html_report


def default_style(fn: Callable) -> Callable:
    """Draw with matplotlib's default style, not the one pesummary sets on import (it thickens the legend lines)."""
    import functools

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        import matplotlib
        matplotlib.use("Agg")
        with matplotlib.rc_context(matplotlib.rcParamsDefault):
            return fn(*args, **kwargs)
    return wrapper

C_KMS = 299792.458
MIN_SKY_SAMPLES = 200


# ---------------------------------------------------------------------------
# registry
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
    print(f"[counterpart] {msg}", flush=True)


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
            s = {c: ps[c] for c in ("luminosity_distance", "ra", "dec", "theta_jn", "iota", "mass_1", "mass_2")
                 if c in ps.dtype.names}
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


def get_counterpart(event: str, ra_deg: Optional[float] = None, dec_deg: Optional[float] = None) -> Counterpart:
    """The registered counterpart of `event`, at another sky position when `ra_deg`, `dec_deg` are given."""
    if event not in COUNTERPARTS:
        raise ValueError(f"no counterpart registered for {event}; supported: {', '.join(COUNTERPARTS)}")
    if (ra_deg is None) != (dec_deg is None):
        raise ValueError("give both the right ascension and the declination of the position")
    cp = COUNTERPARTS[event]
    if ra_deg is not None:
        cp = replace(cp, ra_deg=float(ra_deg), dec_deg=float(dec_deg),
                     transient=f"the position RA {ra_deg:g}°, Dec {dec_deg:+g}°", published=None)
    return cp


def load_event_samples(event: str, cp: Counterpart, pe_file: Optional[str | Path] = None,
                       pe_labels: Optional[Iterable[str]] = None, cache_dir: str | Path = ".cache_gwosc",
                       pe_cache: Optional[str | Path] = None, log_cb: Callable[[str], None] = None
                       ) -> tuple[Path, dict[str, dict]]:
    """(PE file, samples by label) of the event: its Zenodo PE file, or the unofficial bundle for GW170817."""
    from . import hubble_constant as hc

    if pe_file is None:
        if cp.pe_event:
            pe_file = _event_pe_file(cp.pe_event, Path(pe_cache).expanduser() if pe_cache else hc.default_pe_cache())
        else:
            from .unofficial_pe import build_unofficial_pe_bundle

            pe_file = build_unofficial_pe_bundle(event, cache_dir=cache_dir, log_cb=log_cb or _log)
            if pe_file is None:
                raise ValueError(f"could not build the PE bundle of {event}")
    pe_file = Path(pe_file)
    return pe_file, _read_pe_labels(pe_file, pe_labels or cp.pe_labels or None)



# ---------------------------------------------------------------------------
# association
# ---------------------------------------------------------------------------
def _unit_vectors(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    ra, dec = np.asarray(ra, float), np.asarray(dec, float)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)


def sky_credible_level(ra: np.ndarray, dec: np.ndarray, ra0_deg: float, dec0_deg: float, max_samples: int = 20000,
                       seed: int = 0) -> float:
    """Searched probability of a position: the credible level of the smallest sky region of the samples (ra, dec in
    radians) that contains it. A Fisher-kernel density estimate on the sphere, its width the median distance to the
    n^(2/3)-th nearest sample (it shrinks as n^(-1/6), as Scott's rule in two dimensions); the level is the fraction
    of samples where the density is higher than at the position. It converges slowly for a sky with fine structure
    (GW190521 and ZTF19abanrhr: 0.55 from 20,000 samples, 0.62 from 50,000, against 0.64 from the LVK sky map),
    hence the sky map first when there is one (`sky_map_searched_probability`)."""
    x = _unit_vectors(ra, dec)
    if len(x) > max_samples:
        x = x[np.random.default_rng(seed).choice(len(x), max_samples, replace=False)]
    n = len(x)
    k = min(n - 1, max(10, int(n ** (2 / 3))))
    probe = x[np.random.default_rng(seed + 1).choice(n, min(n, 500), replace=False)]
    cos_k = np.sort(probe @ x.T, axis=1)[:, -k - 1]               # the k-th neighbour, the sample itself excluded
    sigma = max(float(np.median(np.arccos(np.clip(cos_k, -1, 1)))), np.radians(0.05))
    kappa = 1.0 / sigma ** 2

    def density(y: np.ndarray) -> np.ndarray:
        return np.concatenate([np.exp(kappa * (y[i:i + 500] @ x.T - 1)).sum(axis=1) for i in range(0, len(y), 500)])

    d0 = density(_unit_vectors(np.radians(ra0_deg), np.radians(dec0_deg))[None, :])[0]
    ds = density(x) - 1.0                                       # leave-one-out
    return float(np.mean(ds > d0))


def sky_map_searched_probability(path: str | Path, ra_deg: float, dec_deg: float) -> float:
    """Searched probability of a position in a FITS sky map (ligo.skymap crossmatch)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from ligo.skymap.io import read_sky_map
    from ligo.skymap.postprocess import crossmatch

    m = read_sky_map(str(path), moc=True)
    return float(crossmatch(m, SkyCoord(ra_deg * u.deg, dec_deg * u.deg)).searched_prob)


def _event_sky_map(cp: Counterpart, label: str, pe_name: str) -> Optional[Path]:
    """The LVK sky map of the event for this PE label (Zenodo archive of its catalog), or None."""
    from .skymap3d import fetch_fits

    try:
        return fetch_fits(cp.pe_event, label, pe_file_name=pe_name, log_cb=_log)
    except Exception as e:                     # no network, no archive: the samples give the level
        _log(f"no sky map for {cp.pe_event} ({e}); searched probability from the PE samples")
        return None


def weighted_quantiles(x: np.ndarray, q: Iterable[float], w: Optional[np.ndarray] = None) -> np.ndarray:
    x = np.asarray(x, float)
    w = np.ones_like(x) if w is None else np.asarray(w, float)
    o = np.argsort(x)
    c = np.cumsum(w[o])
    return np.interp(np.asarray(list(q)) * c[-1], c - 0.5 * w[o], x[o])


def viewing_angle_weights(view_deg: np.ndarray, constraint: tuple[float, float]) -> np.ndarray:
    """Weights of the samples for an independent Gaussian constraint (mean, sigma, degrees) on the viewing angle."""
    mu, s = constraint
    return np.exp(-0.5 * ((np.asarray(view_deg, float) - mu) / s) ** 2)


H0_REFERENCES = (("Planck", 67.4), ("SH0ES", 73.0))


@default_style
def _plot_distance(d: np.ndarray, w: Optional[np.ndarray], refs: list[tuple[str, float]], title: str,
                   out_png: Path) -> Path:
    """Distance posterior along the line of sight of the counterpart, with the distances of the host redshift."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
    fig, ax = plt.subplots(figsize=(7.6, 3.8), dpi=150)
    fig.patch.set_facecolor(surface); ax.set_facecolor(surface)
    bins = np.linspace(*np.percentile(d, [0.2, 99.8]), 60)
    ax.hist(d, bins=bins, density=True, histtype="step", lw=2, color="#2a78d6", label="GW only", zorder=3)
    if w is not None:
        ax.hist(d, bins=bins, weights=w, density=True, histtype="step", lw=2, color="#eb6834",
                label="with the viewing-angle constraint", zorder=4)
    for (name, dl), ls in zip(refs, ((0, (4, 2)), (0, (1, 1.5)))):
        ax.axvline(dl, color=ink2, lw=1.2, ls=ls, zorder=2, label=f"host redshift, {name} H$_0$: {dl:.0f} Mpc")
    ax.grid(False)
    ax.set_xlabel("Luminosity distance d$_L$ (Mpc)", color=ink2)
    ax.set_ylabel("Posterior density", color=ink2)
    ax.set_title(title, color=ink, fontsize=10.5, loc="left", pad=40)
    leg = ax.legend(frameon=False, fontsize=8.5, labelcolor=ink2, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2,
                    borderaxespad=0.1, handlelength=2.4)
    for h in leg.legend_handles:
        h.set_linewidth(1.5)               # gwpy, once imported, widens the legend lines
    ax.grid(axis="y", color=grid, lw=0.8, zorder=0)
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


def _q(x, w=None) -> tuple[float, float, float]:
    lo, med, hi = weighted_quantiles(x, (0.05, 0.5, 0.95), w)
    return float(med), float(lo), float(hi)


def viewing_angle(samples: dict) -> Optional[np.ndarray]:
    """Angle between the line of sight and the orbital (or total) angular momentum, folded to 0-90 degrees."""
    for k in ("theta_jn", "iota"):
        if k in samples:
            th = np.degrees(np.asarray(samples[k], float))
            return np.minimum(th, 180 - th)
    return None


def implied_h0(distances: np.ndarray, z_obs: float) -> np.ndarray:
    """H0 that puts each distance sample at the host redshift: d_L(z_obs; H0) = d."""
    from .h0_bright import luminosity_distance

    return float(luminosity_distance(z_obs, 1.0)) / np.asarray(distances, float)


def degeneracy_table(distances: np.ndarray, view: np.ndarray, z_obs: float) -> pd.DataFrame:
    h = implied_h0(distances, z_obs)
    rows = []
    for lo, hi in ((0, 30), (30, 60), (60, 90)):
        m = (view >= lo) & (view <= hi)
        rows.append(dict(viewing_angle=f"{lo}–{hi}°", fraction=float(m.mean()),
                         distance_median=float(np.median(distances[m])) if m.any() else np.nan,
                         h0_median=float(np.median(h[m])) if m.any() else np.nan))
    return pd.DataFrame(rows)


@default_style
def plot_degeneracy(distances: np.ndarray, view: np.ndarray, z_obs: float, title: str, out_png: Path,
                    band: Optional[tuple[float, float]] = None) -> Path:
    """Distance against viewing angle, colored by the H0 each sample implies; `band` (mean, sigma) shades a
    constraint on the viewing angle."""
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
    if band is not None:
        ax.axvspan(max(0, band[0] - band[1]), min(90, band[0] + band[1]), color="#eb6834", alpha=0.12, lw=0,
                   zorder=1)
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



# ---------------------------------------------------------------------------
# mode
# ---------------------------------------------------------------------------
def run_counterpart(
    event: str = "GW170817",
    ra_deg: Optional[float] = None,
    dec_deg: Optional[float] = None,
    pe_labels: Optional[Iterable[str]] = None,
    pe_file: Optional[str | Path] = None,
    cache_dir: str | Path = ".cache_gwosc",
    pe_cache: Optional[str | Path] = None,
    v_recession: Optional[tuple[float, float]] = None,
    v_peculiar: Optional[tuple[float, float]] = None,
    redshift: Optional[tuple[float, float]] = None,
    sky_radius_deg: float = 3.0,
    viewing_angle_constraint: Optional[tuple[float, float]] = None,
    pe_distance_prior: Optional[str] = None,
    sky_map: Optional[str | Path] = None,
    out_report_html: Optional[str | Path] = "counterpart.html",
    out_summary_tsv: Optional[str | Path] = "counterpart.tsv",
    plots_dir: str | Path = "counterpart_plots",
) -> pd.DataFrame:
    """Where the counterpart of `event` sits in its GW posterior: sky credible level, distance along its line of
    sight against the host redshift, viewing angle, with an optional EM constraint on the viewing angle.
    `sky_map`: FITS sky map for the searched probability; by default the event's LVK map when its PE file is read
    from Zenodo, else (or with "none") a kernel estimate on the PE samples."""
    from . import h0_bright as hb
    from .hubble_constant import pe_distance_prior as distance_prior

    cp = get_counterpart(event, ra_deg, dec_deg)
    if sky_radius_deg <= 0:
        raise ValueError("the sky radius must be > 0")
    if viewing_angle_constraint is not None and not (0 <= viewing_angle_constraint[0] <= 90
                                                     and viewing_angle_constraint[1] > 0):
        raise ValueError("the viewing-angle constraint is MEAN SIGMA in degrees, 0 <= MEAN <= 90, SIGMA > 0")
    z_obs, sigma_z, z_desc = hubble_flow_redshift(cp, v_recession, v_peculiar, redshift)
    refs = [(name, float(hb.luminosity_distance(z_obs, h))) for name, h in H0_REFERENCES]
    _log(f"{cp.host}: {z_desc.replace('<sub>', '').replace('</sub>', '')}")
    pe_path, samples = load_event_samples(event, cp, pe_file, pe_labels, cache_dir, pe_cache)

    rows, notes, images, first = [], [], [], None
    for label, s in samples.items():
        row = dict(label=label, n_samples=len(s["luminosity_distance"]))
        if "ra" in s and "dec" in s:
            sep = angular_separation_deg(np.asarray(s["ra"], float), np.asarray(s["dec"], float), cp.ra_deg, cp.dec_deg)
            fixed = float(np.median(sep)) < 0.1
            smap = None
            if not fixed and not (isinstance(sky_map, str) and sky_map.lower() == "none"):
                smap = Path(sky_map) if sky_map else (_event_sky_map(cp, label, pe_path.name)
                                                      if pe_file is None and cp.pe_event else None)
            if fixed:
                row["sky_credible_level"], row["sky_level_from"] = np.nan, "fixed in the samples"
            elif smap is not None:
                row["sky_credible_level"] = sky_map_searched_probability(smap, cp.ra_deg, cp.dec_deg)
                row["sky_level_from"] = f"sky map {smap.name}"
            else:
                row["sky_credible_level"] = sky_credible_level(s["ra"], s["dec"], cp.ra_deg, cp.dec_deg)
                row["sky_level_from"] = f"kernel estimate on {min(len(s['ra']), 20000)} PE samples"
            row["sky_separation_median_deg"] = float(np.median(sep))
        else:
            fixed = False
        try:
            keep, how = sky_conditioned(s, cp, sky_radius_deg)
        except ValueError as e:
            notes.append(f"{label}: {e}; the distance along the line of sight is not estimated.")
            rows.append(row)
            continue
        d = np.asarray(s["luminosity_distance"], float)[keep]
        prior, _ = distance_prior(d, s["prior_desc"], cp.run, kind=pe_distance_prior)
        row["n_line_of_sight"] = int(keep.sum())
        row["distance_median"], row["distance_low_90"], row["distance_high_90"] = _q(d)
        for name, dl in refs:
            row[f"distance_{name}"] = dl
            row[f"percentile_{name}"] = float(np.mean(d < dl))
        view = viewing_angle(s)
        w = None
        if view is not None:
            view = view[keep]
            row["viewing_angle_median"], row["viewing_angle_low_90"], row["viewing_angle_high_90"] = _q(view)
            if viewing_angle_constraint is not None:
                w = viewing_angle_weights(view, viewing_angle_constraint)
                row["ess_constrained"] = float(w.sum() ** 2 / (w * w).sum())
                (row["distance_constrained_median"], row["distance_constrained_low_90"],
                 row["distance_constrained_high_90"]) = _q(d, w)
                (row["viewing_angle_constrained_median"], row["viewing_angle_constrained_low_90"],
                 row["viewing_angle_constrained_high_90"]) = _q(view, w)
        elif viewing_angle_constraint is not None:
            notes.append(f"{label}: no viewing angle in the samples; the constraint is not applied.")
        rows.append(row)
        sky_txt = (f"the sky position is fixed to {cp.transient} in the samples" if fixed else
                   f"{cp.transient} is at the {100 * row['sky_credible_level']:.0f}% credible level of the sky "
                   f"posterior ({row['sky_level_from']})" if "sky_credible_level" in row
                   else "the samples have no sky position")
        notes.append(f"{label}: {sky_txt}; {how}. Along this line of sight d<sub>L</sub> = {row['distance_median']:.1f} "
                     f"Mpc (90%: {row['distance_low_90']:.1f}–{row['distance_high_90']:.1f}); the host redshift puts "
                     + " and ".join(f"it at {dl:.1f} Mpc for H<sub>0</sub> = {h:g} ({name}, percentile "
                                    f"{100 * row[f'percentile_{name}']:.0f})"
                                    for (name, dl), (_, h) in zip(refs, H0_REFERENCES)) + ".")
        if first is None:
            first = (label, d, view, w)

    table = pd.DataFrame(rows)
    if out_summary_tsv:
        Path(out_summary_tsv).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_summary_tsv, sep="\t", index=False, float_format="%.6g")
    if out_report_html:
        plots = Path(plots_dir)
        if first:
            label, d, view, w = first
            images.append(_plot_distance(d, w, refs, f"{event}: distance toward {cp.transient} "
                                         f"({label.split(':')[-1]})", plots / f"distance_{event}.png"))
            if view is not None:
                images.append(plot_degeneracy(d, view, z_obs, f"{event} ({label}): the distance–inclination degeneracy",
                                              plots / f"distance_inclination_{event}.png", band=viewing_angle_constraint))
                dtab = degeneracy_table(d, view, z_obs)
                parts = "; ".join(f"{r.viewing_angle}: {100 * r.fraction:.0f}% of the samples, d<sub>L</sub> ≈ "
                                  f"{r.distance_median:.3g} Mpc, H<sub>0</sub> ≈ {r.h0_median:.0f}"
                                  for r in dtab.itertuples() if r.fraction > 0)
                notes.append("<b>The distance is degenerate with the inclination of the orbit</b>: an inclined binary "
                             "is fainter than a face-on one at the same distance, so the samples follow a band where "
                             f"the distance falls as the viewing angle grows. By viewing angle ({label}): {parts} "
                             "(H<sub>0</sub> that each distance implies with the host redshift).")
        paras = [f"Counterpart of {event}: {cp.transient}, host {cp.host} ({cp.ref}). Host: {z_desc}. The GW "
                 f"posterior is read from {pe_path.name}; flat ΛCDM, Ω<sub>m</sub> = {hb.OM0}."]
        if cp.caveat:
            paras.append(f"<b>Caveat.</b> {cp.caveat}")
        paras += notes
        if viewing_angle_constraint is not None and first and first[3] is not None:
            r0 = next(r for r in rows if r["label"] == first[0])
            paras.append(f"With the viewing angle constrained to {viewing_angle_constraint[0]:g} ± "
                         f"{viewing_angle_constraint[1]:g}° (independent of the GW signal, e.g. from the jet): "
                         f"θ = {r0['viewing_angle_constrained_median']:.0f}° (90%: "
                         f"{r0['viewing_angle_constrained_low_90']:.0f}–{r0['viewing_angle_constrained_high_90']:.0f}), "
                         f"d<sub>L</sub> = {r0['distance_constrained_median']:.1f} Mpc (90%: "
                         f"{r0['distance_constrained_low_90']:.1f}–{r0['distance_constrained_high_90']:.1f}); "
                         f"{r0['ess_constrained']:.0f} effective samples. The H<sub>0</sub> posterior with this "
                         "constraint: <code>hubble_constant --method bright --viewing-angle MEAN SIGMA</code>.")
        Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
        write_simple_html_report(out_report_html, title=f"Counterpart of {event}", paragraphs=paras, images=images,
                                 tables=[("Counterpart in the GW posterior",
                                          table.to_html(index=False, float_format=lambda x: f"{x:.3g}", na_rep=""))])
        _log(f"report written to {out_report_html}")
    return table

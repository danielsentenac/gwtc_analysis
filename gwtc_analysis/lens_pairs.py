"""Strongly lensed pairs among the catalog events: the posterior-overlap screen.

A galaxy on the line of sight can split a signal into two (or four) images with the same detector-frame masses,
spins and sky position, arriving days to months apart with different amplitudes. This mode tests every pair of
the selected events for that consistency, as the first (tier-1) stage of the LVK lensing searches does with
Posterior Overlap (Haris et al. 2018, arXiv:1807.07062; GWTC-4.0 lensing, arXiv:2512.16347, Sec. 3.1.1 and
App. A.1). It ranks pairs; it does not confirm lensing, which needs joint parameter estimation of the two signals.

For a pair (1, 2) the overlap Bayes factor is

    B_overlap = ∫ p1(θ) p2(θ) / π(θ) dθ,

with p_i the PE posteriors and π the PE prior, over θ = (detector-frame m1, m2, spin magnitudes a1, a2, cosines
of the spin tilts, sky position). It is factorized into an intrinsic part (masses and spins, a Gaussian kernel
density in unbounded coordinates, evaluated by Monte Carlo on the samples of each event, both ways, averaged)
and a sky part (HEALPix histograms of the samples, smoothed). The PE priors are uniform in the component masses
(bilby UniformInComponents chirp mass and mass ratio), in the spin magnitudes, isotropic in the tilts and on
the sky; π here is that prior over a fixed common mass box, so log B_overlap carries one constant offset common
to every pair, irrelevant to the ranking.

The time delay enters through the ratio of its densities for lensed and unrelated pairs, R = p_L(Δt) / p_U(Δt),
in log Δt. p_U is the distribution of the delays of all the pairs of the selection (it follows the runs, their
gaps and the duty cycle). p_L comes from galaxy lenses: singular isothermal spheres with the SDSS early-type
velocity-dispersion function (Choi et al. 2007), weighted by their lensing cross-section along the line of sight
of each event (source redshift: the catalog median of each event), source positions uniform in the Einstein
disk, both images detectable with a magnification bias μ_-^{3/2} on the fainter one. Both densities keep 1% uniform
in log Δt (1.4 minutes to 27 years), which bounds R where both vanish. The ranking statistic is
B_overlap R, as in the LVK search.

The false-positive probability of a pair is the fraction of the pairs of the selection with a higher statistic:
the selection itself is the background, which holds as long as almost every pair is unrelated. The LVK uses
simulated unlensed events instead. A pair flagged at FPP < 1% is a candidate for joint analysis, not a detection.
"""
from __future__ import annotations

import itertools
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from .report import write_simple_html_report

from .cli_defaults import LENS_DEFAULT_CATALOGS as DEFAULT_CATALOGS  # noqa: E402
A_MAX = 0.99                          # spin-magnitude prior bound of the LVK PE
MASS_BOX = (1.0, 1000.0)             # detector-frame component masses of the reference prior (the PE constraint)
FLOOR_LOG10 = -30.0                   # log10 B_overlap of pairs whose chirp masses do not overlap (as the LVK)
FPP_THRESHOLD = 0.01
SKY_NSIDE = 64
SKY_FWHM_DEG = 2.0
_COLUMNS = ("mass_1", "mass_2", "a_1", "a_2", "cos_tilt_1", "cos_tilt_2", "ra", "dec", "geocent_time")
# velocity-dispersion function of early-type galaxies, SDSS (Choi, Park & Vogeley 2007, ApJ 658, 884)
VDF_SIGMA_STAR, VDF_ALPHA, VDF_BETA = 161.0, 2.32, 2.67
C_KMS = 299792.458
_DAY_S = 86400.0
DELAY_UNIFORM = 0.01                  # fraction of the delay densities uniform in log10 Δt over DELAY_RANGE_LOG10
DELAY_RANGE_LOG10 = (-3.0, 4.0)       # 1.4 minutes to 27 years


def _log(msg: str) -> None:
    print(f"[lens_pairs] {msg}", flush=True)


# ---------------------------------------------------------------------------
# events
# ---------------------------------------------------------------------------
def select_events(catalogs: Iterable[str], far_threshold: float = 1.0, min_mass: float = 3.0,
                  exclude: Iterable[str] = ()) -> pd.DataFrame:
    """BBH events of the catalogs with PE (GWOSC masses), FAR at most the threshold, both source masses >= min_mass.

    An event in several catalogs is kept once: from an update (GWTC-4.1) over the catalog it re-analyses, whatever
    the order on the command line."""
    from . import catalog_registry as reg
    from . import gw_stat as gw
    from .hubble_constant import _full_name

    exclude = set(exclude)
    best: dict[str, dict] = {}
    keys = reg.expand_all(catalogs)
    keys = [k for k in keys if not reg.is_update(k)] + [k for k in keys if reg.is_update(k)]   # updates last: they win
    for key in keys:
        cat = reg.CATALOGS[key]
        raw = gw.fetch_gwtc_events(cat.gwosc_list)
        for v in (raw.get("events") or {}).values():
            g, far = v.get("GPS"), v.get("far")
            m1, m2 = v.get("mass_1_source"), v.get("mass_2_source")
            if g is None or far is None or m1 is None or m2 is None or float(far) > far_threshold:
                continue
            if min(float(m1), float(m2)) < min_mass:
                continue
            name, common = _full_name(g), v.get("commonName") or ""
            if name in exclude or common in exclude:
                continue
            run = next((r.name for r in reg.OBSERVING_RUNS.values() if r.start_gps <= float(g) <= r.end_gps), "")
            best[name] = dict(event=name, common_name=common, catalog=key, run=run, gps=float(g),
                              far_per_yr=float(far), mass_1_source=float(m1), mass_2_source=float(m2),
                              redshift=float(v.get("redshift") or np.nan))
    cols = ["event", "common_name", "catalog", "run", "gps", "far_per_yr", "mass_1_source", "mass_2_source",
            "redshift"]
    return pd.DataFrame(list(best.values()), columns=cols).sort_values("gps").reset_index(drop=True)


# ---------------------------------------------------------------------------
# PE samples
# ---------------------------------------------------------------------------
def _extract(src: Path, out: Path) -> None:
    """Keep the overlap parameters of the preferred analysis (IMRPhenomXPHM, as hubble_constant) of a PE file."""
    import h5py

    from .hubble_constant import _choose_label

    with h5py.File(src, "r") as f:
        lab = _choose_label(f)
        if lab is None:
            raise ValueError(f"no IMRPhenomXPHM analysis in {src.name}")
        ps = f[lab]["posterior_samples"]
        names = ps.dtype.names or tuple(ps.keys())
        cols = {}
        for c in _COLUMNS:
            if c in names:
                cols[c] = np.asarray(ps[c], float)
            elif c.startswith("cos_tilt") and c.replace("cos_", "") in names:
                cols[c] = np.cos(np.asarray(ps[c.replace("cos_", "")], float))
            else:
                raise ValueError(f"{src.name} [{lab}]: no {c}")
        prior = ""
        try:
            v = f[lab]["priors"]["analytic"]["chirp_mass"][()]
            v = v[0] if hasattr(v, "__len__") and not isinstance(v, (bytes, str)) and len(v) else v
            prior = v.decode() if isinstance(v, bytes) else str(v)
        except (KeyError, TypeError, ValueError):
            pass
    part = out.with_suffix(".part")
    out.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(part, "w") as o:
        for c, v in cols.items():
            o.create_dataset(c, data=v, compression="gzip")
        o.attrs.update(label=lab, source=src.name, chirp_mass_prior=prior[:300])
    part.replace(out)


def _s3_index(cache: Path) -> dict[str, str]:
    """Event name -> PE object of the gwtc bucket (MinIO mirror of the Zenodo releases), cached for a day."""
    import time

    idx_file = cache / "index" / "s3_pe_objects.json"
    if idx_file.exists() and time.time() - idx_file.stat().st_mtime < 86400:
        return json.loads(idx_file.read_text())
    import re

    from minio import Minio

    client = Minio("minio-dev.odahub.fr", secure=True)
    out: dict[str, str] = {}
    for o in client.list_objects("gwtc", recursive=True):
        n = o.object_name
        m = re.search(r"(GW\d{6}_\d{6})", n)
        if m and "PEDataRelease" in n and n.endswith((".h5", ".hdf5")):
            prev = out.get(m.group(1))
            if prev is None or ("nocosmo" in n and "nocosmo" not in prev):
                out[m.group(1)] = n
    idx_file.parent.mkdir(parents=True, exist_ok=True)
    idx_file.write_text(json.dumps(out))
    return out


def _download(row, cache: Path, data_repo: str, s3_index: Optional[dict]) -> tuple[Path, bool]:
    """The event's PE file, and whether this call downloaded it (a file already in the shared cache belongs to
    whichever mode fetched it and is kept)."""
    from .hubble_constant import _full_name

    names = [row["event"]] + [_full_name(row["gps"] + d) for d in (-2, -1, 1, 2)]
    if data_repo == "s3":
        from minio import Minio

        from .parameters_estimation import _download_s3_with_progress

        obj = next((s3_index[n] for n in names if n in s3_index), None)
        if obj is None:
            raise ValueError(f"no PE file for {row['event']} on the S3 mirror")
        dest = cache / "files" / Path(obj).name
        if dest.exists() and dest.stat().st_size > 0:
            return dest, False
        tmp = dest.with_suffix(dest.suffix + ".part")
        _download_s3_with_progress(Minio("minio-dev.odahub.fr", secure=True), bucket="gwtc", object_name=obj,
                                   out_path=tmp, desc=row["event"])
        tmp.replace(dest)
        return dest, True
    from .spin_population import _pe_file

    before = set((cache / "files").glob("*")) if (cache / "files").exists() else set()
    src = _pe_file(row, cache)
    return src, src not in before


def fetch_samples(events: pd.DataFrame, cache: Path, data_repo: str = "zenodo", workers: int = 4,
                  keep_files: bool = False) -> dict[str, Path]:
    """Extract (once) the overlap parameters of each event; restartable."""
    out_dir = cache / "samples_lens"
    out = {r["event"]: out_dir / f"{r['event']}.h5" for _, r in events.iterrows()}
    todo = [r for _, r in events.iterrows() if not out[r["event"]].exists()]
    if not todo:
        return out
    s3_index = _s3_index(cache) if data_repo == "s3" else None
    if data_repo != "s3":       # the Zenodo PE index, built once here: the download threads only read it
        from .parameters_estimation import build_zenodo_pe_index

        build_zenodo_pe_index(cache_dir=str(cache / "index"), force_refresh=False)
    _log(f"extracting the PE samples of {len(todo)} events ({data_repo})")

    def one(r):
        src, downloaded = _download(r, cache, data_repo, s3_index)
        _extract(src, out[r["event"]])
        if downloaded and not keep_files:
            src.unlink(missing_ok=True)
        return r["event"]

    errors = []
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(one, r): r["event"] for r in todo}
        for i, fu in enumerate(as_completed(futs), 1):
            try:
                _log(f"[{i}/{len(todo)}] {fu.result()}: extracted")
            except Exception as e:
                errors.append(f"{futs[fu]}: {type(e).__name__}: {e}")
    if errors:
        raise ValueError("PE extraction failed (run again to resume): " + "; ".join(errors))
    return out


def load_samples(path: Path, n: int, rng) -> dict:
    import h5py

    with h5py.File(path, "r") as h:
        d = {c: h[c][:] for c in _COLUMNS}
        d["chirp_mass_prior"] = str(h.attrs.get("chirp_mass_prior", ""))
    idx = rng.permutation(len(d["mass_1"]))[:n]
    for c in _COLUMNS:
        d[c] = d[c][idx]
    return d


# ---------------------------------------------------------------------------
# posterior overlap
# ---------------------------------------------------------------------------
def _logit(x):
    return np.log(x) - np.log1p(-x)


def intrinsic_coords(d: dict) -> tuple[np.ndarray, np.ndarray]:
    """Unbounded coordinates y = (ln Mc, logit q, logit a1/a_max, logit a2/a_max, artanh cos t1, artanh cos t2) of
    the samples, and ln of the reference PE prior density in y (uniform in the component masses over MASS_BOX,
    uniform spin magnitudes, isotropic tilts)."""
    m1 = np.maximum(d["mass_1"], d["mass_2"])
    m2 = np.minimum(d["mass_1"], d["mass_2"])
    q = np.clip(m2 / m1, 1e-6, 1 - 1e-9)
    mc = (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2
    u1 = np.clip(d["a_1"] / A_MAX, 1e-9, 1 - 1e-9)
    u2 = np.clip(d["a_2"] / A_MAX, 1e-9, 1 - 1e-9)
    c1 = np.clip(d["cos_tilt_1"], -1 + 1e-9, 1 - 1e-9)
    c2 = np.clip(d["cos_tilt_2"], -1 + 1e-9, 1 - 1e-9)
    y = np.column_stack([np.log(mc), _logit(q), _logit(u1), _logit(u2), np.arctanh(c1), np.arctanh(c2)])
    lo, hi = MASS_BOX
    ln_pi_m1m2 = np.log(2.0 / (hi - lo) ** 2)                  # uniform on m2 <= m1 in the box
    # |d(m1, m2) / d(Mc, q)| = Mc (1 + q)^(2/5) q^(-6/5); then d Mc / d ln Mc = Mc, d q / d logit q = q (1 - q)
    ln_pi = (ln_pi_m1m2 + np.log(mc) + 0.4 * np.log1p(q) - 1.2 * np.log(q) + np.log(mc) + np.log(q * (1 - q))
             + np.log(u1 * (1 - u1)) + np.log(u2 * (1 - u2))           # a = a_max u, uniform on [0, a_max]
             + np.log((1 - c1 ** 2) / 2) + np.log((1 - c2 ** 2) / 2))  # cos tilt uniform on [-1, 1]
    return y, ln_pi


def _kde_logpdf(train: np.ndarray, points: np.ndarray) -> np.ndarray:
    from scipy.stats import gaussian_kde

    return gaussian_kde(train.T).logpdf(points.T)


def log_overlap_intrinsic(ya: np.ndarray, pia: np.ndarray, yb: np.ndarray, pib: np.ndarray) -> float:
    """ln ∫ pa pb / π over the intrinsic coordinates: the Monte Carlo average of p̂b / π on the samples of a and
    of p̂a / π on those of b, averaged in B."""
    from scipy.special import logsumexp

    ab = logsumexp(_kde_logpdf(yb, ya) - pia) - np.log(len(ya))
    ba = logsumexp(_kde_logpdf(ya, yb) - pib) - np.log(len(yb))
    return float(np.logaddexp(ab, ba) - np.log(2))


def sky_density(d: dict, nside: int = SKY_NSIDE, fwhm_deg: float = SKY_FWHM_DEG) -> np.ndarray:
    """Probability per HEALPix pixel of the sky samples, smoothed with a Gaussian beam, floored and normalized."""
    import healpy as hp

    pix = hp.ang2pix(nside, np.pi / 2 - d["dec"], d["ra"] % (2 * np.pi))
    m = np.bincount(pix, minlength=hp.nside2npix(nside)).astype(float)
    m = hp.smoothing(m / m.sum(), fwhm=np.radians(fwhm_deg))
    m = np.maximum(m, 1e-12)
    return m / m.sum()


def log_overlap_sky(pa: np.ndarray, pb: np.ndarray) -> float:
    """ln ∫ pa pb / π over the sky, π = 1 / 4π: 4π Σ Pa Pb / A_pix = N_pix Σ Pa Pb."""
    return float(np.log(len(pa) * np.sum(pa * pb)))


def _chirp_mass_det(d: dict) -> np.ndarray:
    return (d["mass_1"] * d["mass_2"]) ** 0.6 / (d["mass_1"] + d["mass_2"]) ** 0.2


# ---------------------------------------------------------------------------
# time delays
# ---------------------------------------------------------------------------
def sis_time_delay_days(sigma_kms, z_l, z_s, y, cosmo=None) -> np.ndarray:
    """Delay between the two images of a singular isothermal sphere, days:
    Δt = (1 + z_l) / c · D_l D_s / D_ls · 2 y θ_E², θ_E = 4π (σ / c)² D_ls / D_s."""
    from astropy.cosmology import Planck18

    cosmo = cosmo or Planck18
    z_l, z_s = np.asarray(z_l, float), np.asarray(z_s, float)
    dm_l = cosmo.comoving_distance(z_l).to_value("Mpc")
    dm_s = cosmo.comoving_distance(z_s).to_value("Mpc")
    d_l, d_s, d_ls = dm_l / (1 + z_l), dm_s / (1 + z_s), (dm_s - dm_l) / (1 + z_s)
    theta_e = 4 * np.pi * (np.asarray(sigma_kms) / C_KMS) ** 2 * d_ls / d_s
    mpc_km = 3.0856775814913673e19
    return (1 + z_l) / C_KMS * d_l * d_s / d_ls * mpc_km * 2 * np.asarray(y) * theta_e ** 2 / _DAY_S


def lensed_delays(z_sources: np.ndarray, n: int = 200_000, seed: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Galaxy-lens delays (days) and their weights, for sources at the given redshifts (equally likely)."""
    from astropy.cosmology import Planck18
    from scipy.special import gammaln

    rng = np.random.default_rng(seed)
    z_s = rng.choice(np.asarray(z_sources, float), n)
    z_l = z_s * rng.random(n)                                    # uniform in (0, z_s): weight z_s below
    # sigma from phi(sigma) sigma^4 (the SIS cross-section ∝ θ_E² ∝ σ⁴), by inverse CDF on a grid
    sg = np.linspace(20, 600, 4000)
    x = sg / VDF_SIGMA_STAR
    pdf = np.exp(VDF_ALPHA * np.log(x) - x ** VDF_BETA - gammaln(VDF_ALPHA / VDF_BETA)) / sg * sg ** 4
    cdf = np.cumsum(pdf); cdf /= cdf[-1]
    sigma = np.interp(rng.random(n), cdf, sg)
    y = np.sqrt(rng.random(n))                                   # uniform in the Einstein disk
    dm_l = Planck18.comoving_distance(z_l).to_value("Mpc")
    dm_s = Planck18.comoving_distance(z_s).to_value("Mpc")
    d_l, d_s, d_ls = dm_l / (1 + z_l), dm_s / (1 + z_s), (dm_s - dm_l) / (1 + z_s)
    # optical depth per z_l: (D_l D_ls / D_s)^2 (1 + z_l)^2 / E(z_l); detection of both images: mu_-^(3/2)
    mu_minus = (1 - y) / y
    w = z_s * (d_l * d_ls / d_s) ** 2 * (1 + z_l) ** 2 / Planck18.efunc(z_l) * mu_minus ** 1.5
    dt = sis_time_delay_days(sigma, z_l, z_s, y)
    ok = (dt > 0) & np.isfinite(w) & (w > 0)
    return dt[ok], w[ok]


def window_overlap_days(delays_days: np.ndarray, segments_gps: list[tuple[float, float]]) -> np.ndarray:
    """W(Δt) = ∫ I(t) I(t + Δt) dt in days, I the indicator of the observing segments: the time over which two
    signals Δt apart can both arrive inside the observations."""
    dt = np.asarray(delays_days, float) * _DAY_S
    out = np.zeros_like(dt)
    for a1, b1 in segments_gps:
        for a2, b2 in segments_gps:
            out += np.clip(np.minimum(b1, b2 - dt) - np.maximum(a1, a2 - dt), 0.0, None)
    return out / _DAY_S


class DelayRatio:
    """ln R(Δt) = ln p_L(log10 Δt) - ln p_U(log10 Δt), both Gaussian kernel densities.

    With the observing `segments` (GPS), the lensed delays are weighted by W(Δt), the time over which both images
    can arrive inside the observations: a galaxy lens with a delay longer than the runs gives no observed pair. The
    unlensed density, from the delays of the pairs themselves, already includes it."""

    def __init__(self, pair_delays_days: np.ndarray, z_sources: np.ndarray, seed: int = 1,
                 segments: Optional[list[tuple[float, float]]] = None):
        from scipy.stats import gaussian_kde

        u = np.log10(np.maximum(pair_delays_days, 1e-3))
        self.unlensed = gaussian_kde(u)
        dt, w = lensed_delays(z_sources, seed=seed)
        if segments:
            w = w * window_overlap_days(dt, segments)
            ok = w > 0
            dt, w = dt[ok], w[ok]
        rng = np.random.default_rng(seed)
        keep = rng.choice(len(dt), 20_000, p=w / w.sum())         # resampled: a kde of 20 000 points is fast
        self.lensed_log10 = np.log10(dt[keep])
        self.lensed = gaussian_kde(self.lensed_log10)
        self.median_lensed_days = float(10 ** np.median(self.lensed_log10))

    def _ln_mix(self, kde, u):
        """A kernel density with a small part uniform in log10 Δt: bounds R where both densities vanish (minutes,
        years), which a ratio of two kernel tails does not."""
        lo, hi = DELAY_RANGE_LOG10
        uni = np.where((u >= lo) & (u <= hi), 1.0 / (hi - lo), 1e-300)
        return np.logaddexp(np.log1p(-DELAY_UNIFORM) + kde.logpdf(u), np.log(DELAY_UNIFORM * uni))

    def ln_p_lensed(self, u):
        return self._ln_mix(self.lensed, u)

    def ln_p_unlensed(self, u):
        return self._ln_mix(self.unlensed, u)

    def ln_r(self, delays_days: np.ndarray) -> np.ndarray:
        u = np.log10(np.maximum(np.asarray(delays_days, float), 1e-3))
        return self.ln_p_lensed(u) - self.ln_p_unlensed(u)


# ---------------------------------------------------------------------------
# the search
# ---------------------------------------------------------------------------
def search_pairs(samples: dict[str, dict], events: pd.DataFrame, seed: int = 1,
                 progress_every: int = 200) -> tuple[pd.DataFrame, DelayRatio]:
    """All pairs: ln B_intrinsic, ln B_sky, log10 B_overlap, log10 R, log10 statistic, FPP."""
    names = list(events["event"])
    t = {n: float(np.median(samples[n]["geocent_time"])) for n in names}
    coords = {n: intrinsic_coords(samples[n]) for n in names}
    sky = {n: sky_density(samples[n]) for n in names}
    lnmc = {n: np.quantile(np.log(_chirp_mass_det(samples[n])), [0.001, 0.999]) for n in names}
    pairs = list(itertools.combinations(names, 2))
    delays = np.array([abs(t[b] - t[a]) / _DAY_S for a, b in pairs])
    z = events["redshift"].to_numpy(float)
    segments = None
    if "run" in events:
        from .catalog_registry import OBSERVING_RUNS

        runs = sorted({r for r in events["run"] if r in OBSERVING_RUNS})
        segments = [(OBSERVING_RUNS[r].start_gps, OBSERVING_RUNS[r].end_gps) for r in runs] or None
    ratio = DelayRatio(delays, z[np.isfinite(z)], seed=seed, segments=segments)
    ln_r = ratio.ln_r(delays)
    rows = []
    margin = np.log(1.05)
    for i, ((a, b), dt, lr) in enumerate(zip(pairs, delays, ln_r)):
        overlap = lnmc[a][0] - margin <= lnmc[b][1] and lnmc[b][0] - margin <= lnmc[a][1]
        if overlap:
            lbi = log_overlap_intrinsic(*coords[a], *coords[b])
            lbs = log_overlap_sky(sky[a], sky[b])
            l10 = (lbi + lbs) / np.log(10)
        else:
            lbi, lbs, l10 = np.nan, np.nan, FLOOR_LOG10
        l10 = max(l10, FLOOR_LOG10)
        rows.append(dict(event_1=a, event_2=b, delay_days=dt, ln_B_intrinsic=lbi, ln_B_sky=lbs,
                         log10_B_overlap=l10, log10_R=lr / np.log(10), log10_stat=l10 + lr / np.log(10)))
        if progress_every and (i + 1) % progress_every == 0:
            _log(f"{i + 1}/{len(pairs)} pairs")
    df = pd.DataFrame(rows)
    s = df["log10_stat"].to_numpy()
    order = np.sort(s)
    df["fpp"] = (len(s) - np.searchsorted(order, s, side="left")) / len(s)
    return df.sort_values("log10_stat", ascending=False).reset_index(drop=True), ratio


def _plots(df: pd.DataFrame, ratio: DelayRatio, plots_dir: Path, threshold: float) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plots_dir.mkdir(parents=True, exist_ok=True)
    out = []
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    live = df[df["log10_B_overlap"] > FLOOR_LOG10]
    ax[0].hist(live["log10_stat"], bins=50, color="#2a78d6", alpha=0.85)
    cut = df.loc[df["fpp"] <= threshold, "log10_stat"]
    if len(cut):
        ax[0].axvline(cut.min(), color="#eb6834", ls="--", label=f"FPP ≤ {threshold:g}: {len(cut)} pairs")
        ax[0].legend(frameon=False)
    ax[0].set_xlabel("log10 (B_overlap R)")
    ax[0].set_ylabel("pairs")
    ax[0].set_title(f"{len(live)} pairs with overlapping chirp masses (of {len(df)})")
    u = np.linspace(-2, 3.5, 300)
    ax[1].plot(u, np.exp(ratio.ln_p_lensed(u)), color="#2a78d6", label="galaxy-lensed pairs (SIS model)")
    ax[1].plot(u, np.exp(ratio.ln_p_unlensed(u)), color="#8a8984", label="all pairs of the selection")
    ax2 = ax[1].twinx()
    ax2.plot(u, (ratio.ln_p_lensed(u) - ratio.ln_p_unlensed(u)) / np.log(10), color="#eb6834", ls="--",
             label="log10 R")
    ax2.set_ylabel("log10 R", color="#eb6834")
    ax[1].set_xlabel("log10 (time delay / days)")
    ax[1].set_ylabel("density")
    ax[1].legend(frameon=False, fontsize=8, loc="upper left")
    ax[1].set_title("Time-delay densities")
    fig.tight_layout()
    p = plots_dir / "lens_pairs.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    out.append(p)
    return out


def run_lens_pairs(
    catalogs: Iterable[str] = DEFAULT_CATALOGS,
    far_threshold: float = 1.0,
    min_mass: float = 3.0,
    exclude: Iterable[str] = (),
    pe_cache: Optional[str | Path] = None,
    pe_samples: int = 2000,
    data_repo: str = "zenodo",
    workers: int = 4,
    keep_pe_files: bool = False,
    fpp_threshold: float = FPP_THRESHOLD,
    seed: int = 1,
    out_report_html: Optional[str | Path] = "lens_pairs.html",
    out_summary_tsv: Optional[str | Path] = "lens_pairs.tsv",
    plots_dir: str | Path = "lens_pairs_plots",
) -> pd.DataFrame:
    from .catalog_registry import catalog_name, expand_all
    from .hubble_constant import default_pe_cache

    cats = expand_all(catalogs)
    events = select_events(cats, far_threshold, min_mass, exclude)
    if len(events) < 2:
        raise ValueError(f"{len(events)} event(s) selected: nothing to pair")
    npairs = len(events) * (len(events) - 1) // 2
    _log(f"{len(events)} BBH events ({', '.join(catalog_name(c) for c in cats)}, FAR <= {far_threshold:g}/yr, "
         f"masses >= {min_mass:g} M_sun): {npairs} pairs")
    cache = Path(pe_cache).expanduser() if pe_cache else default_pe_cache()
    paths = fetch_samples(events, cache, data_repo=data_repo, workers=workers, keep_files=keep_pe_files)
    rng = np.random.default_rng(seed)
    samples = {n: load_samples(paths[n], pe_samples, rng) for n in events["event"]}
    df, ratio = search_pairs(samples, events, seed=seed)
    common = dict(zip(events["event"], events["common_name"]))
    flagged = df[df["fpp"] <= fpp_threshold]
    floor = int((df["log10_B_overlap"] <= FLOOR_LOG10).sum())
    _log(f"{len(df)} pairs, {floor} without chirp-mass overlap; {len(flagged)} with FPP <= {fpp_threshold:g}; "
         f"top: {df.iloc[0]['event_1']} + {df.iloc[0]['event_2']}, log10 stat {df.iloc[0]['log10_stat']:.2f}")
    if out_summary_tsv:
        Path(out_summary_tsv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_summary_tsv, sep="\t", index=False, float_format="%.5g")
        events.to_csv(Path(out_summary_tsv).with_suffix(".events.tsv"), sep="\t", index=False)
    if out_report_html:
        imgs = _plots(df, ratio, Path(plots_dir), fpp_threshold)
        top = df.head(max(20, len(flagged))).copy()
        top.insert(2, "names", [f"{common.get(a, a)} + {common.get(b, b)}" for a, b in zip(top.event_1, top.event_2)])
        paras = [
            f"<b>Lensed-pair screen</b> of {len(events)} binary black holes "
            f"({', '.join(catalog_name(c) for c in cats)}; FAR ≤ {far_threshold:g}/yr, both source masses ≥ "
            f"{min_mass:g} M<sub>☉</sub>): {len(df)} pairs, {floor} of them with chirp masses that cannot "
            f"overlap (log<sub>10</sub> B<sub>overlap</sub> set to {FLOOR_LOG10:g}). {len(flagged)} pairs have a "
            f"false-positive probability ≤ {fpp_threshold:g}.",
            "Statistic: the posterior-overlap Bayes factor B<sub>overlap</sub> = ∫ p<sub>1</sub> p<sub>2</sub> / π "
            "over the detector-frame masses, spin magnitudes and tilts (kernel density) and the sky position "
            "(HEALPix), times R = p<sub>L</sub>(Δt) / p<sub>U</sub>(Δt), the time-delay density of galaxy-lensed "
            f"pairs (singular isothermal spheres, median delay {ratio.median_lensed_days:.0f} days for these "
            "sources) over that of all the pairs of the selection (Haris et al. 2018; LVK GWTC-4.0 lensing, "
            "arXiv:2512.16347). B<sub>overlap</sub> carries a constant offset from the prior normalization, the same "
            "for every pair.",
            "Significance: the FPP of a pair is the fraction of the pairs of the selection with a higher statistic; "
            "the selection is its own background, valid while almost every pair is unrelated (the LVK uses "
            "simulated unlensed events). This is a screen, as the LVK tier-1 stage: a flagged pair is a candidate "
            "for joint parameter estimation of the two signals, not a lensing detection. Expect about "
            f"{fpp_threshold:g} × {len(df)} ≈ {fpp_threshold * len(df):.0f} flagged pairs from unrelated events alone.",
        ]
        Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
        write_simple_html_report(
            out_report_html, title="Lensed-pair screen (posterior overlap)", paragraphs=paras, images=imgs,
            tables=[("Highest-ranked pairs", top.drop(columns=["ln_B_intrinsic", "ln_B_sky"]).to_html(
                index=False, float_format=lambda x: f"{x:.3g}"))])
        _log(f"report written to {out_report_html}")
    return df

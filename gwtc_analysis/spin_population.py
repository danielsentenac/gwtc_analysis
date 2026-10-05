"""Population of BBH effective spins: the chi_eff distribution and its correlation with the mass ratio.

Hierarchical inference on the BBH events of the catalog (the selection of `hubble_constant`), with the
LVK search-sensitivity injections for the selection effects:

    p(chi_eff | q) = N(chi_eff; mu_0 + alpha (q - 0.5), sigma), truncated to [-1, 1]

(Callister et al. 2021, arXiv:2106.00521; alpha = 0 is the Gaussian chi_eff model). The other spin degrees
of freedom follow isotropic spins of magnitude uniform on [0, 1], reweighted in chi_eff only: the spin factor
of a PE sample or an injection is N(chi_eff) / pi_iso(chi_eff | q), with pi_iso the chi_eff density of such
spins, computed exactly by convolution (a cos(tilt) has density -ln|s| / 2 on [-1, 1]). The PE spin priors
are isotropic and uniform in magnitude, so this factor is all the spin part of the PE weight.

The masses and the redshifts follow a fixed population, that of the `rates` mode: Power Law + Peak (GWTC-3)
and R ∝ (1 + z)^2.9, Planck15. The likelihood is marginalized over the rate (scale-free):

    ln L = Σ_events ln <w f>_PE  -  N ln <w f>_injections,

with w the mass-redshift weight (population over PE prior, or over injection draw) and f the spin factor.
The injections' draw density is read in full (their spins are not drawn isotropically), with the isotropic
spin part of the population in it, as in `hubble_constant`. A posterior point is kept only when the injections
give at least 4 N effective samples and every event at least 10 effective PE samples. Sampler: emcee.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from .report import write_simple_html_report

KAPPA = 2.9                       # redshift evolution of the rates mode
_trapz = getattr(np, "trapezoid", None) or np.trapz
PRIORS = dict(mu0=(-1.0, 1.0), sigma=(0.005, 1.0), alpha=(-5.0, 5.0))
MIN_NEFF_PE = 10                  # as icarogw and the GWTC-4.0 cosmology paper


def _log(msg: str) -> None:
    print(f"[spin_population] {msg}", flush=True)


# ---------------------------------------------------------------------------
# isotropic chi_eff density
# ---------------------------------------------------------------------------
def _s_cdf(s):
    """CDF of s = a cos(tilt), a ~ U(0, 1), cos(tilt) ~ U(-1, 1): density -ln|s| / 2."""
    s = np.clip(np.asarray(s, float), -1.0, 1.0)
    x = np.abs(s)
    with np.errstate(divide="ignore", invalid="ignore"):
        half = np.where(x > 0, 0.5 * x * (1 - np.log(x)), 0.0)
    return 0.5 + np.sign(s) * half


def iso_chi_eff_density(q: float, chi: np.ndarray, n: int = 4001) -> np.ndarray:
    """pi_iso(chi_eff | q) for isotropic spins of magnitude uniform on [0, 1] (exact bin masses, convolved)."""
    edges = np.linspace(-1, 1, n + 1)
    w = np.diff(_s_cdf(edges))                                   # masses of s1 in its bins
    ds = edges[1] - edges[0]
    c = 0.5 * (edges[1:] + edges[:-1])
    # chi (1 + q) = s1 + q s2: the density of t = q s2 on the same grid, then convolution
    if q < 1e-6:
        conv = w / ds
        grid_x = c
    else:
        w2 = np.diff(_s_cdf(edges / q))                         # P(q s2 in bin)
        conv = np.convolve(w, w2) / ds                         # density of s1 + q s2 on a grid of step ds
        grid_x = np.arange(len(conv)) * ds + 2 * c[0]
    x = np.asarray(chi, float) * (1 + q)
    return np.interp(x, grid_x, conv, left=0.0, right=0.0) * (1 + q)


class IsoTable:
    """pi_iso(chi_eff | q) tabulated on a (q, chi_eff) grid, bilinear interpolation."""

    def __init__(self, nq: int = 50, nchi: int = 801):
        self.q = np.linspace(0.02, 1.0, nq)
        self.chi = np.linspace(-1, 1, nchi)
        self.tab = np.array([iso_chi_eff_density(q, self.chi) for q in self.q])

    def __call__(self, chi, q):
        chi, q = np.asarray(chi, float), np.clip(np.asarray(q, float), self.q[0], self.q[-1])
        iq = np.clip(np.searchsorted(self.q, q) - 1, 0, len(self.q) - 2)
        fq = (q - self.q[iq]) / (self.q[iq + 1] - self.q[iq])
        ic = np.clip(np.searchsorted(self.chi, chi) - 1, 0, len(self.chi) - 2)
        fc = np.clip((chi - self.chi[ic]) / (self.chi[ic + 1] - self.chi[ic]), 0, 1)
        t = self.tab
        lo = t[iq, ic] * (1 - fc) + t[iq, ic + 1] * fc
        hi = t[iq + 1, ic] * (1 - fc) + t[iq + 1, ic + 1] * fc
        return lo * (1 - fq) + hi * fq


# ---------------------------------------------------------------------------
# population model
# ---------------------------------------------------------------------------
def trunc_normal(x, mu, sigma):
    from scipy.special import erf

    norm = 0.5 * (erf((1 - mu) / (np.sqrt(2) * sigma)) - erf((-1 - mu) / (np.sqrt(2) * sigma)))
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (np.sqrt(2 * np.pi) * sigma * norm)


def mass_redshift_weight(m1_det, m2_det, dl, prior):
    """Fixed population (Power Law + Peak, R ∝ (1+z)^2.9, Planck15) over a detector-frame density `prior`."""
    from astropy.cosmology import Planck15

    from .catalogs import _ln_power_law_peak

    z_g = np.linspace(0, 5, 5001)
    d_g = Planck15.luminosity_distance(z_g).value
    z = np.interp(dl, d_g, z_g)
    ddl = np.interp(z, z_g, np.gradient(d_g, z_g))
    dvc = np.interp(z, z_g, Planck15.differential_comoving_volume(z_g).value)
    m1s, m2s = m1_det / (1 + z), m2_det / (1 + z)
    hi, lo = np.maximum(m1s, m2s), np.minimum(m1s, m2s)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        p = np.exp(_ln_power_law_peak(hi, lo)) * dvc * (1 + z) ** (KAPPA - 1) / ((1 + z) ** 2 * ddl) / prior
    return np.where(np.isfinite(p), p, 0.0)


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def _pe_file(row, cache: Path) -> Path:
    """The event's PE file in `cache`/files, else downloaded from Zenodo; names rounded by +-2 s are accepted,
    as in `hubble_constant`."""
    from . import hubble_constant as hc
    from . import parameters_estimation as pe

    names = [row["event"]] + [hc._full_name(row["gps"] + d) for d in (-2, -1, 1, 2)]
    files = cache / "files"
    for nm in names:
        local = sorted(files.glob(f"*{nm}*PEDataRelease*.h*5"))
        if local:
            return Path(hc._pick_pe_file([dict(filename=p.name, path=p) for p in local])["path"])
    index = pe.build_zenodo_pe_index(cache_dir=str(cache / "index"), force_refresh=False)
    cands = next((index[nm] for nm in names + [row.get("common_name") or ""] if nm in index), None)
    if not cands:
        raise ValueError(f"no Zenodo PE release found for {row['event']}")
    entry = hc._pick_pe_file(cands)
    files.mkdir(parents=True, exist_ok=True)
    dest = files / entry["filename"]
    tmp = dest.with_suffix(dest.suffix + ".part")
    _log(f"downloading the PE file of {row['event']}")
    pe._download_http_with_progress(entry["url"], tmp, desc=row["event"])
    tmp.replace(dest)
    return dest


def _event_samples(row, cache: Path, n: int, rng) -> dict:
    """mass_1, mass_2, luminosity_distance, chi_eff, mass_ratio and the PE distance-prior density of an event."""
    import h5py

    from . import hubble_constant as hc

    name = row["event"]
    out_dir = cache / "samples_spin"
    out = out_dir / f"{name}.h5"
    if not out.exists():
        src = _pe_file(row, cache)
        out_dir.mkdir(parents=True, exist_ok=True)
        with h5py.File(src, "r") as f:
            lab = hc._choose_label(f)
            if lab is None:
                raise ValueError(f"{name}: no IMRPhenomXPHM analysis in {src.name}")
            ps = f[lab]["posterior_samples"]
            desc = ""
            try:
                v = f[lab]["priors"]["analytic"]["luminosity_distance"][()]
                v = v[0] if hasattr(v, "__len__") and not isinstance(v, (bytes, str)) and len(v) else v
                desc = v.decode() if isinstance(v, bytes) else str(v)
            except (KeyError, TypeError, ValueError):
                pass
            cols = {c: np.asarray(ps[c], float) for c in ("mass_1", "mass_2", "luminosity_distance", "chi_eff")}
        part = out.with_suffix(".part")
        with h5py.File(part, "w") as o:
            for c, v in cols.items():
                o.create_dataset(c, data=v, compression="gzip")
            o.attrs.update(label=lab, prior_desc=desc[:300], source=src.name)
        part.replace(out)
    with h5py.File(out, "r") as h:
        d = {c: h[c][:] for c in ("mass_1", "mass_2", "luminosity_distance", "chi_eff")}
        d["label"], d["prior_desc"] = str(h.attrs["label"]), str(h.attrs["prior_desc"])
    idx = rng.permutation(len(d["chi_eff"]))[:n]
    for c in ("mass_1", "mass_2", "luminosity_distance", "chi_eff"):
        d[c] = d[c][idx]
    d["mass_ratio"] = np.minimum(d["mass_1"], d["mass_2"]) / np.maximum(d["mass_1"], d["mass_2"])
    return d


def prepare(release: str, far_threshold: float, snr_threshold: float, min_mass: float, exclude: Iterable[str],
            pe_cache: Path, pe_samples: int, sensitivity_file=None, seed: int = 1,
            max_injections: Optional[int] = 250_000) -> dict:
    from . import hubble_constant as hc

    rng = np.random.default_rng(seed)
    events = hc.select_h0_events(release, far_threshold, min_mass, exclude)
    _log(f"{len(events)} BBH events (release {release}, FAR <= {far_threshold:g}/yr, masses >= {min_mass:g} M_sun)")
    iso = IsoTable()
    ev = []
    for _, r in events.iterrows():
        d = _event_samples(r, pe_cache, pe_samples, rng)
        prior, _ = hc.pe_distance_prior(d["luminosity_distance"], d["prior_desc"], r["run"])
        w = mass_redshift_weight(d["mass_1"], d["mass_2"], d["luminosity_distance"], prior)
        ev.append(dict(event=r["event"], w=w, chi=d["chi_eff"], q=d["mass_ratio"],
                       iso=np.maximum(iso(d["chi_eff"], d["mass_ratio"]), 1e-12)))
    inj = hc.detector_frame_injections(hc.h0_sensitivity_path(sensitivity_file, release), far_threshold, snr_threshold)
    wi = mass_redshift_weight(inj["mass_1"], inj["mass_2"], inj["luminosity_distance"], inj["prior"])
    keep = np.flatnonzero(wi > 0)
    ntotal = inj["ntotal"]
    if max_injections and len(keep) > max_injections:        # a random subset: the detectable fraction is unbiased
        frac = max_injections / len(keep)
        keep = rng.choice(keep, max_injections, replace=False)
        ntotal = ntotal * frac
    injd = dict(w=wi[keep], chi=inj["chi_eff"][keep], q=inj["mass_ratio"][keep],
                iso=np.maximum(iso(inj["chi_eff"][keep], inj["mass_ratio"][keep]), 1e-12), ntotal=ntotal)
    _log(f"{len(keep)} found injections with a population weight used")
    return dict(events=ev, inj=injd)


# ---------------------------------------------------------------------------
# likelihood and sampling
# ---------------------------------------------------------------------------
class Likelihood:
    def __init__(self, data: dict, correlated: bool = True):
        ev = data["events"]
        self.n = len(ev)
        self.correlated = correlated
        self.seg = np.cumsum([0] + [len(e["chi"]) for e in ev])[:-1]
        self.counts = np.array([len(e["chi"]) for e in ev])
        self.w = np.concatenate([e["w"] / e["iso"] for e in ev])
        self.chi = np.concatenate([e["chi"] for e in ev])
        self.q = np.concatenate([e["q"] for e in ev])
        inj = data["inj"]
        self.iw, self.ichi, self.iq, self.nt = inj["w"] / inj["iso"], inj["chi"], inj["q"], inj["ntotal"]
        self.iw_mass = inj["w"]                     # population mass-redshift weights: averages over the population

    def params(self, theta):
        mu0, sigma = theta[0], theta[1]
        alpha = theta[2] if self.correlated else 0.0
        return mu0, sigma, alpha

    def ln_prob(self, theta) -> float:
        mu0, sigma, alpha = self.params(theta)
        names = ("mu0", "sigma", "alpha") if self.correlated else ("mu0", "sigma")
        for nm, v in zip(names, theta):
            lo, hi = PRIORS[nm]
            if not lo <= v <= hi:
                return -np.inf
        f = self.w * trunc_normal(self.chi, mu0 + alpha * (self.q - 0.5), sigma)
        s1 = np.add.reduceat(f, self.seg)
        per = s1 / self.counts
        fi = self.iw * trunc_normal(self.ichi, mu0 + alpha * (self.iq - 0.5), sigma)
        s = fi.sum()
        if s <= 0 or not np.all(per > 0):
            return -np.inf
        if s * s / (fi * fi).sum() < 4 * self.n:                       # effective injections
            return -np.inf
        if np.min(s1 * s1 / np.add.reduceat(f * f, self.seg)) < MIN_NEFF_PE:   # effective PE samples per event
            return -np.inf
        return float(np.log(per).sum() - self.n * np.log(s / self.nt))

    def neff(self, theta) -> tuple[float, float]:
        mu0, sigma, alpha = self.params(theta)
        fi = self.iw * trunc_normal(self.ichi, mu0 + alpha * (self.iq - 0.5), sigma)
        f = self.w * trunc_normal(self.chi, mu0 + alpha * (self.q - 0.5), sigma)
        s1 = np.add.reduceat(f, self.seg); s2 = np.add.reduceat(f * f, self.seg)
        return float(fi.sum() ** 2 / (fi * fi).sum()), float(np.min(s1 * s1 / s2))


def sample(lik: Likelihood, n_walkers: int = 24, n_steps: int = 1500, burn: int = 500, seed: int = 1) -> pd.DataFrame:
    import emcee

    rng = np.random.default_rng(seed)
    ndim = 3 if lik.correlated else 2
    start = np.array([0.05, 0.1, 0.0][:ndim])
    p0 = start + 1e-2 * rng.normal(size=(n_walkers, ndim))
    p0[:, 1] = np.abs(p0[:, 1])
    s = emcee.EnsembleSampler(n_walkers, ndim, lik.ln_prob)
    s.run_mcmc(p0, n_steps, progress=False)
    chain = s.get_chain(discard=burn, flat=True)
    cols = ["mu0", "sigma", "alpha"][:ndim]
    df = pd.DataFrame(chain, columns=cols)
    df.attrs["acceptance"] = float(np.mean(s.acceptance_fraction))
    try:
        df.attrs["autocorr"] = float(np.max(s.get_autocorr_time(discard=burn, quiet=True)))
    except Exception:
        df.attrs["autocorr"] = float("nan")
    return df


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def predictive(post: pd.DataFrame, q: float, chi=np.linspace(-1, 1, 401), n: int = 2000, seed: int = 1):
    d = post.sample(min(n, len(post)), random_state=seed)
    a = d["alpha"].to_numpy() if "alpha" in d else np.zeros(len(d))
    mu = d["mu0"].to_numpy() + a * (q - 0.5)
    return chi, trunc_normal(chi[None, :], mu[:, None], d["sigma"].to_numpy()[:, None])


def negative_fraction(post: pd.DataFrame, lik: "Likelihood", n: int = 300, seed: int = 1) -> np.ndarray:
    """Fraction of the population with chi_eff < 0, averaged over its mass ratios (the injections' population
    weights), for posterior draws."""
    from scipy.special import erf

    d = post.sample(min(n, len(post)), random_state=seed)
    w = lik.iw_mass / lik.iw_mass.sum()
    out = []
    for _, r in d.iterrows():
        mu = r["mu0"] + (r["alpha"] if "alpha" in r else 0.0) * (lik.iq - 0.5)
        s = r["sigma"] * np.sqrt(2)
        lo, hi = erf((-1 - mu) / s), erf((1 - mu) / s)
        out.append(float(np.sum(w * (erf(-mu / s) - lo) / (hi - lo))))
    return np.array(out)


def plot(post: pd.DataFrame, out_png: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.3), dpi=150)
    fig.patch.set_facecolor(surface)
    ax = axes[0]
    for q, col in ((1.0, "#2a78d6"), (0.5, "#eb6834"), (0.25, "#1a9e77")):
        chi, p = predictive(post, q)
        lo, med, hi = np.percentile(p, [5, 50, 95], axis=0)
        ax.fill_between(chi, lo, hi, color=col, alpha=0.2, lw=0)
        ax.plot(chi, med, color=col, lw=2, label=f"q = {q:g}")
    ax.axvline(0, color=ink2, lw=0.8, ls=":")
    ax.set_xlim(-0.6, 0.8)
    ax.set_xlabel("Effective spin χ$_{\\rm eff}$", color=ink2)
    ax.set_ylabel("Population density", color=ink2)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=ink2, title="mass ratio", title_fontsize=8.5)
    ax.set_title("Population χ$_{\\rm eff}$ distribution (median, 90%)", color=ink, fontsize=10.5, loc="left")
    ax = axes[1]
    if "alpha" in post:
        ax.scatter(post["alpha"], post["mu0"], s=2, color="#2a78d6", alpha=0.25, lw=0)
        ax.axvline(0, color=ink2, lw=0.8, ls=":")
        ax.set_xlabel("α: slope of the mean χ$_{\\rm eff}$ with q", color=ink2)
        ax.set_ylabel("μ$_0$: mean χ$_{\\rm eff}$ at q = 0.5", color=ink2)
        ax.set_title(f"P(α < 0) = {np.mean(post['alpha'] < 0):.3f}", color=ink, fontsize=10.5, loc="left")
    else:
        ax.scatter(post["mu0"], post["sigma"], s=2, color="#2a78d6", alpha=0.25, lw=0)
        ax.set_xlabel("μ$_0$", color=ink2); ax.set_ylabel("σ", color=ink2)
    for a in axes:
        a.set_facecolor(surface); a.grid(color=grid, lw=0.7)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            a.spines[s].set_color(grid)
        a.tick_params(colors=ink2, labelsize=9)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def run_spin_population(
    sensitivity_release: Optional[str] = None,
    sensitivity_file=None,
    far_threshold: float = 0.25,
    snr_threshold: float = 10.0,
    min_mass: float = 3.0,
    exclude: Optional[Iterable[str]] = None,
    pe_cache: Optional[str | Path] = None,
    pe_samples: int = 5000,
    max_injections: Optional[int] = 250_000,
    correlated: bool = True,
    n_walkers: int = 24,
    n_steps: int = 1500,
    out_report_html: Optional[str | Path] = "spin_population.html",
    out_summary_tsv: Optional[str | Path] = "spin_population.tsv",
    plots_dir: str | Path = "spin_population_plots",
    data: Optional[dict] = None,
) -> pd.DataFrame:
    from . import hubble_constant as hc

    release = sensitivity_release or hc.H0_DEFAULT_RELEASE
    if data is None:
        cache = Path(pe_cache).expanduser() if pe_cache else hc.default_pe_cache()
        data = prepare(release, far_threshold, snr_threshold, min_mass,
                       hc.H0_DEFAULT_EXCLUDE if exclude is None else exclude, cache, pe_samples, sensitivity_file,
                       max_injections=max_injections)
    lik = Likelihood(data, correlated=correlated)
    _log(f"sampling ({'mu0, sigma, alpha' if correlated else 'mu0, sigma'}): {n_walkers} walkers x {n_steps} steps")
    post = sample(lik, n_walkers=n_walkers, n_steps=n_steps, burn=n_steps // 3)
    med = post.median().to_numpy()
    neff_inj, neff_pe = lik.neff(med)
    rows = [dict(parameter=c, median=float(post[c].median()), low_90=float(post[c].quantile(0.05)),
                 high_90=float(post[c].quantile(0.95))) for c in post.columns]
    neg = np.percentile(negative_fraction(post, lik), [5, 50, 95])
    table = pd.DataFrame(rows)
    p_alpha = float(np.mean(post["alpha"] < 0)) if "alpha" in post else None
    _log("; ".join(f"{r['parameter']} = {r['median']:.3f} [{r['low_90']:.3f}, {r['high_90']:.3f}]" for r in rows)
         + (f"; P(alpha < 0) = {p_alpha:.3f}" if p_alpha is not None else ""))
    if out_summary_tsv:
        Path(out_summary_tsv).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_summary_tsv, sep="\t", index=False, float_format="%.4g")
        post.to_csv(Path(out_summary_tsv).with_suffix(".posterior.tsv"), sep="\t", index=False, float_format="%.5g")
    if out_report_html:
        img = plot(post, Path(plots_dir) / "spin_population.png")
        t = table.set_index("parameter")
        f = lambda k: f"{t.loc[k, 'median']:.2f} (90%: {t.loc[k, 'low_90']:.2f} to {t.loc[k, 'high_90']:.2f})"
        paras = [
            f"<b>χ<sub>eff</sub> population</b> of {lik.n} BBH events: mean μ<sub>0</sub> = {f('mu0')} at q = 0.5, width "
            f"σ = {f('sigma')}" + (f", slope with the mass ratio α = {f('alpha')}, P(α &lt; 0) = {p_alpha:.3f}" if p_alpha
                                    is not None else "") + f". A fraction {neg[1]:.2f} (90%: {neg[0]:.2f}–{neg[2]:.2f}) of the binaries, averaged "
            "over the population's mass ratios, have χ<sub>eff</sub> &lt; 0 (GWTC-4.0: 0.24–0.42, arXiv:2508.18083).",
            "Model: χ<sub>eff</sub> | q ~ N(μ<sub>0</sub> + α (q − 0.5), σ) on [−1, 1] (Callister et al. 2021; α = 0 is "
            "the Gaussian model), the other spin degrees of freedom isotropic. α &lt; 0 means that unequal-mass "
            "binaries have larger effective spins, the correlation Callister et al. found at 98.7% credibility; "
            "GWTC-4.0 also finds evidence for a χ<sub>eff</sub>–q correlation (LVK, arXiv:2508.18083). Masses and "
            "redshifts follow the fixed population of the rates mode (Power Law + Peak, R ∝ (1+z)<sup>2.9</sup>).",
            f"Selection: {len(lik.iw)} found injections of {release}; at the median, {neff_inj:.0f} effective injections "
            f"(threshold 4N = {4 * lik.n}) and at least {neff_pe:.0f} effective PE samples per event (threshold "
            f"{MIN_NEFF_PE}). emcee: "
            f"{n_walkers} walkers, {n_steps} steps, acceptance {post.attrs.get('acceptance', float('nan')):.2f}, "
            f"autocorrelation time {post.attrs.get('autocorr', float('nan')):.0f} steps.",
        ]
        Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
        write_simple_html_report(out_report_html, title="BBH spin population", paragraphs=paras, images=[img],
                                 tables=[("Posterior", table.to_html(index=False, float_format=lambda x: f"{x:.3g}"))])
        _log(f"report written to {out_report_html}")
    return table

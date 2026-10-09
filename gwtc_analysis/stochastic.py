"""Predicted gravitational-wave background of the unresolved compact binaries, against the O3 upper limit.

The energy density of the background, per logarithmic frequency and in units of the critical density, is
(e.g. Phinney 2001; Regimbau 2011)

    Omega_GW(f) = f / (rho_c c^2) * Integral dz  R(z) <dE/df_s>(f (1 + z)) / [(1 + z) H(z)],

with R(z) the merger rate per comoving volume and source-frame time, <dE/df_s> the energy spectrum of one
merger averaged over the population, and rho_c c^2 = 3 H0^2 c^2 / (8 pi G).

- BBH: R(z) = R(0.2) psi(z) / psi(0.2), with psi the Madau-Dickinson shape and the mass model (Power Law +
  Peak or Multi Peak, recognized from the posterior's parameters) of each posterior draw of a spectral-siren run (`hubble_constant`), and R(0.2), where the catalog
  measures the rate best, from the `rates` mode. Beyond the farthest detected events (z_h) the fitted shape
  is the prior's: by default (`high_z="sfr"`) the rate follows the star-formation history there, joined
  continuously at z_h; `high_z="posterior"` keeps the fitted shape up to z = 10 (a sensitivity check: most
  of Omega_GW(25 Hz) then comes from z > z_h). The energy spectrum is the phenomenological
  inspiral-merger-ringdown one of Ajith et al. 2008 (arXiv:0710.2335), non-spinning; its merger and
  ringdown overestimate the radiated energy (4.2 M_sun c^2 for 36 + 29 M_sun, against ~3.0), which affects
  Omega_GW(25 Hz) at the ~20% level (merger and ringdown are 21% of it).
- BNS and NSBH: the local rates of the `rates` mode, a rate following the star-formation history (Madau &
  Dickinson 2014, no delay time) and the inspiral spectrum up to the innermost stable circular orbit.

The cosmology is Planck15, as the volumes of the `rates` mode. The upper limit from the data through April 2025
is Omega_GW(25 Hz) <= 2.0e-9 (95%) for a power law of index 2/3 (LVK 2026, arXiv:2608.23477), which predicts
Omega_CBC(25 Hz) = 6.3 (+5.0 / -2.2) e-10 from GWTC-5.0; O3 alone gave 3.4e-9 (arXiv:2101.12130).
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .report import write_simple_html_report

G, C, MSUN = 6.67430e-11, 2.99792458e8, 1.98847e30
MPC, YR = 3.0856775814913673e22, 3.15576e7
LIMIT = dict(omega_25=2.0e-9, alpha=2 / 3, fmin=20.0, fmax=90.6, ref="LVK 2026, data through April 2025, arXiv:2608.23477")
O3_LIMIT = dict(omega_25=3.4e-9, ref="LVK 2021, arXiv:2101.12130")
LVK_PREDICTION = dict(median=6.3e-10, plus=5.0e-10, minus=2.2e-10, ref="GWTC-5.0, arXiv:2608.23477")
Z_HORIZON = 1.0                   # farthest detected BBH (median distance), when the work directory does not tell
# Ajith et al. 2008, Table I: f = (a eta^2 + b eta + c) / (pi M)
_AJITH = dict(merg=(2.9740e-1, 4.4810e-2, 9.5560e-2), ring=(5.9411e-1, 8.9794e-2, 1.9111e-1),
              sigma=(5.0801e-1, 7.7515e-2, 2.2369e-2), cut=(8.4845e-1, 1.2848e-1, 2.7299e-1))
_trapz = getattr(np, "trapezoid", None) or np.trapz


def _log(msg: str) -> None:
    print(f"[stochastic] {msg}", flush=True)


# ---------------------------------------------------------------------------
# energy spectra (source frame, J/Hz)
# ---------------------------------------------------------------------------
def _ajith_freqs(m_tot_kg: np.ndarray, eta: np.ndarray) -> dict:
    tm = G * m_tot_kg / C ** 3
    return {k: (a * eta ** 2 + b * eta + c) / (np.pi * tm) for k, (a, b, c) in _AJITH.items()}


def dE_df_imr(f: np.ndarray, m1: np.ndarray, m2: np.ndarray) -> np.ndarray:
    """dE/df of non-spinning binary black holes (Ajith et al. 2008), shape (n_binaries, n_f)."""
    m1, m2 = np.atleast_1d(m1)[:, None] * MSUN, np.atleast_1d(m2)[:, None] * MSUN
    m, eta = m1 + m2, m1 * m2 / (m1 + m2) ** 2
    mc = m * eta ** 0.6
    fr = _ajith_freqs(m, eta)
    f = np.asarray(f, float)[None, :]
    pre = (np.pi * G) ** (2 / 3) * mc ** (5 / 3) / 3
    insp = f ** (-1 / 3)
    merg = f ** (2 / 3) / fr["merg"]
    ring = (f / (1 + ((f - fr["ring"]) / (fr["sigma"] / 2)) ** 2)) ** 2 / (fr["merg"] * fr["ring"] ** (4 / 3))
    out = np.where(f < fr["merg"], insp, np.where(f < fr["ring"], merg, np.where(f < fr["cut"], ring, 0.0)))
    return pre * out


def dE_df_inspiral(f: np.ndarray, m1: np.ndarray, m2: np.ndarray) -> np.ndarray:
    """Newtonian inspiral dE/df up to the innermost stable circular orbit, shape (n_binaries, n_f)."""
    m1, m2 = np.atleast_1d(m1)[:, None] * MSUN, np.atleast_1d(m2)[:, None] * MSUN
    m = m1 + m2
    mc = (m1 * m2) ** 0.6 / m ** 0.2
    f_isco = C ** 3 / (6 ** 1.5 * np.pi * G * m)
    f = np.asarray(f, float)[None, :]
    return np.where(f < f_isco, (np.pi * G) ** (2 / 3) * mc ** (5 / 3) / 3 * f ** (-1 / 3), 0.0)


def radiated_energy(m1: float, m2: float, spectrum=dE_df_imr) -> float:
    """Total energy of one binary (M_sun c^2), integrating its spectrum."""
    f = np.geomspace(1.0, 2e4, 20000)
    return float(_trapz(spectrum(f, m1, m2)[0], f) / (MSUN * C ** 2))


# ---------------------------------------------------------------------------
# populations
# ---------------------------------------------------------------------------
def sample_plp(rng, n: int, alpha, beta, mmin, mmax, delta_m, mu_g, sigma_g, lambda_peak) -> tuple[np.ndarray, np.ndarray]:
    """(m1, m2) source-frame samples of the Power Law + Peak model (icarogw's PowerLawPeak with the
    m1m2_conditioned_lowpass smoothing, as in the spectral siren):

    p(m1) ∝ [(1 - λ) PL(m1; -α, mmin, mmax) + λ N(m1; μ, σ)] S(m1),  p(m2 | m1) ∝ m2^β S(m2) on [mmin, m1].
    """
    g = np.linspace(mmin, max(mmax, mu_g + 6 * sigma_g), 4000)
    pl = np.where(g <= mmax, g ** (-alpha), 0.0)
    pl /= _trapz(pl, g)
    peak = np.exp(-0.5 * ((g - mu_g) / sigma_g) ** 2) / (sigma_g * np.sqrt(2 * np.pi))
    p1 = ((1 - lambda_peak) * pl + lambda_peak * peak) * _smooth(g, mmin, delta_m)
    c1 = np.cumsum(p1); c1 /= c1[-1]
    m1 = np.interp(rng.uniform(size=n), c1, g)
    return m1, _sample_m2(rng, m1, beta, mmin, delta_m)


def _sample_m2(rng, m1, beta, mmin, delta_m) -> np.ndarray:
    """m2 | m1 ∝ m2^β S(m2) on [mmin, m1] (icarogw's m1m2_conditioned_lowpass, for both mass models)."""
    u = rng.uniform(size=len(m1))
    m2 = np.empty(len(m1))
    for i, a in enumerate(m1):
        x = np.linspace(mmin, a, 200)
        w = x ** beta * _smooth(x, mmin, delta_m)
        cw = np.cumsum(w)
        m2[i] = np.interp(u[i], cw / cw[-1], x) if cw[-1] > 0 else a
    return m2


def _truncated_gaussian(g, mu, sigma, lo, hi):
    """Gaussian density truncated to [lo, hi] and normalized there (icarogw's TruncatedGaussian)."""
    from scipy.special import erf

    norm = 0.5 * (erf((hi - mu) / (np.sqrt(2) * sigma)) - erf((lo - mu) / (np.sqrt(2) * sigma)))
    p = np.exp(-0.5 * ((g - mu) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    return np.where((g >= lo) & (g <= hi), p / max(norm, 1e-300), 0.0)


def sample_mltp(rng, n: int, alpha, beta, mmin, mmax, delta_m, mu_g_low, sigma_g_low, lambda_g_low,
                mu_g_high, sigma_g_high, lambda_g) -> tuple[np.ndarray, np.ndarray]:
    """(m1, m2) source-frame samples of the Multi Peak model (icarogw's massprior_MultiPeak, i.e.
    PowerLawTwoGaussians, with the m1m2_conditioned_lowpass smoothing, as in the spectral siren):

    p(m1) ∝ [(1 - λ) PL(m1; -α, mmin, mmax) + λ λ_low N_low(m1) + λ (1 - λ_low) N_high(m1)] S(m1),
    each Gaussian truncated to [mmin, μ + 5σ]; p(m2 | m1) ∝ m2^β S(m2) on [mmin, m1].
    """
    top = max(mmax, mu_g_low + 5 * sigma_g_low, mu_g_high + 5 * sigma_g_high)
    g = np.linspace(mmin, top, 6000)
    pl = np.where(g <= mmax, g ** (-alpha), 0.0)
    pl /= _trapz(pl, g)
    low = _truncated_gaussian(g, mu_g_low, sigma_g_low, mmin, mu_g_low + 5 * sigma_g_low)
    high = _truncated_gaussian(g, mu_g_high, sigma_g_high, mmin, mu_g_high + 5 * sigma_g_high)
    p1 = ((1 - lambda_g) * pl + lambda_g * lambda_g_low * low + lambda_g * (1 - lambda_g_low) * high) \
        * _smooth(g, mmin, delta_m)
    c1 = np.cumsum(p1); c1 /= c1[-1]
    m1 = np.interp(rng.uniform(size=n), c1, g)
    return m1, _sample_m2(rng, m1, beta, mmin, delta_m)


# the BBH mass models: the posterior columns each one needs, its name, and its sampler
_COMMON = ("alpha", "beta", "mmin", "mmax", "delta_m", "gamma", "kappa", "zp")
MASS_MODELS = {
    "plp": dict(name="Power Law + Peak", params=("mu_g", "sigma_g", "lambda_peak")),
    "mltp": dict(name="Multi Peak", params=("mu_g_low", "sigma_g_low", "lambda_g_low", "mu_g_high", "sigma_g_high",
                                             "lambda_g")),
}


def mass_model_of(post: pd.DataFrame) -> str:
    """The BBH mass model of a spectral-siren posterior, from its parameters (plp or mltp)."""
    cols = set(post.columns)
    for key, m in MASS_MODELS.items():
        if set(_COMMON) | set(m["params"]) <= cols:
            return key
    raise ValueError("the posterior has the parameters of neither the Power Law + Peak nor the Multi Peak model "
                     f"(columns: {sorted(cols)}): a hubble_constant run with --mass-model plp or mltp is needed")


def sample_masses(rng, n: int, row, model: str) -> tuple[np.ndarray, np.ndarray]:
    base = (row["alpha"], row["beta"], row["mmin"], row["mmax"], row["delta_m"])
    if model == "mltp":
        return sample_mltp(rng, n, *base, row["mu_g_low"], row["sigma_g_low"], row["lambda_g_low"],
                           row["mu_g_high"], row["sigma_g_high"], row["lambda_g"])
    return sample_plp(rng, n, *base, row["mu_g"], row["sigma_g"], row["lambda_peak"])


def _smooth(m, mmin, delta):
    m = np.asarray(m, float)
    out = np.where(m >= mmin + delta, 1.0, 0.0)
    mid = (m > mmin) & (m < mmin + delta)
    x = m[mid] - mmin
    with np.errstate(over="ignore"):
        out[mid] = 1.0 / (np.exp(delta / x + delta / (x - delta)) + 1.0)
    return out


def sample_bns(rng, n):
    return rng.uniform(1.0, 2.5, n), rng.uniform(1.0, 2.5, n)


def sample_nsbh(rng, n):
    u = rng.uniform(size=n)
    a, lo, hi = 2.35, 2.5, 40.0                                  # m_BH ∝ m^-2.35 on [2.5, 40] (rates mode)
    mbh = (lo ** (1 - a) + u * (hi ** (1 - a) - lo ** (1 - a))) ** (1 / (1 - a))
    return mbh, rng.uniform(1.0, 2.5, n)


# ---------------------------------------------------------------------------
# background
# ---------------------------------------------------------------------------
def _cosmology():
    from astropy.cosmology import Planck15

    h0 = Planck15.H0.to("1/s").value
    return h0, Planck15


def omega_gw(f: np.ndarray, rate_z, dEdf_s: callable, zmax: float = 10.0, nz: int = 400) -> np.ndarray:
    """Omega_GW(f) for R(z) = rate_z(z) (Gpc^-3 yr^-1) and <dE/df_s>(f_s) = dEdf_s(f_s) (J/Hz)."""
    h0, cosmo = _cosmology()
    z = np.linspace(0, zmax, nz)
    hz = h0 * cosmo.efunc(z)
    rho_c2 = 3 * h0 ** 2 * C ** 2 / (8 * np.pi * G)
    r_si = np.asarray(rate_z(z), float) / (1e3 * MPC) ** 3 / YR                 # m^-3 s^-1
    fs = np.asarray(f, float)[:, None] * (1 + z)[None, :]
    integrand = r_si[None, :] * dEdf_s(fs) / ((1 + z) * hz)[None, :]
    return np.asarray(f, float) / rho_c2 * _trapz(integrand, z, axis=1)


def mean_spectrum(m1, m2, spectrum, fgrid=np.geomspace(1.0, 3e4, 600)):
    """<dE/df_s> over a population, as an interpolating function of f_s (log-log)."""
    e = spectrum(fgrid, m1, m2).mean(axis=0)
    loge = np.log(np.where(e > 0, e, 1e-300))
    return lambda fs: np.where(fs <= fgrid[-1], np.exp(np.interp(np.log(fs), np.log(fgrid), loge)), 0.0)


def read_rates(path: str | Path) -> dict:
    """Median and 90% interval of the BNS, NSBH and BBH(z = 0.2) rates of a `rates` TSV."""
    df = pd.read_csv(path, sep="\t")
    out = {}
    for pop in ("BNS", "NSBH"):
        r = df[df["population"] == pop].iloc[0]
        out[pop] = (float(r["rate_05"]), float(r["rate_median"]), float(r["rate_95"]))
    bbh = df[(df["population"] == "BBH") & (df["z_ref"] > 0)]
    if bbh.empty:
        raise ValueError(f"{path}: no evolving BBH rate (z_ref > 0)")
    r = bbh.iloc[0]
    out["BBH"] = (float(r["rate_05"]), float(r["rate_median"]), float(r["rate_95"]))
    out["BBH_z_ref"] = float(r["z_ref"])
    return out


def _lognormal_draws(rng, q, n):
    """Draws matching a median and a 90% interval (log-normal, average width)."""
    lo, med, hi = q
    s = (np.log(hi) - np.log(lo)) / (2 * 1.645)
    return np.exp(rng.normal(np.log(med), s, n))


def bbh_rate(z, gamma, kappa, zp, r_ref, z_ref, high_z="sfr", z_h=Z_HORIZON):
    """BBH rate (Gpc^-3 yr^-1): the fitted shape normalized at z_ref; beyond z_h, the star-formation shape."""
    from .rate_evolution import MD14, md_shape

    z = np.asarray(z, float)
    psi = md_shape(z, gamma, kappa, zp)[0] / md_shape(z_ref, gamma, kappa, zp)[0, 0]
    if high_z == "sfr":
        join = md_shape(z_h, gamma, kappa, zp)[0, 0] / md_shape(z_ref, gamma, kappa, zp)[0, 0]
        sfr = md_shape(z, **MD14)[0] / md_shape(z_h, **MD14)[0, 0]
        psi = np.where(z <= z_h, psi, join * sfr)
    elif high_z != "posterior":
        raise ValueError(f"high_z must be 'sfr' or 'posterior', not {high_z!r}")
    return r_ref * psi


def predict(spectral_post: pd.DataFrame, rates: dict, f: np.ndarray, n_draws: int = 200, n_masses: int = 1500,
            seed: int = 1, high_z: str = "sfr", z_h: float = Z_HORIZON) -> dict:
    """Omega_GW(f) draws for BBH, BNS, NSBH and their sum."""
    from .rate_evolution import MD14, md_shape

    rng = np.random.default_rng(seed)
    model = mass_model_of(spectral_post)
    d = spectral_post.sample(min(n_draws, len(spectral_post)), random_state=seed).reset_index(drop=True)
    n = len(d)
    zr = rates["BBH_z_ref"]
    r_bbh = _lognormal_draws(rng, rates["BBH"], n)
    out = {"BBH": np.empty((n, len(f)))}
    for i, row in d.iterrows():
        m1, m2 = sample_masses(rng, n_masses, row, model)
        spec = mean_spectrum(m1, m2, dE_df_imr)
        g, k, p = row["gamma"], row["kappa"], row["zp"]
        out["BBH"][i] = omega_gw(f, lambda z: bbh_rate(z, g, k, p, r_bbh[i], zr, high_z, z_h), spec)
    sfr = lambda z: md_shape(z, **MD14)[0]
    for pop, sampler in (("BNS", sample_bns), ("NSBH", sample_nsbh)):
        spec = mean_spectrum(*sampler(rng, 20000), dE_df_inspiral)
        unit = omega_gw(f, sfr, spec)                                   # for R(0) = 1 Gpc^-3 yr^-1
        out[pop] = _lognormal_draws(rng, rates[pop], n)[:, None] * unit[None, :]
    out["total"] = out["BBH"] + out["BNS"] + out["NSBH"]
    return out


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def plot(f: np.ndarray, om: dict, out_png: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
    cols = {"total": "#0b0b0b", "BBH": "#2a78d6", "BNS": "#eb6834", "NSBH": "#1a9e77"}
    fig, ax = plt.subplots(figsize=(10.5, 4.6), dpi=150)
    fig.patch.set_facecolor(surface); ax.set_facecolor(surface)
    for pop in ("BBH", "BNS", "NSBH", "total"):
        lo, med, hi = np.percentile(om[pop], [5, 50, 95], axis=0)
        if pop != "total":
            ax.fill_between(f, lo, hi, color=cols[pop], alpha=0.18, lw=0, zorder=1)
        ax.plot(f, med, color=cols[pop], lw=2.2 if pop == "total" else 1.6, zorder=3,
                label=f"{pop}: Ω(25 Hz) = {np.percentile(om[pop][:, np.argmin(abs(f - 25))], 50):.1e}")
    lim = LIMIT
    fl = np.geomspace(lim["fmin"], lim["fmax"], 50)
    ax.plot(fl, lim["omega_25"] * (fl / 25) ** lim["alpha"], color="#b00020", lw=2, zorder=4,
            label=f"upper limit, O1–O4 (95%, index 2/3): {lim['omega_25']:.1e} at 25 Hz")
    ax.plot(fl, O3_LIMIT["omega_25"] * (fl / 25) ** lim["alpha"], color="#b00020", lw=1, ls=":", zorder=4,
            label=f"O3 limit: {O3_LIMIT['omega_25']:.1e}")
    lv = LVK_PREDICTION
    ax.errorbar([25], [lv["median"]], yerr=[[lv["minus"]], [lv["plus"]]], fmt="s", color="#7a5195", ms=5, capsize=3,
                zorder=5, label=f"LVK prediction (GWTC-5.0): {lv['median']:.1e}")
    ax.axvline(25, color=ink2, lw=0.7, ls=":", zorder=0)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(f[0], f[-1])
    ax.set_ylim(1e-12, 3e-8)
    ax.set_xlabel("Frequency (Hz)", color=ink2)
    ax.set_ylabel("Ω$_{\\rm GW}$(f)", color=ink2)
    ax.grid(color=grid, lw=0.7, zorder=0, which="both")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(grid)
    ax.tick_params(colors=ink2, labelsize=9)
    ax.legend(frameon=False, fontsize=8, labelcolor=ink2, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.set_title("Predicted background of unresolved compact binaries (median, 90% bands)", color=ink, fontsize=10.5,
                 loc="left")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def horizon_redshift(workdir: Path) -> Optional[float]:
    """Redshift of the farthest event of a hubble_constant work directory (largest median distance)."""
    import h5py

    inputs = Path(workdir) / "inputs.h5"
    if not inputs.exists():
        return None
    from astropy import units as u
    from astropy.cosmology import Planck15, z_at_value

    with h5py.File(inputs, "r") as h:
        dmax = max((float(np.median(h[k]["dl"][:])) for k in h if not k.startswith("_")), default=None)
    return float(z_at_value(Planck15.luminosity_distance, dmax * u.Mpc)) if dmax else None


def _quantiles_25(om: dict, f: np.ndarray) -> pd.DataFrame:
    i25 = int(np.argmin(abs(f - 25)))
    rows = []
    for pop in ("BBH", "BNS", "NSBH", "total"):
        q = np.percentile(om[pop][:, i25], [5, 50, 95])
        rows.append(dict(population=pop, omega_25Hz_median=q[1], omega_25Hz_05=q[0], omega_25Hz_95=q[2],
                         fraction_of_limit=q[1] / LIMIT["omega_25"]))
    return pd.DataFrame(rows)


def run_stochastic(
    spectral_posterior: str | Path,
    rates_tsv: str | Path,
    n_draws: int = 200,
    high_z: str = "sfr",
    z_horizon: Optional[float] = None,
    fmin: float = 5.0,
    fmax: float = 3000.0,
    out_report_html: Optional[str | Path] = "stochastic.html",
    out_summary_tsv: Optional[str | Path] = "stochastic.tsv",
    plots_dir: str | Path = "stochastic_plots",
    seed: int = 1,
) -> pd.DataFrame:
    """Omega_GW(f) of BBH, BNS and NSBH from a spectral-siren posterior and a rates TSV, against the upper limit."""
    src = Path(spectral_posterior).expanduser()
    path = src
    if src.is_dir():
        path = next((src / n for n in ("posterior_reweighted.tsv", "posterior.tsv") if (src / n).exists()), None)
        if path is None:
            raise ValueError(f"no posterior_reweighted.tsv or posterior.tsv in {spectral_posterior}")
    post = pd.read_csv(path, sep="\t")
    model = mass_model_of(post)
    z_h = z_horizon or (horizon_redshift(src) if src.is_dir() else None) or Z_HORIZON
    rates = read_rates(rates_tsv)
    f = np.geomspace(fmin, fmax, 120)
    _log(f"{n_draws} posterior draws of {path.name} ({MASS_MODELS[model]['name']}); rates from "
         f"{Path(rates_tsv).name}; z_h = {z_h:.2f}, "
         f"high_z = {high_z}")
    om = predict(post, rates, f, n_draws=n_draws, seed=seed, high_z=high_z, z_h=z_h)
    other = "posterior" if high_z == "sfr" else "sfr"
    om_alt = predict(post, rates, np.array([25.0]), n_draws=min(n_draws, 100), seed=seed, high_z=other, z_h=z_h)
    table = _quantiles_25(om, f)
    alt = np.percentile(om_alt["BBH"][:, 0], [5, 50, 95])
    tot = table.set_index("population").loc["total"]
    _log(f"Omega_GW(25 Hz) = {tot['omega_25Hz_median']:.2g} [{tot['omega_25Hz_05']:.2g}, {tot['omega_25Hz_95']:.2g}]; "
         f"{tot['fraction_of_limit']:.2f} of the upper limit; BBH with high_z={other}: {alt[1]:.2g}")
    if out_summary_tsv:
        Path(out_summary_tsv).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_summary_tsv, sep="\t", index=False, float_format="%.4g")
        spec = pd.DataFrame({"f_Hz": f, **{f"{p}_{s}": np.percentile(om[p], q, axis=0)
                                           for p in om for s, q in (("05", 5), ("median", 50), ("95", 95))}})
        spec.to_csv(Path(out_summary_tsv).with_suffix(".spectrum.tsv"), sep="\t", index=False, float_format="%.4g")
    if out_report_html:
        img = plot(f, om, Path(plots_dir) / "omega_gw.png")
        b = table.set_index("population")
        lv = LVK_PREDICTION
        hz = ("beyond the farthest detected events (z<sub>h</sub> = {:.2f}) the rate follows the star-formation "
              "history, joined continuously at z<sub>h</sub>").format(z_h) if high_z == "sfr" else \
             "the fitted shape is kept up to z = 10, beyond the farthest detected events (z<sub>h</sub> = {:.2f})".format(z_h)
        paras = [
            f"<b>Ω<sub>GW</sub>(25 Hz) = {tot['omega_25Hz_median']:.1e}</b> (90%: {tot['omega_25Hz_05']:.1e}–"
            f"{tot['omega_25Hz_95']:.1e}) for the unresolved compact binaries: BBH {b.loc['BBH', 'omega_25Hz_median']:.1e}, "
            f"BNS {b.loc['BNS', 'omega_25Hz_median']:.1e}, NSBH {b.loc['NSBH', 'omega_25Hz_median']:.1e}. That is "
            f"{tot['fraction_of_limit']:.2f} of the upper limit from the data through April 2025 "
            f"({LIMIT['omega_25']:.1e} at 25 Hz, 95%, index 2/3; {LIMIT['ref']}): a factor "
            f"{1 / tot['fraction_of_limit']:.1f} below it. The LVK prediction from GWTC-5.0 is {lv['median']:.1e} "
            f"(+{lv['plus']:.1e} / −{lv['minus']:.1e}).",
            "Ω<sub>GW</sub>(f) = f / (ρ<sub>c</sub>c²) ∫ dz R(z) ⟨dE/df<sub>s</sub>⟩(f(1+z)) / [(1+z) H(z)], Planck15. "
            f"BBH: for each of {len(om['BBH'])} draws of the spectral-siren posterior ({path.name}), the "
            f"{MASS_MODELS[model]['name']} masses and the Madau–Dickinson rate shape, normalized to the BBH rate at "
            f"z = {rates['BBH_z_ref']:g} of the rates mode ({rates['BBH'][1]:.1f} Gpc⁻³ yr⁻¹, 90%: {rates['BBH'][0]:.1f}–"
            f"{rates['BBH'][2]:.1f}); {hz}. Non-spinning inspiral–merger–ringdown spectrum of Ajith et al. 2008.",
            f"<b>Systematics.</b> The high-redshift rate matters most: with high_z = {other}, Ω<sub>BBH</sub>(25 Hz) = "
            f"{alt[1]:.1e} (90%: {alt[0]:.1e}–{alt[2]:.1e}), since most of the background comes from beyond the "
            "detected events, where the fitted shape is only the prior's. The Ajith et al. merger and ringdown "
            "overestimate the radiated energy (4.2 M☉c² for a GW150914-like binary, against ~3.0); merger and "
            "ringdown make 21% of Ω<sub>BBH</sub>(25 Hz), so the bias is ~20% at most.",
            f"BNS ({rates['BNS'][1]:.0f} Gpc⁻³ yr⁻¹) and NSBH ({rates['NSBH'][1]:.0f} Gpc⁻³ yr⁻¹): local rates of the "
            "rates mode, a rate following the star-formation history without delay time, the masses of the rates mode "
            "and the inspiral spectrum up to the innermost stable circular orbit. Their uncertainty is dominated by "
            "the few detected events.",
        ]
        Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
        write_simple_html_report(out_report_html, title="Predicted compact-binary background", paragraphs=paras,
                                 images=[img], tables=[("Ω_GW at 25 Hz", table.to_html(
                                     index=False, float_format=lambda x: f"{x:.3g}"))])
        _log(f"report written to {out_report_html}")
    return table

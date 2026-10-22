"""Primordial black holes: upper limits on the fraction of dark matter f_PBH from the merger rates.

A population of primordial black holes (PBHs) of a single mass m, making up a fraction f_PBH of the dark matter,
merges at a rate that grows with f_PBH. Where no merger is seen, or fewer than such a population would give, the
rate limit becomes a limit on f_PBH. Two formation channels:

- **early binaries**, formed in the radiation era and merging today (Raidal et al. 2019, arXiv:1812.01930;
  Hütsi et al. 2021, JCAP 03, 068, arXiv:2012.02786, eq. 2.1):
  R = 1.6e6 Gpc^-3 yr^-1 f^(53/37) (t/t0)^(-34/37) (M/Msun)^(-32/37) eta^(-34/37) S(f), with M = 2m and eta = 1/4.
  S is the suppression of the binaries disrupted by other PBHs and by matter fluctuations:
  `suppression="hutsi"` (default) is S1 * S2 of Hütsi et al. 2021 (eqs. 2.2, 2.3, A.5 and the z = 0 fit 2.7),
  valid at every f; `suppression="lvk"` is the fit S = 2.3e-3 f^-0.65 of the LVK sub-solar papers, which reproduces
  their limits. The fit equals S1 * S2 at f = 1 but keeps only the f^-0.65 of S2: S1 * S2 is 1.3 times larger at
  f = 0.1 and 2.9 times at f = 0.01 (limits about 6 times tighter at 1 Msun), and the fit exceeds 1, which is
  unphysical, below f ~ 1e-4, where the catalog limits lie;
- **late binaries**, formed by capture in PBH clusters (Phukon et al. 2021, arXiv:2105.11449, eq. 2; Clesse &
  García-Bellido 2022, arXiv:2007.06481): R = f * 420 Gpc^-3 yr^-1 (2m)^(10/7) / (m^2)^(5/7), linear in f.

The rate limits come from two kinds of data:

- **0.2-1 Msun**: the sensitive volume-time of the LVK sub-solar-mass searches, O4a and O4b combined, per pipeline
  and chirp-mass bin (LVK 2026, arXiv:2610.20276, data release Zenodo 23346454). No significant candidate: the 90%
  upper limit is R = -ln(1 - 0.9) / VT = 2.30 / VT, as in the LVK papers, and m = 2^(1/5) * chirp mass. Only the
  bins up to m = 1 Msun are used, as in the LVK papers: the searches require a secondary below 1 Msun, so the
  higher chirp-mass bins hold unequal binaries, not the equal-mass binaries of a single PBH mass;
- **above 1 Msun**: the LVK search-sensitivity injections of the rates mode, reweighted to binaries with both
  masses log-uniform in [m e^-w, m e^w], and the GWOSC events with FAR below the threshold whose median chirp mass
  is in the same window, 2^(-1/5) [m e^-w, m e^w] (the chirp mass is well measured, the mass ratio is not). With
  `attribute="window"` (default) every such event may be primordial: the limit is the 90%
  Poisson upper limit for n events, Gamma^-1(0.9; n + 1) / VT (2.30 / VT for none). `attribute="none"` takes n = 0,
  i.e. assumes the observed events are all astrophysical: a tighter and stronger assumption. The early-binary VT
  includes the rate evolution (t/t0)^(-34/37); the late-binary rate is taken constant.

The limits are those of a monochromatic mass function; the published bounds for extended mass functions, and the
population fits of GWTC-4.0 (Andrés-Carcasona et al. 2026, arXiv:2605.15749), are not reproduced here.
"""
from __future__ import annotations

import io
import os
import zipfile
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .report import write_simple_html_report

SSM_RECORD = "23346454"
SSM_URL = f"https://zenodo.org/records/{SSM_RECORD}/files/results-PBH.zip?download=1"
SSM_REF = "LVK 2026, arXiv:2610.20276 (O4a + O4b), data release Zenodo 23346454"
SSM_PIPELINES = {"GstLAL": "gstlal", "MBTA": "mbta", "PyCBC": "pycbc"}
SIGMA_M = 0.004                  # rescaled variance of the matter fluctuations when the binaries form (Hütsi+ 2021)
R_LATE_CLUSTERS = 420.0          # Gpc^-3 yr^-1, Phukon et al. 2021
CL = 0.9
SUPPRESSIONS = ("hutsi", "lvk")
CHANNELS = ("early", "late")


def _log(msg: str) -> None:
    print(f"[primordial_black_holes] {msg}", flush=True)


def cache_dir() -> Path:
    return Path(os.environ.get("GWTC_PBH_CACHE", Path.home() / ".cache_gwtc_analysis" / "primordial_black_holes"))


# ---------------------------------------------------------------------------
# merger rates of a monochromatic PBH population [Gpc^-3 yr^-1], m in Msun
# ---------------------------------------------------------------------------
def suppression_lvk(f) -> np.ndarray:
    """S = 2.3e-3 f^-0.65: the fit of the LVK sub-solar papers (arXiv:2605.05444, 2610.20276)."""
    return 2.3e-3 * np.asarray(f, dtype=float) ** -0.65


def suppression_hutsi(f) -> np.ndarray:
    """S = S1 * S2(z = 0) of Hütsi et al. 2021 for a monochromatic mass function (<m^2>/<m>^2 = 1, M = 2m)."""
    from scipy.special import gamma, hyperu

    f = np.asarray(f, dtype=float)
    nbar = 2.0 * f / (f + SIGMA_M)                                        # eq. 2.3
    u = gamma(29 / 37) / np.sqrt(np.pi) * hyperu(21 / 74, 0.5, 5 * f ** 2 / (6 * SIGMA_M ** 2))
    c = f ** 2 / SIGMA_M ** 2 / (u ** (-74 / 21) - 1.0)                   # eq. A.5
    s1 = 1.42 * (1.0 / (nbar + c) + SIGMA_M ** 2 / f ** 2) ** (-21 / 74) * np.exp(-nbar)   # eq. 2.2
    s2 = np.minimum(1.0, 9.6e-3 * f ** -0.65 * np.exp(0.03 * np.log(f) ** 2))             # eq. 2.7
    return s1 * s2


def rate_early(m, f, suppression: str = "hutsi") -> np.ndarray:
    """Merger rate today of early binaries of PBHs of mass m (equal masses: M = 2m, eta = 1/4)."""
    m, f = np.asarray(m, dtype=float), np.asarray(f, dtype=float)
    s = suppression_lvk(f) if suppression == "lvk" else suppression_hutsi(f)
    return 1.6e6 * f ** (53 / 37) * (2 * m) ** (-32 / 37) * 0.25 ** (-34 / 37) * s


def rate_late(m, f) -> np.ndarray:
    """Merger rate of late binaries formed in PBH clusters (Phukon et al. 2021, eq. 2)."""
    m, f = np.asarray(m, dtype=float), np.asarray(f, dtype=float)
    return f * R_LATE_CLUSTERS * (2 * m) ** (10 / 7) / (m * m) ** (5 / 7)


def f_pbh_limit(m: float, rate_limit: float, channel: str = "early", suppression: str = "hutsi") -> float:
    """The f_PBH at which the PBH merger rate equals `rate_limit` (values above 1 mean no constraint)."""
    if not np.isfinite(rate_limit) or rate_limit <= 0:
        return np.nan
    if channel == "late":
        return float(rate_limit / rate_late(m, 1.0))
    if suppression == "lvk":                       # R ∝ f^(53/37 - 0.65): closed form, as in the LVK script
        return float((rate_limit / rate_early(m, 1.0, "lvk")) ** (1.0 / (53 / 37 - 0.65)))
    from scipy.optimize import brentq

    g = lambda lnf: float(np.log(rate_early(m, np.exp(lnf), "hutsi"))) - np.log(rate_limit)
    lo, hi = np.log(1e-10), np.log(10.0)
    if g(hi) < 0:
        return np.inf
    if g(lo) > 0:
        return float(np.exp(lo))
    return float(np.exp(brentq(g, lo, hi, xtol=1e-6)))


def rate_upper_limit(n: int, vt: float, cl: float = CL) -> float:
    """Poisson upper limit on the rate for n events (uniform prior on the rate): Gamma^-1(cl; n + 1) / VT."""
    from scipy.stats import gamma

    return float(gamma.ppf(cl, n + 1) / vt) if vt > 0 else np.inf


# ---------------------------------------------------------------------------
# 0.2-1 Msun: the LVK sub-solar-mass searches
# ---------------------------------------------------------------------------
def fetch_ssm(refresh: bool = False) -> Path:
    """The PBH inputs of the O4b sub-solar data release (results-PBH.zip, 27 kB), downloaded once."""
    import requests

    from . import mirror

    p = cache_dir() / f"zenodo_{SSM_RECORD}_results-PBH.zip"
    if p.exists() and p.stat().st_size > 0 and not refresh:
        mirror.cached(SSM_URL)
        return p
    p.parent.mkdir(parents=True, exist_ok=True)

    def http(u, dest):
        r = requests.get(u, timeout=300)
        r.raise_for_status()
        Path(dest).write_bytes(r.content)

    mirror.fetch_to(SSM_URL, p, http, log=_log)
    _log(f"downloaded {SSM_URL}")
    return p


def read_ssm_vt(zip_path: Path) -> pd.DataFrame:
    """Chirp-mass-binned VT [Gpc^3 yr] of each pipeline, O4a + O4b combined, one row per bin and pipeline."""
    rows = []
    with zipfile.ZipFile(zip_path) as z:
        names = z.namelist()
        for pipe, key in SSM_PIPELINES.items():
            name = next((n for n in names if n.endswith(f"{key}_pbh_vt_inputs_mc.txt")), None)
            if name is None:
                raise ValueError(f"{zip_path.name}: no chirp-mass VT file for {pipe}")
            t = pd.read_csv(io.BytesIO(z.read(name)), sep=r"\s+")
            for _, r in t.iterrows():
                mc = 0.5 * (r["mc_min"] + r["mc_max"])
                rows.append(dict(pipeline=pipe, mc_min=r["mc_min"], mc_max=r["mc_max"], m_pbh=2 ** 0.2 * mc,
                                 vt_gpc3_yr=float(r["VT_O4ab"]), n_eff=float(r["neff_O4ab"])))
    return pd.DataFrame(rows)


def ssm_limits(vt: pd.DataFrame, suppression: str = "hutsi", m_max: float = 1.01) -> pd.DataFrame:
    """Limits from the sub-solar VT, for the bins with m <= m_max (equal masses inside the searches' range)."""
    out = vt[vt["m_pbh"] <= m_max].reset_index(drop=True)
    out["n_events"] = 0
    out["rate_90"] = [rate_upper_limit(0, v) for v in out["vt_gpc3_yr"]]
    out["f_early"] = [f_pbh_limit(m, r, "early", suppression) for m, r in zip(out["m_pbh"], out["rate_90"])]
    out["f_late"] = [f_pbh_limit(m, r, "late") for m, r in zip(out["m_pbh"], out["rate_90"])]
    out["source"] = "sub-solar searches O4a+O4b"
    return out


# ---------------------------------------------------------------------------
# above 1 Msun: the injections and the events of the rates mode
# ---------------------------------------------------------------------------
def _t_over_t0(z: np.ndarray) -> np.ndarray:
    from astropy.cosmology import Planck15

    zg = np.linspace(0.0, float(np.max(z)) * 1.01 + 1e-3, 2000)
    t = Planck15.age(zg).value
    return np.interp(z, zg, t / t[0])


def window_vt(inj: dict, m: float, w: float, ns_max_mass: float = 2.5, evolve: bool = True) -> tuple[float, float]:
    """(<VT> [Gpc^3 yr], effective sample size) of binaries with both masses log-uniform in [m e^-w, m e^w], m1 >= m2.

    With `evolve`, the early-binary rate evolution (t/t0)^(-34/37) weights the volume, so that VT * R(z = 0) is the
    expected number of detections; spins as in the rates mode (|a| < 0.4 below `ns_max_mass`, else 0.99)."""
    from .catalogs import _ln_iso_spin

    m1, m2, z = inj["mass1_source"], inj["mass2_source"], inj["redshift"]
    lo, hi = m * np.exp(-w), m * np.exp(w)
    inside = (m2 >= lo) & (m1 <= hi) & (m1 >= m2)
    with np.errstate(divide="ignore"):
        ln_mass = np.where(inside, np.log(2.0) - np.log(m1) - np.log(m2) - 2 * np.log(2 * w), -np.inf)
        ln_evol = (-34 / 37) * np.log(_t_over_t0(z)) if evolve else 0.0
        amax1, amax2 = np.where(m1 < ns_max_mass, 0.4, 0.99), np.where(m2 < ns_max_mass, 0.4, 0.99)
        lnp = (ln_mass + _ln_iso_spin(inj["a1"], amax1) + _ln_iso_spin(inj["a2"], amax2)
               + np.log(inj["dvdz"]) - np.log1p(z) + ln_evol)
    with np.errstate(over="ignore", invalid="ignore"):
        x = inj["weights"] * np.exp(lnp - inj["lnpdraw"])
    x = np.where(np.isfinite(x), x, 0.0)
    s = float(x.sum())
    return inj["T_yr"] * s / inj["N_gen"], (s * s / float((x * x).sum()) if s > 0 else 0.0)


def events_in_window(events: pd.DataFrame, m: float, w: float) -> pd.DataFrame:
    """Events whose median source-frame chirp mass is that of an equal-mass PBH binary in the window,
    2^(-1/5) [m e^-w, m e^w]. The chirp mass is the best-measured mass; the mass ratio is not used, so an unequal
    binary of the same chirp mass also counts (conservative for an upper limit). Without a published chirp mass,
    it is computed from the median component masses."""
    mc = pd.to_numeric(events.get("chirp_mass_source"), errors="coerce") if "chirp_mass_source" in events else None
    m1 = pd.to_numeric(events["mass_1_source"], errors="coerce")
    m2 = pd.to_numeric(events["mass_2_source"], errors="coerce")
    mc_comp = (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2
    mc = mc_comp if mc is None else mc.fillna(mc_comp)
    lo, hi = 2 ** -0.2 * m * np.exp(-w), 2 ** -0.2 * m * np.exp(w)
    return events[(mc >= lo) & (mc <= hi)]


def catalog_limits(inj: dict, events: pd.DataFrame, masses: np.ndarray, w: float, attribute: str = "window",
                   suppression: str = "hutsi", ns_max_mass: float = 2.5, min_neff: float = 20.0) -> pd.DataFrame:
    rows = []
    for m in masses:
        ev = events_in_window(events, m, w)
        n = len(ev) if attribute == "window" else 0
        vt_e, neff = window_vt(inj, m, w, ns_max_mass, evolve=True)
        vt_l, _ = window_vt(inj, m, w, ns_max_mass, evolve=False)
        ok = neff >= min_neff
        r_e, r_l = (rate_upper_limit(n, vt_e), rate_upper_limit(n, vt_l)) if ok else (np.nan, np.nan)
        rows.append(dict(pipeline="all (FAR threshold)", m_pbh=m, mass_min=m * np.exp(-w), mass_max=m * np.exp(w),
                         vt_gpc3_yr=vt_e, n_eff=neff, n_events=len(ev), events=" ".join(ev["event"].astype(str)),
                         rate_90=r_e, f_early=f_pbh_limit(m, r_e, "early", suppression) if ok else np.nan,
                         rate_90_late=r_l, f_late=f_pbh_limit(m, r_l, "late") if ok else np.nan,
                         source="GWTC injections and events"))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# plot and report
# ---------------------------------------------------------------------------
def plot(ssm: pd.DataFrame, cat: pd.DataFrame, out_png: Path, suppression: str) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
    blue, orange = "#2a78d6", "#eb6834"
    fig, ax = plt.subplots(figsize=(9.5, 4.8), dpi=150)
    fig.patch.set_facecolor(surface); ax.set_facecolor(surface)
    marks = {"GstLAL": "o", "MBTA": "s", "PyCBC": "^"}
    for pipe, g in (ssm.groupby("pipeline") if len(ssm) else ()):
        g = g.sort_values("m_pbh")
        ax.plot(g["m_pbh"], g["f_early"], color=blue, marker=marks.get(pipe, "o"), ms=4, lw=1.2,
                label=f"early binaries, sub-solar search ({pipe})")
        ax.plot(g["m_pbh"], g["f_late"], color=orange, marker=marks.get(pipe, "o"), ms=4, lw=1.0, ls="--",
                label=f"late binaries, sub-solar search ({pipe})")
    if len(cat):
        c = cat.sort_values("m_pbh")
        ax.plot(c["m_pbh"], c["f_early"], color=blue, lw=2.2, label="early binaries, GWTC injections and events")
        ax.plot(c["m_pbh"], c["f_late"], color=orange, lw=1.8, ls="--", label="late binaries, GWTC injections and events")
    ax.axhline(1.0, color=ink2, lw=0.8, ls=":")
    ax.text(0.21, 1.15, "f_PBH = 1: all of the dark matter", fontsize=8, color=ink2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.2, max(200.0, float(cat["m_pbh"].max()) * 1.2 if len(cat) else 1.0))
    ax.set_ylim(1e-5, 30)
    ax.set_xlabel("PBH mass m (M$_\\odot$), monochromatic", color=ink2)
    ax.set_ylabel("upper limit on f$_{\\rm PBH}$ (90%)", color=ink2)
    ax.grid(color=grid, lw=0.7, which="both", zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(grid)
    ax.tick_params(colors=ink2, labelsize=9)
    ax.legend(frameon=False, fontsize=7.5, labelcolor=ink2, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.set_title(f"Fraction of dark matter in primordial black holes (suppression: {suppression})", color=ink,
                 fontsize=10.5, loc="left")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def _at(df: pd.DataFrame, col: str, m: float) -> tuple[float, float]:
    """Range over pipelines of `col` at the bin whose mass is closest to m."""
    if not len(df):
        return np.nan, np.nan
    sel = df.loc[df.groupby("pipeline")["m_pbh"].transform(lambda s: (s - m).abs() == (s - m).abs().min()).astype(bool)]
    return float(sel[col].min()), float(sel[col].max())


def run_primordial_black_holes(
    masses: Optional[list[float]] = None,
    mass_window: float = 0.2,
    attribute: str = "window",
    suppression: str = "hutsi",
    far_threshold: float = 1.0,
    sensitivity_file: Optional[str | Path] = None,
    sensitivity_release: Optional[str] = None,
    sub_solar: bool = True,
    out_limits_tsv: Optional[str | Path] = "pbh_limits.tsv",
    out_report_html: Optional[str | Path] = "pbh_limits.html",
    plots_dir: str | Path = "pbh_plots",
    ns_max_mass: float = 2.5,
) -> pd.DataFrame:
    """Upper limits on f_PBH against the PBH mass, from the sub-solar searches and from the catalog."""
    from .catalogs import RATES_EVENT_LISTS, _load_found_injections, _rates_events, _rates_sensitivity_path
    from .cli_defaults import RATES_DEFAULT_RELEASE

    if suppression not in SUPPRESSIONS:
        raise ValueError(f"suppression must be one of {SUPPRESSIONS}")
    if attribute not in ("window", "none"):
        raise ValueError("attribute must be 'window' or 'none'")
    if not 0 < mass_window < 1:
        raise ValueError("--mass-window must be in (0, 1)")
    w = float(mass_window)
    if masses is None:
        masses = np.geomspace(np.exp(w) * 1.001, 150.0, 22)
    masses = np.asarray(sorted(masses), dtype=float)
    if masses.min() * np.exp(-w) < 1.0:
        raise ValueError(f"the catalog part needs m e^-w >= 1 Msun (the injections start at 1 Msun): "
                         f"m >= {np.exp(w):.3f} for --mass-window {w:g}")
    ssm = ssm_limits(read_ssm_vt(fetch_ssm()), suppression) if sub_solar else pd.DataFrame()
    release = sensitivity_release or RATES_DEFAULT_RELEASE
    inj = _load_found_injections(_rates_sensitivity_path(sensitivity_file, release), far_threshold)
    events = _rates_events(inj["segments"], far_threshold, ns_max_mass, RATES_EVENT_LISTS)
    _log(f"{len(inj['redshift'])} found injections ({release}), {len(events)} events with FAR <= {far_threshold:g}/yr")
    cat = catalog_limits(inj, events, masses, w, attribute, suppression, ns_max_mass)
    table = pd.concat([ssm, cat], ignore_index=True, sort=False)
    if len(ssm):
        e1, e035 = _at(ssm, "f_early", 1.0), _at(ssm, "f_early", 0.35)
        _log(f"sub-solar: early binaries f_PBH <= {e1[0]:.3g}-{e1[1]:.3g} at 1 Msun, {e035[0]:.3g}-{e035[1]:.3g} "
             f"at 0.35 Msun (suppression {suppression}; the LVK quote 0.018-0.032 and 0.10-0.18, with lvk)")
    best = cat.loc[cat["f_early"].idxmin()] if cat["f_early"].notna().any() else None
    if best is not None:
        _log(f"catalog: tightest early-binary limit f_PBH <= {best['f_early']:.2g} at {best['m_pbh']:.3g} Msun "
             f"({best['n_events']} events in the window)")
    if out_limits_tsv:
        Path(out_limits_tsv).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_limits_tsv, sep="\t", index=False, float_format="%.4g")
    if out_report_html:
        img = plot(ssm, cat, Path(plots_dir) / "f_pbh.png", suppression)
        paras = [
            "Upper limits (90%) on the fraction of dark matter in primordial black holes of a single mass m, from "
            "the merger rates: a PBH population of fraction f merges at a rate that grows with f, which may not "
            "exceed the rate the data allow.",
            f"<b>Early binaries</b> (formed in the radiation era): R = 1.6×10⁶ Gpc⁻³ yr⁻¹ f<sup>53/37</sup> "
            f"(t/t₀)<sup>−34/37</sup> (2m)<sup>−32/37</sup> η<sup>−34/37</sup> S(f) (Hütsi et al. 2021, "
            f"arXiv:2012.02786), with the suppression S = "
            + ("S₁ S₂ of Hütsi et al. 2021, valid down to small f." if suppression == "hutsi" else
               "2.3×10⁻³ f<sup>−0.65</sup>, the fit of the LVK sub-solar papers (equal to S₁ S₂ at f = 1, "
               "about 3 times smaller at f = 0.01, and above 1, unphysical, below f ~ 10⁻⁴).")
            + " <b>Late binaries</b> (captures in PBH clusters, Phukon et al. 2021): R = 420 Gpc⁻³ yr⁻¹ f "
              "(2m)<sup>10/7</sup> / (m²)<sup>5/7</sup>.",
        ]
        if len(ssm):
            paras.append(
                f"<b>0.2–1 M☉</b>: sensitive volume-time of the LVK sub-solar-mass searches, O4a and O4b combined, "
                f"per pipeline and chirp-mass bin ({SSM_REF}); no significant candidate, so R₉₀ = 2.30 / VT. "
                f"Early binaries: f ≤ {e1[0]:.3g}–{e1[1]:.3g} at 1 M☉ and {e035[0]:.3g}–{e035[1]:.3g} at 0.35 M☉ "
                f"(range over the pipelines; the LVK quote 0.018–0.032 and 0.10–0.18 with their fit of the suppression, reproduced with --suppression lvk).")
        att = ("every event whose median chirp mass is in the window, 2<sup>−1/5</sup> [m e<sup>−w</sup>, "
               "m e<sup>w</sup>], may be primordial (the mass ratio, poorly measured, is not used): R₉₀ is the 90% "
               "Poisson upper limit for those n events" if attribute == "window" else
               "the observed events are taken as astrophysical (n = 0): a tighter limit and a stronger assumption")
        paras.append(
            f"<b>Above 1 M☉</b>: the LVK search-sensitivity injections ({release}) reweighted to binaries with both "
            f"masses log-uniform in [m e<sup>−{w:g}</sup>, m e<sup>{w:g}</sup>], FAR ≤ {far_threshold:g}/yr; {att}. "
            "The early-binary volume includes the rate evolution (t/t₀)<sup>−34/37</sup>. "
            + (f"Tightest early-binary limit: f ≤ {best['f_early']:.2g} at {best['m_pbh']:.3g} M☉." if best is not None else ""))
        paras.append(
            "These are limits for a monochromatic mass function and fixed merger-rate models; the published bounds "
            "for extended mass functions and the population fits of GWTC-4.0 (Andrés-Carcasona et al. 2026, "
            "arXiv:2605.15749) are not reproduced. A sub-solar merger need not be primordial: collapsar disks may "
            "fragment into sub-solar neutron stars (Metzger et al. 2024), which would shine.")
        cols = ["source", "pipeline", "m_pbh", "vt_gpc3_yr", "n_eff", "n_events", "rate_90", "f_early", "f_late"]
        Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
        write_simple_html_report(out_report_html, title="Primordial black holes: limits on f_PBH", paragraphs=paras,
                                 images=[img], tables=[("Limits per mass", table[cols].to_html(
                                     index=False, float_format=lambda x: f"{x:.3g}"))])
        _log(f"report written to {out_report_html}")
    return table

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

These limits are those of a monochromatic mass function (equal-mass binaries). With `extended` (default), unequal
masses too:

- **the map F_PBH(m1, m2)** of the LVK script: in each (m1, m2) bin the rate is F dln m1 dln m2 R(m1, m2; F = 1),
  R = 1.6e6 * 2.3e-3 M^(-32/37) eta^(-34/37) (early, suppression fixed) or 420 M^(10/7) / (m1 m2)^(5/7) (late), from
  the (m1, m2) tables of the sub-solar release below m2 = 1.1 Msun and the injections in log bins above;
- **log-normal mass functions** (Hütsi et al. 2021, eq. 3.6, number of PBHs per ln m): the limit on the total f_PBH,
  the unequal-mass rates (eq. 2.1, with S1 for the mass function) summed against the VT of every pair, tested in
  each chirp-mass window against the events seen there (as the monochromatic limits) or pooled with
  `attribute="none"`.

The population fits of GWTC-4.0 (Andrés-Carcasona et al. 2026, arXiv:2605.15749) are not reproduced here.
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


def suppression_hutsi(f, m_over_mean=2.0, m2_ratio=1.0) -> np.ndarray:
    """S = S1 * S2(z = 0) of Hütsi et al. 2021.

    `m_over_mean` is M / <m>, the total mass of the binary over the mean PBH mass, and `m2_ratio` is <m^2>/<m>^2,
    the averages being over the number of PBHs: 2 and 1 for a monochromatic mass function (the defaults)."""
    from scipy.special import gamma, hyperu

    f = np.asarray(f, dtype=float)
    nbar = np.asarray(m_over_mean, dtype=float) * f / (f + SIGMA_M)       # eq. 2.3
    u = gamma(29 / 37) / np.sqrt(np.pi) * hyperu(21 / 74, 0.5, 5 * f ** 2 / (6 * SIGMA_M ** 2))
    c = f ** 2 * m2_ratio / SIGMA_M ** 2 / (u ** (-74 / 21) - 1.0)        # eq. A.5
    s1 = 1.42 * (m2_ratio / (nbar + c) + SIGMA_M ** 2 / f ** 2) ** (-21 / 74) * np.exp(-nbar)   # eq. 2.2
    return s1 * suppression_hutsi_s2(f)


def suppression_hutsi_s2(f) -> np.ndarray:
    """S2(z = 0) of Hütsi et al. 2021 (eq. 2.7), the disruption in PBH clusters: it depends on f only."""
    f = np.asarray(f, dtype=float)
    return np.minimum(1.0, 9.6e-3 * f ** -0.65 * np.exp(0.03 * np.log(f) ** 2))


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


def _injection_terms(inj: dict, ns_max_mass: float = 2.5, evolve: bool = True) -> np.ndarray:
    """Per found injection, weight * p(spins) dV/dz / (1 + z) [(t/t0)^(-34/37)] / p_draw: summed against a density
    in the source masses (m1 >= m2) and multiplied by T / N_gen, it gives VT [Gpc^3 yr] for a normalised density,
    or the expected number of detections for a rate density [Gpc^-3 yr^-1 Msun^-2]. Spins as in the rates mode
    (|a| < 0.4 below `ns_max_mass`, else 0.99)."""
    from .catalogs import _ln_iso_spin

    m1, m2, z = inj["mass1_source"], inj["mass2_source"], inj["redshift"]
    with np.errstate(divide="ignore"):
        ln_evol = (-34 / 37) * np.log(_t_over_t0(z)) if evolve else 0.0
        amax1, amax2 = np.where(m1 < ns_max_mass, 0.4, 0.99), np.where(m2 < ns_max_mass, 0.4, 0.99)
        lnp = (_ln_iso_spin(inj["a1"], amax1) + _ln_iso_spin(inj["a2"], amax2)
               + np.log(inj["dvdz"]) - np.log1p(z) + ln_evol)
    with np.errstate(over="ignore", invalid="ignore"):
        x = inj["weights"] * np.exp(lnp - inj["lnpdraw"])
    return np.where(np.isfinite(x), x, 0.0)


def _vt_of(inj: dict, x: np.ndarray, density: np.ndarray) -> tuple[float, float]:
    """(T / N_gen * sum x * density, effective sample size) for the injection terms x of `_injection_terms`."""
    y = x * density
    s = float(y.sum())
    return inj["T_yr"] * s / inj["N_gen"], (s * s / float((y * y).sum()) if s > 0 else 0.0)


def window_vt(inj: dict, m: float, w: float, ns_max_mass: float = 2.5, evolve: bool = True) -> tuple[float, float]:
    """(<VT> [Gpc^3 yr], effective sample size) of binaries with both masses log-uniform in [m e^-w, m e^w], m1 >= m2.

    With `evolve`, the early-binary rate evolution (t/t0)^(-34/37) weights the volume, so that VT * R(z = 0) is the
    expected number of detections."""
    return region_vt(inj, (m * np.exp(-w), m * np.exp(w)), (m * np.exp(-w), m * np.exp(w)),
                     _injection_terms(inj, ns_max_mass, evolve))


def region_vt(inj: dict, b1: tuple[float, float], b2: tuple[float, float], x: np.ndarray) -> tuple[float, float]:
    """(<VT>, effective sample size) of binaries log-uniform in m1 in b1 and m2 in b2 with m1 >= m2 (a bin of the
    (m1, m2) plane; a diagonal bin, b1 = b2, is the triangle above the diagonal)."""
    m1, m2 = inj["mass1_source"], inj["mass2_source"]
    inside = (m1 >= b1[0]) & (m1 <= b1[1]) & (m2 >= b2[0]) & (m2 <= b2[1]) & (m1 >= m2)
    area = ordered_ln_area(b1, b2)
    return _vt_of(inj, x, np.where(inside, 1.0 / (area * m1 * m2), 0.0)) if area > 0 else (0.0, 0.0)


def ordered_ln_area(b1: tuple[float, float], b2: tuple[float, float]) -> float:
    """Area in (ln m1, ln m2) of the part of the bin b1 x b2 with m1 >= m2."""
    l1, h1, l2, h2 = np.log([b1[0], b1[1], b2[0], b2[1]])
    u2 = l2 + (np.arange(4000) + 0.5) * (h2 - l2) / 4000                 # midpoints
    return float(np.clip(h1 - np.maximum(u2, l1), 0.0, h1 - l1).mean() * (h2 - l2))


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
# unequal masses: the map F_PBH(m1, m2) and log-normal mass functions
# ---------------------------------------------------------------------------
S_LVK_F1 = 2.3e-3                # the LVK suppression at f = 1, fixed in their F_PBH(m1, m2)
SSM_M2_MAX = 1.1                 # Msun, top of the secondary-mass bins of the sub-solar (m1, m2) tables
LOGNORMAL_WIDTHS = (0.5, 1.0, 2.0)


def rate_early_pair_lvk(m1, m2) -> np.ndarray:
    """dR / (d ln m1 d ln m2) of early binaries per unit F_PBH, as in the LVK script (rateEB_m1m2):
    1.6e6 * 2.3e-3 * M^(-32/37) eta^(-34/37), the suppression fixed at its value for f = 1."""
    m1, m2 = np.asarray(m1, dtype=float), np.asarray(m2, dtype=float)
    mt = m1 + m2
    return 1.6e6 * S_LVK_F1 * mt ** (-32 / 37) * (m1 * m2 / mt ** 2) ** (-34 / 37)


def rate_late_pair(m1, m2, f=1.0) -> np.ndarray:
    """dR / (d ln m1 d ln m2) of late binaries (Phukon et al. 2021, eq. 2; LVK rateLB_m1m2), linear in f:
    f * 420 * M^(10/7) / (m1 m2)^(5/7), per unit product of the two mass fractions."""
    m1, m2 = np.asarray(m1, dtype=float), np.asarray(m2, dtype=float)
    return f * R_LATE_CLUSTERS * (m1 + m2) ** (10 / 7) / (m1 * m2) ** (5 / 7)


def read_ssm_vt_m1m2(zip_path: Path) -> pd.DataFrame:
    """(m1, m2)-binned VT [Gpc^3 yr] of each pipeline, O4a + O4b combined; bins without injections get NaN."""
    rows = []
    with zipfile.ZipFile(zip_path) as z:
        names = z.namelist()
        for pipe, key in SSM_PIPELINES.items():
            name = next((n for n in names if n.endswith(f"{key}_pbh_vt_inputs_m1m2.txt")), None)
            if name is None:
                raise ValueError(f"{zip_path.name}: no (m1, m2) VT file for {pipe}")
            # read as the LVK script does: a row without O4a values is shorter, its last columns are NaN
            t = pd.read_csv(io.BytesIO(z.read(name)), sep=r"\s+", engine="python", na_values=["", "NaN", "nan"])
            for _, r in t.iterrows():
                vt = float(r["VT_O4ab"])
                rows.append(dict(pipeline=pipe, m1_min=r["m1_min"], m1_max=r["m1_max"], m2_min=r["m2_min"],
                                 m2_max=r["m2_max"], vt_gpc3_yr=vt if vt > 0 else np.nan,
                                 n_eff=float(r["neff_O4ab"]) if vt > 0 else np.nan))
    return pd.DataFrame(rows)


def ssm_map(vt2: pd.DataFrame) -> pd.DataFrame:
    """F_PBH(m1, m2) upper limits from the sub-solar (m1, m2) tables, exactly as the LVK script: R_90 = 2.30 / VT,
    the bin's rate is F * dln m1 * dln m2 * R(m1, m2; F = 1) at the arithmetic centre of the bin."""
    out = vt2.copy()
    out["m1"], out["m2"] = 0.5 * (out["m1_min"] + out["m1_max"]), 0.5 * (out["m2_min"] + out["m2_max"])
    dln = np.log(out["m1_max"] / out["m1_min"]) * np.log(out["m2_max"] / out["m2_min"])
    out["n_events"] = 0
    out["rate_90"] = -np.log(1 - CL) / out["vt_gpc3_yr"]
    out["rate_90_late"] = out["rate_90"]
    out["F_early"] = out["rate_90"] / (dln * rate_early_pair_lvk(out["m1"], out["m2"]))
    out["F_late"] = out["rate_90"] / (dln * rate_late_pair(out["m1"], out["m2"]))
    out["source"] = "sub-solar searches O4a+O4b"
    return out


def catalog_map(inj: dict, events: pd.DataFrame, edges: np.ndarray, attribute: str = "window",
                ns_max_mass: float = 2.5, min_neff: float = 20.0) -> pd.DataFrame:
    """F_PBH(m1, m2) upper limits from the GWTC injections and events, in log bins of both masses (m1 >= m2), with
    the definition of the LVK map (geometric centre of the bin). An event counts in the bin of its median masses."""
    xe, xl = _injection_terms(inj, ns_max_mass, True), _injection_terms(inj, ns_max_mass, False)
    a = pd.to_numeric(events["mass_1_source"], errors="coerce").to_numpy()
    b = pd.to_numeric(events["mass_2_source"], errors="coerce").to_numpy()
    e1, e2 = np.maximum(a, b), np.minimum(a, b)
    rows = []
    for j in range(len(edges) - 1):
        for i in range(j, len(edges) - 1):
            b1, b2 = (edges[i], edges[i + 1]), (edges[j], edges[j + 1])
            sel = (e1 >= b1[0]) & (e1 < b1[1]) & (e2 >= b2[0]) & (e2 < b2[1])
            n_ev = int(sel.sum())
            n = n_ev if attribute == "window" else 0
            vt_e, neff = region_vt(inj, b1, b2, xe)
            vt_l, _ = region_vt(inj, b1, b2, xl)
            ok = neff >= min_neff
            m1, m2 = np.sqrt(b1[0] * b1[1]), np.sqrt(b2[0] * b2[1])
            dln = np.log(b1[1] / b1[0]) * np.log(b2[1] / b2[0])
            r_e, r_l = (rate_upper_limit(n, vt_e), rate_upper_limit(n, vt_l)) if ok else (np.nan, np.nan)
            rows.append(dict(pipeline="all (FAR threshold)", m1_min=b1[0], m1_max=b1[1], m2_min=b2[0], m2_max=b2[1],
                             m1=m1, m2=m2, vt_gpc3_yr=vt_e, n_eff=neff, n_events=n_ev,
                             events=" ".join(events["event"].astype(str).to_numpy()[sel]),
                             rate_90=r_e, rate_90_late=r_l,
                             F_early=r_e / (dln * float(rate_early_pair_lvk(m1, m2))),
                             F_late=r_l / (dln * float(rate_late_pair(m1, m2))),
                             source="GWTC injections and events"))
    return pd.DataFrame(rows)


def lognormal_psi(m, mc: float, sigma: float) -> np.ndarray:
    """Log-normal mass function of Hütsi et al. 2021 (eq. 3.6): the distribution of the NUMBER of PBHs per ln m,
    normalised to 1; <m> = mc e^(sigma^2 / 2) and <m^2>/<m>^2 = e^(sigma^2). The mass density per ln m is
    m psi / <m>, a log-normal of the same width peaking at mc e^(sigma^2)."""
    lnm = np.log(np.asarray(m, dtype=float))
    return np.exp(-0.5 * ((lnm - np.log(mc)) / sigma) ** 2) / (np.sqrt(2 * np.pi) * sigma)


def grid_exposure(u1_edges, u2_edges, vt: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cells of a (ln m1, ln m2) grid with a constant VT: centres (m1, m2) and exposure 2 VT x (ln-area of the cell
    with m1 >= m2); a cell crossed by the diagonal counts for its part above it."""
    u1, u2 = np.asarray(u1_edges, float), np.asarray(u2_edges, float)
    c1, c2 = np.meshgrid(0.5 * (u1[1:] + u1[:-1]), 0.5 * (u2[1:] + u2[:-1]), indexing="ij")
    h1, h2 = np.meshgrid(0.5 * np.diff(u1), 0.5 * np.diff(u2), indexing="ij")
    s = (np.arange(16) + 0.5) / 16 * 2 - 1                                # 16 sub-rows in u2
    lo1 = c1[..., None] - h1[..., None]
    v2 = c2[..., None] + s * h2[..., None]
    frac = np.clip((c1[..., None] + h1[..., None] - np.maximum(v2, lo1)) / (2 * h1[..., None]), 0, 1).mean(-1)
    keep = frac > 0
    return np.exp(c1[keep]), np.exp(c2[keep]), (2.0 * vt * 4 * h1 * h2 * frac)[keep]


SUB_SOLAR_GROUP = -1             # the window of the sub-solar searches, which have no candidate


def chirp_window(mc, w: float) -> np.ndarray:
    """Index of the chirp-mass window of width 2w in ln Mc, aligned on that of an equal-mass 1 + 1 Msun binary."""
    return np.floor((np.log(np.asarray(mc, float)) - np.log(2 ** -0.2)) / (2 * w) + 0.5).astype(int)


class _ExposurePoints:
    """Points of the (m1, m2) plane with m1 >= m2, their exposure c and their counting window: the expected number
    of detections of a population is sum c * dR/(d ln m1 d ln m2), dR being the full-plane rate density of Hütsi
    et al. (eq. 2.1), whose integral over all (m1, m2) is the total rate (a pair m1 > m2 appears twice, hence the
    factor 2 here).

    Below m2 = 1.1 Msun: the cells of a fine grid inside the sub-solar (m1, m2) bins, with the bin's VT (no rate
    evolution: the sub-solar horizon is small), all in one window (no candidate); above: the GWTC injections, in
    chirp-mass windows of width 2w."""

    def __init__(self, m1, m2, c_early, c_late, group=None):
        self.m1, self.m2 = np.asarray(m1, float), np.asarray(m2, float)
        self.c_early, self.c_late = np.asarray(c_early, float), np.asarray(c_late, float)
        self.group = np.full(len(self.m1), SUB_SOLAR_GROUP) if group is None else np.asarray(group, int)

    @classmethod
    def build(cls, inj: Optional[dict], vt2: Optional[pd.DataFrame], ns_max_mass: float = 2.5, grid: int = 16,
              w: float = 0.2):
        parts = []
        if vt2 is not None and len(vt2):
            for _, r in vt2.dropna(subset=["vt_gpc3_yr"]).iterrows():
                u1 = np.linspace(np.log(r["m1_min"]), np.log(r["m1_max"]), grid + 1)
                u2 = np.linspace(np.log(r["m2_min"]), np.log(r["m2_max"]), grid + 1)
                m1, m2, c = grid_exposure(u1, u2, r["vt_gpc3_yr"])
                parts.append((m1, m2, c, c, np.full(len(m1), SUB_SOLAR_GROUP)))
        if inj is not None:
            m1, m2 = inj["mass1_source"], inj["mass2_source"]
            keep = (m1 >= m2) & ((m2 >= SSM_M2_MAX) if vt2 is not None and len(vt2) else True)
            scale = inj["T_yr"] / inj["N_gen"] * 2.0 / (m1 * m2)       # per (ln m1, ln m2), both orderings
            xe, xl = _injection_terms(inj, ns_max_mass, True), _injection_terms(inj, ns_max_mass, False)
            grp = chirp_window((m1 * m2) ** 0.6 / (m1 + m2) ** 0.2, w)
            parts.append((m1[keep], m2[keep], (scale * xe)[keep], (scale * xl)[keep], grp[keep]))
        if not parts:
            return cls([], [], [], [])
        return cls(*(np.concatenate([p[k] for p in parts]) for k in range(5)))


F_GRID = np.geomspace(1e-10, 10.0, 67)          # 6 points per decade: scan of the expected counts against f


def _early_counts(a: np.ndarray, mo: np.ndarray, m2_ratio: float, f: float, suppression: str) -> float:
    """Expected detections of early binaries, f^(53/37) sum a S(f, M/<m>)."""
    if suppression == "lvk":
        return float(f ** (53 / 37) * S_LVK_F1 * f ** -0.65 * a.sum())
    return float(f ** (53 / 37) * (a * suppression_hutsi(f, mo, m2_ratio)).sum())


def _first_excluded(windows: list, count, f_grid: np.ndarray = F_GRID) -> tuple[float, float]:
    """(f_limit, f_reallowed) for windows [(data, n90), ...] and count(data, f): f_limit is the smallest f at which
    some window expects more than its n90 (refined between grid points), f_reallowed the first grid f above it at
    which no window does any more (NaN if none): for a broad mass function the suppression S1 can make the early-
    binary count fall with f, the heavy pairs being disrupted, so that larger f may be allowed again."""
    from scipy.optimize import brentq

    if not windows:
        return np.nan, np.nan
    ratio = np.array([[count(d, f) / n90 for f in f_grid] for d, n90 in windows])     # (window, f)
    excl = (ratio >= 1).any(axis=0)
    if not excl.any():
        return np.inf, np.nan
    i = int(np.argmax(excl))
    if i == 0:
        f_lim = float(f_grid[0])
    else:
        lo, hi = np.log(f_grid[i - 1]), np.log(f_grid[i])
        roots = [brentq(lambda u: np.log(count(d, np.exp(u)) / n90), lo, hi, xtol=1e-6)
                 for (d, n90), r in zip(windows, ratio) if r[i] >= 1 and r[i - 1] < 1]
        f_lim = float(np.exp(min(roots)))
    after = np.where(~excl[i:])[0]
    return f_lim, (float(f_grid[i + after[0]]) if len(after) else np.nan)


def lognormal_limit(points: _ExposurePoints, mc: float, sigma: float, n_by_group: Optional[dict] = None,
                    suppression: str = "hutsi", pooled: bool = False, details: bool = False):
    """(f_early, f_late): the smallest f_PBH of a log-normal population (`lognormal_psi`) at which, in some counting
    window, the expected number of detections exceeds the 90% Poisson upper limit for the n events seen there,
    Gamma^-1(0.9; n + 1) (2.30 for none). With `pooled`, all the windows form one, with n = 0. With `details`, also
    f_early_reallowed (see `_first_excluded`).

    Early binaries: Hütsi et al. eq. 2.1, dR/(d ln m1 d ln m2) = 1.6e6 f^(53/37) M^(-32/37) eta^(-34/37) S
    psi(m1) psi(m2) m1 m2 / <m>^2, with S1 S2 for the mass function (eqs. 2.2, 2.3, A.5: M/<m> and <m^2>/<m>^2)
    or the LVK fit 2.3e-3 f^-0.65. Late binaries: f * 420 M^(10/7) / (m1 m2)^(5/7) times the product of the mass
    fractions per ln m, m psi(m) / <m> (the LVK rateLB_m1m2), linear in f."""
    from scipy.stats import gamma as gamma_dist

    if not len(points.m1):
        return (np.nan, np.nan, np.nan) if details else (np.nan, np.nan)
    n_by_group = {} if pooled or n_by_group is None else n_by_group
    m1, m2 = points.m1, points.m2
    mean = mc * np.exp(sigma ** 2 / 2)
    m2_ratio = np.exp(sigma ** 2)
    p1, p2 = lognormal_psi(m1, mc, sigma), lognormal_psi(m2, mc, sigma)
    mt = m1 + m2
    a_early = points.c_early * 1.6e6 * mt ** (-32 / 37) * (m1 * m2 / mt ** 2) ** (-34 / 37) * p1 * p2 * m1 * m2 / mean ** 2
    a_late = points.c_late * rate_late_pair(m1, m2) * (m1 * p1 / mean) * (m2 * p2 / mean)
    keep = (a_early > 1e-12 * a_early.max()) | (a_late > 1e-12 * a_late.max())
    groups = np.zeros(len(m1), int) if pooled else points.group
    early, late = [], []
    for k in np.unique(groups[keep]):
        sel = keep & (groups == k)
        n90 = float(gamma_dist.ppf(CL, n_by_group.get(int(k), 0) + 1))
        if a_late[sel].sum() > 0:
            late.append(n90 / float(a_late[sel].sum()))
        if a_early[sel].sum() > 0:
            early.append(((a_early[sel], (mt / mean)[sel]), n90))
    f_e, f_re = _first_excluded(early, lambda d, f: _early_counts(d[0], d[1], m2_ratio, f, suppression))
    f_l = min(late) if late else np.nan
    return (f_e, f_l, f_re) if details else (f_e, f_l)


def lognormal_limits(inj: Optional[dict], vt2: Optional[pd.DataFrame], events: pd.DataFrame, mcs, sigmas,
                     attribute: str = "window", suppression: str = "hutsi", ns_max_mass: float = 2.5,
                     w: float = 0.2) -> pd.DataFrame:
    """Upper limits on the total f_PBH of a log-normal mass function, per pipeline of the sub-solar searches (their
    (m1, m2) tables below m2 = 1.1 Msun, combined with the GWTC injections above). With `attribute="window"` the
    population is tested in each chirp-mass window of width 2w against the events seen there, as the monochromatic
    limits are, and the tightest window sets the limit; with "none", all the detections are pooled with n = 0."""
    mc_ev = pd.to_numeric(events.get("chirp_mass_source"), errors="coerce") if "chirp_mass_source" in events else None
    e1 = pd.to_numeric(events["mass_1_source"], errors="coerce")
    e2 = pd.to_numeric(events["mass_2_source"], errors="coerce")
    mc_comp = (e1 * e2) ** 0.6 / (e1 + e2) ** 0.2
    mc_ev = (mc_comp if mc_ev is None else mc_ev.fillna(mc_comp)).dropna()
    n_by_group = pd.Series(chirp_window(mc_ev.to_numpy(), w)).value_counts().to_dict() if len(mc_ev) else {}
    pipes = list(vt2["pipeline"].unique()) if vt2 is not None and len(vt2) else [None]
    rows = []
    for pipe in pipes:
        sub = vt2[vt2["pipeline"] == pipe] if pipe is not None else None
        pts = _ExposurePoints.build(inj, sub, ns_max_mass, w=w)
        for sigma in sigmas:
            for mc in mcs:
                f_e, f_l, f_re = lognormal_limit(pts, mc, sigma, n_by_group, suppression,
                                                 pooled=attribute == "none", details=True)
                rows.append(dict(pipeline=pipe or "all (FAR threshold)", m_c=mc, sigma=sigma,
                                 mean_m=mc * np.exp(sigma ** 2 / 2), peak_mass_density=mc * np.exp(sigma ** 2),
                                 f_early=f_e, f_early_reallowed=f_re, f_late=f_l))
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
    ax.axhspan(1.0, 30.0, color=grid, alpha=0.6, lw=0, zorder=0.5)
    ax.axhline(1.0, color=ink2, lw=0.8, ls=":")
    ax.text(0.21, 0.78, "f_PBH = 1: all of the dark matter", fontsize=8, color=ink2, va="top")
    ax.text(0.21, 12.0, "no constraint (f_PBH > 1)", fontsize=8, color=ink2, style="italic", va="center")
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


def _style(ax, ink2, grid):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(grid)
    ax.tick_params(colors=ink2, labelsize=9)


def plot_map(maps: pd.DataFrame, out_png: Path) -> Path:
    """F_PBH(m1, m2) upper limits, early and late binaries: per bin, the weakest of the sub-solar pipelines."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Rectangle

    surface, ink, ink2, grid = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
    keys = ["source", "m1_min", "m1_max", "m2_min", "m2_max"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), dpi=150)
    fig.patch.set_facecolor(surface)
    for ax, col, name in zip(axes, ("F_early", "F_late"), ("early binaries", "late binaries")):
        ax.set_facecolor(surface)
        b = maps.groupby(keys, as_index=False)[col].max(min_count=1).dropna(subset=[col])
        b = b[b[col] > 0]
        if len(b):
            norm = LogNorm(vmin=float(b[col].min()), vmax=float(b[col].max()))
            cmap = plt.get_cmap("viridis")
            for _, r in b.iterrows():
                ax.add_patch(Rectangle((r["m1_min"], r["m2_min"]), r["m1_max"] - r["m1_min"], r["m2_max"] - r["m2_min"],
                                       facecolor=cmap(norm(r[col])), edgecolor=surface, lw=0.4))
            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            cb = fig.colorbar(sm, ax=ax, pad=0.02)
            cb.set_label("upper limit on F$_{\\rm PBH}$(m$_1$, m$_2$) (90%)", color=ink2, fontsize=8.5)
            cb.ax.tick_params(colors=ink2, labelsize=8)
            hi = float(max(b["m1_max"].max(), b["m2_max"].max()))
        else:
            hi = 200.0
        ax.plot([0.19, hi], [0.19, hi], color=ink2, lw=0.7, ls=":")
        ax.axhline(SSM_M2_MAX, color=ink2, lw=0.6, ls="--")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(0.19, hi); ax.set_ylim(0.19, hi)
        ax.set_xlabel("m$_1$ (M$_\\odot$)", color=ink2); ax.set_ylabel("m$_2$ (M$_\\odot$)", color=ink2)
        ax.set_title(name, color=ink, fontsize=10, loc="left")
        _style(ax, ink2, grid)
    fig.suptitle("Unequal masses: F$_{\\rm PBH}$(m$_1$, m$_2$), as defined by the LVK (below the dashed line: "
                 "sub-solar searches, weakest pipeline; above: GWTC)", color=ink, fontsize=10, x=0.01, ha="left")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def plot_lognormal(ln: pd.DataFrame, mono: pd.DataFrame, out_png: Path, suppression: str) -> Path:
    """Upper limit on the total f_PBH of a log-normal mass function against m_c, per width; the band spans the
    sub-solar pipelines. The monochromatic early-binary limit for reference."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
    colors = ["#2a78d6", "#1b9e77", "#8e44ad", "#c0392b", "#d68910"]
    fig, ax = plt.subplots(figsize=(9.5, 4.8), dpi=150)
    fig.patch.set_facecolor(surface); ax.set_facecolor(surface)
    if len(mono):
        g = mono.groupby("m_pbh", as_index=False)["f_early"].max()
        ax.plot(g["m_pbh"], g["f_early"], color="#b9b8b4", lw=1.6, label="monochromatic, early binaries (weakest)")
    for k, (sigma, g) in enumerate(ln.groupby("sigma")):
        c = colors[k % len(colors)]
        for col, ls, name in (("f_early", "-", "early"), ("f_late", "--", "late")):
            r = g.groupby("m_c")[col].agg(["min", "max"]).replace(np.inf, np.nan)
            ax.plot(r.index, r["max"], color=c, lw=1.8 if ls == "-" else 1.3, ls=ls,
                    label=f"σ = {sigma:g}, {name} binaries")
            if (r["max"] > r["min"] * 1.001).any():
                ax.fill_between(r.index, r["min"], r["max"], color=c, alpha=0.15, lw=0)
        isl = g[g["f_early_reallowed"].notna()] if "f_early_reallowed" in g else g.iloc[:0]
        if len(isl):
            ax.plot(isl["m_c"], isl["f_early"], ls="none", marker="o", mfc=surface, mec=c, ms=6, mew=1.4,
                    label=f"σ = {sigma:g}: allowed again at larger f")
    ax.axhspan(1.0, 30.0, color=grid, alpha=0.6, lw=0, zorder=0.5)
    ax.axhline(1.0, color=ink2, lw=0.8, ls=":")
    ax.text(0.21, 12.0, "no constraint (f_PBH > 1)", fontsize=8, color=ink2, style="italic", va="center")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.2, 200.0); ax.set_ylim(1e-5, 30)
    ax.set_xlabel("m$_c$ (M$_\\odot$): peak of the log-normal number distribution", color=ink2)
    ax.set_ylabel("upper limit on total f$_{\\rm PBH}$ (90%)", color=ink2)
    ax.grid(color=grid, lw=0.7, which="both", zorder=0)
    _style(ax, ink2, grid)
    ax.legend(frameon=False, fontsize=7.5, labelcolor=ink2, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.set_title(f"Log-normal mass functions (suppression: {suppression})", color=ink, fontsize=10.5, loc="left")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=surface)
    plt.close(fig)
    return out_png


def extra_table_path(out_limits_tsv: str | Path, name: str) -> Path:
    """`pbh_limits.tsv` -> `pbh_limits.<name>.tsv`: the extra tables sit next to the limits."""
    p = Path(out_limits_tsv)
    return p.with_name(f"{p.stem}.{name}{p.suffix or '.tsv'}")


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
    extended: bool = True,
    lognormal_widths: Optional[list[float]] = None,
    lognormal_masses: Optional[list[float]] = None,
    map_edges: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Upper limits on f_PBH against the PBH mass, from the sub-solar searches and from the catalog; with
    `extended`, also the unequal-mass map F_PBH(m1, m2) and the limits for log-normal mass functions."""
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
    widths = [float(s) for s in (lognormal_widths if lognormal_widths else LOGNORMAL_WIDTHS)]
    if extended and min(widths) <= 0:
        raise ValueError("--lognormal-widths must be > 0")
    mcs = np.asarray(sorted(lognormal_masses), float) if lognormal_masses else np.geomspace(0.3, 150.0, 24)
    if extended and mcs.min() <= 0:
        raise ValueError("--lognormal-masses must be > 0")
    zip_path = fetch_ssm() if sub_solar else None
    ssm = ssm_limits(read_ssm_vt(zip_path), suppression) if sub_solar else pd.DataFrame()
    release = sensitivity_release or RATES_DEFAULT_RELEASE
    inj = _load_found_injections(_rates_sensitivity_path(sensitivity_file, release), far_threshold)
    events = _rates_events(inj["segments"], far_threshold, ns_max_mass, RATES_EVENT_LISTS)
    _log(f"{len(inj['redshift'])} found injections ({release}), {len(events)} events with FAR <= {far_threshold:g}/yr")
    cat = catalog_limits(inj, events, masses, w, attribute, suppression, ns_max_mass)
    table = pd.concat([ssm, cat], ignore_index=True, sort=False)
    maps, ln = pd.DataFrame(), pd.DataFrame()
    if extended:
        vt2 = read_ssm_vt_m1m2(zip_path) if sub_solar else None
        edges = np.asarray(map_edges, float) if map_edges is not None else np.geomspace(SSM_M2_MAX, 160.0, 12)
        maps = pd.concat([ssm_map(vt2) if vt2 is not None else pd.DataFrame(),
                          catalog_map(inj, events, edges, attribute, ns_max_mass)], ignore_index=True, sort=False)
        ln = lognormal_limits(inj, vt2, events, mcs, widths, attribute, suppression, ns_max_mass, w)
        for sigma, g in ln.groupby("sigma"):
            ok = g[np.isfinite(g["f_early"])]
            if len(ok):
                b = ok.loc[ok["f_early"].idxmin()]
                _log(f"log-normal sigma = {sigma:g}: tightest early-binary limit f_PBH <= {b['f_early']:.2g} "
                     f"at m_c = {b['m_c']:.3g} Msun")
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
        if extended:
            maps.to_csv(extra_table_path(out_limits_tsv, "map"), sep="\t", index=False, float_format="%.4g")
            ln.to_csv(extra_table_path(out_limits_tsv, "lognormal"), sep="\t", index=False, float_format="%.4g")
    if out_report_html:
        img = plot(ssm, cat, Path(plots_dir) / "f_pbh.png", suppression)
        imgs = [img]
        if extended:
            imgs += [plot_map(maps, Path(plots_dir) / "f_pbh_map.png"),
                     plot_lognormal(ln, table, Path(plots_dir) / "f_pbh_lognormal.png", suppression)]
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
            "The limits above are for a monochromatic mass function (all PBHs of mass m, equal-mass binaries) and "
            "fixed merger-rate models; the population fits of GWTC-4.0 (Andrés-Carcasona et al. 2026, "
            "arXiv:2605.15749) are not reproduced. A sub-solar merger need not be primordial: collapsar disks may "
            "fragment into sub-solar neutron stars (Metzger et al. 2024), which would shine.")
        if extended:
            paras.append(
                "<b>Unequal masses, the map</b>: one limit per (m₁, m₂) bin on F<sub>PBH</sub>, defined as in the "
                "LVK script: the rate in the bin is F Δln m₁ Δln m₂ R(m₁, m₂; F = 1), with "
                "R = 1.6×10⁶ × 2.3×10⁻³ M<sup>−32/37</sup> η<sup>−34/37</sup> for early binaries (the suppression "
                "fixed at its value for f = 1, whatever --suppression) and 420 M<sup>10/7</sup> / "
                "(m₁m₂)<sup>5/7</sup> for late binaries. Each bin is tested alone: F is not a fraction of the dark "
                "matter but, for a mass function ψ, about f<sup>53/37</sup> ψ(m₁) ψ(m₂) m₁m₂/⟨m⟩² (early binaries). "
                "Below m₂ = 1.1 M☉: the sub-solar (m₁, m₂) tables (the plot shows the weakest pipeline); above: "
                "the GWTC injections in log bins, the events counted in the bin of their median masses.")
            paras.append(
                "<b>Unequal masses, log-normal mass functions</b>: the total f<sub>PBH</sub> at which a log-normal "
                "population, ψ(m) ∝ exp(−ln²(m/m<sub>c</sub>)/2σ²) for the number of PBHs per ln m (Hütsi et al. "
                "2021, eq. 3.6), gives the 90% upper limit on the number of detections, summed over all (m₁, m₂): "
                "the sub-solar (m₁, m₂) tables below m₂ = 1.1 M☉ (one limit per pipeline, the band of the plot) and "
                "the GWTC injections above. Early binaries: Hütsi et al. eq. 2.1, with S₁ for the mass function "
                "(M/⟨m⟩ and ⟨m²⟩/⟨m⟩², eqs. 2.2, 2.3, A.5) or the LVK fit; late binaries: the rate above times the "
                "mass fractions m ψ(m)/⟨m⟩ of the two black holes. Pairs with m₂ < 1.1 M☉ and m₁ > 11 M☉, "
                "outside the sub-solar tables, are left out (conservative). "
                + (f"As for the monochromatic limits, the population is tested in each chirp-mass window of width "
                   f"2w = {2 * w:g} in ln M<sub>c</sub> against the events seen there (every one may be primordial), "
                   "the sub-solar searches forming one more window with no candidate; the tightest window sets the "
                   "limit, which tends to the monochromatic one as σ → 0." if attribute == "window" else
                   "All the detections are pooled and the observed events taken as astrophysical (n = 0)."))
        cols = ["source", "pipeline", "m_pbh", "vt_gpc3_yr", "n_eff", "n_events", "rate_90", "f_early", "f_late"]
        tables = [("Limits per mass", table[cols].to_html(index=False, float_format=lambda x: f"{x:.3g}"))]
        if extended and len(ln):
            lt = (ln.groupby(["sigma", "m_c"], as_index=False)
                  .agg(f_early_min=("f_early", "min"), f_early_max=("f_early", "max"),
                       f_late_min=("f_late", "min"), f_late_max=("f_late", "max")))
            tables.append(("Log-normal mass functions (range over the sub-solar pipelines)",
                           lt.to_html(index=False, float_format=lambda x: f"{x:.3g}")))
        Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
        write_simple_html_report(out_report_html, title="Primordial black holes: limits on f_PBH", paragraphs=paras,
                                 images=imgs, tables=tables)
        _log(f"report written to {out_report_html}")
    return table

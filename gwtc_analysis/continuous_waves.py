"""Continuous waves: the published LVK upper limits, and the boson clouds of the scalar_cloud mode.

No continuous wave has been detected. The mode gathers the upper limits of the LVK data releases (LIGO DCC, public)
and shows them as the papers do:

- **known pulsars, O4a** (LVK 2025, ApJ 983, 99, arXiv:2501.01495; data LIGO-G2402400): the 95% upper limits on
  the strain amplitude h0 of 45 pulsars (Bayesian pipeline), against their spin-down limits; the spin-down limits,
  distances and ellipticity limits come from the paper's table (its LaTeX source). The narrowband search of 16 of
  them (and 2 post-glitch analyses), which lets the frequency wander around twice the rotation frequency, from the
  same release: less sensitive than the targeted one, robust to a mismatch between the GW and EM timings;
- **all-sky, O4a** (LVK 2026, arXiv:2603.14168; data LIGO-P2500416): the sensitivity to unknown isolated neutron
  stars, h0 against frequency, for FrequencyHough, PowerFlux (population average) and SOAP;
- **boson clouds, O3** (LVK 2022, PRD 105, 102001, arXiv:2111.15507; data LIGO-P2100343): the strain upper limits
  of the all-sky search for the annihilation signal of scalar clouds around spinning black holes (20-610 Hz), the
  region of black-hole mass and boson mass it excludes (for a distance, an initial spin and an age of the cloud),
  and the largest distance at which clouds are excluded, against the boson mass.

A cloud of bosons of mass m radiates at f = 2 m c² / h (483.6 Hz for 10⁻¹² eV) while it annihilates; the
`--boson-mass` band (by default the 90% interval of GW190728 with the superradiance prior, from scalar_cloud) is
drawn on the boson plots, to compare the inspiral dephasing with the continuous-wave limits on the same particle.
"""
from __future__ import annotations

import io
import os
import re
import zipfile
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

H_PLANCK_EV_S = 4.135667696e-15
EV_TO_HZ = 2.0 / H_PLANCK_EV_S          # GW frequency of the annihilation signal per eV of boson mass (2 m c^2 / h)
# GW190728 with the superradiance prior (tau_d = 1e6 yr): 5%, 50%, 95% of m_phi, per-sample source-frame mass
GW190728_BOSON_MASS_EV = (0.70e-12, 0.89e-12, 1.17e-12)

DCC = "https://dcc.ligo.org/public"
FILES = {
    "known_pulsars": f"{DCC}/0197/G2402400/001/bayes_ul_O4a.txt",
    "narrowband": f"{DCC}/0197/G2402400/001/narrowband_ul_O4a.txt",
    "allsky_fh": f"{DCC}/0201/P2500416/006/f_vs_h0min_o4a_FH.txt",
    "allsky_powerflux": f"{DCC}/0201/P2500416/006/uls_powerflux_popave.txt",
    "allsky_soap": f"{DCC}/0201/P2500416/006/uls_soap_noheader.txt",
    "boson_o3": f"{DCC}/0178/P2100343/008/O3_boson_cloud_data_release.zip",
    "known_pulsars_tex": "https://arxiv.org/e-print/2501.01495",
}
REFS = {
    "known_pulsars": "LVK 2025, ApJ 983, 99 (arXiv:2501.01495), data LIGO-G2402400",
    "allsky": "LVK 2026, arXiv:2603.14168, data LIGO-P2500416",
    "boson": "LVK 2022, PRD 105, 102001 (arXiv:2111.15507), data LIGO-P2100343",
}


def _log(msg: str) -> None:
    print(f"[continuous_waves] {msg}", flush=True)


def cache_dir() -> Path:
    return Path(os.environ.get("GWTC_CW_CACHE", Path.home() / ".cache_gwtc_analysis" / "continuous_waves"))


def fetch(name: str, refresh: bool = False) -> Path:
    """A file of the data releases, downloaded once."""
    import requests

    url = FILES[name]
    p = cache_dir() / (url.rsplit("/", 1)[-1] if name != "known_pulsars_tex" else "2501.01495.tar.gz")
    if p.exists() and not refresh:
        return p
    p.parent.mkdir(parents=True, exist_ok=True)
    r = requests.get(url, timeout=300)
    r.raise_for_status()
    p.with_suffix(".part").write_bytes(r.content)
    p.with_suffix(".part").replace(p)
    _log(f"downloaded {url}")
    return p


def boson_frequency_hz(m_ev) -> np.ndarray:
    """GW frequency of a cloud of bosons of mass m (eV): f = 2 m c^2 / h."""
    return np.asarray(m_ev, float) * EV_TO_HZ


# ---------------------------------------------------------------------------
# known pulsars
# ---------------------------------------------------------------------------
def _scinum(s: str) -> float:
    """\\scinum{a}{b} or \\scinum{a(c)}{b} -> a x 10^b (the first value); plain numbers as they are."""
    m = re.search(r"\\scinum\{\s*([-\d.]*)(?:\([^)]*\))?\s*\}\{\s*([-\d]+)\s*\}", s)
    if m:
        return float(m.group(1)) * 10.0 ** int(m.group(2)) if m.group(1) else np.nan
    m = re.search(r"-?\d+(?:\.\d*)?", s.replace("\\textminus", "-"))
    return float(m.group(0)) if m else np.nan


def _clean_name(s: str) -> str:
    s = re.sub(r"\\tablenotemark\{[^}]*\}", "", s).replace("\\textminus", "-").replace("$-$", "-")
    m = re.search(r"J\d{4}[+-]\d{2,4}[A-Za-z]?", s)
    return m.group(0) if m else s.strip()


def _known_pulsars_table(tex_tar: Path) -> pd.DataFrame:
    """The Bayesian rows of the results table of the paper: name, f_rot, distance, spin-down limit, h0, ε, ratio."""
    import tarfile

    with tarfile.open(tex_tar) as tf:
        text = tf.extractfile("Bayes-table.tex").read().decode("utf-8", "replace")
    rows = []
    for line in text.splitlines():
        if not line.lstrip().startswith("\\multirow{3}") or "Bayesian" not in line:
            continue
        c = [x.strip() for x in line.split("&")]
        mr = lambda s: re.sub(r"^\\multirow\{3\}\{\*\}\{(.*)\}\s*$", r"\1", s)
        rows.append(dict(pulsar=_clean_name(mr(c[0])), f_rot_hz=_scinum(mr(c[1])), distance_kpc=_scinum(mr(c[3])),
                         h0_spindown=_scinum(mr(c[4])), h0_ul_paper=_scinum(c[6]), ellipticity_ul=_scinum(c[7]),
                         ratio_to_spindown=_scinum(c[9])))
    return pd.DataFrame(rows)


def known_pulsars(refresh: bool = False) -> pd.DataFrame:
    """The 45 pulsars of the O4a targeted search: the DCC upper limits (Bayesian) with the paper's table."""
    d = pd.read_csv(fetch("known_pulsars", refresh), sep=r"\s+", header=None, names=["pulsar", "f_rot_hz", "h0_ul"])
    d["f_gw_hz"] = 2.0 * d["f_rot_hz"]
    try:
        t = _known_pulsars_table(fetch("known_pulsars_tex", refresh))
        d = d.merge(t.drop(columns=["f_rot_hz"]), on="pulsar", how="left")
    except Exception as e:                       # the table adds the spin-down limits; the plot works without
        _log(f"WARN: paper table not read ({type(e).__name__}: {e}); no spin-down limits")
    return d


def narrowband(kp: Optional[pd.DataFrame] = None, refresh: bool = False) -> pd.DataFrame:
    """The narrowband upper limits of O4a: pulsar, post-glitch analysis or not, GW frequency, h0, and the ratio to
    the spin-down limit when `kp` (known_pulsars) has it."""
    d = pd.read_csv(fetch("narrowband", refresh), sep=r"\s+", header=None, names=["label", "f_gw_hz", "h0_ul"])
    d["post_glitch"] = d["label"].str.endswith("pg")
    d["pulsar"] = d["label"].str.replace(r"pg$", "", regex=True)
    if kp is not None and "h0_spindown" in kp:
        d = d.merge(kp[["pulsar", "h0_spindown"]], on="pulsar", how="left")
        d["ratio_to_spindown"] = d["h0_ul"] / d["h0_spindown"]
    return d.drop(columns=["label"])


# ---------------------------------------------------------------------------
# all-sky and boson clouds
# ---------------------------------------------------------------------------
def allsky_o4a(refresh: bool = False) -> dict[str, pd.DataFrame]:
    out = {}
    fh = np.loadtxt(fetch("allsky_fh", refresh))
    out["FrequencyHough"] = pd.DataFrame(dict(f_hz=fh[:, 0], h0_ul=fh[:, 1]))
    pf = np.loadtxt(fetch("allsky_powerflux", refresh))
    out["PowerFlux (population average)"] = pd.DataFrame(dict(f_hz=pf[:, 0], h0_ul=pf[:, 1]))
    so = np.loadtxt(fetch("allsky_soap", refresh))
    out["SOAP"] = pd.DataFrame(dict(f_hz=so[:, 0], h0_ul=10.0 ** so[:, 1]))          # log10 in the file
    return out


def boson_o3(refresh: bool = False) -> dict:
    """{'ul': (f, h0 max, h0 min per 1 Hz band), 'mass': {(d_kpc, chi, age_yr): (M_BH, m_b) excluded points},
    'distance': {(M_BH max, age_yr): (m_b, d_max kpc)}}."""
    z = zipfile.ZipFile(fetch("boson_o3", refresh))
    out = {"mass": {}, "distance": {}}
    for n in z.namelist():
        base = n.rsplit("/", 1)[-1]
        if not base.endswith(".txt"):
            continue
        a = np.loadtxt(io.BytesIO(z.read(n)))
        if base == "boson_O3_UL_1Hz.txt":
            out["ul"] = pd.DataFrame(dict(f_hz=a[:, 0], h0_ul_max=a[:, 1], h0_ul_min=a[:, 2]))
        elif m := re.match(r"boson_O3_mass_exclusion_d(\d+)kpc_chi([\d.]+)_tage([\de.]+)yr\.txt", base):
            age = float(m.group(3).replace("1e4.5", "31622.8"))
            out["mass"][(float(m.group(1)), float(m.group(2)), age)] = pd.DataFrame(dict(m_bh=a[:, 0], m_b_ev=a[:, 1]))
        elif m := re.match(r"boson_O3_distance_exclusion_Mbhmax(\d+)_tage([\de.]+)yr\.txt", base):
            out["distance"][(float(m.group(1)), float(m.group(2)))] = pd.DataFrame(dict(m_b_ev=a[:, 0], d_kpc=a[:, 1]))
    return out


def distance_reach(boson: dict, m_lo: float, m_hi: float) -> pd.DataFrame:
    """Largest distance (kpc) at which O3 excludes clouds, over the boson-mass band, per black-hole mass cap and age."""
    rows = []
    for (mbh, age), d in sorted(boson["distance"].items()):
        sel = d[(d["m_b_ev"] >= m_lo) & (d["m_b_ev"] <= m_hi)]
        rows.append(dict(m_bh_max=mbh, age_yr=age, d_max_kpc_in_band=float(sel["d_kpc"].max()) if len(sel) else 0.0,
                         d_max_kpc_min_in_band=float(sel["d_kpc"].min()) if len(sel) else 0.0))
    return pd.DataFrame(rows)


def mass_exclusion_in_band(boson: dict, m_lo: float, m_hi: float) -> pd.DataFrame:
    """Black-hole masses excluded by O3 for boson masses in the band, per distance, spin and age."""
    rows = []
    for (d, chi, age), t in sorted(boson["mass"].items()):
        sel = t[(t["m_b_ev"] >= m_lo) & (t["m_b_ev"] <= m_hi)]
        rows.append(dict(distance_kpc=d, chi=chi, age_yr=age,
                         m_bh_min=float(sel["m_bh"].min()) if len(sel) else np.nan,
                         m_bh_max=float(sel["m_bh"].max()) if len(sel) else np.nan))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------
INK, INK2, GRID = "#1f2937", "#4b5563", "#e5e7eb"
BLUE, ORANGE, GREY, RED = "#2a78d6", "#e8743b", "#9ca3af", "#d64550"


def _line_handler() -> dict:
    """An explicit legend handler for lines: a module imported by the CLI registers another one for Line2D."""
    from matplotlib.legend_handler import HandlerLine2D
    from matplotlib.lines import Line2D

    return {Line2D: HandlerLine2D()}


def _style(ax):
    ax.grid(True, color=GRID, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=8)


def plot_known_pulsars(kp: pd.DataFrame, allsky: dict, out_png: Path, nb: Optional[pd.DataFrame] = None) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    for (name, c), color in zip(allsky.items(), ("#c7d2fe", "#bae6fd", "#ddd6fe")):
        ax.plot(c["f_hz"], c["h0_ul"], color=color, lw=1.0, label=f"all-sky O4a, {name}")
    if "h0_spindown" in kp:
        sd = kp.dropna(subset=["h0_spindown"])
        ax.scatter(sd["f_gw_hz"], sd["h0_spindown"], marker="v", s=22, color=GREY, label="spin-down limit", zorder=3)
        for _, r in sd.iterrows():
            ax.plot([r["f_gw_hz"]] * 2, [r["h0_ul"], r["h0_spindown"]], color=GREY, lw=0.6, zorder=2)
        below = sd["h0_ul"] < sd["h0_spindown"]
        ax.scatter(sd.loc[below, "f_gw_hz"], sd.loc[below, "h0_ul"], marker="*", s=70, color=ORANGE, zorder=4,
                   label=f"upper limit below the spin-down limit ({int(below.sum())})")
    ax.scatter(kp["f_gw_hz"], kp["h0_ul"], s=14, color=BLUE, zorder=3, label="95% upper limit, known pulsar (O4a)")
    if nb is not None and len(nb):
        ax.scatter(nb["f_gw_hz"], nb["h0_ul"], s=26, marker="D", facecolors="none", edgecolors="#7c3aed", lw=0.9,
                   zorder=3, label=f"narrowband search ({nb['pulsar'].nunique()} pulsars)")
    for name in ("J0534+2200", "J0835-4510", "J0537-6910", "J0437-4715"):
        r = kp[kp["pulsar"] == name]
        if len(r):
            ax.annotate({"J0534+2200": "Crab", "J0835-4510": "Vela"}.get(name, name),
                        (r["f_gw_hz"].iloc[0], r["h0_ul"].iloc[0]), xytext=(4, -10), textcoords="offset points",
                        fontsize=7, color=INK2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10, 2100)
    ax.set_ylim(2e-27, 2e-23)
    ax.set_xlabel("gravitational-wave frequency (Hz)", color=INK2, fontsize=9)
    ax.set_ylabel("strain amplitude h₀", color=INK2, fontsize=9)
    _style(ax)
    ax.legend(frameon=False, fontsize=7, labelcolor=INK2, loc="upper right", handler_map=_line_handler())
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png


def _band(ax, band, label, horizontal=True):
    lo, mid, hi = band
    if horizontal:
        ax.axhspan(lo, hi, color=RED, alpha=0.15, lw=0)
        ax.axhline(mid, color=RED, lw=1.2, label=label)
    else:
        ax.axvspan(lo, hi, color=RED, alpha=0.15, lw=0)
        ax.axvline(mid, color=RED, lw=1.2, label=label)


def plot_boson_exclusion(boson: dict, band, band_label: str, out_png: Path, chi: float = 0.9) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dists = sorted({k[0] for k in boson["mass"]})
    fig, axes = plt.subplots(1, len(dists), figsize=(5.0 * len(dists), 4.0), sharey=True)
    axes = np.atleast_1d(axes)
    colors = ("#1e3a8a", "#2a78d6", "#7dd3fc", "#bae6fd")
    for ax, d in zip(axes, dists):
        # from the youngest clouds to the oldest (which are excluded at lower boson masses: heavier bosons have
        # annihilated by then)
        keys = sorted([k for k in boson["mass"] if k[0] == d and k[1] == chi], key=lambda k: k[2])
        for k, col in zip(keys, colors[::-1][-len(keys):]):
            t = boson["mass"][k]
            ax.scatter(t["m_bh"], t["m_b_ev"], s=1.5, color=col, lw=0, rasterized=True,
                       label=f"excluded, cloud age {k[2]:.0e} yr".replace("e+0", "e"))
        _band(ax, band, band_label)
        ax.set_yscale("log")
        ax.set_xlabel(r"black-hole mass (M$_\odot$)", color=INK2, fontsize=9)
        ax.set_title(f"O3, black hole at {d:g} kpc, initial spin {chi:g}", color=INK, fontsize=9, loc="left")
        _style(ax)
        ax.legend(frameon=False, fontsize=7, labelcolor=INK2, markerscale=6, loc="upper right",
                  handler_map=_line_handler())
    axes[0].set_ylabel("boson mass (eV)", color=INK2, fontsize=9)
    sec = axes[-1].secondary_yaxis("right", functions=(lambda m: m * EV_TO_HZ, lambda f: f / EV_TO_HZ))
    sec.set_ylabel("CW frequency 2mc²/h (Hz)", color=INK2, fontsize=9)
    sec.set_yticks([50, 100, 200, 300, 500])
    sec.set_yticklabels(["50", "100", "200", "300", "500"])
    sec.tick_params(colors=INK2, labelsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png


def plot_boson_distance(boson: dict, band, band_label: str, out_png: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.5, 4.0))
    keys = sorted(boson["distance"], key=lambda k: (k[0], -k[1]))
    styles = {50.0: "-", 100.0: "--"}
    cmap = plt.get_cmap("viridis")
    ages = sorted({k[1] for k in keys})
    for k in keys:
        t = boson["distance"][k]
        ax.plot(t["m_b_ev"], t["d_kpc"], styles.get(k[0], ":"), color=cmap(ages.index(k[1]) / max(len(ages) - 1, 1)),
                lw=1.2, label=rf"M$_{{BH}}$ ≤ {k[0]:g} M$_\odot$, age {k[1]:.0e} yr".replace("e+0", "e"))
    _band(ax, band, band_label, horizontal=False)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("boson mass (eV)", color=INK2, fontsize=9)
    ax.set_ylabel("largest excluded distance (kpc)", color=INK2, fontsize=9)
    sec = ax.secondary_xaxis("top", functions=(lambda m: m * EV_TO_HZ, lambda f: f / EV_TO_HZ))
    sec.set_xlabel("CW frequency (Hz)", color=INK2, fontsize=9)
    sec.set_xticks([30, 50, 100, 200, 300, 500])
    sec.set_xticklabels(["30", "50", "100", "200", "300", "500"])
    sec.tick_params(colors=INK2, labelsize=8)
    _style(ax)
    ax.legend(frameon=False, fontsize=7, labelcolor=INK2, loc="center left", bbox_to_anchor=(1.01, 0.5),
              handlelength=2.5, handler_map=_line_handler())
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
# mode
# ---------------------------------------------------------------------------
def boson_paragraph(lo, hi, f_lo, f_hi, boson_label, res, excl, far, default_band: bool = True) -> str:
    """The boson-cloud paragraph of the report, for the band [lo, hi] eV radiating at [f_lo, f_hi] Hz."""
    inside = 20 <= f_lo and f_hi <= 610
    text = (f"<b>Boson clouds, O3</b> ({REFS['boson']}). A cloud of bosons of mass m radiates at 2mc²/h; the band "
            f"{lo:.2g}–{hi:.2g} eV ({boson_label}) radiates at {f_lo:.0f}–{f_hi:.0f} Hz, "
            + ("inside" if inside else "<b>not entirely inside</b>") + " the 20–610 Hz of the search"
            + (f", where its h0 upper limits are {res['boson_ul_in_band'][0]:.2g}–{res['boson_ul_in_band'][1]:.2g}"
               if res["boson_ul_in_band"] else ": the search has no upper limit there") + ".")
    light = [f"{r.m_bh_min:.0f}–{r.m_bh_max:.0f} M☉ at {r.distance_kpc:g} kpc"
             for r in excl[(excl["chi"] == 0.9) & (excl["age_yr"] == 1e3)].itertuples() if np.isfinite(r.m_bh_min)]
    if light:
        text += (" There the excluded black holes are the light ones, " + ", ".join(light)
                 + " (spin 0.9, cloud age 10³ yr): bosons this heavy are superradiant only around light black holes"
                 + (", the range of GW190728's own components." if default_band else "."))
    else:
        text += " No black-hole mass is excluded in this band (spin 0.9, cloud age 10³ yr)."
    reach = [f"{r.d_max_kpc_in_band:.2g} kpc (M_BH ≤ {r.m_bh_max:g} M☉, {r.age_yr:.0e} yr)"
             for r in far.itertuples() if r.d_max_kpc_in_band > 0]
    if reach:
        text += " Clouds younger than 10⁵ yr are excluded out to " + ", ".join(reach) + "."
    return text + (" These limits concern clouds around isolated black holes of the Galaxy, with assumptions on "
                   "their mass, spin, age and distance; they constrain the boson only together with a population of "
                   "black holes.")


def run_continuous_waves(boson_mass_ev=GW190728_BOSON_MASS_EV, boson_label: str = "GW190728 (scalar_cloud)",
                         out_report_html: Optional[str] = "continuous_waves.html",
                         out_summary_tsv: Optional[str] = "continuous_waves.tsv",
                         plots_dir: str = "continuous_waves_plots", refresh: bool = False) -> dict:
    """The CW upper limits of the LVK data releases, with a boson-mass band on the boson-cloud plots."""
    from .report import write_simple_html_report

    lo, hi = float(min(boson_mass_ev)), float(max(boson_mass_ev))
    band = (lo, float(np.median(boson_mass_ev)) if len(boson_mass_ev) != 3 else float(boson_mass_ev[1]), hi)
    plots = Path(plots_dir)
    plots.mkdir(parents=True, exist_ok=True)
    kp, sky, bos = known_pulsars(refresh), allsky_o4a(refresh), boson_o3(refresh)
    nb = narrowband(kp, refresh)
    reach = distance_reach(bos, lo, hi)
    excl = mass_exclusion_in_band(bos, lo, hi)
    f_lo, f_hi = boson_frequency_hz([lo, hi])
    ul_band = bos["ul"][(bos["ul"]["f_hz"] >= f_lo) & (bos["ul"]["f_hz"] <= f_hi)]

    res = dict(n_pulsars=len(kp))
    if "ratio_to_spindown" in kp:
        sd = kp.dropna(subset=["h0_spindown"])
        res["n_below_spindown"] = int((sd["h0_ul"] < sd["h0_spindown"]).sum())
    best_h0 = kp.loc[kp["h0_ul"].idxmin()]
    res.update(best_h0=float(best_h0["h0_ul"]), best_h0_pulsar=best_h0["pulsar"])
    if "ellipticity_ul" in kp and kp["ellipticity_ul"].notna().any():
        be = kp.loc[kp["ellipticity_ul"].idxmin()]
        res.update(best_ellipticity=float(be["ellipticity_ul"]), best_ellipticity_pulsar=be["pulsar"])
    nb_main = nb[~nb["post_glitch"]]
    res["n_narrowband"] = int(nb_main["pulsar"].nunique())
    if "ratio_to_spindown" in nb:
        res["n_narrowband_below_spindown"] = int((nb_main["ratio_to_spindown"] < 1).sum())
        crab = nb_main[nb_main["pulsar"] == "J0534+2200"]
        if len(crab):
            res["crab_narrowband_ratio"] = float(crab["ratio_to_spindown"].iloc[0])
    allmin = {k: (float(v["h0_ul"].min()), float(v.loc[v["h0_ul"].idxmin(), "f_hz"])) for k, v in sky.items()}
    res.update(boson_mass_band_ev=(lo, hi), cw_band_hz=(float(f_lo), float(f_hi)),
               boson_ul_in_band=(float(ul_band["h0_ul_max"].min()), float(ul_band["h0_ul_max"].max()))
               if len(ul_band) else None)

    rows = [dict(item="known pulsars (O4a)", value=f"{len(kp)}", ref=REFS["known_pulsars"])]
    if "n_below_spindown" in res:
        rows.append(dict(item="upper limit below the spin-down limit (Bayesian pipeline; any of the three: 29)",
                         value=f"{res['n_below_spindown']}", ref=REFS["known_pulsars"]))
    rows.append(dict(item="lowest h0 upper limit", value=f"{res['best_h0']:.2g} ({res['best_h0_pulsar']})",
                     ref=REFS["known_pulsars"]))
    if "best_ellipticity" in res:
        rows.append(dict(item="lowest ellipticity upper limit",
                         value=f"{res['best_ellipticity']:.2g} ({res['best_ellipticity_pulsar']})",
                         ref=REFS["known_pulsars"]))
    if "n_narrowband_below_spindown" in res:
        rows.append(dict(item="narrowband search: pulsars below the spin-down limit",
                         value=f"{res['n_narrowband_below_spindown']} of {res['n_narrowband']}"
                               + (f" (Crab at {100 * res['crab_narrowband_ratio']:.0f}% of it)"
                                  if "crab_narrowband_ratio" in res else ""), ref=REFS["known_pulsars"]))
    for k, (h, f) in allmin.items():
        rows.append(dict(item=f"all-sky O4a {k}, lowest upper limit", value=f"{h:.2g} at {f:.0f} Hz", ref=REFS["allsky"]))
    rows.append(dict(item=f"boson mass band ({boson_label})", value=f"{lo:.3g}–{hi:.3g} eV, CW at {f_lo:.0f}–{f_hi:.0f} Hz",
                     ref=""))
    if res["boson_ul_in_band"]:
        rows.append(dict(item="O3 boson-cloud h0 upper limits in that band (max CR per 1 Hz)",
                         value=f"{res['boson_ul_in_band'][0]:.2g}–{res['boson_ul_in_band'][1]:.2g}", ref=REFS["boson"]))
    table = pd.DataFrame(rows)
    if out_summary_tsv:
        table.to_csv(out_summary_tsv, sep="\t", index=False)
        kp.to_csv(Path(out_summary_tsv).with_suffix(".pulsars.tsv"), sep="\t", index=False, float_format="%.4g")
        nb.to_csv(Path(out_summary_tsv).with_suffix(".narrowband.tsv"), sep="\t", index=False, float_format="%.4g")
        reach.to_csv(Path(out_summary_tsv).with_suffix(".boson_reach.tsv"), sep="\t", index=False, float_format="%.4g")
    if out_report_html:
        import matplotlib

        matplotlib.use("Agg")
        with matplotlib.rc_context(matplotlib.rcParamsDefault):   # other modes (ligo.skymap) change the style
            images = [plot_known_pulsars(kp, sky, plots / "cw_known_pulsars.png", nb),
                      plot_boson_exclusion(bos, band, boson_label, plots / "cw_boson_exclusion.png"),
                      plot_boson_distance(bos, band, boson_label, plots / "cw_boson_distance.png")]
        far = reach[reach["age_yr"] <= 1e5]
        paras = [
            "No continuous wave has been detected; these are the upper limits of the LVK data releases.",
            f"<b>Known pulsars, O4a</b> ({REFS['known_pulsars']}): {len(kp)} pulsars"
            + (f", {res['n_below_spindown']} with an upper limit below the spin-down limit with the Bayesian pipeline "
               "(29 with any of the paper's three)" if "n_below_spindown" in res else "") + f"; lowest h0 {res['best_h0']:.2g} ({res['best_h0_pulsar']})"
            + (f", lowest ellipticity {res['best_ellipticity']:.2g} ({res['best_ellipticity_pulsar']})."
               if "best_ellipticity" in res else "."),
            f"<b>All-sky, O4a</b> ({REFS['allsky']}): " + "; ".join(f"{k} {h:.2g} at {f:.0f} Hz"
                                                                     for k, (h, f) in allmin.items()) + ".",
            boson_paragraph(lo, hi, f_lo, f_hi, boson_label, res, excl, far,
                            default_band=tuple(boson_mass_ev) == GW190728_BOSON_MASS_EV),
        ]
        write_simple_html_report(out_report_html, title="Continuous waves: LVK upper limits", paragraphs=paras,
                                 images=[str(i) for i in images],
                                 tables=[("Summary", table.to_html(index=False)),
                                         ("O3 boson clouds: excluded black-hole masses in the band",
                                          excl.to_html(index=False, float_format=lambda x: f"{x:.3g}", na_rep="none"))])
        _log(f"report written to {out_report_html}")
    return res

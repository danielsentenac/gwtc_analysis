"""Fitted Power Law + Peak and Multi Peak mass models over the selection-corrected primary-mass distribution.

Reads docs/scripts/data/mass_distribution_selection_corrected.tsv (written by
mass_distribution_selection_corrected.py) and docs/scripts/data/mass_models_fit.tsv (written by mass_models_fit.py
from the hubble_constant runs), and writes docs/img/mass_models_plp_mltp.png. The models are fitted without their
rate (scale-free likelihood): each curve is scaled to the rate of the histogram between 7 and 100 M_sun.
Run from anywhere:
    python docs/scripts/mass_models_plp_mltp.py
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "img" / "mass_models_plp_mltp.png"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BARS = "#c9c7c1"
MODELS = (("plp", "Power Law + Peak", "#eb6834"), ("mltp", "Multi Peak", "#2a78d6"))
NORM = (7.0, 100.0)

h = np.genfromtxt(HERE / "data" / "mass_distribution_selection_corrected.tsv", delimiter="\t", names=True)
f = np.genfromtxt(HERE / "data" / "mass_models_fit.tsv", delimiter="\t", names=True)
centers = np.sqrt(h["m1_lo"] * h["m1_hi"])
dln = np.log(h["m1_hi"] / h["m1_lo"])
inside = (centers >= NORM[0]) & (centers <= NORM[1]) & np.isfinite(h["rate"])
hist_rate = np.sum(h["rate"][inside] * dln[inside])
m1, lnm = f["m1"], np.log(f["m1"])
win = (m1 >= NORM[0]) & (m1 <= NORM[1])

fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), dpi=150, sharey=True)
fig.patch.set_facecolor(SURFACE)
ok = np.isfinite(h["rate"]) & (h["rate"] > 0)
for ax, (key, name, color) in zip(axes, MODELS):
    ax.set_facecolor(SURFACE)
    ax.bar(centers[ok], h["rate"][ok], width=(h["m1_hi"] - h["m1_lo"])[ok] * 0.9, color=BARS, linewidth=0,
           zorder=2, label="Events ÷ ⟨VT⟩")
    ax.errorbar(centers[ok], h["rate"][ok], yerr=[(h["rate"] - h["rate_lo90"])[ok], (h["rate_hi90"] - h["rate"])[ok]],
                fmt="none", ecolor="#9a9891", elinewidth=1, zorder=2)
    scale = hist_rate / np.trapezoid(f[f"{key}_median"][win], lnm[win])
    lo, hi = scale * f[f"{key}_lo90"], scale * f[f"{key}_hi90"]
    ax.fill_between(m1, np.maximum(lo, 1e-6), np.maximum(hi, 1e-6), color=color, alpha=0.15, linewidth=0, zorder=3,
                    label="90% band")
    draws = [c for c in f.dtype.names if c.startswith(f"{key}_draw")]
    for i, c in enumerate(draws):
        d = f[c] * hist_rate / np.trapezoid(f[c][win], lnm[win])
        ax.plot(m1[d > 0], d[d > 0], color=color, lw=0.8, alpha=0.55, zorder=4,
                label=f"{len(draws)} posterior draws" if i == 0 else None)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(2, 250)
    ax.set_ylim(1e-2, 1e3)
    ax.set_xticks([2, 5, 10, 20, 35, 50, 100, 200])
    ax.set_xticklabels(["2", "5", "10", "20", "35", "50", "100", "200"])
    for x in (10, 35):
        ax.axvline(x, color=INK2, lw=0.8, ls=(0, (3, 3)), zorder=1)
    ax.grid(axis="y", which="major", color=GRID, lw=0.8, zorder=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.set_xlabel("Primary mass $m_1$ (source frame, M$_\\odot$)", color=INK2, fontsize=10)
    ax.set_title(name, color=INK, fontsize=10.5, loc="left")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2, loc="upper right")
axes[0].set_ylabel("Merger rate dR/d ln m$_1$ (Gpc$^{-3}$ yr$^{-1}$)", color=INK2, fontsize=10)
fig.text(0.01, 0.01, "Models: hubble_constant fits of the 137 GWTC-4.0 BBHs (all injections, 5000 PE samples), shape "
         "only, each curve scaled to the histogram's rate between 7 and 100 M$_\\odot$.\nHistogram: GWTC-1 to "
         "GWTC-5.0 events divided by ⟨VT⟩, as on the home page (m$_2$ uniform in [1, m$_1$], no redshift evolution), "
         "not a population fit. Dashed: 10 and 35 M$_\\odot$.", fontsize=7, color=INK2)
fig.tight_layout(rect=(0, 0.07, 1, 1))
fig.savefig(OUT, facecolor=SURFACE)
print("wrote", OUT)

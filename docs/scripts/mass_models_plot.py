"""Fitted Power Law + Peak, Multi Peak and FullPop-4.0 mass models over the selection-corrected primary-mass distribution.

Reads docs/scripts/data/mass_distribution_selection_corrected.tsv (written by
mass_distribution_selection_corrected.py) and docs/scripts/data/mass_models_fit.tsv (written by mass_models_fit.py
from the hubble_constant runs), and writes docs/img/mass_models.png. The models are fitted without their
rate (scale-free likelihood): each curve is scaled to the rate of the histogram between 7 and 100 M_sun. Each
model is drawn at the posterior median of its parameters, with the 90% band of the posterior at each mass.
Run from anywhere:
    python docs/scripts/mass_models_plot.py
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "img" / "mass_models.png"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BARS = "#c9c7c1"
MODELS = (("plp", "Power Law + Peak", "#eb6834", "-"), ("mltp", "Multi Peak", "#2a78d6", "-"),
          ("fullpop4", "FullPop-4.0", "#1baf7a", (0, (5, 2))))
NORM = (7.0, 100.0)

h = np.genfromtxt(HERE / "data" / "mass_distribution_selection_corrected.tsv", delimiter="\t", names=True)
f = np.genfromtxt(HERE / "data" / "mass_models_fit.tsv", delimiter="\t", names=True)
centers = np.sqrt(h["m1_lo"] * h["m1_hi"])
dln = np.log(h["m1_hi"] / h["m1_lo"])
inside = (centers >= NORM[0]) & (centers <= NORM[1]) & np.isfinite(h["rate"])
hist_rate = np.sum(h["rate"][inside] * dln[inside])
m1, lnm = f["m1"], np.log(f["m1"])
win = (m1 >= NORM[0]) & (m1 <= NORM[1])

fig, ax = plt.subplots(figsize=(8.2, 4.6), dpi=150)
fig.patch.set_facecolor(SURFACE)
ax.set_facecolor(SURFACE)
ok = np.isfinite(h["rate"]) & (h["rate"] > 0)
ax.bar(centers[ok], h["rate"][ok], width=(h["m1_hi"] - h["m1_lo"])[ok] * 0.9, color=BARS, linewidth=0, zorder=2,
       label="Events ÷ ⟨VT⟩ (GWTC-1 to GWTC-5.0)")
ax.errorbar(centers[ok], h["rate"][ok], yerr=[(h["rate"] - h["rate_lo90"])[ok], (h["rate_hi90"] - h["rate"])[ok]],
            fmt="none", ecolor="#9a9891", elinewidth=1, zorder=2)
for key, name, color, style in MODELS:
    curve = f[f"{key}_medpar"]
    scale = hist_rate / np.trapezoid(curve[win], lnm[win])
    lo, hi = (hist_rate / np.trapezoid(f[f"{key}_median"][win], lnm[win]) * f[f"{key}_{s}"] for s in ("lo90", "hi90"))
    ax.fill_between(m1, np.maximum(lo, 1e-6), np.maximum(hi, 1e-6), color=color, alpha=0.12, linewidth=0, zorder=3)
    pos = curve > 0
    ax.plot(m1[pos], scale * curve[pos], color=color, lw=2, ls=style, zorder=4, label=f"{name} fit, with its 90% band")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1, 250)
ax.set_ylim(1e-2, 1e3)
ax.set_xticks([1, 2, 5, 10, 20, 35, 50, 100, 200])
ax.set_xticklabels(["1", "2", "5", "10", "20", "35", "50", "100", "200"])
for x, t in ((1.4, "neutron stars"), (10, "~10 M$_\\odot$ peak"), (35, "~35 M$_\\odot$")):
    ax.axvline(x, color=INK2, lw=0.8, ls=(0, (3, 3)), zorder=1)
    ax.text(x * 1.04, 0.97, t, transform=ax.get_xaxis_transform(), color=INK2, fontsize=8.5, va="top")
ax.grid(axis="y", which="major", color=GRID, lw=0.8, zorder=0)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ax.spines[sp].set_color(GRID)
ax.tick_params(colors=INK2, labelsize=9)
ax.set_xlabel("Primary mass $m_1$ (source frame, M$_\\odot$)", color=INK2, fontsize=10)
ax.set_ylabel("Merger rate dR/d ln m$_1$ (Gpc$^{-3}$ yr$^{-1}$)", color=INK2, fontsize=10)
ax.set_title("Fitted mass models (GWTC-4.0 spectral sirens) over the selection-corrected masses",
             color=INK, fontsize=10.5, loc="left")
ax.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper right", bbox_to_anchor=(1, 0.9))
fig.text(0.01, 0.01, "Curves: the models at the posterior median of each parameter (Power Law + Peak, Multi Peak: 137 "
         "GWTC-4.0 BBHs; FullPop-4.0: 141 events,\nwith neutron stars; all injections, 5000 PE samples); bands: 90% of "
         "the posterior at each mass. Shape only, scaled to the histogram's rate\nbetween 7 and 100 M$_\\odot$. "
         "Histogram: as on the home page, not a population fit.", fontsize=7, color=INK2)
fig.tight_layout(rect=(0, 0.08, 1, 1))
fig.savefig(OUT, facecolor=SURFACE)
print("wrote", OUT)

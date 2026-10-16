"""Primary-mass (source frame) distribution of the GWOSC confident catalogs: docs/img/mass_distribution.png.

The figure of the home page of the docs, as detected (not corrected for selection effects); its corrected companion
is mass_distribution_selection_corrected.py. Run from anywhere: python docs/scripts/mass_distribution.py
"""
import json, urllib.request
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
CATALOGS = ("GWTC-1-confident", "GWTC-2.1-confident", "GWTC-3-confident", "GWTC-4.0", "GWTC-5.0")

events = {}
for cat in CATALOGS:
    for v in json.load(urllib.request.urlopen(f"https://gwosc.org/eventapi/jsonfull/{cat}/"))["events"].values():
        m, lo, hi = v.get("mass_1_source"), v.get("mass_1_source_lower"), v.get("mass_1_source_upper")
        if m is None:
            continue
        events[v["commonName"]] = (float(m), float(lo or 0.0), float(hi or 0.0), cat)   # later catalogs win

m1 = np.array([e[0] for e in events.values()])
print(f"{len(m1)} events with m1_source; range {m1.min():.2f}-{m1.max():.1f} Msun")

# Stacked uncertainty: each event as a split normal in ln m, from its 90% interval
grid = np.logspace(0, np.log10(250), 1200)
lg = np.log(grid)
dens = np.zeros_like(grid)
for m, lo, hi, _ in events.values():
    s_lo = max(np.log(m) - np.log(max(m + lo, 0.5)), 0.02) / 1.645   # lower is negative in GWOSC
    s_hi = max(np.log(m + hi) - np.log(m), 0.02) / 1.645
    s = np.where(lg < np.log(m), s_lo, s_hi)
    dens += np.exp(-0.5 * ((lg - np.log(m)) / s) ** 2) * 2 / (np.sqrt(2 * np.pi) * (s_lo + s_hi))

bins = np.logspace(0, np.log10(250), 41)
dlog = np.log(bins[1]) - np.log(bins[0])
counts, _ = np.histogram(m1, bins)

fig, ax = plt.subplots(figsize=(8.2, 4.4), dpi=150)
fig.patch.set_facecolor(SURFACE); ax.set_facecolor(SURFACE)
centers = np.sqrt(bins[1:] * bins[:-1])
widths = np.diff(bins)
ax.bar(centers, counts, width=widths * 0.9, color=BLUE, alpha=0.85, linewidth=0, label="Events (median mass)", zorder=2)
ax.plot(grid, dens * dlog, color=ORANGE, lw=2, label="Sum of event uncertainties", zorder=3)
ax.set_xscale("log")
ax.set_xlim(1, 250)
ax.set_xticks([1, 2, 5, 10, 20, 35, 50, 100, 200])
ax.set_xticklabels(["1", "2", "5", "10", "20", "35", "50", "100", "200"])
ax.set_xlabel("Primary mass $m_1$ (source frame, M$_\\odot$)", color=INK2)
ax.set_ylabel("Events per bin", color=INK2)
ax.grid(axis="y", color=GRID, lw=0.8, zorder=0)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(GRID)
ax.tick_params(colors=INK2, labelsize=9)
ymax = counts.max()
for x, text in ((1.45, "neutron stars"), (10, "~10 M$_\\odot$ peak"), (35, "~35 M$_\\odot$ peak")):
    ax.annotate(text, xy=(x, ymax * 1.02), ha="center", va="bottom", fontsize=9, color=INK2)
    ax.axvline(x, color=INK2, lw=0.8, ls=(0, (3, 3)), zorder=1)
ax.set_ylim(0, ymax * 1.18)
ax.set_title(f"Primary masses of {len(m1)} confident GW events (GWTC-1 to GWTC-5.0, GWOSC)", color=INK, fontsize=11, loc="left")
ax.legend(frameon=False, fontsize=9, loc="upper right", bbox_to_anchor=(1.0, 0.86), labelcolor=INK2)
fig.text(0.01, 0.01, "Observed distribution, not corrected for selection effects: heavy binaries are over-represented.",
         fontsize=8, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(Path(__file__).resolve().parents[1] / "img" / "mass_distribution.png", facecolor=SURFACE)

for lo_, hi_ in ((1, 3), (3, 6), (6, 15), (15, 25), (25, 45), (45, 80), (80, 250)):
    print(f"  m1 {lo_:>3}-{hi_:<3} Msun: {int(((m1 >= lo_) & (m1 < hi_)).sum()):3d} events")

"""Primary-mass distribution in the detector frame: docs/img/mass_distribution_detector_frame.png.

The first figure of the home page of the docs: the masses as the waveform measures them, m1_det = (1 + z) m1_source,
before any cosmology is assumed. Its companions are mass_distribution.py (source frame) and
mass_distribution_selection_corrected.py. Run from anywhere: python docs/scripts/mass_distribution_detector_frame.py
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

# GWOSC lists source-frame masses and the redshift (Planck15); (1 + z) undoes the conversion. The product of the
# medians approximates the median detector-frame mass; the interval ignores the redshift uncertainty.
events = {}
for cat in CATALOGS:
    for v in json.load(urllib.request.urlopen(f"https://gwosc.org/eventapi/jsonfull/{cat}/"))["events"].values():
        m, lo, hi, z = (v.get(k) for k in ("mass_1_source", "mass_1_source_lower", "mass_1_source_upper", "redshift"))
        if m is None or z is None:
            continue
        events[v["commonName"]] = (float(m), float(lo or 0.0), float(hi or 0.0), float(z))   # later catalogs win

m1_src = np.array([e[0] for e in events.values()])
z = np.array([e[3] for e in events.values()])
m1 = m1_src * (1 + z)
print(f"{len(m1)} events with m1_source and z; detector frame {m1.min():.2f}-{m1.max():.1f} Msun; "
      f"median z {np.median(z):.2f}")

# Stacked uncertainty: each event as a split normal in ln m, from its 90% interval scaled by (1 + z)
grid = np.logspace(0, np.log10(400), 1200)
lg = np.log(grid)
dens = np.zeros_like(grid)
for m, lo, hi, zz in events.values():
    m, lo, hi = m * (1 + zz), lo * (1 + zz), hi * (1 + zz)
    s_lo = max(np.log(m) - np.log(max(m + lo, 0.5)), 0.02) / 1.645   # lower is negative in GWOSC
    s_hi = max(np.log(m + hi) - np.log(m), 0.02) / 1.645
    s = np.where(lg < np.log(m), s_lo, s_hi)
    dens += np.exp(-0.5 * ((lg - np.log(m)) / s) ** 2) * 2 / (np.sqrt(2 * np.pi) * (s_lo + s_hi))

bins = np.logspace(0, np.log10(400), 44)
dlog = np.log(bins[1]) - np.log(bins[0])
counts, _ = np.histogram(m1, bins)
counts_src, _ = np.histogram(m1_src, bins)

fig, ax = plt.subplots(figsize=(8.2, 4.4), dpi=150)
fig.patch.set_facecolor(SURFACE); ax.set_facecolor(SURFACE)
centers = np.sqrt(bins[1:] * bins[:-1])
widths = np.diff(bins)
ax.bar(centers, counts, width=widths * 0.9, color=BLUE, alpha=0.85, linewidth=0, label="Events (median mass)", zorder=2)
ax.plot(grid, dens * dlog, color=ORANGE, lw=2, label="Sum of event uncertainties", zorder=3)
ax.stairs(counts_src, bins, color=INK2, lw=1, ls=(0, (2, 2)), label="Same events, source frame", zorder=3)
ax.set_xscale("log")
ax.set_xlim(1, 400)
ax.set_xticks([1, 2, 5, 10, 20, 35, 50, 100, 200, 400])
ax.set_xticklabels(["1", "2", "5", "10", "20", "35", "50", "100", "200", "400"])
ax.set_xlabel("Primary mass $m_1(1+z)$ (detector frame, M$_\\odot$)", color=INK2)
ax.set_ylabel("Events per bin", color=INK2)
ax.grid(axis="y", color=GRID, lw=0.8, zorder=0)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(GRID)
ax.tick_params(colors=INK2, labelsize=9)
ymax = max(counts.max(), counts_src.max())
for x, text in ((1.45, "neutron stars"), (10, "10 M$_\\odot$ (source)"), (35, "35 M$_\\odot$ (source)")):
    ax.annotate(text, xy=(x, ymax * 1.02), ha="center", va="bottom", fontsize=9, color=INK2)
    ax.axvline(x, color=INK2, lw=0.8, ls=(0, (3, 3)), zorder=1)
ax.set_ylim(0, ymax * 1.18)
ax.set_title(f"Primary masses of {len(m1)} confident GW events, detector frame (GWTC-1 to GWTC-5.0, GWOSC)",
             color=INK, fontsize=11, loc="left")
ax.legend(frameon=False, fontsize=9, loc="upper right", bbox_to_anchor=(1.0, 0.86), labelcolor=INK2)
fig.text(0.01, 0.01, "As the waveform measures them, before any cosmology. The redshift shifts and widens "
         "the 35 M$_\\odot$ peak.", fontsize=8, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(Path(__file__).resolve().parents[1] / "img" / "mass_distribution_detector_frame.png", facecolor=SURFACE)

for lo_, hi_ in ((1, 3), (3, 6), (6, 15), (15, 25), (25, 45), (45, 80), (80, 400)):
    print(f"  m1(1+z) {lo_:>3}-{hi_:<3} Msun: {int(((m1 >= lo_) & (m1 < hi_)).sum()):3d} events")

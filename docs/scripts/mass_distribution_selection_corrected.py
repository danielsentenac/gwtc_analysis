"""Primary-mass distribution of the GW events, corrected for selection effects, with gwtc_analysis.

Same figure as docs/img/mass_distribution.png (median-mass histogram + sum of the event uncertainties), with each
bin and each point of the curve divided by the sensitive volume-time <VT> of binaries of that primary mass:
dR/dln m1 = N_bin / (<VT>_bin dln m1). The events and the injections are those of the `rates` mode (FAR at most
1/yr, inside the periods covered by the injections, GWTC-1 to GWTC-5 with the gwtc5 sensitivity release, the
semi-analytic O1-O2 injections included), and <VT> is computed with its `_sensitive_vt`.

Writes docs/img/mass_distribution_selection_corrected.png. Run from anywhere:
    python docs/scripts/mass_distribution_selection_corrected.py
The first run downloads the injection file of the release (about 1 GB) into the Zenodo cache; it needs about
3 GB of memory.
"""
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # the repository
from gwtc_analysis import catalogs as cat                      # noqa: E402
from gwtc_analysis import catalog_registry as reg              # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "img"
CATALOGS = ["GWTC-1", "GWTC-2.1", "GWTC-3", "GWTC-4", "GWTC-5"]
FAR, SNR, NS_MAX = 1.0, 10.0, 2.5
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"

# --- injections and events, as the rates mode selects them --------------------------------------------------------
runs = cat.catalog_runs(CATALOGS)
inj = cat._load_found_injections(cat._rates_sensitivity_path(None, "gwtc5", semi_analytic=True), FAR, SNR, runs=runs)
events = cat._rates_events(inj["segments"], FAR, NS_MAX, reg.gwosc_lists(CATALOGS))
events = events[events["class"] != "unknown"]
print(f"{len(inj['redshift'])} found injections; {len(events)} events (FAR <= {FAR}/yr, inside the injection periods)")

# 90% bounds of m1 (GWOSC gives offsets), for the uncertainty curve
bounds = {}
for lst in reg.gwosc_lists(CATALOGS):
    for v in json.load(urllib.request.urlopen(f"https://gwosc.org/eventapi/jsonfull/{lst}/"))["events"].values():
        if v.get("mass_1_source") is not None:
            bounds[v["commonName"]] = (float(v["mass_1_source"]), float(v.get("mass_1_source_lower") or 0.0),
                                       float(v.get("mass_1_source_upper") or 0.0))

# --- <VT> per bin of m1, as in rates' mass-distribution plot: m1 uniform in ln m1 within the bin, m2 | m1 ~ U[1, m1],
#     isotropic spins up to 0.99, no redshift evolution -------------------------------------------------------------
m1i, m2i = inj["mass1_source"], inj["mass2_source"]
with np.errstate(divide="ignore", invalid="ignore"):
    ln_pair = np.where((m2i >= 1) & (m2i <= m1i), -np.log(np.maximum(m1i - 1, 1e-9)), -np.inf)
edges = np.logspace(0, np.log10(250), 41)                     # the bins of the observed figure
dln = np.log(edges[1]) - np.log(edges[0])
vt = np.zeros(len(edges) - 1)
neff = np.zeros(len(edges) - 1)
for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
    with np.errstate(divide="ignore", invalid="ignore"):
        ln_m1 = np.where((m1i >= lo) & (m1i < hi), -np.log(m1i) - np.log(np.log(hi / lo)), -np.inf)
    vt[i], neff[i] = cat._sensitive_vt(inj, ln_m1 + ln_pair, 0.99, 0.99, kappa=0.0)
ok = (vt > 0) & (neff >= 20)                                   # bins with a reliable <VT>
centers = np.sqrt(edges[1:] * edges[:-1])

# --- corrected histogram and corrected sum of the event uncertainties ---------------------------------------------
em1 = np.maximum(events["mass_1_source"].astype(float), events["mass_2_source"].astype(float)).to_numpy()
counts, _ = np.histogram(em1, edges)
rate_bins = np.where(ok, counts / np.where(ok, vt, 1) / dln, np.nan)
lo90, hi90 = (np.array([cat._rate_quantiles(int(n), v)[k] if o and n else np.nan for n, v, o in zip(counts, vt, ok)]) / dln
              for k in (0, 2))

grid = np.logspace(0, np.log10(250), 1200)
lg = np.log(grid)
dens = np.zeros_like(grid)
missing = [n for n in events["event"] if n not in bounds]
for name in events["event"]:
    if name not in bounds:
        continue
    m, lo, hi = bounds[name]
    s_lo = max(np.log(m) - np.log(max(m + lo, 0.5)), 0.02) / 1.645
    s_hi = max(np.log(m + hi) - np.log(m), 0.02) / 1.645
    s = np.where(lg < np.log(m), s_lo, s_hi)
    dens += np.exp(-0.5 * ((lg - np.log(m)) / s) ** 2) * 2 / (np.sqrt(2 * np.pi) * (s_lo + s_hi))
# divide by <VT>(m1), interpolated in ln m1 between the reliable bins
vt_grid = np.exp(np.interp(lg, np.log(centers[ok]), np.log(vt[ok]), left=np.nan, right=np.nan))
curve = dens / vt_grid

# --- figure --------------------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.2, 4.6), dpi=150)
fig.patch.set_facecolor(SURFACE)
ax.set_facecolor(SURFACE)
ax.bar(centers, rate_bins, width=np.diff(edges) * 0.9, color=BLUE, alpha=0.85, linewidth=0, zorder=2,
       label="Events (median mass) ÷ ⟨VT⟩ of the bin")
ax.errorbar(centers, rate_bins, yerr=[rate_bins - lo90, hi90 - rate_bins], fmt="none", ecolor=BLUE, elinewidth=1,
            alpha=0.6, zorder=2)
ax.plot(grid, curve, color=ORANGE, lw=2, zorder=3, label="Sum of event uncertainties ÷ ⟨VT⟩")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1, 250)
ax.set_xticks([1, 2, 5, 10, 20, 35, 50, 100, 200])
ax.set_xticklabels(["1", "2", "5", "10", "20", "35", "50", "100", "200"])
for x, t in ((1.4, "neutron stars"), (10, "~10 M$_\\odot$ peak"), (35, "~35 M$_\\odot$ peak")):
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
ax.set_title(f"Primary masses of {len(em1)} GW events, corrected for selection effects (GWTC-1 to GWTC-5.0)",
             color=INK, fontsize=10.5, loc="left")
ax.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper right", bbox_to_anchor=(1, 0.9))
fig.text(0.01, 0.01, "⟨VT⟩ from the LVK injections (gwtc5 release, FAR ≤ 1/yr): m$_2$ uniform in [1, m$_1$], isotropic "
         "spins,\nno redshift evolution. 90% Poisson intervals; bins with fewer than 20 effective injections left out.",
         fontsize=7, color=INK2)
fig.tight_layout(rect=(0, 0.06, 1, 1))
out = OUT / "mass_distribution_selection_corrected.png"
fig.savefig(out, facecolor=SURFACE)
print("wrote", out, "| events without m1 bounds:", missing)
# the corrected histogram, for the overlay of the fitted mass models (mass_models_plot.py)
data = Path(__file__).resolve().parent / "data" / "mass_distribution_selection_corrected.tsv"
data.parent.mkdir(exist_ok=True)
np.savetxt(data, np.column_stack([edges[:-1], edges[1:], counts, rate_bins, lo90, hi90]), delimiter="\t",
           header="m1_lo\tm1_hi\tevents\trate\trate_lo90\trate_hi90", comments="", fmt="%.6g")
print("wrote", data)
for c, n, v, r, ne in zip(centers, counts, vt, rate_bins, neff):
    if n:
        print(f"m1 ~ {c:6.1f}: N = {n:3d}  <VT> = {v:9.3g} Gpc3 yr  n_eff = {ne:7.0f}  dR/dln m1 = {r:8.3g}")

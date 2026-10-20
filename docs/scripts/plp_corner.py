"""Corner plot of H0 and the population parameters most correlated with it, for a hubble_constant spectral siren.

Reads the runs' posterior (posterior.tsv) of a finished `hubble_constant --method spectral` work directory and
weights each sample by exp(ln L_target - ln L_runs), from the log-likelihoods stored by the reweighting stage
(reweight/chunk*.npz): the posterior with all the injections and the target PE samples, without the duplicates of
the resampled posterior_reweighted.tsv. Writes docs/img/plp_corner.png (or --out).

Runs with any Python that has numpy, pandas, matplotlib and corner (the icarogw environment has them):
    python docs/scripts/plp_corner.py WORKDIR [--out FILE]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import corner

PARAMS = {"H0": r"$H_0$", "mu_g": r"$\mu_g$", "alpha": r"$\alpha$", "lambda_peak": r"$\lambda$",
          "delta_m": r"$\delta_m$", "gamma": r"$\gamma$"}
UNITS = {"H0": "km/s/Mpc", "mu_g": "M$_\\odot$", "delta_m": "M$_\\odot$"}
SURFACE, INK, INK2, BLUE = "#fcfcfb", "#0b0b0b", "#52514e", "#2a78d6"

ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument("workdir")
ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "img" / "plp_corner.png"))
a = ap.parse_args()
wd = Path(a.workdir)

post = pd.read_csv(wd / "posterior.tsv", sep="\t")
chunks = sorted((wd / "reweight").glob("chunk*.npz"))
weights = np.ones(len(post))
if chunks:
    d = [np.load(c) for c in chunks]
    idx = np.concatenate([x["idx"] for x in d])
    dl = np.concatenate([x["target"] - x["runs"] for x in d])
    weights = np.zeros(len(post))
    weights[idx] = np.exp(dl - np.nanmax(dl))
    weights = np.nan_to_num(weights)
ess = weights.sum() ** 2 / (weights ** 2).sum()
print(f"{len(post)} samples, effective sample size {ess:.0f}")

x = post[list(PARAMS)].to_numpy()
fig = corner.corner(x, weights=weights, labels=[PARAMS[k] + (f" ({UNITS[k]})" if k in UNITS else "")
                                                for k in PARAMS],
                    titles=list(PARAMS.values()), range=[0.995] * len(PARAMS),
                    color=BLUE, smooth=1.2, smooth1d=1.0, bins=30, levels=(0.68, 0.9), plot_datapoints=False,
                    plot_density=False, fill_contours=True, quantiles=(0.05, 0.5, 0.95), show_titles=True,
                    title_quantiles=(0.05, 0.5, 0.95), title_fmt=".3g", title_kwargs={"fontsize": 10},
                    label_kwargs={"fontsize": 12, "color": INK2}, hist_kwargs={"linewidth": 1.5})
fig.set_size_inches(9, 9)
fig.patch.set_facecolor(SURFACE)
for ax in fig.get_axes():
    ax.set_facecolor(SURFACE)
    ax.tick_params(labelsize=9, colors=INK2)
fig.suptitle("Power Law + Peak spectral siren: H$_0$ and the parameters most correlated with it\n"
             "(68% and 90% regions; titles: median and 90% interval)", color=INK, fontsize=12, x=0.98, y=0.98,
             ha="right")
fig.savefig(a.out, dpi=150, facecolor=SURFACE)
print("wrote", a.out)

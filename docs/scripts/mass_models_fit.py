"""Primary-mass distribution of fitted hubble_constant mass models (Power Law + Peak, Multi Peak), with icarogw.

For each finished `hubble_constant --method spectral` work directory given, evaluates the source-frame
primary-mass density p(m1 | Λ) of its mass model for each posterior sample Λ (posterior_reweighted.tsv, else
posterior.tsv), with icarogw's own model (the mass prior smoothed at m_min, as in the likelihood), and writes the
median and the 90% band of m1·p(m1) (the density per ln m1), and the model at the posterior median of each
parameter, on a grid of m1:
docs/scripts/data/mass_models_fit.tsv, read by mass_models_plp_mltp.py.

Runs with the Python of the icarogw environment, from anywhere:
    /path/to/icarogw-env/bin/python docs/scripts/mass_models_fit.py WORKDIR_PLP WORKDIR_MLTP [--out FILE]
"""
import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np

MODELS = {"plp": ("massprior_PowerLawPeak", "Power Law + Peak"), "mltp": ("massprior_MultiPeak", "Multi Peak")}
GRID = np.logspace(np.log10(2.0), np.log10(250.0), 600)

ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument("workdirs", nargs="+", help="finished hubble_constant spectral-siren work directories")
ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "data" / "mass_models_fit.tsv"))
ap.add_argument("--max-samples", type=int, default=2000, help="posterior samples used per model")
a = ap.parse_args()

# icarogw reads `config` (CUPY=False: CPU mode) from the import path; the work directories hold one
sys.path.insert(0, str(Path(a.workdirs[0]).resolve()))
warnings.filterwarnings("ignore")
import pandas as pd                                    # noqa: E402
import icarogw.wrappers as w                           # noqa: E402

cols = {"m1": GRID}
for wd in map(Path, a.workdirs):
    model = json.loads((wd / "run_settings.json").read_text())["mass_model"]
    if model not in MODELS:
        raise SystemExit(f"{wd}: mass model {model!r} not handled (plp or mltp)")
    post = pd.read_csv(wd / ("posterior_reweighted.tsv" if (wd / "posterior_reweighted.tsv").exists()
                             else "posterior.tsv"), sep="\t")
    if len(post) > a.max_samples:
        post = post.sample(a.max_samples, random_state=1)
    wrapper = getattr(w, MODELS[model][0])()
    dens = np.empty((len(post), GRID.size))
    for i, (_, row) in enumerate(post.iterrows()):
        wrapper.update(**{k: float(row[k]) for k in wrapper.population_parameters})
        p1 = w.LowpassSmoothedProb(wrapper.prior, float(row["delta_m"]))
        dens[i] = GRID * np.nan_to_num(np.asarray(p1.pdf(GRID), dtype=float))
    q = np.percentile(dens, [5, 50, 95], axis=0)
    for name, v in zip(("lo90", "median", "hi90"), q):
        cols[f"{model}_{name}"] = v
    med = post.median()     # the model at the posterior median of each parameter: a typical fit
    wrapper.update(**{k: float(med[k]) for k in wrapper.population_parameters})
    p1 = w.LowpassSmoothedProb(wrapper.prior, float(med["delta_m"]))
    cols[f"{model}_medpar"] = GRID * np.nan_to_num(np.asarray(p1.pdf(GRID), dtype=float))
    print(f"{wd.name}: {MODELS[model][1]}, {len(post)} samples; median m1 p(m1) peaks at m1 = "
          f"{GRID[np.argmax(q[1])]:.1f}")

out = Path(a.out)
out.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(cols).to_csv(out, sep="\t", index=False, float_format="%.6g")
print("wrote", out)

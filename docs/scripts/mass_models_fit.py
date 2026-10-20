"""Primary-mass distribution of fitted hubble_constant mass models (Power Law + Peak, Multi Peak, FullPop-4.0),
with icarogw.

For each finished `hubble_constant --method spectral` work directory given, evaluates the source-frame
primary-mass density p(m1 | Λ) of its mass model for each posterior sample Λ (posterior_reweighted.tsv, else
posterior.tsv), with icarogw's own model, and writes the median and the 90% band of m1·p(m1) (the density per
ln m1), and the model at the posterior median of each parameter, on a grid of m1:
docs/scripts/data/mass_models_fit.tsv, read by mass_models_plp_mltp.py.

- Power Law + Peak and Multi Peak: the primary-mass prior smoothed at m_min, as in the likelihood.
- FullPop-4.0: a joint density p(m1) p(m2) f(m1, m2) with a pairing function f; p(m1) is its integral over m2 < m1.

Runs with the Python of the icarogw environment, from anywhere:
    /path/to/icarogw-env/bin/python docs/scripts/mass_models_fit.py WORKDIR [WORKDIR ...] [--out FILE]
"""
import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np

MODELS = {"plp": ("massprior_PowerLawPeak", "Power Law + Peak"), "mltp": ("massprior_MultiPeak", "Multi Peak"),
          "fullpop4": ("m1m2_paired_massratio_bplmulti_dip", "FullPop-4.0")}
GRID = np.logspace(0.0, np.log10(250.0), 700)

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

M1, M2 = (x.ravel() for x in np.meshgrid(GRID, GRID, indexing="ij"))


def density(model: str, wrapper, row) -> np.ndarray:
    """m1 p(m1) on GRID for the parameters in `row`."""
    wrapper.update(**{k: float(row[k]) for k in wrapper.population_parameters})
    if model == "fullpop4":
        base = np.nan_to_num(np.asarray(wrapper.prior.pdf_base.pdf(GRID), dtype=float))
        pair = np.asarray(wrapper.prior.pairing_function(M1, M2), dtype=float).reshape(GRID.size, GRID.size)
        p1 = base * np.trapezoid(pair * base[None, :], GRID, axis=1)
    else:
        p1 = np.nan_to_num(np.asarray(w.LowpassSmoothedProb(wrapper.prior, float(row["delta_m"])).pdf(GRID),
                                      dtype=float))
    return GRID * p1 / np.trapezoid(p1, GRID)


cols = {"m1": GRID}
for wd in map(Path, a.workdirs):
    model = json.loads((wd / "run_settings.json").read_text())["mass_model"]
    if model not in MODELS:
        raise SystemExit(f"{wd}: mass model {model!r} not handled ({', '.join(MODELS)})")
    post = pd.read_csv(wd / ("posterior_reweighted.tsv" if (wd / "posterior_reweighted.tsv").exists()
                             else "posterior.tsv"), sep="\t")
    if len(post) > a.max_samples:
        post = post.sample(a.max_samples, random_state=1)
    wrapper = getattr(w, MODELS[model][0])()
    dens = np.array([density(model, wrapper, row) for _, row in post.iterrows()])
    q = np.percentile(dens, [5, 50, 95], axis=0)
    for name, v in zip(("lo90", "median", "hi90"), q):
        cols[f"{model}_{name}"] = v
    cols[f"{model}_medpar"] = density(model, wrapper, post.median())   # a typical fit
    print(f"{wd.name}: {MODELS[model][1]}, {len(post)} samples; median m1 p(m1) peaks at m1 = "
          f"{GRID[np.argmax(q[1])]:.1f}")

out = Path(a.out)
out.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(cols).to_csv(out, sep="\t", index=False, float_format="%.6g")
print("wrote", out)

# icarogw and gwcosmo

The LVK uses two independent codes for cosmology with gravitational waves:
[icarogw](https://github.com/icarogw-developers/icarogw)
(Mastrogiovanni et al. 2024 [\[57\]](../references.md#ref-57)) and
[gwcosmo](https://git.ligo.org/lscsoft/gwcosmo)
(Gray et al. 2020 [\[24\]](../references.md#ref-24); Gray et al. 2022 [\[25\]](../references.md#ref-25);
Gray et al. 2023 [\[26\]](../references.md#ref-26)). The GWTC-4.0 and GWTC-5.0 cosmology analyses
(LVK 2026 [\[29\]](../references.md#ref-29); LVK 2026 [\[30\]](../references.md#ref-30)) run **both**,
and their published results are an equally weighted mixture of the posterior samples of the two codes
(50% from each), "to incorporate any residual (small) systematic uncertainty associated with differences
in the numerical implementation of the likelihood". `gwtc_analysis` uses icarogw.

## What they share

- **The statistical framework:** the same hierarchical likelihood
  ([Hubble constant](spectral-siren.md#the-hierarchical-likelihood)), with the selection function
  estimated by reweighting the same LVK search-sensitivity injections, under the same criterion of more
  than 4N effective injections (Farr 2019 [\[37\]](../references.md#ref-37)).
- **The methods:** bright sirens (GW170817), dark sirens with a galaxy catalog, and spectral sirens, which
  are the dark-siren method with an empty catalog. The GWTC-4.0 paper shows the spectral-siren results of
  the two codes separately, and they are consistent.
- **The population models** of the papers (Power Law + Peak, Multi Peak, FullPop-4.0) and their priors.
- **The sampler** of the papers: nessai, nested sampling with normalizing flows
  (Williams et al. 2021 [\[61\]](../references.md#ref-61)), run through bilby. `gwtc_analysis` uses
  dynesty (Speagle 2020 [\[60\]](../references.md#ref-60)).

## Where they differ

| | icarogw | gwcosmo |
|---|---|---|
| Per-event term | Monte Carlo sum of p_pop / π_PE over the PE samples, rejected when the effective number of samples is below 10 | kernel density estimate of the redshift distribution of each event in each HEALPix sky pixel, built from the population-weighted PE samples |
| Trade-off | simple and unbiased, but noisy: hence the effective-sample-size checks | more numerically stable; "susceptible to systematic uncertainties if re-weighted sample sizes are too small" (GWTC-4.0 paper) |
| Galaxy catalogs | its own galaxy-catalog classes | precomputed line-of-sight redshift priors per sky pixel (GLADE 2.4 or GLADE+) |
| Interface | a Python library: the analysis is written as a driver (in `gwtc_analysis`, `h0_icarogw.py`) | command-line programs (`gwcosmo_dark_siren_posterior`, `gwcosmo_bright_siren_posterior`) taking JSON dictionaries of PE files and skymaps, an injection file and a dictionary of parameters, each fixed, sampled or on a grid |
| Hardware | CPU or GPU (CuPy) | CPU (multiprocess) |
| Distribution | GitHub; Python ≥ 3.12 | PyPI (`pip install gwcosmo`, version 3.1.0) and git.ligo.org; pinned older dependencies (`numpy<=1.24.2`, `healpy==1.17.3`, `ligo.skymap<=1.0.7`), so it needs an environment of its own |

The two numerical strategies for the per-event term are the main difference. It matters most for events
whose PE samples overlap the population model poorly: in the spectral-siren analyses, the lightest
binary black holes (GW190924 and some O4a events), which come close to icarogw's threshold on the
effective number of PE samples ([Numerical stability](spectral-siren.md#numerical-stability)).

## Consequences for the reproduction

`gwtc_analysis` reproduces the published results with icarogw alone, while the papers combine both
codes. The two codes agree within their Monte Carlo uncertainties, so the difference is small, but it
adds to the other reasons why an exact match is not expected: the number of PE samples per event, the
sampler (dynesty against nessai), and the sampling noise of the runs
([Hubble constant](spectral-siren.md#injection-subsets-and-reweighting)).

gwcosmo could become a second engine of the `hubble_constant` mode, to reproduce the papers' combination.
It would need its own environment, the translation of the prepared events, PE files, injections and
priors into gwcosmo's input format, and the merging of the samples of the two codes. It is not
implemented.

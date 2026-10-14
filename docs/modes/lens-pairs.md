# lens_pairs

A **screen for strongly lensed pairs** among the catalog events. A galaxy close to the line of sight can split a
gravitational wave into two or four images that reach the Earth days to months apart: the network records
the same merger several times, as separate events. This mode tests every pair of the selected binary black
holes for that, as the first (tier-1) stage of the LVK lensing searches does.

```bash
gwtc_analysis lens_pairs                                  # the O4a events (GWTC-4.0), as the LVK search
gwtc_analysis lens_pairs --catalogs GWTC-2.1 GWTC-3 GWTC-4 # pairs across O1-O4a
gwtc_analysis lens_pairs --catalogs ALL --data-repo s3     # every run, PE files from the faster S3 mirror
```

It ranks pairs. It does not confirm lensing, which needs a joint parameter estimation of the two signals
(the LVK tier-2 and tier-3 analyses, GOLUM and Hanabi).

## Why lensed copies look alike

A lens changes neither the frequencies nor the frequency evolution of the signal (geometric optics), so the
**detector-frame masses and the spins** of the images are the same, and the images are about an arcsecond
apart, so their **sky position** is the same within the GW localization. What differs is the amplitude
(each image has its own magnification μ, so its apparent distance is off by √μ), the arrival time, and for a
saddle-point image a phase shift of π/2. See [How the redshift is obtained](../science/redshift.md) for the
(1 + z) degeneracy that hides the true masses.

## Method

**Posterior overlap** (Haris et al. 2018 [\[108\]](../references.md#ref-108)). For two events with
posteriors \(p_1, p_2\) and the PE prior \(\pi\),

\[
\mathcal{B}_\text{overlap} = \int \frac{p_1(\theta)\, p_2(\theta)}{\pi(\theta)}\, d\theta ,
\]

large when the two posteriors sit at the same place compared with the prior volume. θ is the set the LVK
uses (GWTC-4.0 lensing [\[109\]](../references.md#ref-109), App. A.1): the detector-frame component masses,
the spin magnitudes, the cosines of the spin tilts and the sky position.

- The integral is factorized into an **intrinsic** part (masses and spins) and a **sky** part.
- *Intrinsic:* a Gaussian kernel density of each event in the unbounded coordinates
  (ln 𝓜, logit q, logit a₁/0.99, logit a₂/0.99, artanh cos θ₁, artanh cos θ₂), and the Monte Carlo average of
  \(\hat p_2/\pi\) on the samples of event 1 and of \(\hat p_1/\pi\) on those of event 2, averaged.
- *Sky:* HEALPix histograms (nside 64) of the samples, smoothed by a 2° beam; \(\mathcal{B}_\text{sky} =
  4\pi \sum P_1 P_2 / A_\text{pix}\).
- *Prior:* the LVK PE priors are uniform in the component masses (bilby `UniformInComponents` chirp mass and
  mass ratio), uniform in the spin magnitudes up to 0.99, isotropic in the tilts and on the sky. π is that
  prior over a fixed mass box (1–1000 M☉), so \(\log \mathcal{B}_\text{overlap}\) carries a constant offset,
  the same for every pair, which changes no ranking.
- Pairs whose detector-frame chirp masses cannot overlap (99.8 % intervals disjoint, 5 % margin) get
  log₁₀ 𝓑 = −30, as in the LVK results.
- PE samples: 2000 per event (`--pe-samples`) of the IMRPhenomXPHM analysis, as `hubble_constant`.

**Time delay.** The statistic is \(\mathcal{B}_\text{overlap}\,R\), with \(R = p_L(\Delta t)/p_U(\Delta t)\)
(Haris et al. 2018; the LVK \(\mathcal{R}^L_U\)):

- \(p_U\): the distribution of the delays of all the pairs of the selection. It follows the runs, the gaps
  between them and the duty cycle.
- \(p_L\): galaxy lenses as singular isothermal spheres with the velocity-dispersion function of SDSS
  early-type galaxies (Choi et al. 2007 [\[111\]](../references.md#ref-111)), weighted by their
  cross-section (∝ σ⁴) and optical depth along the line of sight of each event (source redshift: the catalog
  median), source position uniform in the Einstein disk, and both images detectable (magnification bias
  \(\mu_-^{3/2}\) on the fainter one). The delay of the two images is
  \(\Delta t = (1+z_l)\,\frac{D_l D_s}{c\,D_{ls}}\,2y\,\theta_E^2\), \(\theta_E = 4\pi(\sigma/c)^2 D_{ls}/D_s\).
  Typical delays are days to a few months; pairs a year or more apart are disfavoured.
- Both densities keep 1 % uniform in log Δt (1.4 minutes to 27 years): where both vanish, a ratio of kernel
  tails would be arbitrary; this bounds |log₁₀ R| to about 2.

**False-positive probability.** The FPP of a pair is the fraction of the pairs of the selection with a higher
statistic: the selection is its own background, valid while almost every pair is unrelated (the expected
lensed fraction is ~10⁻³). The LVK computes it on simulated unlensed events instead. With a 1 % threshold,
unrelated events alone flag about 1 % of the pairs, so **a flagged pair is a candidate for joint analysis,
not a detection**; a lensed pair would show as an outlier at the top of the ranking.

## Validation: the O4a search of the LVK

The GWTC-4.0 lensing paper [\[109\]](../references.md#ref-109) tested every pair of the 84 O4a binary black holes
with FAR < 1 /yr, **3486 pairs**, and its data release ([Zenodo 18163632](https://zenodo.org/records/18163632)) gives the
posterior-overlap results of each pair. `gwtc_analysis lens_pairs` (default catalog GWTC-4) selects the same
84 events.

Run of 14 October 2026 (PE files from the S3 mirror, 2000 samples per event), compared pair by pair with the
LVK results (superevent names mapped to the catalog events through their GraceDB times):

| | LVK (posterior overlap) | `lens_pairs` |
|---|---|---|
| Events, pairs | 84, 3486 | 84, 3486 (the same) |
| Pairs without chirp-mass overlap (log₁₀ 𝓑 = −30) | 2270 | 1803 (1777 in common) |
| Pairs flagged at FPP ≤ 1 % | 38 | 34, **28 of them among the LVK's 38** |
| Highest-ranked pair | S231113bw + S231223bg (GW231113_200417 + GW231223_202619) | the same |
| The LVK's 5 highest pairs | ranks 1–5 | ranks 1, 2, 4, 5, 6 |

- **Overlap Bayes factor.** Over the 429 pairs with LVK log₁₀ 𝓑_overlap > 0, ours equals the LVK's plus a
  constant 1.1 (the normalization of the reference prior, the same for every pair), with a scatter of 0.2
  (robust standard deviation); Spearman rank correlation 0.89. Further down the tail the two differ more, where
  neither matters (pairs the LVK puts at −30 are at a median −3.5 here).
- **Time-delay factor.** log₁₀ R against the LVK's: Spearman 0.998, Pearson 0.99, median offset −0.14: the
  isothermal-lens model reproduces their galaxy-lens delay distribution.
- **Differences.** The FPP here is the rank among the pairs themselves, the LVK's comes from simulated unlensed
  events, which is why 34 pairs are flagged rather than 38. The largest single discrepancy, S231018by +
  S231118an (LVK rank 6, rank 84 here), comes from the sky term (ln 𝓑_sky = −0.8): our sky densities are
  smoothed histograms of 2000 samples of an earlier PE release version. The LVK passed 105 pairs (posterior
  overlap or Phazap) to tier 2; none survived tier 3.

## Outputs

- `lens_pairs.tsv`: every pair, ranked: the two events, the delay (days), ln 𝓑 intrinsic and sky,
  log₁₀ 𝓑_overlap, log₁₀ R, log₁₀ of the statistic and the FPP; `lens_pairs.events.tsv`: the events.
- `lens_pairs.html`: the distribution of the statistic with the FPP threshold, the time-delay densities and R,
  and the highest-ranked pairs.

## Options

| Option | Default | |
|---|---|---|
| `--catalogs` | `GWTC-4` | Catalogs whose events are paired; several, or `ALL`. With an update (GWTC-4.1) its events replace those of the catalog it re-analyses. |
| `--far-threshold` | 1 | FAR at most this, per year (the LVK lensing selection). |
| `--min-mass` | 3 | Both source-frame masses above it: binary black holes. |
| `--data-repo` | `zenodo` | `s3`: the MinIO mirror of the releases, much faster; it may hold an earlier version of a release than the latest on Zenodo. |
| `--pe-samples` | 2000 | PE samples per event for the kernel densities. |
| `--fpp-threshold` | 0.01 | Pairs flagged in the report. |

## Limits

- Tier 1 only: no joint parameter estimation, and no phase consistency test (the LVK runs Phazap beside
  posterior overlap; a type-II image's π/2 phase is not used here).
- Only events with PE: LVK marginal candidates have none. Faint second images below the detection threshold
  need a targeted search of the strain data (the LVK sub-threshold search), outside this mode.
- Events of other catalogs (IAS, OGC), which the O1–O4a search of Barsode et al.
  [\[110\]](../references.md#ref-110) includes, are not read.
- The lens model is a single isothermal sphere: no ellipticity (four images), no clusters (delays of years).

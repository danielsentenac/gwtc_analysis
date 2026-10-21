# GWTC Analysis

**GWTC Analysis** (`gwtc_analysis`) is a command-line suite for exploring the public
**Gravitational-Wave Transient Catalogs** (GWTC) of the **LIGO–Virgo–KAGRA** (LVK) collaboration,
from GWTC-1 (2015) to GWTC-5.0 (the O4b observing run, 2025).

<!-- CATALOG_COVERAGE_BEGIN -->
> **Catalogs in version 0.8.2**: GWTC-1, GWTC-2.1, GWTC-3, GWTC-4.0, GWTC-5.0, and the update GWTC-4.1 (of GWTC-4.0), observing runs O1 to O4b. **Latest catalog: GWTC-5.0 (O4b)**. Registry checked against GWOSC and Zenodo on 2026-10-05; catalogs published later need a newer version of the package (`gwtc_analysis check_catalogs` tells whether GWOSC has published one).
<!-- CATALOG_COVERAGE_END -->

It reads the catalogs directly from their public releases (GWOSC, Zenodo, an S3 mirror or Galaxy
collections) and produces TSV tables, plots and self-contained HTML reports.

With it, a user can visualize the content of the catalogs (the masses, distances, spins and sky localizations of
the events, and the posteriors and strain data of each one), and illustrate and test the published results by
running the analyses behind them on the same public data: the Hubble constant from spectral, dark and bright
sirens and their combination, the localization of a transient event and of its electromagnetic counterpart on the
sky and in distance, the merger rates of each source class, the spin population of binary black holes, and tests
of the black-hole area law, of strong lensing and of the neutron-star equation of state. Each analysis is compared
with the published result it reproduces, and the comparison is documented. The same analyses extend to the
catalogs that follow, once added to the package registry (`check_catalogs` reports a new GWOSC release,
[Adding a new catalog](data-sources.md#adding-a-new-catalog)), and to new datasets, such as deeper galaxy
catalogs for the dark sirens (DES, Rubin), so that the published results can be redone and extended as the data
grow.

## What it does

| Mode | Purpose |
|---|---|
| [`catalog_statistics`](modes/catalog-statistics.md) | Catalog-wide tables and plots: masses, distances, SNR, detector networks, sky-localization areas, remnants (radiated energy, final spin) |
| [`event_selection`](modes/event-selection.md) | Events selected by mass, distance and χ_eff ranges, or by class: neutron stars, lower mass gap, hierarchical-merger candidates |
| [`search_skymaps`](modes/search-skymaps.md) | Events whose sky localization contains a given sky position; with `--counterparts`, gamma-ray bursts and neutrinos coincident with the events |
| [`parameters_estimation`](modes/parameters-estimation.md) | Posterior plots of one event, whitened strain overlays, q-transforms, matched-filter SNR, 3D sky map and host-galaxy candidates |
| [`build_unofficial_pe`](modes/unofficial-pe.md) | A PESummary-compatible bundle for GW170817, rebuilt from public GWTC-1 products |
| [`rates`](modes/rates.md) | BNS, NSBH and BBH merger rates, R = N / ⟨VT⟩, from the catalogs and the LVK sensitivity injections |
| [`hubble_constant`](modes/hubble-constant.md) | The Hubble constant by `--method`: the compact-binary mass spectrum, BBH only or all binaries with FullPop-4.0 (spectral siren) or with a galaxy catalog (dark siren), with icarogw; an event with an identified host ([bright siren](modes/hubble-constant-bright.md): GW170817, the candidate GW190521 flare); the joint posterior of independent results |
| [`galaxy_catalog`](modes/galaxy-catalog.md) | The galaxy catalog of the dark siren: GLADE+ K band, or any catalog read in chunks, turned into icarogw's line-of-sight redshift prior; local or Slurm (cluster) execution |
| [`counterpart`](modes/counterpart.md) | Where an electromagnetic counterpart sits in the GW posterior: searched probability on the sky, distance along its line of sight against the host redshift, viewing angle and an independent constraint on it |
| [`area_law`](modes/area-law.md) | Hawking's area law with GW250114: initial areas from the inspiral against the remnant area from the ringdown, reproducing the published 4.4σ |
| [`scalar_cloud`](modes/scalar-cloud.md) | A binary black hole in an ultralight scalar-field cloud: the IMRPhenomXPHM_Scalar waveform, best fit against vacuum calibrated on off-source noise, and the Bayes factor (GW190728) |
| [`continuous_waves`](modes/continuous-waves.md) | Continuous waves: the LVK upper limits on known pulsars, unknown neutron stars and boson clouds around black holes, with the boson mass of GW190728 at its continuous-wave frequency |
| [`stochastic`](modes/stochastic.md) | The background of unresolved compact binaries, Ω_GW(f), from the spectral-siren masses and rate evolution and the `rates` local rates, against the upper limits |
| [`primordial_black_holes`](modes/primordial-black-holes.md) | Upper limits on the fraction of dark matter in primordial black holes of one mass, 0.35–150 M☉: the LVK sub-solar searches (O4a + O4b) and the catalog injections and events |
| [`neutron_star_eos`](modes/neutron-star-eos.md) | The neutron-star equation of state from GW170817 and GW190425 jointly: Λ₁.₄ and R₁.₄ |
| [`spin_population`](modes/spin-population.md) | The BBH effective-spin distribution and its correlation with the mass ratio, corrected for selection effects |
| [`lens_pairs`](modes/lens-pairs.md) | Strongly lensed pairs: every pair of BBH events ranked by posterior overlap and time delay, the tier-1 screen of the LVK lensing searches |
| `zenodo_releases` | The versions of the Zenodo catalog releases ([Data sources](data-sources.md#zenodo-release-versions)) |
| `check_catalogs` | Compares the catalogs and Zenodo records published by GWOSC with the package registry, and reports the new ones ([Adding a new catalog](data-sources.md#adding-a-new-catalog)) |

## Highlights

- **Hubble constant from black holes alone.** The `hubble_constant` mode reproduces the
  *Power Law + Peak* spectral-siren measurement of the GWTC-4.0 cosmology paper [\[29\]](references.md#ref-29):
  H₀ = 107.6 (+45.6 / −35.2) km/s/Mpc against the published 105.5 (+46.4 / −35.8), with the same
  137 binary black holes, PE samples, injections and priors. With FullPop-4.0, the paper's fiducial population of
  all 141 compact binaries, it gives 75.1 (+25.8 / −19.2) against the published 72.9 (+21.9 / −18.8). See
  [Hubble constant (spectral siren)](science/spectral-siren.md) and the [mass models](modes/hubble-constant.md#mass-models).
- **Dark sirens with a galaxy catalog.** The `galaxy_catalog` mode builds the GLADE+ K-band catalog of the
  GWTC-4.0 analysis (1.0 million galaxies, about 45 min as Slurm jobs on a cluster); with it,
  `hubble_constant --method dark` gives H₀ = 114.6 (+41.5 / −33.8) km/s/Mpc against the published
  115.4 (+40.1 / −33.8) with Power Law + Peak, 86.8 (+39.8 / −26.1) against 86.3 (+41.3 / −26.3) with Multi Peak,
  and 82.0 (+21.0 / −16.4) against 82.0 (+21.5 / −16.0) with FullPop-4.0 (the icarogw values of the data release),
  with 5000 PE samples per event. Combined with GW170817, the FullPop-4.0 dark siren gives 76.7 (+13.4 / −8.9)
  against the paper's headline 75.4 (+12.8 / −9.1), which also includes gwcosmo. Deeper catalogs (DES, Rubin) go
  through the same pipeline. See [galaxy_catalog](modes/galaxy-catalog.md) and the
  [validations](modes/hubble-constant.md#validation-gwtc-40-power-law-peak-glade-k-band).
- **Strongly lensed pairs.** `lens_pairs` repeats the first stage of the LVK lensing search on the 84 O4a binary
  black holes (3486 pairs): the same top pair, and 28 of the LVK's 38 pairs flagged at FPP < 1 % (ours ≤ 1 %), compared pair
  by pair with their data release. See [lens_pairs](modes/lens-pairs.md).
- **Merger rates** of the three source classes, corrected for selection effects with the LVK
  search-sensitivity injections. See [Merger rates](science/merger-rates.md).
- **Robust access to the releases:** Zenodo versions resolved at run time, and missing PSDs of
  official releases filled from the events' discovery releases. See
  [Known issues in the public releases](known-data-issues.md).

![Primary masses of the confident GW events](img/mass_distribution.png)

![Primary masses corrected for selection effects](img/mass_distribution_selection_corrected.png)

*Top: the primary masses as detected. Bottom: the same events divided by the sensitive volume-time ⟨VT⟩ of each
mass, from the LVK injections, as the [rates](modes/rates.md) mode computes it: a merger rate per ln m₁. Heavy
binaries are seen much farther, so the 35 M☉ peak of the detections is mostly a selection effect, while the
10 M☉ peak remains. This is not a population fit: m₂ uniform in [1, m₁], no redshift evolution. The two figures are made by the
scripts of [`docs/scripts`](https://github.com/danielsentenac/gwtc_analysis/tree/main/docs/scripts).*

## Where to go next

- New users: [Installation](installation.md), then [Quick start](quickstart.md).
- The physics and statistics behind the rates and H₀ modes:
  [What a GW signal measures](science/gw-signals.md) and the following pages.
- Every option of every mode: [CLI reference](cli-reference.md).
- All the papers the tool and its methods rely on: [References](references.md).

The source code is on [GitHub](https://github.com/danielsentenac/gwtc_analysis) under the MIT
license.

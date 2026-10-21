# Changelog

## Unreleased

### New

- **primordial_black_holes**: upper limits (90%) on the fraction of dark matter in primordial black holes of one
  mass, for early binaries (Hütsi et al. 2021 rate, suppression S₁ S₂ or, with `--suppression lvk`, the fit of the
  LVK sub-solar papers) and late binaries in clusters (Phukon et al. 2021). From 0.35 to 1 M☉, the sensitive
  volume-time of the LVK sub-solar searches, O4a + O4b (Zenodo 23346454): with `--suppression lvk` it reproduces
  the published 0.018–0.032 at 1 M☉ and 0.10–0.18 at 0.35 M☉. Above 1 M☉, the search-sensitivity injections of
  `rates` and the events matched by chirp mass: ~10⁻³ at 10–40 M☉, ~10⁻⁴ near 100 M☉. The events of `rates` now
  carry the GWOSC chirp mass.

### Changes

- **scalar_cloud**: the best-fit statistic is renamed the **deviance** 𝒟 (the likelihood-ratio test statistic),
  formerly Λ, which the population analyses use for their hyperparameters. Only the text changes: the result
  files keep the key `Lambda`.
- Documentation: the FullPop-4.0 dark siren with GLADE+ K band, 82.0 (+21.0 / −16.4) km/s/Mpc against 82.0
  (+21.5 / −16.0) in the GWTC-4.0 release (icarogw; 79.6 in the paper, with gwcosmo), and its combination with
  GW170817, 76.7 (+13.4 / −8.9) against the paper's headline 75.4 (+12.8 / −9.1). The spectral sirens are also
  compared with the icarogw release values. Results of diagnostics made outside the `gwtc_analysis` command are
  labelled as such.

### Fixes

- **hubble_constant --method joint**: the number of events of a spectral or dark siren counted the events with a
  common name twice (160 instead of 141 for the FullPop-4.0 spectral siren); it is now the number of rows of its
  `events.tsv`. The posteriors were not affected.

## 0.8.2 (2026-10-16)

### Fixes

- **hubble_constant report**: it said "N BBH events" with `--mass-model fullpop4`, which keeps neutron-star
  binaries too; it now gives the number per class (BBH, NSBH, BNS) when the events are not all BBHs.
- **Help texts**: `search_skymaps --skymap-label` gives the full fallback order of the maps and says that `any`
  acts as `Mixed` with `--counterparts`; `--data-repo galaxy` says that a staged collection must hold the maps.
- **Docstrings** of `hubble_constant`, `h0_icarogw` and `stochastic` brought up to date (dark siren,
  FullPop-4.0, `--data-repo`, the 2025 upper limit).

## 0.8.1 (2026-10-16)

### Fixes

- **Containers with an arbitrary uid** (`docker run -u`, Galaxy, Apptainer): building the GW170817 PE bundle
  (`build_unofficial_pe`, and the modes run on GW170817: bright-siren `hubble_constant`, `counterpart`,
  `neutron_star_eos`, `parameters_estimation`) failed with
  `OSError: No username set in the environment`. pesummary records the user name with `getpass.getuser()`,
  which raises when the uid has no passwd entry; `USER` is now set when the name is unknown.
- **Docker image**: the caches under `HOME=/tmp` were writable by root only, so under `docker run -u` pesummary
  stopped with `PermissionError: /tmp/.cache/pesummary` and matplotlib warned at every start. They are now
  writable by any uid, and the image sets `USER` and `MPLCONFIGDIR`. The image now also includes `pycbc`:
  without it, the GW170817 runs could not fit the time, phase and polarization of the maximum-likelihood sample
  to the strain, so the waveform overlay was not coherent (the results were unchanged).
- **Galaxy wrapper**: exports `USER`, which the container of the tools-iuc CI lacks (bright-siren test).
- **Network failures are retried** (new module `net`): a connection error or a timeout (such as a DNS failure),
  or a status sent by an overloaded or throttling server (429, 500, 502, 503, 504), is retried twice, after the
  delay the server asks for (Retry-After, at most 60 s), else after 5 s then 15 s. Other HTTP errors (404, 403)
  are raised at once. This covers the downloads of public files (DCC, Zenodo, arXiv) and of PE files (Zenodo,
  Galaxy), and the Zenodo, Galaxy, AMON and HEASARC queries. The GWOSC and VizieR queries already retried.
  In the tools-iuc CI, a DNS failure on dcc.ligo.org and a 503 from the Zenodo API each failed one test.
- **Backup of the public files**: when a public file (DCC, Zenodo, arXiv) still cannot be downloaded after the
  retries, its copy in the public usegalaxy.org GWTC history is read instead, whatever `--data-repo` is. The
  history now holds the 21 public files the modes use (12 from DCC, 8 from Zenodo, 1 from arXiv), each named
  after its mirror key (`public/dcc/...`). The report lists such a file as "Galaxy history (<host> failed)".
- **Docs**: the Docker example runs as you (`-u $(id -u):$(id -g)`), so the output files belong to you.

## 0.8.0 (2026-10-15)

### Hubble constant: one mode, four methods

`hubble_constant` now holds every way of measuring H₀, chosen with `--method`:

- `spectral` (default): the BBH mass spectrum, as before.
- `dark`: the mass spectrum and a galaxy catalog (`--galaxy-catalog`, built by the new `galaxy_catalog` mode).
  Later stages take the method from the work directory, without `--galaxy-catalog` again.
- `bright`: an event with an identified host (GW170817 and NGC 4993, the candidate GW190521 flare), in seconds.
  New `--viewing-angle MEAN SIGMA`: an independent constraint on the viewing angle (e.g. the radio jet of
  GW170817), applied as weights on the PE samples. The result is a work directory (`posterior_grid.tsv`,
  `bright.json`).
- `joint`: the product of independent results (`--inputs DIR DIR ...`: work directories of the other methods, or
  posterior TSVs), with the table of each input and of the joint posterior. The inputs are checked to be
  independent: at most one spectral or dark siren, no event in two inputs.

The options of another method are refused. The default work directory and outputs are named after the method:
`hubble_constant_<method>`, `hubble_constant_<method>.html`, `hubble_constant_<method>.tsv` (they were
`hubble_constant_run`, `hubble_constant.html`, `hubble_constant.tsv`).

### New mode: `counterpart`

Where an electromagnetic counterpart sits in the GW posterior of an event: the searched probability of its position
(from the LVK sky map, else the PE samples), the distance along its line of sight against the distance of the host
redshift for Planck and SH0ES H₀, the viewing angle and the distance–inclination degeneracy, and an optional
viewing-angle constraint. `--ra`/`--dec` test another position.

### New mode: `lens_pairs`

A screen for strongly lensed pairs: every pair of the BBH events of the chosen LVK catalogs (`--catalogs`,
default GWTC-4) is ranked by the posterior-overlap Bayes factor (detector-frame masses, spins, sky position)
times the ratio of time-delay densities of galaxy-lensed and unrelated pairs (isothermal lenses), with a
false-positive probability from the ranking of all the pairs. On the 84 O4a events it reproduces the LVK
GWTC-4.0 search: the same 3486 pairs, the same top pair, 28 of the LVK's 38 flagged pairs. `--data-repo s3`
reads the PE files from the MinIO mirror.

### New mode: `continuous_waves`

The upper limits of the LVK continuous-wave data releases (LIGO DCC): the 45 known pulsars of O4a against their
spin-down limits, the all-sky O4a sensitivity, and the O3 boson-cloud exclusions, with the GW190728 boson mass of
`scalar_cloud` (0.70–1.17 × 10⁻¹² eV) at its continuous-wave frequency (340–570 Hz). The `scalar_cloud` report gives
that band. The narrowband upper limits of O4a (16 pulsars) are on the pulsar plot.

### Multi-messenger coincidences

`search_skymaps --counterparts grb neutrino`: the gamma-ray bursts (Fermi-GBM, Swift) and neutrinos (IceCube alert
tracks from June 2019, point-source events of O1 and O2) within a time window of each GW merger, scored by their
overlap with the GW map, their chance probability and RAVEN's coincidence FAR; `--skymap FILE` for given maps (a
GraceDB alert). Over the 280 catalog maps and GW170817, the only coincidence is GRB 170817A; for a gamma-ray burst, the speed of gravity it implies (the distance along its line of sight from the 3D
sky map): −2.8 × 10⁻¹⁵ to +7.2 × 10⁻¹⁶ for GW170817, as published.
`hubble_constant --method bright` and `counterpart` take any event of the GWTC PE releases with its counterpart
given on the command line (`--ra`, `--dec`, `--redshift`, `--host`, `--counterpart-class`), e.g. the afterglow
and host of a coincident GRB.

### New mode: `scalar_cloud`

Is a binary black hole better fitted with an ultralight scalar cloud around it? The IMRPhenomXPHM_Scalar waveform of
Roy et al. (arXiv:2510.17967), IMRPhenomXPHM plus the dephasing of the angular momentum lost to the cloud, against
the vacuum IMRPhenomXPHM of the catalog, on GWTC PE settings and GWOSC strain. Stages `prepare`, `fit`,
`background`, `sample` (nested sampling with the paper's agnostic and superradiance priors: the Bayes factor),
`seeded` (the best-fit statistic Λ with the cloud fits started in the cloud mode found by the sampling, on the event
and on off-source injections: its p-value), `report`; `--executor slurm` for the cluster stages. For GW190728:
Λ = 8.98, p ≈ 0.01 against 150 vacuum injections (and 4.7% of the injections of the cloud best fit at
or below it), ln B = 1.34 with the superradiance prior (the paper finds 3.5), −0.82 with the agnostic prior up to 10⁹ g/cm³
(0.34 up to 3 × 10⁸, the paper finds 0.4); the
improvement is concentrated near 30 Hz, mostly in L1.

### Removed: `bright_siren`

Replaced by `hubble_constant --method bright` (H₀), `hubble_constant --method joint` (combination with the spectral
or dark siren, formerly `--spectral-posterior`) and `counterpart` (degeneracy plot and table). `--src-name` is now
`--event`; `--plots-dir` is replaced by `<workdir>/plots`. The module `bright_siren.py` is split into `h0_bright.py`,
`h0_joint.py` and `counterpart.py`.

### Dark sirens

- `galaxy_catalog` mode: an icarogw line-of-sight galaxy catalog (GLADE+ K band by default, or any catalog read in
  chunks), local or as Slurm jobs; GLADE+ selection options (`--glade-types`, `--glade-redshift`, `--glade-sigmaz`,
  `--sigmaz`, `--where`, `--zmin`, `--ptype`).
- `hubble_constant` dark siren with icarogw's `CBC_catalog_vanilla_rate`; Slurm executor (`--executor slurm`,
  chains of batch jobs); likelihood thresholds `--neff-pe`, `--neff-inj`; settings files (`--settings`); the
  resolved options recorded in `<workdir>/options_<mode>.json`.
- Validation: the GWTC-4.0 Power Law + Peak dark siren with GLADE+, 114.6 (+41.5 / −33.8) km/s/Mpc with 5000 PE
  samples per event, against 115.4 (+40.1 / −33.8) published; Multi Peak, 86.8 (+39.8 / −26.1) against
  86.3 (+41.3 / −26.3); with GW170817, a median of 77.5 against the paper's 75.4 (76.7 with
  `--population fullpop4`).
- `hubble_constant --method bright --population fullpop4`: the FullPop-4.0 mass distribution and the
  Madau–Dickinson rate at the GWTC-4.0 medians, in the numerator and the selection term, as the GWTC-4.0
  reanalysis of GW170817 (78.5 against the paper's 77.8; 79.9 by default).
- `--pe-distance-prior dl2|comoving|source-frame` (`hubble_constant --method bright`, `counterpart`): the PE distance
  prior, instead of the recorded one or the catalog default.
- `--catalogs` outside the observing runs of the injection release is refused (it was silently restricted).

### Data sources

- `--data-repo zenodo|s3|galaxy` in every mode that downloads files: Zenodo or the public source (LIGO DCC), the S3
  (MinIO) mirror, or the public usegalaxy.org history. The mirrors serve the same files, those of the selected Zenodo
  release version: a file a mirror lacks, or a failed transfer, comes from the public source.
- `--zenodo-version` applies with any `--data-repo`.
- The PE file depends only on the options: one rule per mode (nocosmo files up to GWTC-3 for the population and
  H₀ modes, mixed_cosmo for `parameters_estimation`, the combined file from GWTC-4.0).
- Every report lists the PE files and public data files it used, and ends with the options of its run, defaults
  included.

### Galaxy

- The Galaxy wrapper is generated from the command line (`galaxy/gen_galaxy_wrapper.py`): 18 analyses (every mode
  except `galaxy_catalog`, whose icarogw is not on conda; `hubble_constant` as its bright and joint methods), with
  10 tests that pass under planemo.

### FullPop-4.0

- `hubble_constant --mass-model fullpop4`: FullPop-4.0 spectral and dark sirens (all compact binaries). GWTC-4.0
  spectral siren: 74.4 (+25.4 / −18.9) km/s/Mpc against the published 72.9 (+21.9 / −18.8).
- `rates --mass-model fullpop4`: one rate for all compact binaries, split into BNS, NSBH and BBH.

### Other changes

- `gwtc_analysis -h` and `gwtc_analysis MODE -h` answer in about 0.2 s instead of 10–15 s: the parser takes its
  defaults from the light `cli_defaults` module, and each mode imports its analysis code only when it runs.
- `rates`: Multi Peak mass model, BBH rate over the draws of a `hubble_constant` posterior
  (`--population-posterior`), evolving BBH rate also at z = 0.
- `stochastic`: Multi Peak mass model, read from the `hubble_constant` posterior.
- `parameters_estimation`: 3D sky map from the release FITS archive, with GLADE+ host-galaxy candidates.
- `search_skymaps`, `skymap3d`: one sky map per event, the requested waveform label or a shared fallback.
- New modes `neutron_star_eos` (joint equation of state from GW170817 and GW190425) and `spin_population` (BBH
  effective-spin population).
- Documentation: Waveform models page, the stochastic background of compact binaries, dark sirens, the bright siren
  and the counterpart pages.

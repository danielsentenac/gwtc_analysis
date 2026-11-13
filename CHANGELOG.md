# Changelog

## Unreleased

### New

- **primordial_black_holes**: unequal masses. The LVK map F_PBH(m₁, m₂), one upper limit per pair of masses, as
  defined by the script of the sub-solar release (its (m₁, m₂) tables below m₂ = 1.1 M☉, reproduced exactly; the
  search-sensitivity injections in log bins above), and the limits on the total f_PBH of log-normal mass functions
  (Hütsi et al. 2021, eq. 3.6; S₁ for the mass function), summed over all pairs and tested in each chirp-mass window
  as the monochromatic limits, with the f at which a broad mass function is allowed again when the suppression makes
  the count fall. New options `--lognormal-widths` (0.5 1 2), `--lognormal-masses`, `--no-extended`; tables
  `<out-limits>.map.tsv` and `<out-limits>.lognormal.tsv`, plots `f_pbh_map.png` and `f_pbh_lognormal.png`. The
  default run takes about 5 min instead of 2.
- **primordial_black_holes**: the region f_PBH > 1 of the plot is shaded as "no constraint".
- **event_selection**: `--area-max` keeps the events whose 90% sky area (deg²) is at most a value, computed from
  the skymaps of the events passing the other cuts (`--area-cred` for another level, `--data-repo`,
  `--zenodo-version` as catalog_statistics); the area is added to the TSV.
- **event_selection**: cuts on the mass ratio (`--q-min`, `--q-max`, q = m₂/m₁ of the GWOSC medians), the network
  SNR (`--snr-min`), the significance (`--p-astro-min`, `--far-max` per year) and the observing run (`--runs O3a
  O4a`, from the GPS time); each adds its column to the TSV. GW241011-like hierarchical candidates in one command:
  `--preset hierarchical --pisn-gap-min 15 --q-max 0.35`.
- **event_selection**: `--detectors H L V` keeps the events by their detector network (the one GWOSC lists for
  the event, asked for the events passing the other cuts and added to the TSV as `detectors`);
  `--detectors-match` applies the list as `all` (default: the network holds all of them), `exact` (these and no
  other) or `any` (at least one).

### Changes

- **scalar_cloud** docs: the GW190728 section follows a new full run with the paper's 1000 live points for the
  three samplings: 𝒟 = 9.78 (was 8.98), p = 0.013 against 150 vacuum injections, ln B = 1.06 ± 0.24 with the
  superradiance prior and −1.21 ± 0.25 with the agnostic one (were 1.34 and −0.82); boson mass 0.67–1.17 × 10⁻¹² eV,
  continuous waves at 324–568 Hz; 3.3% of the 150 injections of the cloud best fit at or below the event's 𝒟
  (was 4.7%). The matched-filter figure replaces the whitened strain overlay.
- **continuous_waves**: the default boson band is the one of the new scalar_cloud run of GW190728,
  0.67–0.86–1.17 × 10⁻¹² eV (was 0.70–0.89–1.17): 324–566 Hz, light black holes of 8–30 M☉ at 1 kpc, young clouds
  excluded out to 16–36 kpc.
- **scalar_cloud**: the fits table and the figures of the report use the vacuum fit the deviance is computed
  against (the one refitted with the cloud fits), so that they give the ln L of the "Best fit" line; with the
  seeded fits the table showed the first vacuum fit, 0.03 lower for GW190728.
- **scalar_cloud**: the report draws a matched-filter figure (`plots/matched_filter.png`: the Q-transform of the
  whitened data with the frequency track of each fit, and the matched-filter SNR of the catalog, vacuum and cloud
  fits over the full band and over 20–45 Hz) instead of the whitened strain overlay, on which the fits could not be
  told apart, and the cumulative ln L of the cloud over vacuum against frequency per detector
  (`plots/lnl_by_frequency.png`, with the table by band in the report).
- **catalog_statistics**: the dashed lines of the remnants figure are labelled as reference values (E_rad = 4.8% of
  M and a final spin of 0.686 for equal masses without spin).
- **scalar_cloud**: the mean density of the cloud keeps the bar of the paper's ρ̄_φ in the code: the sampling
  parameter and result key `rho_phi` is now `rhobar_phi` (`rhobar` in the profile points, `rhobar_g_cm3` in the
  report, `rhobar_cgs` in the functions). Work directories and sampling results written before are still read,
  their old keys renamed.
- **search_skymaps, catalog_statistics `--include-area`, event_selection `--area-max`**: the credible level at a
  position and the sky areas are computed on the multi-order pixels of the skymaps instead of the maps flattened to
  their finest resolution (new `skymap_moc.py`). Peak memory over ALL catalogs: 0.8 GB for search_skymaps (was
  15.8 GB) and 0.3 GB for the two others (was 13 GB), 3 to 6 times faster. The areas are unchanged; the credible
  level at a position changes by less than 0.01 percentage point (inside a coarse pixel it is now the level of the
  whole pixel).
- **parameters_estimation**: the PE file is read one analysis at a time (new `pe_lean.py`). pesummary reads the
  samples, noise spectra, sky maps and prior samples of every label of a file, at several times their size, and
  builds the samples of all the labels again at each access; the mode now reads the labels and settings from a
  copy of the file without them, and the samples (the table of the file, as pesummary returns it), the sky map and
  the noise spectra of a label when they are first used. The tables over all the labels (remnant, multipoles) read
  only their columns. Same samples, spectra and maps; a file with another layout is read whole, as before.
- **parameters_estimation, counterpart**: the credible volumes and the galaxies' probability densities of the 3D
  sky map are computed by `skymap3d.volume_stats`, on the voxels and with the formulas of ligo.skymap's
  `crossmatch` (which holds eight arrays of pixels × 999 distance shells, 1.3 GB for 20,000 pixels) but keeping one
  array and its sort order: the same values to 10⁻¹¹. The 3D plot rasterizes one column of the map at a time, and
  counterpart reads the searched probability of the position without the distance columns. Peak memory: 4.7 → 1.7
  GB for GW200115 (ten labels), 3.5 → 1.95 GB for GW190814 with GLADE+, 3.1 → 1.5 GB for a GWTC-4 event,
  1.9 → 0.7 GB for counterpart on GW190521.
- **Galaxy wrapper**: regenerated for this version: the new event_selection options, `--area-max` (missing so
  far), credible levels as numbers between 0 and 1, and without the `--pe-distance-prior` that counterpart refuses.
- **Tests**: a GitHub workflow runs the unit tests (without those that need the network) on each push to `main`
  and on pull requests.
- **Shared caches**: one writer per file. A run that needs a file another run (or thread) is downloading into the
  same cache waits for it and uses it (`net.file_lock`, on `.<name>.lock` next to the file while it is written) instead of overwriting
  its partial file: the PE files, the skymap tarballs and the public files of the mirrors. The reduced sample files
  (hubble_constant, spin_population, lens_pairs, the nside-limited PE copy, the cached strain) are written to one
  partial file per process. The plots of the hits of search_skymaps are written in `--plots-dir` directly, not in
  the working directory first (two runs started in one directory lost each other's plots).
- **GLADE+**: when the CDS VizieR fails, the queries (galaxy_catalog, the host-galaxy candidates of parameters_estimation) go
  to its CfA mirror, and stay there for the rest of the run.

### Fixes

- **catalog_statistics, event_selection** with pandas 3: when GWOSC lists the heavier component second, the 90%
  bounds of the two masses were not swapped with them (both got the bounds of the heavier one). Found by the new
  test workflow; pandas 2 was not affected.
- **primordial_black_holes**: `--no-sub-solar` no longer crashes at the plot (the empty sub-solar table had no
  `pipeline` column).
- **search_skymaps `--prob`, catalog_statistics and event_selection `--area-cred`**: a value outside (0, 1] is
  refused (`--prob 1.5` was accepted and reported every event as inside).
- **scalar_cloud**: the prepare stage read the analysis settings under the names of the GWTC-4 PE files only and
  ended in a `KeyError: segment_length` for a GWTC-3 event (GW191109_010717); the bilby names (`duration`,
  `start_time`, `sampling_frequency`) are read too, and a label without settings or noise spectra (`C01:Mixed`) is
  refused with the labels to choose from.
- **counterpart**: `-h` no longer lists `--pe-distance-prior`, which the mode refuses.
- **Docs, bright siren**: the GRB example for GW190425 needs `--sky-radius 10` (it stopped with "only 119 samples
  within 3.0°").
- **hubble_constant**: where `squeue` exists, the seed lock of a run started outside Slurm on another host was
  taken for the lock of a finished job and removed, so the same seed could run twice at once. `test_seed_lock`,
  which failed inside a Slurm job (the lock then also holds `SLURM_JOB_ID`), covers it.

## 0.9.0 (2026-10-21)

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
- Documentation: primordial black holes and sub-solar masses on their own science page (formation channels, early
  binaries and their suppression, late binaries, the LVK sub-solar bounds of O4a + O4b); the three approaches to
  waveform models (numerical relativity and surrogates, effective-one-body, phenomenological) compared; the
  search pipelines (GstLAL, PyCBC, MBTA, cWB) and the ranking statistic in the glossary; references 115–136.

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

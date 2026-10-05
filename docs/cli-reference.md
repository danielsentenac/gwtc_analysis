# CLI reference

Every option of every mode, generated from `gwtc_analysis/cli.py` by
`python gwtc_analysis/gen_readme_cli_tables.py`. Each mode also has its own help:
`gwtc_analysis <MODE> -h`.

## `catalog_statistics`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--catalogs` | `` | Catalog keys, space-separated (e.g. GWTC-1 GWTC-2.1 GWTC-3 GWTC-4 GWTC-4.1 GWTC-5). ALL takes them all except the updates (GWTC-4.1, update of GWTC-4), which are used only when named. |
| `--out-events` | `catalogs_statistics.tsv` | Output TSV path (per-event table). |
| `--out-report` | `catalogs_statistics.html` | Output HTML report path. |
| `--include-detectors` | `False` | Include detector network via GWOSC v2 calls. |
| `--include-area` | `False` | Compute sky localization area Axx if skymaps are available. |
| `--area-cred` | `0.9` | Credible level for sky area: 0.9→A90, 0.5→A50, 0.95→A95. |
| `--plots-dir` | `cat_plots` | Directory for plots (default: cat_plots). |
| `--data-repo` | `zenodo` | Where to read data from: galaxy \| zenodo \| s3. |
| `--zenodo-version` | `` | With --data-repo zenodo, read an older Zenodo release version of a catalog instead of the latest (e.g. --zenodo-version GWTC-3=v2 GWTC-4=v1). Versions are numbered from the oldest (v1); list them with the zenodo_releases mode. |

## `rates`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--out-rates` | `merger_rates.tsv` | Output TSV of rates per population. |
| `--out-events` | `merger_rates_events.tsv` | Output TSV of the events counted. |
| `--out-report` | `merger_rates.html` | Output HTML report path. |
| `--plots-dir` | `rates_plots` | Directory for plots (default: rates_plots). |
| `--sensitivity-release` | `gwtc5` | LVK search-sensitivity release retrieved automatically from Zenodo: gwtc5 = GWTC-5.0 cumulative, real O3 + O4a + O4b injections (~900 MB); gwtc4 = GWTC-4.0 cumulative, real O3 + O4a injections (~400 MB). |
| `--sensitivity-file` | `` | Local LVK injection HDF file to use instead of --sensitivity-release. |
| `--far-threshold` | `1.0` | FAR threshold [1/yr] for both injections and events. |
| `--ns-max-mass` | `2.5` | Maximum neutron-star mass [Msun] separating NS from BH. |
| `--bbh-kappa` | `2.9` | BBH rate evolution R ∝ (1+z)^kappa. |
| `--bbh-z-ref` | `0.2` | Redshift at which the evolving BBH rate is reported. |
| `--catalogs` | `` | Catalog keys (GWTC-1 GWTC-2.1 GWTC-3 GWTC-4 GWTC-4.1 GWTC-5, or ALL): events and injections are restricted to their observing runs (GWTC-1: O1-O2, GWTC-2.1: O3a, GWTC-3: O3b, GWTC-4: O4a, GWTC-4.1: O4a, GWTC-5: O4b). Default: the runs of the real-injection mixture (O3 onward). |
| `--snr-threshold` | `10.0` | Network SNR threshold for the semi-analytic O1+O2 injections (with GWTC-1). |

## `hubble_constant`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--stages` | `['prepare', 'sample', 'combine', 'reweight', 'report']` | Stages to run (default: all). |
| `--workdir` | `hubble_constant_run` | Work directory (inputs, runs, posterior). |
| `--out-report` | `hubble_constant.html` | Output HTML report path. |
| `--out-summary` | `hubble_constant.tsv` | Output TSV of the posterior quantiles. |
| `--sensitivity-release` | `gwtc4` | LVK search-sensitivity release (and matching catalogs and runs): gwtc4 = GWTC-4.0 cumulative, semi-analytic O1+O2 + real O3+O4a injections; gwtc5 = GWTC-5.0 cumulative, semi-analytic O1+O2 + real O3+O4a+O4b injections. |
| `--catalogs` | `` | Catalog keys (GWTC-1 GWTC-2.1 GWTC-3 GWTC-4 GWTC-4.1 GWTC-5, or ALL): events and injections are restricted to their observing runs. Default: all the runs of --sensitivity-release (gwtc4: O1-O4a; gwtc5: O1-O4b; gwtc4: the published analysis). |
| `--sensitivity-file` | `` | Local LVK injection mixture file (semi-analytic O1+O2 + real) instead of the release's. |
| `--far-threshold` | `0.25` | FAR threshold [1/yr] for the events and the real injections. |
| `--snr-threshold` | `10.0` | Network SNR threshold for the semi-analytic O1+O2 injections. |
| `--min-mass` | `3.0` | Minimum source-frame mass [Msun] of both components (potential neutron stars excluded). |
| `--exclude` | `['GW231123_135430', 'GW200105_162426']` | Events left out. |
| `--pe-cache` | `` | PE cache directory (files/, samples/, index/); default ~/.cache_gwtc_analysis/pe_catalog or $GWTC_PE_CACHE. |
| `--keep-pe-files` | `False` | Keep the full PE files after extraction. |
| `--mass-model` | `plp` | BBH primary-mass model: plp = Power Law + Peak; mltp = Multi Peak. Use one work directory per model. |
| `--seeds` | `[1]` | One sampler run per seed. |
| `--parallel` | `1` | Seeds run at the same time on this machine (each with --npool processes; logs in <workdir>/logs). |
| `--nlive` | `100` | dynesty live points per run. |
| `--npool` | `4` | Worker processes per run: random walks of one seed run at the same time. |
| `--naccept` | `60` | dynesty accepted steps per MCMC walk. |
| `--pe-samples` | `1500` | PE samples per event. |
| `--inj-fraction` | `auto` | Fraction of the found injections used by the sampler runs: 'auto' (a probe chooses the fastest reliable subset, the posterior being then reweighted to all the injections), or a number in (0, 1], 1 = all the injections, as in the paper. |
| `--min-ess-fraction` | `0.5` | With --inj-fraction auto: smallest predicted effective-sample-size fraction accepted for the reweighting to all the injections. |
| `--probe-points` | `30` | With --inj-fraction auto: finite-likelihood prior points used by the probe. |
| `--reweight-pe-samples` | `` | PE samples per event of the reweighting target (default: those of the runs). |
| `--icarogw-python` | `` | Python interpreter of the icarogw environment (default: the current one). |

## `bright_siren`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--src-name` | `GW170817` | Event with an identified host galaxy. |
| `--pe-label` | `` | PE label(s) to use (default: all the labels of the PE file, LowSpin first). |
| `--pe-file` | `` | PE file to read instead of the event's bundle. |
| `--cache-dir` | `.cache_gwosc` | Cache root of the unofficial PE bundle (as in build_unofficial_pe). |
| `--v-recession` | `` | Recession velocity of the host and its uncertainty, km/s (default for GW170817: 3327 72, the NGC 4993 group in the CMB frame). |
| `--v-peculiar` | `` | Peculiar velocity of the host and its uncertainty, km/s (default for GW170817: 310 150). |
| `--redshift` | `` | Hubble-flow redshift of the host and its uncertainty, instead of the velocities (default for GW190521: 0.438 0.0015). |
| `--selection` | `auto` | Selection term: euclidean (GW-limited, nearby sources: beta ∝ H0^3), injections (LVK sensitivity injections of the event's run), auto (euclidean below z = 0.05). |
| `--sensitivity-release` | `` | Injections of the selection term (default: gwtc4). |
| `--sensitivity-file` | `` | Local sensitivity file instead of the release. |
| `--far-threshold` | `0.25` | Found injections: FAR below this, per year. |
| `--snr-threshold` | `10.0` | Found semi-analytic O1+O2 injections: network SNR above this. |
| `--pe-cache` | `` | PE cache of the events read from Zenodo (default: that of hubble_constant). |
| `--sky-radius` | `3.0` | For samples not fixed to the counterpart's position: keep those within this angle (deg). |
| `--spectral-posterior` | `` | Spectral-siren H0 posterior to combine with: a hubble_constant work directory or a posterior TSV with an H0 column. |
| `--h0-range` | `[10.0, 200.0]` | Flat H0 prior range, km/s/Mpc (that of the spectral siren by default). |
| `--out-report` | `bright_siren.html` | Output HTML report path. |
| `--out-summary` | `bright_siren.tsv` | Output TSV of the H0 summary (the posterior grid goes to <name>.posterior.tsv). |
| `--plots-dir` | `bright_siren_plots` | Directory for the plots. |

## `area_law`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--src-name` | `GW250114` | Event (only GW250114). |
| `--cache-dir` | `` | Where the release is extracted (default: the Zenodo cache). |
| `--with-imr` | `False` | Also show the area change of the full-signal PE (NR fits: a consistency check, not a test). |
| `--out-report` | `area_law.html` | Output HTML report path. |
| `--out-summary` | `area_law.tsv` | Output TSV of the comparison with the paper (scans in <name>.truncation.tsv and <name>.ringdown.tsv). |
| `--plots-dir` | `area_law_plots` | Directory for the plots. |

## `stochastic`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--spectral-posterior` | `` | hubble_constant work directory (--mass-model plp) or its posterior TSV. |
| `--rates` | `` | TSV written by the rates mode (--out-rates). |
| `--high-z` | `sfr` | BBH rate beyond the farthest detected events: the star-formation history (sfr), or the fitted shape, which there is the prior's (posterior). |
| `--z-horizon` | `` | Redshift of the farthest detected events (default: from the work directory, else 1). |
| `--n-draws` | `200` | Posterior draws. |
| `--out-report` | `stochastic.html` | Output HTML report path. |
| `--out-summary` | `stochastic.tsv` | Output TSV of Omega_GW(25 Hz) (the spectrum goes to <name>.spectrum.tsv). |
| `--plots-dir` | `stochastic_plots` | Directory for the plots. |

## `neutron_star_eos`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--events` | `['GW170817', 'GW190425']` | Binary neutron stars to combine. |
| `--spin-prior` | `low` | PE analyses with the low-spin (\|chi\| <= 0.05) or high-spin (\|chi\| <= 0.89) prior. |
| `--lambda-max` | `5000.0` | Upper bound of the uniform PE priors on Lambda_1, Lambda_2. |
| `--cache-dir` | `.cache_gwosc` | Cache root of the GW170817 bundle. |
| `--pe-cache` | `` | PE cache of the Zenodo files (default: that of hubble_constant). |
| `--out-report` | `neutron_star_eos.html` | Output HTML report path. |
| `--out-summary` | `neutron_star_eos.tsv` | Output TSV of Lambda_1.4 and R_1.4 (the posteriors go to <name>.posterior.tsv). |
| `--plots-dir` | `neutron_star_eos_plots` | Directory for the plots. |

## `event_selection`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--catalogs` | `` | Catalog keys, space-separated (e.g. GWTC-1 GWTC-2.1 GWTC-3 GWTC-4 GWTC-4.1 GWTC-5). ALL takes them all except the updates (GWTC-4.1, update of GWTC-4), which are used only when named. |
| `--out-selection` | `event_selection.tsv` | Output TSV path for the selected events. |
| `--m1-min` | `` | Minimum primary mass (source frame). |
| `--m1-max` | `` | Maximum primary mass (source frame). |
| `--m2-min` | `` | Minimum secondary mass (source frame). |
| `--m2-max` | `` | Maximum secondary mass (source frame). |
| `--dl-min` | `` | Minimum luminosity distance (Mpc). |
| `--dl-max` | `` | Maximum luminosity distance (Mpc). |
| `--chi-eff-min` | `` | Minimum effective spin chi_eff. |
| `--chi-eff-max` | `` | Maximum effective spin chi_eff. |
| `--preset` | `` | Class of sources (the cuts apply on top): neutron-stars (a component below --ns-max-mass), mass-gap (a component in --mass-gap), hierarchical (primary above --pisn-gap-min, or chi_eff < 0 at 90%%: earlier-generation black holes). |
| `--ns-max-mass` | `3.0` | Maximum neutron-star mass (M_sun). |
| `--mass-gap` | `[3.0, 5.0]` | Lower mass gap between neutron stars and black holes (M_sun). |
| `--pisn-gap-min` | `50.0` | Lower edge of the pair-instability mass gap (M_sun; ~45-65 in the literature). |
| `--out-plot` | `` | Optional PNG of the selected events among all the events of the catalogs (m2 and D_L against m1). |

## `search_skymaps`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--catalogs` | `` | Catalog keys, space-separated (e.g. GWTC-1 GWTC-2.1 GWTC-3 GWTC-4 GWTC-4.1 GWTC-5). ALL takes them all except the updates (GWTC-4.1, update of GWTC-4), which are used only when named. |
| `--ra-deg` | `` | Right ascension (deg). |
| `--dec-deg` | `` | Declination (deg). |
| `--prob` | `0.9` | Credible-level threshold (0–1). Common values: 0.9, 0.5, 0.95. |
| `--skymap-label` | `Mixed` | Label selector used to filter skymap (default: Mixed). |
| `--out-events` | `search_skymaps.tsv` | Output TSV file (default: search_skymaps.tsv). |
| `--out-report` | `search_skymaps.html` | Optional output HTML report path for hits. |
| `--plots-dir` | `sky_plots` | Directory for hit plots (default: sky_plots). |
| `--data-repo` | `zenodo` | Where to read data from: galaxy \| zenodo \| s3. |
| `--zenodo-version` | `` | With --data-repo zenodo, read an older Zenodo release version of a catalog instead of the latest (e.g. --zenodo-version GWTC-3=v2 GWTC-4=v1). Versions are numbered from the oldest (v1); list them with the zenodo_releases mode. |

## `parameters_estimation`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--out-report` | `parameters_estimation.html` | Output HTML report path. |
| `--src-name` | `` | Source event name (e.g. GW231223_032836). |
| `--data-repo` | `zenodo` | Where to read data from: galaxy \| zenodo \| s3. |
| `--zenodo-version` | `` | With --data-repo zenodo, read an older Zenodo release version of a catalog instead of the latest (e.g. --zenodo-version GWTC-3=v2 GWTC-4=v1). Versions are numbered from the oldest (v1); list them with the zenodo_releases mode. |
| `--pe-vars` | `` | Extra posterior sample variables to plot (space-separated). Example: --pe-vars chi_eff chi_p luminosity_distance. |
| `--pe-pairs` | `` | Extra 2D posterior pairs to plot as 'x:y' tokens. Example: --pe-pairs mass_1_source:mass_2_source chi_eff:chi_p. |
| `--plots-dir` | `pe_plots` | Directory for output PE plots (default: pe_plots). |
| `--start` | `0.2` | Default seconds before GPS time for overlay and q-transform windows. |
| `--stop` | `0.1` | Default seconds after GPS time for overlay and q-transform windows. |
| `--fmin` | `20.0` | Default low frequency bound (Hz) used for overlay filtering and q-transform range. |
| `--fmax` | `300.0` | Default high frequency bound (Hz) used for overlay filtering and q-transform range. |
| `--overlay-start` | `` | Override seconds before GPS time for the whitened overlay window. |
| `--overlay-stop` | `` | Override seconds after GPS time for the whitened overlay window. |
| `--overlay-fmin` | `` | Override low frequency bound (Hz) for overlay whitening/bandpass. |
| `--overlay-fmax` | `` | Override high frequency bound (Hz) for overlay whitening/bandpass. |
| `--q-start` | `` | Override seconds before GPS time for the q-transform window. |
| `--q-stop` | `` | Override seconds after GPS time for the q-transform window. |
| `--q-fmin` | `` | Override low frequency bound (Hz) for the q-transform. |
| `--q-fmax` | `` | Override high frequency bound (Hz) for the q-transform. |
| `--q-fscale` | `log` | Frequency axis scaling for q-transform plots (default: log). |
| `--pe-label` | `` | PE label used to select posterior samples and metadata. If omitted and --waveform-engine is provided, the tool selects the closest PE label by substring match in the PE label. If both are omitted: the Mixed label for the posteriors, and for the strain overlay the IMRPhenomXPHM label when the Mixed one has no PSD. |
| `--waveform-engine` | `` | Waveform engine used to generate a time-domain waveform for strain overlay. If omitted, a sensible default engine is used for overlays. |

## `build_unofficial_pe`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--src-name` | `` | Source event name (e.g. GW170817). |
| `--cache-dir` | `.cache_gwosc` | Cache root where unofficial_pe/<bundle>.h5 will be written. |
| `--force` | `False` | Force rebuilding the unofficial bundle even if a cached copy already exists and is up to date. |

## `check_catalogs`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--out-json` | `` | Optional JSON file with the full report. |
| `--sample-events` | `3` | Events of each new list whose PE links are used to find its Zenodo records. |

## `zenodo_releases`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--catalogs` | `['ALL']` | Catalog keys, space-separated (e.g. GWTC-3 GWTC-4). ALL key takes them all. |

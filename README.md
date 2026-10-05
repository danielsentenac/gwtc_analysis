# GWTC Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23077807.svg)](https://doi.org/10.5281/zenodo.23077807)
[![conda-forge](https://img.shields.io/conda/vn/conda-forge/gwtc_analysis.svg)](https://anaconda.org/conda-forge/gwtc_analysis)

## Overview

**GWTC Analysis** is a command-line analysis suite for exploring publicly released
**Gravitational-Wave Transient Catalogs (GWTC)** from the **LIGO–Virgo–KAGRA (LVK) Collaboration**.

The tool provides:

- Search of gravitational-wave sky localizations around a given sky position
- Visualization of parameter-estimation results for individual events
- Selection of events based on physical constraints (masses, distance, χ_eff), or by class of sources: neutron stars, the lower mass gap (3–5 M☉), hierarchical-merger candidates (pair-instability gap, negative χ_eff)
- Global catalog statistics, including detector-network participation, sky-localization performance, and remnants (radiated energy, final-spin estimate)
- BNS, NSBH and BBH merger-rate estimates (R = N / ⟨VT⟩) from the catalogs and the LVK search-sensitivity injections
- Hubble-constant estimate from the BBH mass spectrum (spectral siren, with icarogw)
- BBH effective-spin population and its correlation with the mass ratio
- Neutron-star equation of state from GW170817 and GW190425 jointly (Λ₁.₄, R₁.₄)
- Predicted background of unresolved compact binaries, Ω_GW(f), against the stochastic upper limits
- Test of Hawking's area law with GW250114 (inspiral against ringdown, reproducing the published 4.4σ)
- Hubble-constant estimate from events with an identified host (bright siren: GW170817, the candidate GW190521 flare), alone or combined with the spectral siren

All gravitational-wave data products are retrieved from the **Gravitational Wave Open Science Center (GWOSC)**,
or from supported alternative repositories (Zenodo / S3 / Galaxy collections).

📖 **Documentation:** https://danielsentenac.github.io/gwtc_analysis/ (user guide, methods of the rates and Hubble-constant modes, references).

<!-- CATALOG_COVERAGE_BEGIN -->
> **Catalogs in version 0.7.0**: GWTC-1, GWTC-2.1, GWTC-3, GWTC-4.0, GWTC-5.0, and the update GWTC-4.1 (of GWTC-4.0), observing runs O1 to O4b. **Latest catalog: GWTC-5.0 (O4b)**. Registry checked against GWOSC and Zenodo on 2026-10-05; catalogs published later need a newer version of the package (`gwtc_analysis check_catalogs` tells whether GWOSC has published one).
<!-- CATALOG_COVERAGE_END -->

---

## Containerized Distribution (Docker)

`gwtc_analysis` is distributed as a ready-to-use **Docker container** named **`gwtc-tool`**.

- **Docker image name:** `gwtc-tool`
- **Docker Hub repository:** https://hub.docker.com/r/danielsentenac/gwtc-tool/

Using the Docker image is recommended for reproducibility, portability, and integration with workflow systems
(e.g. CI pipelines, computing clusters).

---

## Astrophysical Sources

The GWTC catalogs contain compact binary merger events involving:

- **Binary Black Holes (BBH)**
- **Binary Neutron Stars (BNS)**
- **Neutron Star – Black Hole systems (NSBH)**

These mergers are detected by the **LVK detector network**:
**H1 (Hanford), L1 (Livingston), V1 (Virgo), K1 (KAGRA)**.

---

## Supported GW Catalog Names

Catalog identifiers are **case-sensitive**. All of them are **confident catalogs** (see below the table):

| Key | Catalog | Observing runs of its events | Events | PE and skymaps on Zenodo |
|---|---|---|---|---|
| `GWTC-1` | GWTC-1 | O1, O2 | 11 | the GWTC-2.1 release, [Zenodo 6513631](https://zenodo.org/records/6513631), which re-analysed O1–O2 |
| `GWTC-2.1` | GWTC-2.1 | O3a, plus 10 O1–O2 events re-analysed | 54 | [Zenodo 6513631](https://zenodo.org/records/6513631) |
| `GWTC-3` | GWTC-3 | O3b | 35 | [Zenodo 22685054](https://zenodo.org/records/22685054) |
| `GWTC-4` | GWTC-4.0 | O4a, plus GW230518 from the engineering run ER15 | 129 | [Zenodo 17602505](https://zenodo.org/records/17602505) |
| `GWTC-4.1` | GWTC-4.1, update of GWTC-4.0 | O4a, plus two events from ER15 | 140 | [Zenodo 20275769](https://zenodo.org/records/20275769) |
| `GWTC-5` | GWTC-5.0 | O4b, plus 5 events of 6–8 April 2024, just before O4b | 161 | [Zenodo 20348005](https://zenodo.org/records/20348005) (part 1, with the skymaps) and [20348006](https://zenodo.org/records/20348006) (part 2) |
| `ALL` | all the catalogs above except the update GWTC-4.1 | O1 to O4b | | |

`GWTC-4.1` is an **update** of GWTC-4.0: the same O4a data re-analysed, with the 129 events of
GWTC-4.0 and 11 new ones. It is used only when named (`--catalogs GWTC-4.1`), in place of `GWTC-4` for
the O4a events: `ALL` and the defaults of every mode keep GWTC-4.0, the catalog of the published
analyses. Its PE files are read only on request (`--zenodo-version GWTC-4.1=latest`).

All the keys are confident catalogs: every event of their GWOSC lists has p_astro ≥ 0.5 (the
re-analysed O1–O2 events of GWTC-2.1 carry no p_astro value). Event counts of the GWOSC lists in
September 2026; the Zenodo version read is resolved at run time (see below).

`GWTC-5` resolves to the GWOSC `GWTC-5.0` endpoint and to the Zenodo records below.

All the catalogs, runs and injection releases are described in `gwtc_analysis/catalog_registry.py`.
`gwtc_analysis check_catalogs` compares it with what GWOSC and Zenodo publish, and drafts the registry
entry of any new catalog (see [Adding a new catalog](https://danielsentenac.github.io/gwtc_analysis/data-sources/#adding-a-new-catalog)).

---

## Command-Line Interface (CLI)

```text
usage: gwtc_analysis [-h] MODE ...

positional arguments:
  MODE
    catalog_statistics
    rates
    hubble_constant
    bright_siren
    area_law
    stochastic
    neutron_star_eos
    spin_population
    event_selection
    search_skymaps
    parameters_estimation
    build_unofficial_pe
    check_catalogs
    zenodo_releases
```

Each mode has its own help:

```bash
gwtc_analysis <MODE> -h
```

The `gwtc_analysis` command is installed with the package (conda-forge, PyPI or Docker); from a source checkout that is not installed, `python -m gwtc_analysis.cli` is equivalent.

---

## General Units and Ranges

- Right Ascension: degrees [0, 360)
- Declination: degrees [-90, +90]
- Probability threshold: [0, 1]
- Masses: solar masses (M☉)
- Distances: megaparsecs (Mpc)

---

## Data repositories

The GWTC catalogs (Parameter Estimation and Skymaps) can be directly downloaded from different supports:

- The Zenodo portal (official catalogs PE/skymaps tarballs): GWTC-2.1 (https://zenodo.org/records/6513631), GWTC-3 (https://zenodo.org/records/22685054), GWTC-4.0 (https://zenodo.org/records/17602505) and GWTC-5.0, which is split across two records: https://zenodo.org/records/20348005 (part 1, plus the archived skymaps tarball) and https://zenodo.org/records/20348006 (part 2). These records are only the starting point: the Zenodo version used is resolved at run time (see [Zenodo release versions](#zenodo-release-versions)).
- A s3 Minio bucket called gwtc on  https://minio-dev.odahub.fr
- Galaxy collections under the name GWTC at https://usegalaxy.org. With `--data-repo galaxy`, parameter-estimation files are first looked up in locally staged Galaxy collections (`./galaxy_inputs/<CATALOG>-PE`) and, if none is found, downloaded over HTTP from the public usegalaxy.org **"GWTC" published history** (anonymous, no API key — works while the history stays published).

### Zenodo release versions

The Zenodo releases are versioned (for instance GWTC-3 has v1, v2 and v3). With `--data-repo zenodo`, each catalog uses its **latest** version by default; the version listings are fetched from the Zenodo API and cached for one day in `~/.cache_gwtc_analysis/zenodo`, so a new release is picked up automatically.

To read an older version, pass `--zenodo-version CATALOG=VERSION` (modes `catalog_statistics`, `search_skymaps`, `parameters_estimation`). Versions are numbered from the oldest (`v1`); `latest` is also accepted:

```bash
gwtc_analysis zenodo_releases --catalogs GWTC-3 GWTC-4   # list the versions
gwtc_analysis search_skymaps --catalogs GWTC-3 --ra-deg 40 --dec-deg -30 --zenodo-version GWTC-3=v2
gwtc_analysis parameters_estimation --src-name GW200105_162426 --zenodo-version GWTC-3=v2
```

Skymap tarballs are cached per Zenodo record (`.cache_gwosc/zenodo_<record>_<file>`), and the PE index is rebuilt when the selected records change. If zenodo.org is unreachable, the cached version listing is used; with no cache at all, the latest version falls back to the record listed above.

---

## Usage

The tool runs as a Python package installed from conda-forge or PyPI, or as a Docker tool:

- [conda package](https://anaconda.org/conda-forge/gwtc_analysis) on conda-forge: `conda install -c conda-forge gwtc_analysis`
- [PyPI package](https://pypi.org/project/gwtc_analysis/): `pip install gwtc_analysis` (PyPI displays it as `gwtc-analysis`, the same project; light dependencies only: the PE, strain and skymap modes also need the GW software stack, e.g. an IGWN conda environment)
- [docker image](https://hub.docker.com/r/danielsentenac/gwtc-tool/): `docker pull danielsentenac/gwtc-tool`

Galaxy, like the S3 bucket, is only a data repository (`--data-repo galaxy`, below).

### Inputs
- Catalog selections are passed as parameters separated by space
- Data repositories accept `--data-repo` to choose where data products are read from:
	- `galaxy`: read inputs from locally staged Galaxy collections, falling back to the public usegalaxy.org "GWTC" published history over HTTP when no staged file is found
	- `zenodo`: official releases from Zenodo
	- `s3`: S3-compatible bucket

### Outputs
- TSV tables
- HTML reports
- Plot images

---

## CLI options (auto-generated)

The tables below are generated directly from `cli.py` to stay aligned with the real CLI.

To regenerate locally (from the repository root):

```bash
python gwtc_analysis/gen_readme_cli_tables.py
```

<!-- CLI_TABLES_BEGIN -->

### `catalog_statistics`

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

### `rates`

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

### `hubble_constant`

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

### `bright_siren`

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

### `area_law`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--src-name` | `GW250114` | Event (only GW250114). |
| `--cache-dir` | `` | Where the release is extracted (default: the Zenodo cache). |
| `--with-imr` | `False` | Also show the area change of the full-signal PE (NR fits: a consistency check, not a test). |
| `--out-report` | `area_law.html` | Output HTML report path. |
| `--out-summary` | `area_law.tsv` | Output TSV of the comparison with the paper (scans in <name>.truncation.tsv and <name>.ringdown.tsv). |
| `--plots-dir` | `area_law_plots` | Directory for the plots. |

### `stochastic`

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

### `neutron_star_eos`

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

### `spin_population`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--no-correlation` | `False` | Gaussian chi_eff model without the q slope. |
| `--sensitivity-release` | `` | Injections and catalogs (default: gwtc4). |
| `--sensitivity-file` | `` | Local sensitivity file instead of the release. |
| `--far-threshold` | `0.25` | Events and found injections: FAR below this. |
| `--snr-threshold` | `10.0` | Found semi-analytic O1+O2 injections. |
| `--min-mass` | `3.0` | Both source-frame masses above it. |
| `--exclude` | `['GW231123_135430', 'GW200105_162426']` | Events left out. |
| `--pe-cache` | `` | PE cache (default: that of hubble_constant). |
| `--pe-samples` | `5000` | PE samples per event. |
| `--max-injections` | `250000` | Random subset of the found injections (0: all). |
| `--walkers` | `20` | emcee walkers. |
| `--steps` | `1200` | emcee steps (the first third is burn-in). |
| `--out-report` | `spin_population.html` | Output HTML report path. |
| `--out-summary` | `spin_population.tsv` | Output TSV of the posterior quantiles (the samples go to <name>.posterior.tsv). |
| `--plots-dir` | `spin_population_plots` | Directory for the plots. |

### `event_selection`

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

### `search_skymaps`

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

### `parameters_estimation`

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
| `--no-skymap-3d` | `True` | Skip the 3D sky map (FITS map of the Zenodo skymap archive: credible volume, distance by direction, host-galaxy candidates). |
| `--galaxies` | `glade` | Galaxy catalog cross-matched with the 3D sky map: 'glade' (GLADE+ via VizieR, default), a CSV/TSV file with ra, dec and dist [Mpc] or z columns, or 'none'. |
| `--galaxy-max-area` | `100.0` | Largest 90 percent credible area (deg2) for which GLADE+ is queried (default: 100). |

### `build_unofficial_pe`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--src-name` | `` | Source event name (e.g. GW170817). |
| `--cache-dir` | `.cache_gwosc` | Cache root where unofficial_pe/<bundle>.h5 will be written. |
| `--force` | `False` | Force rebuilding the unofficial bundle even if a cached copy already exists and is up to date. |

### `check_catalogs`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--out-json` | `` | Optional JSON file with the full report. |
| `--sample-events` | `3` | Events of each new list whose PE links are used to find its Zenodo records. |

### `zenodo_releases`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--catalogs` | `['ALL']` | Catalog keys, space-separated (e.g. GWTC-3 GWTC-4). ALL key takes them all. |

<!-- CLI_TABLES_END -->

### `rates`: Merger Rates From The Catalogs

`rates` estimates the BNS, NSBH and BBH merger rates (per Gpc³ per year) as R = N / ⟨VT⟩:

- **N**: GWOSC candidates (confident and marginal lists of GWTC-2.1, GWTC-3, GWTC-4.0, GWTC-5.0) inside the observing periods covered by the injections, with FAR below `--far-threshold` (default 1/yr), classified by their median source-frame masses (neutron stars below `--ns-max-mass`, default 2.5 M☉).
- **⟨VT⟩**: the sensitive volume-time, from the LVK search-sensitivity injections (simulated signals added to the real data and searched by the real pipelines), reweighted to each population by importance sampling. The injection file is retrieved automatically from Zenodo (latest version of the record) and cached in `~/.cache_gwtc_analysis/zenodo`. Choose the release with `--sensitivity-release`:
  - `gwtc5` (default): GWTC-5.0 cumulative, real O3 + O4a + O4b injections ([Zenodo 19500052](https://zenodo.org/records/19500052), ~900 MB);
  - `gwtc4`: GWTC-4.0 cumulative, real O3 + O4a injections ([Zenodo 16740128](https://zenodo.org/records/16740128), ~400 MB).

  `--sensitivity-file` uses a local injection file instead (same LVK mixture format).
- **Populations** (fixed shapes): BNS with both masses uniform in [1, 2.5] M☉; NSBH with the black hole ∝ m^-2.35 on [2.5, 40] M☉; BBH with the GWTC-3 *Power Law + Peak* model, reported with R ∝ (1+z)^κ at z = 0.2 (`--bbh-kappa`, `--bbh-z-ref`) and without evolution.
- **Intervals**: 90% Poisson (Jeffreys prior). The LVK population papers fit the population shapes together with the rates, so their intervals are wider and model-dependent.

- **Catalogs** (`--catalogs`, as in the other modes): the events **and** the injections are restricted to the observing runs of the selected catalogs (GWTC-1: O1–O2, GWTC-2.1: O3a, GWTC-3: O3b, GWTC-4: O4a, GWTC-5: O4b; `ALL` for O1 to O4b), so that counts and ⟨VT⟩ describe the same observing time. GWTC-1 uses the release's mixture with semi-analytic O1+O2 injections (`--snr-threshold`, 10). For O4b alone (`--catalogs GWTC-5`): 104 BBH candidates, BBH rate 24.8 [21.0, 29.0] Gpc⁻³ yr⁻¹ at z = 0.2.

Outputs: `--out-rates` (TSV per population), `--out-events` (TSV of the events counted), `--out-report` (HTML report with the observed and the selection-corrected primary-mass distributions).

```bash
gwtc_analysis rates                               # GWTC-5.0 injections (O3 + O4a + O4b)
gwtc_analysis rates --sensitivity-release gwtc4   # GWTC-4.0 injections (O3 + O4a)
```

| Release | Candidates | BNS | NSBH | BBH at z = 0.2 |
|---|---|---|---|---|
| `gwtc5` (O3–O4b, 2.5 yr) | 259 | 26 [4, 87] | 33 [13, 67] | 25 [23, 28] |
| `gwtc4` (O3–O4a, 1.7 yr) | 155 | 43 [6, 143] | 54 [21, 109] | 26 [23, 30] |

Rates in Gpc⁻³ yr⁻¹, median [90%]. They are consistent with the LVK population papers: [GWTC-5.0](https://arxiv.org/abs/2605.27226) (BBH 27.5–49.4 at z = 0.2 for masses 2.5–200 M☉), [GWTC-4.0](https://arxiv.org/abs/2508.18083) (z = 0: BNS 7.6–250, NSBH 9.1–84, BBH 14–26) and [GWTC-3](https://arxiv.org/abs/2111.03634).

### `hubble_constant`: Hubble Constant From The BBH Mass Spectrum

`hubble_constant` measures H₀ with the **spectral-siren** method. Each event gives its luminosity distance D_L and its detector-frame masses m_det = m_src (1+z), but not its redshift. The redshift comes from the population: for a trial H₀, every D_L gives a z and every m_det a source-frame mass, and only the right H₀ makes the events near and far fall on one distance-independent mass distribution. The *Power Law + Peak* (PLP) mass model and the Madau–Dickinson rate evolution are therefore fitted together with H₀ (flat ΛCDM, Ω_m = 0.3065), with the hierarchical likelihood of [icarogw](https://github.com/icarogw-developers/icarogw) and bilby/dynesty. Selection effects are corrected with the LVK search-sensitivity injections.

The default setup reproduces the spectral-siren measurements of the [GWTC-4.0 cosmology paper](https://arxiv.org/abs/2509.04348) (published version v3): H₀ = 105.5 (+46.4 / −35.8) km/s/Mpc with the *Power Law + Peak* mass model (`--mass-model plp`, the default) and 72.3 (+42.5 / −25.6) km/s/Mpc with the *Multi Peak* model (`--mass-model mltp`: a power law and two Gaussian peaks, found near 9 and 27 M☉ in the paper and near 10 and 30 M☉ in the reproduction run). Use one work directory per mass model. The setup:

- **Catalogs** (`--catalogs`): events and injections restricted to the observing runs of the selected catalogs; by default all the runs of `--sensitivity-release` (`gwtc4`: O1–O4a, the published analysis; `gwtc5`: O1–O4b, compared with the GWTC-5.0 cosmology paper, MLTP 71.0 (+21.0 / −17.5) km/s/Mpc).
- **Events**: BBHs of O1–O4a from the GWOSC confident and marginal lists, lowest FAR at most 0.25/yr (`--far-threshold`; the published FARs are rounded, so they are compared inclusively), both source-frame masses above 3 M☉ (`--min-mass`), GW231123 and GW200105 left out (`--exclude`): 137 events, as in the paper.
- **PE samples**: from the Zenodo PE releases, `C01:IMRPhenomXPHM` up to O3 and `C00:IMRPhenomXPHM-SpinTaylor` in O4a, reduced to (m1_det, m2_det, D_L). The PE distance prior of each event is read from its file (D_L² up to O3, uniform in source-frame comoving volume in O4a) and divided out.
- **Injections** (`--sensitivity-release`): `gwtc4` (default), the GWTC-4.0 semi-analytic O1+O2 + real O3+O4a mixture ([Zenodo 16740128](https://zenodo.org/records/16740128)), found when the semi-analytic SNR exceeds 10 (`--snr-threshold`) or the lowest search FAR is below the threshold; `gwtc5` adds GWTC-5.0 and O4b ([Zenodo 19500052](https://zenodo.org/records/19500052)) and is not yet validated against a published result.
- **Priors**: those of the paper (Tables 3 and 6), H₀ uniform in [10, 200] km/s/Mpc.

The work is split into stages (`--stages`, all by default), sharing `--workdir`:

| Stage | Does | Cost |
|---|---|---|
| `prepare` | selects the events, downloads their PE files (restartable; only the extracted samples are kept in `--pe-cache` unless `--keep-pe-files`), prepares the injections → `inputs.h5`, `events.tsv` | ~35 GB of downloads the first time |
| `sample` | with `--inj-fraction auto` (default), first a probe (a few minutes) that chooses the fastest reliable subset of the injections; then one dynesty run per `--seeds` value (`--mass-model`, `--nlive`, `--npool`, `--naccept`, `--pe-samples`), `--parallel` of them at a time on this machine (logs in `<workdir>/logs`); resumable from its checkpoint | hours per run |
| `combine` | merges the runs → `posterior.tsv`, `corner.png`, `summary.json`, with the effective numbers of injections and PE samples over the posterior (icarogw's stability criteria) | minutes |
| `reweight` | when the runs used a subset of the injections, reweights their posterior to all of them (importance weights exp(ln L_all − ln L_runs)) → `posterior_reweighted.tsv`, with the effective sample size in `summary.json` | minutes to an hour |
| `report` | `--out-report` (HTML) and `--out-summary` (TSV of the posterior quantiles) | seconds |

**icarogw.** `gwtc_analysis/h0_icarogw.py` is a driver of icarogw, not a modified copy: icarogw is used as installed, through its public API. icarogw provides the hierarchical likelihood (PE and injection reweighting, selection term, scale-free rate marginalisation, effective-sample-size checks), the population models (`massprior_PowerLawPeak` with the `m1m2_conditioned_lowpass` smoothing, `rateevolution_Madau`, `FlatLambdaCDM_wrap`, combined by `CBC_vanilla_rate`) and the detector-frame conversion for each trial H₀. The driver reads `inputs.h5` into icarogw's `posterior_samples` and `injections` objects, chooses the model components and the priors (Tables 3 and 6 of the paper), runs bilby/dynesty, merges the runs and computes the diagnostics with icarogw's own methods. The analysis choices made here, outside icarogw, are the input preparation in `hubble_constant.py` (event selection, PE distance prior read from each file, injection draw density carried to the detector frame with the spin part divided out and the mixture weights applied) and three settings of the driver: at least 10 effective PE samples per event (the paper's choice; the default of icarogw's likelihood class is 20), at least 4 × N_events effective injections (icarogw's default), and the `--inj-fraction` subset of the injections.

Only the `sample` and `combine` stages need icarogw; `prepare` and `report` run in the gwtc_analysis environment. icarogw needs Python ≥ 3.12, so it usually has its own environment; pass its interpreter with `--icarogw-python` (default: the interpreter running gwtc_analysis; the mode stops with an error before sampling if icarogw or bilby cannot be imported there). The stages run `h0_icarogw.py` with it, in CPU mode (a `config.py` with `CUPY=False` in the work directory) and with the environment's `lib/` on `LD_LIBRARY_PATH` (for its `libstdc++`). icarogw is not on PyPI; an installation that works:

```bash
conda create -n icarogw python=3.12
conda activate icarogw
export TMPDIR=~/tmp                                                  # the torch wheels are large
pip install torch --index-url https://download.pytorch.org/whl/cpu   # CPU torch first, not the CUDA build
pip install git+https://github.com/icarogw-developers/icarogw.git
```

If other packages in that environment need an older numpy (e.g. ligo.skymap), pin it (`numpy==2.1.1 scipy==1.14.1` worked).

```bash
# prepare in the gwtc_analysis environment, then 4 runs, 2 at a time with 2 processes each, and the report
gwtc_analysis hubble_constant --stages prepare
gwtc_analysis hubble_constant --stages sample combine report \
    --icarogw-python ~/.conda/envs/icarogw/bin/python --seeds 1 2 3 4 --parallel 2 --npool 2
```

**Seeds.** Each seed is an independent dynesty run (`result/<model>_seed<N>_result.json`, with `<model>` = `plp` or `mltp`); `combine` merges all the finished ones, weighted by their evidence. All runs sample the same likelihood: the PE samples are shuffled once in `prepare` and the injection subset is drawn with a fixed seed, and `run_settings.json` refuses runs with other `--mass-model`, `--nlive`, `--pe-samples` or `--inj-fraction` values in the same work directory. Launching again resumes the interrupted runs from their checkpoint and skips the finished ones. A lock file (`result/<model>_seed<N>.lock`) prevents the same seed from running twice at once; interrupting the launcher (Ctrl-C) stops its runs after they write their checkpoint.

**`--npool` and `--parallel`.** Nearly all the time of a run goes into likelihood evaluations. At each iteration dynesty replaces the live point of lowest likelihood L_min by a new point with L > L_min, found by a random walk (about `--naccept` accepted steps, one likelihood evaluation per step) from another live point. With `--npool N`, bilby starts N worker processes, each holding a copy of the likelihood, and dynesty runs N such walks at the same time; their new points replace the next N worst points. The walks all start from the same L_min, so some of their points are no longer good enough when used: N workers give less than N times the speed. `--npool` makes one seed finish sooner without changing its result; `--parallel` runs several seeds at the same time. The machine then runs `--parallel` × `--npool` processes, which should not exceed its number of CPUs, and each of them holds the likelihood data in memory. For the same CPUs, several seeds with few workers each use the machine better than one seed with many workers.

| Machine | Suggested settings |
|---|---|
| 4 CPUs, 8 GB (laptop) | `--parallel 1 --npool 4`, or `--parallel 2 --npool 2` if memory allows |
| 8 CPUs, 16 GB | `--parallel 2 --npool 4` |

**How many seeds?** The seeds do not change the physics: they set how precisely the sampler describes the posterior. They serve two purposes:

1. **Checking that the runs agree** (at least 2 seeds). The evidences ln Z of the runs should agree within their quoted errors (about 0.4), and so should their H₀ intervals. Runs that disagree beyond their errors are not fixed by more seeds but by more live points (`--nlive`).
2. **Precision of the quoted numbers.** One run of 100 live points gives about 560 posterior samples, so its median wanders. In the 10 runs of the reproduction above, the per-run H₀ medians range from 111.7 to 126.1 km/s/Mpc (standard deviation 4.8), and the ln Z values have a standard deviation of 0.32, consistent with their errors. Combining N runs divides the scatter by about √N:

| Seeds (100 live points) | Uncertainty on the H₀ median | Relative to the posterior width (±40) |
|---|---|---|
| 1 | ±4.8 km/s/Mpc | 12% |
| 4 | ±2.4 km/s/Mpc | 6% |
| 10 | ±1.5 km/s/Mpc | 4% |

The PLP posterior is broad, so 3–5 seeds give the result to two significant digits; 10 seeds allow a comparison with a published value at the level of a few km/s/Mpc. The error is a fixed fraction of the posterior width, so the same numbers of seeds hold for narrower posteriors. Fewer runs with more live points are equivalent: 10 runs of 100 live points give about as many samples as one run of about 1000, but small runs can be spread over machines and interrupted, and each explores less carefully, which makes the agreement check more important.

| Purpose | Settings |
|---|---|
| Quick look | 2 seeds |
| Result to report | 4–5 seeds, or 2 seeds with `--nlive 500` |
| Precise comparison with a paper | about 10 seeds |

Runs can be spread over several machines that share the work directory: start `python gwtc_analysis/h0_icarogw.py run --workdir DIR --seed N` with the icarogw interpreter on each, then run the `combine` and `report` stages once.

**Without icarogw on the local machine.** `h0_icarogw.py` only needs numpy, h5py, icarogw and bilby, so the sampling can run on another machine that has icarogw (e.g. a computing cluster):

1. locally: `hubble_constant --stages prepare --workdir DIR`, then copy `DIR/inputs.h5` (about 45 MB) and `gwtc_analysis/h0_icarogw.py` to a work directory on the remote machine;
2. remotely, with the icarogw interpreter (and `LD_LIBRARY_PATH=<env>/lib` if needed): `python h0_icarogw.py run --workdir RDIR --seed N` for each seed, then `python h0_icarogw.py combine --workdir RDIR`;
3. locally: copy `RDIR/summary.json`, `RDIR/posterior.tsv` and `RDIR/corner.png` (a few MB) back into `DIR`, which still holds `events.tsv`, and run `hubble_constant --stages report --workdir DIR`.

With 10 runs of 100 live points and 1500 PE samples per event, sampled with 10% of the found injections and 136 events, then reweighted to all the injections and the paper's 137 events, the result is H₀ = 105.8 (+44.7 / −33.2) km/s/Mpc [68%], 90%: 54.4–175.7, and μ_g = 28.7 (+3.8 / −4.6) M☉ (PLP), against 105.5 (+46.4 / −35.8), 90%: 50.5–176.1, and 28.3 (+4.1 / −4.4) M☉ in the published paper; MLTP gives 78.6 (+38.0 / −26.5) against 72.3 (+42.5 / −25.6). Without the reweighting, the 10% subset alone shifts H₀ by about +13 km/s/Mpc for PLP (119.3) and +10 for MLTP (89.1): a subset estimates the detectable fraction without bias, but its Monte Carlo noise, raised to the power N = 136 in the likelihood, tilts the posterior. This is why `--inj-fraction auto` (the default) probes the subsets before sampling and reweights afterwards; `--inj-fraction 1` samples with all the injections, as the paper does. The number of PE samples per event also matters (3000 instead of 1500 raises H₀ by about 5 km/s/Mpc). H₀ is strongly anti-correlated with μ_g, and the upper part of its interval depends on the prior bound of 200 km/s/Mpc. The lower peak of MLTP, near 9–10 M☉ (σ ≈ 0.8 M☉ in the run), is a much sharper feature than the ~28–30 M☉ bump, hence its tighter and lower H₀; FullPop-4.0 (72.9 km/s/Mpc in the paper) is not implemented.

### `bright_siren`: Hubble Constant From An Identified Host

`bright_siren` measures H₀ with the **bright-siren** method: the luminosity distance of an event from the GW signal alone, at the sky position of its counterpart, against the Hubble-flow redshift of its host, with d_L = (c/H₀) D(z) in flat ΛCDM (Ω_m = 0.3065). The likelihood marginalizes over the true redshift with a population uniform in comoving volume, divides out the PE distance prior recorded in the file, and corrects for selection effects. It runs in seconds.

```bash
gwtc_analysis bright_siren                                            # GW170817 and NGC 4993
gwtc_analysis bright_siren --spectral-posterior hubble_constant_run   # combined with the spectral siren
gwtc_analysis bright_siren --src-name GW190521                        # candidate AGN flare ZTF19abanrhr
```

| Event | Host | Selection term | H₀ (km/s/Mpc): maximum a posteriori, 68% |
|---|---|---|---|
| GW170817, LowSpin | NGC 4993: v_H = 3327 − 310 = 3017 ± 166 km/s ([LVK 2017](https://arxiv.org/abs/1710.05835)) | euclidean, β ∝ H₀³ | 69.8 (+23.6 / −8.2) |
| GW170817, HighSpin | | | 69.4 (+14.8 / −7.4) |
| GW170817, published (2017 distance) | | | 70.0 (+12.0 / −8.0) |
| GW170817 × spectral siren (Power Law + Peak) | | | 71.7 (+22.3 / −8.0) |
| GW190521, IMRPhenomXPHM | **candidate** AGN flare at z = 0.438 ([Graham et al. 2020](https://arxiv.org/abs/2006.14122)) | O3a injections, BBH population | 30.4 (+61.0 / −7.0), bimodal |

- **Selection term** (`--selection`, `auto` by default):
  - **nearby events** (z < 0.05): GW-limited detection, β ∝ H₀³. This cancels the volume factor and gives the LVK 2017 formula p(H₀) ∝ Σᵢ N(v_H; H₀ dᵢ, σ).
  - **distant events:** β is computed from the LVK sensitivity injections of the event's run, reweighted at each trial H₀ to the population of its class.
- **GW170817:** the upper side of its interval is wider than published because the public GWTC-1 distance has a longer low-distance tail, from the distance–inclination degeneracy. Its samples come from the [unofficial bundle](#build_unofficial_pe-unofficial-bundle-workflow), whose sky position is fixed to AT2017gfo.
- **GW190521:** the association with ZTF19abanrhr is not established; [Ashton et al. 2021](https://arxiv.org/abs/2009.12346) find odds of 1 to 12. The report gives H₀ *if* the flare is the counterpart.
- **Validation:** mock bright sirens (simulated detections at z ≈ 0.5 with a known H₀) check the analysis. It recovers H₀ = 70 to within its uncertainty with the selection term, against a 3σ bias without it, and passes a P–P coverage test.
- **Combination:** `--spectral-posterior` takes a `hubble_constant` work directory or a TSV with an `H0` column; the posteriors multiply (different events, same flat prior 10–200 km/s/Mpc).

Method, validation and options: [bright_siren](https://danielsentenac.github.io/gwtc_analysis/modes/bright-siren/) in the documentation.

### `area_law`: Hawking's Area Law With GW250114

`area_law` tests Hawking's area law, A_f ≥ A₁ + A₂, with GW250114 (network SNR 80), from the LVK data release of its discovery paper ([Zenodo 16877102](https://zenodo.org/records/16877102), [arXiv:2509.08054](https://arxiv.org/abs/2509.08054)). The Kerr horizon area is A = 8π (GM/c²)² (1 + √(1 − χ²)).

- **The two sides come from different parts of the signal:** the initial areas from parameter estimation on the data truncated before the peak, and the remnant area from fits of the ringdown quasinormal modes, which use the Kerr spectrum alone. The full-signal PE, whose remnant comes from fits to numerical relativity, obeys the law by construction; `--with-imr` shows it only for comparison.
- **Result:** A_f > A_i at 4.45σ (published 4.4σ) for a truncation 40 M before the peak; at least 3.36σ over all the truncations (3.4σ); above 5σ from −10 M (−10 M); 3.63σ with the overtone from 6 M_f (3.6σ). (A_f − A_i)/A_i = 0.73 (90%: 0.41–1.16).
- **Only GW250114** has a public release of inspiral-only and ringdown-only analyses.

```bash
gwtc_analysis area_law               # downloads 114 MB once
```

### `stochastic`: Predicted Compact-Binary Background

`stochastic` predicts the gravitational-wave background of the compact binaries too faint to be detected one by one, Ω_GW(f) = f/(ρ_c c²) ∫ dz R(z) ⟨dE/df_s⟩(f(1+z)) / [(1+z) H(z)], and compares it with the upper limits of the stochastic searches.

- **BBH:** the Power Law + Peak masses and the Madau–Dickinson rate shape of each draw of a spectral-siren posterior (`hubble_constant`), normalized to the BBH rate at z = 0.2 of the `rates` mode, with the inspiral–merger–ringdown spectrum of Ajith et al. 2008. Beyond the farthest detected events (z ≈ 1) the rate follows the star-formation history (`--high-z sfr`, default) or the fitted shape, there the prior's (`--high-z posterior`).
- **BNS, NSBH:** the local rates of `rates`, the star-formation history, the inspiral spectrum up to the ISCO.
- **Result (GWTC-4.0):** Ω_GW(25 Hz) = 7.5 × 10⁻¹⁰ (90%: 5.5–11.7 × 10⁻¹⁰), BBH 6.0 × 10⁻¹⁰. The LVK prediction from GWTC-5.0 is 6.3 (+5.0 / −2.2) × 10⁻¹⁰, and the upper limit from the data through April 2025 is 2.0 × 10⁻⁹ ([arXiv:2608.23477](https://arxiv.org/abs/2608.23477)): a factor 2.7 below detection.

```bash
gwtc_analysis rates --out-rates merger_rates.tsv
gwtc_analysis stochastic --spectral-posterior hubble_constant_run --rates merger_rates.tsv
```

The `hubble_constant` report now also plots the fitted merger-rate evolution R(z)/R(0) against the star-formation history, and the `parameters_estimation` report summarizes the remnant (final mass and spin, radiated energy, peak luminosity) of every label.

### `neutron_star_eos`: Joint Neutron-Star Equation Of State

`neutron_star_eos` constrains the equation of state of neutron-star matter with GW170817 and GW190425 together. All neutron stars share one equation of state, so the tidal deformabilities of the four stars follow one curve, Λ(m) = Λ₁.₄ (m / 1.4 M☉)⁻⁶ ([De et al. 2018](https://arxiv.org/abs/1804.08583)). The likelihood of Λ₁.₄ is built from each event's PE samples of the measured combination Λ̃, with the PE priors on Λ divided out exactly; the events multiply. The radius follows from Λ₁.₄ = 2.88 × 10⁻⁶ (R₁.₄/km)^7.5 ([Annala et al. 2018](https://arxiv.org/abs/1711.02644)).

| Analysis (low-spin priors) | Λ₁.₄, median (90%) | R₁.₄ (km) |
|---|---|---|
| GW170817 | 216 (80–628) | 11.2 (9.8–12.9) |
| GW190425 | 628 (109–1783) | 12.9 (10.2–14.9) |
| **joint** | **241 (97–631)** | **11.4 (10.1–12.9)** |
| GW170817, LVK 2018 | 190 (70–580) | |

```bash
gwtc_analysis neutron_star_eos
```

### `spin_population`: BBH Effective-Spin Population

`spin_population` infers the distribution of the effective spin χ_eff of the binary black holes, and its correlation with the mass ratio q, hierarchically over the 137 events of the `hubble_constant` selection, with the LVK search-sensitivity injections for the selection effects. The model is χ_eff | q ~ N(μ₀ + α (q − 0.5), σ) on [−1, 1] ([Callister et al. 2021](https://arxiv.org/abs/2106.00521); `--no-correlation` for α = 0), with the other spin components isotropic and the masses and redshifts of the `rates` population.

| Parameter (GWTC-4.0 setup) | Median (90%) |
|---|---|
| μ₀, mean χ_eff at q = 0.5 | 0.18 (0.11–0.23) |
| σ | 0.077 (0.067–0.097) |
| α, slope with q | −0.44 (−0.58 to −0.24), P(α < 0) > 0.999 |
| fraction with χ_eff < 0 | 0.35 (0.26–0.43); GWTC-4.0: 0.24–0.42 |

Unequal-mass binaries have larger effective spins, as Callister et al. found; [GWTC-4.0](https://arxiv.org/abs/2508.18083) also finds evidence for a χ_eff–q correlation. A run takes about 30 minutes.

```bash
gwtc_analysis spin_population
```

### `parameters_estimation`: Shared Defaults And Overrides

The `parameters_estimation` workflow now separates shared plotting defaults from
per-product overrides.

- Shared defaults apply to both the whitened overlay and the q-transform:
  `--start`, `--stop`, `--fmin`, `--fmax`
- Overlay-only overrides affect only the whitened waveform overlay:
  `--overlay-start`, `--overlay-stop`, `--overlay-fmin`, `--overlay-fmax`
- Q-transform-only overrides affect only the time-frequency panel:
  `--q-start`, `--q-stop`, `--q-fmin`, `--q-fmax`, `--q-fscale {linear,log}`

If an override is omitted, the corresponding shared default is used.

When the posterior is BNS-like (median `chirp_mass` < 5 M☉), the workflow
automatically switches the overlay and q-transform windows to a BNS profile
(longer windows, wider frequency range) for any parameter you did **not** set
explicitly on the CLI. An explicit `--overlay-*` / `--q-*` value always wins.

### `parameters_estimation`: Higher multipoles and precession

From GWTC-4.0 on, the PE samples store the network SNR of the (3,3), (4,4) and (2,1) multipoles beyond the (2,2) one ([Mills & Fairhurst 2021](https://arxiv.org/abs/2007.04313)) and the precession SNR ρ_p ([Fairhurst et al. 2020](https://arxiv.org/abs/1908.05707)). The report summarizes them for every label of the file:

- the median and 90% interval of each SNR, written to `<event>_multipoles_precession.tsv`, and a plot of the selected label's posteriors;
- a comparison with noise alone: without the effect, ρ follows a Rayleigh distribution, with P(ρ > 2.1) = 11% and P(ρ > 3) = 1%. Medians ≥ 3 count as clear evidence and 2.1–3 as a hint;
- a note when the waveform models disagree.

GW231123, for example, shows a clear (4,4) multipole with the default `C00:Mixed` label, but its (3,3) SNR ranges from 2.1 to 10.5 across the waveform models. The GWTC-1 to GWTC-3 files do not store these SNRs, and the report says so.

### `parameters_estimation`: Matched-filter SNR

For each detector present in the PE file, the workflow also produces a
**matched-filter SNR** time series `|ρ(t)|`: the maximum-likelihood projected
waveform is matched-filtered against the detector strain, and the peak should
sit at the coalescence time and rise to the detector's recovered SNR.

- The strain is conditioned the canonical PyCBC way before filtering (high-pass
  at 15 Hz, resampled to a 2048 Hz grid, edges cropped of filter transients), so
  the off-source `|ρ(t)|` has unit-scale RMS (~0.7). A normalization guard warns
  if it strays from that range.
- **Short (BBH-like) signals only.** A single maximum-likelihood template cannot
  coherently recover a long BNS inspiral — over the many thousands of inspiral
  cycles, small parameter/phase differences accumulate and the SNR is lost
  (reliable BNS recovery requires a template bank, not just more strain). When the
  template is longer than the available conditioned data, the matched-filter SNR
  is **skipped with a warning**; all other plots are still produced. This is *not*
  fixable by fetching a longer strain segment.

### `build_unofficial_pe`: Unofficial Bundle Workflow

For supported special-case events such as `GW170817`, `build_unofficial_pe`
creates a PESummary-compatible `PEDataRelease` bundle from locally cached source
products such as posterior samples, PSDs, and skymaps.

- `GW170817` is reconstructed from public GWTC-1 releases only, downloaded into
  `~/.gwcache` the first time they are needed:

  | Product | DCC release | File |
  |---|---|---|
  | Posterior samples (IMRPhenomPv2_NRTidal, low and high spin) | [LIGO-P1800370](https://dcc.ligo.org/LIGO-P1800370/public) | `GW170817_GWTC-1.hdf5` |
  | PSDs (H1, L1, V1) | [LIGO-P1900011](https://dcc.ligo.org/LIGO-P1900011/public) | `GWTC1_GW170817_PSDs.dat` |
  | Calibration uncertainty envelopes | [LIGO-P1900040](https://dcc.ligo.org/LIGO-P1900040/public) | `GWTC1_GW170817_CalEnv/GWTC1_GW170817_{H,L,V}_CalEnv.txt` |
  | Skymap | [LIGO-P1800381](https://dcc.ligo.org/LIGO-P1800381/public) | `GW170817_skymap.fits.gz` |

- The public samples carry no polarization, phase, coalescence time or
  likelihood. The builder draws `psi` and `phase` from their priors, sets
  `geocent_time` to the trigger time and derives a synthetic ranking
  `log_likelihood`. For the sample it ranks first (the "maxL" sample used by the
  strain overlay), it then fits `geocent_time`, `phase` and `psi` to the GWOSC
  strain by maximizing the coherent network likelihood with the public PSDs,
  keeping the intrinsic parameters and sky position. The build log reports the
  recovered matched-filter SNRs (about H1 18, L1 25, network 31 for GW170817).
  All other samples keep their prior-drawn extrinsic values.
- Detector arrival times (`H1_time`, `L1_time`, `V1_time`) are derived from
  `geocent_time`, `ra` and `dec` with LAL detector delays.
- If a source file is missing and cannot be downloaded, the bundle is not built
  and strain overlays will not proceed from this special-case path. The strain
  fit needs network access to GWOSC; if it fails, the bundle is still built and
  the log warns that the overlay will not be coherent.
- The built bundle is cached with a recipe fingerprint
  (`<bundle>.recipe.json`); it is rebuilt when the recipe or a source file changes.

- Use `gwtc_analysis build_unofficial_pe --src-name GW170817` to
  build or reuse the cached unofficial bundle explicitly.
- Use `--force` to rebuild the bundle even if the cached output is up to date.
- `parameters_estimation` keeps its transparent fallback for supported special
  cases, but the explicit builder is the recommended way to prepare an uploadable
  bundle for S3 or local inspection.

### `parameters_estimation`: Missing PSDs and skymaps in official releases

A few official PE files ship without noise PSDs (empty `psds` groups) or
skymaps. The strain overlay needs a PSD to whiten the data, so the policy is:
if a PE label has no PSD (or no skymap), the event is looked up in a registry of
public supplementary releases (`gwtc_analysis/pe_supplements.py`), usually the
data release of the event's discovery paper, and the missing product is taken
from there and attached to every label that lacks one. Only the PSD group is
read (over HTTP range requests where the server allows it), and the result is
cached in `~/.gwcache/psd_supplements`. The log names the source and its caveat.

| Event | Missing in | Supplementary source | Caveat |
|---|---|---|---|
| `GW230529_181500` | GWTC-4.0 / 4.1 (PSDs and skymaps) | Discovery release, [Zenodo 10845779](https://zenodo.org/records/10845779) | L1 PSD identical in all 15 discovery runs |
| `GW190425_081805` | GWTC-2.1 (PSDs) | Discovery release, [LIGO-P2000026](https://dcc.ligo.org/LIGO-P2000026/public) | PSDs of the earlier LALInference analysis, not the GWTC-2.1 bilby ones |
| `GW200105_162426` | GWTC-3 v1/v2 (PSDs; fixed in v3) | Discovery release, [LIGO-P2100143](https://dcc.ligo.org/LIGO-P2100143/public) | Agree with the GWTC-3 v3 PSDs to ~2% median, mostly differing on lines |

Events with missing PSDs and no registered supplement are reported in the log,
and the overlay then whitens with a PSD estimated from the strain.

### Choosing `--pe-label` and `--waveform-engine`

The **parameter estimation** workflow distinguishes between **which PE label is used** (to read posteriors and metadata from the PE file) and **which waveform engine is requested** (to synthesize a time-domain signal for strain overlays).

**`--pe-label` (PE samples / posteriors)**

* Selects the PE label used to read posterior samples and associated metadata (e.g. `C00:Mixed`, `C00:IMRPhenomXPHM-SpinTaylor`, `C00:SEOBNRv5PHM`).
* If explicitly provided, this choice takes priority.

**`--waveform-engine` (waveform engine for strain overlay)**

* Selects the waveform generator used to build the time-domain waveform for strain overlays (engine name, e.g. `IMRPhenomXPHM`).
* This is an engine name and does not have to exactly match a PE label stored in the file.

---

#### Automatic label selection rules

1. **If `--pe-label` is explicitly provided**
   * That label is used for posterior plots.
   * The same label is also used as the source of PSDs and maximum-likelihood parameters for strain overlays.

2. **If `--pe-label` is *not* provided but `--waveform-engine` *is***
   * The tool selects the PE label whose waveform string best matches the requested engine by **substring match** on the PE label string (no hardcoded waveform mappings).
   * The selected label is then used consistently for posteriors, PSDs/detector lists, and maximum-likelihood parameters.

   Example:
   * `--waveform-engine IMRPhenomXPHM`
   * → selects `C00:IMRPhenomXPHM-SpinTaylor` if present in the PE file

3. **If neither option is provided**
   * The `Mixed` PE label is used for the posteriors (the plain `Cxx:Mixed` one when there are variants such as `Mixed:NSBH:*`).
   * The strain overlay needs a PSD, which the `Mixed` labels of the Zenodo releases do not carry: it then uses the `IMRPhenomXPHM` label (e.g. `C00:IMRPhenomXPHM-SpinTaylor`), else the first label with a PSD.
   * Files without a `Mixed` label (the GW170817 bundle) use their first label.

---

#### Waveform synthesis fallback

If the requested waveform engine cannot be instantiated (e.g. unsupported parameter range), the tool:

* logs a clear warning
* falls back to an alternative engine when possible
* explicitly reports both the requested and the actually used engine in the logs and plot titles

This ensures robustness while keeping model choices transparent.

### `parameters_estimation`: Tidal deformability of neutron stars

As the binary spirals in, each neutron star is distorted by the tidal field of
its companion. How easily it deforms is measured by the dimensionless tidal
deformability

```text
Λ = (2/3) k₂ (R c² / G m)⁵
```

where k₂ is the Love number and R the radius. Λ is 0 for a black hole and of
order 100–1000 for a 1.4 M☉ neutron star. Through the R⁵ factor it depends
steeply on the equation of state (EOS).

The deformation drains orbital energy and speeds up the end of the inspiral.
This adds a phase term at 5PN order, mostly above a few hundred Hz. The data
therefore constrain mainly one mass-weighted combination,

```text
Λ̃ = (16/13) [(m₁ + 12 m₂) m₁⁴ Λ₁ + (m₂ + 12 m₁) m₂⁴ Λ₂] / (m₁ + m₂)⁵
```

(`lambda_tilde`), while the second combination, δΛ̃ (`delta_lambda`), is
essentially unconstrained.

**Tidal samples are only in the tidal-waveform labels.** The default `Mixed`
label of BNS and NSBH events carries no tidal parameters, so select a tidal
label explicitly with `--pe-label`. Those labels provide `lambda_1`, `lambda_2`,
`lambda_tilde` and `delta_lambda`:

```bash
gwtc_analysis parameters_estimation --src-name GW170817 \
  --pe-label C02:IMRPhenomPv2_NRTidal-LowSpin \
  --pe-vars lambda_tilde delta_lambda lambda_1 lambda_2 \
  --pe-pairs lambda_1:lambda_2 chi_eff:lambda_tilde mass_ratio:lambda_tilde
```

The log lists the labels available in the PE file, and the report flags any
requested variable that the selected label does not have.

| Event | Tidal labels | Λ̃: median, 90% upper bound |
|---|---|---|
| `GW170817` (BNS, [unofficial bundle](#build_unofficial_pe-unofficial-bundle-workflow)) | `C02:IMRPhenomPv2_NRTidal-LowSpin`, `-HighSpin` | 406, ≤ 793 (low spin); 328, ≤ 746 (high spin) |
| `GW190425_081805` (BNS) | `C01:IMRPhenomPv2_NRTidal:LowSpin`, `:HighSpin` | 398, ≤ 1247 (low spin); 986, ≤ 2063 (high spin) |
| `GW200105_162426`, `GW200115_042309` (NSBH) | `C01:IMRPhenomNSBH:*`, `C01:SEOBNRv4_ROM_NRTidalv2_NSBH:*` (GW200115 also `C01:Mixed:NSBH:*`) | Λ₂ uninformative (see below) |
| `GW230529_181500` (NSBH, mass-gap primary) | `C00:IMRPhenomNSBH`, `C00:SEOBNRv4_ROM_NRTidalv2_NSBH`, `C00:IMRPhenomPv2_NRTidalv2` | Λ₂ uninformative |

Values are computed from the posterior samples the tool reads. Only GW170817 has
enough high-frequency SNR for a real measurement: its Λ̃ rules out the stiffest
EOSs and points to radii of about 11–13 km for 1.4 M☉ stars. GW190425 was seen
essentially by one detector at lower SNR, so its bound is weak.

**NSBH events: read `lambda_2`, not `lambda_tilde`.** The NSBH models fix
Λ₁ = 0 for the black hole, and the heavy primary dominates the mass weighting,
so `lambda_tilde` comes out small (median ~20–150) *whatever the neutron star
is*. The informative quantity is `lambda_2`, and for these events it stays close
to its flat 0–5000 prior (median ~2500): the neutron star was swallowed without
measurable tidal disruption. `C00:IMRPhenomPv2_NRTidalv2` for GW230529 instead
treats both objects as neutron stars (Λ₁ ≠ 0).

**Spin and tides are correlated.** The tidal effect does not depend on spin
directly, but the measurement does:

- The aligned spin (`chi_eff`, entering at 1.5PN) is degenerate with the mass
  ratio. Λ depends steeply on mass, so a wider spin prior spreads q and hence
  Λ₁ and Λ₂. This is why LVK publishes a **LowSpin** prior (|χ| ≤ 0.05, like the
  Galactic double neutron stars) and a **HighSpin** prior (|χ| ≤ 0.89). Compare
  the two labels, and look at the `chi_eff:lambda_tilde` and
  `mass_ratio:lambda_tilde` pair plots.
- A spinning neutron star also has a spin-induced quadrupole moment (2PN), which
  depends on the EOS (κ = 1 for a black hole, ~2–14 for a neutron star). The
  waveform models tie it to Λ through quasi-universal Love–Q relations instead of
  fitting it separately.
- Tidal torques do not spin up the stars: viscosity is far too low for tidal
  locking, so the models assume irrotational stars.

BNS signals are long, so a run with the strain overlays takes several minutes
(about 7 min for GW170817).

---

## Testing

```bash
gwtc_analysis search_skymaps --catalogs GWTC-4 --ra-deg 265.0 --dec-deg -46.0 --prob 0.6 --data-repo s3
gwtc_analysis event_selection --catalogs GWTC-4
gwtc_analysis event_selection --catalogs ALL --preset hierarchical --out-plot hierarchical.png
gwtc_analysis catalog_statistics --catalogs GWTC-4 --data-repo s3
gwtc_analysis catalog_statistics --catalogs GWTC-5 --data-repo zenodo
gwtc_analysis build_unofficial_pe --src-name GW170817
gwtc_analysis bright_siren
gwtc_analysis area_law
gwtc_analysis parameters_estimation --src-name GW231223_032836 --data-repo zenodo
gwtc_analysis parameters_estimation --src-name GW170817 --overlay-start 0.2 --overlay-stop 0.2 --overlay-fmax 1000 --q-start 2 --q-stop 2 --q-fmax 1000 --q-fscale log
```

---

## LIGO–Virgo–KAGRA (LVK)

- LIGO Scientific Collaboration: https://www.ligo.org
- Virgo Collaboration: https://www.virgo-gw.eu
- KAGRA Collaboration: https://gwcenter.icrr.u-tokyo.ac.jp/en

---

## Software Stack

- **GWpy** – detector strain handling and time-series analysis: https://gwpy.github.io
- **GWOSC** – public access to gravitational-wave data and metadata: https://gwosc.org
- **pesummary** – parameter-estimation posteriors handling and visualization: https://pesummary.readthedocs.io
- **ligo.skymap** – sky-localization map I/O and plotting: https://lscsoft.docs.ligo.org/ligo.skymap
- **icarogw** – hierarchical population likelihood and spectral-siren cosmology (`hubble_constant` mode): https://github.com/icarogw-developers/icarogw

---

## Acknowledgements

`gwtc_analysis` is developed as part of **ACME**, the [Astrophysics Centre for Multi-messenger studies in Europe](https://www.acme-astro.eu/).

This software is part of a project that has received funding from the European Union's Horizon Europe Research and innovation programme under Grant Agreement No 101131928.

<p>
  <a href="https://www.acme-astro.eu/"><img src="https://raw.githubusercontent.com/danielsentenac/gwtc_analysis/main/docs/img/acme_logo_dark.jpg" alt="ACME – Astrophysics Centre for Multimessenger studies in Europe" height="100"></a>
  &nbsp;&nbsp;
  <img src="https://raw.githubusercontent.com/danielsentenac/gwtc_analysis/main/docs/img/funded_by_eu.png" alt="Funded by the European Union" height="100">
</p>

Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or of the European Research Executive Agency (REA). Neither the European Union nor the granting authority can be held responsible for them.

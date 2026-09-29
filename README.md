# GWTC Analysis

## Overview

**GWTC Analysis** is a command-line analysis suite for exploring publicly released
**Gravitational-Wave Transient Catalogs (GWTC)** from the **LIGO–Virgo–KAGRA (LVK) Collaboration**.

The tool provides:

- Search of gravitational-wave sky localizations around a given sky position
- Visualization of parameter-estimation results for individual events
- Selection of events based on physical constraints (masses, distance)
- Global catalog statistics, including detector-network participation and sky-localization performance
- BNS, NSBH and BBH merger-rate estimates (R = N / ⟨VT⟩) from the catalogs and the LVK search-sensitivity injections
- Hubble-constant estimate from the BBH mass spectrum (spectral siren, with icarogw)

All gravitational-wave data products are retrieved from the **Gravitational Wave Open Science Center (GWOSC)**,
or from supported alternative repositories (Zenodo / S3 / Galaxy collections).

📖 **Documentation:** https://danielsentenac.github.io/gwtc_analysis/ (user guide, methods of the rates and Hubble-constant modes, references).

> **New:** the **GWTC-5.0** catalog (O4b observing run) is now available and fully supported — use the `GWTC-5` catalog key.

---

## Containerized Distribution (Docker)

`gwtc_analysis` is distributed as a ready-to-use **Docker container** named **`gwtc-tool`**.

- **Docker image name:** `gwtc-tool`
- **Docker Hub repository:** https://hub.docker.com/r/danielsentenac/gwtc-tool/

Using the Docker image is recommended for reproducibility, portability, and integration with workflow systems
(e.g. Galaxy, CI pipelines).

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
| `GWTC-5` | GWTC-5.0 | O4b, plus 5 events of 6–8 April 2024, just before O4b | 161 | [Zenodo 20348005](https://zenodo.org/records/20348005) (part 1, with the skymaps) and [20348006](https://zenodo.org/records/20348006) (part 2) |
| `ALL` | all the catalogs above | O1 to O4b | | |

All the keys are confident catalogs: every event of their GWOSC lists has p_astro ≥ 0.5 (the
re-analysed O1–O2 events of GWTC-2.1 carry no p_astro value). Event counts of the GWOSC lists in
September 2026; the Zenodo version read is resolved at run time (see below).

`GWTC-5` resolves to the GWOSC `GWTC-5.0` endpoint and to the Zenodo records below.

---

## Command-Line Interface (CLI)

```text
usage: gwtc_analysis [-h] MODE ...

positional arguments:
  MODE
    catalog_statistics
    rates
    hubble_constant
    event_selection
    search_skymaps
    parameters_estimation
    build_unofficial_pe
    zenodo_releases
```

Each mode has its own help:

```bash
python -m gwtc_analysis.cli <MODE> -h
```

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
python -m gwtc_analysis.cli zenodo_releases --catalogs GWTC-3 GWTC-4   # list the versions
python -m gwtc_analysis.cli search_skymaps --catalogs GWTC-3 --ra-deg 40 --dec-deg -30 --zenodo-version GWTC-3=v2
python -m gwtc_analysis.cli parameters_estimation --src-name GW200105_162426 --zenodo-version GWTC-3=v2
```

Skymap tarballs are cached per Zenodo record (`.cache_gwosc/zenodo_<record>_<file>`), and the PE index is rebuilt when the selected records change. If zenodo.org is unreachable, the cached version listing is used; with no cache at all, the latest version falls back to the record listed above.

---

## Usage

This tool is designed to run either on your laptop as a docker image or conda package, or on several user-friendly platforms:

- [docker image](https://hub.docker.com/repository/docker/danielsentenac/gwtc-tool)
- [conda package](https://anaconda.org/conda-forge/gwtc_analysis) on conda-forge: `conda install -c conda-forge gwtc_analysis`
- [Galaxy tool](https://usegalaxy.org)
- [MMODA-LIGO-VIRGO-KAGRA service](https://www.astro.unige.ch/mmoda/)

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
| `--catalogs` | `` | Catalog keys, space-separated (e.g. GWTC-1 GWTC-2.1 GWTC-3 GWTC-4 GWTC-5). ALL key takes them all. |
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

### `hubble_constant`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--stages` | `['prepare', 'sample', 'combine', 'report']` | Stages to run (default: all). |
| `--workdir` | `hubble_constant_run` | Work directory (inputs, runs, posterior). |
| `--out-report` | `hubble_constant.html` | Output HTML report path. |
| `--out-summary` | `hubble_constant.tsv` | Output TSV of the posterior quantiles. |
| `--sensitivity-release` | `gwtc4` | LVK search-sensitivity release (and matching catalogs and runs): gwtc4 = GWTC-4.0 cumulative, semi-analytic O1+O2 + real O3+O4a injections; gwtc5 = GWTC-5.0 cumulative, semi-analytic O1+O2 + real O3+O4a+O4b injections. |
| `--sensitivity-file` | `` | Local LVK injection mixture file (semi-analytic O1+O2 + real) instead of the release's. |
| `--far-threshold` | `0.25` | FAR threshold [1/yr] for the events and the real injections. |
| `--snr-threshold` | `10.0` | Network SNR threshold for the semi-analytic O1+O2 injections. |
| `--min-mass` | `3.0` | Minimum source-frame mass [Msun] of both components (potential neutron stars excluded). |
| `--exclude` | `['GW231123_135430', 'GW200105_162426']` | Events left out. |
| `--pe-cache` | `` | PE cache directory (files/, samples/, index/); default ~/.cache_gwtc_analysis/pe_catalog or $GWTC_PE_CACHE. |
| `--keep-pe-files` | `False` | Keep the full PE files after extraction. |
| `--seeds` | `[1]` | One sampler run per seed. |
| `--parallel` | `1` | Seeds run at the same time on this machine (each with --npool processes; logs in <workdir>/logs). |
| `--nlive` | `100` | dynesty live points per run. |
| `--npool` | `4` | Worker processes per run: random walks of one seed run at the same time. |
| `--naccept` | `60` | dynesty accepted steps per MCMC walk. |
| `--pe-samples` | `1500` | PE samples per event. |
| `--inj-fraction` | `0.1` | Random fraction of the found injections used (1 = all; the result stays unbiased). |
| `--icarogw-python` | `` | Python interpreter of the icarogw environment (default: the current one). |

### `event_selection`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--catalogs` | `` | Catalog keys, space-separated (e.g. GWTC-1 GWTC-2.1 GWTC-3 GWTC-4 GWTC-5). ALL key takes them all. |
| `--out-selection` | `event_selection.tsv` | Output TSV path for the selected events. |
| `--m1-min` | `` | Minimum primary mass (source frame). |
| `--m1-max` | `` | Maximum primary mass (source frame). |
| `--m2-min` | `` | Minimum secondary mass (source frame). |
| `--m2-max` | `` | Maximum secondary mass (source frame). |
| `--dl-min` | `` | Minimum luminosity distance (Mpc). |
| `--dl-max` | `` | Maximum luminosity distance (Mpc). |

### `search_skymaps`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--catalogs` | `` | Catalog keys, space-separated (e.g. GWTC-1 GWTC-2.1 GWTC-3 GWTC-4 GWTC-5). ALL key takes them all. |
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
| `--pe-label` | `` | PE label used to select posterior samples and metadata. If omitted and --waveform-engine is provided, the tool selects the closest PE label by substring match in the PE label. If both are omitted, defaults to Mixed. |
| `--waveform-engine` | `` | Waveform engine used to generate a time-domain waveform for strain overlay. If omitted, a sensible default engine is used for overlays. |

### `build_unofficial_pe`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--src-name` | `` | Source event name (e.g. GW170817). |
| `--cache-dir` | `.cache_gwosc` | Cache root where unofficial_pe/<bundle>.h5 will be written. |
| `--force` | `False` | Force rebuilding the unofficial bundle even if a cached copy already exists and is up to date. |

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

Outputs: `--out-rates` (TSV per population), `--out-events` (TSV of the events counted), `--out-report` (HTML report with the observed and the selection-corrected primary-mass distributions).

```bash
python -m gwtc_analysis.cli rates                               # GWTC-5.0 injections (O3 + O4a + O4b)
python -m gwtc_analysis.cli rates --sensitivity-release gwtc4   # GWTC-4.0 injections (O3 + O4a)
```

| Release | Candidates | BNS | NSBH | BBH at z = 0.2 |
|---|---|---|---|---|
| `gwtc5` (O3–O4b, 2.5 yr) | 258 | 26 [4, 87] | 33 [13, 67] | 25 [23, 28] |
| `gwtc4` (O3–O4a, 1.7 yr) | 154 | 43 [6, 143] | 54 [21, 109] | 26 [23, 30] |

Rates in Gpc⁻³ yr⁻¹, median [90%]. They are consistent with the LVK population papers: [GWTC-5.0](https://arxiv.org/abs/2605.27226) (BBH 27.5–49.4 at z = 0.2 for masses 2.5–200 M☉), [GWTC-4.0](https://arxiv.org/abs/2508.18083) (z = 0: BNS 7.6–250, NSBH 9.1–84, BBH 14–26) and [GWTC-3](https://arxiv.org/abs/2111.03634).

### `hubble_constant`: Hubble Constant From The BBH Mass Spectrum

`hubble_constant` measures H₀ with the **spectral-siren** method. Each event gives its luminosity distance D_L and its detector-frame masses m_det = m_src (1+z), but not its redshift. The redshift comes from the population: for a trial H₀, every D_L gives a z and every m_det a source-frame mass, and only the right H₀ makes the events near and far fall on one distance-independent mass distribution. The *Power Law + Peak* (PLP) mass model and the Madau–Dickinson rate evolution are therefore fitted together with H₀ (flat ΛCDM, Ω_m = 0.3065), with the hierarchical likelihood of [icarogw](https://github.com/icarogw-developers/icarogw) and bilby/dynesty. Selection effects are corrected with the LVK search-sensitivity injections.

The default setup reproduces the PLP measurement of the [GWTC-4.0 cosmology paper](https://arxiv.org/abs/2509.04348), H₀ = 112.7 (+51.0 / −35.9) km/s/Mpc:

- **Events**: BBHs of O1–O4a from the GWOSC confident and marginal lists, lowest FAR below 0.25/yr (`--far-threshold`), both source-frame masses above 3 M☉ (`--min-mass`), GW231123 and GW200105 left out (`--exclude`): 136 events.
- **PE samples**: from the Zenodo PE releases, `C01:IMRPhenomXPHM` up to O3 and `C00:IMRPhenomXPHM-SpinTaylor` in O4a, reduced to (m1_det, m2_det, D_L). The PE distance prior of each event is read from its file (D_L² up to O3, uniform in source-frame comoving volume in O4a) and divided out.
- **Injections** (`--sensitivity-release`): `gwtc4` (default), the GWTC-4.0 semi-analytic O1+O2 + real O3+O4a mixture ([Zenodo 16740128](https://zenodo.org/records/16740128)), found when the semi-analytic SNR exceeds 10 (`--snr-threshold`) or the lowest search FAR is below the threshold; `gwtc5` adds GWTC-5.0 and O4b ([Zenodo 19500052](https://zenodo.org/records/19500052)) and is not yet validated against a published result.
- **Priors**: those of the paper (Tables 3 and 6), H₀ uniform in [10, 200] km/s/Mpc.

The work is split into stages (`--stages`, all by default), sharing `--workdir`:

| Stage | Does | Cost |
|---|---|---|
| `prepare` | selects the events, downloads their PE files (restartable; only the extracted samples are kept in `--pe-cache` unless `--keep-pe-files`), prepares the injections → `inputs.h5`, `events.tsv` | ~35 GB of downloads the first time |
| `sample` | one dynesty run per `--seeds` value (`--nlive`, `--npool`, `--naccept`, `--pe-samples`, `--inj-fraction`), `--parallel` of them at a time on this machine (logs in `<workdir>/logs`); resumable from its checkpoint | hours per run |
| `combine` | merges the runs → `posterior.tsv`, `corner.png`, `summary.json`, with the effective numbers of injections and PE samples over the posterior (icarogw's stability criteria) | minutes |
| `report` | `--out-report` (HTML) and `--out-summary` (TSV of the posterior quantiles) | seconds |

**icarogw.** `gwtc_analysis/h0_icarogw.py` is a driver of icarogw, not a modified copy: icarogw is used as installed, through its public API. icarogw provides the hierarchical likelihood (PE and injection reweighting, selection term, scale-free rate marginalisation, effective-sample-size checks), the population models (`massprior_PowerLawPeak` with the `m1m2_conditioned_lowpass` smoothing, `rateevolution_Madau`, `FlatLambdaCDM_wrap`, combined by `CBC_vanilla_rate`) and the detector-frame conversion for each trial H₀. The driver reads `inputs.h5` into icarogw's `posterior_samples` and `injections` objects, chooses the model components and the priors (Tables 3 and 6 of the paper), runs bilby/dynesty, merges the runs and computes the diagnostics with icarogw's own methods. The analysis choices made here, outside icarogw, are the input preparation in `hubble_constant.py` (event selection, PE distance prior read from each file, injection draw density carried to the detector frame with the spin part divided out and the mixture weights applied) and three settings of the driver: at least 10 effective PE samples per event (icarogw's default is 20), at least 4 × N_events effective injections (icarogw's default), and the `--inj-fraction` subset of the injections.

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
python -m gwtc_analysis.cli hubble_constant --stages prepare
python -m gwtc_analysis.cli hubble_constant --stages sample combine report \
    --icarogw-python ~/.conda/envs/icarogw/bin/python --seeds 1 2 3 4 --parallel 2 --npool 2
```

**Seeds.** Each seed is an independent dynesty run (`result/plp_seed<N>_result.json`); `combine` merges all the finished ones, weighted by their evidence. All runs sample the same likelihood: the PE samples are shuffled once in `prepare` and the injection subset is drawn with a fixed seed, and `run_settings.json` refuses runs with other `--nlive`, `--pe-samples` or `--inj-fraction` values in the same work directory. Launching again resumes the interrupted runs from their checkpoint and skips the finished ones. A lock file (`result/plp_seed<N>.lock`) prevents the same seed from running twice at once; interrupting the launcher (Ctrl-C) stops its runs after they write their checkpoint.

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

With 10 runs of 100 live points, 1500 PE samples per event and 10% of the found injections, the result is H₀ = 119.3 (+46.1 / −34.9) km/s/Mpc [68%], 90%: 62.9–186.1, and the mass peak μ_g = 27.8 (+4.0 / −4.7) M☉, against 28.6 (+3.9 / −4.9) M☉ in the paper. H₀ is strongly anti-correlated with μ_g, and the upper part of its interval depends on the prior bound of 200 km/s/Mpc. The richer mass models of the paper give tighter results (MLTP 77.1, FullPop-4.0 76.4 km/s/Mpc); only PLP is implemented here.

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

- Use `python -m gwtc_analysis.cli build_unofficial_pe --src-name GW170817` to
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
   * The default `Mixed` PE label is used.

---

#### Waveform synthesis fallback

If the requested waveform engine cannot be instantiated (e.g. unsupported parameter range), the tool:

* logs a clear warning
* falls back to an alternative engine when possible
* explicitly reports both the requested and the actually used engine in the logs and plot titles

This ensures robustness while keeping model choices transparent.

---

## Testing

```bash
python -m gwtc_analysis.cli search_skymaps --catalogs GWTC-4 --ra-deg 265.0 --dec-deg -46.0 --prob 0.6 --data-repo s3
python -m gwtc_analysis.cli event_selection --catalogs GWTC-4
python -m gwtc_analysis.cli catalog_statistics --catalogs GWTC-4 --data-repo s3
python -m gwtc_analysis.cli catalog_statistics --catalogs GWTC-5 --data-repo zenodo
python -m gwtc_analysis.cli build_unofficial_pe --src-name GW170817
python -m gwtc_analysis.cli parameters_estimation --src-name GW231223_032836 --data-repo zenodo
python -m gwtc_analysis.cli parameters_estimation --src-name GW170817 --overlay-start 0.2 --overlay-stop 0.2 --overlay-fmax 1000 --q-start 2 --q-stop 2 --q-fmax 1000 --q-fscale log
```

---

## LIGO–Virgo–KAGRA (LVK)

- LIGO Scientific Collaboration: https://www.ligo.org
- Virgo Collaboration: https://www.virgo-gw.eu
- KAGRA Collaboration: https://gwcenter.icrr.u-tokyo.ac.jp/en

---

## Software Stack

- **GWpy** – detector strain handling and time-series analysis: https://gwpy.github.io
- **GWOSC** – public access to gravitational-wave data and metadata: https://www.gw-openscience.org
- **pesummary** – parameter-estimation posteriors handling and visualization: https://pesummary.readthedocs.io
- **ligo.skymap** – sky-localization map I/O and plotting: https://lscsoft.docs.ligo.org/ligo.skymap

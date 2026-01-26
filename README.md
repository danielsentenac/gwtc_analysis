# GWTC Analysis

## Overview

**GWTC Analysis** is a command-line analysis suite for exploring publicly released
**Gravitational-Wave Transient Catalogs (GWTC)** from the **LIGO–Virgo–KAGRA (LVK) Collaboration**.

The tool provides:

- Search of gravitational-wave sky localizations around a given sky position
- Visualization of parameter-estimation results for individual events
- Selection of events based on physical constraints (masses, distance)
- Global catalog statistics, including detector-network participation and sky-localization performance

All gravitational-wave data products are retrieved from the **Gravitational Wave Open Science Center (GWOSC)**,
or from supported alternative repositories (Zenodo / S3 / Galaxy collections).

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

Catalog identifiers are **case-sensitive**:

| Catalog name | Description |
|---|---|
| `GWTC-1` | Confident subset of GWTC-1 |
| `GWTC-2.1` | Confident subset of GWTC-2.1 |
| `GWTC-3` | Confident subset of GWTC-3 |
| `GWTC-4` | GWTC-4 public release |
| `ALL` | Expands to all catalogs above |

---

## Command-Line Interface (CLI)

```text
usage: gwtc_analysis [-h] MODE ...

positional arguments:
  MODE
    catalog_statistics
    event_selection
    search_skymaps
    parameters_estimation
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

- The Zenodo portal (official catalogs PE/skymaps tarballs) at https://zenodo.org/records/8177023|17014085|6513631.
- A s3 Minio bucket called gwtc on  https://minio-dev.odahub.fr
- Galaxy collections under the name GWTC at https://usegalaxy.org

---

## Usage

This tool is designed to run either on your laptop as a docker image or conda package, or on several user-friendly platforms:

- [docker image](https://hub.docker.com/repository/docker/danielsentenac/gwtc-tool)
- [conda package](https://anaconda.org/channels/danielsentenac/packages/gwtc_analysis/overview)
- [Galaxy tool](https://usegalaxy.org)
- [MMODA-LIGO-VIRGO-KAGRA service](https://www.astro.unige.ch/mmoda/)

### Inputs
- Catalog selections are passed as parameters separated by space
- Data repositories accept `--data-repo` to choose where data products are read from:
	- `galaxy`: read inputs from Galaxy collections
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
python tools/gen_readme_cli_tables.py
```

<!-- CLI_TABLES_BEGIN -->

### `catalog_statistics`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--catalogs` | `` | Catalog keys, space-separated (e.g. GWTC-1 GWTC-2.1 GWTC-3 GWTC-4). ALL key takes them all. |
| `--out-events` | `catalogs_statistics.tsv` | Output TSV path (per-event table). |
| `--out-report` | `catalogs_statistics.html` | Output HTML report path. |
| `--include-detectors` | `False` | Include detector network via GWOSC v2 calls. |
| `--include-area` | `False` | Compute sky localization area Axx if skymaps are available. |
| `--area-cred` | `0.9` | Credible level for sky area: 0.9→A90, 0.5→A50, 0.95→A95. |
| `--plots-dir` | `cat_plots` | Directory for plots (default: cat_plots). |
| `--data-repo` | `zenodo` | Where to read data from: galaxy \| zenodo \| s3. |

### `event_selection`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--catalogs` | `` | Catalog keys, space-separated (e.g. GWTC-1 GWTC-2.1 GWTC-3 GWTC-4). ALL key takes them all. |
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
| `--catalogs` | `` | Catalog keys, space-separated (e.g. GWTC-1 GWTC-2.1 GWTC-3 GWTC-4). ALL key takes them all. |
| `--ra-deg` | `` | Right ascension (deg). |
| `--dec-deg` | `` | Declination (deg). |
| `--prob` | `0.9` | Credible-level threshold (0–1). Common values: 0.9, 0.5, 0.95. |
| `--skymap_label` | `Mixed` | Label selector used to filter skymap (default: Mixed). |
| `--out-events` | `search_skymaps.tsv` | Output TSV file (default: search_skymaps.tsv). |
| `--out-report` | `search_skymaps.html` | Optional output HTML report path for hits. |
| `--plots-dir` | `sky_plots` | Directory for hit plots (default: sky_plots). |
| `--data-repo` | `zenodo` | Where to read data from: galaxy \| zenodo \| s3. |

### `parameters_estimation`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--out-report` | `parameters_estimation.html` | Output HTML report path. |
| `--src-name` | `` | Source event name (e.g. GW231223_032836). |
| `--data-repo` | `zenodo` | Where to read data from: galaxy \| zenodo \| s3. |
| `--pe-vars` | `` | Extra posterior sample variables to plot (space-separated). Example: --pe-vars chi_eff chi_p luminosity_distance. |
| `--pe-pairs` | `` | Extra 2D posterior pairs to plot as 'x:y' tokens. Example: --pe-pairs mass_1_source:mass_2_source chi_eff:chi_p. |
| `--plots-dir` | `pe_plots` | Directory for output PE plots (default: pe_plots). |
| `--start` | `0.2` | Seconds before GPS time for strain window. |
| `--stop` | `0.1` | Seconds after GPS time for strain window. |
| `--fs-low` | `20.0` | Bandpass low frequency (Hz). |
| `--fs-high` | `300.0` | Bandpass high frequency (Hz). |
| `--pe-label` | `` | PE label used to select posterior samples and metadata. If omitted and --waveform-engine is provided, the tool selects the closest PE label by substring match in the PE label. If both are omitted, defaults to Mixed. |
| `--waveform-engine` | `` | Waveform engine used to generate a time-domain waveform for strain overlay. If omitted, a sensible default engine is used for overlays. |

<!-- CLI_TABLES_END -->

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
python -m gwtc_analysis.cli parameters_estimation --src-name GW231223_032836 --data-repo zenodo
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

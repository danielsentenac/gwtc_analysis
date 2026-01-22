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

Many modes accept `--data-repo` to choose where data products are read from:

- `local`: use local files only (no downloads).
- `galaxy`: read inputs from Galaxy collections (collections may expand to lists of file paths).
- `zenodo`: download official releases from Zenodo; cached under `.cache_gwosc/`.
- `s3`: read skymaps from an S3-compatible object store (optimized for large catalogs).

---

## Galaxy usage

This tool is designed to run smoothly inside **Galaxy**.

### Inputs
- Catalog selections are passed as parameters.
- Skymaps may be provided as:
  - Galaxy collections (directories or file lists)
  - Individual FITS / FITS.gz files

### Outputs
- TSV tables are produced as standard Galaxy datasets.
- HTML reports are single-file (no external assets) and render in the Galaxy UI.
- Plot images are written under the directory given by `--plots-dir`.

### Collections
When a Galaxy collection is provided:
- a single directory collection is used directly when possible
- a multi-file collection may be materialized into a working directory (symlinks preferred; copies as fallback)

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
| `--catalogs` | `` | Catalog keys. Space-separated; commas also accepted (e.g. GWTC-4,GWTC-3). |
| `--out-events` | `catalogs_statistics.tsv` | Output TSV path (per-event table). |
| `--out-report` | `catalogs_statistics.html` | Output HTML report path. |
| `--include-detectors` | `False` | Include detector network via GWOSC v2 calls. |
| `--include-area` | `False` | Compute sky localization area Axx if skymaps are available. |
| `--area-cred` | `0.9` | Credible level for sky area: 0.9→A90, 0.5→A50, 0.95→A95. |
| `--plots-dir` | `cat_plots` | Directory for plots (default: cat_plots). |
| `--data-repo` | `zenodo` | Where to read data from: galaxy \| zenodo \| s3 |

### `event_selection`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--catalogs` | `` | Catalogs space-separated: GWTC-1, GWTC-2.1, GWTC-3, GWTC-4, ALL |
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
| `--catalogs` | `` | Catalogs space-separated: GWTC-1, GWTC-2.1, GWTC-3, GWTC-4, ALL |
| `--ra-deg` | `` | Right ascension (deg). |
| `--dec-deg` | `` | Declination (deg). |
| `--prob` | `0.9` | Credible-level threshold (0–1). Common values: 0.9, 0.5, 0.95. |
| `--waveform` | `Mixed` | Waveform/approximant selector used to filter skymap filenames (default: Mixed). Use 'any' to disable filtering. |
| `--out-events` | `search_skymaps.tsv` | Output TSV file (default: search_skymaps.tsv) |
| `--out-report` | `search_skymaps.html` | Optional output HTML report path for hits. |
| `--plots-dir` | `sky_plots` | Directory for hit plots (default: sky_plots). |
| `--data-repo` | `zenodo` | Where to read data from: galaxy \| zenodo \| s3 |

### `parameters_estimation`

| Option | Default | Description |
|---|---:|---|
| `-h, --help` | `` | show this help message and exit |
| `--out-report` | `parameters_estimation.html` | Output HTML report path. |
| `--src-name` | `` | Source event name (e.g. GW231223_032836). |
| `--data-repo` | `zenodo` | Where to read data from: galaxy \| zenodo \| s3 |
| `--pe-vars` | `` | Extra posterior sample variables to plot (space-separated). Example: --pe-vars chi_eff chi_p luminosity_distance |
| `--pe-pairs` | `` | Extra 2D posterior pairs to plot as 'x:y' tokens. Example: --pe-pairs mass_1_source:mass_2_source chi_eff:chi_p |
| `--plots-dir` | `pe_plots` | Directory for output PE plots (default: pe_plots). |
| `--start` | `0.2` | Seconds before GPS time for strain window. |
| `--stop` | `0.1` | Seconds after GPS time for strain window. |
| `--fs-low` | `20.0` | Bandpass low frequency (Hz). |
| `--fs-high` | `300.0` | Bandpass high frequency (Hz). |
| `--sample-method` | `Mixed` | Posterior sample label/model selector (e.g. Mixed, IMRPhenomXPHM, SEOBNRv4PHM). |
| `--strain-approximant` | `IMRPhenomXPHM` | Waveform model used to generate time-domain waveform for strain overlay. |

<!-- CLI_TABLES_END -->

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

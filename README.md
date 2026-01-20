# GWTC Analysis

## Overview

**GWTC Analysis** is a command-line analysis suite for exploring publicly released  
**Gravitational-Wave Transient Catalogs (GWTC)** from the  
**LIGO–Virgo–KAGRA (LVK) Collaboration**.

The tool provides:

- Search of gravitational-wave sky localizations around a given sky position
- Visualization of parameter-estimation results for individual events
- Selection of events based on physical constraints (masses, distance)
- Global catalog statistics, including detector-network participation and sky-localization performance

All gravitational-wave data products are retrieved from the  
**Gravitational Wave Open Science Center (GWOSC)**.

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

The tool supports **only** the following catalog identifiers (case-sensitive):

| Catalog name | Description |
|--------------|-------------|
| **GWTC-1-confident** | Confident subset of GWTC-1 |
| **GWTC-2.1-confident** | Confident subset of GWTC-2.1 |
| **GWTC-3-confident** | Confident subset of GWTC-3 |
| **GWTC-4.0** | GWTC-4.0 public release |

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

---

## General Units and Ranges

- Right Ascension: degrees [0, 360)
- Declination: degrees [-90, +90]
- Probability threshold: [0, 1]
- Masses: solar masses (M☉)
- Distances: megaparsecs (Mpc)

---

## Mode: `search_skymaps`

Given a sky position (RA/Dec in degrees) and a list of sky-localization FITS files,
report which events contain that position above the requested credible level.

### Input Options

| Option | Description | Default |
|------|------------|---------|
| `--catalogs` | Catalog keys (space-separated; commas also accepted). | None |
| `--ra-deg` | Right ascension (degrees). | None |
| `--dec-deg` | Declination (degrees). | None |
| `--prob` | Credible-level threshold (0–1). | 0.9 |
| `--skymaps` | List of skymap FITS / FITS.gz files. | None |
| `--out-events` | Optional output TSV path for hits. | None |
| `--out-report` | Optional output HTML report path. | None |
| `--plots-dir` | Directory for hit plots. | Current directory |
| `--events-json` | Offline mode: local jsonfull-like events JSON. | None |

---

## Mode: `event_selection`

Select events by simple cuts on source-frame masses and luminosity distance.

### Input Options

| Option | Description | Default |
|------|------------|---------|
| `--catalogs` | Catalog keys (space-separated; commas also accepted). | None |
| `--out-selection` | Output TSV path for selected events. | None |
| `--m1-min` | Minimum primary mass (source frame). | None |
| `--m1-max` | Maximum primary mass (source frame). | None |
| `--m2-min` | Minimum secondary mass (source frame). | None |
| `--m2-max` | Maximum secondary mass (source frame). | None |
| `--dl-min` | Minimum luminosity distance (Mpc). | None |
| `--dl-max` | Maximum luminosity distance (Mpc). | None |
| `--events-json` | Offline mode: local jsonfull-like events JSON. | None |

---

## Mode: `catalog_statistics`

Fetch events for one or more catalogs and compute derived statistics.

### Input Options

| Option | Description | Default |
|------|------------|---------|
| `--catalogs` | Catalog keys (space-separated; commas also accepted). | None |
| `--out-events` | Output TSV path (per-event table). | None |
| `--out-report` | Output HTML report path. | None |
| `--include-detectors` | Include detector network via GWOSC v2 calls. | False |
| `--include-area` | Compute sky-localization area Axx. | False |
| `--area-cred` | Credible level for sky area (e.g. 0.9 → A90). | 0.9 |
| `--skymaps-gwtc21` | GWTC-2.1 skymaps (files or collections). | None |
| `--skymaps-gwtc3` | GWTC-3 skymaps (files or collections). | None |
| `--skymaps-gwtc4` | GWTC-4.0 skymaps (files or collections). | None |
| `--events-json` | Offline mode: local jsonfull-like events JSON. | None |

---

## Mode: `parameters_estimation`

Generate parameter-estimation plots for a single event.

### Input Options

| Option | Description | Default |
|------|------------|---------|
| `--src-name` | Source event name (e.g. `GW231223_032836`). | None |
| `--plots-dir` | Directory for output plots. | Current directory |
| `--start` | Seconds before GPS time for strain window. | 0.2 |
| `--stop` | Seconds after GPS time for strain window. | 0.1 |
| `--fs-low` | Bandpass low frequency (Hz). | 20.0 |
| `--fs-high` | Bandpass high frequency (Hz). | 300.0 |
| `--sample-method` | Posterior sample selector. | Mixed |
| `--strain-approximant` | Waveform model for strain overlay. | IMRPhenomXPHM |

### Detector Strain Data

Strain time-series data are retrieved from GWOSC using **GWpy**.
Plots are generated for all available detectors: **H1, L1, V1, K1**.

---

## LIGO–Virgo–KAGRA (LVK)

The LIGO–Virgo–KAGRA (LVK) Collaboration operates the global network of ground-based
gravitational-wave detectors and produces the GWTC catalogs used by this tool.

- LIGO Scientific Collaboration: https://www.ligo.org
- Virgo Collaboration: https://www.virgo-gw.eu
- KAGRA Collaboration: https://gwcenter.icrr.u-tokyo.ac.jp/en

---

## Software Stack

This tool relies on open-source software developed and maintained by the LVK collaboration and the broader community:

- **GWpy** – detector strain handling and time-series analysis  
  https://gwpy.github.io
- **GWOSC** – public access to gravitational-wave data and metadata  
  https://www.gw-openscience.org
- **pesummary** – parameter-estimation posteriors handling and visualization  
  https://pesummary.readthedocs.io
- **ligo.skymap** – sky-localization map I/O and plotting  
  https://lscsoft.docs.ligo.org/ligo.skymap

---

## Testing

```bash
python -m gwtc_analysis.cli search_skymaps --catalogs GWTC-4.0 --ra-deg 265.0 --dec-deg -46.0 --skymaps data/*.fits
python -m gwtc_analysis.cli event_selection --catalogs GWTC-4.0 --out-selection selected.tsv
python -m gwtc_analysis.cli catalog_statistics --catalogs GWTC-4.0 --out-events events.tsv --out-report report.html
python -m gwtc_analysis.cli parameters_estimation --src-name GW231223_032836 --plots-dir plots
```

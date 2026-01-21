**GWTC Analysis** is a command-line analysis suite for exploring publicly released
**Gravitational-Wave Transient Catalogs (GWTC)** from the
**LIGO–Virgo–KAGRA (LVK) Collaboration**.

It provides:

- Search of gravitational-wave sky localizations around a given sky position
- Visualization of parameter-estimation results for individual events
- Selection of events based on physical constraints (masses, distance)
- Global catalog statistics, including detector-network participation and sky-localization performance

All gravitational-wave data products are retrieved from the
**Gravitational Wave Open Science Center (GWOSC)**.

---

### Containerized Distribution (Docker)

`gwtc_analysis` is also distributed as a ready-to-use Docker container:

- **Image:** `danielsentenac/gwtc-tool`
- **Docker Hub:** https://hub.docker.com/r/danielsentenac/gwtc-tool/

---

### Astrophysical Sources

The GWTC catalogs contain compact binary mergers:

- Binary Black Holes (BBH)
- Binary Neutron Stars (BNS)
- Neutron Star – Black Hole systems (NSBH)

Detected by the LVK detector network:
**H1 (Hanford), L1 (Livingston), V1 (Virgo), K1 (KAGRA)**.

---

### Supported GW Catalog Names

The following catalog identifiers are supported (case-sensitive):

- `GWTC-2.1`
- `GWTC-3`
- `GWTC-4`

---

### Command-Line Interface

```text
gwtc_analysis [-h] MODE ...
```

Available modes:

- `search_skymaps`
- `event_selection`
- `catalog_statistics`
- `parameters_estimation`

---

### General Units

- Right Ascension: degrees [0, 360)
- Declination: degrees [-90, +90]
- Probability / credible level: [0, 1]
- Masses: solar masses (M☉)
- Distances: megaparsecs (Mpc)

---

### Mode: search_skymaps

Search sky-localization maps for events consistent with a given sky position.

### Input options

- `--catalogs` (required): Catalog keys (space-separated or comma-separated)
- `--ra-deg` (required): Right ascension in degrees
- `--dec-deg` (required): Declination in degrees
- `--prob`: Credible-level threshold in [0, 1] (default: 0.9)
- `--skymaps` (required): List of skymap FITS / FITS.gz files
- `--out-events`: Optional output TSV file
- `--out-report`: Optional output HTML report
- `--plots-dir`: Directory for plots (default: current directory)
- `--events-json`: Offline mode using a local events JSON

---

## Mode: event_selection

Select events using cuts on source-frame masses and luminosity distance.

### Input options

- `--catalogs` (required): Catalog keys
- `--out-selection` (required): Output TSV file
- `--m1-min`: Minimum primary mass (M☉)
- `--m1-max`: Maximum primary mass (M☉)
- `--m2-min`: Minimum secondary mass (M☉)
- `--m2-max`: Maximum secondary mass (M☉)
- `--dl-min`: Minimum luminosity distance (Mpc)
- `--dl-max`: Maximum luminosity distance (Mpc)
- `--events-json`: Offline mode using a local events JSON

---

### Mode: catalog_statistics

Compute catalog-wide statistics and summary plots.

### Input options

- `--catalogs` (required): Catalog keys
- `--out-events` (required): Output TSV table
- `--out-report` (required): Output HTML report
- `--include-detectors`: Include detector network information
- `--include-area`: Compute sky-localization area
- `--area-cred`: Credible level for sky area (default: 0.9)
- `--skymaps-gwtc21`: GWTC-2.1 skymaps
- `--skymaps-gwtc3`: GWTC-3 skymaps
- `--skymaps-gwtc4`: GWTC-4.0 skymaps
- `--events-json`: Offline mode using a local events JSON

---

### Mode: parameters_estimation

Generate parameter-estimation plots for a single event.

### Input options

- `--src-name` (required): Source event name (e.g. GW231223_032836)
- `--plots-dir`: Output directory (default: current directory)
- `--start`: Seconds before GPS time for strain window (default: 0.2)
- `--stop`: Seconds after GPS time for strain window (default: 0.1)
- `--fs-low`: Bandpass low frequency in Hz (default: 20.0)
- `--fs-high`: Bandpass high frequency in Hz (default: 300.0)
- `--sample-method`: Posterior sample selector (default: Mixed)
- `--strain-approximant`: Waveform model for strain overlay (default: IMRPhenomXPHM)

Strain data are retrieved using **GWpy** for all available detectors:
**H1, L1, V1, K1**.

---

### Software Stack

- GWpy: https://gwpy.github.io
- GWOSC: https://www.gw-openscience.org
- pesummary: https://pesummary.readthedocs.io
- ligo.skymap: https://lscsoft.docs.ligo.org/ligo.skymap

### Testing

- `python -m gwtc_analysis.cli search_skymaps --catalogs GWTC-4.0 --out-events hits.tsv --ra-deg 265.0 --dec-deg -46.0 --skymaps /data`
- `python -m gwtc_analysis.cli event_selection --catalogs GWTC-4.0 --out-selection selected.tsv`
- `python -m gwtc_analysis.cli catalog_statistics --catalogs GWTC-4.0 --out-events events.tsv --out-report report.html`
- `python -m gwtc_analysis.cli parameters_estimation --src-name GW231223_032836 --plots-dir plots`


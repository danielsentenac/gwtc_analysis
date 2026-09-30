# GWTC Analysis

**GWTC Analysis** (`gwtc_analysis`) is a command-line suite for exploring the public
**Gravitational-Wave Transient Catalogs** (GWTC) of the **LIGO–Virgo–KAGRA** (LVK) collaboration,
from GWTC-1 (2015) to GWTC-5.0 (the O4b observing run, 2025).

It reads the catalogs directly from their public releases (GWOSC, Zenodo, an S3 mirror or Galaxy
collections) and produces TSV tables, plots and self-contained HTML reports.

## What it does

| Mode | Purpose |
|---|---|
| [`catalog_statistics`](modes/catalog-statistics.md) | Catalog-wide tables and plots: masses, distances, SNR, detector networks, sky-localization areas |
| [`event_selection`](modes/event-selection.md) | Events selected by mass and distance ranges |
| [`search_skymaps`](modes/search-skymaps.md) | Events whose sky localization contains a given sky position |
| [`parameters_estimation`](modes/parameters-estimation.md) | Posterior plots of one event, whitened strain overlays, q-transforms, matched-filter SNR |
| [`build_unofficial_pe`](modes/unofficial-pe.md) | A PESummary-compatible bundle for GW170817, rebuilt from public GWTC-1 products |
| [`rates`](modes/rates.md) | BNS, NSBH and BBH merger rates, R = N / ⟨VT⟩, from the catalogs and the LVK sensitivity injections |
| [`hubble_constant`](modes/hubble-constant.md) | The Hubble constant from the binary-black-hole mass spectrum (spectral siren), with icarogw |
| `zenodo_releases` | The versions of the Zenodo catalog releases ([Data sources](data-sources.md#zenodo-release-versions)) |

## Highlights

- **Hubble constant from black holes alone.** The `hubble_constant` mode reproduces the
  *Power Law + Peak* spectral-siren measurement of the GWTC-4.0 cosmology paper [\[27\]](references.md#ref-27):
  H₀ = 106.5 (+45.0 / −34.0) km/s/Mpc against the published 105.5 (+46.4 / −35.8), with the same
  136 binary black holes, PE samples, injections and priors. See
  [Hubble constant (spectral siren)](science/spectral-siren.md).
- **Merger rates** of the three source classes, corrected for selection effects with the LVK
  search-sensitivity injections. See [Merger rates](science/merger-rates.md).
- **Robust access to the releases:** Zenodo versions resolved at run time, and missing PSDs of
  official releases filled from the events' discovery releases. See
  [Known issues in the public releases](known-data-issues.md).

![Primary masses of the confident GW events](img/mass_distribution.png)

## Where to go next

- New users: [Installation](installation.md), then [Quick start](quickstart.md).
- The physics and statistics behind the rates and H₀ modes:
  [What a GW signal measures](science/gw-signals.md) and the following pages.
- Every option of every mode: [CLI reference](cli-reference.md).
- All the papers the tool and its methods rely on: [References](references.md).

The source code is on [GitHub](https://github.com/danielsentenac/gwtc_analysis) under the MIT
license.

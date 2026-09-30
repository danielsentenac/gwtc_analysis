# catalog_statistics

Catalog-wide statistics of one or several catalogs: a table with one row per event, plots, and an
HTML report.

```bash
gwtc_analysis catalog_statistics --catalogs GWTC-4 GWTC-5 --include-detectors --include-area
```

**Inputs.** Event metadata from the GWOSC event lists (names, times, false-alarm rates, median
source-frame masses, distances, network SNR); with `--include-detectors`, the detector network of each
event from the GWOSC v2 API; with `--include-area`, the sky-localization area at the credible level
`--area-cred` (A90 by default), computed from the skymaps of `--data-repo`.

**Outputs.**

- `--out-events` (TSV): one row per event;
- `--plots-dir`: the primary against secondary mass diagram colored by SNR, histograms of the main
  parameters, the fraction of each source type (BBH, BNS, NSBH), the detector networks, and the
  cumulative distribution of the sky-localization areas;
- `--out-report`: the HTML report with the tables and plots.

## Example

```bash
gwtc_analysis catalog_statistics --catalogs GWTC-4 --include-detectors --include-area
```

GWTC-4.0 (O4a), 86 events with a PE release; every figure below was produced by this command.

![Primary against secondary source-frame mass, colored by network SNR](../img/modes/catstat_m1_m2_snr.png)

*Median source-frame masses, colored by network SNR. The heaviest event is GW231123; the two points
near the horizontal axis are the NSBH candidates.*

![Histograms of total mass, luminosity distance and network SNR](../img/modes/catstat_histograms.png)

| Source types | Detector networks |
|---|---|
| ![Source types](../img/modes/catstat_source_types.png) | ![Detector networks](../img/modes/catstat_network.png) |

*84 BBHs and 2 NSBHs; 77 events seen by the two LIGO detectors, 9 by one only (Virgo did not observe
during O4a).*

![Cumulative distribution of the 90% sky areas](../img/modes/catstat_area_cdf.png){ width="560" }

*90% sky areas from the Zenodo PE skymaps: median 2 065 deg² with two detectors, about 27 600 deg²
with a single detector, which localizes an event only to a large part of the sky.*

All options: [CLI reference](../cli-reference.md#catalog_statistics).

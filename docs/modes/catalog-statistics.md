# catalog_statistics

Catalog-wide statistics of one or several catalogs: a table with one row per event, plots, and an
HTML report.

```bash
python -m gwtc_analysis.cli catalog_statistics --catalogs GWTC-4 GWTC-5 --include-detectors --include-area
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

All options: [CLI reference](../cli-reference.md#catalog_statistics).

# event_selection

Events of the chosen catalogs whose median source-frame masses and luminosity distance fall in given
ranges. Every bound is optional.

```bash
# heavy binaries within 2 Gpc
python -m gwtc_analysis.cli event_selection --catalogs ALL --m1-min 50 --dl-max 2000

# potential neutron-star companions
python -m gwtc_analysis.cli event_selection --catalogs ALL --m2-max 3
```

| Option | Selects on |
|---|---|
| `--m1-min`, `--m1-max` | primary mass, source frame (M☉) |
| `--m2-min`, `--m2-max` | secondary mass, source frame (M☉) |
| `--dl-min`, `--dl-max` | luminosity distance (Mpc) |

The selection is written to `--out-selection` (TSV). It uses the GWOSC confident lists; for analyses
that need the marginal candidates as well, see [Known issues](../known-data-issues.md#gw200105-only-in-the-marginal-list).

All options: [CLI reference](../cli-reference.md#event_selection).

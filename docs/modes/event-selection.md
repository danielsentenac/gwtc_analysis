# event_selection

Events of the chosen catalogs whose median source-frame masses and luminosity distance fall in given
ranges. Every bound is optional.

```bash
# heavy binaries within 2 Gpc
gwtc_analysis event_selection --catalogs ALL --m1-min 50 --dl-max 2000

# potential neutron-star companions
gwtc_analysis event_selection --catalogs ALL --m2-max 3
```

| Option | Selects on |
|---|---|
| `--m1-min`, `--m1-max` | primary mass, source frame (M☉) |
| `--m2-min`, `--m2-max` | secondary mass, source frame (M☉) |
| `--dl-min`, `--dl-max` | luminosity distance (Mpc) |

The selection is written to `--out-selection` (TSV). It uses the GWOSC confident lists; for analyses
that need the marginal candidates as well, see [Known issues](../known-data-issues.md#gw200105-only-in-the-marginal-list).

## Example

```bash
gwtc_analysis event_selection --catalogs ALL --m1-min 50 --dl-max 2000
```

Heavy binaries within 2 Gpc; the TSV (`event_selection.tsv`) lists:

| event_id | catalog_key | mass_1_source | mass_2_source | luminosity_distance |
|---|---|---|---|---|
| GW191109_010717-v1 | GWTC-3 | 65.0 | 47.0 | 1290 |
| GW240519_012815-v1 | GWTC-5 | 65.0 | 39.0 | 1740 |
| GW241127_061008-v1 | GWTC-5 | 63.8 | 20.8 | 1080 |
| GW241225_082815-v1 | GWTC-5 | 55.7 | 42.2 | 1880 |

All options: [CLI reference](../cli-reference.md#event_selection).

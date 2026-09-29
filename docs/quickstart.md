# Quick start

Every mode has its own help:

```bash
python -m gwtc_analysis.cli -h
python -m gwtc_analysis.cli <MODE> -h
```

Catalogs are named by case-sensitive keys. All of them are **confident catalogs** (see below):

| Key | Catalog | Observing runs of its events | Events | PE and skymaps on Zenodo |
|---|---|---|---|---|
| `GWTC-1` | GWTC-1 | O1, O2 | 11 | the GWTC-2.1 release, [Zenodo 6513631](https://zenodo.org/records/6513631), which re-analysed O1–O2 |
| `GWTC-2.1` | GWTC-2.1 | O3a, plus 10 O1–O2 events re-analysed | 54 | [Zenodo 6513631](https://zenodo.org/records/6513631) |
| `GWTC-3` | GWTC-3 | O3b | 35 | [Zenodo 22685054](https://zenodo.org/records/22685054) |
| `GWTC-4` | GWTC-4.0 | O4a, plus GW230518 from the engineering run ER15 | 129 | [Zenodo 17602505](https://zenodo.org/records/17602505) |
| `GWTC-5` | GWTC-5.0 | O4b, plus 5 events of 6–8 April 2024, just before O4b | 161 | [Zenodo 20348005](https://zenodo.org/records/20348005) (part 1, with the skymaps) and [20348006](https://zenodo.org/records/20348006) (part 2) |
| `ALL` | all the catalogs above | O1 to O4b | | |

All the keys are **confident** catalogs: every event of their GWOSC lists has a probability of
astrophysical origin p_astro ≥ 0.5 (the re-analysed O1–O2 events of GWTC-2.1 carry no p_astro value).
The marginal candidates are only read by the `rates` and `hubble_constant` modes. Event counts of the
GWOSC lists in September 2026; the Zenodo records are the starting point, the version read being
resolved at run time ([Zenodo release versions](data-sources.md#zenodo-release-versions)); the run
dates are in [Data sources](data-sources.md#observing-runs).

Units: right ascension in degrees [0, 360), declination in degrees [−90, 90], masses in solar
masses (M☉), distances in megaparsecs (Mpc), probabilities in [0, 1].

## Examples

```bash
# catalog statistics of GWTC-4.0 from Zenodo, with sky-localization areas
python -m gwtc_analysis.cli catalog_statistics --catalogs GWTC-4 --include-area

# events with a primary mass above 50 Msun
python -m gwtc_analysis.cli event_selection --catalogs ALL --m1-min 50

# events whose 90% sky region contains a given position
python -m gwtc_analysis.cli search_skymaps --catalogs GWTC-4 --ra-deg 265.0 --dec-deg -46.0 --prob 0.9

# posteriors, strain overlay and q-transform of one event
python -m gwtc_analysis.cli parameters_estimation --src-name GW231223_032836

# GW170817, from public GWTC-1 products
python -m gwtc_analysis.cli parameters_estimation --src-name GW170817 \
    --overlay-start 0.2 --overlay-stop 0.2 --overlay-fmax 1000 --q-start 2 --q-stop 2 --q-fmax 1000

# merger rates, with the GWTC-5.0 injections downloaded automatically
python -m gwtc_analysis.cli rates

# Hubble constant: prepare the inputs, then 4 sampler runs, 2 at a time
python -m gwtc_analysis.cli hubble_constant --stages prepare
python -m gwtc_analysis.cli hubble_constant --stages sample combine report \
    --icarogw-python ~/.conda/envs/icarogw/bin/python --seeds 1 2 3 4 --parallel 2 --npool 2
```

## Outputs

Every mode writes:

- **TSV tables** (one row per event, per population or per parameter);
- **plots** (PNG, in a `--plots-dir` directory);
- an **HTML report** (`--out-report`), self-contained with its images embedded, suitable for Galaxy.

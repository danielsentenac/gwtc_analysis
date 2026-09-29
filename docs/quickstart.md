# Quick start

Every mode has its own help:

```bash
python -m gwtc_analysis.cli -h
python -m gwtc_analysis.cli <MODE> -h
```

Catalogs are named by case-sensitive keys:

| Key | Catalog |
|---|---|
| `GWTC-1` | Confident events of GWTC-1 (O1, O2) |
| `GWTC-2.1` | Confident events of GWTC-2.1 (O3a) |
| `GWTC-3` | Confident events of GWTC-3 (O3b) |
| `GWTC-4` | GWTC-4.0 (O4a) |
| `GWTC-5` | GWTC-5.0 (O4b) |
| `ALL` | All of the above |

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

# galaxy_catalog

The **galaxy catalog of a dark siren**: galaxies turned into the line-of-sight redshift prior that
[`hubble_constant --galaxy-catalog`](hubble-constant.md#dark-sirens-with-a-galaxy-catalog) uses, with the
pixelated catalog pipeline of [icarogw](https://github.com/icarogw-developers/icarogw)
[\[57\]](../references.md#ref-57) (the functions reviewed for the LVK analyses).

A dark siren adds the galaxies to a spectral siren: run the spectral siren of the same events first
([hubble_constant quick start](hubble-constant.md#quick-start)), then the two steps of the dark siren:

```
1. galaxy_catalog   galaxies → catalog_<band>_nside<N>_eps<ε>.hdf5      (once per catalog, any machine or cluster)
2. hubble_constant  --galaxy-catalog <that file>: prepare, sample, ...  (once per event set and mass model)
```

The catalog does not depend on the events or on the mass model: one catalog serves every `hubble_constant` work
directory. Its defaults are the GLADE+ K-band catalog of the GWTC-4.0 analysis [\[29\]](../references.md#ref-29),
which reproduces the published dark siren
([validation](hubble-constant.md#validation-gwtc-40-power-law-peak-glade-k-band)).

## Quick start

**GLADE+ K band on one machine** (about an hour on 4 CPUs, 4 GB of catalog file):

```bash
gwtc_analysis galaxy_catalog --workdir glade_k --jobs 4 --icarogw-python ~/.conda/envs/icarogw/bin/python
# then the dark siren
gwtc_analysis hubble_constant --workdir h0_dark_plp --stages prepare \
    --galaxy-catalog glade_k/catalog_K-glade+_nside64_eps1.hdf5
gwtc_analysis hubble_constant --workdir h0_dark_plp --stages sample combine reweight report \
    --icarogw-python ~/.conda/envs/icarogw/bin/python --seeds 1 2 3 4 --npool 4
```

**On a Slurm cluster**: download the galaxies, then the icarogw stages as batch jobs (about 45 min on CC-IN2P3):

```bash
gwtc_analysis galaxy_catalog --workdir /sps/.../glade_k --stages galaxies
gwtc_analysis galaxy_catalog --workdir /sps/.../glade_k --stages shard pixels gather prepare init interpolate finish summary \
    --executor slurm --jobs 32 --icarogw-python /sps/.../venv/bin/python \
    --slurm-option=--partition=htc --slurm-option=--licenses=sps --slurm-option=--mem=4G \
    --slurm-option=--cpus-per-task=1 --slurm-option=--time=06:00:00 --slurm-assembly-option=--mem=16G --submit
gwtc_analysis galaxy_catalog --workdir /sps/.../glade_k --stages report        # when the jobs have finished
```

**Another catalog** (DES, Rubin, …): give the file and map its columns ([details](#deeper-catalogs-and-clusters)):

```bash
gwtc_analysis galaxy_catalog --workdir des_r --input-catalog des_y6_gold/ --band r-upglade \
    --columns ra=RA dec=DEC z=DNF_Z sigmaz=DNF_ZSIGMA m=SOF_CM_MAG_CORRECTED_R \
    --where "EXT_XGB == 3 and SOF_CM_MAG_CORRECTED_R < 23.9" --nside 128 --nside-mthr 128 --zmin 0.05 --zcut 0.35 ...
```

## What the catalog is

For each event, the dark siren needs the redshift prior along every line of sight
[\[29\]](../references.md#ref-29) [\[26\]](../references.md#ref-26):

- the **in-catalog** part, the galaxies of the pixel brighter than its apparent-magnitude threshold, each a
  redshift probability (from its redshift and uncertainty) weighted by its luminosity to the power ε;
- the **out-of-catalog** part, the galaxies missed because they are fainter than the threshold, from the
  Schechter luminosity function of the band, uniform in comoving volume.

The likelihood evaluates both at every PE sample (sky pixel, redshift). This mode precomputes the in-catalog part
as an interpolant in redshift for every HEALPix pixel, and the threshold map that sets the out-of-catalog part.

## Options

Every option once, grouped by the question it answers. The column "GWTC-4.0" says where the default comes from:
**paper** = given in the paper, **open** = not given there (our choice, which can be changed).

### Which galaxies

| Option | Default | GWTC-4.0 | Meaning |
|---|---|---|---|
| `--source` | `glade-kband` | paper | GLADE+ galaxies with a Ks magnitude, downloaded from VizieR (VII/291) |
| `--glade-types` | `G` | paper ([checked](#glade-k-band)) | GLADE+ object types: `G` galaxies, `G,Q` with the quasars |
| `--glade-redshift` | `zcmb` | paper | `zcmb` (CMB frame, peculiar velocities corrected below z = 0.05) or `zhelio` (heliocentric) |
| `--glade-sigmaz` | `quadrature` | open | redshift error: measurement and peculiar-velocity errors in quadrature, `measurement` (`e_zhelio`) or `peculiar` (`e_z`) |
| `--sigmaz`, `--sigmaz-relative` | none | open | a constant redshift error instead (per 1 + z with `--sigmaz-relative`); also for a catalog without an error column |
| `--where` | none | — | a cut, as a pandas query on the columns: GLADE+ VizieR columns (`RAJ2000`, `DEJ2000`, `Kmag`, `zhelio`, `zcmb`, `f_zcmb`, `e_z`, `e_zhelio`), or those of `--input-catalog` |
| `--input-catalog` | none | — | another catalog: Parquet (file or HATS/LSDB tree), FITS, HDF5 or CSV, read in chunks |
| `--columns` | — | — | with `--input-catalog`: `ra=… dec=… z=… m=… [sigmaz=…]` |
| `--band` | `K-glade+` | paper | the icarogw band of the magnitude, which sets the Schechter function (K-glade+: M* = −23.39, α = −1.09, M from −27 to −19 [\[92\]](../references.md#ref-92)) |
| `--input-format`, `--angle-unit` | from the name, `deg` | — | format of `--input-catalog`, unit of its angles |
| `--galaxies` | none | — | a standard galaxy file already made (skips the `galaxies` stage) |

The galaxy file records its selection: changing it needs another work directory (the stage refuses to reuse a file
made with another selection).

### Completeness and weights

| Option | Default | GWTC-4.0 | Meaning |
|---|---|---|---|
| `--nside` | 64 | paper | HEALPix resolution of the catalog (0.84 deg² pixels) |
| `--nside-mthr` | 32 | paper | resolution of the apparent-magnitude threshold map (3.35 deg² pixels) |
| `--mthr-percentile` | 50 | paper | the threshold of a pixel: this percentile of its galaxies' magnitudes (the median) |
| `--epsilon` | 1 | paper | luminosity weight L^ε of the galaxies (0: every galaxy equally likely) |

### Redshift prior

| Option | Default | GWTC-4.0 | Meaning |
|---|---|---|---|
| `--ptype` | `gaussian` | open (the paper tests both Gaussian forms: "negligible differences") | each galaxy's redshift probability: `gaussian` = Gaussian likelihood × uniform-in-comoving-volume prior; `gaussian_nocom` = the Gaussian itself; `uniform` = uniform in volume within ±`--numsigma` σ |
| `--numsigma` | 3 | open | width of each galaxy's redshift probability, in σ |
| `--zmin`, `--zcut` | 0, 0.5 | paper (from 0), open (0.5) | redshift range of the in-catalog term. Galaxies below `--zmin` are left out; outside the range every galaxy counts as missed (completeness correction alone). `--zmin` needs a logarithmic grid, which then starts at it |
| `--nintegration` | `logspace:0.0001:5000` | open | the redshift grid: `logspace:ZMIN:N`, one logarithmic grid up to `--zcut` (5000 points: a step of 0.17%); an integer, icarogw's adaptive grid ([why not by default](#deeper-catalogs-and-clusters)) |

The published dark siren is compared in the `hubble_constant` report only when the catalog has all the settings
marked **paper**; the **open** ones do not prevent it.

### Execution

| Option | Default | Meaning |
|---|---|---|
| `--workdir` | `galaxy_catalog_run` | the work directory (restartable: finished stages and chunks are skipped) |
| `--stages` | all | stages to run ([below](#stages)) |
| `--jobs` | 4 | chunks of the parallel stages: processes on this machine, or array tasks on Slurm |
| `--icarogw-python` | the current interpreter | Python of the icarogw environment, which runs the icarogw stages |
| `--executor` | `local` | `local` or `slurm`: one batch script per stage, chained by dependencies |
| `--slurm-option` | none | an `#SBATCH` option of every job, repeatable |
| `--slurm-assembly-option` | none | an `#SBATCH` option of the jobs that hold the whole catalog (`init`, `finish`, `summary`), typically more memory |
| `--slurm-env-setup`, `--submit` | none, no | shell lines run first in each job; submit at once |
| `--nshards` | 1024 | files of the first pixelation pass (by pixel range) |
| `--out-report` | `galaxy_catalog.html` | the HTML report |
| `--settings` | none | a file of option values ([below](#settings-file)) |

Changing any catalog setting needs another work directory: the mode refuses to continue a catalog built with other
settings.

### Settings file

`--settings FILE` reads option values from a JSON file (or YAML, with PyYAML), the keys being the long option names;
options given on the command line override it. The resolved options of each invocation are written to
`<workdir>/options_galaxy_catalog.json`.

```json
{"glade-sigmaz": "measurement", "ptype": "gaussian_nocom", "jobs": 32, "executor": "slurm",
 "icarogw-python": "/sps/.../venv/bin/python",
 "slurm-option": ["--partition=htc", "--mem=4G", "--cpus-per-task=1", "--time=06:00:00", "--licenses=sps"]}
```

```bash
gwtc_analysis galaxy_catalog --workdir glade_k_nocom --settings glade_nocom.json --submit
```

All options with their help text: [CLI reference](../cli-reference.md#galaxy_catalog).

## GLADE+ K band

The galaxies of GLADE+ with a Ks magnitude (2MASS, Vega), downloaded from VizieR VII/291 in declination bands
(each cached, so an interrupted download resumes): 1,004,455 galaxies with a redshift and its uncertainty
(median z = 0.082, median σ_z = 0.015, median Ks = 13.5).

**The selection of the paper.** GLADE+ has 1,155,738 entries with a Ks magnitude, the "approximately 1.16 million
sources" of the paper; 133,158 of them have no redshift (mostly faint, median Ks = 14.1) and 18,106 are quasars
(median z = 0.76), which leaves the 1,004,455 galaxies. The threshold map shows the paper made the same choice:
with the galaxies that have a redshift, 4.3% of the nside-32 pixels are empty, as the paper's "~5%"; with every Ks
entry only 1.5% would be. The thresholds of 20, 40, 60 and 80% of the sky, 13.37, 13.47, 13.54 and 13.63, agree
with the labels of its Figure 3 (13.3, 13.5, 13.6, 13.7) to their rounding.

**Validated.** With this catalog, the GWTC-4.0 dark siren (Power Law + Peak, 137 BBH, 5000 PE samples per event)
gives H₀ = 114.6 (+41.5 / −33.8) km/s/Mpc against 115.4 (+40.1 / −33.8) in the GWTC-4.0 release
([hubble_constant](hubble-constant.md#validation-gwtc-40-power-law-peak-glade-k-band)). Built on CC-IN2P3 in about
45 min: 505,575 galaxies enter the in-catalog term, 94% of the sky has a threshold, median threshold Ks = 13.5;
luminosity-weighted completeness 0.71 at z = 0.04, 0.25 at 0.08, 0.04 at 0.13.

## Stages

Each stage is restartable: a finished stage or chunk leaves a marker in `<workdir>/done`.

| Stage | Does | Parallel |
|---|---|---|
| `galaxies` | the standard galaxy file (`ra`, `dec` in rad, `z`, `sigmaz`, `m`) | — |
| `shard` | one pass over the galaxies, spread by pixel range over `--nshards` files | — |
| `pixels` | one file per pixel (icarogw's layout), NaNs marked | `--jobs` chunks |
| `gather` | the list of filled pixels | — |
| `prepare` | threshold (percentile of the magnitudes in the coarse pixel) and redshift grid of each pixel | `--jobs` chunks |
| `init` | common redshift grid and threshold map, in the catalog file | — |
| `interpolate` | in-catalog line-of-sight interpolant of each pixel | `--jobs` chunks |
| `finish` | the catalog file read by the likelihood | — |
| `summary` | threshold map and completeness | — |
| `report` | HTML report | — |

`galaxies` and `report` run in the gwtc_analysis environment; the others in the icarogw environment
(`--icarogw-python`), through the self-contained `dark_catalog_icarogw.py`. Most of the time goes into `prepare`
and `interpolate`.

## Deeper catalogs and clusters

The pipeline is built so that catalogs much larger than GLADE+ (DES Y6 Gold: 350 million galaxies after cuts;
Rubin: billions) never need to fit in memory, and run on a cluster:

- **Input in chunks.** `--input-catalog` reads Parquet (a file, or a directory such as a HATS/LSDB partition tree),
  FITS, HDF5 or CSV in chunks, with the column mapping `--columns`, the cut `--where` applied to each chunk, and the
  icarogw band of the magnitude (`--band`); the column names of the example above are those of one DES release:
  use those of the catalog at hand.
- **Pixelation in shards.** `shard` streams the galaxies once; the per-pixel work is then split by pixel.
- **Batch jobs.** `--executor slurm` writes one sbatch script per stage in `<workdir>/slurm` (array jobs of `--jobs`
  tasks for `pixels`, `prepare` and `interpolate`) and `submit.sh`, which chains them with `afterok`
  dependencies; `--submit` runs it. The runner is copied next to the scripts, so the jobs only need the icarogw
  environment. On CC-IN2P3, every job must give its time limit, CPU count and memory, and jobs reading or writing
  `/sps` declare `--licenses=sps`.
- **Redshift grid.** The default is one logarithmic grid (`--nintegration logspace:ZMIN:N`), icarogw's fixed-grid
  path. icarogw's adaptive grid (an integer: points per galaxy) refines around every galaxy over ±`numsigma` σ: with
  the spectroscopic redshifts of GLADE+ (σ_z down to 1.5 × 10⁻⁴), the merged grid grows to tens of thousands of
  redshifts, and its serial merge (`init`) slowed to a few pixels per second on CC-IN2P3 before it was stopped. The
  grid must resolve the narrowest galaxy redshift probabilities: N points from ZMIN to `--zcut` give a relative step
  ln(zcut/ZMIN)/N (0.17% for the default).
- **Memory of the result.** The catalog file holds the interpolants on (redshifts × sky pixels); the likelihood
  loads them (in single precision). Its size is in the report; choose `--nside` and the redshift grid with it in
  mind.

The DES Y6 settings of the GWTC-5.0 analysis [\[30\]](../references.md#ref-30) [\[93\]](../references.md#ref-93)
are r band, nside 128 for the catalog and the threshold map, in-catalog part between z = 0.05 and 0.35
(`--zmin 0.05 --zcut 0.35`), and a Schechter function M* = −20.9, α = −1.01, M from −24.29 to −16.33. icarogw's
`r-upglade` band uses different Schechter parameters and expects per-galaxy K-corrections; a DES analysis needs its
own band definition, which is not in this version.

## Outputs

- `<workdir>/catalog_<band>_nside<N>_eps<ε>.hdf5`: the catalog for `hubble_constant --galaxy-catalog`;
- `<workdir>/summary.json`, `plots/mthr_map.png`, `plots/completeness.png` and `--out-report` (HTML);
- `<workdir>/catalog_settings.json`, `options_galaxy_catalog.json`: the settings of the catalog and of the run;
- `<workdir>/pixels/`: the per-pixel files (several GB for deep catalogs; they can be removed once the catalog
  file is finished).

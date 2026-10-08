# galaxy_catalog

The **galaxy catalog of the dark-siren analysis**: galaxies turned into the line-of-sight redshift prior that
[`hubble_constant --galaxy-catalog`](hubble-constant.md#dark-sirens-with-a-galaxy-catalog) uses, with the
pixelated catalog pipeline of [icarogw](https://github.com/icarogw-developers/icarogw)
[\[57\]](../references.md#ref-57) (the functions reviewed for the LVK analyses).

```bash
# GLADE+ K band, the GWTC-4.0 setup, 4 parallel chunks on this machine
gwtc_analysis galaxy_catalog --workdir glade_k --jobs 4 --icarogw-python ~/.conda/envs/icarogw/bin/python

# then the dark siren
gwtc_analysis hubble_constant --workdir h0_dark_plp --galaxy-catalog glade_k/catalog_K-glade+_nside64_eps1.hdf5 ...
```

## What the catalog is

For each event, the dark siren needs the redshift prior along every line of sight
[\[29\]](../references.md#ref-29) [\[26\]](../references.md#ref-26):

- the **in-catalog** part, the galaxies of the pixel brighter than its apparent-magnitude threshold, each a
  redshift likelihood (Gaussian, from its redshift and uncertainty) weighted by its luminosity to the power ε;
- the **out-of-catalog** part, the galaxies missed because they are fainter than the threshold, from the
  Schechter luminosity function of the band, uniform in comoving volume.

The likelihood evaluates both at every PE sample (pixel, redshift). This mode precomputes the in-catalog part as
an interpolant in redshift for every HEALPix pixel, and the threshold map that sets the out-of-catalog part.

## Defaults: the GWTC-4.0 analysis

| Setting | Option | Default | Source |
|---|---|---|---|
| galaxies | `--source glade-kband` | GLADE+ with a Ks magnitude | [\[91\]](../references.md#ref-91), [\[29\]](../references.md#ref-29) §3.2 |
| GLADE+ selection | `--glade-types`, `--glade-redshift`, `--glade-sigmaz`, `--sigmaz`, `--where` | galaxies (G), `zcmb`, measurement and peculiar-velocity errors in quadrature | entries: [\[29\]](../references.md#ref-29) §3.2 (checked below); errors: not given in the paper |
| band | `--band` | `K-glade+`: M* = −23.39, α = −1.09, M from −27 to −19 | [\[92\]](../references.md#ref-92), [\[29\]](../references.md#ref-29) §3.2 |
| catalog pixels | `--nside` | 64 | [\[29\]](../references.md#ref-29) §3.2 |
| threshold | `--nside-mthr`, `--mthr-percentile` | median magnitude in nside-32 pixels | [\[29\]](../references.md#ref-29) §2.2 |
| luminosity weight | `--epsilon` | 1 (ε = 0: every galaxy equally likely) | [\[29\]](../references.md#ref-29) §2.2 |
| galaxy redshift probability | `--ptype`, `--numsigma` | `gaussian`: Gaussian likelihood × uniform-in-comoving-volume prior, ±3σ | [\[29\]](../references.md#ref-29) §2.2 tests it and `gaussian_nocom` ("negligible differences") without saying which is the main result |
| redshift grid | `--nintegration`, `--zmin`, `--zcut` | logarithmic, 5000 points from z = 10⁻⁴ to 0.5 (step 0.17%); in-catalog part from z = 0 to 0.5 | icarogw's fixed-grid path (not given in the paper) |

**GLADE+ K band.** The galaxies of GLADE+ with a Ks magnitude (2MASS, Vega), downloaded from VizieR VII/291 in
declination bands (each cached, so an interrupted download resumes): 1,004,455 galaxies with a redshift and its
uncertainty (median z = 0.082, median σ_z = 0.015, median Ks = 13.5), as in the GWTC-5.0 analysis (9.9 × 10⁵
after cuts). The redshift is `zcmb`, corrected for peculiar velocities below z = 0.05; its uncertainty adds the
measurement error (2MPZ photometric for most) and the peculiar-velocity error in quadrature. Quasars (`Type` Q)
are left out.

**Validated.** With this catalog, the GWTC-4.0 dark siren (Power Law + Peak, 137 BBH) gives H₀ = 120.7
(+39.9 / −37.6) km/s/Mpc against 115.4 (+40.1 / −33.8) in the GWTC-4.0 release, 0.14σ apart
([hubble_constant](hubble-constant.md#validation-gwtc-40-power-law-peak-glade-k-band)). Built on CC-IN2P3 in about
45 min: 505,575 galaxies enter the in-catalog term, 94% of the sky has a threshold, median threshold Ks = 13.5;
luminosity-weighted completeness 0.71 at z = 0.04, 0.25 at 0.08, 0.04 at 0.13.

**The selection of the paper.** GLADE+ has 1,155,738 entries with a Ks magnitude, the "approximately 1.16 million
sources" of the paper; 133,158 of them have no redshift (mostly faint, median Ks = 14.1) and 18,106 are quasars
(median z = 0.76), which leaves the 1,004,455 galaxies above. The threshold map shows the paper made the same
choice: with the galaxies that have a redshift, 4.3% of the nside-32 pixels are empty, as the paper's "~5%";
with every Ks entry only 1.5% would be. The thresholds of 20, 40, 60 and 80% of the sky, 13.37, 13.47, 13.54 and
13.63, agree with the labels of its Figure 3 (13.3, 13.5, 13.6, 13.7) to their rounding.

## Settings the paper leaves open

The GWTC-4.0 paper gives the catalog, the band, the pixels, the threshold and the luminosity weight, but not the
redshift error of each galaxy, the form of its redshift probability, the redshift grid, nor the Monte Carlo settings
of the likelihood ([hubble_constant](hubble-constant.md#likelihood-thresholds) `--pe-samples`, `--neff-pe`,
`--neff-inj`). They can all be set:

| Option | Values |
|---|---|
| `--glade-types` | `G` (default), `G,Q` (with the quasars) |
| `--glade-redshift` | `zcmb` (default; CMB frame, peculiar velocities corrected below z = 0.05), `zhelio` (heliocentric) |
| `--glade-sigmaz` | `quadrature` (default), `measurement` (`e_zhelio`), `peculiar` (`e_z`); or a constant with `--sigmaz` (per 1 + z with `--sigmaz-relative`) |
| `--where` | a cut on the GLADE+ VizieR columns, e.g. `"f_zcmb == 1"` or `"Kmag < 13.5"` |
| `--ptype` | `gaussian` (default), `gaussian_nocom` (the Gaussian taken as the redshift posterior), `uniform` (uniform in comoving volume within ±`--numsigma` σ) |
| `--zmin`, `--zcut` | the redshift range of the in-catalog term: galaxies below `--zmin` are left out and, outside the range, icarogw counts every galaxy as missed (completeness correction alone). `--zmin` needs a logarithmic grid, which then starts at it |

A changed GLADE+ selection needs a new work directory (the galaxy file records its selection, and the stage refuses
an existing file made with another one); so does any other change of the catalog settings. The report of
`hubble_constant` compares with the published dark siren only when the catalog has the settings the paper gives
(band, ε, nside, threshold map, redshift range from 0, galaxy selection); `--ptype` and the grid, which the paper
leaves open, do not prevent it.

**Settings file.** `--settings FILE` reads the option values from a JSON (or YAML) file, keys as the long option
names; options on the command line override it. The resolved options of each run are written to
`<workdir>/options_galaxy_catalog.json`.

```json
{"glade-types": "G", "glade-sigmaz": "measurement", "ptype": "gaussian_nocom",
 "nintegration": "logspace:0.0001:5000", "zcut": 0.5}
```

```bash
gwtc_analysis galaxy_catalog --workdir glade_k_nocom --settings glade_nocom.json --jobs 4
```

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

The icarogw stages run in the icarogw environment (`--icarogw-python`), through the self-contained
`dark_catalog_icarogw.py`. GLADE+ takes about an hour on 4 CPUs, most of it in `prepare` and `interpolate`.

## Deeper catalogs and clusters

The pipeline is built so that catalogs much larger than GLADE+ (DES Y6 Gold: 350 million galaxies after cuts;
Rubin: billions) never need to fit in memory, and run on a cluster:

- **Input in chunks.** `--input-catalog` reads Parquet (a file, or a directory such as a HATS/LSDB partition tree),
  FITS, HDF5 or CSV in chunks, with a column mapping, a quality cut applied to each chunk, and the icarogw band of
  the magnitude:

  ```bash
  gwtc_analysis galaxy_catalog --stages galaxies --workdir des_r \
      --input-catalog des_y6_gold/ --columns ra=RA dec=DEC z=DNF_Z sigmaz=DNF_ZSIGMA m=SOF_CM_MAG_CORRECTED_R \
      --band r-upglade --where "EXT_XGB == 3 and SOF_CM_MAG_CORRECTED_R < 23.9"
  ```

  (the column names are an example: use those of the catalog at hand).
- **Pixelation in shards.** `shard` streams the galaxies once; the per-pixel work is then split by pixel.
- **Batch jobs.** `--executor slurm` writes one sbatch script per stage in `<workdir>/slurm` (array jobs of `--jobs`
  tasks for `pixels`, `prepare` and `interpolate`) and `submit.sh`, which chains them with `afterok`
  dependencies; `--submit` runs it. The runner is copied next to the scripts, so the jobs only need the icarogw
  environment:

  ```bash
  gwtc_analysis galaxy_catalog --workdir /sps/.../des_r --galaxies /sps/.../des_r/galaxies.h5 \
      --band r-upglade --nside 128 --nside-mthr 128 --nintegration logspace:0.001:2000 --zcut 0.35 \
      --executor slurm --jobs 200 --slurm-option=--partition=htc --slurm-option=--mem=4G --slurm-option=--cpus-per-task=1 \
      --slurm-option=--time=24:00:00 --slurm-option=--licenses=sps --slurm-assembly-option=--mem=32G \
      --slurm-env-setup "source ~/miniforge3/etc/profile.d/conda.sh; conda activate icarogw" \
      --icarogw-python ~/miniforge3/envs/icarogw/bin/python --submit
  ```

  `--slurm-assembly-option` adds options to the single jobs that hold the whole catalog in memory (`init`, `finish`,
  `summary`), typically more memory. On CC-IN2P3, every job must give its time limit, CPU count and memory, and jobs reading or writing `/sps` declare `--licenses=sps`.
  The `report` stage then runs anywhere with access to the work directory.
- **Redshift grid.** The default is one logarithmic grid (`--nintegration logspace:ZMIN:N`), icarogw's fixed-grid
  path. icarogw's adaptive grid (an integer: points per galaxy) refines around every galaxy over ±`numsigma` σ: with
  the spectroscopic redshifts of GLADE+ (σ_z down to 1.5 × 10⁻⁴), the merged grid grows to tens of thousands of
  redshifts, and its serial merge (`init`) slowed to a few pixels per second on CC-IN2P3 before it was stopped. The
  grid must resolve the narrowest galaxy redshift likelihoods: N points from ZMIN to `--zcut` give a relative step
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
- `<workdir>/pixels/`: the per-pixel files (several GB for deep catalogs; they can be removed once the catalog
  file is finished).

All options: [CLI reference](../cli-reference.md#galaxy_catalog).

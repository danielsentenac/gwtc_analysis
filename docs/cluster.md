# Running on a cluster

Most modes run in seconds to minutes on a laptop. Three modes have stages that take hundreds of CPU-hours:
`hubble_constant` (spectral and dark sirens), `galaxy_catalog` and `scalar_cloud`. For those stages, gwtc_analysis
writes the batch jobs itself and can submit them, for **Slurm**. This page gives:

- the Slurm parametrization of gwtc_analysis ([below](#slurm-parametrization));
- a complete setup on **CC-IN2P3**, the computing centre of Virgo, which runs Slurm ([below](#the-cc-in2p3-case-virgo));
- how to run on the **LIGO** clusters, which run HTCondor instead ([below](#ligo-clusters-htcondor)).

## When to use a cluster

Measured costs of the heavy stages (on CC-IN2P3):

| Task | Cost | Where |
|---|---|---|
| `hubble_constant --method spectral`, Power Law + Peak, 137 events | 8 runs of 4–5 h on 8 CPUs each (about 300 CPU-hours) | cluster |
| `hubble_constant --method spectral`, Multi Peak | 8 runs of 7–10 h on 8 CPUs each | cluster |
| `hubble_constant --method dark`, GLADE+ K band | runs of 8–10 h (Power Law + Peak) and 15–25 h (Multi Peak) on 8 CPUs each | cluster |
| reweighting of a spectral or dark siren to all the injections and 5000 PE samples | 16–32 chunks of 2–9 min | cluster, in the same chain |
| `galaxy_catalog`, GLADE+ K band | about 1 h on 4 CPUs on one machine; 45 min as 32 array jobs | either; a cluster for deeper catalogs (DES, Rubin) |
| `scalar_cloud`, GW190728 | about 2 CPU-hours per fit; 150 background injections of 1–2 CPU-hours; nested sampling of each cloud model about a day on 64 CPUs with 500 live points | cluster |
| `hubble_constant --method bright` and `--method joint`, the other modes | seconds to minutes | one machine |

On one machine, a spectral siren at 8 CPUs, one seed after another, takes about 1.5 days (Power Law + Peak) to
3 days (Multi Peak); a dark siren 3 to 8 days. The same commands then run with `--executor local` (the default):
`--parallel` runs several seeds at once, each with `--npool` processes.

## Slurm parametrization

### What gwtc_analysis does for you

With `--executor slurm`, a mode writes one sbatch script per step in its work directory (`<workdir>/slurm`, and for
`scalar_cloud` also `slurm_sample` and `slurm_seeded`), plus `submit.sh`, which submits them in order, each after
the previous one succeeded (`afterok`). Parallel steps are array jobs. With `--submit`, the mode runs `submit.sh`,
so one command queues the whole computation; without it, the scripts are only written (no Slurm is needed for
that). The runner script is copied next to the sbatch scripts, so the jobs only need the analysis Python
environment (icarogw for `hubble_constant` and `galaxy_catalog`), not gwtc_analysis itself. Each job sets
`HDF5_USE_FILE_LOCKING=FALSE` and puts the environment's `lib/` on `LD_LIBRARY_PATH`.

You still do these steps yourself:

- **before the jobs:** the `prepare` stage (and the `galaxies` stage of `galaxy_catalog`), which downloads the PE
  files or the galaxies. It runs anywhere the downloads are reachable, interactively or as a job of its own;
- **after the jobs:** the `report` stage, anywhere with access to the work directory;
- **follow-up:** monitoring the jobs (`squeue`, `sacct`, the logs in `<workdir>/logs`) and handling a run that does
  not finish: nested sampling occasionally enters a slow tail, and the chain can continue without that run
  ([how](modes/hubble-constant.md#on-a-slurm-cluster)).

### Options

| Option | Modes | What it sets |
|---|---|---|
| `--executor slurm` | the three | write sbatch scripts instead of running local processes |
| `--submit` | the three | submit the chain at once |
| `--slurm-option=OPT` | the three | an `#SBATCH` line of every job, repeatable: partition, memory, time limit, account, licences |
| `--slurm-assembly-option=OPT` | `galaxy_catalog` | an `#SBATCH` line of the jobs holding the whole catalog in memory (`init`, `finish`, `summary`), typically `--mem=16G` |
| `--slurm-env-setup "LINES"` | the three | shell lines run first in each job (module loads, conda activation) |
| `--icarogw-python PATH` | `hubble_constant`, `galaxy_catalog` | the Python of the icarogw environment, used by the jobs |
| `--python PATH` | `scalar_cloud` | the Python of the jobs (default: the one running gwtc_analysis) |
| `--seeds 1 2 ...`, `--npool N` | `hubble_constant` | one array task per seed, N CPUs each |
| `--reweight-jobs N`, `--reweight-pe-samples 5000` | `hubble_constant` | chunks of the reweighting (default 16); PE samples per event of its target |
| `--jobs N` | `galaxy_catalog` | array tasks of the chunked stages |
| `--npool N` | `scalar_cloud` | CPUs of each sampling job (default 32) |
| `--settings FILE` | `hubble_constant`, `galaxy_catalog` | option values from a JSON (or YAML) file, keys as the long option names; the command line overrides them |

The modes also set some `#SBATCH` options per job, written after those of `--slurm-option`, so that they win:

| Mode | Options set per job |
|---|---|
| `hubble_constant` | `--cpus-per-task`: `--npool` for the probe and the runs, 1 for the other jobs. Memory and time come from `--slurm-option` |
| `scalar_cloud` | CPUs, memory and time of every job: 1 CPU, 3–4 GB, 6 h to 1 day for the fits; `--npool` CPUs, max(8, 1.5 × `--npool`) GB and 4 days for the samplings |
| `galaxy_catalog` | none: give CPUs, memory and time with `--slurm-option` (and `--slurm-assembly-option`) |

### How each mode splits the work

**`hubble_constant`** (spectral and dark sirens), scripts in `<workdir>/slurm`
([details](modes/hubble-constant.md#on-a-slurm-cluster)):

| Job | Split | CPUs |
|---|---|---|
| `probe`, then `plan` | one job each, only with `--inj-fraction auto` and no runs yet: measures which injection subset is fast and accurate enough, and writes the choice to `plan.json` | `--npool`, then 1 |
| `run` | an array, **one task per seed**: independent nested-sampling runs | `--npool` each |
| `combine` | one job: the runs merged | 1 |
| `reweight`, `reweight_merge` | an array of `--reweight-jobs` chunks over the posterior samples, then one merge job: the posterior reweighted to all the injections and to `--reweight-pe-samples` PE samples per event | 1 each |

The runs sample with a subset of the injections and 1500 PE samples per event, which keeps each likelihood
evaluation fast; the reweighting corrects both subsets. For a comparison with a published result, add
`--reweight-pe-samples 5000`: with 1500 samples, the choice of the samples alone can move H₀ by up to about
8 km/s/Mpc ([comparison](science/spectral-siren.md#injection-subsets-and-reweighting)). This reweighting holds all
the PE samples in memory: 24 GB per chunk was enough for 137 events.

**`galaxy_catalog`**, scripts in `<workdir>/slurm` ([details](modes/galaxy-catalog.md#deeper-catalogs-and-clusters)):
the stages split by sky pixel. `pixels`, `prepare` and `interpolate` are arrays of `--jobs` chunks, each handling a
range of pixels. `shard`, `gather`, `init`, `finish` and `summary` are single jobs. Each finished stage or chunk
leaves a marker, so a resubmitted chain skips the work already done.

**`scalar_cloud`**, three chains ([details](modes/scalar-cloud.md)):

| Directory | Jobs |
|---|---|
| `slurm/` | the vacuum fit, an array of 28 profile points, the polish, an array of background segments |
| `slurm_sample/` | one nested-sampling task per model |
| `slurm_seeded/` | the seeds, then an array over the event and the injections |

`--stages sample` and `--stages seeded` are submitted separately; `seeded` starts from the superradiance sampling,
so it is submitted once that sampling has finished.

## The CC-IN2P3 case (Virgo)

[CC-IN2P3](https://doc.cc.in2p3.fr/) runs Slurm. Its rules for batch jobs:

- every job states its time limit, CPU count and memory (`--time`, `--cpus-per-task`, `--mem`);
- jobs reading or writing the `/sps` storage declare the licence `--licenses=sps`;
- the jobs run in the `htc` partition; the Slurm account (`virgo` for Virgo members) is set at submission.

The commands below are bash: if your login shell is tcsh, start `bash` first or put them in a bash script.

### Environment

Work directories, the Python environment and the PE cache go on `/sps`, which the jobs read and write
(`/sps/virgo/USERS/<user>` for Virgo members). One environment holds gwtc_analysis and icarogw (Python ≥ 3.12,
[installation](installation.md#icarogw-for-the-hubble_constant-mode)):

```bash
W=/sps/virgo/USERS/$USER/gwtc
python3.12 -m venv $W/venv
$W/venv/bin/pip install gwtc_analysis git+https://github.com/icarogw-developers/icarogw.git
export GWTC_PE_CACHE=$W/pe_cache          # the PE files, on /sps rather than in the home directory
```

### Settings files

The site options are best kept in a settings file per mode, passed with `--settings`. For `hubble_constant`:

```json
{"executor": "slurm", "icarogw-python": "/sps/virgo/USERS/<user>/gwtc/venv/bin/python",
 "seeds": [1, 2, 3, 4, 5, 6, 7, 8], "npool": 8, "reweight-pe-samples": 5000,
 "slurm-option": ["--partition=htc", "--licenses=sps", "--mem=24G", "--time=3-00:00:00"]}
```

`--cpus-per-task` is set by the mode from `--npool`. 16 GB is enough for the runs; 24 GB is for the reweighting
to 5000 PE samples, which uses the same option. For `galaxy_catalog`:

```json
{"executor": "slurm", "jobs": 32, "icarogw-python": "/sps/virgo/USERS/<user>/gwtc/venv/bin/python",
 "slurm-option": ["--partition=htc", "--licenses=sps", "--mem=4G", "--cpus-per-task=1", "--time=06:00:00"],
 "slurm-assembly-option": ["--mem=16G"]}
```

`scalar_cloud` has no settings file; its jobs set their own CPUs, memory and time, and need only the partition and
the licence: `--slurm-option=--partition=htc --slurm-option=--licenses=sps`.

### A spectral siren, step by step

```bash
W=/sps/virgo/USERS/$USER/gwtc; export GWTC_PE_CACHE=$W/pe_cache
# 1. inputs (PE samples, injections): downloads, once per work directory
$W/venv/bin/gwtc_analysis hubble_constant --method spectral --mass-model plp --workdir $W/h0_spectral_plp \
    --stages prepare --data-repo s3
# 2. the chain: probe, plan, 8 runs, combine, reweighting (all the injections, 5000 PE samples)
$W/venv/bin/gwtc_analysis hubble_constant --method spectral --mass-model plp --workdir $W/h0_spectral_plp \
    --stages sample combine reweight --settings h0_cc.json --submit
# 3. follow the jobs
squeue -u $USER
sacct -X -S today -o JobID%16,JobName%20,State%12,Elapsed
# 4. once the reweight_merge job is COMPLETED
$W/venv/bin/gwtc_analysis hubble_constant --method spectral --workdir $W/h0_spectral_plp --stages report
```

`--data-repo s3` downloads from the S3 mirror, faster than Zenodo ([data sources](data-sources.md)). The
`prepare` stage can also run as a job of its own (1 CPU and 16 GB were enough for 137 events), with the same
`--partition`, `--licenses` and `--time` lines. The dark siren is the same with `--method dark` and the catalog
built by `galaxy_catalog` (`--galaxy-catalog`); the joint posterior with GW170817 takes seconds and runs
interactively ([joint posterior](modes/hubble-constant.md#joint-posterior)).

## LIGO clusters (HTCondor)

The LIGO and IGWN clusters (the LDAS clusters at Caltech and the observatories, the IGWN computing grid) run
HTCondor, not Slurm ([IGWN computing guide](https://computing.docs.ligo.org/guide/htcondor/)), and gwtc_analysis
writes Slurm jobs only. Every IGWN job must give an `accounting_group` tag, or it is rejected
([accounting](https://computing.docs.ligo.org/guide/htcondor/accounting/)). The environment and the `prepare`
and `report` stages are the same as above. Two ways to run the heavy stages:

- **One HTCondor job running the local executor:** the same command with `--executor local`, `--seeds`,
  `--parallel` and `--npool`, in one job with `request_cpus` = `--parallel` × `--npool`. This is simple, but the whole
  computation must fit in one job on one machine.
- **One HTCondor job per task, chained by a DAG**, like the Slurm chain. Run the command once with
  `--executor slurm` and without `--submit`: this writes the scripts (no Slurm is needed) with the exact runner
  commands, which the HTCondor submit files then call. For `hubble_constant`, the runner is
  `<workdir>/slurm/h0_icarogw.py`. A fixed `--inj-fraction` replaces the `probe` and `plan` jobs (the probe chose 0.1
  for Power Law + Peak and Multi Peak, 0.5 for FullPop-4.0).

An example for a spectral siren, on a cluster with a shared file system:

```
# run.sub: one job per seed
universe         = vanilla
executable       = /path/to/icarogw-env/bin/python
arguments        = /path/to/h0_spectral_plp/slurm/h0_icarogw.py run --workdir /path/to/h0_spectral_plp --seed $(seed) --mass-model plp --prior-set gwtc4 --nlive 100 --npool 8 --naccept 60 --pe-samples 1500 --inj-fraction 0.1
environment      = "HDF5_USE_FILE_LOCKING=FALSE LD_LIBRARY_PATH=/path/to/icarogw-env/lib"
request_cpus     = 8
request_memory   = 16 GB
accounting_group = <your accounting tag>
log              = logs/h0.log
output           = logs/run_$(seed).out
error            = logs/run_$(seed).err
queue seed in (1 2 3 4 5 6 7 8)
```

`combine.sub` and `merge.sub` are the same with 1 CPU and the `combine` and `reweight-merge` commands. `reweight.sub`
runs `reweight --chunk $(chunk) --nchunks 16 --target-inj-fraction 1.0 --target-pe-samples 5000 --skip-if-done`, with
`request_memory = 24 GB` and `queue chunk from seq 0 15 |`. The DAG chains them:

```
# h0.dag, submitted with: condor_submit_dag h0.dag
JOB run run.sub
JOB combine combine.sub
JOB reweight reweight.sub
JOB merge merge.sub
PARENT run CHILD combine
PARENT combine CHILD reweight
PARENT reweight CHILD merge
```

Then run `--stages report` as on Slurm. The `galaxy_catalog` and `scalar_cloud` chains translate the same way: one
submit file per script in `<workdir>/slurm`, with an array becoming a `queue` over its chunks. On the IGWN grid
without a shared file system, the input and output files must also be transferred with the jobs; these examples
do not cover that.

# Running on a cluster

Most modes run in seconds to minutes on a laptop. Three modes have stages that take hundreds of CPU-hours. For
those stages, gwtc_analysis writes Slurm batch jobs itself and can submit them: `hubble_constant` (spectral and
dark sirens), `galaxy_catalog` and `scalar_cloud`. The other modes have no batch executor and run on one machine.

## When to use a cluster

Measured costs of the heavy stages (CC-IN2P3, `htc` partition):

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
3 days (Multi Peak); a dark siren 3 to 8 days. The same commands then run with `--executor local`
(the default): `--parallel` runs several seeds at once, each with `--npool` processes.

## What gwtc_analysis does for you

With `--executor slurm`, a mode writes one sbatch script per step in its work directory. Each step depends on the
previous one (`afterok`), and parallel steps are array jobs. With `--submit`, it also submits the chain, so a
single command queues the whole computation. The runner script is copied next to the sbatch scripts, so the
jobs only need the analysis Python environment (icarogw for `hubble_constant` and `galaxy_catalog`), not
gwtc_analysis itself.

You choose:

- **the size of the split:** the seeds (`--seeds`), the CPUs of each sampling job (`--npool`), and the number of
  chunks (`--jobs`, `--reweight-jobs`);
- **the cluster options**, passed through unchanged as `#SBATCH` lines (`--slurm-option`, repeatable): partition,
  memory, time limit, licences. On CC-IN2P3, every job must give its time limit, CPU count and memory, and jobs
  using `/sps` declare `--licenses=sps`;
- **the environment** of the jobs: `--icarogw-python` (or `--python` for `scalar_cloud`), and `--slurm-env-setup`
  for shell lines run first in each job (module loads, conda activation).

You still do these steps yourself:

- **before the jobs:** the `prepare` stage (and the `galaxies` stage of `galaxy_catalog`), where the PE files and
  the internet are available. Its output (`inputs.h5`, about 50 MB) can be copied to the cluster;
- **after the jobs:** the `report` stage, anywhere with access to the work directory;
- **follow-up:** monitoring the jobs (`squeue`, `sacct`, the logs in `<workdir>/logs`) and handling a run that does
  not finish. Nested sampling occasionally enters a slow tail; the chain can continue without that run
  ([how](modes/hubble-constant.md#on-a-slurm-cluster)).

## How each mode splits the work

**`hubble_constant`** (spectral and dark sirens), scripts in `<workdir>/slurm`
([details](modes/hubble-constant.md#on-a-slurm-cluster)):

| Job | Split | CPUs |
|---|---|---|
| `probe`, then `plan` | one job each, only with `--inj-fraction auto` and no runs yet: measures which injection subset is fast and accurate enough, and writes the choice to `plan.json` | `--npool`, then 1 |
| `run` | an array, **one task per seed**: independent nested-sampling runs | `--npool` each (`--cpus-per-task` set for you) |
| `combine` | one job: the runs merged | 1 |
| `reweight`, `reweight_merge` | an array of `--reweight-jobs` chunks (default 16) over the posterior samples, then one merge job: the posterior reweighted to all the injections and to `--reweight-pe-samples` PE samples per event | 1 each |

The runs sample with a subset of the injections and 1500 PE samples per event, which keeps each likelihood
evaluation fast. The reweighting then corrects both subsets. For a comparison with a published result, add
`--reweight-pe-samples 5000`: with 1500 samples, the choice of the samples alone can move H₀ by up to about
8 km/s/Mpc ([comparison](science/spectral-siren.md#injection-subsets-and-reweighting)). This reweighting holds all
the PE samples in memory: give its jobs more memory (24 GB per chunk was enough for 137 events).

```bash
gwtc_analysis hubble_constant --method spectral --workdir /sps/.../h0_spectral_plp --stages prepare
gwtc_analysis hubble_constant --method spectral --workdir /sps/.../h0_spectral_plp \
    --stages sample combine reweight --seeds 1 2 3 4 5 6 7 8 --npool 8 --reweight-pe-samples 5000 \
    --icarogw-python /sps/.../venv/bin/python --executor slurm \
    --slurm-option=--partition=htc --slurm-option=--mem=24G --slurm-option=--time=3-00:00:00 \
    --slurm-option=--licenses=sps --submit
gwtc_analysis hubble_constant --method spectral --workdir /sps/.../h0_spectral_plp --stages report   # afterwards
```

**`galaxy_catalog`**, scripts in `<workdir>/slurm` ([details](modes/galaxy-catalog.md#deeper-catalogs-and-clusters)):
the stages split by sky pixel. `pixels`, `prepare` and `interpolate` are arrays of `--jobs` chunks, each handling a
range of pixels. `shard`, `gather`, `init`, `finish` and `summary` are single jobs. The last three hold the whole
catalog in memory and take their own options (`--slurm-assembly-option`, for example `--mem=16G`). Each finished
stage or chunk leaves a marker, so a resubmitted chain skips the work already done.

**`scalar_cloud`**, three chains ([details](modes/scalar-cloud.md)):

| Directory | Jobs |
|---|---|
| `slurm/` | the vacuum fit, an array of 28 profile points, the polish, an array of background segments (1 CPU, 3 GB each) |
| `slurm_sample/` | one nested-sampling task per model, `--npool` CPUs each (default 32), at most 4 days |
| `slurm_seeded/` | the seeds, then an array over the event and the injections |

`--stages sample` and `--stages seeded` are submitted separately; `seeded` starts from the superradiance sampling,
so it is submitted once that sampling has finished.

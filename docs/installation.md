# Installation

`gwtc_analysis` runs on Linux and macOS with Python 3.10 or later. It is distributed as a **Python
package on conda-forge and on PyPI**, and as a **Docker image**:

| Distribution | Where | Install |
|---|---|---|
| Conda package `gwtc_analysis` | [conda-forge/gwtc_analysis](https://anaconda.org/conda-forge/gwtc_analysis) | `conda install -c conda-forge gwtc_analysis` |
| PyPI package `gwtc_analysis` | [pypi.org/project/gwtc_analysis](https://pypi.org/project/gwtc_analysis/) | `pip install gwtc_analysis` |
| Docker image `gwtc-tool` | [Docker Hub](https://hub.docker.com/r/danielsentenac/gwtc-tool/) | `docker pull danielsentenac/gwtc-tool` |
| Source | [GitHub](https://github.com/danielsentenac/gwtc_analysis) | `pip install -e .` |

The Docker image is the recommended choice for reproducibility and for workflow systems (CI
pipelines, computing clusters).

## From Docker

```bash
docker pull danielsentenac/gwtc-tool
docker run --rm -v "$PWD":/work -w /work danielsentenac/gwtc-tool gwtc_analysis -h
```

## From conda-forge

```bash
conda install -c conda-forge gwtc_analysis
gwtc_analysis -h
```

## From PyPI

```bash
pip install gwtc_analysis
gwtc_analysis -h
```

PyPI displays the project as `gwtc-analysis`, the spelling of its first registration; package names are
normalized, so `gwtc_analysis` and `gwtc-analysis` are the same project, and the package, its files, the
import name and the command are all `gwtc_analysis`.

The PyPI package declares only the light dependencies (numpy, pandas, matplotlib, minio, requests). The
modes that read PE files, strain and skymaps also need the gravitational-wave software stack (GWpy,
PESummary, PyCBC, ligo.skymap, LALSuite, h5py, astropy): install it first, for instance in an IGWN conda
environment, or use the conda-forge package or the Docker image.

## From source

The PE and strain modes need the gravitational-wave software stack (GWpy, PESummary, PyCBC,
ligo.skymap, LALSuite), which is easiest to get from the IGWN conda environments:

```bash
git clone https://github.com/danielsentenac/gwtc_analysis
cd gwtc_analysis
conda activate igwn          # or any environment with the GW stack
pip install -e .
gwtc_analysis -h
```

Installing the package (from conda-forge, PyPI or Docker) provides the `gwtc_analysis` command used
throughout this documentation. From a source checkout that is not installed, `python -m gwtc_analysis.cli`
is equivalent.

## icarogw, for the `hubble_constant` mode

The `sample` and `combine` stages of [`hubble_constant`](modes/hubble-constant.md) need
[icarogw](https://github.com/icarogw-developers/icarogw) [\[57\]](references.md#ref-57) and bilby [\[58\]](references.md#ref-58). icarogw requires Python ≥ 3.12
and is not on PyPI, so it usually lives in an environment of its own:

```bash
conda create -n icarogw python=3.12
conda activate icarogw
export TMPDIR=~/tmp                                                  # the torch wheels are large
pip install torch --index-url https://download.pytorch.org/whl/cpu   # CPU torch first, not the CUDA build
pip install git+https://github.com/icarogw-developers/icarogw.git
```

The mode is then pointed to that interpreter with `--icarogw-python ~/.conda/envs/icarogw/bin/python`.
It runs icarogw in CPU mode and puts the environment's `lib/` on `LD_LIBRARY_PATH` itself.

!!! tip "numpy version"
    If other packages in the icarogw environment need an older numpy (for instance ligo.skymap), pin
    it: `pip install numpy==2.1.1 scipy==1.14.1` worked.

## Caches

Downloads are cached so that each file is fetched once:

| Directory | Content |
|---|---|
| `~/.cache_gwtc_analysis/zenodo` | Zenodo version listings (one day) and sensitivity-injection files |
| `~/.cache_gwtc_analysis/pe_catalog` | PE samples extracted for `hubble_constant` (`$GWTC_PE_CACHE` or `--pe-cache` to move it) |
| `.cache_gwosc/` | Skymap tarballs and the PE index, per Zenodo record |
| `~/.gwcache` | Public GWTC-1 products for GW170817, supplementary PSDs |

# Installation

`gwtc_analysis` runs on Linux and macOS with Python 3.10 or later, and on several platforms.

| Distribution | Where |
|---|---|
| Docker image `gwtc-tool` | [Docker Hub](https://hub.docker.com/r/danielsentenac/gwtc-tool/) |
| Conda package | [anaconda.org/danielsentenac](https://anaconda.org/channels/danielsentenac/packages/gwtc_analysis/overview) |
| Galaxy tool | [usegalaxy.org](https://usegalaxy.org) |
| MMODA service | [MMODA LIGO-Virgo-KAGRA](https://www.astro.unige.ch/mmoda/) |
| Source | [GitHub](https://github.com/danielsentenac/gwtc_analysis) |

The Docker image is the recommended choice for reproducibility and for workflow systems (Galaxy,
CI pipelines).

## From source

The PE and strain modes need the gravitational-wave software stack (GWpy, PESummary, PyCBC,
ligo.skymap, LALSuite), which is easiest to get from the IGWN conda environments:

```bash
git clone https://github.com/danielsentenac/gwtc_analysis
cd gwtc_analysis
conda activate igwn          # or any environment with the GW stack
pip install -e .
python -m gwtc_analysis.cli -h
```

## icarogw, for the `hubble_constant` mode

The `sample` and `combine` stages of [`hubble_constant`](modes/hubble-constant.md) need
[icarogw](https://github.com/icarogw-developers/icarogw) and bilby. icarogw requires Python ≥ 3.12
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

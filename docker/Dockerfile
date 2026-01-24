FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# ---- system deps ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl bzip2 git \
    build-essential \
 && rm -rf /var/lib/apt/lists/*

# ---- install Miniforge (conda-forge) ----
ARG MINIFORGE_VERSION=24.7.1-0
RUN curl -L -o /tmp/miniforge.sh \
      "https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-Linux-x86_64.sh" \
 && bash /tmp/miniforge.sh -b -p /opt/conda \
 && rm -f /tmp/miniforge.sh

ENV PATH="/opt/conda/bin:$PATH"

# Keep conda quiet and fast-ish
RUN conda config --set channel_priority strict \
 && conda config --add channels conda-forge

# ---- create env with your runtime deps ----
# Adjust versions as you like; these are typical for your stack.
RUN conda create -y -n gwtc python=3.10 \
    numpy pandas matplotlib requests tqdm \
    minio \
    gwpy gwosc pesummary ligo.skymap astropy h5py \
 && conda clean -afy

# Always use env python without "conda activate"
ENV PATH="/opt/conda/envs/gwtc/bin:$PATH"

# Galaxy/Planemo friendliness: writable runtime dirs
ENV HOME=/tmp
ENV XDG_CACHE_HOME=/tmp/.cache
RUN mkdir -p /tmp/.cache

# ---- install your package into that env ----
WORKDIR /app
COPY . /app

# If your package metadata is correct, this will pull remaining python deps
RUN pip install --no-cache-dir -U pip setuptools wheel \
 && pip install --no-cache-dir .

# Smoke test (optional but recommended)
RUN python -c "import pandas, minio, gwpy, gwosc, pesummary, ligo.skymap; import gwtc_analysis; print('deps OK')"

# IMPORTANT: do not override entrypoint (Galaxy provides command)
ENTRYPOINT []
CMD ["python", "-m", "gwtc_analysis.cli", "--help"]


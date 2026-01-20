FROM mambaorg/micromamba:1.5.10

WORKDIR /opt/tool

# env spec lives in docker/
COPY docker/environment.yml /tmp/environment.yml

RUN micromamba create -y -n gwtc -f /tmp/environment.yml \
 && micromamba clean -a -y

# Copy the repository sources into the image
COPY . /opt/tool

# Install gwtc_analysis into the conda env
RUN micromamba run -n gwtc python -m pip install -e .

# Make the env binaries (python, pip, etc.) available without activation
ENV PATH=/opt/conda/envs/gwtc/bin:$PATH

ENTRYPOINT ["micromamba", "run", "-n", "gwtc"]

CMD ["python", "-c", "import gwtc_analysis; print('gwtc_analysis ok')"]


FROM mambaorg/micromamba:1.5.10

WORKDIR /opt/tool

COPY environment.yml /tmp/environment.yml

RUN micromamba create -y -n gwtc -f /tmp/environment.yml \
 && micromamba clean -a -y

# Make the env binaries (python, pip, etc.) available without activation
ENV PATH=/opt/conda/envs/gwtc/bin:$PATH

ENTRYPOINT ["micromamba", "run", "-n", "gwtc"]

CMD ["python", "-c", "import pesummary; print('pesummary ok')"]


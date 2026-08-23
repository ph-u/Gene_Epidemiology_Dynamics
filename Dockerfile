FROM mambaorg/micromamba
#FROM ghcr.io/mamba-org/micromamba:git-fddee42-cuda12.2.2-ubuntu20.04
RUN micromamba install --yes unzip && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root
COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && micromamba clean --all --yes
#RUN Rscript -e "BiocManager::install('xlsx2dfs')"
COPY raw/params.csv /tmp/params.csv
COPY src /src
CMD ["R", "--version"]
# podman build --tag gdyn . && podman run gdyn

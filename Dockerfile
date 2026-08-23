FROM mambaorg/micromamba
#FROM ghcr.io/mamba-org/micromamba:git-fddee42-cuda12.2.2-ubuntu20.04
RUN micromamba install --yes unzip && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root
COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH
#RUN Rscript -e "BiocManager::install('xlsx2dfs')"
COPY raw /raw
COPY src /src
RUN mkdir -p /data
RUN Rscript -e "devtools::install_github('GaelBn/BRREWABC')"
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
WORKDIR /src
CMD ["R", "--version"]
# podman build --tag gdyn . && podman run gdyn

FROM mambaorg/micromamba:2.0.5
#FROM ghcr.io/mamba-org/micromamba:git-fddee42-cuda12.2.2-ubuntu20.04
ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root
#RUN micromamba install --yes unzip && micromamba clean --all --yes
COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH
RUN Rscript -e "remotes::install_github('GaelBn/BRREWABC@<commit-sha', upgrade='never')"

ARG CACHEBUST
RUN echo "$CACHEBUST"

#RUN Rscript -e "BiocManager::install('xlsx2dfs')"
COPY raw /raw
COPY src /src
RUN mkdir -p /data
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ENV PATH=/src:${PATH}
WORKDIR /src
CMD ["R", "--version"]
# podman build --tag gdyn . && podman run gdyn

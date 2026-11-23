FROM mambaorg/micromamba:2.0.5
ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root

# ── System dependencies (via conda) ───────────────────────────────────────────
COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH

# ── Julia 1.11 ────────────────────────────────────────────────────────────────
# Installed directly from julialang.org to guarantee the pinned version.
# conda-forge's Julia package is typically one minor version behind.
ENV JULIA_VERSION=1.11.1
RUN curl -fsSL \
      "https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" \
    | tar xz -C /opt \
    && ln -s /opt/julia-${JULIA_VERSION}/bin/julia /usr/local/bin/julia \
    && julia --version

# ── Julia depot in a fixed, user-independent path ─────────────────────────────
# ROOT CAUSE OF THE SINGULARITY BUG:
# Docker builds run as root, so Pkg.instantiate() writes packages to /root/.julia/.
# Singularity runs containers as the HOST USER (e.g. ph17), so inside the container
# ~ resolves to that user's NFS home (/nfs/users/nfs_p/ph17) — not /root/.
# The host home has no Julia packages, so every `using` statement fails immediately.
#
# Fix: redirect the entire Julia depot to a fixed path that any user can read.
# JULIA_DEPOT_PATH is inherited by Singularity from the Docker ENV, so the same
# path is used at build time and at runtime regardless of who runs the container.
ENV JULIA_DEPOT_PATH=/opt/julia-depot

# ── Julia packages ────────────────────────────────────────────────────────────
# Instantiate from Project.toml so the exact versions from [compat] are locked.
# Precompiling here means zero compilation time at first run on the cluster.
COPY src/Project.toml /src/Project.toml
RUN mkdir -p /opt/julia-depot && \
    julia --project=/src -e ' \
    using Pkg; \
    Pkg.instantiate(); \
    Pkg.precompile(); \
    println("Julia packages installed OK — depot: " * ENV["JULIA_DEPOT_PATH"]) \
    '

# ── Repository contents ───────────────────────────────────────────────────────
COPY raw /raw
COPY src /src
# /data is created empty; mount real data at run time:
#   docker run -v /host/data:/data gdyn julia --project=/src --threads N /src/model.jl /data/input.csv 1
RUN mkdir -p /data

# Thread control: let the caller set JULIA_NUM_THREADS at runtime.
# OMP/BLAS flags prevent accidental over-subscription with BLAS threads.
ENV JULIA_NUM_THREADS=auto
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ENV PATH=/src:${PATH}
WORKDIR /src

# Smoke test
CMD ["julia", "--version"]
# Build:  podman build --tag gdyn .
# Run:    podman run --rm -it --threads 10 gdyn \
#           julia --project=/src --threads 10 /src/model.jl ../data/input.csv 1

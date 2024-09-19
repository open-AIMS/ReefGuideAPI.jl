# See https://hub.docker.com/_/julia for valid versions.
ARG JULIA_VERSION="1.10.4"

#------------------------------------------------------------------------------
# internal-base build target: julia with OS updates and an empty @reefguide
# Julia environment prepared for use. NOT intended for standalone use.
#------------------------------------------------------------------------------
FROM julia:${JULIA_VERSION}-bookworm AS internal-base

# Record the actual base image used from the FROM command as label in the compiled image
ARG BASE_IMAGE=$BASE_IMAGE
LABEL org.opencontainers.image.base.name=${BASE_IMAGE}

# Update all pre-installed OS packages (to get security updates)
# and add a few extra utilities
RUN --mount=target=/var/lib/apt/lists,type=cache,sharing=locked \
    --mount=target=/var/cache/apt,type=cache,sharing=locked \
    apt-get update \
    && apt-get -y upgrade \
    && apt-get install --no-install-recommends -y \
    git \
    less \
    nano \
    && apt-get clean \
    && apt-get autoremove --purge \
    && rm -rf /var/lib/apt/lists/*

# Tweak the JULIA_DEPOT_PATH setting so that our shared environments will end up
# in a user-agnostic location, not in ~/.julia => /root/.julia which is the default.
# See https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH
# This allows apps derived from this image to drop privileges and run as non-root
# user accounts, but still activate environments configured by this dockerfile.
ENV JULIA_DEPOT_PATH="/usr/local/share/julia"

# Prepare an empty @reefguide Julia environment for derived images to use - this is created in the shared depot path
RUN mkdir -p "${JULIA_DEPOT_PATH}" && \
    chmod 0755 "${JULIA_DEPOT_PATH}" && \
    julia -e 'using Pkg; Pkg.activate("reefguide", shared=true)'

# Ensure the @reefguide environment is in the load path for Julia, so that apps derived
# from this image can access any packages installed to there.
# (See https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_LOAD_PATH)
ENV JULIA_LOAD_PATH="@:@reefguide:@v#.#:@stdlib"

# Run Julia commands by default as the container launches.
# Derived applications should override the command.
ENTRYPOINT ["julia", "--project=@reefguide"]


#------------------------------------------------------------------------------
# reefguide-base build target: ReefGuideApi.jl preinstalled as a non-dev package.
# Use `--target=reefguide-base` on your `docker build command to build *just* this.
#------------------------------------------------------------------------------
FROM internal-base AS reefguide-base

# What version of ReefGuideApi.jl from package registry to install in reefguide-base
# TODO this doesn't make sense yet
ARG REEFGUIDE_VERSION="0.11.0"

# Which julia package registry version to install
ENV REEFGUIDE_VERSION=$REEFGUIDE_VERSION

# Try to coerce Julia to build across multiple targets
ENV JULIA_CPU_TARGET=x86_64;haswell;skylake;skylake-avx512;tigerlake

# Install ReefGuideApi.jl into the @reefguide shared environment as an unregistered package.
# - Allow the package source and version to be overridden at build-time.
# - Include citation information for ReefGuideApi.jl in the image labels.
RUN mkdir -p "${JULIA_DEPOT_PATH}" && \
    chmod 0755 "${JULIA_DEPOT_PATH}" && \
    julia --project=@reefguide -e "using Pkg; Pkg.add(name=\"ReefGuideApi\", version=\"${REEFGUIDE_VERSION}\"); Pkg.instantiate(); Pkg.precompile(); using ReefGuideApi;"
LABEL au.gov.aims.reefguideapi.source="https://github.com/open-AIMS/ReefGuideApi.jl/releases/tag/v${REEFGUIDE_VERSION}" \
    au.gov.aims.reefguideapi.version="${REEFGUIDE_VERSION}" \
    au.gov.aims.reefguideapi.vendor="Australian Institute of Marine Science" \
    au.gov.aims.reefguideapi.licenses=MIT

#------------------------------------------------------------------------------
# reefguide-dev build target: installs directly from source files in this repo.
#------------------------------------------------------------------------------
FROM internal-base AS reefguide-dev

ENV REEFGUIDE_ENV_DIR="${JULIA_DEPOT_PATH}/environments/reefguide" \
    REEFGUIDE_SRC_DIR="/usr/local/src/reefguide"

# Install the versioned .toml file(s) into the shared reefguide environment and use
# those to set up the ReefGuideApi source code as a development package in the
# shared @reefguide environment, pre-installing and precompiling dependencies.
# This should *hugely* speeds up the build when other ReefGuideApi files have changed,
# as this very slow step won't need to be run again.

WORKDIR "${REEFGUIDE_SRC_DIR}"
COPY ./Project.toml ./Project.toml
# TODO include Manifest.toml
COPY ./Manifest.toml ./Manifest.toml
RUN julia --project=@reefguide -e 'using Pkg;  Pkg.instantiate(verbose=true)'

# Install the ReefGuideApi source code and configure it as a development
# package in the @reefguide shared environment.
# Should be v speedy if the .toml file is unchanged, because all the
# dependencies *should* already be installed.
COPY . .
RUN julia --project=@reefguide \
    -e  'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.precompile(); using ReefGuideApi;'
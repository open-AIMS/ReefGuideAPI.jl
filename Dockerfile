# See https://hub.docker.com/_/julia for valid versions.
ARG JULIA_VERSION="1.11.5"

#------------------------------------------------------------------------------
# internal-base build target: julia with OS updates and an empty @reefguide
# Julia environment prepared for use. NOT intended for standalone use.
#------------------------------------------------------------------------------
FROM julia:${JULIA_VERSION}-bookworm AS internal-base

# Record the actual base image used from the FROM command as label in the compiled image
ARG BASE_IMAGE="julia:${JULIA_VERSION}-bookworm"
LABEL org.opencontainers.image.base.name=${BASE_IMAGE}


# Update all pre-installed OS packages (to get security updates)
# and add a few extra utilities -
# Installs build essentials for sysimage
RUN --mount=target=/var/lib/apt/lists,type=cache,sharing=locked \
    --mount=target=/var/cache/apt,type=cache,sharing=locked \
    apt-get update \
    && apt-get -y upgrade \
    && apt-get install --no-install-recommends -y \
    build-essential \
    git \
    less \
    nano \
    gdal-bin \
    libgdal-dev \
    libfftw3-dev \
    && apt-get clean \
    && apt-get autoremove --purge \
    && rm -rf /var/lib/apt/lists/*


# Tweak the JULIA_DEPOT_PATH setting so that our shared environments will end up
# in a user-agnostic location, not in ~/.julia => /root/.julia which is the default.
# See https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH
# This allows apps derived from this image to drop privileges and run as non-root
# user accounts, but still activate environments configured by this dockerfile.
ENV JULIA_DEPOT_PATH="/usr/local/share/julia"
ENV JULIA_PKG_USE_CLI_GIT=true

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
# reefguide-base build target: ReefGuideAPI.jl preinstalled as a non-dev package.
# Use `--target=reefguide-base` on your `docker build command to build *just* this.
#------------------------------------------------------------------------------

# TODO enable and update this once the package is available from registry
#    FROM internal-base AS reefguide-base
#
#    # What version of ReefGuideAPI.jl from package registry to install in reefguide-base
#    # TODO this doesn't make sense yet
#    ARG REEFGUIDE_VERSION="0.1.0"
#
#    # Which julia package registry version to install
#    ENV REEFGUIDE_VERSION=$REEFGUIDE_VERSION
#
#    # Try to coerce Julia to build across multiple targets
#    ENV JULIA_CPU_TARGET=x86_64;haswell;skylake;skylake-avx512;tigerlake
#
#    # Install ReefGuideAPI.jl into the @reefguide shared environment as an unregistered package.
#    # - Allow the package source and version to be overridden at build-time.
#    # - Include citation information for ReefGuideAPI.jl in the image labels.
#    RUN mkdir -p "${JULIA_DEPOT_PATH}" && \
#        chmod 0755 "${JULIA_DEPOT_PATH}" && \
#        julia --project=@reefguide -e "using Pkg; Pkg.add(name=\"ReefGuideAPI\", version=\"${REEFGUIDE_VERSION}\"); Pkg.instantiate(); Pkg.precompile(); using ReefGuideAPI;"
#    LABEL au.gov.aims.reefguideapi.source="https://github.com/open-AIMS/ReefGuideAPI.jl/releases/tag/v${REEFGUIDE_VERSION}" \
#        au.gov.aims.reefguideapi.version="${REEFGUIDE_VERSION}" \
#        au.gov.aims.reefguideapi.vendor="Australian Institute of Marine Science" \
#        au.gov.aims.reefguideapi.licenses=MIT
#
#    # Expect to include the prepped data at /data/reefguide and the config at
#    # /data/.config.toml
#    VOLUME ["/data/reefguide"]
#
#    EXPOSE 8000
#
#    # Run Julia commands by default as the container launches.
#    # Derived applications should override the command.
#    ENTRYPOINT ["julia", "--project=@reefguide", "-t", "1,auto", "-e", "using ReefGuideAPI; ReefGuideAPI.start_server(\"/data/reefguide/config.toml\")"]

#------------------------------------------------------------------------------
# reefguide-src build target: installs directly from source files in this repo.
#------------------------------------------------------------------------------
FROM internal-base AS reefguide-src

ENV REEFGUIDE_ENV_DIR="${JULIA_DEPOT_PATH}/environments/reefguide" \
    REEFGUIDE_SRC_DIR="/usr/local/src/reefguide" \
    JULIA_PKG_USE_CLI_GIT=true

# Try to coerce Julia to build across multiple targets
ENV JULIA_CPU_TARGET=x86_64;haswell;skylake;skylake-avx512;tigerlake

# Install the versioned .toml file(s) into the shared reefguide environment and use
# those to set up the ReefGuideAPI source code as a development package in the
# shared @reefguide environment, pre-installing and precompiling dependencies.
WORKDIR "${REEFGUIDE_SRC_DIR}"

# Copy project and manifest - includes Manifest-v1.11 etc
COPY Project.toml Manifest*.toml ./

# Then fire up a julia execution just to dump out the version
RUN echo $(julia --project=.  -e 'using Pkg; println(Pkg.dependencies()[Base.UUID("856f044c-d86e-5d09-b602-aeab76dc8ba7")].version)') | cut -d '+' -f 1 >> mkl.dep

# Compile MKL_jll first - this improves build time significantly - unsure exactly why
RUN MKL_VERSION=$(cat mkl.dep) julia -e 'using Pkg; Pkg.add(PackageSpec(name="MKL_jll", version=ENV["MKL_VERSION"])); Pkg.precompile()'

# Install package compiler
RUN julia --project=@reefguide \
    -e 'using Pkg; Pkg.add("PackageCompiler")'

# Add custom fork of GeometryOps
# TODO: Remove once fix that resolves precompilation issues gets released
#       (blocks system image generation)
RUN julia --project=@reefguide \
    -e 'using Pkg; Pkg.add(url="https://github.com/ConnectedSystems/GeometryOps.jl", rev="main");'

# Install the ReefGuideAPI source code and configure it as a development
# package in the @reefguide shared environment.
# Should be v speedy if the .toml file is unchanged, because all the
# dependencies *should* already be installed.
COPY ./src src

# Dev and precompile reefguide API
RUN julia --project=@reefguide \
    -e  'using Pkg; Pkg.develop(PackageSpec(path=pwd())); using ReefGuideAPI; Pkg.precompile();'

# Build custom sys image
RUN julia --project=@reefguide \
    -e  'using ReefGuideAPI; using PackageCompiler; create_sysimage(["ReefGuideAPI"]; sysimage_path="ReefGuideSysImage.so", cpu_target="x86_64;haswell;skylake;skylake-avx512;tigerlake")'

# Expect to include the prepped data at /data/reefguide and the config at
# /data/.config.toml
VOLUME ["/data/reefguide"]

# Server port
EXPOSE 8000

# By default, drops the user into a  julia shell with ReefGuideAPI activated (from custom sys image)
ENTRYPOINT ["julia", "--project=@reefguide", "-t", "auto,1", "-J", "ReefGuideSysImage.so", "-e"]

# Derived applications should override the command e.g. to start worker use
CMD ["using ReefGuideAPI; ReefGuideAPI.start_worker()"]

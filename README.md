# Reef Guidance API

API for supporting reef suitability assessments.

[![Documentation](https://img.shields.io/badge/docs-dev-blue)](https://open-aims.github.io/ReefGuideAPI.jl/dev/)

## Table of Contents

- [Quick Start](#quick-start)
- [Setup](#setup)
- [Configuration](#configuration)
- [JWT Authentication](#jwt-auth-configuration)
- [API Usage](#api-usage)
- [Development](#development-setup)
- [Performance Notes](#performance-notes)
- [Reef Edge Alignment](#reef-edge-alignment-for-site-searching)
- [Docker Usage](#docker-build-and-run)
  - [Dockerfile Configuration](#dockerfile-configuration)
  - [Publishing Releases](#publish-release)
  - [Data Requirements](#mounting-files-and-required-data)
  - [Building and Running](#building-and-running)
- [Troubleshooting](#troubleshooting)

## Quick Start

Start the API server with a configuration file:

```julia
using ReefGuideAPI

# To enable debug messages:
# ENV["JULIA_DEBUG"] = "ReefGuideAPI"

# If multiple threads are available, a parallel server will be spun up
ReefGuideAPI.start_server(".config.toml")
```

### Using Docker

```bash
# Run with Docker
docker run -p 8000:8000 -v ./data:/data/reefguide reefguide

# Or with Docker Compose
docker compose up reefguide-src
```

## Setup

Initialize the project the usual way:

```julia
]instantiate
```

A TOML file should be defined indicating location of the MPA dataset.
These are currently the files/data created in Step/Script 1a in https://github.com/open-AIMS/GBR-reef-guidance-assessment

```toml
[prepped_data]
PREPPED_DATA_DIR = "C:/some_path_to_data/MPA/"

[server_config]
TIFF_CACHE_DIR = "<some location to cache geotiffs>"
REGIONAL_CACHE_DIR = "<some location to cache regional datasets>"
DEBUG_MODE = "false"  # Optional, disables file caching and displays debug logs
COG_THREADS = "2"  # Optional, Number of threads to use when creating COGs (defaults to 1)
TILE_SIZE = "256"  # Optional, tile block size to use (defaults to 256)
```

By convention, this file is named `.config.toml` (note the leading `.`).

## JWT Auth Configuration

The API can be additionally configured to expect a valid JWT in the `Authorisation: Bearer <token>` header format.

Add the following to `.config.toml`:

```toml
[jwt_auth]
# Enable JWT auth : bool true/false
JWT_ENABLED = true
# Which iss to validate for the JWTs?
JWT_ISS = "https://issuer.com"
# WKT JWKS endpoint where public key can be retrieved
WKT_ENDPOINT = "https://issuer.com/api/.well-known/jwks.json"
```

Pay attention to the issuer and wkt endpoints. The first should exactly match the expected JWT issuer claim. The second should be web-resolvable and return a WKT JSON which provides the public key.

### Authentication Notes

- When JWT authentication is enabled, most API routes require a valid token
- The health check endpoint (`/health`) remains accessible without authentication
- Tokens must be included in the Authorisation header format: `Authorization: Bearer <token>`
- For testing purposes, you can disable authentication by setting `JWT_ENABLED = false`

## Publishing a new version

1. Bump the version in `Project.toml`
1. Identify the previous version
2. Identify the intended version bump (new tag)
3. Create tag

```bash
git tag v1.x.y -a
# Then fill in description
```

4. Push to remote

```bash
git push origin --tags
```

5. On GitHub - draft the new release targeting this as the new version

## API Usage

### Dynamic COG Generation

Example URL:

```
http://127.0.0.1:8000/assess/Cairns-Cooktown/slopes?Depth=-9.0:0.0&Slope=0.0:40.0
```

### Simple Slippy Tiles

Example URL:

```
http://127.0.0.1:8000/tile/8/231/139?region=Cairns-Cooktown&rtype=slopes&Depth=-9.0:0.0&Slope=0.0:40.0&Rugosity=0.0:3.0
```

## Development Setup

The steps below assumes you are in the project root.

Create a sandbox environment:

```bash
$ mkdir sandbox
$ cd sandbox
$ julia --project=.
(sandbox) julia> ]add Revise Infiltrator
(sandbox) julia> ]dev ..
```

Copy the quickstart to a file (e.g., `dev_server.jl`).

Create the `.config.toml` file and save to the sandbox directory.

Assuming VS Code is configured to default to the sandbox environment and start the
Julia REPL at project root:

```julia
;cd sandbox
include("dev_server.jl")
```

Note that the server now caches the initially loaded spatial data in between server
launches to reduce downtime. It will be necessary to restart the Julia session to reload
spatial data.

## Performance Notes

The config setting `COG_THREADS` controls how many threads should be requested when writing
out COGs. Ideally this will be set to at least 2 (preferably 4).
Higher values do seem to reduce write times but with diminishing returns (tested up to 8).
Locally, write times with four threads configured range from 10 to 15 seconds.

## Reef Edge Alignment for Site Searching

`identify_edge_aligned_sites()` can be used to identify potential sites that only align with
the nearest reef edge (or specified rotations away from this angle).
This method works by identifying the closest edge of reef polygon geometries that have been
converted into lines.

The following processing is required before use:

- Reef polygons should be simplified (`GO.simplify()`) and buffered to avoid matching possibly inaccurate reef edges.
- Simplified reef polygons should be provided as vertex-vertex lines with `polygon_to_lines()`.
- Require raster of target pixels to search, and their indices (currently a vector of `CartesianIndices` for identifying search pixels). Use `findall(bool_search_raster)` to return pixel indices.
- Raster of search pixels should be masked by reef polygons or simplified reef polygons.
- The target region name should be specified in GBRMPA format.
  - E.g. "Townsville/Whitsunday Management Area" rather than "Townsville-Whitsunday".

### Parquet Assessment Additional Setup

- A parquet GeoDataFrame must be loaded and filtered for unsuitable pixels based on user criteria thresholds using a Dict and `within_thresholds()`.
- `lons` and `lats` columns (FLoat64) must be added to the GeoDataFrame.
  - E.g. `valid_pixels.lons = first.(GI.coordinates.(valid_pixels.geometry))`
    The column used for masking should be the same as the column specified as geometry_col in
    `identify_edge_aligned_sites` (default = `:geometry`).

## Docker Build and Run

The ReefGuideAPI.jl package has an associated `Dockerfile` and build/publish process. This means you can run an instance of the ReefGuideAPI.jl package without needing to compile/build it with a local `Julia` installation. You will be able to view the latest published versions of the Docker image on the repository packages page.

### Dockerfile Configuration

The Julia version is specified by the `JULIA_VERSION` arg. The full version is specified
to maintain build stability, but should be bumped to the latest version of Julia when a
release is published.

#### A Note About MKL_jll

Due to how Julia (particularly v1.11) handles precompilation, it significantly reduces the build time by explicitly installing MKL_jll before installing any explicit project dependencies.

For this reason, the Dockerfile extracts the MKL_jll version from the Manifest file using Pkg.dependency(), precompiles this in an anonymous project, then compiles the main dependencies. This cuts the build time from around 15 minutes down to around 6-7.

### Publish Release

To publish a new version of the Docker image:

1. Bump version in Project.toml and [PublishDockerImage.yml](../../.github/workflows/PublishDockerImage.yml)
2. Create PR, merge to `main` branch
3. Publish Release on GitHub (this triggers the PublishDockerImage workflow)

### Mounting Files and Required Data

The `ReefGuideAPI.jl` package requires:

- a `config.toml` file and
- a set of input data files

Please include these in a folder called `data` in your working directory.

When running the commands below, it is assumed you have `data` available locally with the required files.

**Note**: Due to how Docker excludes `.` files, we have named the config file `config.toml` in the data folder (not `.config.toml`). This is required to launch the server.

### Building and Running

#### Build from Source Files

Using Docker:

```bash
docker build . --target reefguide-src -t reefguide
```

Using Docker Compose:

```bash
docker compose build reefguide-src
# Or to build and run in one command:
docker compose up --build reefguide-src
```

#### Run Server with Mounted Files

Using Docker:

```bash
docker run -p 8000:8000 -v ./data:/data/reefguide reefguide
```

Using Docker Compose:

```bash
docker compose up reefguide-src
```

#### Run Interactive Shell with Mounted Files

This will start a Julia shell where `ReefGuideAPI` is compiled and ready for use:

```bash
docker run --rm --interactive --entrypoint="julia" --tty -v ./data:/data/reefguide reefguide
```

Then in the Julia REPL:

```julia
using ReefGuideAPI
ReefGuideAPI.start_server("/data/reefguide/config.toml")
```

## Troubleshooting

### Common Issues

#### Configuration Problems

**Issue**: Server fails to start with configuration error

- **Solution**: Double-check that your config file has the correct path format for your OS (Windows uses backslashes, Unix uses forward slashes)
- **Solution**: Ensure all required directories specified in the config actually exist on your system

**Issue**: Data not found

- **Solution**: Verify the `PREPPED_DATA_DIR` path is correct and contains the required MPA dataset files
- **Solution**: When using Docker, check that the volume mount path (`-v ./data:/data/reefguide`) is correct and the `data` directory contains your files

#### Docker-Specific Issues

**Issue**: Docker container exits immediately after starting

- **Solution**: Check Docker logs with `docker logs <container-id>` to see the specific error
- **Solution**: Verify that your `config.toml` is named correctly (not `.config.toml`) in the data directory

**Issue**: "Permission denied" errors when accessing data directory

- **Solution**: Check file permissions on your data directory and ensure the Docker user has read/write access

#### Performance Issues

**Issue**: Slow COG generation

- **Solution**: Increase the `COG_THREADS` value in your config (recommend 2-4)
- **Solution**: Ensure the server has sufficient memory for processing (minimum 4GB recommended)

#### Authentication Issues

**Issue**: JWT token rejected as invalid

- **Solution**: Verify the token's issuer (`iss` claim) matches exactly what's in your config
- **Solution**: Check that the token is not expired
- **Solution**: Ensure the `WKT_ENDPOINT` is accessible from the server

### Getting Help

If you encounter issues not covered here, please:

1. Check the Julia REPL output for specific error messages
2. Enable debug logs with `ENV["JULIA_DEBUG"] = "ReefGuideAPI"`
3. Open an issue on the GitHub repository with a detailed description of the problem

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

1. Bump the version in `Project.toml` and commit
2. Identify the previous version
3. Identify the intended version bump (new tag)
4. Create tag

```bash
git tag v1.x.y -a
# Then fill in description
```

5. Push to remote

```bash
git push origin --tags
```

6. On GitHub - draft the new release targeting this as the new version

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

## Job worker component

The ReefGuide Job System is a job processing framework written in Julia. The system is designed to poll for available jobs from the reefguide-web-api, claim them, process them using appropriate handlers, and report results back to the API.

### System Architecture

The job system is built around these core components:

- **WorkerService**: The main orchestrator that manages the job processing lifecycle
- **JobHandler**: A wrapper around the job registry which manages dispatching jobs to the proper registry handler, and reporting the results
- **API Client**: For communication with the ReefGuide API, automatically managing logins/token refreshing
- **Job Registry**: A centralized registry that maps job types to their handlers, including input and output payload configurations

### Configuration

#### Environment Variables

The worker requires the following environment variables to be set:

| Variable                        | Description                                                | Example                                   |
| ------------------------------- | ---------------------------------------------------------- | ----------------------------------------- |
| `API_ENDPOINT`                  | Base URL of the ReefGuide API                              | `"https://api.reefguide.example.com/api"` |
| `JOB_TYPES`                     | Comma-separated list of job types the worker should handle | `"CRITERIA_POLYGONS"`                     |
| `USERNAME`                      | API authentication username                                | `"worker-service"`                        |
| `PASSWORD`                      | API authentication password                                | `"secure-password"`                       |
| `POLL_INTERVAL_MS`              | (Optional) Polling interval in milliseconds                | `2000` (default)                          |
| `IDLE_TIMEOUT_MS`               | (Optional) Idle timeout in milliseconds                    | `300000` (default)                        |
| `AWS_REGION`                    | AWS region for S3 operations                               | `"ap-southeast-2"`                        |
| `ECS_CONTAINER_METADATA_URI_V4` | ECS task metadata endpoint                                 | Automatically set in ECS                  |

#### Using DotEnv for Local Development

For local development, you can use the DotEnv.jl package to load environment variables from a `.env` file:

```julia
using DotEnv

# Load environment variables from .env file
DotEnv.load!()

# Create and start worker
using ReefGuideAPI
ReefGuideAPI.start_worker()
```

Example `.env` file:

```
API_ENDPOINT=http://localhost:8000
JOB_TYPES=CRITERIA_POLYGONS
USERNAME=local-dev
PASSWORD=local-password
POLL_INTERVAL_MS=5000
IDLE_TIMEOUT_MS=600000
AWS_REGION=ap-southeast-2
```

### Worker Operation

#### Polling Loop

The worker operates in a continuous loop:

1. **Poll** for available jobs matching the configured job types
2. If a job is found, **claim** it to get an assignment
3. **Process** the job with the appropriate handler
4. **Report** the result (success/failure + payload)
5. **Sleep** for the configured poll interval
6. Check for **idle timeout** and shut down if idle too long

The loop includes error handling to ensure that failures in one job don't affect the processing of subsequent jobs.

#### Idle Timeout

The worker will automatically shut down after being idle (no jobs processed) for the configured timeout period. This helps manage resources in cloud environments where workers can be dynamically scaled.

### Adding New Job Types

To add a new job type:

1. Define a new value in the `JobType` enum in `Jobs.jl`
2. Create input and output type definitions that extend `AbstractJobInput` and `AbstractJobOutput`
3. Implement a handler that extends `AbstractJobHandler`
4. Register the handler during application initialization

NOTE: the JobType, input and output should correspond to the typed interfaces in the reefguide-web-api project. This just provides Julia structs around the existing types. The API will reject improperly formed types.

Example:

```julia
# 1. Add to enum
@enum JobType begin
    CRITERIA_POLYGONS
    NEW_JOB_TYPE  # New job type
end

# 2. Define input/output types
struct NewJobInput <: AbstractJobInput
    # Define input fields
    parameter1::String
    parameter2::Int
end

struct NewJobOutput <: AbstractJobOutput
    # Define output fields
    result::String
end

# 3. Implement handler
struct NewJobHandler <: AbstractJobHandler end

function handle_job(
    ::NewJobHandler, input::NewJobInput, storage_uri::String
)::NewJobOutput
    # Implement job processing logic
    result = process_data(input.parameter1, input.parameter2)

    # Return result
    return NewJobOutput(result)
end

# 4. Register in __init__
function __init__()
    # Register existing handlers
    register_job_handler!(
        CRITERIA_POLYGONS,
        CriteriaPolygonsHandler(),
        CriteriaPolygonsInput,
        CriteriaPolygonsOutput
    )

    # Register new handler
    register_job_handler!(
        NEW_JOB_TYPE,
        NewJobHandler(),
        NewJobInput,
        NewJobOutput
    )
end
```

### Launching the Worker

In a Julia REPL, you can launch the worker using `ReefGuideAPI.start_worker()`.

e.g.

#### From the Command Line

```bash
# Set environment variables
export API_ENDPOINT=https://api.reefguide.example.com
export JOB_TYPES=CRITERIA_POLYGONS
export USERNAME=worker-service
export PASSWORD=secure-password

# Start the Julia REPL with the ReefGuide module
cd sandbox
julia --project=.

```

```
pkg > dev ..
pkg > instantiate
julia > using ReefGuideAPI
julia > ReefGuideAPI.start_worker()
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

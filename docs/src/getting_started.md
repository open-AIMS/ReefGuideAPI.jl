# Getting Started


## Setup

Initialize the project the usual way:

```julia
]instantiate
```

A TOML file should be defined indicating location of the MPA dataset.
These are currently the files/data created in Step/Script 1a in https://github.com/open-AIMS/GBR-reef-guidance-assessment

Other settings configure the server

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

::: info

By convention, this file is named `.config.toml` (note the leading `.`).

:::

## Quickstart

```julia
using ReefGuideAPI

# To enable debug messages:
# ENV["JULIA_DEBUG"] = "ReefGuideAPI"

# If multiple threads are available, a parallel server will be spun up
ReefGuideAPI.start_server(".config.toml")
```

::: tip

For best performance, start with at least one interactive thread:

```bash
julia --project=. --threads=4,1
```

:::

::: info

The config setting `COG_THREADS` controls how many threads should be requested when writing
out COGs. Ideally this will be set to at least 2 (preferably 4).
Higher values do seem to reduce write times but with diminishing returns (tested up to 8).
Locally, write times with four threads configured range from 10 to 15 seconds.

:::

In its current state, the main page displays a simple form for dev/testing purposes.

## Dynamic COG generation

Example URL:

```code
http://127.0.0.1:8000/assess/Cairns-Cooktown/slopes?criteria_names=Depth,Slope&lb=-9.0,0.0&ub=-2.0,40.0
```

## Simple Slippy Tiles

Example URL:

```code
http://127.0.0.1:8000/tile/8/231/139?region=Cairns-Cooktown&rtype=slopes&criteria_names=Depth,Slope,Rugosity&lb=-9.0,0.0,0.0&ub=-2.0,40.0,0.0
```

## Development setup

The steps below assumes you are in the project root.

Create a sandbox enviroment:

```bash
$ mkdir sandbox
$ cd sandbox
$ julia --project=.
(sandbox) julia> ]add Revise Infiltrator Chairmarks
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

## JWT Auth configuration

The API can be configured to expect a valid JWT in the `Authorization: Bearer <token>`
header format.

Add the following to your `.config.toml`:

```toml
[jwt_auth]
# Enable JWT auth : bool true/false
JWT_ENABLED = true
# Which iss to validate for the JWTs?
JWT_ISS = "https://issuer.com"
# WKT JWKS endpoint where public key can be retrieved
WKT_ENDPOINT = "https://https://issuer.com/api/.well-known/jwks.json"
```

Pay attention to the issuer and wkt endpoints. The first should exactly match the expected
JWT issuer claim. The second should be web-resolvable and return a WKT JSON which provides
the public key.

### Auth TODOs

- ensure health check route is not authorised
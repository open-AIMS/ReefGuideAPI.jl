# Reef Guidance API

API for supporting reef suitability assessments.

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
CACHE_DIR = "<some location to cache geotiffs>"
DEBUG_MODE = "false"  # optional, disables file caching
```

By convention, this file is named `.config.toml` (note the leading `.`).

## Quickstart

```julia
using ReefGuideAPI

# If multiple threads are available, a parallel server will be spun up
ReefGuideAPI.start_server(".config.toml")

# An example URL to query:
# http://127.0.0.1:8000/assess/Cairns-Cooktown/slopes?criteria_names=Depth,Slope&lb=-9.0,0.0&ub=-2.0,40.0
```

In its current state, the main page displays a simple form for dev/testing purposes.

## Development setup

The steps below assumes you are in the project root.

Create a sandbox enviroment:

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

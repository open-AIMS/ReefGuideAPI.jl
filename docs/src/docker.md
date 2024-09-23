# Docker build and run

The ReefGuideAPI.jl package has an associated `Dockerfile` and build/publish process. This means you can run an instance of the ReefGuideAPI.jl package without needing to compile/build it with a local `Julia` installation. You will be able to view the latest published versions of the Docker image on the repository packages page.

## Mounting files and required data

As mentioned in [Getting Started](@ref getting_started), the `ReefGuideAPI.jl` package currently requires

- a `config.toml` file and
- a set of input data files

Please include these in a folder called `data` in your working directory.

When running the below commands, it is assumed you have `data` available locally with the required files.

**Note**: Due to how Docker excludes `.` files, we have named the config file `config.toml` in the data folder. This is required to launch the server.

## To build from src files using Docker

```bash
docker build . --target reefguide-src -t reefguide
```

## To build from src files using Docker Compose

```bash
docker compose up --build reefguide-src
```

## To run with mounted files (launch server) using Docker

```bash
docker run -p 8000:8000 -v ./data:/data/reefguide reefguide
```

## To run with mounted files (launch server) using Docker Compose

```bash
docker compose up reefguide-src
```

## To run with mounted files (interactive shell) using Docker

This will start a Julia shell where `ReefGuideAPI` is compiled and ready for use e.g.

```julia
using ReefGuideAPI
ReefGuideAPI.start_server("/data/reefguide/config.toml")
```

```bash
docker run --rm --interactive --entrypoint="julia" --tty -v ./data:/data/reefguide reefguide
```

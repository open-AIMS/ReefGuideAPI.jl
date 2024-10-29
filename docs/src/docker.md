# Docker build and run

The ReefGuideAPI.jl package has an associated `Dockerfile` and build/publish process. This
means you can run an instance of the ReefGuideAPI.jl package without needing to
compile/build it with a local `Julia` installation. You will be able to view the latest
published versions of the Docker image on the repository packages page.

## A note about MKL_jll

Due to how Julia (particularly v1.11) handles precompilation, it significantly reduces the build time by explicitly installed MKL_jll before installing any of explicit project dependencies.

For this reason, the project includes a few helpers to optimise this process

- `get_mkl.sh`: this script extracts the `uuid` field of the `[[deps.MKL_jll]]` element in the `Manifest.toml` - it is written into;
- `MKL_jll.dep`: which caches the output of the above script
- `Dockerfile`: the dockerfile includes some explicit steps to preinstall this dependency against the `uuid`

```Dockerfile
# Get the MKL_jll dependency hash
COPY MKL_jll.dep .

# Compile MKL_jll first - this improves build time significantly - unsure exactly why
RUN export MKL_JLL_HASH=$(cat MKL_jll.dep); julia -e 'using Pkg; Pkg.add(PackageSpec(name="MKL_jll", uuid=ENV["MKL_JLL_HASH"])); Pkg.precompile()'
```

Before running a docker build, I recommend running the script

```bash
./get_mkl.sh
```

However PRs will automatically perform this sync, if out of date.

### Github action to sync MKL_jll

There is an action which ensures the MKL_jll dep file is up to date, and pushes if not. See `.github/workflows/CheckMklDep.yml`.

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

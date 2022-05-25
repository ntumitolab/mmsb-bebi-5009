FROM ghcr.io/sosiristseng/docker-juliabook:latest

WORKDIR /usr/src/app

# Build Julia environment
COPY Project.toml Manifest.toml build.jl ./
# COPY src ./src
RUN julia --threads=auto --color=yes --project=@. -e 'import Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile()' && julia --threads=auto --color=yes build.jl

# Julia command
CMD ["julia"]

# If using sysimage
# CMD ["julia", "-Jsysimage.so"]

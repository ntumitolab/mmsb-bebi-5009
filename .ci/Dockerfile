FROM registry.gitlab.com/sosiristseng/docker-juliabook:1.7.3.1

ENV LC_ALL C.UTF-8

WORKDIR /usr/src/app

# Python deps (jupyter-book and matplotlib)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Build Julia environment
COPY Project.toml Manifest.toml build.jl ./
# COPY src/ src
RUN julia --threads=auto --color=yes --project=@. -e 'import Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile()' && julia --threads=auto --color=yes build.jl

# Julia command
CMD ["julia"]

# If using sysimage
# CMD ["julia", "-Jsysimage.so"]

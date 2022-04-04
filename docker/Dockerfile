FROM python:3.10-slim

WORKDIR /usr/src/app

RUN apt-get update && apt-get install -y build-essential && rm -rf /var/lib/apt/lists/*

# Julia
ENV JULIA_PATH /usr/local/julia
ENV PATH $JULIA_PATH/bin:$PATH
COPY --from=julia:1.7.2 ${JULIA_PATH} /usr/local/

# Python deps (jupyter-book and matplotlib)
COPY requirements.txt .
RUN pip install --no-cache-dir -U pip wheel setuptools \
    && pip install --no-cache-dir -r requirements.txt

# Build Julia environment
COPY Project.toml Manifest.toml build_kernel.jl ./
RUN julia --threads=auto --color=yes --project=@. -e 'import Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile()' && julia --threads=auto --color=yes build_kernel.jl

# Julia command
# CMD ["julia"]

# If using sysimage
CMD ["julia", "--sysimage=sysimage.so"]

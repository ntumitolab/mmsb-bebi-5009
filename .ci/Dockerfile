FROM python:3.10-slim

WORKDIR /usr/src/app

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
RUN apt-get update && apt-get install -y build-essential git gzip brotli parallel --no-install-recommends && rm -rf /var/lib/apt/lists/*

# Julia
ENV JULIA_PATH /usr/local/julia
ENV PATH $JULIA_PATH/bin:$PATH
COPY --from=julia:1.7.2 ${JULIA_PATH} /usr/local/

# Python deps (jupyter-book and matplotlib)
COPY requirements.txt .
COPY .ci/requirements.txt .ci/requirements.txt
RUN pip install --no-cache-dir -U pip wheel setuptools \
    && pip install --no-cache-dir -r .ci/requirements.txt -r requirements.txt

# Build Julia environment
COPY Project.toml Manifest.toml build.jl ./
# COPY src ./src
RUN julia --threads=auto --color=yes --project=@. -e 'import Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile()' && julia --threads=auto --color=yes build.jl

# Julia command
CMD ["julia"]

# If using sysimage
# CMD ["julia", "-Jsysimage.so"]

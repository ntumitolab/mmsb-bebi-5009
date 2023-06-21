FROM ghcr.io/sosiristseng/docker-jupyterbook:latest as base

# Supress warnings from parallel
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# Julia config
ENV JULIA_CI true
ENV JULIA_NUM_THREADS "auto"
# Let PythonCall use built-in python
ENV JULIA_CONDAPKG_BACKEND "Null"
# Headless Plots.jl
ENV GKSwstype 140
ENV JULIA_CPU_TARGET "generic"
ENV JULIA_PATH /usr/local/julia/
ENV JULIA_DEPOT_PATH /srv/juliapkg/
ENV PATH ${JULIA_PATH}/bin:${PATH}
COPY --from=julia:1.9.1 ${JULIA_PATH} ${JULIA_PATH}

FROM base

WORKDIR /app

# Python dependencies. e.g. matplotlib
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Julia environment
COPY Project.toml Manifest.toml ./
COPY src/ src
RUN julia --project="" --color=yes -e 'import Pkg; Pkg.add(["Literate", "PrettyTables"])' && \
    julia --project=@. --color=yes -e 'import Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile()'

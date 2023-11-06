FROM julia:1.9.3 as julia
FROM python:3.12.0-slim as base

# Julia config
ENV JULIA_CI true
ENV JULIA_NUM_THREADS 'auto'
# Let PythonCall use built-in python
ENV JULIA_CONDAPKG_BACKEND 'Null'
ENV JULIA_CPU_TARGET 'generic;haswell,clone_all'
ENV GKSwstype '100'
ENV JULIA_PATH /usr/local/julia/
ENV JULIA_DEPOT_PATH /srv/juliapkg/
ENV PATH ${JULIA_PATH}/bin:${PATH}
COPY --from=julia ${JULIA_PATH} ${JULIA_PATH}

FROM base

WORKDIR /work

# Python dependencies. e.g. matplotlib
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Julia environment
COPY Project.toml Manifest.toml ./
COPY src/ src
RUN julia --color=yes -e 'using Pkg; Pkg.add(["PrettyTables", "Literate"]); Pkg.activate("."); Pkg.instantiate()'

FROM python:3.11.4-slim as base

ENV JULIA_CI true
ENV JULIA_NUM_THREADS "auto"
# Let PythonCall use built-in python
ENV JULIA_CONDAPKG_BACKEND "Null"
# Headless Plots.jl
ENV GKSwstype 140
ENV JULIA_PATH /usr/local/julia/
ENV PATH ${JULIA_PATH}/bin:${PATH}
ENV JULIA_CPU_TARGET 'generic;znver3,clone_all'
COPY --from=julia:1.9.1 ${JULIA_PATH} ${JULIA_PATH}

FROM base

# Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Julia environment
COPY Project.toml Manifest.toml ./
COPY src/ src
RUN julia --project="" --color=yes -e 'import Pkg; Pkg.add(["Literate", "PrettyTables"])' && \
    julia --project=@. --color=yes -e 'import Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile()'

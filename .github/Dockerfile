FROM python:3.11.3-slim as base

# Julia
ENV JULIA_CI true
ENV JULIA_CONDAPKG_BACKEND "Null"
ENV JULIA_NUM_THREADS "auto"
ENV JULIA_PATH /usr/local/julia/
ENV JULIA_DEPOT_PATH /srv/juliapkg/
ENV PATH ${JULIA_PATH}/bin:${PATH}
COPY --from=julia:1.9.0-rc2 ${JULIA_PATH} ${JULIA_PATH}

FROM base

WORKDIR /app

RUN pip install --no-cache-dir -U nbconvert matplotlib

# Julia environment
COPY Project.toml Manifest.toml ./
COPY src/ src
RUN julia --color=yes --project="" -e 'import Pkg; Pkg.add("IJulia"); using IJulia; installkernel("Julia", "--project=@.")' && \
    julia --color=yes --project=@. -e 'import Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile()'

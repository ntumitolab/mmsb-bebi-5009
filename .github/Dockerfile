FROM jupyter/base-notebook:python-3.10.10 as base

# Julia
USER root
ENV JULIA_CI true
ENV JULIA_NUM_THREADS "auto"
ENV JULIA_CONDAPKG_BACKEND "System"
ENV JULIA_DEPOT_PATH /srv/juliapkg/
RUN mkdir -p ${JULIA_DEPOT_PATH} && chown ${NB_UID}:${NB_UID} ${JULIA_DEPOT_PATH}
ENV JULIA_PATH /usr/local/julia/
ENV PATH ${JULIA_PATH}/bin:${PATH}
COPY --from=julia:1.8.5 ${JULIA_PATH} ${JULIA_PATH}

FROM base as runtime
USER ${NB_UID}
COPY --chown=${NB_UID}:${NB_UID} Project.toml Manifest.toml ./
COPY --chown=${NB_UID}:${NB_UID} src/ src
RUN julia --color=yes --project="" -e 'import Pkg; Pkg.add("IJulia"); Pkg.build("IJulia")' && \
    julia --color=yes --project=@. -e 'import Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile()' && \
    julia --color=yes --project="" -e 'import Pkg; Pkg.resolve(); Pkg.precompile()'

FROM runtime
COPY --chown=${NB_UID}:${NB_UID} docs/ docs

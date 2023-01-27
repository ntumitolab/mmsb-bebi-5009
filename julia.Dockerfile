FROM julia:1.8.5 AS julia

FROM ghcr.io/sosiristseng/docker-jupyterbook:0.13.1.7

# Julia
ENV JULIA_PATH /usr/local/julia/
ENV PATH ${JULIA_PATH}/bin:${PATH}
COPY --from=julia ${JULIA_PATH} ${JULIA_PATH}

# Python dependencies. e.g. matplotlib
# RUN pip install --no-cache-dir matplotlib

# Julia environment
COPY Project.toml Manifest.toml ./
COPY src/ src
RUN julia --color=yes --project="" -e 'import Pkg; Pkg.add("IJulia"); using IJulia; installkernel("Julia", "--project=@.", env=Dict("JULIA_NUM_THREADS"=>"auto"))' && \
    julia --color=yes --project=@. -e 'import Pkg; Pkg.instantiate()'

CMD ["jupyter-book"]

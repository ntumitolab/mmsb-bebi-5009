FROM python:3.11.1-slim

# Julia
ENV JULIA_CI true
ENV JULIA_PATH /usr/local/julia/
ENV JULIA_DEPOT_PATH /srv/juliapkg/
ENV PATH ${JULIA_PATH}/bin:${PATH}
COPY --from=julia:1.8.5 ${JULIA_PATH} ${JULIA_PATH}

# Python dependencies. e.g. matplotlib
RUN pip install --no-cache-dir nbconvert

# Julia environment
COPY Project.toml Manifest.toml ./
COPY src/ src
RUN julia --color=yes --project="" -e 'import Pkg; Pkg.add("IJulia"); using IJulia; installkernel("Julia", "--project=@.", env=Dict("JULIA_NUM_THREADS"=>"auto"))' && \
    julia --color=yes --project=@. -e 'import Pkg; Pkg.instantiate()'

CMD ["julia"]

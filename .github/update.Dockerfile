FROM python:3.11.1-slim

# Julia
ENV JULIA_CI true
ENV JULIA_PKG_PRECOMPILE_AUTO 0
ENV JULIA_NUM_THREADS auto
ENV JULIA_PATH /usr/local/julia/
ENV JULIA_DEPOT_PATH /srv/juliapkg/
ENV PATH ${JULIA_PATH}/bin:${PATH}
COPY --from=julia:1.8.5 ${JULIA_PATH} ${JULIA_PATH}

CMD ["julia"]

FROM python:3.10-slim

WORKDIR /usr/src/app

# Install apt packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    git openssh-client curl parallel brotli \
    # build-essential \
    && rm -rf /var/lib/apt/lists/*

# Julia version, Download and extact Julia
ENV JL_MAJOR="1.7"
ENV JL_MINOR="2"
ENV JULIA_URL="https://julialang-s3.julialang.org/bin/linux/x64/${JL_MAJOR}/julia-${JL_MAJOR}.${JL_MINOR}-linux-x86_64.tar.gz"
RUN curl -fL ${JULIA_URL} | tar -xzf - -C /usr/local --strip-components 1 \
    && julia --version

# Python deps
COPY requirements.txt .
RUN pip install --no-cache-dir -U pip wheel setuptools && pip install --no-cache-dir -r requirements.txt

# Build Julia environment
COPY Project.toml Manifest.toml ./
RUN julia --threads=auto --color=yes -e 'import Pkg; Pkg.add("IJulia")' \
    && julia --threads=auto --color=yes --project=@. -e 'import Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile()'

# Julia sysimage (optional)
# RUN julia --threads=auto --color=yes -e 'import Pkg; Pkg.add("PackageCompiler"); \
#     create_sysimage(["Example"]; project= ".", sysimage_path="sysimage.so"); \
#     using IJulia; installkernel("Julia", "-J sysimage.so") \
#     Pkg.remove("PackageCompiler"); Pkg.gc()'

# Julia command
CMD ["julia"]

# If using sysimage
# CMD ["julia", "-J", "sysimage.so"]

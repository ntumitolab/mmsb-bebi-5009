FROM ghcr.io/sosiristseng/docker-juliabook:1.7.3.12

ENV LC_ALL C.UTF-8

WORKDIR /usr/src/app

# Python deps (matplotlib)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Build Julia environment
COPY Project.toml Manifest.toml ./
# COPY src/ src
RUN julia --threads=auto --color=yes --project=@. -e 'import Pkg; Pkg.instantiate()'

# Julia command
CMD ["jupyter-book"]

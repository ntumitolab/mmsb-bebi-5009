# AGENTS.md

Course website for *Mathematical Modeling in Systems Biology* (BEBI-5009). Content is authored as Julia/Literate sources under `docs/`, executed into Jupyter notebooks, then built as a Jupyter Book site deployed to GitHub Pages.

## Setup

```sh
julia --project=@. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

- Run Julia from the repo root (or a subdirectory) with `--project=@.`; it finds the root `Project.toml`.
- The Julia version is pinned by `Manifest.toml` (`julia_version`, currently 1.13.2). CI derives its Julia version from that field.
- `Startup` is a local package (`Startup/`, wired via `[sources]` and `[workspace]` in the root `Project.toml`). It exports plotting helpers `meshgrid`, `get_gradient`, `normalize_gradient` used by `docs/figs/ch4-*` and `ch5/ch7` vector-field figures.
- Headless plotting: set `GKSwstype=100` (`.github/ci.jl` does this for workers).
- Python side is only `requirements.txt` (`matplotlib==3.11.2`), installed by CI with `uv`; the site builder itself is npm-installed `jupyter-book`, not a Python dependency.

## Build and verify

```sh
# Run one source directly for fast feedback (does not generate a notebook)
GKSwstype=100 julia --project=@. docs/intro-03-diffeq.jl

# Execute all new/changed docs and refresh the notebook cache
GKSwstype=100 julia --project=@. -p 4 .github/ci.jl
cp -rf .cache/docs/. docs/

# Build the site (requires notebooks generated above)
cd docs && npm install -g jupyter-book && jupyter-book build --html
```

CI (`.github/workflows/ci.yml`) runs the same pipeline on `ubuntu-26.04-arm` with `LITERATE_PROC=4`, then builds with `jupyter-book build --html` and deploys only on pushes to `main`.

- `.github/ci.jl` caches by SHA256 per source file under `.cache/docs/**/<name>.sha`; unchanged files are not re-executed. Delete the matching `.sha` or `.cache/` to force a rerun.
- The cache copy step is required: notebooks are written to `.cache/docs/...`, and the site build only sees them after copying back into `docs/`.
- `docs/myst.yml` is the table of contents. `intro-*.ipynb` and `figs/*.ipynb` are included by glob; new pages outside those patterns must be added explicitly (as `hw-01.ipynb` is).

## Source conventions

- Most files under `docs/` are Literate.jl sources, not runnable Python/Jupyter content:
  - `docs/*.jl`: markdown in `#=== ... ===#` blocks, code in between.
  - `docs/figs/*.jl`: markdown as `#` comments, code blocks separated by `#---`.
- Both forms are valid Julia and can be executed directly.
- Textbook figures mostly build Catalyst.jl/ModelingToolkit.jl ODE models and use the GR backend (`Plots.gr(linewidth=1.5)`); keep the GR backend and explicit `Plots.gr(...)` resets when adding figures.
- `docs/figs/ch1.ipynb`, `ch2.ipynb`, and `ch3.ipynb` are checked-in notebooks with no `.jl` source. On a cache miss, `ci.jl` deletes each `.ipynb`, generates and runs a `.jl` sibling, then copies a regenerated notebook back, so avoid running the full CI script if you only need to edit those notebooks.
- Generated `.ipynb` files under `docs/` are not gitignored (the `*.ipynb` rule is commented out); do not commit build outputs. `*.png` and `*.csv` are ignored.
- `Manifest.toml` is intentionally committed and updated weekly by `.github/workflows/update-manifest.yml`; do not hand-edit it.

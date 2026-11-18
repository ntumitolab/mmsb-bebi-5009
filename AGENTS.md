# AGENTS.md

Course website for *Mathematical Modeling in Systems Biology* (BEBI-5009). Content is authored as Julia/Literate sources under `docs/`, executed into Jupyter notebooks, then built as a Jupyter Book site deployed to GitHub Pages.

## Setup

- To install project dependencies

```sh
julia --project=@. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```
- Prefer `repld` skill to test and run Julia scripts.
- `Startup` is a local package (in `Startup/`). It exports helper functions.
- Headless plotting: set `GKSwstype=100`.

## Source conventions

- Content under `docs/` is authored as Markdown or Literate.jl scripts:
  - `docs/*.jl`: markdown in `#=== ... ===#` blocks, code in between.
  - `docs/figs/*.jl`: markdown as `#` comment lines; `#---` starts a new chunk and preserves the current chunk type.
  - Inside code, comments must start with `##`; Literate treats a single leading `#` as markdown. Both forms are valid Julia and can be executed directly.
- Do not hand-edit `Manifest.toml`. It should be updated by the following command

```sh
julia --project=@. -e 'using Pkg; Pkg.update()'
```

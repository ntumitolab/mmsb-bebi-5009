# AGENTS.md

Course website for *Mathematical Modeling in Systems Biology* (BEBI-5009). Content is authored as Julia/Literate sources under `docs/`, executed into Jupyter notebooks, then built as a Jupyter Book site deployed to GitHub Pages.

## Layout

- `Startup/` — the inner model package.
- `docs/` — Literate notebook sources (`*.jl`) plus `index.md` (model description) and `myst.yml` for the site TOC.
- CI (`.github/workflows/ci.yml`) runs notebooks (by `.github/ci.jl`) and builds the site by `jupyter-book`.
- `extras/` — **read-only archive** of superseded model variants. Never edited, executed, or built by CI. Do not "fix" it; keep it as historical reference.

## Workflow

- Use the repld skill (if available) to run Julia scripts interactively.
- There is no unit test suite — the notebooks *are* the tests. After any change, run the full notebook pipeline: `julia --project=@. --color=yes -p 4 .github/ci.jl`.
- **Cache gotcha:** `ci.jl` caches executed notebooks into `.cache/` keyed by the notebook source hash *only*. After editing `Startup/src/`, delete the relevant `.cache/docs/*.sha` files (or all of `.cache/`) or the changed behavior will not be exercised. CI gets this right by including `Startup/**` and the manifests in its cache key.
- Do not hand-edit `Manifest.toml`; it is updated by Pkg — manually with `julia --project=@. --color=yes -e 'using Pkg; Pkg.update()'`, or by the weekly `update-manifest.yml` bot PR, whose command additionally runs `Pkg.add(["JSON", "Literate", "MarkdownTables", "SHA", "Tables"])` before `Pkg.update()`. Keep those packages as direct deps of the root `Project.toml`: `.github/ci.jl` does `using Tables, MarkdownTables` at top level and `using Literate, JSON` on the workers.

## Conventions

- Keep things in one function unless composable or reusable.
- Do not extract single-use helpers preemptively. Inline the logic at the call site unless the helper is reused, hides a genuinely complex boundary, or has a clear independent name that improves the caller.
- Before adding complexity for a speculative or vanishingly unlikely race or security edge case, explain the concrete failure mode, likelihood, and complexity cost to the user and get their buy-in. Do not silently expand scope for theoretical robustness.
- Notebook sources are Literate files: markdown blocks are `#===` … `===#`, code cells are separated by `#---`. `index.md` is static markdown (never executed) and feeds the site TOC alongside the built notebooks (`myst.yml`).
- **Literate gotcha:** comments inside a code cell must start with `##` (rendered as `# ` in the output). Any line starting with a single `# ` is parsed as markdown and splits the code chunk there — if that lands inside a multiline expression, Literate cuts the block in half and execution errors.
- Plots use GR; headless runs set `ENV["GKSwstype"] = "100"` (ci.jl does this for workers).

# Mathematical Modeling in Systems Biology

## Commands

### Install Julia dependencies without updating

Requires `julia` to be installed.

```bash
julia --project=. --color=yes -e 'using Pkg; Pkg.instantiate()'
```

### Update Julia dependencies

Requires `julia` to be installed.

```bash
julia --project=. --color=yes -e 'using Pkg; Pkg.update()'
```

### Run all the notebooks locally

Requires
- Julia dependencies installed
- Jupyter `nbconvert`
- GNU `parallel`

```bash
find . -type f -name '*.ipynb' -print0 | parallel -0 -j$(nproc) jupyter nbconvert --to notebook --ExecutePreprocessor.timeout=600 --execute --inplace {}
```

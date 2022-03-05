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

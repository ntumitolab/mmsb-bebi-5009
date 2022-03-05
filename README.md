# Mathematical Modeling in Systems Biology

![GitHub repo size](https://img.shields.io/github/repo-size/NTUMitoLab/mmsb-bebi-5009?style=flat-square)

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

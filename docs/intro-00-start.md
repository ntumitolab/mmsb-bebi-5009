# Start using Julia

See [Julia in Visual Studio Code](https://code.visualstudio.com/docs/languages/julia) for details.

## Tools

- Install [Julia](https://julialang.org/downloads/) in the official website or via the [Windows store](https://www.microsoft.com/zh-tw/p/julia/9njnww8pvkmn?rtc=1&activetab=pivot:overviewtab).
- Install [VS Code](https://code.visualstudio.com/).
- Install the [Julia extension for VS Code](https://www.julia-vscode.org/). See also [Julia extension's website]().
- (Optional) Install the [Jupyter extension for VS Code](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter) for opening and running Jupyter notebooks in VS Code.

## Running Julia in VS Code

Julia extension will be automatically activated upon opening a Julia file (`*.jl`).

You can open the command palette (`Ctrl + Shift + P` in windows) and search "Julia" for available commands and keyboard shortcuts. The most used on is `Shift + Enter`: to execute the current selected line.

## Package management in Julia

[Pkg in Julia docs](https://docs.julialang.org/en/v1/stdlib/Pkg/)

In regular environments `Project.toml` and `Manifest.toml` describe the dependencies.

### Install packages

In the Julia script

```julia
using Pkg

# Using the current directory as the project folder
Pkg.activate(".")

# Or: Using the directory where the script resides as the project folder
Pkg.activate(@__DIR__)

# Add packages
Pkg.add("Plots")
```

### See installed packages

```julia
using Pkg
Pkg.status()
```

### Remove packages

```julia
using Pkg
Pkg.remove("Plots")
```

### Update installed packages

```julia
using Pkg
Pkg.update()
```

### Create / Use a environment

```julia

using Pkg

# Activate environment in the foldername directory
Pkg.activate("foldername")

# Or activate the current working directory
# current_project() is a little bit misleading
# since it actually looks for an available Project.toml file
Pkg.activate(".")

# Reproduce the environment
Pkg.instantiate()

# If the above failed, try this
Pkg.resolve()
```

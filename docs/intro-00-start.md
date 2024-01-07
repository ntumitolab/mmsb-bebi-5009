# Start using Julia

See [Julia in Visual Studio Code](https://code.visualstudio.com/docs/languages/julia) for details.

## Installation

- [Install Julia](https://julialang.org/downloads/) in the official website or via the [Windows store](https://www.microsoft.com/zh-tw/p/julia/9njnww8pvkmn?rtc=1&activetab=pivot:overviewtab).
- [Install VS Code](https://code.visualstudio.com/).
- [Install Julia extension for VS Code](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia). See alos [Julia extension's website](https://www.julia-vscode.org/).
- (Optional) [Install Jupyter extension for VS Code](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter) for opening and running Jupyter notebooks in VS Code.

## Running Julia in VS Code

Julia extension will be automatically activated upon opening a Julia file (`*.jl`).

You can open the command pallete (`Ctrl + Shift + P` in windows) and search "Julia" for available commands and keyboard shortcuts. The most used on is `Shift + Enter`: to execute the current selected line.

## Package management in Julia

[Pkg in Julia docs](https://docs.julialang.org/en/v1/stdlib/Pkg/)

In regular environments `Project.toml` and `Manifest.toml` describe the dependencies.

### Install packages

In the Julia script

```julia
using Pkg

# Function form
Pkg.add("Plots")

# Or use Pkg's special strings
pkg"add Plots"
```


In the Julia REPL:

```julia-repl
] add Plots
```

### See installed packages

```julia
using Pkg
Pkg.status()

# Or
pkg"st" # pkg"status"

```

In the Julia REPL:

```julia-repl
] st
```

### Remove packages

```julia
using Pkg
Pkg.remove("Plots")

# Or
pkg"rm Plots"
```

```julia-repl
] rm Plots
```

### Update installed packages

```julia
using Pkg
Pkg.update()

# Or
pkg"up" # pkg"update"
```

### Create / Use a environment

```julia

using Pkg

# Activate environment in the foldername directory
Pkg.activate("foldername")

# Or activate the current working directory
# current_project() is a little bit misleading
# since it actually looks for available Project.toml file
Pkg.activate(Base.current_project())

# Install the packages according to the environment
Pkg.instantiate()

# If the above failed, try this
Pkg.resolve()
```

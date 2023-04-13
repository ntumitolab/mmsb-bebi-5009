# # Prob 2.4.6

using DifferentialEquations
using Plots
Plots.default(linewidth=2)
import DisplayAs.SVG

fig = ODEProblem((u, p, t) -> p * (1. - u), 0., 10., 1.) |> solve |> plot

fig |> SVG

# ## Runtime information
import Pkg
Pkg.status()

#---
import InteractiveUtils
InteractiveUtils.versioninfo()

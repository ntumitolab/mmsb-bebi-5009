# # Prob 2.4.6
using OrdinaryDiffEq
using Plots
Plots.default(linewidth=2)

#---
ODEProblem((u, p, t) -> p * (1. - u), 0., 10., 1.) |> solve |> plot

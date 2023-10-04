# # Prob 2.4.6

using DifferentialEquations
using Plots
Plots.default(linewidth=2)

fig = ODEProblem((u, p, t) -> p * (1. - u), 0., 10., 1.) |> solve |> plot

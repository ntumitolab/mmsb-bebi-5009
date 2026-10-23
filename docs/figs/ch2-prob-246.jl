# # Problem 2.4.6
using OrdinaryDiffEq
using Plots

sol246 = ODEProblem((u, p, t) -> p * (1.0 - u), 0.0, 10.0, 1.0) |> solve
plot(sol246, title="Problem 2.4.6", xlabel="Time", ylabel="Concentration")

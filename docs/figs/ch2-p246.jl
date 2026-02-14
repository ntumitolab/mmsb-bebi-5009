# # Prob 2.4.6
using OrdinaryDiffEq
using Plots
Plots.gr(framestyle = :box, linewidth=1.5)
#---
sol = ODEProblem((u, p, t) -> p * (1.0 - u), 0.0, 10.0, 1.0) |> solve

# Fig 2.46
plot(sol, xlabel="Time", ylabel="Concentration", title="Prob. 2.4.6", legend=nothing)

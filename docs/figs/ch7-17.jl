# # Fig 7.17
# Goodwin oscillator model: https://en.wikipedia.org/wiki/Goodwin_model_(biology)
using Startup: meshgrid, get_gradient, normalize_gradient
using DifferentialEquations
using DiffEqCallbacks
using DifferentialEquations
using Catalyst
using Plots
Plots.gr(linewidth=1.5)

@time "Build system" rn717 = @reaction_network begin
    (a / (k^n + Z^n), b), 0 <--> X
    (α * X, β), 0 <--> Y
    (γ * Y, δ), 0 <--> Z
end

#---
alg = FBDF()
up717 = Dict(:a => 360.0, :k => 1.368, :b => 1.0, :α => 1.0, :β => 0.6, :γ => 1.0, :δ => 0.8, :n => 12.0, :X => 0.0, :Y => 0.0, :Z => 0.0)

tend = 35.0
@time "Build problem" prob717 = ODEProblem(rn717, up717, (0.0, tend))
@time "Solve problem" sol717 = solve(prob717, alg)
plot(sol717, idxs=[:X, :Y, :Z], labels=["X" "Y" "Z"], title="Fig 7.17 (A)", xlabel="Time", ylabel="Concentration")

#---
plot(sol717, idxs=(:X, :Y, :Z), labels=false, title="Fig 7.17 (B)", size=(600, 600), xlabel="X", ylabel="Y", zlabel="Z")

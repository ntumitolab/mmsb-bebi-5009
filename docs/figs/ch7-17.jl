#===
# Fig 7.17

Goodwin oscillator model: https://en.wikipedia.org/wiki/Goodwin_model_(biology)
===#
using ComponentArrays
using SimpleUnPack
using OrdinaryDiffEq
using Plots
Plots.default(linewidth=2)

#---
function model717!(D, u, p, t)
    @unpack a, k, b, α, β, γ, δ, n = p
    @unpack X, Y, Z = u
    D.X = a / (k^n + Z^n) - b * X
    D.Y = α * X - β * Y
    D.Z = γ * Y - δ * Z
    nothing
end

#---
ps717 = ComponentArray(
    a = 360.0,
    k = 1.368,
    b = 1.0,
    α = 1.0,
    β = 0.6,
    γ = 1.0,
    δ = 0.8,
    n = 12.0
)

u0 = ComponentArray(
    X = 0.0,
    Y = 0.0,
    Z = 0.0
)

tend = 35.0
prob717 = ODEProblem(model717!, u0, tend, ps717)

#---
@time sol = solve(prob717)

#---
plot(sol, title="Fig 7.17 (A)", xlabel="Time", ylabel="Concentration", labels=["X" "Y" "Z"])

#---
plot(sol, idxs=(1, 2, 3), title="Fig 7.17 (B)", legend=false, size=(600, 600))

#===
# Fig 4.18

Continuation diagram

See also [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl)
===#
using OrdinaryDiffEq
using SteadyStateDiffEq
using ComponentArrays: ComponentArray
using SimpleUnPack
using CairoMakie

# Model
function model418!(D, u, p, t)
    @unpack k1, k2, k3, k4, k5, n = p
    @unpack A, B = u
    D.A = k1 / (1 + B^n) - (k3 + k5) * A
    D.B = k2 + k5 * A - k4 * B
    return nothing
end

#---
ps418 = ComponentArray(
    k1 = 0.0,
    k2 = 5.0,
    k3 = 5.0,
    k4 = 5.0,
    k5 = 2.0,
    n = 4.0
)
u0418 = ComponentArray(
    A = 0.0,
    B = 0.0
)

prob = SteadyStateProblem(model418!, u0418, ps418)

function ainf(k1val)
    sol = solve(remake(prob, p=ComponentArray(ps418; k1=k1val)))
    return sol.u[1]
end

#---
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "K1",
    ylabel = "Steady state [A]",
    title = "Fig 4.18\nContinuation diagram"
)
lines!(ax, 0.0..1000.0, k1 -> ainf(k1), color=:blue)
fig

#===
# Fig 2.09

Metabolic network simulation
===#
using OrdinaryDiffEq
using ComponentArrays
using SimpleUnPack
using Plots
Plots.default(linewidth=2)

#---
function model209!(du, u, p, t)
    @unpack k1, k2, k3, k4, k5 = p
    @unpack A, B, C, D = u
    v1 = k1
    v2 = k2 * A
    v3 = k3 * A * B
    v4 = k4 * C
    v5 = k5 * D
    du.A = v1 - v2 - v3
    du.B = v2 - v3
    du.C = v3 - v4
    du.D = v3 - v5
    nothing
end

# Setup problem
ps = ComponentArray(
    k1 = 3.0,
    k2 = 2.0,
    k3 = 2.5,
    k4 = 3.0,
    k5 = 4.0
)
u0 = ComponentArray(
    A = 0.0,
    B = 0.0,
    C = 0.0,
    D = 0.0
)
tend = 10.0

prob = ODEProblem(model209!, u0, tend, ps)

# Solve the problem
@time sol = solve(prob, Tsit5())

# Visualize the results
plot(sol, legend=:bottomright, title="Fig 2.9",
    xlims=(0., 4.), ylims=(0., 1.),
    xlabel="Time (sec)", ylabel="Concentration (mM)",
    labels=["A" "B" "C" "D"]
)

#===
# Fig 2.09

Metabolic network simulation
===#
using OrdinaryDiffEq
import ComponentArrays as CA
using SimpleUnPack
using CairoMakie

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
ps = CA.ComponentArray(
    k1=3.0,
    k2=2.0,
    k3=2.5,
    k4=3.0,
    k5=4.0
)
u0 = CA.ComponentArray(
    A=0.0,
    B=0.0,
    C=0.0,
    D=0.0
)
tend = 10.0
prob = ODEProblem(model209!, u0, tend, ps)

# Solve the problem
@time sol = solve(prob, Tsit5())

# Visualize the results
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel="Time",
    ylabel="Concentration",
    title="Fig 2.9\nMetabolic Network Simulation"
)
lines!(ax, 0 .. tend, t -> sol(t).A, label="A")
lines!(ax, 0 .. tend, t -> sol(t).B, label="B")
lines!(ax, 0 .. tend, t -> sol(t).C, label="C")
lines!(ax, 0 .. tend, t -> sol(t).D, label="D")
axislegend(ax, position=:rb)

fig

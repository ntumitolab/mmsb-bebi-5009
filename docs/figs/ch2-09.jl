# # Fig 2.09
# Metabolic network simulation
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
Plots.gr(framestyle = :box, linewidth=1.5)
#---
function prob209(; tend = 10.0)
    @parameters k1=3.0 k2=2.0 k3=2.5 k4=3.0 k5=4.0
    @variables a(t)=0.0 b(t)=0.0 c(t)=0.0 d(t)=0.0
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
    k1=3.0,
    k2=2.0,
    k3=2.5,
    k4=3.0,
    k5=4.0
)
u0 = ComponentArray(
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

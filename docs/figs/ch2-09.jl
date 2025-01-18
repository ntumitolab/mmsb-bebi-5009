#===
# Fig 2.09

Metabolic network simulation

Using `Catalyst.jl` to simulate a metabolic network.
===#
using OrdinaryDiffEq
using Catalyst
using ModelingToolkit
using Plots
Plots.default(linewidth=2)

#---
rn = @reaction_network begin
    k1, 0 --> A
    k2, A --> B
    k3, A + B --> C + D
    k4, C --> 0
    k5, D --> 0
end

# Showing the differential equations in the reaction network
osys = convert(ODESystem, rn) |> complete
equations(osys)

# Solve the problem
ps = [:k1=>3., :k2=>2., :k3=>2.5, :k4=>3., :k5=>4.]
u0 = [:A=>0., :B=>0., :C=>0., :D=>0.]
tend = 10.
prob = ODEProblem(osys, u0, tend, ps)
sol = solve(prob)

# Visual
plot(sol, legend=:bottomright, title="Fig 2.9",
    xlims=(0., 4.), ylims=(0., 1.),
    xlabel="Time (sec)", ylabel="Concentration (mM)"
)

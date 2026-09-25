# # Fig 4.18
# Continuation diagram.
# See also [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl) for a more robust continuation diagram.
using DifferentialEquations
using SteadyStateDiffEq
using Catalyst
using Plots
Plots.gr(linewidth=1.5)

@time "Build system" model415 = @reaction_network begin
    (k1 / (1 + B^n), k3), 0 <--> A
    k5, A --> B
    (k2, k4), 0 <--> B
end

#---
up418 = Dict(:k1 => 0.0, :k2 => 5.0, :k3 => 5.0, :k4 => 5.0, :k5 => 2.0, :n => 4.0, :A => 0.0, :B => 0.0)
@time "Build problem" prob = SteadyStateProblem(model415, up418)
k1s = 0.0:10.0:1000.0
@time sols = [solve(remake(prob, p=[:k1 => k1]), DynamicSS(FBDF()))[:A] for k1 in k1s]
plot(k1s, sols, xlabel="Parameter k1", ylabel="Steady state [A]", title="Fig. 4.18", label=false)

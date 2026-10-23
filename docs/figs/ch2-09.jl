# Fig 2.09
# Metabolic network simulation
using OrdinaryDiffEq
using Catalyst
using Plots

@time "Build system" rn209 = @reaction_network model209 begin
    k1, 0 --> A
    k2, A --> B
    k3, A + B --> C + D
    k4, C --> 0
    k5, D --> 0
end

#---
up209 = (:k1 => 3.0, :k2 => 2.0, :k3 => 2.5, :k4 => 3.0, :k5 => 4.0, :A => 0.0, :B => 0.0, :C => 0.0, :D => 0.0)
tend = 10.0
@time "Build problem" prob209 = ODEProblem(rn209, up209, tend)
@time "Solve problem" sol = solve(prob209, Tsit5())
plot(sol, title="Fig 2.9", xlabel="Time", ylabel="Concentration")

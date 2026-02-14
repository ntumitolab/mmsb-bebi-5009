#===
# Fig 2.09

Metabolic network simulation
===#
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
    v2 = k2 * a
    v3 = k3 * a * b
    v4 = k4 * c
    v5 = k5 * d
    eqs = [
        D(a) ~ v1 - v2 - v3
        D(b) ~ v2 - v3
        D(c) ~ v3 - v4
        D(d) ~ v3 - v5
    ]
    @mtkcompile sys = ODESystem(eqs, t)
    return ODEProblem(sys, [], tend)
end

@time prob = prob209()

# Solve the problem
@time sol = solve(prob, Tsit5())

# Visualize the results
plot(sol, xlabel="Time", ylabel="Concentration", title="Fig 2.9\nMetabolic Network Simulation", legend=:bottomright)

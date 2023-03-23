# # Chapter 2
# ## Fig 2.04 Exponential decay

using Plots
Plots.default(linewidth=2)

plot(title= "Figure 2.4 Exponential decay")
for k in 1:3
    plot!(t -> 3 * exp(-k*t), 0., 5., label = "exp(-$(k)t)")
end

plot!(xlim = (0, 5), ylim=(0, 3.2), xlabel="Time", ylabel="Concentration")

#===
## Fig 2.09 Metabolic network simulation

Using `Catalyst.jl` to simulate a metabolic network.
===#

using DifferentialEquations
using Catalyst
using ModelingToolkit
using Plots
Plots.default(linewidth=2)

rn = @reaction_network begin
    k1, 0 --> A
    k2, A --> B
    k3, A + B --> C + D
    k4, C --> 0
    k5, D --> 0
end

#---
osys = convert(ODESystem, rn)
for eq in osys.eqs
    println(eq)
end

#---
ps = [:k1 => 3., :k2 => 2., :k3 => 2.5, :k4 => 3., :k5 => 4. ]
u0 = [:A=>0., :B=>0., :C=>0., :D=>0.]
tend = 10.
sol = solve(ODEProblem(rn, u0, tend, ps))

plot(sol, legend=:bottomright, title="Figure 2.09 Metabolic network",
    xlims=(0., 4.), ylims=(0., 1.),
    xlabel="Time (sec)", ylabel="Concentration (mM)")

# ## Figure 2.11
# Model reduction of ODE metabolic networks.

using DifferentialEquations
using Catalyst
using ModelingToolkit
using Plots
Plots.default(linewidth=2)

rn211 = @reaction_network begin
    k0, 0 --> A
    (k1, km1), A <--> B
    k2, B --> 0
end

#--
@unpack k0, k1, km1, k2, A, B = rn211
ps1 = [k0=>0., k1=>9., km1=>12., k2=>2.]
u0 = [A=>0., B=>10.]
tend = 3.0
sol211 = solve(ODEProblem(rn211, u0, tend, ps1))

#--
plot(
    sol211,
    xlabel="Time (AU)",
    ylabel="Concentration (AU)",
    title="Fig. 2.11 (Full model)"
)

# ## Figure 2.12 : Rapid equilibrium assumption

function make_212(;name)
    @parameters k0 k1 km1 k2
    @variables t
    @variables A(t) B(t) C(t)
    D = Differential(t)
    eqs = [
        C ~ A + B
        B ~ C * k1 / (km1 + k1)
        D(C) ~ k0 - k2 * B
    ]
    sys = ODESystem(eqs; name)
    structural_simplify(sys)
end

#--
@named model212 = make_212()

#--
states(model212)

#--
observed(model212)

#--
parameters(model212)

#--
independent_variables(model212)

#--
@unpack k0, k1, km1, k2, C, A, B = model212
ps1 = [k0=>0., k1=>9., km1=>12., k2=>2.]
u0 = [C=>10.]
tend = 3.
prob = ODEProblem(model212, u0, tend, ps1)
sol212 = solve(prob)

#--
plot(sol211, line=(:dash, 1), label=["A (full solution)" "B (full solution)"])
plot!(sol212, idxs=[A, B], label=["A (rapid equilibrium)" "B (rapid equilibrium)"])
plot!(
    title="Fig. 2.12 (Rapid equilibrium model)",
    xlabel="Time (AU)",
    ylabel="Concentration (AU)"
)

#===
## Figure 2.13: Rapid equilibrium (take 2)

When another set of parameters is not suitable for rapid equilibrium assumption.
===#

ps2 = [k0=>9., k1=>20., km1=>12., k2=>2.]
u0 = [A=>8., B=>4.]
tend = 3.0

sol213full = solve(ODEProblem(rn211, u0, tend, ps2))
sol213re = solve(ODEProblem(model212, [C => sum(last.(u0))], tend, ps2))

plot(sol213full, line=(:dash, 1), label=["A (full solution)" "B (full solution)"])
plot!(sol213re, idxs=[A, B], label=["A (rapid equilibrium)" "B (rapid equilibrium)"])
plot!(
    title="Fig. 2.13 (Rapid equilibrium model)",
    xlabel="Time (AU)",
    ylabel="Concentration (AU)"
)

#===
## Figure 2.14 : QSSA
Quasi-steady state assumption on species A
===#

function make_214(;name)
    @parameters k0 k1 km1 k2
    @variables t
    @variables A(t) B(t)
    D = Differential(t)
    eqs = [
        A ~ (k0 + km1 * B)/k1
        D(B) ~ k1 * A - (km1 + k2) * B
    ]
    sys = ODESystem(eqs; name)
    structural_simplify(sys)
end

@named model214 = make_214()

# Initial condations could also be expressed in ModelingToolkit symbols

sol214 = solve(ODEProblem(model214, [B => (k1 * sum(last.(u0)) - k0) / (k1 + km1)], tend, ps2))

plot(sol213full, line=(:dash))
plot!(sol214, idxs=[A, B], label=["A (QSSA)" "B (QSSA)"])
plot!(
    xlabel="Time (AU)",
    ylabel="Concentration (AU)",
    title="Figure 2.14: Ref vs QSSA",
    xlims=(0.0, tend)
)

# ## Problem 2.4.6

using DifferentialEquations
using Plots
Plots.default(linewidth=2)

# Using pipe operator |>
ODEProblem((u, p, t) -> p * (1. - u), 0., 10., 1.) |> solve |> plot

# ## Runtime information

import InteractiveUtils
InteractiveUtils.versioninfo()

#---

import Pkg
Pkg.status()

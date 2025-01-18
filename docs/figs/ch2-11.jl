#===
# Fig 2.11, 2.12, 2.13, 2.14

Model reduction of ODE metabolic networks.
===#
using OrdinaryDiffEq
using Catalyst
using ModelingToolkit
using Plots
Plots.default(linewidth=2)

#---
rn211 = @reaction_network begin
    k0, 0 --> A
    (k1, km1), A <--> B
    k2, B --> 0
end

#---
@unpack k0, k1, km1, k2, A, B = rn211
ps1 = [k0 => 0.0, k1 => 9.0, km1 => 12.0, k2 => 2.0]
u0 = [A => 0.0, B => 10.0]
tend = 3.0
sol211 = solve(ODEProblem(rn211, u0, tend, ps1))

# Fig 2.11
plot(
    sol211,
    xlabel="Time (AU)",
    ylabel="Concentration (AU)",
    title="Fig. 2.11 (Full model)"
)


# ## Figure 2.12 : Rapid equilibrium assumption

function make_212(; name)
    @parameters k0 k1 km1 k2
    @independent_variables t
    @variables A(t) B(t) C(t)
    D = Differential(t)
    eqs = [
        C ~ A + B
        B ~ C * k1 / (km1 + k1)
        D(C) ~ k0 - k2 * B
    ]
    return ODESystem(eqs, t; name)
end

#---
@mtkbuild model212 = make_212()

#---
unknowns(model212)

#---
observed(model212)

#---
parameters(model212)

#---
independent_variables(model212)

#---
@unpack k0, k1, km1, k2, C, A, B = model212
ps1 = [k0 => 0.0, k1 => 9.0, km1 => 12.0, k2 => 2.0]
u0 = [C => 10.0]
tend = 3.0
prob = ODEProblem(model212, u0, tend, ps1)
sol212 = solve(prob)

#---
fig = plot(sol211, line=(:dash, 1), label=["A (full solution)" "B (full solution)"])
plot!(fig, sol212, idxs=[A, B], label=["A (rapid equilibrium)" "B (rapid equilibrium)"])
plot!(fig,
    title="Fig. 2.12 (Rapid equilibrium model)",
    xlabel="Time (AU)",
    ylabel="Concentration (AU)"
)

fig

#===
## Figure 2.13: Rapid equilibrium (take 2)

When another set of parameters is not suitable for rapid equilibrium assumption.
===#

ps2 = [k0 => 9.0, k1 => 20.0, km1 => 12.0, k2 => 2.0]
u0 = [A => 8.0, B => 4.0]
tend = 3.0

sol213full = solve(ODEProblem(rn211, u0, tend, ps2))
sol213re = solve(ODEProblem(model212, [C => sum(last.(u0))], tend, ps2))

fig = plot(sol213full, line=(:dash, 1), label=["A (full solution)" "B (full solution)"])
plot!(fig, sol213re, idxs=[A, B], label=["A (rapid equilibrium)" "B (rapid equilibrium)"])
plot!(fig,
    title="Fig. 2.13 (Rapid equilibrium model)",
    xlabel="Time (AU)",
    ylabel="Concentration (AU)"
)

fig

#===
## Figure 2.14 : QSSA

Quasi-steady state assumption on species A
===#

function make_214(; name)
    @parameters k0 k1 km1 k2
    @independent_variables t
    @variables A(t) B(t)
    D = Differential(t)
    eqs = [
        A ~ (k0 + km1 * B) / k1
        D(B) ~ k1 * A - (km1 + k2) * B
    ]
    return ODESystem(eqs, t; name)
end

#---
@mtkbuild model214 = make_214()

# Initial conditions can also be represented in symbols
sol214 = solve(ODEProblem(model214, [B => (k1 * sum(last.(u0)) - k0) / (k1 + km1)], tend, ps2))

fig = plot(sol213full, line=(:dash))
plot!(fig, sol214, idxs=[A, B], label=["A (QSSA)" "B (QSSA)"])
plot!(fig,
    xlabel="Time (AU)",
    ylabel="Concentration (AU)",
    title="Figure 2.14: Ref vs QSSA",
    xlims=(0.0, tend)
)

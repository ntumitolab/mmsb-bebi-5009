#===
# Fig 2.11, 2.12, 2.13, 2.14

Model reduction of ODE metabolic networks.
===#
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots

#---
function model211(; tend = 3.0)
    @parameters k0=0.0 k1=9.0 km1=12.0 k2=2.0
    @variables a(t)=0.0 b(t)=10.0
    v0 = k0
    v1 = k1 * a - km1 * b
    v2 = k2 * b
    eqs = [
        D(a) ~ v0 - v1
        D(b) ~ v1 - v2
    ]
    @mtkcompile sys = ODESystem(eqs, t)
    return ODEProblem(sys, [], tend)
end

@time prob = model211()

#---
@time sol211 = solve(prob, Tsit5())

# Fig 2.11
plot(sol211, xlabel="Time", ylabel="Concentration", title="Fig. 2.11\nFull model")

# ## Figure 2.12
# Rapid equilibrium assumption
function model212(; tend = 3.0)
    @parameters k0=0.0 k1=9.0 km1=12.0 k2=2.0
    @variables a(t) b(t) c(t)=10.0
    v0 = k0
    v2 = k2 * b
    eqs = [
        a ~ c * km1 / (km1 + k1)
        b ~ c * k1 / (km1 + k1)
        D(c) ~ v0 - v2
    ]
    @mtkcompile sys = ODESystem(eqs, t)
    return ODEProblem(sys, [], tend)
end

@time prob2 = model212()

#---
@time sol212 = solve(prob2, Tsit5())
#---
@unpack a, b, c = prob2.f.sys
plot(sol212, idxs=[a, b], xlabel="Time", ylabel="Concentration", title="Fig. 2.12\nRapid equilibrium model")
plot!(sol211, label=["a (full solution)" "b (full solution)"], linestyle=:dash)

#===
## Figure 2.13

Rapid equilibrium (take 2)

When another set of parameters is not suitable for rapid equilibrium assumption.
===#
@unpack k0, k1, km1, k2 = prob2.f.sys
prob213full = remake(prob, p=[k0 => 9, k1 => 20, km1 => 12, k2 => 2], u0=[a => 8, b => 4])
prob213re = remake(prob2, p=[k0 => 9, k1 => 20, km1 => 12, k2 => 2], u0=[c => 12])

@time sol213full = solve(prob213full, Tsit5())
@time sol213re = solve(prob213re, Tsit5())
#---
plot(sol213re, idxs=[a, b], xlabel="Time", ylabel="Concentration", title="Fig. 2.13\nRapid equilibrium model (take 2)")
plot!(sol213full, label=["a (full solution)" "b (full solution)"], linestyle=:dash)

#===
## Figure 2.14 : QSSA

Quasi-steady state assumption on species A
===#
function model214(u0; tend = 3.0)
    @parameters k0=9 k1=20 km1=12 k2=2
    @variables a(t) b(t)
    eqs = [
        a ~ (k0 + km1 * b) / k1 ## Quasi-steady state assumption
        D(b) ~ k1 * a - (km1 + k2) * b
    ]
    @mtkcompile sys = ODESystem(eqs, t)
    return ODEProblem(sys, [b => (k1 * u0 - k0) / (k1 + km1)], tend)
end

@time prob214 = model214(12)

# Solve QSSA model
@time sol214 = solve(prob214, Tsit5())

#---
plot(sol214, idxs=[a, b], xlabel="Time", ylabel="Concentration", title="Fig. 2.14\nQSSA model")
plot!(sol213full, label=["a (full solution)" "b (full solution)"], linestyle=:dash)

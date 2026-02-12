#===
# Fig 2.11, 2.12, 2.13, 2.14

Model reduction of ODE metabolic networks.
===#
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
Plots.pythonplot()

#---
function prob211(; tend = 3.0)
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

@time prob = prob211()

#---
@time sol211 = solve(prob, Tsit5())

# Fig 2.11
plot(sol211, xlabel="Time", ylabel="Concentration", title="Fig. 2.11\nFull model")

# ## Figure 2.12
# Rapid equilibrium assumption
function prob212(; tend = 3.0)
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

@time prob2 = prob212()

#---
@time sol212 = solve(prob2, Tsit5())
#---
@unpack a, b = prob2.f.sys
plot(sol212, idxs=[a, b], xlabel="Time", ylabel="Concentration", title="Fig. 2.12\nRapid equilibrium model")
plot!(sol211, label=["a (full solution)" "b (full solution)"], linestyle=:dash)

#===
## Figure 2.13

Rapid equilibrium (take 2)

When another set of parameters is not suitable for rapid equilibrium assumption.
===#
@unpack k0, k1, km1, k2, c = prob2.f.sys
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
_a214(u, p, t) = (p.k0 + p.km1 * u.B) / p.k1
_u0214(u0, p) = (p.k1 * u0 - p.k0) / (p.k1 + p.km1)
function model214!(du, u, p, t)
    @unpack k0, k1, km1, k2 = p
    @unpack B = u
    A = _a214(u, p, t)
    du.B = k1 * A - (km1 + k2) * B
    nothing
end

ps214 = ps213
u0214 = ComponentArray(B=_u0214(sum(u0), ps214))

# Solve QSSA model
@time sol214 = solve(ODEProblem(model214!, u0214, tend, ps214), Tsit5())

#---
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel="Time",
    ylabel="Concentration",
    title="Fig. 2.14\nQSSA vs Full model"
)
lines!(ax, 0 .. tend, t -> sol213full(t).A, label="A (full solution)", linestyle=:dash)
lines!(ax, 0 .. tend, t -> sol213full(t).B, label="B (full solution)", linestyle=:dash)
lines!(ax, 0 .. tend, t -> _a214(sol214(t), ps214, t), label="A (QSSA)")
lines!(ax, 0 .. tend, t -> sol214(t).B, label="B (QSSA)")
axislegend(ax, position=:rt)
fig

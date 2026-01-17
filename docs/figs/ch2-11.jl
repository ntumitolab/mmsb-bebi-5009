#===
# Fig 2.11, 2.12, 2.13, 2.14

Model reduction of ODE metabolic networks.
===#
using OrdinaryDiffEq
import ComponentArrays as CA
using SimpleUnPack
using CairoMakie

#---
function model211!(du, u, p, t)
    @unpack k0, k1, km1, k2 = p
    @unpack A, B = u
    v0 = k0
    v1 = k1 * A - km1 * B
    v2 = k2 * B
    du.A = v0 - v1
    du.B = v1 - v2
    nothing
end

#---
ps211 = CA.ComponentArray(
    k0=0.0,
    k1=9.0,
    km1=12.0,
    k2=2.0
)

u0 = CA.ComponentArray(
    A=0.0,
    B=10.0
)

tend = 3.0
prob211 = ODEProblem(model211!, u0, tend, ps211)

#---
@time sol211 = solve(prob211, Tsit5())

# Fig 2.11
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel="Time",
    ylabel="Concentration",
    title="Fig. 2.11\nFull model"
)
lines!(ax, 0 .. tend, t -> sol211(t).A, label="A")
lines!(ax, 0 .. tend, t -> sol211(t).B, label="B")
axislegend(ax, position=:rt)
fig

# ## Figure 2.12
# Rapid equilibrium assumption
_a212(u, p, t) = u.C * p.km1 / (p.km1 + p.k1)
_b212(u, p, t) = u.C * p.k1 / (p.km1 + p.k1)
function model212!(du, u, p, t)
    @unpack k0, k2 = p
    @unpack C = u
    A = _a212(u, p, t)
    B = C - A
    v0 = k0
    v2 = k2 * B
    du.C = v0 - v2
    nothing
end

#---
tend = 3.0
u0212 = CA.ComponentArray(C=sum(u0))
prob212 = ODEProblem(model212!, u0212, tend, ps211)

#---
@time sol212 = solve(prob212, Tsit5())
#---
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel="Time",
    ylabel="Concentration",
    title="Fig. 2.12\nRapid equilibrium model"
)

lines!(ax, 0 .. tend, t -> sol211(t).A, label="A (full solution)", linestyle=:dash)
lines!(ax, 0 .. tend, t -> sol211(t).B, label="B (full solution)", linestyle=:dash)
lines!(ax, 0 .. tend, t -> _a212(sol212(t), ps211, t), label="A (rapid equilibrium)")
lines!(ax, 0 .. tend, t -> _b212(sol212(t), ps211, t), label="B (rapid equilibrium)")
axislegend(ax, position=:rt)
fig

#===
## Figure 2.13

Rapid equilibrium (take 2)

When another set of parameters is not suitable for rapid equilibrium assumption.
===#

ps213 = CA.ComponentArray(k0=9.0, k1=20.0, km1=12.0, k2=2.0)
u0 = CA.ComponentArray(A=8.0, B=4.0)
tend = 3.0

@time sol213full = solve(ODEProblem(model211!, u0, tend, ps213), Tsit5())
@time sol213re = solve(ODEProblem(model212!, CA.ComponentArray(C=sum(u0)), tend, ps213), Tsit5())

#---
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel="Time",
    ylabel="Concentration",
    title="Fig. 2.13\nFull vs Rapid equilibrium"
)
lines!(ax, 0 .. tend, t -> sol213full(t).A, label="A (full solution)", linestyle=:dash)
lines!(ax, 0 .. tend, t -> sol213full(t).B, label="B (full solution)", linestyle=:dash)
lines!(ax, 0 .. tend, t -> _a212(sol213re(t), ps213, t), label="A (rapid equilibrium)")
lines!(ax, 0 .. tend, t -> _b212(sol213re(t), ps213, t), label="B (rapid equilibrium)")
axislegend(ax, position=:rt)
fig

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
u0214 = CA.ComponentArray(B=_u0214(sum(u0), ps214))

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

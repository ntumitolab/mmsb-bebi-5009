#===
# Fig 2.11, 2.12, 2.13, 2.14

Model reduction of ODE metabolic networks.
===#
using OrdinaryDiffEq
using ComponentArrays
using SimpleUnPack
using Plots
Plots.default(linewidth=2)

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
ps211 = ComponentArray(
    k0 = 0.0,
    k1 = 9.0,
    km1 = 12.0,
    k2 = 2.0
)

u0 = ComponentArray(
    A = 0.0,
    B = 10.0
)

tend = 3.0
prob211 = ODEProblem(model211!, u0, tend, ps211)

#---
@time sol211 = solve(prob211)

# Fig 2.11
plot(
    sol211,
    xlabel="Time (AU)",
    ylabel="Concentration (AU)",
    title="Fig. 2.11 (Full model)",
    labels=["A" "B"]
)

# ## Figure 2.12
# Rapid equilibrium assumption
_a212(u, p, t) =  u.C * p.km1 / (p.km1 + p.k1)
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
prob212 = ODEProblem(model212!, ComponentArray(C = sum(u0)), tend, ps211)

#---
@time sol212 = solve(prob212)

#---
fig = plot(sol211, line=(:dash, 1), label=["A (full solution)" "B (full solution)"])
ts = sol212.t
cs = sol212[1, :]
as = _a212.(sol212.u, Ref(ps212), ts)
bs = cs .- as
plot!(fig, ts, [as bs], label=["A (rapid equilibrium)" "B (rapid equilibrium)"])

plot!(fig,
    title="Fig. 2.12 (Rapid equilibrium model)",
    xlabel="Time (AU)",
    ylabel="Concentration (AU)"
)

fig

#===
## Figure 2.13

Rapid equilibrium (take 2)

When another set of parameters is not suitable for rapid equilibrium assumption.
===#

ps213 =ComponentArray(k0 = 9.0, k1 = 20.0, km1 = 12.0, k2 = 2.0)
u0 = ComponentArray(A = 8.0, B = 4.0)
tend = 3.0

@time sol213full = solve(ODEProblem(model211!, u0, tend, ps213))
@time sol213re = solve(ODEProblem(model212!, ComponentArray(C = sum(u0)), tend, ps213))

#---

fig = plot(sol213full, line=(:dash, 1), label=["A (full solution)" "B (full solution)"])
ts = sol213re.t
cs = sol213re[1, :]
as = _a212.(sol213re.u, Ref(ps213), ts)
bs = cs .- as
plot!(fig, ts, [as bs], label=["A (rapid equilibrium)" "B (rapid equilibrium)"])
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
u0214= ComponentArray(B = _u0214(sum(u0), ps214))

# Solve QSSA model
@time sol214 = solve(ODEProblem(model214!, u0214, tend, ps214))

fig = plot(sol213full, line=(:dash), label=["A (full solution)" "B (full solution)"])
ts = sol214.t
as = _a214.(sol214.u, Ref(ps214), ts)
bs = sol214[1, :]
plot!(fig, ts, [as bs], label=["A (QSSA)" "B (QSSA)"])
plot!(fig,
    xlabel="Time (AU)",
    ylabel="Concentration (AU)",
    title="Figure 2.14: Ref vs QSSA",
    xlims=(0.0, tend)
)

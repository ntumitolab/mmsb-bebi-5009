#===
# Fig 3.03

Michaelis-Menten kinetics
===#
using OrdinaryDiffEq
using ComponentArrays
using SimpleUnPack
using Plots
Plots.default(linewidth=2)

# Enzyme kinetics full model
_e303(u, p, t) = p.ET - u.ES
function model303!(du, u, p, t)
    @unpack k1, km1, k2 = p
    @unpack S, ES, P = u
    E = _e303(u, p, t)
    v1 = k1 * S * E - km1 * ES
    v2 = k2 * ES
    du.S = -v1
    du.ES = v1 - v2
    du.P = v2
    nothing
end

#---
ps303 = ComponentArray(
    k1 = 30.0,
    km1 = 1.0,
    k2 = 10.0,
    ET = 1.0
)
u0 = ComponentArray(
    S = 5.0,
    ES = 0.0,
    P = 0.0
)
tend = 1.0
prob303 = ODEProblem(model303!, u0, tend, ps303)

#---
@time sol = solve(prob303, Tsit5())

#---
pl303 = plot(sol, xlabel="Time (AU)", ylabel="Concentration (AU)", legend=:right, title="Fig 3.03", labels=["S (full)" "ES (full)" "P (full)"])
e303_t = t -> _e303(sol(t), ps303, t)
plot!(pl303, e303_t, 0, tend, label="E (full)")

# ## QSSA of ES complex
_s303mm(u, p, t) = p.S0 - u.P
function model303mm!(du, u, p, t)
    @unpack k1, km1, k2, ET = p
    @unpack P = u
    S = _s303mm(u, p, t)
    ES = (k1 * S * ET) / (km1 + k2 + k1 * S)
    v2 = k2 * ES
    du.P = v2
    nothing
end

ps303mm = ComponentArray(ps303; S0=5.0)
prob303mm = ODEProblem(model303mm!, ComponentArray(P=0.0), tend, ps303mm)

#---
@time sol303mm = solve(prob303mm, Tsit5())

#---
plot(sol, idxs=[1, 3], label=["S (full)" "P (full)"], line=(:dash))
p303mm_t = t -> sol303mm(t).P
s303mm_t = t -> _s303mm(sol303mm(t), ps303mm, t)

plot!([s303mm_t p303mm_t], 0, tend, label=["S (MM)" "P (MM)"], title="Fig. 3.03", xlabel="Time (AU)", ylabel="Concentration (AU)", xlims=(0., tend), ylims=(0., 5.), legend=:right)

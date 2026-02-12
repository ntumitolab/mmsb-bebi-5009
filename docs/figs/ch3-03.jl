#===
# Fig 3.03

Michaelis-Menten kinetics
===#
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots

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
    k1=30.0,
    km1=1.0,
    k2=10.0,
    ET=1.0
)
u0 = ComponentArray(
    S=5.0,
    ES=0.0,
    P=0.0
)
tend = 1.0
prob303 = ODEProblem(model303!, u0, tend, ps303)

#---
@time sol = solve(prob303, Tsit5())

#---
fig303 = Figure()
ax303 = Axis(
    fig303[1, 1],
    xlabel="Time",
    ylabel="Concentration",
    title="Fig. 3.03 (Full model)"
)
lines!(ax303, 0 .. tend, t -> sol(t).S, label="S")
lines!(ax303, 0 .. tend, t -> sol(t).ES, label="ES")
lines!(ax303, 0 .. tend, t -> sol(t).P, label="P")
lines!(ax303, 0 .. tend, t -> _e303(sol(t), ps303, t), label="E")
axislegend(ax303, position=:rt)

fig303

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
fig303mm = Figure()
ax303mm = Axis(
    fig303mm[1, 1],
    xlabel="Time",
    ylabel="Concentration",
    title="Fig. 3.03 (QSSA)"
)
lines!(ax303mm, 0 .. tend, t -> sol(t).S, label="S (full)", linestyle=:dash)
lines!(ax303mm, 0 .. tend, t -> sol(t).P, label="P (full)", linestyle=:dash)
lines!(ax303mm, 0 .. tend, t -> _s303mm(sol303mm(t), ps303mm, t), label="S (MM)")
lines!(ax303mm, 0 .. tend, t -> sol303mm(t).P, label="P (MM)")
axislegend(ax303mm, position=:rc)
fig303mm

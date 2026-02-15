# # Fig 3.03
# Michaelis-Menten kinetics
using OrdinaryDiffEq
using ComponentArrays: ComponentArray
using SimpleUnPack
using CairoMakie

# Enzyme kinetics full model
function model303(u, p, t)
    @unpack k1, km1, k2 = p
    @unpack S, ES, P = u
    E = p.ET - ES
    v1 = k1 * S * E - km1 * ES
    v2 = k2 * ES
    return (; dS=-v1, dES=v1-v2, dP=v2, E)
end

function model303!(du, u, p, t)
    @unpack dS, dES, dP = model303(u, p, t)
    du.S = dS
    du.ES = dES
    du.P = dP
    nothing
end

#---
ps303 = ComponentArray(k1=30.0, km1=1.0, k2=10.0, ET=1.0)
u0 = ComponentArray(S=5.0, ES=0.0, P=0.0)
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
lines!(ax303, 0 .. tend, t -> model303(sol(t), ps303, t).E, label="E")
axislegend(ax303, position=:rc)
fig303

# ## QSSA of ES complex
function model303mm(u, p, t)
    @unpack k1, km1, k2, ET, S0 = p
    @unpack P = u
    S = S0 - P
    ES = (k1 * S * ET) / (km1 + k2 + k1 * S)
    E = ET - ES
    return (; dP=k2 * ES, S, E, ES)
end
function model303mm!(du, u, p, t)
    @unpack dP = model303mm(u, p, t)
    du.P = dP
    nothing
end

ps303mm = ComponentArray(ps303; S0=5.0)
u0303mm = ComponentArray(P=0.0)
prob303mm = ODEProblem(model303mm!, u0303mm, tend, ps303mm)

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
lines!(ax303mm, 0 .. tend, t -> model303mm(sol303mm(t), ps303mm, t).S, label="S (MM)")
lines!(ax303mm, 0 .. tend, t -> sol303mm(t).P, label="P (MM)")
axislegend(ax303mm, position=:rc)
fig303mm

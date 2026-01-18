#===
# Fig 6.16

Model of apoptosis signalling pathway
===#
using ComponentArrays: ComponentArray
using SimpleUnPack
using OrdinaryDiffEq
using DiffEqCallbacks
using CairoMakie

#---
function model616!(D, u, p, t)
    @unpack k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, I = p
    @unpack C8, C8s, C3, C3s, BAR, IAP, C8sBAR, C3sIAP = u
    v1 = k1 - k2 * C8
    v2 = k3 * (C3s + I) * C8
    v3 = k4 * C8s
    v4 = k5 * C8s * BAR - k6 * C8sBAR
    v5 = k7 - k8 * C3
    v6 = k9 * C8s * C3
    v7 = k10 * C3s
    v8 = k11 * C3s * IAP - k12 * C3sIAP
    v9 = k13 - k14 * BAR
    v10 = k15 - (k16 + k17 * C3s) * IAP
    v11 = k18 * C8sBAR
    v12 = k19 * C3sIAP
    D.C8 = v1 - v2
    D.C8s = v2 - v3 - v4
    D.C3 = v5 - v6
    D.C3s = v6 - v7 - v8
    D.BAR = v9 - v4
    D.IAP = v10 - v8
    D.C8sBAR = v4 - v11
    D.C3sIAP = v8 - v12
    nothing
end

ps616 = ComponentArray(
    k1 = 507.0,
    k2 = 3.9e-3,
    k3 = 1e-5,
    k4 = 5.8e-3,
    k5 = 5e-4,
    k6 = 0.21,
    k7 = 81.9,
    k8 = 3.9e-3,
    k9 = 5.8e-6,
    k10 = 5.8e-3,
    k11 = 5e-4,
    k12 = 0.21,
    k13 = 40.0,
    k14 = 1e-3,
    k15 = 464.0,
    k16 = 1.16e-2,
    k17 = 3e-4,
    k18 = 1.16e-2,
    k19 = 1.73e-2,
    I = 0.0
)

u0616 = ComponentArray(
    C8 = 1.3E5,
    C8s = 0.0,
    C3 = 0.21E5,
    C3s = 0.0,
    BAR = 0.4E5,
    IAP = 0.4E5,
    C8sBAR = 0.0,
    C3sIAP = 0.0
)

# Event: increase I at t=100, decrease I at t=1200
affect_i1!(integrator) = integrator.p.I = 200.0
affect_i2!(integrator) = integrator.p.I = 0.0
event_i1 = PresetTimeCallback([100.0], affect_i1!)
event_i2 = PresetTimeCallback([1200.0], affect_i2!)
cbs = CallbackSet(event_i1, event_i2)
tend = 1800.0
prob = ODEProblem(model616!, u0616, tend, ps616)

#---
@time sol = solve(prob, TRBDF2(), callback=cbs)

#---
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time", ylabel="Concentration", title="Fig 6.16\nApoptosis signalling pathway")
lines!(ax, 0..tend, t -> sol(t).C8s, label="C8s")
lines!(ax, 0..tend, t -> sol(t).C3s, label="C3s")
lines!(ax, 0..tend, t -> 100 * 200 * (100<=t<1200), label="I (×100)")
axislegend(ax, position=:rc)
fig

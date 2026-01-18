#===
# Fig 6.14

Model of E. coli chemotaxis signalling pathway
===#
using ComponentArrays: ComponentArray
using SimpleUnPack
using OrdinaryDiffEq
using DiffEqCallbacks
using CairoMakie

#---
function model614!(D, u, p, t)
    hil(x, k) = x / (k + x)
    @unpack k1, k2, k3, km1, km2, km3, k4, km4, k5, km5, KM1, KM2, L, R, Atotal, Btotal = p
    @unpack Am, AmL, AL, BP = u
    A = Atotal - Am - AmL - AL
    B = Btotal - BP
    v1 = k1 * BP * hil(Am, KM1)
    v2 = k2 * BP * hil(AmL, KM2)
    v3 = km1 * R
    v4 = km2 * R
    v5 = k3 * L * Am - km3 * AmL
    v6 = k4 * L * A - km4 * AL
    v7 = k5 * Am * B - km5 * BP
    D.Am = v3 - v1 - v5
    D.AmL = v4 - v2 + v5
    D.AL = v2 + v6
    D.BP = v7
    nothing
end

#---
ps614 = ComponentArray(
    k1 = 200.0,
    k2 = 1.0,
    k3 = 1.0,
    km1 = 1.0,
    km2 = 1.0,
    km3 = 1.0,
    k4 = 1.0,
    km4 = 1.0,
    k5 = 0.05,
    km5 = 0.005,
    KM1 = 1.0,
    KM2 = 1.0,
    L = 20.0,
    R = 1.0,
    Atotal = 2.0,
    Btotal = 1.0
)

u0614 = ComponentArray(
    Am = 0.0360,
    AmL = 1.5593,
    AL = 0.3504,
    BP = 0.2644
)

# Events to increase L
affect_L1!(integrator) = integrator.p.L = 40.0
affect_L2!(integrator) = integrator.p.L = 80.0
event_L1 = PresetTimeCallback([10.0], affect_L1!)
event_L2 = PresetTimeCallback([30.0], affect_L2!)
cbs = CallbackSet(event_L1, event_L2)
tend = 50.0
prob614 = ODEProblem(model614!, u0614, tend, ps614)

#---
@time sol614 = solve(prob614, Tsit5(), callback=cbs)

#---
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time", ylabel="Concentration", title="Fig 6.14\nE. coli chemotaxis signalling pathway")
lines!(ax, 0..tend, t -> sol614(t).Am, color=:blue, label="Am")
axislegend(ax, position=:rc)
fig

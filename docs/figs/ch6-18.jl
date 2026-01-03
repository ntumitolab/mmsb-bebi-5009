#===
# Fig 6.18

Model of calcium-induced calcium release in hepatocytes
===#
using ComponentArrays
using SimpleUnPack
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
Plots.default(linewidth=2)

#---
hil(x, k) = x / (k + x)
hil(x, k, n) = hil(x^n, k^n)
_r618(u, p, t) = p.Rtot - u.RI - u.RIC - u.RICC
function model618!(D, u, p, t)
    @unpack k1, km1, k2, km2, k3, km3, vr, γ0, γ1, p1, p2, Cer, I, Rtot = p
    @unpack RI, RIC, RICC, C = u
    R = _r618(u, p, t)
    v1 = k1 * I * R - km1 * RI
    v2 = k2 * C * RI - km2 * RIC
    v3 = k3 * C * RIC - km3 * RICC
    v4 = vr * (γ0 + γ1 * RIC) * (Cer - C) - p1 * hil(C, p2, 4)
    D.RI = v1 - v2
    D.RIC = v2 - v3
    D.RICC = v3
    D.C = v4
    nothing
end

#---
ps618 = ComponentArray(
    k1 = 12.0,
    km1 = 8.0,
    k2 = 15.0,
    km2 = 1.65,
    k3 = 1.8,
    km3 = 0.21,
    vr = 0.185,
    γ0 = 0.1,
    γ1 = 20.5,
    p1 = 8.5,
    p2 = 0.065,
    Cer = 8.37,
    I = 0.0,
    Rtot = 1.0
)

u0618 = ComponentArray(
    C = 0.0,
    RIC = 0.0,
    RICC = 0.0,
    RI = 0.0
)

# ## Fig 6.18 (A)
ps618a = ComponentArray(ps618; I = 2.0)
tend = 25.0
prob618a = ODEProblem(model618!, u0618, tend, ps618a)
@time sol618a = solve(prob618a)

#---
plot(sol618a, idxs=[2, 3, 1], title="Fig 6.18 (A)", xlabel="Time", ylabel="Abundance", legend=:topright, labels=["RIC" "RICC" "C"])

# ## Fig 6.18 (B)
affect_i1!(integrator) = integrator.p.I = 0.7
affect_i2!(integrator) = integrator.p.I = 1.2
affect_i3!(integrator) = integrator.p.I = 4.0
event_1 = PresetTimeCallback(20.0, affect_i1!)
event_2 = PresetTimeCallback(60.0, affect_i2!)
event_3 = PresetTimeCallback(90.0, affect_i3!)
cbs = CallbackSet(event_1, event_2, event_3)

tend = 120.
prob618b = ODEProblem(model618!, u0618, tend, ps618)

#---
@time sol618b = solve(prob618b, callback=cbs)

#---
plot(sol618b, idxs=[1], title="Fig 6.18 (B)", xlabel="Time", ylabel="Ca concentration", legend=false, ylim=(0, 2.5))

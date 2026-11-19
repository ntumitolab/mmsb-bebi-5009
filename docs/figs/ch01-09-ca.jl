# # Fig 1.9
# Hodgkin-Huxley model for action potentials in neurons.
using ComponentArrays
using SimpleUnPack
using DiffEqCallbacks
using OrdinaryDiffEq
using Plots
Plots.gr(linewidth=1.5)

#---
exprel(x) = x / expm1(x)

## The algebraic variables of the MTK version (`iNa`, `iK`, `iLeak`, and the
## rate functions) are only intermediates, so they are computed locally.
## `iStim` is a discrete parameter, like `@discretes iStim(t)` in MTK.
function hh_ode!(du, u, p, t)
    @unpack v, m, h, n = u
    ma = exprel(-0.10 * (v + 35))
    mb = 4.0 * exp(-(v + 60) / 18.0)
    ha = 0.07 * exp(-(v + 60) / 20)
    hb = 1 / (exp(-(v + 30) / 10) + 1)
    na = 0.1 * exprel(-0.1 * (v + 50))
    nb = 0.125 * exp(-(v + 60) / 80)
    iNa = p.G_N_BAR * (v - p.E_N) * (m^3) * h
    iK = p.G_K_BAR * (v - p.E_K) * (n^4)
    iLeak = p.G_LEAK * (v - p.E_LEAK)
    du.v = -(iNa + iK + iLeak + p.iStim) / p.C_M
    du.m = -(ma + mb) * m + ma
    du.h = -(ha + hb) * h + ha
    du.n = -(na + nb) * n + na
    return nothing
end

## Electrical stimulation events, replacing the `SymbolicDiscreteCallback`s.
## `PresetTimeCallback` sets `p.iStim` at each event time and makes the
## solver end its step there, so the current is piecewise constant between
## events.
affect_stim_1_on!(integrator) = integrator.p.iStim = -6.5
affect_stim_2_on!(integrator) = integrator.p.iStim = -7.5
affect_stim_off!(integrator) = integrator.p.iStim = 0.0
cbs = CallbackSet(
    PresetTimeCallback(20.0, affect_stim_1_on!),
    PresetTimeCallback(21.0, affect_stim_off!),
    PresetTimeCallback(60.0, affect_stim_2_on!),
    PresetTimeCallback(61.0, affect_stim_off!),
)

#---
p = ComponentArray(E_N=55.0, E_K=-72.0, E_LEAK=-49.0, G_N_BAR=120.0, G_K_BAR=36.0, G_LEAK=0.3, C_M=1.0, iStim=0.0)
u0 = ComponentArray(v=-59.8977, m=0.0536, h=0.5925, n=0.3192)

tend = 100.0
@time "Build problem" prob = ODEProblem(hh_ode!, u0, (0.0, tend), p)
@time "Solve problem" sol = solve(prob, FBDF(), callback=cbs)

# `idxs=[:v]` selects the membrane potential by name.
plot(sol, idxs=[:v], xlabel="Time (ms)", ylabel="Membrane potential (mV)", title="Fig 1.9")

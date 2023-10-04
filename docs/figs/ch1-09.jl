#===

# Fig 1.09
Hodgkin-Huxley model

===#

using DifferentialEquations
using LabelledArrays
using SimpleUnPack
using Plots
Plots.default(linewidth=2)

# Convenience functions
hil(x, k) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)
exprel(x) = x / expm1(x)

# Stimulation current
_istim(t) = ifelse(20 <= t <= 21, -6.6, 0.0) + ifelse(60 <= t <= 61, -6.9, 0.0)

# HH Neuron model
function hh!(D, u, p, t)
    @unpack G_N_BAR, E_N, G_K_BAR, E_K, G_LEAK, E_LEAK, C_M = p
    @unpack v, m, h, n = u
    mα = exprel(-0.10 * (v + 35))
    mβ  = 4.0 * exp(-(v + 60) / 18.0)
    hα  = 0.07 * exp(- ( v + 60) / 20)
    hβ  = 1 / (exp(-(v+30)/10) + 1)
    nα  = 0.1 * exprel(-0.1 * (v+50))
    nβ  = 0.125 * exp( -(v+60) / 80)
    iNa = G_N_BAR * (v - E_N) * (m^3) * h
    iK  = G_K_BAR * (v - E_K) * (n^4)
    iLeak = G_LEAK * (v - E_LEAK)
    iStim = _istim(t)
    D.v = -(iNa + iK + iLeak + iStim) / C_M
    D.m = -(mα + mβ) * m + mα
    D.h = -(hα + hβ) * h + hα
    D.n = -(nα + nβ) * n + nα
    return nothing
end

#---
ps = (
    E_N = 55.0,       ## Reversal potential of Na (mV)
    E_K = -72.0,      ## Reversal potential of K (mV)
    E_LEAK = -49.0,   ## Reversal potential of leaky channels (mV)
    G_N_BAR = 120.0,  ## Max. Na channel conductance (mS/cm^2)
    G_K_BAR = 36.0,   ## Max. K channel conductance (mS/cm^2)
    G_LEAK = 0.30,    ## Max. leak channel conductance (mS/cm^2)
    C_M = 1.0         ## membrane capacitance (uF/cm^2))
)
u0 = LVector(v=-59.8977, m=0.0536, h=0.5925, n=0.3192)
tend = 100.0

#---
prob = ODEProblem(hh!, u0, tend, ps)

#---
sol = solve(prob, tstops=[20., 60.])

#---
p1 = plot(sol, idxs=[:v], ylabel="Membrane potential (mV)", xlabel="", legend=false, title="Fig 1.9")
p2 = plot(sol, idxs = [:m, :h, :n], xlabel="")
p3 = plot(_istim, sol.t, xlabel="Time (ms)", ylabel="Current", labels="Stimulation current")
fig = plot(p1, p2, p3, layout=(3, 1), size=(600, 900), leftmargin=5*Plots.mm)

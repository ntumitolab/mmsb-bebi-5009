#===
# Fig 1.09

Hodgkin-Huxley model
===#
using OrdinaryDiffEq
using ComponentArrays
using DiffEqCallbacks
using SimpleUnPack
using Plots
Plots.default(linewidth=2)

# Convenience functions
hil(x, k) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)
exprel(x) = x / expm1(x)
_mα(u, p, t) = exprel(-0.10 * (u.v + 35))
_mβ(u, p, t) = 4.0 * exp(-(u.v + 60) / 18.0)
_hα(u, p, t) = 0.07 * exp(-(u.v + 60) / 20)
_hβ(u, p, t) = 1 / (exp(-(u.v + 30) / 10) + 1)
_nα(u, p, t) = 0.1 * exprel(-0.1 * (u.v + 50))
_nβ(u, p, t) = 0.125 * exp(-(u.v + 60) / 80)
_iNa(u, p, t) = p.G_N_BAR * (u.v - p.E_N) * (u.m^3) * u.h
_iK(u, p, t) = p.G_K_BAR * (u.v - p.E_K) * (u.n^4)
_iLeak(u, p, t) = p.G_LEAK * (u.v - p.E_LEAK)

# HH Neuron model
function hh_neuron!(du, u, p, t)
    @unpack E_N, E_K, E_LEAK, G_N_BAR, G_K_BAR, G_LEAK, C_M, iStim = p
    @unpack v, m, h, n = u
    mα = _mα(u, p, t)
    mβ = _mβ(u, p, t)
    hα = _hα(u, p, t)
    hβ = _hβ(u, p, t)
    nα = _nα(u, p, t)
    nβ = _nβ(u, p, t)
    iNa = _iNa(u, p, t)
    iK = _iK(u, p, t)
    iLeak = _iLeak(u, p, t)
    du.v = -(iNa + iK + iLeak + iStim) / C_M
    du.m = -(mα + mβ) * m + mα
    du.h = -(hα + hβ) * h + hα
    du.n = -(nα + nβ) * n + nα
    nothing
end

# Problem definition
ps = ComponentArray(
    E_N = 55.0,       ## Reversal potential of Na (mV)
    E_K = -72.0,      ## Reversal potential of K (mV)
    E_LEAK = -49.0,   ## Reversal potential of leaky channels (mV)
    G_N_BAR = 120.0,  ## Max. Na channel conductance (mS/cm^2)
    G_K_BAR = 36.0,   ## Max. K channel conductance (mS/cm^2)
    G_LEAK = 0.30,    ## Max. leak channel conductance (mS/cm^2)
    C_M = 1.0,        ## membrane capacitance (uF/cm^2))
    iStim = 0.0       ## stimulation current
)

u0 = ComponentArray(
    v = -59.8977,
    m = 0.0536,
    h = 0.5925,
    n = 0.3192,
)

tspan = (0.0, 100.0)

prob = ODEProblem(hh_neuron!, u0, tspan, ps)

# Callbacks for external current and solve the problem
affect_stim_on1!(integrator) = integrator.p.iStim = -6.6
affect_stim_off1!(integrator) = integrator.p.iStim = 0.0
affect_stim_on2!(integrator) = integrator.p.iStim = -6.9
affect_stim_off2!(integrator) = integrator.p.iStim = 0.0
cb_stim_on1 = PresetTimeCallback(20.0, affect_stim_on1!)
cb_stim_off1 = PresetTimeCallback(21.0, affect_stim_off1!)
cb_stim_on2 = PresetTimeCallback(60.0, affect_stim_on2!)
cb_stim_off2 = PresetTimeCallback(61.0, affect_stim_off2!)
cbs = CallbackSet(cb_stim_on1, cb_stim_off1, cb_stim_on2, cb_stim_off2)
@time sol = solve(prob, Tsit5(), callback=cbs)

# Visualization
plot(t->sol(t).v, 0, 100, xlabel="Time (ms)", ylabel="Membrane potential (mV)", title="Fig 1.9", label=false, legend=:topleft)

#===
# Fig 1.09

Hodgkin-Huxley model
===#
using OrdinaryDiffEq
using ComponentArrays: ComponentArray
using DiffEqCallbacks
using SimpleUnPack
using CairoMakie

# HH Neuron model
function hh_neuron(u, p, t)
    exprel(x) = x / expm1(x)
    @unpack E_N, E_K, E_LEAK, G_N_BAR, G_K_BAR, G_LEAK, C_M, iStim = p
    @unpack v, m, h, n = u
    ma = exprel(-0.10 * (v + 35))
    mb = 4.0 * exp(-(v + 60) / 18.0)
    ha = 0.07 * exp(-(v + 60) / 20)
    hb = 1 / (exp(-(v + 30) / 10) + 1)
    na = 0.1 * exprel(-0.1 * (v + 50))
    nb = 0.125 * exp(-(v + 60) / 80)
    iNa = G_N_BAR * (v - E_N) * (m^3) * h
    iK = G_K_BAR * (v - E_K) * (n^4)
    iLeak = G_LEAK * (v - E_LEAK)
    return (;
        dv = -(iNa + iK + iLeak + iStim) / C_M,
        dm = -(ma + mb) * m + ma,
        dh = -(ha + hb) * h + ha,
        dn = -(na + nb) * n + na,
        ma, mb, ha, hb, na, nb, iNa, iK, iLeak
    )
end

# Inplace version of the HH neuron model
function hh_neuron!(du, u, p, t)
    @unpack dv, dm, dh, dn = hh_neuron(u, p, t)
    du.v = dv
    du.m = dm
    du.h = dh
    du.n = dn
    nothing
end

# Problem definition
ps = ComponentArray(
    E_N=55.0,       ## Reversal potential of Na (mV)
    E_K=-72.0,      ## Reversal potential of K (mV)
    E_LEAK=-49.0,   ## Reversal potential of leaky channels (mV)
    G_N_BAR=120.0,  ## Max. Na channel conductance (mS/cm^2)
    G_K_BAR=36.0,   ## Max. K channel conductance (mS/cm^2)
    G_LEAK=0.30,    ## Max. leak channel conductance (mS/cm^2)
    C_M=1.0,        ## membrane capacitance (uF/cm^2))
    iStim=0.0       ## stimulation current
)

u0 = ComponentArray(
    v=-59.8977,
    m=0.0536,
    h=0.5925,
    n=0.3192,
)

tend = 100.0
prob = ODEProblem(hh_neuron!, u0, tend, ps)

# Callbacks for external current and solve the problem
cbs = let
    affect_stim_on1!(integrator) = integrator.p.iStim = -6.6
    affect_stim_off1!(integrator) = integrator.p.iStim = 0.0
    affect_stim_on2!(integrator) = integrator.p.iStim = -6.9
    affect_stim_off2!(integrator) = integrator.p.iStim = 0.0
    cb_stim_on1 = PresetTimeCallback(20.0, affect_stim_on1!)
    cb_stim_off1 = PresetTimeCallback(21.0, affect_stim_off1!)
    cb_stim_on2 = PresetTimeCallback(60.0, affect_stim_on2!)
    cb_stim_off2 = PresetTimeCallback(61.0, affect_stim_off2!)
    cbs = CallbackSet(cb_stim_on1, cb_stim_off1, cb_stim_on2, cb_stim_off2)
end

@time sol = solve(prob, Tsit5(), callback=cbs)

# Visualization
fig, ax, sp = lines(0 .. tend, t -> sol(t).v, axis=(title="Fig 1.9\nHodgkin-Huxley Neuron", xlabel="Time (ms)", ylabel="Membrane potential (mV)"))

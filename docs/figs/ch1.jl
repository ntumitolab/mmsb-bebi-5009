# # Chapter 1
# ## Fig 1.7
# For Figures 1.7, 7.13, 7.14, 7.15.
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using Plots
Plots.gr(linewidth=1.5)

# Collins toggle switch model
function collins_sys(; name=:collins, simplify=true)
    @parameters a1=3 a2=2.5 β=4 γ=4
    @variables i1(t) i2(t) s1(t)=0.075 s2(t)=2.5
    hil(x, k) = x / (x + k)
    hil(x, k, n) = hil(x^n, k^n)
    eqs = [
        D(s1) ~ a1 * hil(1 + i2, s2, β) - s1,
        D(s2) ~ a2 * hil(1 + i1, s1, γ) - s2,
        i1 ~ 10 * (t > 30) * (t < 40),
        i2 ~ 10 * (t > 10) * (t < 20)
    ]
    sys = ODESystem(eqs, t; name)
    return simplify ? mtkcompile(sys) : sys
end

tend = 50.0
@time "Build system" sys = collins_sys()
@time "Build problem" prob = ODEProblem(sys, [], tend)
@time sol = solve(prob, Tsit5(), tstops=[10.0, 20.0, 30.0, 40.0])

# Visualization
plot(sol, title="Fig. 1.7", xlabel="Time", ylabel="Concentration")

#---
plot(sol, idxs=[sys.i1, sys.i2], labels=["i1" "i2"], xlabel="Time", ylabel="Signal strength")

# ## Fig 1.09
# Hodgkin-Huxley model
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
Plots.gr(linewidth=1.5)

# HH Neuron model
function hh_sys(; name=:hh, simplify=true)
    exprel(x) = x / expm1(x)
    @discretes iStim(t)=0.0
    @parameters E_N=55 E_K=-72 E_LEAK=-49 G_N_BAR=120 G_K_BAR=36 G_LEAK=0.3 C_M=1
    @variables v(t)=-59.8977 m(t)=0.0536 h(t)=0.5925 n(t)=0.3192 iNa(t) iK(t) iLeak(t) ma(t) mb(t) ha(t) hb(t) na(t) nb(t)

    ## Electrical stimulation events
    stim_on_1 = ModelingToolkit.SymbolicDiscreteCallback([20.0] => [iStim ~ -6.6], discrete_parameters = iStim, iv = t)
    stim_off_1 = ModelingToolkit.SymbolicDiscreteCallback([21.0] => [iStim ~ 0.0], discrete_parameters = iStim, iv = t)
    stim_on_2 = ModelingToolkit.SymbolicDiscreteCallback([60.0] => [iStim ~ -6.9], discrete_parameters = iStim, iv = t)
    stim_off_2 = ModelingToolkit.SymbolicDiscreteCallback([61.0] => [iStim ~ 0.0], discrete_parameters = iStim, iv = t)

    eqs = [
        ma ~ exprel(-0.10 * (v + 35)),
        mb ~ 4.0 * exp(-(v + 60) / 18.0),
        ha ~ 0.07 * exp(-(v + 60) / 20),
        hb ~ 1 / (exp(-(v + 30) / 10) + 1),
        na ~ 0.1 * exprel(-0.1 * (v + 50)),
        nb ~ 0.125 * exp(-(v + 60) / 80),
        iNa ~ G_N_BAR * (v - E_N) * (m^3) * h,
        iK ~ G_K_BAR * (v - E_K) * (n^4),
        iLeak ~ G_LEAK * (v - E_LEAK),
        D(v) ~ -(iNa + iK + iLeak + iStim) / C_M,
        D(m) ~ -(ma + mb) * m + ma,
        D(h) ~ -(ha + hb) * h + ha,
        D(n) ~ -(na + nb) * n + na
    ]
    sys = ODESystem(eqs, t; name, discrete_events=[stim_on_1, stim_off_1, stim_on_2, stim_off_2])
    return simplify ? mtkcompile(sys) : sys
end

tend = 100.0
@time "Build system" sys = hh_sys()
@time "Build problem" prob = ODEProblem(sys, [], tend)
@time "Solve problem" sol = solve(prob, FBDF())

#---
plot(sol, idxs=[sys.v], xlabel="Time (ms)", ylabel="Membrane potential (mV)", title="Fig 1.9")

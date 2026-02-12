#===
# Fig 1.09

Hodgkin-Huxley model
===#
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
Plots.pythonplot()

#---
function prob109(; tend = 100.0)
    ## Convenience functions
    hil(x, k) = x / (x + k)
    hil(x, k, n) = hil(x^n, k^n)
    exprel(x) = x / expm1(x)
    @parameters E_N=55 E_K=-72 E_LEAK=-49 G_N_BAR=120 G_K_BAR=36 G_LEAK=0.30 C_M=1
    @variables v(t)=-59.8977 m(t)=0.0536 h(t)=0.5925 n(t)=0.3192
    @variables mα(t) mβ(t) hα(t) hβ(t) nα(t) nβ(t) iNa(t) iK(t) iLeak(t) iStim(t)
    eqs = [
        mα ~ exprel(-0.10 * (v + 35))
        mβ ~ 4.0 * exp(-(v + 60) / 18)
        hα ~ 0.07 * exp(-(v + 60) / 20)
        hβ ~ 1 / (exp(-(v + 30) / 10) + 1)
        nα ~ 0.1 * exprel(-0.1 * (v + 50))
        nβ ~ 0.125 * exp(-(v + 60) / 80)
        iNa ~ G_N_BAR * (v - E_N) * (m^3) * h
        iK ~ G_K_BAR * (v - E_K) * (n^4)
        iLeak ~ G_LEAK * (v - E_LEAK)
        iStim ~ -6.6 * (20 < t) * (t < 21) - 8.0 * (60 < t) * (t < 61)
        D(v) ~ -(iNa + iK + iLeak + iStim) / C_M
        D(m) ~ -(mα + mβ) * m + mα
        D(h) ~ -(hα + hβ) * h + hα
        D(n) ~ -(nα + nβ) * n + nα
    ]
    @mtkcompile sys = ODESystem(eqs, t)
    return ODEProblem(sys, [], tend)
end

@time prob = prob109()

#---
@time sol = solve(prob, KenCarp47(), tstops=[20, 21, 60, 61])

# Visualization
plot(sol, idxs= prob.f.sys.v, xlabel="Time", ylabel="Voltage (mV)", title="Fig. 1.9\nHodgkin-Huxley Model")

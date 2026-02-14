#===
# Fig 1.7

For Figures 1.7, 7.13, 7.14, 7.15.
===#
using OrdinaryDiffEq
using DiffEqCallbacks
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
Plots.gr(framestyle = :box, linewidth=1.5)

#---
function prob107(; tend = 50.0)
    hil(x, k) = x / (x + k)
    hil(x, k, n) = hil(x^n, k^n)
    @parameters a1=3.0 a2=2.5 β=4.0 γ=4.0
    @variables s1(t)=0.075 s2(t)=2.5 i1(t) i2(t)
    eqs = [
        D(s1) ~ a1 * hil(1 + i2, s2, β) - s1
        D(s2) ~ a2 * hil(1 + i1, s1, γ) - s2
        i2 ~ 10 * (10 < t) * (t < 20)
        i1 ~ 10 * (30 < t) * (t < 40)
    ]
    @mtkcompile sys = ODESystem(eqs, t)
    return ODEProblem(sys, [], tend)
end

@time prob = prob107()

#---
@time sol = solve(prob, Tsit5(); tstops=10:10:40)

#---
plot(sol, xlabel="Time", ylabel="Concentration", title="Fig. 1.7")

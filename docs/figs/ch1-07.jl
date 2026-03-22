# # Fig 1.7
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

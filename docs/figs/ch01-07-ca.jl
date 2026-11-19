# # Fig 1.7 A
# Collins toggle switch model for Figures 1.7, 7.13, 7.14, and 7.15.
using ComponentArrays
using DiffEqCallbacks
using OrdinaryDiffEq
using Plots
using Startup: hil
Plots.gr(linewidth=1.5)

#---
## External input signals, replacing the algebraic equations
## `i1 ~ 10 * (t > 30) * (t < 40)` and `i2 ~ 10 * (t > 10) * (t < 20)`.
## A signal is ON on the half-open interval [t_on, t_off).
signal(t, t_on, t_off) = 10 * (t_on <= t < t_off)
i1(t) = signal(t, 30, 40)
i2(t) = signal(t, 10, 20)

#---
## `u` and `du` are ComponentArrays with components `s1` and `s2`, and the
## input signals `i1` and `i2` are kept in the parameter array `p`, so the
## equations read just like the MTK version.
function collins_ode!(du, u, p, t)
    du.s1 = p.a1 * hil(1 + p.i2, u.s2, p.β) - u.s1
    du.s2 = p.a2 * hil(1 + p.i1, u.s1, p.γ) - u.s2
    return nothing
end

## The signals jump at t = 10, 20, 30, and 40. `PresetTimeCallback` updates
## `p.i1` and `p.i2` at each jump time and splits the integration there.
## Passing `tstops` to `solve` is not enough for implicit methods like `FBDF`:
## the signal would be sampled at the step endpoints only, where it is zero.
update_signals!(integrator) = (integrator.p.i1, integrator.p.i2) = (i1(integrator.t), i2(integrator.t))
cbs = PresetTimeCallback([10.0, 20.0, 30.0, 40.0], update_signals!)

#---
p = ComponentArray(a1=3.0, a2=2.5, β=4.0, γ=4.0, i1=0.0, i2=0.0)
u0 = ComponentArray(s1=0.075, s2=2.5)

tend = 50.0
@time "Build problem" prob = ODEProblem(collins_ode!, u0, (0.0, tend), p)
@time "Solve problem" sol = solve(prob, FBDF(), callback=cbs)

#---
## Component names work as indices into the solution too:
## `idxs=[:s1, :s2]` selects the state variables by name, without knowing
## their order in `u`.
plot(sol, idxs=[:s1, :s2], labels=["s1" "s2"], title="Fig. 1.7", xlabel="Time", ylabel="Concentration")

#---
## Signal strength. The callback saves the solution at the jump times, so
## `steppost` draws exact rectangles.
plot(sol.t, [i1.(sol.t) i2.(sol.t)], st=:steppost, labels=["i1" "i2"], xlabel="Time", ylabel="Signal strength")

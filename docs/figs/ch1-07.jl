#===
# Fig 1.7

For Figures 1.7, 7.13, 7.14, 7.15.
===#
using OrdinaryDiffEq
using DiffEqCallbacks
using ComponentArrays
using Plots
Plots.default(linewidth=2)

# Convenience functions
hil(x, k) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)

# Collins toggle switch model
function collins!(du, u, p, t)
    a1, a2, β, γ, i1, i2 = p
    s1, s2 = u
    du[1] = a1 * hil(1 + i2, s2, β) - s1
    du[2] = a2 * hil(1 + i1, s1, γ) - s2
    nothing
end

# Setup the problem
tspan = (0.0, 50.0)
ps = [3.0, 2.5, 4.0, 4.0, 0.0, 0.0] # a1, a2, β, γ, i1, i2
u0 = [0.075, 2.5] # s1, s2
prob = ODEProblem(collins!, u0, tspan, ps)

# Callbacks
affect_i2_on!(integrator) = integrator.p[6] = 10.0
affect_i2_off!(integrator) = integrator.p[6] = 0.0
affect_i1_on!(integrator) = integrator.p[5] = 10.0
affect_i1_off!(integrator) = integrator.p[5] = 0.0
cb_i2_on = PresetTimeCallback(10.0, affect_i2_on!)
cb_i2_off = PresetTimeCallback(20.0, affect_i2_off!)
cb_i1_on = PresetTimeCallback(30.0, affect_i1_on!)
cb_i1_off = PresetTimeCallback(40.0, affect_i1_off!)
cbs = CallbackSet(cb_i2_on, cb_i2_off, cb_i1_on, cb_i1_off)

# Solve the problem
@time sol = solve(prob, callback=cbs)

# Visual
plot(sol, legend=:right, xlabel = "Time", ylabel="Concentration", title="Fig 1.7", labels=["s1" "s2"])

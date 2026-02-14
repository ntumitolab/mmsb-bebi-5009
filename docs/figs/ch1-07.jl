# # Fig 1.7
# For Figures 1.7, 7.13, 7.14, 7.15.
using OrdinaryDiffEq
using DiffEqCallbacks
using ComponentArrays: ComponentArray
using SimpleUnPack
using CairoMakie

# Collins toggle switch model
function collins!(du, u, p, t)
    hil(x, k) = x / (x + k)
    hil(x, k, n) = hil(x^n, k^n)
    @unpack a1, a2, β, γ, i1, i2 = p
    @unpack s1, s2 = u
    du[1] = a1 * hil(1 + i2, s2, β) - s1
    du[2] = a2 * hil(1 + i1, s1, γ) - s2
    nothing
end

# Setup the problem
tend = 50.0
ps = ComponentArray(a1=3.0, a2=2.5, β=4.0, γ=4.0, i1=0.0, i2=0.0)
u0 = ComponentArray(s1=0.075, s2=2.5)
prob = ODEProblem(collins!, u0, (0.0, tend), ps)

# Callbacks and solve the problem
cbs = let
    affect_i2_on!(integrator) = integrator.p.i2 = 10.0
    affect_i2_off!(integrator) = integrator.p.i2 = 0.0
    affect_i1_on!(integrator) = integrator.p.i1 = 10.0
    affect_i1_off!(integrator) = integrator.p.i1 = 0.0
    cb_i2_on = PresetTimeCallback(10.0, affect_i2_on!)
    cb_i2_off = PresetTimeCallback(20.0, affect_i2_off!)
    cb_i1_on = PresetTimeCallback(30.0, affect_i1_on!)
    cb_i1_off = PresetTimeCallback(40.0, affect_i1_off!)
    cbs = CallbackSet(cb_i2_on, cb_i2_off, cb_i1_on, cb_i1_off)
end

@time sol = solve(prob, Tsit5(), callback=cbs)

# Visualization
ts = 0:0.1:tend
data = Array(sol(ts; idxs=[1, 2]))
fig, ax, sp = series(ts, data; labels=["s1", "s2"], axis=(title="Fig. 1.7", xlabel="Time", ylabel="Concentration"))
axislegend(ax)
fig

@time_imports begin
    using CairoMakie
    using OrdinaryDiffEq
    using DiffEqCallbacks
    import ComponentArrays as CA
    using SimpleUnPack
end

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
tspan = (0.0, 50.0)
ps = CA.ComponentArray(a1=3.0, a2=2.5, β=4.0, γ=4.0, i1=0.0, i2=0.0) # a1, a2, β, γ, i1, i2
u0 = CA.ComponentArray(s1=0.075, s2=2.5) # s1, s2
prob = ODEProblem(collins!, u0, tspan, ps)

# Callbacks and solve the problem
affect_i2_on!(integrator) = integrator.p.i2  = 10.0
affect_i2_off!(integrator) = integrator.p.i2 = 0.0
affect_i1_on!(integrator) = integrator.p.i1 = 10.0
affect_i1_off!(integrator) = integrator.p.i1 = 0.0
cb_i2_on = PresetTimeCallback(10.0, affect_i2_on!)
cb_i2_off = PresetTimeCallback(20.0, affect_i2_off!)
cb_i1_on = PresetTimeCallback(30.0, affect_i1_on!)
cb_i1_off = PresetTimeCallback(40.0, affect_i1_off!)
cbs = CallbackSet(cb_i2_on, cb_i2_off, cb_i1_on, cb_i1_off)

@time sol = solve(prob, Tsit5(), callback=cbs)

# Visualization
f = Figure(size = (800, 500))
ax = Axis(f[1, 1])
lines!(ax, sol, idxs=1, color=Cycled(1))
lines!(ax, sol, idxs=2, color=Cycled(2))
f

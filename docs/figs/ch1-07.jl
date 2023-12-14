#===
# Fig 1.07

Collins toggle switch

For Figures 1.7, 7.13, 7.14, 7.15
===#

using DifferentialEquations
using LabelledArrays
using SimpleUnPack
using Plots
Plots.default(linewidth=2)

# PNG output in Literate.jl
PNG(fig) = display("image/png", fig)

# Convenience functions
hil(x, k) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)
exprel(x) = x / expm1(x)

# Collins switch model
function collins!(D, u, p, t)
    @unpack a1, a2, β, γ, i1, i2 = p
    @unpack s1, s2 = u
    D.s1 = a1 * hil(1 + i2, s2, β) - s1
    D.s2 = a2 * hil(1 + i1, s1, γ) - s2
    return nothing
end

# events
on_10!(integrator) = integrator.p.i2 = 10.
on_20!(integrator) = integrator.p.i2 = 0.
on_30!(integrator) = integrator.p.i1 = 10.
on_40!(integrator) = integrator.p.i1 = 0.

events = CallbackSet(
    PresetTimeCallback(10., on_10!),
    PresetTimeCallback(20., on_20!),
    PresetTimeCallback(30., on_30!),
    PresetTimeCallback(40., on_40!),
)

# Prepare
ps = LVector(a1=3.0, a2=2.5, β=4.0, γ=4.0, i1=0.0, i2=0.0)
u0 = LVector(s1=0.075, s2=2.5)
tend = 50.0

# solve and visualize
prob = ODEProblem(collins!, u0, tend, ps)
sol = solve(prob, callback=events)

fig = plot(sol, legend=:right, xlabel = "Time", ylabel="Concentration", title="Fig 1.7")

fig |> PNG

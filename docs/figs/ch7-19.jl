# # Fig 7.19
# Circadian rhythm model
using OrdinaryDiffEq
using ComponentArrays: ComponentArray
using SimpleUnPack
using CairoMakie
using Peaks

#---
function model719!(D, u, p, t)
    hil(x, k=one(x)) = x / (x + k)
    hil(x, k, n) = hil(x^n, k^n)
    @unpack vs, vm, vd, ks, kt1, kt2, v1, v2, v3, v4, k1, k2, k3, k4, ki, km1, kd, n = p
    @unpack M, P0, P1, P2, PN = u
    rM = vs * hil(ki, PN, n) - vm * hil(M, km1)
    rP01 = v1 * hil(P0, k1) - v2 * hil(P1, k2)
    rP12 = v3 * hil(P1, k3) - v4 * hil(P2, k4)
    rP2N = k1 * P2 - k2 * PN
    rP2 = vd * hil(P2, kd)
    D.M = rM
    D.P0 = ks * M - rP01
    D.P1 = rP01 - rP12
    D.P2 = rP12 - rP2N - rP2
    D.PN = rP2N
    nothing
end

#---
ps719 = ComponentArray(
    vs = 0.76,
    vm = 0.65,
    vd = 0.95,
    ks = 0.38,
    kt1 = 1.9,
    kt2 = 1.3,
    v1 = 3.2,
    v2 = 1.58,
    v3 = 5.0,
    v4 = 2.5,
    k1 = 1.0,
    k2 = 1.0,
    k3 = 2.0,
    k4 = 2.0,
    ki = 1.0,
    km1 = 0.5,
    kd = 0.2,
    n = 4
)

ics719 = ComponentArray(
    M = 1.0,
    P0 = 1.0,
    P1 = 0.0,
    P2 = 0.0,
    PN = 0.0
)

tspan = (-50.0, 200.0)
prob719 = ODEProblem(model719!, ics719, tspan, ps719)

#---
@time sol719 = solve(prob719, KenCarp47())

#---
_total_P(sol) = sol.P0 .+ sol.P1 .+ sol.P2 .+ sol.PN
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time",
    ylabel = "Concentration",
    title = "Fig 7.19 (A)"
)
lines!(ax, 0..tspan[2], t-> sol719(t).M, label = "M")
lines!(ax, 0..tspan[2], t-> sol719(t).PN, label = "Nuclear PER")
lines!(ax, 0..tspan[2], t->_total_P(sol719(t)), label = "Total PER")
axislegend(ax, position = :rt)
fig

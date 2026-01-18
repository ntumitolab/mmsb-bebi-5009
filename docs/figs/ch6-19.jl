# # Fig 6.19
# Sine wave response of g-protein signalling pathway
using OrdinaryDiffEq
using ComponentArrays: ComponentArray
using SimpleUnPack
using CairoMakie

#---
_l619(p, t) = p.lt + (p.lt / p.l_AMP) * cospi(2t / p.l_per)

function model619!(D, u, p, t)
    @unpack kRL, kRLm, kGa, kGd0, kG1, Rtotal, Gtotal, lt, l_per, l_AMP = p
    @unpack RL, Ga, Gd = u
    R = Rtotal - RL
    G = Gtotal - Ga - Gd
    Gbg = Ga + Gd
    L = _l619(p, t)
    v1 = kRL * R * L - kRLm * RL
    v2 = kGa * G * RL
    v3 = kGd0 * Ga
    v4 = kG1 * Gd * Gbg
    D.RL = v1
    D.Ga = v2 - v3
    D.Gd = v3 - v4
    nothing
end

#---
ps619 = ComponentArray(
    kRL = 2e6,
    kRLm = 0.01,
    kGa = 1e-5,
    kGd0 = 0.11,
    kG1 = 1.0,
    Rtotal = 4e3,
    Gtotal = 1e4,
    lt = 1e-9,
    l_per = 200.0,
    l_AMP = 5.0
)

ics619 = ComponentArray(
    RL = 0.0,
    Ga = 0.0,
    Gd = 0.0,
)

tend = 1000.0
prob619 = ODEProblem(model619!, ics619, tend, ps619)
#---
@time sol619 = solve(prob619, KenCarp47())

#---
fig = Figure()
ax1 = Axis(fig[1, 1], xlabel="Time", ylabel="Concentration", title="Fig 6.19 A")
lines!(ax1, 0..tend, t -> sol619(t).Ga, label="Ga")
lines!(ax1, 0..tend, t -> _l619(ps619, t) * 1E12, label=L"L \cdot 10^{12}")
axislegend(ax1, position=:rb)
fig

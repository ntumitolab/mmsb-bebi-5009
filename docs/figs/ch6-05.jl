#===
# Fig 6.5

Model of G-protein signalling pathway
===#
using ComponentArrays: ComponentArray
using SimpleUnPack
using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
using CairoMakie

#---
function model605!(D, u, p, t)
    @unpack kRL, kRLm, kGa, kGd0, kG1, L, Rtotal, Gtotal = p
    @unpack RL, Ga, Gd = u
    R = Rtotal - RL
    G = Gtotal - Ga - Gd
    Gbg = Ga + Gd
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
ps605 = ComponentArray(
    kRL = 2e6,
    kRLm = 0.01,
    kGa = 1e-5,
    kGd0 = 0.11,
    kG1 = 1.0,
    Rtotal = 4e3,
    Gtotal = 1e4,
    L = 0.0
)

u0605 = ComponentArray(
    RL = 0.0,
    Ga = 0.0,
    Gd = 0.0,
)

# Events
affect_L1!(integrator) = integrator.p.L = 1e-9
affect_L2!(integrator) = integrator.p.L = 0.0
event_L1 = PresetTimeCallback([200.0], affect_L1!)
event_L2 = PresetTimeCallback([800.0], affect_L2!)
cbs = CallbackSet(event_L1, event_L2)

# ## Fig 6.5 A
tend = 1200.0
prob605 = ODEProblem(model605!, u0605, tend, ps605)

#---
@time sol = solve(prob605, KenCarp47(), callback=cbs)

#---
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time", ylabel="Concentration", title="Fig 6.05 (A)")
lines!(ax, 0 .. tend, t -> sol(t).RL, label="RL")
lines!(ax, 0 .. tend, t -> sol(t).Ga, label="Ga")
axislegend(ax, position=:rc)
fig

# ## Fig 6.5 B
lrange = range(0, 20 * 1e-9, 101)

prob_func = (prob, i, repeat) -> remake(prob, p=ComponentArray(ps605; L=lrange[i]))
prob605b = SteadyStateProblem(model605!, u0605, ps605)
trajectories = length(lrange)
alg = DynamicSS(KenCarp47())
eprob = EnsembleProblem(prob605b; prob_func)
@time sim = solve(eprob, alg; trajectories);

#---
ga = map(s->s.u.Ga, sim)
rl = map(s->s.u.RL, sim)
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel = "Ligand (nM)",
    ylabel = "Steady-state abundance",
    title = "Fig. 6.5 (B)"
)

lines!(ax, lrange .* 1e9, ga, label="Ga")
lines!(ax, lrange .* 1e9, rl, label="RL")
axislegend(ax, position=:rc)
fig

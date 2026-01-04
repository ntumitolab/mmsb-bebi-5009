#===
# Fig 6.3
Two component pathway
===#
using ComponentArrays
using SimpleUnPack
using OrdinaryDiffEq
using SteadyStateDiffEq
using Plots
Plots.default(linewidth=2)

#---
function model603!(D, u, p, t)
    @unpack k1, km1, k2, k3, L, Rtot, Ptot = p
    @unpack RL, Ps = u
    P = Ptot - Ps
    R = Rtot - RL
    v1 = k1 * L * R - km1 * RL
    v2 = k2 * RL * P
    v3 = k3 * Ps
    D.RL = v1
    D.Ps = v2 - v3
    nothing
end

#---
ps603 = ComponentArray(
    k1 = 5.0,
    km1 = 1.0,
    k2 = 6.0,
    k3 = 2.0,
    L = 0.0,
    Rtot = 3.0,
    Ptot = 8.0
)

u0603 = ComponentArray(
    RL = 0.0,
    Ps = 0.0
)

# Events
affect_L1!(integrator) = integrator.p.L = 3.0
affect_L2!(integrator) = integrator.p.L = 0.0
event_L1 = PresetTimeCallback([1.0], affect_L1!)
event_L2 = PresetTimeCallback([3.0], affect_L2!)
cbs = CallbackSet(event_L1, event_L2)

# ## Fig. 6.3 A
tspan = (0., 10.)
prob603a = ODEProblem(model603!, u0603, tspan, ps603)
@time sol603a = solve(prob603a, KenCarp47(), callback=cbs)

#---
fig = plot(sol603a, idxs=[1, 2], labels= ["RL" "P*"])
plot!(fig, t -> 3 * (1<=t<=3), label="Ligand", line=(:black, :dash), linealpha=0.7)
plot!(fig, title="Fig. 6.3 (A)", xlabel="Time", ylabel="Concentration")

# ## Fig 6.3 B
lrange = 0:0.01:1

prob_func = (prob, i, repeat) -> remake(prob, p=ComponentArray(ps603; L=lrange[i]))
prob603b = SteadyStateProblem(model603!, u0603, ps603)
trajectories = length(lrange)
alg = DynamicSS(KenCarp47())
eprob = EnsembleProblem(prob603b; prob_func)
@time sim = solve(eprob, alg; trajectories);

pstar = map(s->s.u.Ps, sim)
rl = map(s->s.u.RL, sim)
plot(lrange, [pstar rl], label=["P*" "RL"], title="Fig. 6.3 (B)",
xlabel="Ligand", ylabel="Steady-state concentration", xlims=(0, 1), ylims=(0, 8))

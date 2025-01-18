#===
# Fig 6.5

Model of G-protein signalling pathway
===#
using ModelingToolkit
using Catalyst
using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
using Plots
Plots.default(linewidth=2)

#---
rn = @reaction_network begin
    @parameters L(t)
    (kRL * L, kRLm), R <--> RL
    kGa, G + RL --> Ga + Gbg + RL
    kGd0, Ga --> Gd
    kG1, Gd + Gbg --> G
end

#---
setdefaults!(rn, [
    :kRL => 2e6,
    :kRLm => 0.01,
    :kGa => 1e-5,
    :kGd0 => 0.11,
    :kG1 => 1,
    :R => 4e3,
    :RL => 0.,
    :G => 1e4,
    :Ga => 0.,
    :Gd => 0.,
    :Gbg => 0.,
    :L => 0.
])

@unpack L = rn
discrete_events = [[200] => [L~1e-9], [800] => [L~0.0]]
osys = convert(ODESystem, rn; discrete_events, remove_conserved = true) |> complete

# ## Fig 6.5 A
tend = 1200.0
prob = ODEProblem(osys, [], tend)

sol = solve(prob, Rodas5P(), abstol=1e-8, reltol=1e-8)

@unpack RL, Ga = osys
plot(sol, idxs=[RL, Ga], title="Fig 6.05 (A)", xlabel="Time", ylabel="Abundance")

# ## Fig 6.5 B
lrange = range(0, 20 * 1e-9, 101)

prob_func = (prob, i, repeat) -> remake(prob, p=[L => lrange[i]])
prob = SteadyStateProblem(osys, [], [])
trajectories = length(lrange)
alg = DynamicSS(Rodas5())
eprob = EnsembleProblem(prob; prob_func)
sim = solve(eprob, alg; trajectories, abstol=1e-10, reltol=1e-10)

ga = map(s->s[Ga], sim)
rl = map(s->s[RL], sim)
plot(lrange .* 1e9, [ga rl], label=["Ga" "RL"], title="Fig. 6.5 (B)",
xlabel="Ligand (nM)", ylabel="Steady-state abundance", xlims=(0, 20), ylims=(0, 3500))

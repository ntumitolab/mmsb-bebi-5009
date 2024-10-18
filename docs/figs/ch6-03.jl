#===
# Fig 6.3
Two component pathway
===#
using ModelingToolkit
using Catalyst
using OrdinaryDiffEq
using Plots
Plots.default(linewidth=2)

#---
rn = @reaction_network begin
    (k1 * L, km1), R <--> RL
    k2, P + RL --> Ps + RL
    k3, Ps--> P
end

#---
setdefaults!(rn, [
    :R => 3.,
    :RL => 0.,
    :P => 8.0,
    :Ps => 0.,
    :L => 0.,
    :k1 => 5.,
    :km1 => 1.,
    :k2 => 6.,
    :k3 => 2.,
])

@independent_variables t
@unpack L = rn
discrete_events = [[1.0] => [L~3.0], [3.0] => [L~0.0]]

osys = convert(ODESystem, rn; discrete_events, remove_conserved = true) |> structural_simplify

#---
observed(osys)

# ## Fig. 6.3 A
tspan = (0., 10.)
prob = ODEProblem(osys, [], tspan, [])
sol = solve(prob)

#---
@unpack RL, Ps = osys
fig = plot(sol, idxs=[RL, Ps], labels= ["RL" "P*"])
plot!(fig, t -> 3 * (1<=t<=3), label="Ligand", line=(:black, :dash), linealpha=0.7)
plot!(fig, title="Fig. 6.3 (A)", xlabel="Time", ylabel="Concentration")

# ## Fig 6.3 B
@unpack RL, Ps = osys
lrange = 0:0.01:1

prob_func = function(prob, i, repeat)
    remake(prob, p=[L => lrange[i]])
end

prob = ODEProblem(osys, [], Inf, [])
callback = TerminateSteadyState()
trajectories = length(lrange)
alg = Rodas5()
eprob = EnsembleProblem(prob; prob_func)
sim = solve(eprob, alg; save_everystep=false, trajectories, callback)

pstar = map(s->s[Ps][end], sim)
rl = map(s->s[RL][end], sim)
plot(lrange, [pstar rl], label=["P*" "RL"], title="Fig. 6.3 (B)",
xlabel="Ligand", ylabel="Steady-state concentration", xlims=(0, 1), ylims=(0, 8))

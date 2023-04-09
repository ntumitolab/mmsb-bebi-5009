#===
# Fig 6.5

Model of G-protein signalling pathway
===#

using ModelingToolkit
using Catalyst
using DifferentialEquations
using Plots
Plots.default(linewidth=2)

#---

rn = @reaction_network begin
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

osys = convert(ODESystem, rn; remove_conserved = true)

#---
observed(osys)
#---
equations(osys)

#---
@unpack L = osys
idx = findfirst(isequal(L), parameters(osys))

# ## Fig. 6.5 A
cb1 = PresetTimeCallback([200.0], i -> begin i.p[idx] = 1e-9; set_proposed_dt!(i, 0.01) end)
cb2 = PresetTimeCallback([800.0], i -> begin i.p[idx] = 0.0; set_proposed_dt!(i, 0.01) end)
cbs = CallbackSet(cb1, cb2)
prob = ODEProblem(osys, [], (0., 1200.))

#---
sol = solve(prob, Rodas5(), callback=cbs, abstol=1e-8, reltol=1e-8, saveat=0:10:1200)
@unpack RL, Ga = osys

fig = plot(sol, idxs=[RL, Ga], title="Fig 6.05 (A)", xlabel="Time", ylabel="Abundance")

# ## Fig 6.3 (B)

lrange = range(0, 20 * 1e-9, 101)

prob_func = function(prob, i, repeat)
    p = copy(prob.p)
    p[idx] = lrange[i]
    remake(prob, p=p)
end

#---
prob = SteadyStateProblem(osys, [])
#---
eprob = EnsembleProblem(prob; prob_func)
#---
sim = solve(eprob, DynamicSS(Rodas5()); trajectories=length(lrange));
#---

ga = map(s->s[Ga], sim)
rl = map(s->s[RL], sim)
plot(lrange .* 1e9, [ga rl], label=["Ga" "RL"], title="Fig. 6.5 (B)",
xlabel="Ligand (nM)", ylabel="Steady-state abundance", xlims=(0, 20), ylims=(0, 3500))

# ## Runtime information

import InteractiveUtils
InteractiveUtils.versioninfo()

#---

import Pkg
Pkg.status()

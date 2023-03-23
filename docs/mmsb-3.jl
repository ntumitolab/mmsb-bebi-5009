# # Chapter 3
# ## Figure 3.03 Michaelis-Menten kinetics

using DifferentialEquations
using Catalyst
using Plots
Plots.default(linewidth=2)

rn303 = @reaction_network begin
    (k1, km1), S + E <--> ES
    k2, ES --> E + P
end

#---
u0 = [:S=>5., :ES=>0., :E=>1., :P=>0.]
ps = [:k1 => 30., :km1 => 1., :k2 => 10.]
tend = 1.0

prob = ODEProblem(rn303, u0, tend, ps)
sol = solve(prob)
plot(sol, xlabel="Time (AU)", ylabel="Concentration (AU)", legend=:right)

#---

rn303mm = @reaction_network begin
    mm(S, k2 * ET, (km1 + k2) / k1), S ⇒ P ## \Rightarrow
end

u0 = [:S=>5., :P=>0.]
ps = [:k1 => 30., :km1 => 1., :k2 => 10., :ET=>1.]
tend = 1.0
probmm = ODEProblem(rn303mm, u0, tend, ps)
solmm = solve(probmm)

@unpack S, P = rn303
plot(sol, idxs=[S, P], line=(:dash), label=["S (full)" "P (full)"])
plot!(solmm, idxs=[S, P], label=["S (MM)" "P (MM)"])
plot!(
    title="Fig. 3.03",
    xlabel="Time (AU)",
    ylabel="Concentration (AU)",
    xlims=(0., tend),
    ylims=(0., 5.),
    legend=:right)

# ## Fig 3.13 GMA and Michaelis-Menten rate laws

using Plots
Plots.default(linewidth=2)

plot(
    [t -> 2t / (1+t), t -> t^0.4], 0., 4.,
    label = ["MM" "GMA"], title = "Fig 3.13",
    xlabel= "Substrate concentration (AU)",
    ylabel="Reaction rate (AU)"
)

# ## Runtime information

import InteractiveUtils
InteractiveUtils.versioninfo()

#---

import Pkg
Pkg.status()
